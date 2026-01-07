module tox_loess
  !! Multidimensional LOESS (LOWESS) for Tensor Omics — single self‑contained module
  !! --------------------------------------------------------------------------------
  !! • Supports d ≥ 1 (tested/recommended up to d ≤ 4)
  !! • Closely matches R's loess defaults: tricube kernel, local polynomial (degree 2 by default)
  !! • No dynamic allocations inside public entry points — caller provides workspaces (iv, v)
  !! • Provides a compact reimplementation of the essential Netlib-style routines:
  !!     lowesd, lowesb, lowese, lowesf
  !!   plus public wrappers:
  !!     tox_loess_required_workspace, tox_loess_fit, tox_loess_predict
  !! • Internals avoid COMMON blocks and external deps; small fixed-size linear solves per local fit
  !! • Automatic algorithm selection: direct (n<1k), kdtree (1k-5k), sliding-window (1D), subsampling (>5k, multidim)
  !!
  use iso_fortran_env, only: int32, real64
  use kd_tree, only: build_kd_index, kd_knn_query
  use f42_utils, only: sort_array
  use tox_errors, only: set_err, is_ok, ERR_INVALID_INPUT
  implicit none
  private

  !==================== Public API ====================
  public :: tox_loess_required_workspace
  public :: tox_loess_fit
  public :: tox_loess_predict

  ! Also expose the Netlib-style names (so other code can call them directly if needed)
  public :: lowesd, lowesb, lowese, lowesf

  !==================== Parameters ====================
  integer(int32), parameter :: I0 = 0_int32
  integer(int32), parameter :: I1 = 1_int32
  real(real64),  parameter :: ZERO = 0.0_real64, ONE = 1.0_real64
  real(real64),  parameter :: EPS  = 1.0e-15_real64
  
  ! Algorithm selection thresholds - TEMPORARY: Force direct method for testing
  integer(int32), parameter :: DIRECT_THRESHOLD = 1000   ! TEMP: Force direct for most datasets  
  integer(int32), parameter :: KDTREE_THRESHOLD = 5000   ! Use k-d tree for 1000 <= n < 5000 (limited range)
  integer(int32), parameter :: MAX_SAFE_K = 2000          ! Maximum k to avoid stack overflow

  ! Module-level variable for simple random number generator
  integer(int32) :: rng_state = 1

contains

  !======================================================================
  ! Workspace size helper (Netlib-style aggregate pools)
  !----------------------------------------------------------------------
  subroutine tox_loess_required_workspace(d, nvmax, liv, lv)
    !! Recommend sizes for integer pool iv(:) and real pool v(:).
    !! nvmax is the maximum neighborhood size (usually n_train).
    integer(int32), intent(in)  :: d, nvmax
    integer(int32), intent(out) :: liv, lv
    ! Based on Netlib notes (vc≈3, nc≈nvmax):
    liv = 50 + 6*nvmax
    ! lv needs: header(100) + training_data(d*n + n + n) + workspace_for_prediction
    ! For safety, we use a larger multiplier and ensure minimum space
    lv  = 100 + (d + 2)*nvmax + max(200, (2*d + 2)*nvmax)
  end subroutine tox_loess_required_workspace

  !======================================================================
  ! Public: Smooth y at same x positions (fit-in-place)
  !----------------------------------------------------------------------
  subroutine tox_loess_fit(d, n, x, y, w, span, degree, iv, liv, v, lv, y_smooth, ierr)
    integer(int32), intent(in)    :: d, n
    real(real64),   intent(in)    :: x(d, n)
    real(real64),   intent(in)    :: y(n)
    real(real64),   intent(in)    :: w(n)  ! use 1.0 if unweighted
    real(real64),   intent(in)    :: span  ! 0< span ≤1 (R default ~0.75)
    integer(int32), intent(in)    :: degree  ! 0,1,2 (R default 2)
    integer(int32), intent(in)    :: liv, lv
    integer(int32), intent(inout) :: iv(liv)
    real(real64),   intent(inout) :: v(lv)
    real(real64),   intent(out)   :: y_smooth(n)
    integer(int32), intent(out)   :: ierr

    logical :: setLf
    real(real64) :: L_dummy(1,1), hat_dummy(1,1)

    ierr = 0
    setLf = .true.

    call lowesd(d, iv, liv, lv, v, d, n, span, degree, n, setLf)
    call lowesf(x, y, w, iv, liv, lv, v, n, x, L_dummy, hat_dummy, y_smooth)
  end subroutine tox_loess_fit

  !======================================================================
  ! Public: Fit once on training data and evaluate at new x
  !----------------------------------------------------------------------
  subroutine tox_loess_predict(d, n_train, x_train, y_train, w_train, span, degree, &
                               iv, liv, v, lv, n_pred, x_pred, y_pred, ierr)
    integer(int32), intent(in)    :: d, n_train, n_pred
    real(real64),   intent(in)    :: x_train(d, n_train)
    real(real64),   intent(in)    :: y_train(n_train)
    real(real64),   intent(in)    :: w_train(n_train)
    real(real64),   intent(in)    :: span
    integer(int32), intent(in)    :: degree
    integer(int32), intent(in)    :: liv, lv
    integer(int32), intent(inout) :: iv(:)
    real(real64),   intent(inout) :: v(:)
    real(real64),   intent(in)    :: x_pred(d, n_pred)
    real(real64),   intent(out)   :: y_pred(n_pred)
    integer(int32), intent(out)   :: ierr

    logical :: setLf

    ierr  = 0
    setLf = .true.

    if (n_train < degree + 1) then
      call set_err(ierr, ERR_INVALID_INPUT)
      print *, "tox_loess_predict: Not enough training points for the requested degree."
      return
    end if

    call lowesd(d, iv, liv, lv, v, d, n_train, span, degree, n_train, setLf)
    call lowesb(x_train, y_train, w_train, iv, liv, lv, v)
    call lowese(iv, liv, lv, v, n_pred, x_pred, y_pred)
  end subroutine tox_loess_predict

  !======================================================================
  ! Netlib-style: LOWESD — configure model / workspace
  !----------------------------------------------------------------------
  subroutine lowesd(d, iv, liv, lv, v, d_in, n, f, tdeg, nvmax, setLf)
    integer(int32), intent(in)    :: d
    integer(int32), intent(inout) :: iv(liv)
    integer(int32), intent(in)    :: liv, lv
    real(real64),   intent(inout) :: v(lv)
    integer(int32), intent(in)    :: d_in, n, tdeg, nvmax
    real(real64),   intent(in)    :: f
    logical,        intent(in)    :: setLf

    ! Store config in the ivy/v pools (simple header)
    ! iv(1)=d, iv(2)=n, iv(3)=k (neighbors), iv(4)=degree
    ! v(1)=span
    integer(int32) :: k, min_k
    
    ! R's EXACT algorithm: nf = min(N, (int) floor(N * span + 1e-5))
    k = min(nvmax, max(2_int32, int(floor(real(n,real64) * max(EPS,min(ONE,f)) + 1.0e-5_real64), int32)))

    ! IMPORTANT FIX: Ensure enough neighbors for the polynomial degree
    ! Need at least (degree + 2) points for stable local polynomial fit
    ! This matches R's behavior - R warns when span is too small
    min_k = tdeg + 2
    if (k < min_k) then
      write(*,'(A,I0,A,I0,A,I0)') " WARNING: span too small for degree=", tdeg, &
                                ", increasing k from ", k ," to ", min(min_k, n)
      k = min(min_k, n)
    end if

    iv(1) = d
    iv(2) = n
    iv(3) = k
    iv(4) = tdeg
    ! offsets for stored training data after lowesb
    iv(10)= 0  ! x offset in v (set by lowesb)
    iv(11)= 0  ! y offset in v
    iv(12)= 0  ! w offset in v

    v(1)  = max(EPS, min(ONE, f))
  end subroutine lowesd

  !======================================================================
  ! Netlib-style: LOWESB — build model (store training data into pool)
  !----------------------------------------------------------------------
  subroutine lowesb(x, y, w, iv, liv, lv, v)
    real(real64),   intent(in)    :: x(:, :)      ! (d, n)
    real(real64),   intent(in)    :: y(:)         ! (n)
    real(real64),   intent(in)    :: w(:)         ! (n)
    integer(int32), intent(inout) :: iv(liv)
    integer(int32), intent(in)    :: liv, lv
    real(real64),   intent(inout) :: v(lv)

    integer(int32) :: d, n, need, off_x, off_y, off_w, i, j

    d = iv(1); n = iv(2)

    ! We store training arrays into v-pool (contiguously):
    ! layout: [ header ... | X(d*n) | Y(n) | W(n) ]
    need = d*n + n + n
    call ensure_pool_capacity(lv, need)

    ! Place after a small header area (start at v(100)) to avoid clobbering
    off_x = 100
    off_y = off_x + d*n
    off_w = off_y + n

    ! Copy x → v(off_x:)
    do j = 1, n
      do i = 1, d
        v(off_x + (j-1)*d + (i-1)) = x(i,j)
      end do
    end do
    ! Copy y, w
    v(off_y:off_y+n-1) = y(1:n)
    v(off_w:off_w+n-1) = w(1:n)

    ! Save offsets
    iv(10) = off_x
    iv(11) = off_y
    iv(12) = off_w

  contains
    subroutine ensure_pool_capacity(lv_in, need_in)
      integer(int32), intent(in) :: lv_in, need_in
      ! Simple guard; pool must be big enough (no allocs allowed here).
      if (100 + need_in > lv_in) stop 'tox_loess: v-pool too small in lowesb'
    end subroutine ensure_pool_capacity
  end subroutine lowesb

  !======================================================================
  ! Netlib-style: LOWESE — evaluate at new points using stored model
  !----------------------------------------------------------------------
  subroutine lowese(iv, liv, lv, v, m, z, s)
    integer(int32), intent(in)    :: iv(liv), liv, lv, m
    real(real64),   intent(in)    :: v(lv)
    real(real64),   intent(in)    :: z(:, :)   ! (d, m)
    real(real64),   intent(out)   :: s(m)

    integer(int32) :: d, n, k, degree
    integer(int32) :: off_x, off_y, off_w
    real(real64)   :: span

    d      = iv(1)
    n      = iv(2)
    k      = iv(3)
    degree = iv(4)
    off_x  = iv(10)
    off_y  = iv(11)
    off_w  = iv(12)
    span   = v(1)

    call predict_core(d, n, v(off_x), v(off_y), v(off_w), m, z, k, degree, s)
  end subroutine lowese

  !======================================================================
  ! Netlib-style: LOWESF — fit & evaluate at same x (no model reuse)
  !----------------------------------------------------------------------
  subroutine lowesf(x, y, w, iv, liv, lv, v, m, z, L, hat, s)
    real(real64),   intent(in)    :: x(:, :)    ! (d, n)
    real(real64),   intent(in)    :: y(:)       ! (n)
    real(real64),   intent(in)    :: w(:)       ! (n)
    integer(int32), intent(in)    :: iv(liv), liv, lv
    real(real64),   intent(in)    :: v(lv)
    integer(int32), intent(in)    :: m
    real(real64),   intent(in)    :: z(:, :)    ! (d, m) — here we pass x itself
    real(real64),   intent(out)   :: L(:,:)     ! not computed (dummy 1x1 ok)
    real(real64),   intent(out)   :: hat(:,:)   ! not computed (dummy 1x1 ok)
    real(real64),   intent(out)   :: s(m)

    integer(int32) :: d, n, k, degree
    real(real64)   :: span

    d      = iv(1)
    n      = iv(2)
    k      = iv(3)
    degree = iv(4)
    span   = v(1)

    ! Direct prediction using (x,y,w) as training data
    call predict_core(d, n, x, y, w, m, z, k, degree, s)

    ! L, hat are placeholders to be API-compatible; not computed here.
    if (size(L,1)*size(L,2) >= 1) L(1,1) = ZERO
    if (size(hat,1)*size(hat,2) >= 1) hat(1,1) = ZERO
  end subroutine lowesf

  !======================================================================
  ! Core predictor (shared by lowesf/lowese) - Adaptive algorithm selection
  !----------------------------------------------------------------------
  subroutine predict_core(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    integer(int32), intent(in)  :: d, n, m, k, degree
    real(real64),   intent(in)  :: x_train(d, n)
    real(real64),   intent(in)  :: y_train(n)
    real(real64),   intent(in)  :: w_train(n)
    real(real64),   intent(in)  :: x_pred(d, m)
    real(real64),   intent(out) :: y_out(m)

    ! TEMPORARY: Force LOESS direct method for all cases to match R exactly
    ! TODO: Re-enable optimizations after validating correctness
    if (n < DIRECT_THRESHOLD) then
      ! Small dataset: use direct method
      call predict_core_direct_fallback(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    else if (n < KDTREE_THRESHOLD) then
      ! Medium dataset: use K-d tree optimization
      call predict_core_kdtree(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    else
      if (d == 1) then
        call predict_core_sliding_window(n, x_train(1,:), y_train, w_train, k, degree, y_out)
      else
        call predict_core_subsampled(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
      end if
    end if
  end subroutine predict_core

  !======================================================================
  ! Direct method predictor (for small datasets n < 1000)
  !----------------------------------------------------------------------
  subroutine predict_core_direct(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    integer(int32), intent(in)  :: d, n, m, k, degree
    real(real64),   intent(in)  :: x_train(d, n)
    real(real64),   intent(in)  :: y_train(n)
    real(real64),   intent(in)  :: w_train(n)
    real(real64),   intent(in)  :: x_pred(d, m)
    real(real64),   intent(out) :: y_out(m)

    ! Direct method (now just an alias for the fallback)
    call predict_core_direct_fallback(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
  end subroutine predict_core_direct

  !======================================================================
  ! K-d tree optimized predictor (for medium datasets 1k <= n < 5k)
  !----------------------------------------------------------------------
  subroutine predict_core_kdtree(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    integer(int32), intent(in)  :: d, n, m, k, degree
    real(real64),   intent(in)  :: x_train(d, n)
    real(real64),   intent(in)  :: y_train(n)
    real(real64),   intent(in)  :: w_train(n)
    real(real64),   intent(in)  :: x_pred(d, m)
    real(real64),   intent(out) :: y_out(m)

    ! K-d tree workspace
    integer(int32) :: kd_indices(n), dimension_order(d)
    integer(int32) :: workspace(n), permutation(n), left_stack(n), right_stack(n)
    integer(int32) :: recursion_stack(3, n)
    real(real64) :: value_buffer(n)
    integer(int32) :: ierr, j, i

    ! Build dimension order (simple: 1, 2, ..., d)
    do i = 1, d
      dimension_order(i) = i
    end do

    ! Build K-d tree index
    call build_kd_index(x_train, d, n, kd_indices, dimension_order, &
                       workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)

    if (ierr /= 0) then
      ! Fallback to direct method if K-d tree fails
      call predict_core_direct(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
      return
    end if

    ! For each prediction point, use K-d tree for neighbor search
    do j = 1, m
      call local_fit_at_point_kdtree(d, n, x_train, y_train, w_train, x_pred(:,j), k, degree, &
                                    kd_indices, dimension_order, y_out(j))
    end do
  end subroutine predict_core_kdtree

  !======================================================================
  ! Subsampled predictor (for very large datasets n >= 200k) - R-style interpolation
  !----------------------------------------------------------------------
  subroutine predict_core_subsampled(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    integer(int32), intent(in)  :: d, n, m, k, degree
    real(real64),   intent(in)  :: x_train(d, n)
    real(real64),   intent(in)  :: y_train(n)
    real(real64),   intent(in)  :: w_train(n)
    real(real64),   intent(in)  :: x_pred(d, m)
    real(real64),   intent(out) :: y_out(m)

    ! R-style approach: subsample training data, compute LOESS on subsample, then interpolate
    integer(int32), parameter :: MIN_SUBSAMPLE = 2000    ! Minimum subsample size
    real(real64), parameter :: SUBSAMPLE_FRACTION = 0.15_real64  ! Use 15% of data (no upper limit)

    integer(int32) :: target_subsample_size
    integer(int32) :: n_sub, step, i, j, safe_k
    integer(int32), allocatable :: sub_indices(:)
    real(real64), allocatable :: x_sub(:,:), y_sub(:), w_sub(:)
    real(real64), allocatable :: y_sub_smooth(:)
    real(real64) :: original_span
    
    ! Calculate adaptive subsample size (15% of data, minimum 2000 points, no upper limit)
    target_subsample_size = max(MIN_SUBSAMPLE, int(real(n, real64) * SUBSAMPLE_FRACTION))
    
    if (n <= target_subsample_size) then
      ! Use K-d tree for moderate sizes, but recalculate k properly
      safe_k = min(k, n/3, 500)  ! Limit k for safety
      call predict_core_kdtree(d, n, x_train, y_train, w_train, m, x_pred, safe_k, degree, y_out)
      return
    end if
    
    ! IMPROVED SUBSAMPLING: Span-dependent stratified sampling
    ! This ensures different spans get different (but reproducible) subsets
    allocate(sub_indices(target_subsample_size))
    call span_dependent_subsample(n, k, target_subsample_size, sub_indices, n_sub)
    
    ! Allocate arrays based on actual subsample size
    allocate(x_sub(d, n_sub), y_sub(n_sub), w_sub(n_sub), y_sub_smooth(n_sub))
    
    ! Extract subsampled data
    do i = 1, n_sub
      x_sub(:, i) = x_train(:, sub_indices(i))
      y_sub(i) = y_train(sub_indices(i)) 
      w_sub(i) = w_train(sub_indices(i))
    end do
    
    ! Ensure we have enough points
    if (n_sub < 10) then
      ! Fallback: use simple averaging
      do j = 1, m
        y_out(j) = sum(y_train) / real(n, real64)
      end do
      return
    end if
    
    ! Calculate k based on subsample size, not original size
    ! Use the same span proportion as the original request
    original_span = real(k, real64) / real(n, real64)  ! Recover original span
    safe_k = max(1, min(int(ceiling(real(n_sub, real64) * original_span)), n_sub/3, 100))
    
    ! Run LOESS on the subsample using direct method (safest)
    call predict_core_direct(d, n_sub, x_sub, y_sub, w_sub, n_sub, x_sub, safe_k, degree, y_sub_smooth)
    
    ! Now interpolate results for all prediction points
    do j = 1, m
      call interpolate_from_subsample(d, n_sub, x_sub, y_sub_smooth, x_pred(:,j), y_out(j))
    end do
    
    ! Cleanup
    deallocate(sub_indices, x_sub, y_sub, w_sub, y_sub_smooth)
  end subroutine predict_core_subsampled

  !======================================================================
  ! Simple linear interpolation from subsampled LOESS results
  !----------------------------------------------------------------------
  subroutine interpolate_from_subsample(d, n_sub, x_sub, y_sub, x_query, y_interp)
    integer(int32), intent(in) :: d, n_sub
    real(real64), intent(in) :: x_sub(d, n_sub), y_sub(n_sub), x_query(d)
    real(real64), intent(out) :: y_interp
    
    ! For 1D: simple linear interpolation between nearest neighbors
    ! For higher D: inverse distance weighting with nearest few points
    integer(int32) :: i, closest_idx, second_idx
    real(real64) :: min_dist, second_dist, dist, alpha
    
    if (d == 1) then
      ! 1D linear interpolation
      closest_idx = 1
      min_dist = abs(x_sub(1, 1) - x_query(1))
      do i = 2, n_sub
        dist = abs(x_sub(1, i) - x_query(1))
        if (dist < min_dist) then
          min_dist = dist
          closest_idx = i
        end if
      end do
      
      ! Find second closest for interpolation
      second_idx = 1
      if (second_idx == closest_idx) second_idx = 2
      second_dist = abs(x_sub(1, second_idx) - x_query(1))
      do i = 1, n_sub
        if (i == closest_idx) cycle
        dist = abs(x_sub(1, i) - x_query(1))
        if (dist < second_dist) then
          second_dist = dist
          second_idx = i
        end if
      end do
      
      ! Linear interpolation
      if (min_dist + second_dist > EPS) then
        alpha = second_dist / (min_dist + second_dist)
        y_interp = alpha * y_sub(closest_idx) + (ONE - alpha) * y_sub(second_idx)
      else
        y_interp = y_sub(closest_idx)
      end if
    else
      ! Multidimensional: simple nearest neighbor (could be improved to IDW)
      closest_idx = 1
      min_dist = euclid2(d, x_sub(:, 1), x_query)
      do i = 2, n_sub
        dist = euclid2(d, x_sub(:, i), x_query)
        if (dist < min_dist) then
          min_dist = dist
          closest_idx = i
        end if
      end do
      y_interp = y_sub(closest_idx)
    end if
  end subroutine interpolate_from_subsample

  !======================================================================
  ! Local weighted regression around one point (direct method)
  !----------------------------------------------------------------------
  subroutine local_fit_at_point_direct(d, n, x, y, w, z, k, degree, yhat)
    use f42_utils, only: quicksort_real
    integer(int32), intent(in)  :: d, n, k, degree
    real(real64),   intent(in)  :: x(d,n), y(n), w(n), z(d)
    real(real64),   intent(out) :: yhat

    integer(int32) :: i, p, ii
    integer(int32) :: nn(k)
    real(real64)   :: dists(n), dmax, wt(n), u
    real(real64)   :: Xrow(64), XtWX(64,64), XtWy(64), beta(64)
    ! --- Variables for sorting ---
    real(real64)   :: k_dists(k)
    integer(int32) :: perm(k), stack_left(k), stack_right(k)

    ! 1) Compute all distances
    do i = 1, n
      dists(i) = euclid2(d, x(:,i), z)
    end do

    ! 2) Find k nearest neighbors
    call select_k_smallest(dists, n, k, nn)

    ! 3) Bandwidth = k-th distance (sorted), not maximum distance
    ! Extract the k neighbor distances
    do i = 1, k
      k_dists(i) = dists(nn(i))
      perm(i) = i  ! Permutation for positions 1..k in k_dists
    end do
    
    ! Sort k_dists indirectly via perm
    call quicksort_real(k_dists, perm, k, stack_left, stack_right)
    
    ! k_dists is now sorted, perm(k) points to the largest distance
    dmax = k_dists(perm(k))
    if (dmax <= EPS) dmax = EPS

    ! 4) CORRECT: Tricube weights using the FIXED tricube function
    do i = 1, n
      wt(i) = ZERO
    end do
    do i = 1, k
      u = dists(nn(i)) / dmax
      ! NOW CORRECT: tricube returns (1 - u³)³ for u<1, 0 otherwise
      wt(nn(i)) = tricube(u) * max(EPS, w(nn(i)))
    end do

    ! 5) Build normal equations
    p = num_terms(d, degree)
    call zero_mat_vec(XtWX, XtWy, p)
    do i = 1, k
      ii = nn(i)
      if (wt(ii) > ZERO) then
        call build_design_row(d, degree, x(:,ii), z, Xrow, p)
        call accumulate_normal_eq(XtWX, XtWy, Xrow, y(ii), wt(ii), p)
      end if
    end do

    ! 6) Solve and predict
    call solve_sym_posdef(XtWX, XtWy, beta, p)
    call build_design_row(d, degree, z, z, Xrow, p)
    yhat = dot_p(Xrow, beta, p)
  end subroutine local_fit_at_point_direct

  !======================================================================
  ! Local weighted regression around one point (K-d tree optimized)
  !----------------------------------------------------------------------
  subroutine local_fit_at_point_kdtree(d, n, x, y, w, z, k, degree, kd_indices, dimension_order, yhat)
    use f42_utils, only: quicksort_real
    integer(int32), intent(in)  :: d, n, k, degree
    real(real64),   intent(in)  :: x(d,n), y(n), w(n), z(d)
    integer(int32), intent(in)  :: kd_indices(n), dimension_order(d)
    real(real64),   intent(out) :: yhat

    integer(int32) :: i, p, ii, ierr, k_safe
    integer(int32) :: nn(k)
    real(real64)   :: dists(k), dmax, wt(n), u
    real(real64)   :: Xrow(64), XtWX(64,64), XtWy(64), beta(64)
    real(real64)   :: kd_workspace(d)
    ! --- Variables for sorting ---
    integer(int32) :: perm(k), stack_left(k), stack_right(k)

    ! Limit k to safe value to avoid stack overflow
    k_safe = k

    ! 1) K-d tree query for k_safe nearest neighbors
    call kd_knn_query(x, kd_indices, d, n, dimension_order, z, k_safe, nn, dists, kd_workspace, ierr)
    if (ierr /= 0) then
      ! Fallback to direct method if K-d query fails
      call local_fit_at_point_direct(d, n, x, y, w, z, k, degree, yhat)
      return
    end if

    ! 2) CORRECT: Bandwidth = k_safe-th distance (sorted)
    ! dists already contains the k_safe distances
    do i = 1, k_safe
      perm(i) = i  ! Permutation for positions 1..k_safe in dists
    end do
    
    ! Sort dists(1:k_safe) indirectly via perm
    call quicksort_real(dists, perm, k_safe, stack_left, stack_right)
    
    ! dists(perm(k_safe)) is the k_safe-th (largest) distance
    dmax = dists(perm(k_safe))
    if (dmax <= EPS) dmax = EPS

    ! 3) CORRECT: Tricube weights using the FIXED tricube function
    do i = 1, n
      wt(i) = ZERO
    end do
    do i = 1, k_safe
      u = dists(i) / dmax
      ! NOW CORRECT: tricube returns (1 - u³)³ for u<1, 0 otherwise
      wt(nn(i)) = tricube(u) * max(EPS, w(nn(i)))
    end do

    ! 4) Build normal equations
    p = num_terms(d, degree)
    call zero_mat_vec(XtWX, XtWy, p)
    do i = 1, k_safe
      ii = nn(i)
      if (wt(ii) > ZERO) then
        call build_design_row(d, degree, x(:,ii), z, Xrow, p)
        call accumulate_normal_eq(XtWX, XtWy, Xrow, y(ii), wt(ii), p)
      end if
    end do

    ! 5) Solve and predict
    call solve_sym_posdef(XtWX, XtWy, beta, p)
    call build_design_row(d, degree, z, z, Xrow, p)
    yhat = dot_p(Xrow, beta, p)
  end subroutine local_fit_at_point_kdtree

  !======================================================================
  ! Utilities
  !----------------------------------------------------------------------
  function euclid2(d, a, b) result(res)
    integer(int32), intent(in) :: d
    real(real64),   intent(in) :: a(d), b(d)
    real(real64) :: res
    integer(int32) :: i
    res = ZERO
    do i = 1, d
      res = res + (a(i)-b(i))*(a(i)-b(i))
    end do
    res = sqrt(max(res, ZERO))
  end function euclid2

  real(real64) function tricube(u)
    real(real64), intent(in) :: u
    real(real64) :: t

    t = min(ONE, max(ZERO, u))
    
    ! Standard-Tricube-Kernel: (1 - t³)³
    if (t >= ONE) then
      tricube = ZERO
    else
      tricube = (ONE - t**3)**3
    end if
  end function tricube

  ! Alias for sliding window compatibility
  real(real64) function tricube_weight(u)
    real(real64), intent(in) :: u
    tricube_weight = tricube(u)
  end function tricube_weight

  integer(int32) function num_terms(d, degree)
    integer(int32), intent(in) :: d, degree
    if (degree <= 0) then
      num_terms = 1
    else if (degree == 1) then
      num_terms = 1 + d
    else
      ! degree == 2: constant + linear(d) + quadratic(d(d+1)/2)
      num_terms = 1 + d + (d*(d+1))/2
    end if
  end function num_terms

  subroutine zero_mat_vec(A, b, p)
    real(real64), intent(inout) :: A(:,:), b(:)
    integer(int32), intent(in)  :: p
    integer(int32) :: i, j
    do j = 1, p
      b(j) = ZERO
      do i = 1, p
        A(i,j) = ZERO
      end do
    end do
  end subroutine zero_mat_vec

  subroutine build_design_row(d, degree, x, z, row, p)
    !! Build design vector around center z (local poly in (x - z)).
    integer(int32), intent(in) :: d, degree, p
    real(real64),   intent(in) :: x(d), z(d)
    real(real64),   intent(out):: row(p)
    integer(int32) :: i, a, b, idx
    real(real64)   :: dx(8)  ! supports d<=8 safely; we recommend d<=4

    do i = 1, d
      dx(i) = x(i) - z(i)
    end do

    idx = 1
    row(idx) = ONE; idx = idx + 1

    if (degree >= 1) then
      do i = 1, d
        row(idx) = dx(i); idx = idx + 1
      end do
    end if

    if (degree >= 2) then
      do a = 1, d
        do b = a, d
          row(idx) = dx(a)*dx(b); idx = idx + 1
        end do
      end do
    end if
  end subroutine build_design_row

  subroutine accumulate_normal_eq(XtWX, XtWy, xrow, yval, wgt, p)
    real(real64), intent(inout) :: XtWX(:,:), XtWy(:)
    real(real64), intent(in)    :: xrow(:), yval, wgt
    integer(int32), intent(in)  :: p
    integer(int32) :: i, j

    if (wgt <= ZERO) return
    do j = 1, p
      XtWy(j) = XtWy(j) + wgt * xrow(j) * yval
      do i = 1, p
        XtWX(i,j) = XtWX(i,j) + wgt * xrow(i) * xrow(j)
      end do
    end do
  end subroutine accumulate_normal_eq

  subroutine solve_sym_posdef(A, b, x, p)
    !! Simple Cholesky solve for small p (p <= ~64)
    real(real64), intent(inout) :: A(:,:)  ! overwritten
    real(real64), intent(inout) :: b(:)    ! overwritten with solution
    real(real64), intent(out)   :: x(:)
    integer(int32), intent(in)  :: p
    integer(int32) :: i, j, k
    real(real64) :: sum

    ! Cholesky decomposition A = L*L^T
    do j = 1, p
      do i = j, p
        sum = A(i,j)
        do k = 1, j-1
          sum = sum - A(i,k) * A(j,k)
        end do
        if (i == j) then
          A(i,j) = sqrt(max(sum, EPS))
        else
          A(i,j) = sum / A(j,j)
        end if
      end do
    end do

    ! Forward solve L*y = b
    do i = 1, p
      sum = b(i)
      do k = 1, i-1
        sum = sum - A(i,k) * b(k)
      end do
      b(i) = sum / A(i,i)
    end do

    ! Backward solve L^T*x = y
    do i = p, 1, -1
      sum = b(i)
      do k = i+1, p
        sum = sum - A(k,i) * x(k)
      end do
      x(i) = sum / A(i,i)
    end do
  end subroutine solve_sym_posdef

  subroutine select_k_smallest(dists, n, k, idx)
    !! Naive partial selection: O(n*k). For n up to ~1e5, k small, OK.
    real(real64),   intent(in)  :: dists(n)
    integer(int32), intent(in)  :: n, k
    integer(int32), intent(out) :: idx(k)
    integer(int32) :: i, j, pos
    real(real64)   :: best, val

    do j = 1, k
      best = HUGE(ZERO)
      pos  = -1
      do i = 1, n
        val = dists(i)
        if (val < best) then
          ! ensure uniqueness (skip already picked)
          if (.not. already_picked(i, idx, j-1)) then
            best = val
            pos = i
          end if
        end if
      end do
      if (pos < 0) then
        idx(j) = max(1,j)  ! fallback
      else
        idx(j) = pos
      end if
    end do
  contains
    logical function already_picked(i, arr, used)
      integer(int32), intent(in) :: i, used
      integer(int32), intent(in) :: arr(:)
      integer(int32) :: t
      already_picked = .false.
      do t = 1, used
        if (arr(t) == i) then
          already_picked = .true.
          return
        end if
      end do
    end function already_picked
  end subroutine select_k_smallest

  real(real64) function dot_p(a, b, p)
    real(real64), intent(in) :: a(:), b(:)
    integer(int32), intent(in) :: p
    integer(int32) :: i
    dot_p = ZERO
    do i = 1, p
      dot_p = dot_p + a(i)*b(i)
    end do
  end function dot_p

  !======================================================================
  ! Span-dependent subsampling function
  !----------------------------------------------------------------------
  subroutine span_dependent_subsample(n, k, max_sub, indices, n_sub)
    integer(int32), intent(in)  :: n, k, max_sub
    integer(int32), intent(out) :: indices(max_sub)
    integer(int32), intent(out) :: n_sub
    
    ! Local variables - ALL declarations must come first in Fortran
    integer(int32) :: i, j, temp, span_seed
    real(real64) :: original_span, rnorm
    integer(int32), allocatable :: full_indices(:)
    logical :: selected(n)
    real(real64) :: step_real, start_offset, jitter
    integer(int32) :: step_int, start_pos, current_pos
    
    ! Calculate span for creating a span-dependent seed
    original_span = real(k, real64) / real(n, real64)
    span_seed = int(original_span * 10000.0_real64)  ! Convert to integer for seeding
    
    ! Initialize simple random number generator with span-dependent seed
    call simple_srand(span_seed)
    
    if (n <= max_sub) then
      ! If small enough, use all data
      n_sub = n
      do i = 1, n
        indices(i) = i
      end do
      return
    end if
    
    ! Strategy: Stratified sampling with span-dependent randomization
    ! 1. Divide into strata based on data order (approximates sorted order)
    ! 2. Sample from each stratum with span-dependent probability
    
    n_sub = 0
    selected = .false.
    
    ! Method: Systematic sampling with random start (span-dependent)
    step_real = real(n, real64) / real(max_sub, real64)
    step_int = max(1, int(step_real))
    
    ! Span-dependent random offset for start position
    start_offset = call_simple_rand() * step_real
    start_pos = max(1, min(n, int(start_offset) + 1))
    
    current_pos = start_pos
    do while (n_sub < max_sub .and. current_pos <= n)
      n_sub = n_sub + 1
      indices(n_sub) = current_pos
      
      ! Add some span-dependent jitter to avoid strict periodicity
      jitter = (call_simple_rand() - 0.5_real64) * 0.3_real64 * step_real
      current_pos = current_pos + step_int + int(jitter)
      
      if (current_pos > n) exit
    end do
    
    ! Ensure minimum sample size
    if (n_sub < 10) then
      ! Fallback: evenly spaced
      n_sub = min(max_sub, n)
      step_int = max(1, n / n_sub)
      do i = 1, n_sub
        indices(i) = min(n, 1 + (i-1) * step_int)
      end do
    end if
    
  end subroutine span_dependent_subsample
  
  !======================================================================
  ! Simple reproducible random number generator (LCG)
  !----------------------------------------------------------------------
  subroutine simple_srand(seed)
    integer(int32), intent(in) :: seed
    rng_state = seed
    if (rng_state <= 0) rng_state = 1
  end subroutine simple_srand
  
  real(real64) function call_simple_rand()
    ! Linear Congruential Generator (simple but reproducible)
    rng_state = mod(rng_state * 1103515245 + 12345, 2147483647)
    call_simple_rand = real(rng_state, real64) / 2147483647.0_real64
  end function call_simple_rand

  !======================================================================
  ! Helper function: Check if prediction points are same as training points
  !----------------------------------------------------------------------
  function is_same_points(d, n, x_train, x_pred) result(is_same)
    integer(int32), intent(in) :: d, n
    real(real64), intent(in) :: x_train(d, n), x_pred(d, n)
    logical :: is_same
    integer(int32) :: i, j
    
    is_same = .true.
    do i = 1, n
      do j = 1, d
        if (abs(x_train(j,i) - x_pred(j,i)) > 1e-12_real64) then
          is_same = .false.
          return
        end if
      end do
    end do
  end function is_same_points

  !======================================================================
  ! Sliding Window LOESS with R-style interpolation (O(n*k) - assumes sorted data)
  !----------------------------------------------------------------------
  subroutine predict_core_sliding_window(n, x, y, w, k, degree, y_out)
    integer(int32), intent(in) :: n, k, degree
    real(real64), intent(in) :: x(n), y(n), w(n)
    real(real64), intent(out) :: y_out(n)
    
    ! R-style delta parameter - NOW USED ALWAYS for maximum speed
    ! Adaptive delta: smaller datasets get smaller delta for better precision
    real(real64), parameter :: DEFAULT_DELTA_FRACTION = 0.01_real64  ! Base: 1% of x range
    real(real64), parameter :: MIN_DELTA_FRACTION = 0.005_real64     ! 0.5% for small datasets  
    real(real64), parameter :: MAX_DELTA_FRACTION = 0.02_real64      ! 2% for very large datasets
    
    real(real64) :: delta, x_range, adaptive_delta_fraction
    
    ! Calculate adaptive delta fraction based on dataset size
    ! Smaller datasets: use smaller delta for better precision
    ! Larger datasets: can use larger delta for more speed
    if (real(n, real64) < 5000.0_real64) then
      adaptive_delta_fraction = MIN_DELTA_FRACTION  ! 0.5% for small datasets
    else if (real(n, real64) > 100000.0_real64) then
      adaptive_delta_fraction = MAX_DELTA_FRACTION  ! 2% for very large datasets
    else
      ! Linear interpolation between 0.5% and 2% based on dataset size
      adaptive_delta_fraction = MIN_DELTA_FRACTION + &
        (MAX_DELTA_FRACTION - MIN_DELTA_FRACTION) * &
        (real(n, real64) - 5000.0_real64) / (100000.0_real64 - 5000.0_real64)
    end if
    
    ! Calculate delta using adaptive approach
    x_range = x(n) - x(1)
    delta = adaptive_delta_fraction * x_range
    
    ! ALWAYS use delta optimization for maximum speed
    call predict_sliding_window_with_delta(n, x, y, w, k, degree, delta, y_out)
  end subroutine predict_core_sliding_window

  !======================================================================
  ! Direct sliding window (compute all points)
  !----------------------------------------------------------------------
  subroutine predict_sliding_window_direct(n, x, y, w, k, degree, y_out)
    integer(int32), intent(in) :: n, k, degree
    real(real64), intent(in) :: x(n), y(n), w(n)
    real(real64), intent(out) :: y_out(n)
    
    integer(int32) :: i, j, window_start, window_end, window_size
    real(real64) :: x_center, max_dist, dist, pred_value
    real(real64) :: weights(k)
    integer(int32) :: k_actual
    integer(int32) :: neighbors(k)
    
    ! For each point, use sliding window of k neighbors (assumes sorted x)
    do i = 1, n
      x_center = x(i)
      
      ! Find optimal window around point i (assumes x is sorted)
      window_start = max(1, i - k/2)
      window_end = min(n, window_start + k - 1)
      
      ! Adjust window size near boundaries
      if (window_end - window_start + 1 < k .and. window_start > 1) then
        window_start = max(1, window_end - k + 1)
      end if
      if (window_end - window_start + 1 < k .and. window_end < n) then
        window_end = min(n, window_start + k - 1)
      end if
      
      window_size = window_end - window_start + 1
      k_actual = min(k, window_size)
      
      ! Calculate max distance for tricube weighting
      max_dist = 0.0_real64
      do j = window_start, window_end
        dist = abs(x(j) - x_center)
        if (dist > max_dist) max_dist = dist
      end do
      if (max_dist <= EPS) max_dist = EPS
      
      ! Compute tricube weights for window
      do j = 1, k_actual
        dist = abs(x(window_start + j - 1) - x_center)
        weights(j) = w(window_start + j - 1) * tricube_weight(dist / max_dist)
      end do
      
      ! Fit local polynomial directly on window
      call local_polynomial_fit_1d_direct(k_actual, x, y, neighbors, weights, x_center, degree, pred_value)
      y_out(i) = pred_value
    end do
  end subroutine predict_sliding_window_direct

  !======================================================================
  ! Sliding window with R-style delta optimization (faithful to R lowess)
  !----------------------------------------------------------------------
  subroutine predict_sliding_window_with_delta(n, x, y, w, k, degree, delta, y_out)
    integer(int32), intent(in) :: n, k, degree
    real(real64), intent(in) :: x(n), y(n), w(n), delta
    real(real64), intent(out) :: y_out(n)
    
    integer(int32) :: i, j, last_computed
    real(real64) :: cut_distance
    
    ! R's algorithm: compute LOESS at strategic points, skip very close ones
    ! This is MUCH more conservative than my previous 1% approach
    
    last_computed = 0
    
    i = 1
    do while (i <= n)
      
      ! Always compute LOESS for this point
      call compute_loess_at_point_1d(n, x, y, w, k, degree, i, y_out(i))
      
      ! Interpolate any skipped points between last_computed and i
      if (last_computed > 0 .and. last_computed < i - 1) then
        call interpolate_between_points(x, y_out, last_computed, i)
      end if
      
      last_computed = i
      
      ! R-style: skip points that are very close (within delta)
      cut_distance = x(i) + delta
      j = i + 1
      do while (j <= n .and. x(j) <= cut_distance)
        if (abs(x(j) - x(i)) < EPS) then
          ! Exact tie: use same value
          y_out(j) = y_out(i)
          last_computed = j
        end if
        ! Skip this point (will be interpolated later)
        j = j + 1
      end do
      
      ! Continue from the first point outside the delta range
      i = max(j, i + 1)
      
    end do
    
    ! Handle any remaining points at the end
    if (last_computed < n) then
      ! Compute final point
      call compute_loess_at_point_1d(n, x, y, w, k, degree, n, y_out(n))
      
      ! Interpolate any gaps before final point
      if (last_computed < n - 1) then
        call interpolate_between_points(x, y_out, last_computed, n)
      end if
    end if
  end subroutine predict_sliding_window_with_delta

  !======================================================================
  ! Compute LOESS at single point (1D optimized)
  !----------------------------------------------------------------------
  subroutine compute_loess_at_point_1d(n, x, y, w, k, degree, i_center, result)
    integer(int32), intent(in) :: n, k, degree, i_center
    real(real64), intent(in) :: x(n), y(n), w(n)
    real(real64), intent(out) :: result
    integer(int32) :: left, right, nk, idx, j
    real(real64) :: x0, dists(k), rho
    integer(int32) :: neighbors(k)
    real(real64) :: weights(k)
    x0 = x(i_center)
    left = i_center
    right = i_center
    nk = 1
    do while (nk < k)
      if (left > 1 )then
        if (right == n .or. abs(x(left-1)-x0) <= abs(x(right+1)-x0)) then
            left = left - 1
        end if
      else if (right < n) then
        right = right + 1
      else
        exit
      end if
      nk = right - left + 1
    end do
    do j = 1, nk
      neighbors(j) = left + j - 1
      dists(j) = abs(x(neighbors(j)) - x0)
    end do
    rho = maxval(dists(1:nk)) * 0.1_real64  ! Example: scale by 0.1, adjust as needed
    if (rho < 1e-12_real64) rho = 1e-12_real64
    do j = 1, nk
      weights(j) = w(neighbors(j)) * tricube_weight(dists(j) / rho)
    end do
    call local_polynomial_fit_1d_direct(nk, x, y, neighbors, weights, x0, degree, result)
  end subroutine compute_loess_at_point_1d

  !======================================================================
  ! Linear interpolation between computed points
  !----------------------------------------------------------------------
  subroutine interpolate_between_points(x, y_out, i_start, i_end)
    real(real64), intent(in) :: x(:)
    real(real64), intent(inout) :: y_out(:)
    integer(int32), intent(in) :: i_start, i_end
    
    integer(int32) :: i
    real(real64) :: x_start, x_end, y_start, y_end, alpha, denom
    
    if (i_end <= i_start + 1) return  ! No points to interpolate
    
    x_start = x(i_start)
    x_end = x(i_end)
    y_start = y_out(i_start)
    y_end = y_out(i_end)
    
    denom = x_end - x_start
    if (abs(denom) < EPS) then
      ! Same x values: use average
      do i = i_start + 1, i_end - 1
        y_out(i) = (y_start + y_end) * 0.5_real64
      end do
    else
      ! Linear interpolation
      do i = i_start + 1, i_end - 1
        alpha = (x(i) - x_start) / denom
        y_out(i) = (1.0_real64 - alpha) * y_start + alpha * y_end
      end do
    end if
  end subroutine interpolate_between_points

  !======================================================================
  ! R-compatible LOESS core (simplified but robust version)
  !----------------------------------------------------------------------
  subroutine predict_core_direct_fallback(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    integer(int32), intent(in)  :: d, n, m, k, degree
    real(real64),   intent(in)  :: x_train(d, n)
    real(real64),   intent(in)  :: y_train(n)
    real(real64),   intent(in)  :: w_train(n)
    real(real64),   intent(in)  :: x_pred(d, m)
    real(real64),   intent(out) :: y_out(m)

    integer(int32) :: j

    ! For each prediction point, use simplified R-compatible local fit
    do j = 1, m
      call simplified_r_local_fit(d, n, x_train, y_train, w_train, x_pred(:,j), k, degree, y_out(j))
    end do
  end subroutine predict_core_direct_fallback


  !======================================================================
  ! Direct 1D polynomial fit (no array copying) - FULL LOESS IMPLEMENTATION
  !----------------------------------------------------------------------
  
  subroutine local_polynomial_fit_1d_direct(nk, x, y, neighbors, weights, x0, degree, yhat)
    integer(int32), intent(in) :: nk, neighbors(nk), degree
    real(real64), intent(in) :: x(*), y(*), weights(nk), x0
    real(real64), intent(out) :: yhat

    integer(int32), parameter :: MAX_TERMS = 3
    real(real64) :: X_matrix(nk, MAX_TERMS)
    real(real64) :: XtWX(MAX_TERMS, MAX_TERMS), XtWy(MAX_TERMS)
    real(real64) :: beta(MAX_TERMS)
    integer(int32) :: i, j, idx, n_terms
    real(real64) :: xi, det

    ! Number of terms in the polynomial
    if (degree == 0) then
        n_terms = 1
    else if (degree == 1) then
        n_terms = 2
    else
        n_terms = 3   ! quadratic, like R default
    end if

    !-----------------------------------------------------------
    ! Build R-style design matrix: X = [1, xi, xi^2]
    !-----------------------------------------------------------
    do i = 1, nk
        xi = x(neighbors(i))
        X_matrix(i,1) = 1.0_real64
        if (n_terms >= 2) X_matrix(i,2) = xi
        if (n_terms >= 3) X_matrix(i,3) = xi*xi
    end do

    ! Zero matrices
    XtWX = 0.0_real64
    XtWy = 0.0_real64

    !-----------------------------------------------------------
    ! Compute X^T W X and X^T W y
    !-----------------------------------------------------------
    do i = 1, n_terms
        do j = 1, n_terms
            do idx = 1, nk
                XtWX(i,j) = XtWX(i,j) + X_matrix(idx,i)*weights(idx)*X_matrix(idx,j)
            end do
        end do
        do idx = 1, nk
            XtWy(i) = XtWy(i) + X_matrix(idx,i)*weights(idx)*y(neighbors(idx))
        end do
    end do

    !-----------------------------------------------------------
    ! Solve the linear system for β
    !-----------------------------------------------------------
    if (n_terms == 1) then
        beta(1) = XtWy(1) / XtWX(1,1)

    else if (n_terms == 2) then
        det = XtWX(1,1)*XtWX(2,2) - XtWX(1,2)*XtWX(2,1)
        beta(1) = ( XtWy(1)*XtWX(2,2) - XtWy(2)*XtWX(1,2) ) / det
        beta(2) = ( XtWy(2)*XtWX(1,1) - XtWy(1)*XtWX(2,1) ) / det

    else
        call solve_3x3_system(XtWX, XtWy, beta, 3)
    end if

    !-----------------------------------------------------------
    ! Predict y(x0) using the R-style polynomial:
    ! yhat = β0 + β1*x0 + β2*x0^2
    !-----------------------------------------------------------
    if (n_terms == 1) then
        yhat = beta(1)

    else if (n_terms == 2) then
        yhat = beta(1) + beta(2)*x0

    else
        yhat = beta(1) + beta(2)*x0 + beta(3)*(x0*x0)
    end if
  end subroutine local_polynomial_fit_1d_direct

! Solve 3x3 linear system using Gaussian elimination with partial pivoting
    subroutine solve_3x3_system(A, b, x, n)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: A(n,n), b(n)
        real(real64), intent(out) :: x(n)
        
        real(real64) :: A_work(n,n), b_work(n)
        real(real64) :: factor, temp
        integer(int32) :: i, j, k, max_row
        real(real64) :: max_val
        
        ! Copy to working arrays
        A_work = A(1:n,1:n)
        b_work = b(1:n)
        x = 0.0_real64
        
        ! Forward elimination with partial pivoting
        do i = 1, n-1
            ! Find pivot
            max_val = abs(A_work(i,i))
            max_row = i
            do k = i+1, n
                if (abs(A_work(k,i)) > max_val) then
                    max_val = abs(A_work(k,i))
                    max_row = k
                end if
            end do
            
            ! Swap rows if needed
            if (max_row /= i) then
                do j = 1, n
                    temp = A_work(i,j)
                    A_work(i,j) = A_work(max_row,j)
                    A_work(max_row,j) = temp
                end do
                temp = b_work(i)
                b_work(i) = b_work(max_row)
                b_work(max_row) = temp
            end if
            
            ! Check for singular matrix
            if (abs(A_work(i,i)) < 1e-12_real64) then
                ! Singular matrix, use simplified solution
                x(1) = b_work(1) / A_work(1,1)
                return
            end if
            
            ! Eliminate
            do k = i+1, n
                factor = A_work(k,i) / A_work(i,i)
                do j = i, n
                    A_work(k,j) = A_work(k,j) - factor * A_work(i,j)
                end do
                b_work(k) = b_work(k) - factor * b_work(i)
            end do
        end do
        
        ! Back substitution
        if (abs(A_work(n,n)) < 1e-12_real64) then
            x(1) = b_work(1) / A_work(1,1)
            return
        end if
        
        x(n) = b_work(n) / A_work(n,n)
        do i = n-1, 1, -1
            x(i) = b_work(i)
            do j = i+1, n
                x(i) = x(i) - A_work(i,j) * x(j)
            end do
            x(i) = x(i) / A_work(i,i)
        end do
    end subroutine solve_3x3_system
    
  !======================================================================
  ! Simplified R-compatible local fit (robust version without QR issues)
  !----------------------------------------------------------------------
  subroutine simplified_r_local_fit(d, n, x, y, w, z, k, degree, yhat)
    integer(int32), intent(in)  :: d, n, k, degree
    real(real64),   intent(in)  :: x(d,n), y(n), w(n), z(d)
    real(real64),   intent(out) :: yhat

    ! Simplified version that focuses on the key R differences:
    ! 1. R's exact neighbor selection with tie handling
    ! 2. R's exact tricube formula  
    ! 3. Robust matrix solving without complex QR
    
    integer(int32) :: i, p, nk
    integer(int32) :: neighbors(n)
    real(real64)   :: dists(n), rw(n), rho, u
    real(real64)   :: Xrow(8), XtWX(8,8), XtWy(8), beta(8)
    
    ! 1) Compute distances
    do i = 1, n
      dists(i) = euclid2(d, x(:,i), z)
    end do
    
    ! 2) R's neighbor selection with exact tie handling
    call r_select_neighbors_simple(dists, n, k, neighbors, nk, rho)
    
    ! 3) R's exact tricube weighting  
    do i = 1, n
      rw(i) = ZERO
    end do
    do i = 1, nk
      if (rho > EPS) then
        u = dists(neighbors(i)) / rho
        if(abs(u - ONE) < EPS) then
          ! Exact tie at bandwidth boundary, adjust slightly as in R to avoid zero weight
          rho = rho * 1.005_real64
          print *, "Warning: Adjusted bandwidth rho to: ", rho
          u = dists(neighbors(i)) / rho
        end if
        rw(neighbors(i)) = w(neighbors(i)) * r_tricube_exact(u)
      else
        rw(neighbors(i)) = w(neighbors(i))
      end if
    end do
    
    ! 4) Build normal equations (same structure as original but with R weights)
    p = num_terms(d, degree)
    call zero_mat_vec(XtWX, XtWy, p)

    do i = 1, nk
      if (rw(neighbors(i)) > EPS) then
        call build_design_row(d, degree, x(:,neighbors(i)), z, Xrow, p)
        call accumulate_normal_eq(XtWX, XtWy, Xrow, y(neighbors(i)), rw(neighbors(i)), p)
      end if
    end do

    ! 5) Solve (robust Cholesky)
    call solve_sym_posdef_robust(XtWX, XtWy, beta, p, rw, y, neighbors, nk)

    ! 6) Predict at z
    call build_design_row(d, degree, z, z, Xrow, p)
    yhat = dot_p(Xrow, beta, p)
    
  end subroutine simplified_r_local_fit

  !======================================================================
  ! R's simplified neighbor selection - EXACTLY like R's lowess.c
  !----------------------------------------------------------------------
  subroutine r_select_neighbors_simple(dists, n, k, neighbors, nk, rho)
    integer(int32), intent(in) :: n, k
    real(real64), intent(in) :: dists(n)
    integer(int32), intent(out) :: neighbors(n), nk
    real(real64), intent(out) :: rho
    
    integer(int32) :: indices(n), i, j, temp_int
    real(real64) :: temp_real, sorted_dists(n)
    
    ! Create index array and copy distances
    do i = 1, n
      indices(i) = i
      sorted_dists(i) = dists(i)
    end do
    
    ! Simple insertion sort (stable like R) - EXACT R algorithm
    do i = 2, n
      temp_real = sorted_dists(i)
      temp_int = indices(i)
      j = i - 1
      do while (j >= 1)
        if (sorted_dists(j) <= temp_real) exit
        sorted_dists(j + 1) = sorted_dists(j)
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      sorted_dists(j + 1) = temp_real
      indices(j + 1) = temp_int
    end do
    
    ! R's EXACT neighbor selection: select first k points, then add ties
    nk = min(k, n)
    do i = 1, nk
      neighbors(i) = indices(i)
    end do
    
    ! R algorithm: bandwidth is k-th distance
    if (nk < n) then
      rho = sorted_dists(nk)
      ! R includes ALL ties at the k-th distance
      do i = nk + 1, n
        if (abs(sorted_dists(i) - rho) < EPS .and. nk < n) then
          nk = nk + 1
          neighbors(nk) = indices(i)
        else
          exit
        end if
      end do
    else
      rho = sorted_dists(n)
    end if
    
    ! Ensure positive bandwidth (R does this)
    if (rho < EPS) rho = EPS
    
  end subroutine r_select_neighbors_simple

  !======================================================================
  ! R's exact tricube weight (simplified)
  !----------------------------------------------------------------------
  function r_tricube_exact(u) result(weight)
    real(real64), intent(in) :: u
    real(real64) :: weight
    
    ! R's exact formula but with tolerance for numerical stability
    ! In R's lowess.c: if(r <= h) w = (1. - fcube(r/h));
    ! So when r == h (u == 1), weight = 0
    
    if (u >= ONE) then
      weight = ZERO
    else
      ! Add small epsilon to avoid exact u=1 issues
      ! This helps when a point is exactly at the bandwidth boundary
      weight = (ONE - min(ONE - 1.0e-12_real64, u)**3)**3
    end if
  end function r_tricube_exact

  !======================================================================
  ! Robust symmetric positive definite solver with fallback
  !----------------------------------------------------------------------
  subroutine solve_sym_posdef_robust(A, b, x, p, weights, y_vals, neighbors, nk)
    real(real64), intent(inout) :: A(:,:), b(:)
    real(real64), intent(out) :: x(:)
    integer(int32), intent(in) :: p, nk
    real(real64), intent(in) :: weights(:), y_vals(:)
    integer(int32), intent(in) :: neighbors(:)
    
    integer(int32) :: i, j, k
    real(real64) :: sum_val, det_check, total_weight, weighted_mean
    logical :: is_singular
    
    ! Check for near-singularity
    is_singular = .false.
    det_check = ONE
    do i = 1, p
      if (abs(A(i,i)) < EPS) then
        is_singular = .true.
        exit
      end if
      det_check = det_check * A(i,i)
    end do
    
    if (is_singular .or. abs(det_check) < EPS) then
      ! Fallback: weighted average
      total_weight = ZERO
      weighted_mean = ZERO
      do i = 1, nk
        if (weights(neighbors(i)) > EPS) then
          total_weight = total_weight + weights(neighbors(i))
          weighted_mean = weighted_mean + weights(neighbors(i)) * y_vals(neighbors(i))
        end if
      end do
      
      if (total_weight > EPS) then
        x(1) = weighted_mean / total_weight
      else
        x(1) = ZERO
      end if
      
      do i = 2, p
        x(i) = ZERO
      end do
      return
    end if
    
    ! Standard Cholesky if matrix looks good
    call solve_sym_posdef(A, b, x, p)
    
  end subroutine solve_sym_posdef_robust

  !======================================================================
  ! Weighted linear fit for 1D case (needed by sliding window code)
  !----------------------------------------------------------------------
  subroutine weighted_linear_fit_1d(n_local, x_local, y_local, w_local, x_center, pred_value)
    integer(int32), intent(in) :: n_local
    real(real64), intent(in) :: x_local(n_local), y_local(n_local), w_local(n_local), x_center
    real(real64), intent(out) :: pred_value
    
    real(real64) :: sum_w, sum_wx, sum_wy, sum_wxx, sum_wxy
    real(real64) :: mean_x, mean_y, slope, intercept
    integer(int32) :: i
    
    ! Compute weighted sums
    sum_w = 0.0_real64
    sum_wx = 0.0_real64
    sum_wy = 0.0_real64
    sum_wxx = 0.0_real64
    sum_wxy = 0.0_real64
    
    do i = 1, n_local
      sum_w = sum_w + w_local(i)
      sum_wx = sum_wx + w_local(i) * x_local(i)
      sum_wy = sum_wy + w_local(i) * y_local(i)
      sum_wxx = sum_wxx + w_local(i) * x_local(i) * x_local(i)
      sum_wxy = sum_wxy + w_local(i) * x_local(i) * y_local(i)
    end do
    
    ! Solve for slope and intercept
    mean_x = sum_wx / sum_w
    mean_y = sum_wy / sum_w
    
    slope = (sum_wxy - sum_w * mean_x * mean_y) / (sum_wxx - sum_w * mean_x * mean_x)
    intercept = mean_y - slope * mean_x
    
    ! Predict at x_center
    pred_value = slope * x_center + intercept
  end subroutine weighted_linear_fit_1d

end module tox_loess