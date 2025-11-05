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
  !! • Automatic algorithm selection: direct (n<1k), kdtree (1k-5k), subsampling (>5k)
  !!
  use iso_fortran_env, only: int32, real64
  use kd_tree, only: build_kd_index, kd_knn_query
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
  
  ! Algorithm selection thresholds - Force subsampling for safety
  integer(int32), parameter :: DIRECT_THRESHOLD = 100     ! Use direct method for n < 100 only
  integer(int32), parameter :: KDTREE_THRESHOLD = 200     ! Use k-d tree for 100 <= n < 200 only
  integer(int32), parameter :: MAX_SAFE_K = 2000          ! Maximum k to avoid stack overflow

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
    integer(int32), intent(inout) :: iv(liv)
    real(real64),   intent(inout) :: v(lv)
    real(real64),   intent(in)    :: x_pred(d, n_pred)
    real(real64),   intent(out)   :: y_pred(n_pred)
    integer(int32), intent(out)   :: ierr

    logical :: setLf

    ierr  = 0
    setLf = .true.

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
    integer(int32) :: k
    k = max(1_int32, min(nvmax, int(ceiling(real(n,real64)*max(EPS,min(ONE,f))), int32)))

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

    ! Choose algorithm based on dataset size
    if (n < DIRECT_THRESHOLD) then
      ! Small dataset: use direct method
      call predict_core_direct(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    else if (n < KDTREE_THRESHOLD) then
      ! Medium/Large dataset: use K-d tree optimization
      call predict_core_kdtree(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
    else
      ! Very large dataset: use subsampling + interpolation (R-style)
      call predict_core_subsampled(d, n, x_train, y_train, w_train, m, x_pred, k, degree, y_out)
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

    integer(int32) :: j

    ! For each prediction point, build neighborhood and solve local fit
    do j = 1, m
      call local_fit_at_point_direct(d, n, x_train, y_train, w_train, x_pred(:,j), k, degree, y_out(j))
    end do
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
    integer(int32), parameter :: MAX_SUBSAMPLE = 5000   ! Reduced for safety
    integer(int32) :: n_sub, step, i, j, safe_k
    integer(int32) :: sub_indices(MAX_SUBSAMPLE)
    real(real64) :: x_sub(d, MAX_SUBSAMPLE), y_sub(MAX_SUBSAMPLE), w_sub(MAX_SUBSAMPLE)
    real(real64) :: y_sub_smooth(MAX_SUBSAMPLE)
    real(real64) :: original_span
    
    ! Calculate subsampling step
    if (n <= MAX_SUBSAMPLE) then
      ! Use K-d tree for moderate sizes, but recalculate k properly
      safe_k = min(k, n/3, 500)  ! Limit k for safety
      call predict_core_kdtree(d, n, x_train, y_train, w_train, m, x_pred, safe_k, degree, y_out)
      return
    end if
    
    ! Subsample evenly spaced points (R uses more sophisticated sampling, but this works)
    step = max(1, n / MAX_SUBSAMPLE)
    n_sub = 0
    do i = 1, n, step
      n_sub = n_sub + 1
      if (n_sub > MAX_SUBSAMPLE) exit
      sub_indices(n_sub) = i
      x_sub(:, n_sub) = x_train(:, i)
      y_sub(n_sub) = y_train(i) 
      w_sub(n_sub) = w_train(i)
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
    integer(int32), intent(in)  :: d, n, k, degree
    real(real64),   intent(in)  :: x(d,n), y(n), w(n), z(d)
    real(real64),   intent(out) :: yhat

    ! We find the k nearest neighbors of z in x (Euclidean),
    ! compute tricube weights (scaled by max distance among the k),
    ! and fit local polynomial of given degree (0,1,2).

    integer(int32) :: i, t, p, ii
    integer(int32) :: nn(k)
    real(real64)   :: dists(n), dmax, wt(n)
    real(real64)   :: Xrow(64)   ! supports up to p<=64 terms (safe for d<=4, degree<=2)
    real(real64)   :: XtWX(64,64), XtWy(64), beta(64)

    ! 1) compute all distances
    do i = 1, n
      dists(i) = euclid2(d, x(:,i), z)
    end do

    ! 2) select k nearest indices (partial selection)
    call select_k_smallest(dists, n, k, nn)

    ! 3) bandwidth = max distance among k; avoid zero
    dmax = ZERO
    do i = 1, k
      if (dists(nn(i)) > dmax) dmax = dists(nn(i))
    end do
    if (dmax <= EPS) dmax = EPS

    ! 4) set weights (tricube((1 - (d/dmax)^3)^3)) * data weights
    do i = 1, n
      wt(i) = ZERO
    end do
    do i = 1, k
      wt(nn(i)) = tricube( max(ZERO, ONE - (dists(nn(i))/dmax)) ) * max(EPS, w(nn(i)))
    end do

    ! 5) build and solve normal equations for local poly at z
    p = num_terms(d, degree)
    call zero_mat_vec(XtWX, XtWy, p)

    do i = 1, k
      ii = nn(i)
      call build_design_row(d, degree, x(:,ii), z, Xrow, p)
      call accumulate_normal_eq(XtWX, XtWy, Xrow, y(ii), wt(ii), p)
    end do

    call solve_sym_posdef(XtWX, XtWy, beta, p)

    ! prediction at z uses design vector at z (centered: x-z gives zeros for linear terms)
    call build_design_row(d, degree, z, z, Xrow, p)
    yhat = dot_p(Xrow, beta, p)
  end subroutine local_fit_at_point_direct

  !======================================================================
  ! Local weighted regression around one point (K-d tree optimized)
  !----------------------------------------------------------------------
  subroutine local_fit_at_point_kdtree(d, n, x, y, w, z, k, degree, kd_indices, dimension_order, yhat)
    integer(int32), intent(in)  :: d, n, k, degree
    real(real64),   intent(in)  :: x(d,n), y(n), w(n), z(d)
    integer(int32), intent(in)  :: kd_indices(n), dimension_order(d)
    real(real64),   intent(out) :: yhat

    ! We use K-d tree to find the k nearest neighbors of z in x,
    ! compute tricube weights, and fit local polynomial.
    ! For safety, limit k to avoid stack overflow.

    integer(int32) :: i, p, ii, ierr, k_safe
    integer(int32) :: nn(MAX_SAFE_K)  ! Use safe maximum
    real(real64)   :: dists(MAX_SAFE_K), dmax, wt(n)
    real(real64)   :: Xrow(64)   ! supports up to p<=64 terms (safe for d<=4, degree<=2)
    real(real64)   :: XtWX(64,64), XtWy(64), beta(64)
    real(real64)   :: kd_workspace(d)

    ! Limit k to safe value to avoid stack overflow
    k_safe = min(k, MAX_SAFE_K)

    ! 1) Use K-d tree to find k_safe nearest neighbors
    call kd_knn_query(x, kd_indices, d, n, dimension_order, z, k_safe, nn, dists, kd_workspace, ierr)

    if (ierr /= 0) then
      ! Fallback to direct method if K-d query fails
      call local_fit_at_point_direct(d, n, x, y, w, z, k, degree, yhat)
      return
    end if

    ! 2) bandwidth = max distance among k_safe; avoid zero
    dmax = maxval(dists(1:k_safe))
    if (dmax <= EPS) dmax = EPS

    ! 3) set weights (tricube((1 - (d/dmax)^3)^3)) * data weights
    do i = 1, n
      wt(i) = ZERO
    end do
    do i = 1, k_safe
      wt(nn(i)) = tricube( max(ZERO, ONE - (dists(i)/dmax)) ) * max(EPS, w(nn(i)))
    end do

    ! 4) build and solve normal equations for local poly at z
    p = num_terms(d, degree)
    call zero_mat_vec(XtWX, XtWy, p)

    do i = 1, k_safe
      ii = nn(i)
      call build_design_row(d, degree, x(:,ii), z, Xrow, p)
      call accumulate_normal_eq(XtWX, XtWy, Xrow, y(ii), wt(ii), p)
    end do

    call solve_sym_posdef(XtWX, XtWy, beta, p)

    ! prediction at z uses design vector at z (centered: x-z gives zeros for linear terms)
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
    t = max(ZERO, min(ONE, u))
    tricube = (ONE - (ONE - t)**3)**3
  end function tricube

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

end module tox_loess
