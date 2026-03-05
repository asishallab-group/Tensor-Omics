module anwil
  use safeguard
  use kd_tree, only: build_kd_index, kd_knn_query
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
  use f42_utils, only: sort_real, calculate_mean, count_valid, get_median
  use knn_smoothing_nadaraya_watson, only: smooth_vectors_gaussian_adaptive_nw
  implicit none

    interface
    pure subroutine dpotrf(uplo, n, a, lda, info)
      use iso_fortran_env, only: real64, int32
      character(len=1), intent(in) :: uplo
      integer(int32), intent(in)   :: n, lda
      real(real64), intent(inout)  :: a(lda, *)
      integer(int32), intent(out)  :: info
    end subroutine dpotrf

    pure subroutine dpotri(uplo, n, a, lda, info)
      use iso_fortran_env, only: real64, int32
      character(len=1), intent(in) :: uplo
      integer(int32), intent(in)   :: n, lda
      real(real64), intent(inout)  :: a(lda, *)
      integer(int32), intent(out)  :: info
    end subroutine dpotri
  end interface

contains

pure real(real64) function kernel_weight(d2, kernel_type) result(w)
  use iso_fortran_env, only: real64
  implicit none

  real(real64), intent(in) :: d2
  integer,      intent(in) :: kernel_type

  real(real64) :: r

  ! Numerical safety
  if (d2 < 0.0_real64) then
    w = 0.0_real64
    return
  end if

  select case (kernel_type)

  case (1)
    ! Standard Gaussian: exp(-d2)
    w = exp(-d2)

  case (2)
    ! Tricubic: d2 = r^2 / 2  -->  r = sqrt(2*d2)
    r = sqrt(2.0_real64 * d2)
    if (r < 1.0_real64) then
      w = (1.0_real64 - r**3)**3
    else
      w = 0.0_real64
    end if

  case default
    w = 0.0_real64

  end select
end function kernel_weight


subroutine anwil_smooth_sigma( &
    coords, vectors, smoothed, &
    n_coord_dims, n_vector_dims, n_points, k_neighbors, &
    kd_indices, dimension_order, workspace, value_buffer, &
    permutation, left_stack, right_stack, sigma_raw, sd_arr, &
    anisotropy_mode, anisotropy_factor, k_neighbors_sigma, kernel_type, ierr , &
    ! --- ARGUMENTS FOR MODE 3 ---
    U, singular_values, k_intrinsic )

    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err_once, is_ok, ERR_INVALID_INPUT
    use f42_utils, only: sort_real
    use kd_tree, only: build_kd_index, kd_knn_query
    implicit none

    ! Input/Output Arguments
    integer(int32), intent(in)    :: n_coord_dims, n_vector_dims, n_points, k_neighbors
    real(real64),    intent(in)    :: coords(n_coord_dims, n_points), vectors(n_vector_dims, n_points)
    real(real64),    intent(out)   :: smoothed(n_vector_dims, n_points), sigma_raw(1, n_points)
    integer(int32), intent(out)   :: ierr, kd_indices(n_points), dimension_order(n_coord_dims)
    integer(int32), intent(inout) :: workspace(n_points), permutation(n_points)
    integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
    real(real64),    intent(inout) :: value_buffer(n_points), sd_arr(n_vector_dims, n_points)
    integer(int32), intent(in) :: kernel_type ! 1 gaussian, 2 tricubic
    integer(int32), intent(in) :: k_neighbors_sigma
    
    ! Anisotropy Arguments
    integer(int32), intent(in)    :: anisotropy_mode
    real(real64),    intent(in)    :: anisotropy_factor
    integer(int32) :: i, recursion_stack(3, n_points), j 
    integer(int32) :: n_tmp(k_neighbors_sigma) 
    real(real64)   :: d_tmp(k_neighbors_sigma)

    ! Optional arguments for AMANLE mode
    real(real64), intent(in), optional :: U(:,:)
    real(real64), intent(in), optional :: singular_values(:)
    integer(int32), intent(in), optional :: k_intrinsic
    logical :: use_global_median
    real(real64)   :: global_median_val

    use_global_median = .false.   ! Set to .false. to use Nadaraya-Watson local smoothing

  call set_ok(ierr)

  if (anisotropy_mode == 3) then
    if (.not. present(U) .or. .not. present(singular_values) .or. .not. present(k_intrinsic)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if
  end if

  ! 1. Build KD-Tree
  dimension_order = [(i, i=1, n_coord_dims)]
  call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                      workspace, value_buffer, permutation, left_stack, right_stack, &
                      recursion_stack, ierr )
  if (.not. is_ok(ierr)) return

  ! Initialization
  smoothed(:,:) = 0.0_real64

  ! ==========================================================
  ! PHASE 1: RAW SIGMAS (Concurrent)
  !   - Sort by distance (using permutation)
  !   - Detect valid neighbors (exclude NaN/HUGE)
  !   - Calculate REAL median of distances excluding self (j=2..j_last)
  !   - c = 1/sqrt(log(k_eff+1)) 
  ! ==========================================================
  do concurrent (i = 1:n_points)
    block
      integer(int32) :: n_loc(k_neighbors), p_loc(k_neighbors), ierr_loc
      integer(int32) :: j, j_last, k_eff, j1, j2
      real(real64)   :: d_loc(k_neighbors)
      real(real64)   :: c_loc, local_median
      integer(int32) :: l_stack_loc(k_neighbors), r_stack_loc(k_neighbors)

      call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:,i), k_neighbors, n_loc, d_loc, ierr_loc)

      if (.not. is_ok(ierr_loc)) then
        sigma_raw(1,i) = 1.0e-12_real64
        cycle
      end if

      p_loc = [(j, j = 1, k_neighbors)]
      call sort_real(d_loc, p_loc, l_stack_loc, r_stack_loc)

      ! Find last valid index (NaN and HUGE are placed at the end by sort)
      j_last = 0
      do j = 1, k_neighbors
        if (d_loc(p_loc(j)) < 1.0e20_real64) j_last = j
      end do

      ! k_eff = number of real neighbors excluding self (assuming self is at j=1)
      if (j_last >= 2) then
        k_eff = j_last - 1
      else
        k_eff = 0
      end if

      if (k_eff < 1) then
        sigma_raw(1,i) = 1.0e-12_real64
      else
        ! Median of distances in range j=2..j_last (k_eff elements)
        if (mod(k_eff, 2) == 0) then
          ! Even: average of the two central elements
          j1 = 2 + (k_eff/2) - 1      
          j2 = j1 + 1
          local_median = 0.5_real64 * ( d_loc(p_loc(j1)) + d_loc(p_loc(j2)) )
        else
          ! Odd: the central element
          j1 = 2 + (k_eff - 1)/2
          local_median = d_loc(p_loc(j1))
        end if

        ! Kernel constant (including self in count: k_eff + 1)
        c_loc = 1.0_real64 / sqrt(log(real(k_eff + 1, real64)))

        sigma_raw(1,i) = max(c_loc * local_median, 1.0e-12_real64)
      end if

    end block
  end do


  ! ==========================================================
  ! PHASE 2: SIGMA SMOOTHING (Method Selection)
  ! ==========================================================
  
  if (use_global_median) then
      if (n_points > 0) then
          ! 1. Initialize permutation
          do i = 1, n_points
              permutation(i) = i
          end do

          ! 2. Sort based on raw sigmas
          call sort_real(sigma_raw(1, 1:n_points), &
                         permutation(1:n_points), &
                         left_stack, right_stack)

          ! 3. Extract global median and assign to all points
          if (mod(n_points, 2) /= 0) then
              ! Odd case
              i = permutation(n_points/2 + 1)
              global_median_val = sigma_raw(1, i)
          else
              ! Even case
              i = permutation(n_points/2)
              j = permutation(n_points/2 + 1)
              global_median_val = (sigma_raw(1, i) + sigma_raw(1, j)) * 0.5_real64
          end if
          
          sd_arr(1, 1:n_points) = global_median_val
      end if
      
  else
      ! ------------------------------------------------------
      ! LOCAL SMOOTHING: Nadaraya-Watson (NW)
      ! ------------------------------------------------------
        call smooth_vectors_gaussian_adaptive_nw( &
          coords,           & ! Pass X and Y (n_coord_dims = 2)
          sigma_raw,        & ! The noisy sigmas
          sd_arr(1:1, :),   & ! Result goes to first row of sd_arr
          n_coord_dims,     & 
          1,                & ! Smoothing 1 scalar (the sigma)
          n_points,         &
          k_neighbors_sigma,&
          kd_indices, dimension_order, n_tmp, d_tmp, & 
          workspace, value_buffer, permutation, left_stack, right_stack, &
          0.7_real64, ierr, 1)
  end if


  ! ==========================================================
  ! PHASE 3: SMOOTHING WITH ADAPTIVE ANISOTROPY
  ! ==========================================================
  do concurrent (i = 1:n_points)
    block
      integer(int32) :: n_loc(k_neighbors), p_loc(k_neighbors), ierr_loc, j, idx, r
      real(real64)   :: d_loc(k_neighbors), w_loc, wsum_loc, d2
      real(real64)   :: work_loc(n_vector_dims), delta(n_coord_dims)
      real(real64)   :: s_i, s_dot_ur, zeta_r
      real(real64)   :: displacement(n_vector_dims)
      integer(int32) :: l_stack_loc(k_neighbors), r_stack_loc(k_neighbors)

      ! 1. Neighbor Search
      call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:,i), k_neighbors, n_loc, d_loc, ierr_loc)

      ! 2. Local sorting to ensure consistency
      p_loc = [(j, j = 1, k_neighbors)]
      call sort_real(d_loc, p_loc, l_stack_loc, r_stack_loc) 

      ! 3. Determine Local Sigma (using the smoothed sd_arr from Phase 2)
      ! Factor 3.0 used as a general multiplier
      s_i = max(sd_arr(1,i), 1.0e-12_real64) * 3.0_real64 

      ! 4. Weighting Loop
      wsum_loc = 0.0_real64
      work_loc = 0.0_real64

      do j = 1, k_neighbors
        idx = n_loc(p_loc(j))
        if (idx <= 0 .or. idx == i) cycle
        
        delta = coords(:,idx) - coords(:,i)
        displacement = vectors(:, idx) - coords(:, i) 

        select case (anisotropy_mode)
        case (0) ! Isotropic
            d2 = (d_loc(p_loc(j))**2) / (2.0_real64 * s_i**2)
            
        case (3) ! AMANLE MODE (Global SVD + Adaptive Sigma)
              d2 = 0.0_real64
              do r = 1, n_coord_dims
                  ! Project displacement onto singular component r
                  s_dot_ur = dot_product(displacement, U(:, r))
                  
                  if (r <= k_intrinsic) then
                      ! TANGENTIAL: Use s_i (smoothed adaptive sigma)
                      d2 = d2 + (s_dot_ur**2) / (2.0_real64 * s_i**2)
                  else
                      ! NORMAL: Penalize with zeta_r
                      zeta_r = max(singular_values(r), 1.0e-12_real64)
                      ! Normal sigma is s_i scaled by anisotropy/singular values
                      d2 = d2 + (s_dot_ur**2) / (2.0_real64 * s_i**2 * zeta_r)
                  end if
              end do
          end select

        w_loc = kernel_weight(d2, kernel_type)

        work_loc = work_loc + w_loc * vectors(:, idx)
        wsum_loc = wsum_loc + w_loc
      end do

      ! 5. Final Assignment
      if (wsum_loc > 1.0e-18_real64) then
        smoothed(:, i) = work_loc / wsum_loc
      else
        smoothed(:, i) = vectors(:, i)
      end if
    end block
  end do

end subroutine anwil_smooth_sigma


  pure subroutine identity_matrix(A, n)
    integer(int32), intent(in) :: n
    real(real64),   intent(out) :: A(n,n)
    integer(int32) :: ii
    A(:,:) = 0.0_real64
    do ii = 1, n
      A(ii,ii) = 1.0_real64
    end do
  end subroutine identity_matrix

  pure subroutine invert_small_matrix_spd(A, n, ierr_local)
    ! In-place inverse of SPD matrix using Cholesky (dpotrf+dpotri).
    integer(int32), intent(in) :: n
    real(real64),   intent(inout) :: A(n,n)
    integer(int32), intent(out) :: ierr_local
    integer(int32) :: info, r, c

    call set_ok(ierr_local)

    call dpotrf('L', n, A, n, info)
    if (info /= 0) then
      call set_err_once(ierr_local, ERR_INVALID_INPUT)
      return
    end if

    call dpotri('L', n, A, n, info)
    if (info /= 0) then
      call set_err_once(ierr_local, ERR_INVALID_INPUT)
      return
    end if

    ! Symmetrize: copy lower triangle to upper triangle
    do r = 1, n
      do c = r+1, n
        A(r,c) = A(c,r)
      end do
    end do
  end subroutine invert_small_matrix_spd

end module anwil