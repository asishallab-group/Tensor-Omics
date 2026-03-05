module manle_module
  use safeguard
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
  use kd_tree, only: build_kd_index, kd_knn_query
  use anwil, only: anwil_smooth_sigma
  use f42_utils, only: get_approx_diameter, compute_rmse, compute_roughness, compute_coverage, compute_smoothing_score

  implicit none

  interface

  ! External LAPACK routines for QR and SVD
  subroutine dgeqrf(m, n, a, lda, tau, work, lwork, info)
    use iso_fortran_env, only: real64, int32
    integer(int32), intent(in) :: m, n, lda, lwork
    real(real64), intent(inout) :: a(lda, *)
    real(real64), intent(out) :: tau(*)
    real(real64), intent(inout) :: work(*)
    integer(int32), intent(out) :: info
  end subroutine dgeqrf

  subroutine dorgqr(m, n, k, a, lda, tau, work, lwork, info)
    use iso_fortran_env, only: real64, int32
    integer(int32), intent(in) :: m, n, k, lda, lwork
    real(real64), intent(inout) :: a(lda, *)
    real(real64), intent(in) :: tau(*)
    real(real64), intent(inout) :: work(*)
    integer(int32), intent(out) :: info
  end subroutine dorgqr

  subroutine dgesvd(jobu, jobvt, m, n, a, lda, s, u, ldu, vt, ldvt, &
                    work, lwork, info)
    use iso_fortran_env, only: real64, int32
    character(len=1), intent(in) :: jobu, jobvt
    integer(int32), intent(in) :: m, n, lda, ldu, ldvt, lwork
    real(real64), intent(inout) :: a(lda, *)
    real(real64), intent(out) :: s(*)
    real(real64), intent(out) :: u(ldu, *)
    real(real64), intent(out) :: vt(ldvt, *)
    real(real64), intent(inout) :: work(*)
    integer(int32), intent(out) :: info
  end subroutine dgesvd

end interface

contains

  !> Subroutine to center the data by subtracting the mean of each dimension
  subroutine center_module(data, n_points, n_dims, centered_data, ierr)
    integer(int32), intent(in) :: n_points, n_dims
    real(real64), intent(in) :: data(n_dims, n_points)
    real(real64), intent(out) :: centered_data(n_dims, n_points)
    integer(int32), intent(out) :: ierr
    integer(int32) :: i, j
    real(real64) :: mean_value

    call set_ok(ierr)

    do i = 1, n_dims
      mean_value = sum(data(i, :)) / real(n_points, real64)
      do j = 1, n_points
        centered_data(i, j) = data(i, j) - mean_value
      end do
    end do
  end subroutine center_module

  !> Computes the principal subspace using a randomized SVD approach
  subroutine compute_svd_randomized_f42( &
    data, n_points, n_dims, k, p, &
    Omega, Y, Q, B, &
    Stmp, Utmp, tau, &
    S, U, &
    work, lwork, ierr)

    ! Input
    integer(int32), intent(in) :: n_points, n_dims, k, p
    real(real64), intent(in)  :: data(n_dims, n_points)

    ! Work Buffers (pre-allocated externally)
    real(real64), intent(inout) :: Omega(n_points, k+p)
    real(real64), intent(inout) :: Y(n_dims,   k+p)
    real(real64), intent(inout) :: Q(n_dims,   k+p)
    real(real64), intent(inout) :: B(k+p, n_points)
    real(real64), intent(inout) :: Stmp(k+p)
    real(real64), intent(inout) :: Utmp(k+p, k+p)
    real(real64), intent(inout) :: work(lwork)
    integer(int32), intent(in)  :: lwork
    real(real64), intent(inout) :: tau(k+p) 

    ! Output
    real(real64), intent(out) :: S(k)
    real(real64), intent(out) :: U(n_dims, k)
    integer(int32), intent(out) :: ierr

    ! Locals
    integer(int32) :: l, info
    real(real64) :: work_query(1)
    integer(int32) :: lwork_svd

    call set_ok(ierr)

    if (n_points <= 0 .or. n_dims <= 0 .or. k <= 0 .or. p < 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    if (lwork < n_dims * (k + p)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    l = min(k + p, n_dims)

    ! 1) Generate Random Matrix Omega
    call random_number(Omega(:,1:l))

    ! 2) Sample Matrix Y = A * Omega
    Y(:,1:l) = matmul(data, Omega(:,1:l))

    ! 3) Orthonormalize Y to find Basis Q
    call qr_orthonormalize_f42(Y, Q, n_dims, l, tau, work, lwork, info)
    if (info /= 0) then
        ierr = info
        return
    end if

    ! 4) Project data to smaller subspace: B = Q^T * A
    B(1:l,:) = matmul(transpose(Q(:,1:l)), data)

    ! 5) Perform SVD on the smaller matrix B
    ! Workspace query for dgesvd
    call dgesvd('S','S', l, n_points, B, l, &
                Stmp, Utmp, l, &
                B, l, work_query, -1, info)

    if (info /= 0) then
      ierr = info
      return
    end if

    lwork_svd = int(work_query(1))

    if (lwork < lwork_svd) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! Actual SVD execution
    call dgesvd('S','S', l, n_points, B, l, &
                Stmp, Utmp, l, &
                B, l, work, lwork_svd, info)

    if (info /= 0) then
      ierr = info
      return
    end if

    ! 6) Reconstruct the ambient space U
    U = matmul(Q(:,1:l), Utmp(1:l, 1:k))
    S = Stmp(1:k)

  end subroutine compute_svd_randomized_f42

  subroutine qr_orthonormalize_f42(A, Q, m, n, tau, work, lwork, info)
    integer(int32), intent(in) :: m, n, lwork
    real(real64), intent(in)  :: A(m,n)
    real(real64), intent(out) :: Q(m,n)
    real(real64), intent(inout) :: tau(min(m,n))
    real(real64), intent(inout) :: work(lwork)
    integer(int32), intent(out) :: info
    integer(int32) :: qcols

    Q = A
    qcols = min(m, n)

    ! 1) QR factorization
    call dgeqrf(m, qcols, Q, m, tau, work, lwork, info)
    if (info /= 0) return

    ! 2) Build explicit Q matrix
    call dorgqr(m, qcols, qcols, Q, m, tau, work, lwork, info)
  end subroutine qr_orthonormalize_f42

  !> Projects data onto the k-dimensional principal subspace Uk
  subroutine project_to_subspace(data, Uk, k, n_points, n_dims, projected_data, ierr)
    use iso_fortran_env, only: real64, int32
    use tox_errors,     only: set_ok, set_err_once, ERR_INVALID_INPUT
    implicit none

    integer(int32), intent(in) :: n_points, n_dims, k
    real(real64), intent(in)  :: data(n_dims, n_points)
    real(real64), intent(in)  :: Uk(n_dims, k)
    real(real64), intent(out) :: projected_data(n_dims, n_points)
    integer(int32), intent(out) :: ierr
    real(real64) :: coords(k, n_points)

    call set_ok(ierr)

    if (n_points <= 0 .or. n_dims <= 0 .or. k <= 0 .or. k > n_dims) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! 1) Subspace coordinates: C = Uk^T * V
    coords = matmul(transpose(Uk), data)

    ! 2) Ambient space reconstruction: V_proj = Uk * C
    projected_data = matmul(Uk, coords)
  end subroutine project_to_subspace

  !> Wrapper for ANWIL adaptive smoothing
  subroutine apply_anwil_smoothing(coords, vectors, smoothed, n_coord_dims, n_vector_dims, n_points, k_neighbors, k_neighbors_sigma, kernel_type, ierr)
    integer(int32), intent(in) :: n_coord_dims, n_vector_dims, n_points, k_neighbors
    real(real64), intent(inout) :: coords(n_coord_dims, n_points)
    real(real64), intent(inout) :: vectors(n_vector_dims, n_points)
    real(real64), intent(out) :: smoothed(n_vector_dims, n_points)
    integer(int32), intent(out) :: ierr
    integer(int32), intent(in) :: kernel_type
    integer(int32), intent(in) :: k_neighbors_sigma

    integer(int32) :: kd_indices(n_points), dimension_order(n_coord_dims)
    integer(int32) :: workspace(n_points), permutation(n_points)
    integer(int32) :: left_stack(n_points), right_stack(n_points)
    real(real64) :: value_buffer(n_points), sd_arr(n_coord_dims,n_points), sigma_raw(1,n_points)

    call anwil_smooth_sigma( coords, vectors, smoothed, &
                             n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                             kd_indices, dimension_order, &
                             workspace, value_buffer, permutation, left_stack, right_stack, sigma_raw, sd_arr, 0, 1.0_real64, k_neighbors_sigma, kernel_type, ierr )
  end subroutine apply_anwil_smoothing

  !> Applies anisotropic smoothing using the AManLe algorithm and KD-Tree.
  !> The tree is built on manifold geometry, but smoothing operates on data space.
  subroutine apply_amanle_smoothing( &
      coords, vectors, smoothed, &
      n_coord_dims, n_vector_dims, n_points, k_neighbors, &
      k_intrinsic, U, singular_values, sigma_T, sigma_N, ierr, &
      kd_indices, dimension_order, workspace, value_buffer, &
      neighbors, distances, permutation, left_stack, right_stack )

      use iso_fortran_env, only: int32, real64
      use tox_errors, only: set_ok, is_ok, set_err_once, ERR_INVALID_INPUT
      implicit none

      ! Input
      integer(int32), intent(in) :: n_coord_dims, n_vector_dims, n_points, k_neighbors
      real(real64), intent(in) :: coords(n_coord_dims, n_points) 
      real(real64), intent(in) :: vectors(n_vector_dims, n_points) 
      integer(int32), intent(in) :: k_intrinsic                 
      real(real64), intent(in) :: U(n_vector_dims, n_vector_dims) 
      real(real64), intent(in) :: singular_values(n_vector_dims) 
      real(real64), intent(in) :: sigma_T, sigma_N

      ! Workspace
      integer(int32), intent(inout) :: kd_indices(n_points)
      integer(int32), intent(inout) :: dimension_order(n_coord_dims)
      integer(int32), intent(inout) :: workspace(n_points)
      real(real64), intent(inout) :: value_buffer(n_points)
      integer(int32), intent(inout) :: neighbors(k_neighbors)
      real(real64), intent(inout) :: distances(k_neighbors)
      integer(int32), intent(inout) :: permutation(k_neighbors)
      integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
      integer(int32) :: recursion_stack(3, n_points)

      ! Output
      real(real64), intent(out) :: smoothed(n_vector_dims, n_points)
      integer(int32), intent(out) :: ierr

      ! Locals
      real(real64) :: displacement(n_vector_dims)
      real(real64) :: local_wsum
      real(real64) :: d2, s_dot_ur, zeta_r, local_w
      integer(int32) :: i, j, r, idx, l
      
      call set_ok(ierr)

      ! 1. Build KD-Tree on manifold coordinates
      dimension_order = [(l, l=1,n_coord_dims)]
      call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
      if (.not. is_ok(ierr)) return

      smoothed(:,:) = 0.0_real64

      ! 2. Iterate over each point for anisotropic smoothing
      do i = 1, n_points
          call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:, i), k_neighbors, neighbors, distances, ierr)
          if (.not. is_ok(ierr)) return

          local_wsum = 0.0_real64
          
          do j = 2, k_neighbors
              idx = neighbors(j)
              displacement = vectors(:, idx) - coords(:, i)
              d2 = 0.0_real64

              ! Compute Anisotropic Distance D^2
              do r = 1, n_vector_dims
                  s_dot_ur = dot_product(displacement, U(:, r))
                  
                  if (r <= k_intrinsic) then
                      ! TANGENT space: Standard smoothing
                      d2 = d2 + (s_dot_ur**2) / (2.0_real64 * sigma_T**2)
                  else
                      ! NORMAL space: Penalty weighted by singular values
                      zeta_r = singular_values(r)
                      if (zeta_r < 1.0e-12_real64) zeta_r = 1.0e-12_real64 
                      d2 = d2 + (s_dot_ur**2) / (2.0_real64 * sigma_N**2 * zeta_r)
                  end if
              end do
              
              local_w = exp(-d2)
              smoothed(:, i) = smoothed(:, i) + local_w * vectors(:, idx)
              local_wsum = local_wsum + local_w
          end do
          
          if (local_wsum > 0.0_real64) then
              smoothed(:, i) = smoothed(:, i) / local_wsum
          else
              smoothed(:, i) = vectors(:, i)
          end if
      end do
  end subroutine apply_amanle_smoothing

  !> Execution pipeline for standard MANLE algorithm
  subroutine manle_pipeline( &
      data, n_points, n_dims, &
      k_neighbors, n_iters_max, tol, &
      Omega, Y, Q, B, Stmp, Utmp, tau, work, lwork, &
      manifold, svd_line, k_neighbors_sigma, kernel_type, ierr )

    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
    implicit none

    integer(int32), intent(in) :: n_points, n_dims
    integer(int32), intent(in) :: k_neighbors, n_iters_max
    real(real64), intent(in)   :: tol
    real(real64), intent(inout) :: data(n_dims, n_points)
    real(real64), intent(inout) :: Omega(:, :), Y(:, :), Q(:, :), B(:, :)
    real(real64), intent(inout) :: Stmp(:), Utmp(:, :), tau(:)
    real(real64), intent(inout) :: work(lwork)
    integer(int32), intent(in)  :: lwork
    integer(int32), intent(in) :: kernel_type, k_neighbors_sigma
    real(real64), intent(out) :: manifold(n_dims, n_points)
    real(real64), intent(out) :: svd_line(n_dims, n_points)
    integer(int32), intent(out) :: ierr

    real(real64) :: Uk(n_dims, 1)
    real(real64) :: coords(1, n_points)
    real(real64) :: tmp(n_dims, n_points), diff(n_dims, n_points)
    real(real64) :: delta, data_centered(n_dims, n_points)
    integer(int32) :: t

    call set_ok(ierr)

    ! 1) Data Centering
    call center_module(data, n_points, n_dims, data_centered, ierr)
    if (ierr /= 0) return

    ! 2) Compute Initial SVD (randomized)
    call compute_svd_randomized_f42(data_centered, n_points, n_dims, 1, 5, &
         Omega, Y, Q, B, Stmp, Utmp, tau, Stmp(1:1), Uk, work, lwork, ierr )
    if (ierr /= 0) return

    ! 3) Initial projection onto first singular vector
    call project_to_subspace(data_centered, Uk, 1, n_points, n_dims, manifold, ierr)
    if (ierr /= 0) return
    svd_line = manifold

    ! 4) Iterate MANLE refinement
    do t = 1, n_iters_max
       coords(:,:) = matmul(transpose(Uk), manifold)
       call apply_anwil_smoothing(coords, data_centered, tmp, 1, n_dims, n_points, &
                                  k_neighbors, k_neighbors_sigma, kernel_type, ierr )
       if (ierr /= 0) return
       diff  = tmp - manifold
       delta = maxval( sqrt(sum(diff**2, dim=1)) )
       manifold = tmp
       if (delta < tol) exit
    end do
  end subroutine manle_pipeline

  !> Execution pipeline for anisotropic AManLe algorithm
  subroutine amanle_pipeline( &
      data, n_points, n_dims, &
      k_intrinsic, k_neighbors, n_iters_max, tol, &
      work, lwork, &
      manifold, svd_line, k_neighbors_sigma, kernel_type, ierr )

    use iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT, is_ok
    implicit none

    integer(int32), intent(in) :: n_points, n_dims, k_intrinsic, k_neighbors, n_iters_max
    real(real64), intent(in)   :: tol, data(n_dims, n_points)
    real(real64), intent(inout):: work(lwork)
    integer(int32), intent(in) :: lwork, kernel_type, k_neighbors_sigma
    real(real64), intent(out) :: manifold(n_dims, n_points), svd_line(n_dims, n_points)
    integer(int32), intent(out) :: ierr

    real(real64) :: S_full(n_dims)
    real(real64), allocatable :: U_full(:,:), Vt_dummy(:,:)
    integer(int32) :: info, lwork_svd
    real(real64) :: work_query(1)
    real(real64), allocatable :: data_centered(:,:), tmp_svd(:,:), tmp_next(:,:), diff(:,:)
    real(real64) :: delta, dotval

    integer(int32) :: kd_indices(n_points), dimension_order(n_dims)
    integer(int32) :: workspace_knn(n_points), permutation(n_points)
    integer(int32) :: left_stack(n_points), right_stack(n_points)
    real(real64)   :: value_buffer_knn(n_points), sd_arr(n_dims, n_points), sigma_raw(1, n_points)
    integer(int32) :: t, r

    call set_ok(ierr)
    allocate(U_full(n_dims, n_dims), data_centered(n_dims, n_points), tmp_svd(n_dims, n_points), &
             tmp_next(n_dims, n_points), diff(n_dims, n_points), Vt_dummy(1,1), stat=info)
    if (info /= 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    ! 1) Data Centering
    call center_module(data, n_points, n_dims, data_centered, ierr)
    if (.not. is_ok(ierr)) return

    ! 2) Full SVD for initial geometry
    tmp_svd(:,:) = data_centered(:,:)
    call dgesvd('A', 'N', n_dims, n_points, tmp_svd, n_dims, S_full, &
                U_full, n_dims, Vt_dummy, 1, work_query, -1, info)
    lwork_svd = int(work_query(1), kind=int32)

    if (lwork < lwork_svd) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if

    tmp_svd(:,:) = data_centered(:,:)
    call dgesvd('A', 'N', n_dims, n_points, tmp_svd, n_dims, S_full, &
                U_full, n_dims, Vt_dummy, 1, work, lwork_svd, info)

    ! 3) Initial Projection
    svd_line(:,:) = 0.0_real64
    do r = 1, k_intrinsic
      do concurrent (t = 1:n_points)
        dotval = dot_product(data_centered(:,t), U_full(:,r))
        svd_line(:,t) = svd_line(:,t) + dotval * U_full(:,r)
      end do
    end do
    manifold = svd_line

    ! 4) Refinement loop with AManLe anisotropic smoothing
    do t = 1, n_iters_max
      call anwil_smooth_sigma(manifold, data_centered, tmp_next, n_dims, n_dims, n_points, k_neighbors, &
          kd_indices, dimension_order, workspace_knn, value_buffer_knn, permutation, left_stack, right_stack, &
          sigma_raw, sd_arr, 3, 1.0_real64, k_neighbors_sigma, kernel_type, ierr, U_full, S_full, k_intrinsic)

      if (.not. is_ok(ierr)) return
      diff = tmp_next - manifold
      delta = maxval( sqrt(sum(diff**2, dim=1)) )
      manifold = tmp_next
      if (delta < tol) exit
    end do

    if (allocated(U_full)) deallocate(U_full)
    if (allocated(Vt_dummy)) deallocate(Vt_dummy)
  end subroutine amanle_pipeline

  !> Iterative ANWIL optimizer with automated stop conditions and scoring history
  subroutine anwil_iterative( &
      data, n_points, n_dims, k_neighbors, n_iters_max, &
      patience_k, tol_rel, min_iters, k_neighbors_sigma, kernel_type, &
      anisotropy_mode, anisotropy_factor, &
      manifold, ierr, method_flag, w_r, w_e, w_c, &
      history_scores, history_coverage, history_penalty, history_rmse, history_roughness, &
      stop_iter, stop_reason, best_iter_out, best_score_out )

      use iso_fortran_env, only: real64, int32
      use tox_errors,      only: set_ok, is_ok
      implicit none

      ! Inputs/Outputs
      integer, intent(in)        :: method_flag ! 1 = Arithmetic, 2 = Geometric
      real(real64), intent(in)   :: w_r, w_e, w_c ! Weights for Roughness, Error, Coverage
      integer(int32), intent(in) :: n_points, n_dims, k_neighbors, n_iters_max
      integer(int32), intent(in) :: patience_k, min_iters, k_neighbors_sigma, kernel_type
      integer(int32), intent(in) :: anisotropy_mode
      real(real64),   intent(in) :: tol_rel, anisotropy_factor
      real(real64),   intent(in) :: data(n_dims, n_points)
      real(real64),   intent(out):: manifold(n_dims, n_points)
      integer(int32), intent(out):: ierr

      ! History/Metrics
      real(real64),   intent(out):: history_scores(n_iters_max), history_roughness(n_iters_max)
      real(real64),   intent(out):: history_rmse(n_iters_max), history_coverage(n_iters_max), history_penalty(n_iters_max)
      integer(int32), intent(out):: stop_iter, stop_reason, best_iter_out
      real(real64),   intent(out):: best_score_out

      ! Workspace
      real(real64) :: coords_t(n_dims, n_points), smoothed_t(n_dims, n_points), best_state(n_dims, n_points)
      real(real64) :: sigma_raw(1, n_points), sd_arr(1, n_points), value_buffer(n_points)
      integer(int32) :: kd_indices(n_points), dimension_order(n_dims), workspace(n_points)
      integer(int32) :: permutation(n_points), left_stack(n_points), right_stack(n_points), recursion_stack(3, n_points)
      integer(int32) :: nb_fix(k_neighbors, n_points)
      real(real64)   :: dist_fix(k_neighbors, n_points)

      real(real64) :: r_ref, rt, rmse_t, d0, score_t, best_score, pen_t, cov_t
      integer(int32) :: i1, i2, t, no_improve, i, d, ierr_loc
      real(real64), parameter :: eps_score = 1.0d-12

      call set_ok(ierr)
      history_scores = huge(1.0_real64)
      stop_reason = 0; stop_iter = 0; best_iter_out = 0

      ! 1) Init State
      coords_t = data; manifold = data; best_state = data
      best_score = huge(1.0_real64); no_improve = 0

      ! 2) Build Geometric Reference (Diameter and Fixed Neighbors)
      call get_approx_diameter(data, n_points, n_dims, d0, i1, i2)
      dimension_order = [(d, d=1_int32, n_dims)]
      call build_kd_index(data, n_dims, n_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr )
      
      do i = 1, n_points
          call kd_knn_query(data, kd_indices, n_dims, n_points, dimension_order, &
                            data(:,i), k_neighbors, nb_fix(:,i), dist_fix(:,i), ierr_loc)
      end do

      ! 3) Optimization Refinement Loop
      do t = 1, n_iters_max
          stop_iter = t
          call anwil_smooth_sigma(coords_t, manifold, smoothed_t, n_dims, n_dims, n_points, k_neighbors, &
               kd_indices, dimension_order, workspace, value_buffer, permutation, left_stack, right_stack, &
               sigma_raw, sd_arr, anisotropy_mode, anisotropy_factor, k_neighbors_sigma, kernel_type, ierr )

          if (t == 1) then
              call compute_roughness(data, sd_arr(1,:), nb_fix, dist_fix, n_points, n_dims, k_neighbors, r_ref)
              r_ref = max(r_ref, eps_score)
          end if

          ! Compute Scoring Metrics
          call compute_roughness(smoothed_t, sd_arr(1,:), nb_fix, dist_fix, n_points, n_dims, k_neighbors, rt)
          call compute_rmse(data, smoothed_t, n_points, n_dims, rmse_t)
          call compute_coverage(d0, i1, i2, smoothed_t, n_points, n_dims, cov_t, pen_t)
          call compute_smoothing_score(method_flag, w_r, w_e, w_c, rt, r_ref, rmse_t, d0, pen_t, eps_score, score_t)

          history_scores(t) = score_t; history_roughness(t) = rt; history_rmse(t) = rmse_t; history_coverage(t) = cov_t; history_penalty(t) = pen_t

          ! Convergence and Patience Logic
          if (score_t < best_score * (1.0_real64 - tol_rel)) then
              best_score = score_t; best_iter_out = t; best_state = smoothed_t; no_improve = 0
          else
              no_improve = no_improve + 1
          end if

          manifold = smoothed_t
          if (t >= min_iters .and. no_improve >= patience_k) then
              stop_reason = 1; exit ! Patience exceeded
          end if
      end do

      if (stop_reason == 0) stop_reason = 2 ! Max iterations reached
      manifold = best_state; best_score_out = best_score
  end subroutine anwil_iterative

end module manle_module