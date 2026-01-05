module manle_module
  use iso_fortran_env, only: int32, real64
  use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
  use kd_tree, only: build_kd_index, kd_knn_query

  implicit none

  interface

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

  !> Subroutine to center the data
  subroutine center_module(data, n_points, n_dims, centered_data, ierr)
    integer(int32), intent(in) :: n_points, n_dims
    real(real64), intent(in) :: data(n_dims, n_points)
    real(real64), intent(out) :: centered_data(n_dims, n_points)
    integer(int32), intent(out) :: ierr
    integer(int32) :: i, j
    real(real64) :: mean_value

    ! Initialize ierr to 0 (success) using tox_errors
    call set_ok(ierr)

    ! Correct centering logic: subtract the mean of each dimension
    do i = 1, n_dims
      mean_value = sum(data(i, :)) / real(n_points, real64)
      do j = 1, n_points
        centered_data(i, j) = data(i, j) - mean_value
      end do
    end do
  end subroutine center_module

  !> Subroutine to compute the principal subspace using SVD
  subroutine compute_svd_randomized_f42( &
    data, n_points, n_dims, k, p, &
    Omega, Y, Q, B, &
    Stmp, Utmp, tau, &
    S, U, &
    work, lwork, ierr)

    ! ---------------------------
    ! INPUT
    ! ---------------------------
    integer(int32), intent(in) :: n_points, n_dims, k, p
    real(real64), intent(in)  :: data(n_dims, n_points)

    ! ---------------------------
    ! WORK BUFFERS (YA ALOCADOS FUERA)
    ! ---------------------------
    real(real64), intent(inout) :: Omega(n_points, k+p)
    real(real64), intent(inout) :: Y(n_dims,   k+p)
    real(real64), intent(inout) :: Q(n_dims,   k+p)
    real(real64), intent(inout) :: B(k+p, n_points)

    real(real64), intent(inout) :: Stmp(k+p)
    real(real64), intent(inout) :: Utmp(k+p, k+p)

    real(real64), intent(inout) :: work(lwork)
    integer(int32), intent(in)  :: lwork
     real(real64), intent(inout) :: tau(k+p) 

    ! ---------------------------
    ! OUTPUT
    ! ---------------------------
    real(real64), intent(out) :: S(k)
    real(real64), intent(out) :: U(n_dims, k)
    integer(int32), intent(out) :: ierr

    ! ---------------------------
    ! LOCALS
    ! ---------------------------
    integer(int32) :: l, info
    real(real64) :: work_query(1)
    integer(int32) :: lwork_svd

    ! Initialize ierr to 0 (success) using tox_errors
    call set_ok(ierr)

    ! Validate inputs
    if (n_points <= 0 .or. n_dims <= 0 .or. k <= 0 .or. p < 0) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    ! Validate work buffer size
    if (lwork < n_dims * (k + p)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        return
    end if

    l = min(k + p, n_dims)


    ! 1) Random Omega
    call random_number(Omega(:,1:l))

    ! 2) Y = A * Omega
    Y(:,1:l) = matmul(data, Omega(:,1:l))

    ! 3) QR → Q
    call qr_orthonormalize_f42(Y, Q, n_dims, l, tau, work, lwork, info)
    if (info /= 0) then
        ierr = info
        return
    end if

    ! 4) B = Q^T * A
    B(1:l,:) = matmul(transpose(Q(:,1:l)), data)

    ! ==========================
    ! 5) Small SVD of B (robust workspace query)
    ! ==========================
    

    ! Workspace query
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

    ! Actual SVD
    call dgesvd('S','S', l, n_points, B, l, &
                Stmp, Utmp, l, &
                B, l, work, lwork_svd, info)

    if (info /= 0) then
    ierr = info
    return
    end if


    ! 6) Reconstruction of the large U
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

        ! 2) Build explicit Q with only qcols columns
        call dorgqr(m, qcols, qcols, Q, m, tau, work, lwork, info)
    end subroutine qr_orthonormalize_f42


!> Projects data onto the k-dimensional principal subspace U_k
subroutine project_to_subspace(data, Uk, k, n_points, n_dims, projected_data, ierr)
  use iso_fortran_env, only: real64, int32
  use tox_errors,     only: set_ok, set_err_once, ERR_INVALID_INPUT
  implicit none

  integer(int32), intent(in) :: n_points, n_dims, k
  real(real64), intent(in)  :: data(n_dims, n_points)
  real(real64), intent(in)  :: Uk(n_dims, k)      ! ONLY the first k vectors
  real(real64), intent(out) :: projected_data(n_dims, n_points)
  integer(int32), intent(out) :: ierr

  real(real64) :: coords(k, n_points)

  call set_ok(ierr)

  if (n_points <= 0 .or. n_dims <= 0 .or. k <= 0 .or. k > n_dims) then
    call set_err_once(ierr, ERR_INVALID_INPUT)
    return
  end if

  ! 1) Coordinates in the subspace: C = U_k^T * V
  coords = matmul(transpose(Uk), data)

  ! 2) Reconstruction in the original space: V_proj = U_k * C
  projected_data = matmul(Uk, coords)

end subroutine project_to_subspace

  !> Subroutine to apply adaptive smoothing using the ANWIL module
  subroutine apply_anwil_smoothing(coords, vectors, smoothed, n_coord_dims, n_vector_dims, n_points, k_neighbors, ierr)
    use anwil, only: smooth_vectors_gaussian_adaptive
    integer(int32), intent(in) :: n_coord_dims, n_vector_dims, n_points, k_neighbors
    real(real64), intent(inout) :: coords(n_coord_dims, n_points)
    real(real64), intent(inout) :: vectors(n_vector_dims, n_points)
    real(real64), intent(out) :: smoothed(n_vector_dims, n_points)
    integer(int32), intent(out) :: ierr

    ! Buffers for k-d tree and smoothing
    integer(int32) :: kd_indices(n_points), dimension_order(n_coord_dims)
    integer(int32) :: neighbors(k_neighbors), workspace(n_points), permutation(n_points)
    integer(int32) :: left_stack(n_points), right_stack(n_points)
    real(real64) :: distances(k_neighbors), value_buffer(n_points)

    ! Call the smoothing subroutine from the ANWIL module
    call smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
                                          n_coord_dims, n_vector_dims, n_points, k_neighbors, &
                                          kd_indices, dimension_order, neighbors, distances, &
                                          workspace, value_buffer, permutation, left_stack, right_stack, 0, 1.0_real64, ierr)
  end subroutine apply_anwil_smoothing

  !> Subroutine to apply anisotropic smoothing using the AManLe algorithm and KD-Tree.
  !> The KD-Tree is built on 'coords' (manifold geometry) but the smoothing
  !> operates on 'vectors' (data space).
  subroutine apply_amanle_smoothing( &
      coords, vectors, smoothed, &
      n_coord_dims, n_vector_dims, n_points, k_neighbors, &
      k_intrinsic, U, singular_values, sigma_T, sigma_N, ierr, &
      kd_indices, dimension_order, workspace, value_buffer, &
      neighbors, distances, permutation, left_stack, right_stack )

      use iso_fortran_env, only: int32, real64
      use tox_errors, only: set_ok, is_ok, set_err_once, ERR_INVALID_INPUT
      implicit none

      ! --- INPUT ---
      integer(int32), intent(in) :: n_coord_dims, n_vector_dims, n_points, k_neighbors
      real(real64), intent(in) :: coords(n_coord_dims, n_points) ! Coords for KD-Tree (Manifold)
      real(real64), intent(in) :: vectors(n_vector_dims, n_points) ! Vectors to be smoothed (Data_Centered)
      integer(int32), intent(in) :: k_intrinsic                ! Dimensión intrínseca k
      real(real64), intent(in) :: U(n_vector_dims, n_vector_dims) ! Matriz de vectores singulares
      real(real64), intent(in) :: singular_values(n_vector_dims) ! Valores singulares
      real(real64), intent(in) :: sigma_T, sigma_N

      ! --- WORK BUFFERS (KD-Tree and KNN) ---
      integer(int32), intent(inout) :: kd_indices(n_points)
      integer(int32), intent(inout) :: dimension_order(n_coord_dims)
      integer(int32), intent(inout) :: workspace(n_points)
      real(real64), intent(inout) :: value_buffer(n_points)
      
      integer(int32), intent(inout) :: neighbors(k_neighbors)
      real(real64), intent(inout) :: distances(k_neighbors)
      integer(int32), intent(inout) :: permutation(k_neighbors)
      integer(int32), intent(inout) :: left_stack(n_points), right_stack(n_points)
      integer(int32) :: recursion_stack(3, n_points)
      ! --- OUTPUT ---
      real(real64), intent(out) :: smoothed(n_vector_dims, n_points)
      integer(int32), intent(out) :: ierr

      ! --- LOCALS ---
      real(real64) :: displacement(n_vector_dims)
      real(real64) :: local_wsum
      real(real64) :: d2, s_dot_ur, zeta_r, local_w
      integer(int32) :: i, j, r, idx, l
      
      call set_ok(ierr)

      ! 1. BUILD KD-TREE (Rebuilt in each iteration because 'coords' changes)
      dimension_order = [(l, l=1,n_coord_dims)]

      ! The construction must be accessible. If it is in the 'anwil' module, it must have its interface.

      ! Correct the call to build_kd_index
      call build_kd_index(coords, n_coord_dims, n_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
      if (.not. is_ok(ierr)) return

      ! Initialization of outputs
      smoothed(:,:) = 0.0_real64

      ! 2. ITERATION OVER EACH CENTRAL POINT (i)
      do i = 1, n_points

          ! A. QUERY K-NEAREST NEIGHBORS (Using coords - the current geometry of the manifold)

          call kd_knn_query(coords, kd_indices, n_coord_dims, n_points, dimension_order, &
                        coords(:, i), k_neighbors, neighbors, distances, value_buffer, ierr)
          if (.not. is_ok(ierr)) return

          ! Verify that the indices are within bounds before calling kd_knn_query
          if (any(neighbors < 1) .or. any(neighbors > n_points)) then
              call set_err_once(ierr, ERR_INVALID_INPUT)
              print *, "Error: Neighbor indices out of range in kd_knn_query."
              return
          end if

          local_wsum = 0.0_real64
          
          ! Iterate over the K_NEIGHBORS found:
          do j = 2, k_neighbors
              idx = neighbors(j) ! Use the neighbor's index

              ! B. CALCULATE DISPLACEMENT (in the ambient space, NOT the projected one)
              displacement = vectors(:, idx) - vectors(:, i) 
              
              d2 = 0.0_real64

              ! C. CALCULATE ANISOTROPIC DISTANCE D^2
              do r = 1, n_vector_dims ! n_vector_dims = d (ambient dimension)
                  s_dot_ur = dot_product(displacement, U(:, r))
                  
                  if (r <= k_intrinsic) then
                      ! TANGENT TERM (r=1 to k): Standard Gaussian smoothing
                      d2 = d2 + (s_dot_ur**2) / (2.0_real64 * sigma_T**2)
                  else
                      ! NORMAL TERM (r=k+1 to d): Penalty by SVD
                      zeta_r = singular_values(r)
                      
                      ! Regularization: Avoid division by very small singular values
                      if (zeta_r < 1.0e-12_real64) zeta_r = 1.0e-12_real64 
                      
                      d2 = d2 + (s_dot_ur**2) / (2.0_real64 * sigma_N**2 * zeta_r)
                  end if
              end do
              
              ! D. COMPUTE ANISOTROPIC WEIGHT w_ij = exp(-D_ij^2)
              local_w = exp(-d2)
              
              ! E. ACCUMULATE WEIGHTED VECTORS (ANWIL estimator)
              smoothed(:, i) = smoothed(:, i) + local_w * vectors(:, idx)
              local_wsum = local_wsum + local_w
          end do
          
          ! F. NORMALIZE
          if (local_wsum > 0.0_real64) then
              smoothed(:, i) = smoothed(:, i) / local_wsum
          else
              ! Fallback: no valid neighbors, keep the original value
              smoothed(:, i) = vectors(:, i)
          end if
      end do

  end subroutine apply_amanle_smoothing

subroutine manle_pipeline( &
    data, n_points, n_dims, &
    k_neighbors, n_iters_max, tol, &
    Omega, Y, Q, B, Stmp, Utmp, tau, work, lwork, &
    manifold, svd_line, ierr )

  use iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT
  implicit none

  integer(int32), intent(in) :: n_points, n_dims
  integer(int32), intent(in) :: k_neighbors, n_iters_max
  real(real64), intent(in)   :: tol

  real(real64), intent(inout) :: data(n_dims, n_points)

  real(real64), intent(inout) :: Omega(:, :)
  real(real64), intent(inout) :: Y(:, :)
  real(real64), intent(inout) :: Q(:, :)
  real(real64), intent(inout) :: B(:, :)
  real(real64), intent(inout) :: Stmp(:)
  real(real64), intent(inout) :: Utmp(:, :)
  real(real64), intent(inout) :: tau(:)
  real(real64), intent(inout) :: work(lwork)
  integer(int32), intent(in)  :: lwork

  real(real64), intent(out) :: manifold(n_dims, n_points)
  real(real64), intent(out) :: svd_line(n_dims, n_points)  ! Added output for SVD line
  integer(int32), intent(out) :: ierr

  ! -------------------------
  ! LOCALS
  ! -------------------------
  real(real64) :: Uk(n_dims, 1)
  real(real64) :: coords(1, n_points)
  real(real64) :: tmp(n_dims, n_points)
  real(real64) :: diff(n_dims, n_points)
  real(real64) :: delta
  real(real64) :: data_centered(n_dims, n_points)
  integer(int32) :: t, i

  call set_ok(ierr)

  ! -----------------------------------------
  ! 1) CENTRADO (SIEMPRE NECESARIO PARA SVD)
  ! -----------------------------------------
  call center_module(data, n_points, n_dims, data_centered, ierr)
   if (ierr /= 0) return


  ! -----------------------------------------
  ! 2) SVD + primer vector singular (u1)
  ! -----------------------------------------
  call compute_svd_randomized_f42( &
       data_centered, n_points, n_dims, &
       1, 5, &
       Omega, Y, Q, B, &
       Stmp, Utmp, tau, &
       Stmp(1:1), Uk, &  
       work, lwork, ierr )

  if (ierr /= 0) return

  

  ! -----------------------------------------
  ! 3) PROYECCIÓN INICIAL SOBRE span{u1}
  ! -----------------------------------------
  call project_to_subspace(data_centered, Uk, 1, n_points, n_dims, manifold, ierr)
  if (ierr /= 0) return
  svd_line = manifold

  ! -----------------------------------------
  ! 4) ITERACIONES MANLE
  ! -----------------------------------------
  do t = 1, n_iters_max

     coords(:,:) = matmul(transpose(Uk), manifold)

     call apply_anwil_smoothing( &
          coords, data_centered, tmp, &
          1, n_dims, n_points, k_neighbors, ierr )
     if (ierr /= 0) return

     diff  = tmp - manifold
     delta = maxval( sqrt(sum(diff**2, dim=1)) )

     manifold = tmp

     if (delta < tol) exit
  end do



end subroutine manle_pipeline

subroutine amanle_pipeline( &
    data, n_points, n_dims, &
    k_intrinsic, k_neighbors, n_iters_max, tol, &
    work, lwork, &
    manifold, svd_line, ierr )

  use iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT, is_ok
  implicit none

  ! INPUT
  integer(int32), intent(in) :: n_points, n_dims
  integer(int32), intent(in) :: k_intrinsic
  integer(int32), intent(in) :: k_neighbors, n_iters_max
  real(real64), intent(in)   :: tol
  real(real64), intent(in)   :: data(n_dims, n_points)   ! <- fixed originals
  real(real64), intent(inout):: work(lwork)
  integer(int32), intent(in) :: lwork

  ! OUTPUT
  real(real64), intent(out) :: manifold(n_dims, n_points)
  real(real64), intent(out) :: svd_line(n_dims, n_points)
  integer(int32), intent(out) :: ierr

  ! LOCALS
  real(real64) :: S_full(n_dims)
  real(real64), allocatable :: U_full(:,:)
  real(real64), allocatable :: Vt_dummy(:,:)
  integer(int32) :: info, lwork_svd
  real(real64) :: work_query(1)

  real(real64), allocatable :: data_centered(:,:)
  real(real64), allocatable :: tmp_svd(:,:)     ! scratch ONLY for dgesvd overwrite
  real(real64), allocatable :: tmp_next(:,:)    ! next manifold

  real(real64), allocatable :: diff(:,:)
  real(real64) :: delta
  real(real64) :: current_sigma_T, current_sigma_N

  ! KD / KNN buffers (persistentes)
  integer(int32) :: kd_indices(n_points), dimension_order(n_dims)
  integer(int32) :: neighbors(k_neighbors), workspace_knn(n_points), permutation(n_points)
  integer(int32) :: left_stack(n_points), right_stack(n_points)
  real(real64)   :: distances(k_neighbors), value_buffer_knn(n_points)

  integer(int32) :: t, r
  real(real64)   :: dotval

  call set_ok(ierr)

  allocate(U_full(n_dims, n_dims), stat=info)
  allocate(data_centered(n_dims, n_points), stat=info)
  allocate(tmp_svd(n_dims, n_points), stat=info)
  allocate(tmp_next(n_dims, n_points), stat=info)
  allocate(diff(n_dims, n_points), stat=info)

  if (info /= 0) then
    call set_err_once(ierr, ERR_INVALID_INPUT)
    print *, "Error: allocate(U_full) failed."
    return
  end if

  ! No necesitamos V^T, pero dgesvd pide argumento; damos dummy 1x1
  allocate(Vt_dummy(1, 1), stat=info)
  if (info /= 0) then
    call set_err_once(ierr, ERR_INVALID_INPUT)
    print *, "Error: allocate(Vt_dummy) failed."
    if (allocated(U_full)) deallocate(U_full)
    return
  end if

  ! Validación mínima
  if (k_intrinsic < 1 .or. k_intrinsic >= n_dims) then
    call set_err_once(ierr, ERR_INVALID_INPUT)
    print *, "Error: Invalid k_intrinsic."
    return
  end if

  ! 1) CENTRADO (una sola vez)
  call center_module(data, n_points, n_dims, data_centered, ierr)
  if (.not. is_ok(ierr)) then
    print *, "Error in center_module. ierr=", ierr
    return
  end if

  ! 2) SVD COMPLETO (una sola vez) sobre data_centered (= M)
  tmp_svd(:,:) = data_centered(:,:)

  call dgesvd('A', 'N', n_dims, n_points, tmp_svd, n_dims, S_full, &
              U_full, n_dims, Vt_dummy, 1, work_query, -1, info)
  if (info /= 0) then
    ierr = info
    print *, "Error in SVD workspace query. info=", info
    return
  end if
  lwork_svd = int(work_query(1), kind=int32)

  print *, "DBG: lwork passed=", lwork
  print *, "DBG: work_query(1)=", work_query(1)
  print *, "DBG: lwork_svd=", lwork_svd


  if (lwork < lwork_svd) then
    call set_err_once(ierr, ERR_INVALID_INPUT)
    print *, "Error: Insufficient workspace size. Need lwork >=", lwork_svd
    return
  end if

  tmp_svd(:,:) = data_centered(:,:)
  
  call dgesvd('A', 'N', n_dims, n_points, tmp_svd, n_dims, S_full, &
              U_full, n_dims, Vt_dummy, 1, work, lwork_svd, info)
  if (info /= 0) then
    ierr = info
    print *, "Error in SVD computation. info=", info
    return
  end if

  ! 3) PROYECCIÓN INICIAL AL SUBESPACIO TANGENTE: svd_line
  !    svd_line(:,i) = sum_{r=1..k} (v_i · u_r) u_r
  svd_line(:,:) = 0.0_real64
  do r = 1, k_intrinsic
    do concurrent (t = 1:n_points)
      dotval = dot_product(data_centered(:,t), U_full(:,r))
      svd_line(:,t) = svd_line(:,t) + dotval * U_full(:,r)
    end do
  end do

  ! Inicialización del manifold: arranca en la proyección (más fiel a tu texto)
  manifold = svd_line

  ! Sigmas (fijos aquí)
  current_sigma_T = 1.0_real64
  current_sigma_N = 1.0_real64

  ! 4) ITERACIONES: solo evoluciona manifold; targets y SVD quedan fijos
  do t = 1, n_iters_max

     call apply_amanle_smoothing( &
          manifold, data_centered, tmp_next, &
          n_dims, n_dims, n_points, k_neighbors, &
          k_intrinsic, U_full, S_full, current_sigma_T, current_sigma_N, ierr, &
          kd_indices, dimension_order, workspace_knn, value_buffer_knn, &
          neighbors, distances, permutation, left_stack, right_stack )

     if (.not. is_ok(ierr)) then
        print *, "Error in apply_amanle_smoothing. ierr=", ierr
        return
     end if

     diff  = tmp_next - manifold
     delta = maxval( sqrt(sum(diff**2, dim=1)) )
     print *, "Delta: ", delta

     manifold = tmp_next

     if (delta < tol) then
        print *, "Convergence achieved."
        exit
     end if
  end do

  print *, "amanle_pipeline completed successfully."

  if (allocated(U_full))    deallocate(U_full)
  if (allocated(Vt_dummy))  deallocate(Vt_dummy)

end subroutine amanle_pipeline


subroutine anwil_iterative( &
    data, n_points, n_dims, &
    k_neighbors, n_iters_max, tol, &
    tmp, manifold, ierr )

  use iso_fortran_env, only: real64, int32
  use tox_errors, only: set_ok
  implicit none

  ! INPUT
  integer(int32), intent(in) :: n_points, n_dims
  integer(int32), intent(in) :: k_neighbors, n_iters_max
  real(real64),   intent(in) :: tol
  real(real64),   intent(in) :: data(n_dims, n_points)   ! (x1, x2, ..., xn)

  ! WORK
  real(real64), intent(inout) :: tmp(n_dims, n_points)

  ! OUTPUT
  real(real64), intent(out) :: manifold(n_dims, n_points)
  integer(int32), intent(out) :: ierr

  ! LOCALS
  real(real64) :: coords(n_dims, n_points)
  real(real64) :: smoothed(n_dims, n_points)
  real(real64) :: diff(n_dims, n_points)
  real(real64) :: delta
  integer :: t

  call set_ok(ierr)

  ! -------------------------
  ! 1) INICIALIZAR COORDENADAS
  ! -------------------------
  coords(:,:) = data(:,:)
  manifold(:,:) = data(:,:)

  ! -------------------------
  ! 2) ITERACIONES MANLE (SUAVIZADO GLOBAL)
  ! -------------------------
  do t = 1, n_iters_max

     call apply_anwil_smoothing( &
          coords, manifold, smoothed, &
          n_dims, n_dims, n_points, k_neighbors, ierr )
     if (ierr /= 0) return

     diff  = smoothed - manifold
     delta = maxval(sqrt(sum(diff**2, dim=1)))

     manifold = smoothed

     if (delta < tol) exit
  end do

end subroutine anwil_iterative

end module manle_module
