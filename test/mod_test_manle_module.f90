!> @file mod_test_manle_module.f90
!> Unit test suite for ManLe module (manle_module.F90)
!> @details Unit tests for all public subroutines in manle_module.

module mod_test_manle_module
  use asserts
  use manle_module
  use anwil
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(7)

    all_tests(1) = test_case("test_center_module", test_center_module)
    all_tests(2) = test_case("test_compute_svd_randomized_f42", test_compute_svd_randomized_f42)
    all_tests(3) = test_case("test_qr_orthonormalize_f42", test_qr_orthonormalize_f42)
    all_tests(4) = test_case("test_compute_svd_randomized_f42_large", test_compute_svd_randomized_f42_large)
    all_tests(5) = test_case("test_project_to_subspace", test_project_to_subspace)
    all_tests(6) = test_case("test_project_to_subspace_k1", test_project_to_subspace_k1)
    all_tests(7) = test_case("test_smooth_vectors_gaussian_adaptive_mode0", test_smooth_vectors_gaussian_adaptive_mode0)

  end function get_all_tests

  subroutine run_all_tests_manle()
    type(test_case) :: all_tests(7)
    integer(int32) :: i

    all_tests = get_all_tests()

    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All ManLe module tests passed successfully."
  end subroutine run_all_tests_manle

  subroutine run_named_tests_manle(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(7)
    integer(int32) :: i, j
    logical :: found
    
    all_tests = get_all_tests()
    
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_manle

  !-------------------
  ! Test center_module
  !-------------------
  subroutine test_center_module()
    integer(int32), parameter :: n_points = 5, n_dims = 3
    real(real64) :: data(n_dims, n_points), centered_data(n_dims, n_points)
    integer(int32) :: ierr

    ! Initialize data
    data = reshape([1.0, 2.0, 3.0, 4.0, 5.0, &
                    6.0, 7.0, 8.0, 9.0, 10.0, &
                    11.0, 12.0, 13.0, 14.0, 15.0], [n_dims, n_points])

    ! Call center_module
    call center_module(data, n_points, n_dims, centered_data, ierr)

    ! Check for no errors
    call assert_equal_int(ierr, 0, "center_module: no error")

    ! Check that the mean of each dimension is approximately zero
    call assert_true(all(abs(sum(centered_data, dim=2)) < 1.0e-12), "center_module: mean centered")

    ! Check specific values
    call assert_true(abs(centered_data(1, 1) - (-6.0)) < 1.0e-12, "center_module: value check (1,1)")
    call assert_true(abs(centered_data(2, 5) - ( 6.0)) < 1.0e-12, "center_module: value check (2,5)")

    ! Additional test: constant matrix
    data = 5.0d0
    call center_module(data, n_points, n_dims, centered_data, ierr)
    call assert_equal_int(ierr, 0, "center_module: constant matrix no error")
    call assert_true(all(abs(centered_data) < 1.0e-12), "center_module: constant matrix centered to zero")

    ! Additional test: zero matrix
    data = 0.0d0
    call center_module(data, n_points, n_dims, centered_data, ierr)
    call assert_equal_int(ierr, 0, "center_module: zero matrix no error")
    call assert_true(all(abs(centered_data) < 1.0e-12), "center_module: zero matrix remains zero")

  end subroutine test_center_module

  subroutine test_compute_svd_randomized_f42()
    use iso_fortran_env, only: real64, int32
    implicit none

    integer(int32), parameter :: n_dims = 6
    integer(int32), parameter :: n_points = 20
    integer(int32), parameter :: k = 2
    integer(int32), parameter :: p = 4
    integer(int32), parameter :: l = k + p

    real(real64) :: M(n_dims, n_points)
    real(real64) :: Omega(n_points, l)
    real(real64) :: Y(n_dims, l)
    real(real64) :: Q(n_dims, l)
    real(real64) :: B(l, n_points)
    real(real64) :: Stmp(l)
    real(real64) :: Utmp(l, l)
    real(real64) :: tau(l)
    real(real64) :: work(512)

    real(real64) :: S(k)
    real(real64) :: U(n_dims, k)

    real(real64) :: Icheck(k,k)
    real(real64) :: v(n_dims), vproj(n_dims)
    real(real64) :: S_full(n_dims), U_full(n_dims, n_dims), VT_full(n_points, n_points)
    real(real64) :: M_copy(n_dims, n_points)


    integer(int32) :: i, j, ierr

    ! -------------------------
    ! Test Data
    ! -------------------------
    integer :: seed(8)
    seed = [123456, 654321, 111111, 222222, 333333, 444444, 555555, 666666]
    call random_seed(put=seed)

    call random_number(M)


    ! -------------------------
    ! Randomized SVD
    ! -------------------------
    call compute_svd_randomized_f42( &
      M, n_points, n_dims, k, p, &
      Omega, Y, Q, B, &
      Stmp, Utmp, tau, &
      S, U, &
      work, size(work), ierr)

    call assert_equal_int(ierr, 0, "compute_svd_randomized_f42: ierr == 0")

    ! -------------------------
    ! Test 1: Singular values are positive
    ! -------------------------
    do i = 1, k
      call assert_true(S(i) >= 0.0d0, "randomized_svd: singular value >= 0")
    end do

    ! -------------------------
    ! Test 2: U orthonormal → UᵀU = I
    ! -------------------------
    Icheck = matmul(transpose(U), U)

    do i = 1, k
      do j = 1, k
        if (i == j) then
          call assert_true(abs(Icheck(i,j) - 1.0d0) < 1.0d-10, &
              "randomized_svd: diag(U^T U) ~ 1")
        else
          call assert_true(abs(Icheck(i,j)) < 1.0d-10, &
              "randomized_svd: offdiag(U^T U) ~ 0")
        end if
      end do
    end do

    ! -------------------------
    ! Test 3: Projection lies in span(U)
    ! -------------------------
    v = M(:,1)
    vproj = 0.0d0
    do i = 1, k
      vproj = vproj + dot_product(v, U(:,i)) * U(:,i)
    end do

    call assert_true(norm2(vproj) > 0.0d0, &
        "randomized_svd: projected vector non-zero")

    ! -------------------------
    ! Test 4 (DEBUG): Compare subspace with full SVD + prints
    ! -------------------------

    M_copy = M   ! dgesvd sobreescribe la matriz

    call dgesvd('S', 'S', n_dims, n_points, M_copy, n_dims, &
                S_full, U_full, n_dims, VT_full, n_points, &
                work, size(work), ierr)

    call assert_equal_int(ierr, 0, "full_svd: ierr == 0")

    ! Producto entre subespacios
    Icheck = matmul(transpose(U), U_full(:,1:k))

    ! Ahora sí validamos con tolerancia razonable
    Icheck = abs(Icheck)

    do i = 1, k
      do j = 1, k
        if (i == j) then
          call assert_true(abs(Icheck(i,j) - 1.0d0) < 1.0d-3, &
              "randomized_svd: principal subspace matches full SVD (diag)")
        else
          call assert_true(Icheck(i,j) < 1.0d-3, &
              "randomized_svd: principal subspace matches full SVD (offdiag)")
        end if
      end do
    end do



  end subroutine test_compute_svd_randomized_f42

  subroutine test_compute_svd_randomized_f42_large()
  use iso_fortran_env, only: real64, int32
  implicit none

  integer(int32), parameter :: n_dims   = 50
  integer(int32), parameter :: n_points = 200
  integer(int32), parameter :: k        = 2
  integer(int32), parameter :: p        = 5
  integer(int32), parameter :: l        = min(k + p, n_dims)

  real(real64) :: M(n_dims, n_points)
  real(real64) :: M_copy(n_dims, n_points)

  real(real64) :: Omega(n_points, l)
  real(real64) :: Y(n_dims, l)
  real(real64) :: Q(n_dims, l)
  real(real64) :: B(l, n_points)
  real(real64) :: Stmp(l)
  real(real64) :: Utmp(l, l)
  real(real64) :: tau(l)
  real(real64) :: work(8000)

  real(real64) :: S(k)
  real(real64) :: U(n_dims, k)

  real(real64) :: S_full(n_dims)
  real(real64) :: U_full(n_dims, n_dims)
  real(real64) :: VT_full(n_points, n_points)

  real(real64) :: Icheck(k, k)

  integer(int32) :: i, j, ierr

  ! -------------------------
  ! Datos de prueba
  ! -------------------------
  integer :: seed(8)
  seed = [123456, 654321, 111111, 222222, 333333, 444444, 555555, 666666]
  call random_seed(put=seed)

  call random_number(M)

  M_copy = M

  ! -------------------------
  ! Randomized SVD
  ! -------------------------
  call compute_svd_randomized_f42( &
       M, n_points, n_dims, k, p, &
       Omega, Y, Q, B, &
       Stmp, Utmp, tau, &
       S, U, &
       work, size(work), ierr)

  call assert_equal_int(ierr, 0, "large_randomized_svd: ierr == 0")

  ! -------------------------
  ! Full SVD para comparación
  ! -------------------------
  call dgesvd('S', 'S', n_dims, n_points, M_copy, n_dims, &
              S_full, U_full, n_dims, VT_full, n_points, &
              work, size(work), ierr)
  call assert_equal_int(ierr, 0, "large_full_svd: ierr == 0")

  ! -------------------------
  ! DEBUG: imprime comparación
  ! -------------------------
  print *, "========================================"
  print *, "DEBUG Large Randomized SVD vs Full SVD"
  print *, "n_dims=", n_dims, " n_points=", n_points, " k=", k, " p=", p, " l=", l
  print *, "Singular values (randomized, first k):"
  do i = 1, k
    print *, "  S_rand(", i, ") = ", S(i)
  end do
  print *, "Singular values (full, first k):"
  do i = 1, k
    print *, "  S_full(", i, ") = ", S_full(i)
  end do

  Icheck = matmul(transpose(U), U_full(:,1:k))

  print *, "U_rand^T * U_full(:,1:k) (raw):"
  do i = 1, k
    print *, Icheck(i,1:k)
  end do

  Icheck = abs(Icheck)
  print *, "U_rand^T * U_full(:,1:k) (abs):"
  do i = 1, k
    print *, Icheck(i,1:k)
  end do
  print *, "========================================"

  ! -------------------------
  ! Checks razonables (aproximado, no exacto)
  ! -------------------------

  ! 1) valores singulares no negativos
  do i = 1, k
    call assert_true(S(i) >= 0.0_real64, "large_randomized_svd: S(i) >= 0")
  end do

  ! después de Icheck = abs(Icheck)

    ! Eje 1: muy bien alineado
    call assert_true(Icheck(1,1) > 0.99_real64, &
        "large_randomized_svd: first axis strongly aligned")
    call assert_true(Icheck(1,2) < 0.1_real64, &
        "large_randomized_svd: limited leakage into second axis")


end subroutine test_compute_svd_randomized_f42_large


  subroutine test_qr_orthonormalize_f42()
  use iso_fortran_env, only: real64, int32
  implicit none

  integer(int32), parameter :: m = 5, n = 3
  real(real64) :: A(m,n), Q(m,n)
  real(real64) :: tau(n)
  real(real64) :: work(128)
  real(real64) :: Icheck(n,n)
  integer(int32) :: i, j, info

  ! -------------------------
  ! Test 1: General matrix
  ! -------------------------
  A = reshape([ &
    1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, &
    2.0d0, 1.0d0, 2.0d0, 1.0d0, 2.0d0, &
    3.0d0, 2.0d0, 1.0d0, 2.0d0, 1.0d0  &
  ], [m,n])

  call qr_orthonormalize_f42(A, Q, m, n, tau, work, size(work), info)
  call assert_equal_int(info, 0, "qr_orthonormalize_f42: info == 0")

  Icheck = matmul(transpose(Q), Q)
  do i = 1, n
    do j = 1, n
      if (i == j) then
        call assert_true(abs(Icheck(i,j) - 1.0d0) < 1.0d-10, &
            "qr_orthonormalize_f42: diag(Q^T Q) ~ 1")
      else
        call assert_true(abs(Icheck(i,j)) < 1.0d-10, &
            "qr_orthonormalize_f42: offdiag(Q^T Q) ~ 0")
      end if
    end do
  end do

  ! -------------------------
  ! Test 2: Rank-deficient matrix
  ! -------------------------
  A = 0.0d0
  A(1,1) = 1.0d0
  call qr_orthonormalize_f42(A, Q, m, n, tau, work, size(work), info)
  call assert_equal_int(info, 0, "qr_orthonormalize_f42: rank-deficient matrix info == 0")

  Icheck = matmul(transpose(Q), Q)
  do i = 1, n
    do j = 1, n
      if (i == j) then
        call assert_true(abs(Icheck(i,j) - 1.0d0) < 1.0d-10, &
            "qr_orthonormalize_f42: rank-deficient diag(Q^T Q) ~ 1")
      else
        call assert_true(abs(Icheck(i,j)) < 1.0d-10, &
            "qr_orthonormalize_f42: rank-deficient offdiag(Q^T Q) ~ 0")
      end if
    end do
  end do

end subroutine test_qr_orthonormalize_f42
subroutine test_project_to_subspace()
  use iso_fortran_env, only: real64, int32
  implicit none

  integer(int32), parameter :: n_points = 4, n_dims = 3, k = 2
  real(real64) :: data(n_dims, n_points)
  real(real64) :: Uk(n_dims, k)
  real(real64) :: projected_data(n_dims, n_points)
  real(real64) :: residual(n_dims, n_points)
  real(real64) :: check(k, n_points)
  integer(int32) :: ierr, i, j

  print *, "========================================"
  print *, "DEBUG test_project_to_subspace"
  print *, "n_dims=", n_dims, " n_points=", n_points, " k=", k

  ! -------------------------
  ! Datos de prueba
  ! -------------------------
  data = reshape([ &
      1.0d0, 2.0d0, 3.0d0, 4.0d0, &
      5.0d0, 6.0d0, 7.0d0, 8.0d0, &
      9.0d0,10.0d0,11.0d0,12.0d0  &
    ], [n_dims, n_points])

  print *, "DATA:"
  do j = 1, n_points
    print *, "v(",j,") =", data(:,j)
  end do

  ! ✅ Base ortonormal exacta
  Uk = reshape([ &
  1.0d0, 0.0d0, 0.0d0, &  ! col 1
  0.0d0, 1.0d0, 0.0d0   &  ! col 2
], [n_dims, k])


  print *, "Uk:"
  do i = 1, k
    print *, "u(",i,") =", Uk(:,i)
  end do

  ! -------------------------
  ! PROYECCIÓN
  ! -------------------------
  call project_to_subspace(data, Uk, k, n_points, n_dims, projected_data, ierr)
  print *, "ierr =", ierr

  print *, "PROJECTED:"
  do j = 1, n_points
    print *, "p(",j,") =", projected_data(:,j)
  end do

  ! -------------------------
  ! RESIDUO
  ! -------------------------
  residual = data - projected_data

  print *, "RESIDUAL:"
  do j = 1, n_points
    print *, "r(",j,") =", residual(:,j)
  end do

  ! -------------------------
  ! CHECK = Uk^T * residual
  ! -------------------------
  check = matmul(transpose(Uk), residual)

  print *, "CHECK = Uk^T * residual:"
  do j = 1, n_points
    print *, "check(",j,") =", check(:,j)
  end do

  print *, "========================================"

  ! -------------------------
  ! AHORA SÍ EL ASSERT
  ! -------------------------
  do j = 1, n_points
    call assert_true(all(abs(check(:,j)) < 1.0d-12), &
        "project_to_subspace: residual orthogonal to Uk")
  end do

end subroutine test_project_to_subspace


subroutine test_project_to_subspace_k1()
  use iso_fortran_env, only: real64, int32
  implicit none

  integer(int32), parameter :: n_points = 4, n_dims = 3, k = 1
  real(real64) :: data(n_dims, n_points)
  real(real64) :: Uk(n_dims, k)
  real(real64) :: projected_data(n_dims, n_points)
  real(real64) :: expected_projected_data(n_dims, n_points)
  integer(int32) :: ierr

  integer(int32) :: i, j

  print *, "========================================"
  print *, "DEBUG test_project_to_subspace_k1"
  print *, "n_dims=", n_dims, " n_points=", n_points, " k=", k

  ! -------------------------
  ! Datos
  ! -------------------------
  data = reshape([ &
      1.0d0, 2.0d0, 3.0d0, 4.0d0, &
      5.0d0, 6.0d0, 7.0d0, 8.0d0, &
      9.0d0,10.0d0,11.0d0,12.0d0  &
    ], [n_dims, n_points])

  print *, "DATA:"
  do j = 1, n_points
    print *, "v(", j, ") =", data(:,j)
  end do

  ! u1 = eje X
  Uk(:,1) = [1.0d0, 0.0d0, 0.0d0]

  print *, "Uk:"
  do i = 1, k
    print *, "u(", i, ") =", Uk(:,i)
  end do

  ! Proyección esperada: (x,0,0)
  expected_projected_data = reshape([ &
    1.0d0, 0.0d0, 0.0d0, &  ! col 1
    4.0d0, 0.0d0, 0.0d0, &  ! col 2
    7.0d0, 0.0d0, 0.0d0, &  ! col 3
    10.0d0, 0.0d0, 0.0d0   &  ! col 4
    ], [n_dims, n_points])


  ! -------------------------
  ! Llamada a project_to_subspace
  ! -------------------------
  call project_to_subspace(data, Uk, k, n_points, n_dims, projected_data, ierr)

  print *, "ierr =", ierr

  print *, "PROJECTED:"
  do j = 1, n_points
    print *, "p(", j, ") =", projected_data(:,j)
  end do

  print *, "EXPECTED:"
  do j = 1, n_points
    print *, "e(", j, ") =", expected_projected_data(:,j)
  end do

  print *, "DIFF:"
  do j = 1, n_points
    print *, "d(", j, ") =", projected_data(:,j) - expected_projected_data(:,j)
  end do

  print *, "========================================"

  call assert_equal_int(ierr, 0, "project_to_subspace_k1: no error")

  call assert_true(all(abs(projected_data - expected_projected_data) < 1.0d-12), &
      "project_to_subspace_k1: projected data matches expected")

end subroutine test_project_to_subspace_k1

subroutine test_smooth_vectors_gaussian_adaptive_mode0()
    use iso_fortran_env, only: real64, int32
    implicit none

    ! Dimension parameters
    integer(int32), parameter :: n_points = 5, n_dims = 2
    
    ! Arguments for the subroutine
    real(real64)    :: coords(n_dims, n_points)
    real(real64)    :: vecs(n_dims, n_points)
    real(real64)    :: smoothed(n_dims, n_points)
    real(real64)    :: sigma_raw(1, n_points)      ! Required by the new signature
    integer(int32)  :: kd_indices(n_points)
    integer(int32)  :: dimension_order(n_dims)
    integer(int32)  :: workspace(n_points)
    real(real64)    :: value_buffer(n_points)
    integer(int32)  :: permutation(n_points)
    integer(int32)  :: left_stack(n_points)
    integer(int32)  :: right_stack(n_points)
    real(real64)    :: sd_arr(n_dims, n_points)    ! Adjusted to (n_vector_dims, n_points)
    integer(int32)  :: ierr
    
    ! Control parameters
    integer(int32)  :: k_neighbors = 5
    integer(int32)  :: k_neighbors_sigma = 3       
    integer(int32)  :: kernel_type = 1             ! 1: Gaussian
    integer(int32)  :: anisotropy_mode = 0         ! Mode 0: Isotropic
    real(real64)    :: anisotropy_factor = 1.0_real64

    ! Test data and expected results
    real(real64)    :: expected_smoothed(n_dims, n_points)

    ! Initialize coordinates and vectors
    ! Using reshape to ensure specific mapping [dim, point]
    coords = reshape([1.0, 1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0, 5.0], [n_dims, n_points])
    vecs = reshape([10.0, 5.0, 20.0, 10.0, 30.0, 15.0, 40.0, 20.0, 50.0, 25.0], [n_dims, n_points])

    ! Expected smoothed values from original test logic
    expected_smoothed = reshape([23.861929093451685, 25.231420749740611, &
                                 21.907633336401815, 22.135772559620364, &
                                 16.722814001864784, 27.451268882624056, &
                                 18.932113720189815, 24.585696583374318, &
                                 17.731420749740607, 23.995634110957173], [n_dims, n_points])

    ! Call the updated subroutine
    ! Optional arguments (U, singular_values, k_intrinsic) are omitted for Mode 0
    call anwil_smooth_sigma( &
        coords              = coords, &
        vectors             = vecs, &
        smoothed            = smoothed, &
        n_coord_dims        = n_dims, &
        n_vector_dims       = n_dims, &
        n_points            = n_points, &
        k_neighbors         = k_neighbors, &
        kd_indices          = kd_indices, &
        dimension_order     = dimension_order, &
        workspace           = workspace, &
        value_buffer        = value_buffer, &
        permutation         = permutation, &
        left_stack          = left_stack, &
        right_stack         = right_stack, &
        sigma_raw           = sigma_raw, &
        sd_arr              = sd_arr, &
        anisotropy_mode     = anisotropy_mode, &
        anisotropy_factor   = anisotropy_factor, &
        k_neighbors_sigma   = k_neighbors_sigma, &
        kernel_type         = kernel_type, &
        ierr                = ierr)

    ! Assertions
    call assert_equal_int(ierr, 0, "test_smooth_vectors_mode0: ierr should be 0")

    call assert_true(all(abs(smoothed - expected_smoothed) < 1.0e-6), &
        "test_smooth_vectors_mode0: result should match expected values")

  end subroutine test_smooth_vectors_gaussian_adaptive_mode0

end module mod_test_manle_module
