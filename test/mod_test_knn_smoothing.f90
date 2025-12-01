!> @file mod_test_knn_smoothing.f90
!> Unit test suite for KNN smoothing (knn_smoothing.F90)
!> @details Unit tests for all public subroutines in knn_smoothing except kd_request_neighbors.

module mod_test_knn_smoothing
  use asserts
  use knn_smoothing
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
    type(test_case) :: all_tests(13)

    all_tests(1) = test_case("test_validate_k_params", test_validate_k_params)
    all_tests(2) = test_case("test_compute_sigma_for_k", test_compute_sigma_for_k)
    all_tests(3) = test_case("test_smooth_point_with_k", test_smooth_point_with_k)
    all_tests(4) = test_case("test_compute_second_derivative", test_compute_second_derivative)
    all_tests(5) = test_case("test_compute_curvature_for_k", test_compute_curvature_for_k)
    all_tests(6) = test_case("test_compute_normalized_roughness", test_compute_normalized_roughness)
    all_tests(7) = test_case("test_apply_elbow_rule", test_apply_elbow_rule)
    all_tests(8) = test_case("test_smooth_multiple_k", test_smooth_multiple_k)
    all_tests(9) = test_case("test_compute_global_roughness", test_compute_global_roughness)
    all_tests(10) = test_case("test_smooth_vectors_gaussian_adaptive", test_smooth_vectors_gaussian_adaptive)
    all_tests(11) = test_case("test_smooth_multiple_k_extra", test_smooth_multiple_k)
    all_tests(12) = test_case("test_compute_global_roughness_extra", test_compute_global_roughness_extra)
    all_tests(13) = test_case("test_smooth_vectors_gaussian_adaptive_extra", test_smooth_vectors_gaussian_adaptive_extra)
  end function get_all_tests

  subroutine run_all_tests_knn_smoothing()
    type(test_case) :: all_tests(13)
    integer(int32) :: i

    all_tests = get_all_tests()

    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All KNN smoothing tests passed successfully."
  end subroutine run_all_tests_knn_smoothing

  subroutine run_named_tests_knn_smoothing(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(13)
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
  end subroutine run_named_tests_knn_smoothing

  !-------------------
  ! Test validate_k_params
  !-------------------
  subroutine test_validate_k_params()
    integer(int32) :: ierr, k_grid(3)
    ! Caso válido
    k_grid = (/2, 3, 4/)
    call validate_k_params(2, 4, k_grid, 3, 10, ierr)
    call assert_equal_int(ierr, 0, 'validate_k_params: valid input')
    ! k_grid no creciente
    k_grid = (/2, 2, 4/)
    call validate_k_params(2, 4, k_grid, 3, 10, ierr)
    call assert_true(ierr /= 0, 'validate_k_params: non-increasing k_grid')
    ! k_grid(1) /= k_min
    k_grid = (/1, 2, 3/)
    call validate_k_params(2, 4, k_grid, 3, 10, ierr)
    call assert_true(ierr /= 0, 'validate_k_params: k_grid(1) /= k_min')
    ! k_min > k_max
    call validate_k_params(2, 1, (/2/), 1, 10, ierr)
    call assert_true(ierr /= 0, 'validate_k_params: k_min > k_max')
    ! n_k < 1
    call validate_k_params(2, 4, (/2,3,4/), 0, 10, ierr)
    call assert_true(ierr /= 0, 'validate_k_params: n_k < 1')
    ! k_grid fuera de rango
    k_grid = (/2, 3, 11/)
    call validate_k_params(2, 10, k_grid, 3, 10, ierr)
    call assert_true(ierr /= 0, 'validate_k_params: k_grid(kk) > k_max')
  end subroutine

  !-------------------
  ! Test compute_sigma_for_k
  !-------------------
  subroutine test_compute_sigma_for_k()
    real(real64) :: distances(3), sigma
    ! Caso normal
    distances = (/1.0, 2.0, 3.0/)
    call compute_sigma_for_k(distances, 2, 1.5_real64, sigma)
    call assert_true(abs(sigma - 3.0) < 1e-12, 'compute_sigma_for_k: normal case')
    ! sigma_factor = 0
    call compute_sigma_for_k(distances, 2, 0.0_real64, sigma)
    call assert_true(sigma > 0.0, 'compute_sigma_for_k: fallback for small sigma')
    ! k fuera de rango
    call compute_sigma_for_k(distances, 0, 1.0_real64, sigma)
    call assert_true(sigma == 0.0, 'compute_sigma_for_k: k=0')
    ! distancias negativas
    distances = (/-1.0, -2.0, -3.0/)
    call compute_sigma_for_k(distances, 2, 1.0_real64, sigma)
    call assert_true(sigma == 2.0, 'compute_sigma_for_k: negative distances (abs)')
  end subroutine

  !-------------------
  ! Test smooth_point_with_k
  !-------------------
  subroutine test_smooth_point_with_k()
    real(real64) :: vectors(1,3), distances(3), result(1)
    integer(int32) :: neighbors(3)
    ! Promedio simple
    vectors = reshape((/1.0, 2.0, 3.0/), (/1,3/))
    distances = (/0.0, 1.0, 2.0/)
    neighbors = (/1,2,3/)
    call smooth_point_with_k(neighbors, distances, vectors, 1, 3, 3, 1.0_real64, result)
    call assert_true(abs(result(1) - sum(vectors)/3.0) < 1.0, 'smooth_point_with_k: basic average')
    ! k=0
    call smooth_point_with_k(neighbors, distances, vectors, 1, 3, 0, 1.0_real64, result)
    call assert_true(all(result == 0.0), 'smooth_point_with_k: k=0 returns zero')
    ! k=1
    call smooth_point_with_k(neighbors, distances, vectors, 1, 3, 1, 1.0_real64, result)
    call assert_true(abs(result(1) - vectors(1,1)) < 1e-12, 'smooth_point_with_k: k=1')
    ! pesos extremos (sigma muy pequeño)
    call smooth_point_with_k(neighbors, distances, vectors, 1, 3, 3, 1e-12_real64, result)
    call assert_true(abs(result(1) - vectors(1,1)) < 1.0, 'smooth_point_with_k: small sigma')
    ! vecinos repetidos
    neighbors = (/2,2,2/)
    call smooth_point_with_k(neighbors, distances, vectors, 1, 3, 3, 1.0_real64, result)
    call assert_true(abs(result(1) - vectors(1,2)) < 1e-12, 'smooth_point_with_k: repeated neighbors')
  end subroutine

  !-------------------
  ! Test compute_second_derivative
  !-------------------
  subroutine test_compute_second_derivative()
    real(real64) :: current(2), prev(2), prev2(2), a(2)
    ! Caso lineal
    current = (/3.0, 5.0/); prev = (/2.0, 3.0/); prev2 = (/1.0, 1.0/)
    call compute_second_derivative(current, prev, prev2, a)
    call assert_true(all(abs(a - (/0.0, 0.0/)) < 1e-12), 'compute_second_derivative: linear')
    ! Curvatura positiva
    current = (/2.0, 8.0/); prev = (/1.0, 3.0/); prev2 = (/0.0, 0.0/)
    call compute_second_derivative(current, prev, prev2, a)
    call assert_true(all(a >= 0.0), 'compute_second_derivative: positive curvature')
    ! Curvatura negativa
    current = (/0.0, 0.0/); prev = (/1.0, 3.0/); prev2 = (/2.0, 0.0/)
    call compute_second_derivative(current, prev, prev2, a)
    call assert_true(any(a < 0.0) .and. all(a <= 0.0), 'compute_second_derivative: negative curvature')
    ! Arrays de tamaño 1
    call compute_second_derivative([1.0_real64], [0.0_real64], [1.0_real64], a(1:1))
    call assert_true(abs(a(1) - 2.0) < 1e-12, 'compute_second_derivative: size 1')
    ! Arrays iguales
    call compute_second_derivative([1.0_real64, 1.0_real64], [1.0_real64, 1.0_real64], [1.0_real64, 1.0_real64], a)
    call assert_true(all(abs(a) < 1e-12), 'compute_second_derivative: all equal')
  end subroutine

  !-------------------
  ! Test compute_curvature_for_k
  !-------------------
    subroutine test_compute_curvature_for_k()
        real(real64) :: t_hat_k(2), t_prev(2), t_prev2(2), C

        ! Case 1: curvatura positiva (convexa)
        t_prev2 = (/0.0, 0.0/)
        t_prev  = (/1.0, 2.0/)
        t_hat_k = (/3.0, 5.0/)
        call compute_curvature_for_k(t_hat_k, t_prev, t_prev2, C)
        call assert_true(C > 0.0_real64, 'convex curvature')

        ! Case 2: lineal → curvatura cero
        t_prev2 = (/0.0, 0.0/)
        t_prev  = (/1.0, 2.0/)
        t_hat_k = (/2.0, 4.0/)
        call compute_curvature_for_k(t_hat_k, t_prev, t_prev2, C)
        call assert_true(abs(C) < 1e-12_real64, 'linear curvature = 0')

        ! Case 3: cóncava
        t_prev2 = (/0.0, 0.0/)
        t_prev  = (/2.0, 4.0/)
        t_hat_k = (/1.0, 2.0/)
        call compute_curvature_for_k(t_hat_k, t_prev, t_prev2, C)
        call assert_true(C > 0.0_real64, 'concave curvature')

        ! Case 4: tamaño 1, con curvatura real
        call compute_curvature_for_k([4.0_real64], [1.0_real64], [0.0_real64], C)
        call assert_true(abs(C - 4.0_real64) < 1e-12_real64, 'size=1 curvature')

    end subroutine




  !-------------------
  ! Test compute_normalized_roughness
  !-------------------
  subroutine test_compute_normalized_roughness()
    real(real64) :: C(3), Rtilde(3)
    ! Caso normal
    C = (/2.0, 1.0, 0.5/)
    call compute_normalized_roughness(C, 3, Rtilde)
    call assert_true(abs(Rtilde(1) - 1.0) < 1e-12, 'compute_normalized_roughness: first is 1')
    call assert_true(Rtilde(2) == 0.5 .and. Rtilde(3) == 0.25, 'compute_normalized_roughness: ratios')
    ! Todos ceros
    C = 0.0
    call compute_normalized_roughness(C, 3, Rtilde)
    call assert_true(all(Rtilde == 1.0), 'compute_normalized_roughness: all zero')
    ! C(1) negativo
    C = (/-2.0, -1.0, -0.5/)
    call compute_normalized_roughness(C, 3, Rtilde)
    call assert_true(Rtilde(1) == 1.0, 'compute_normalized_roughness: negative first')
    ! Arrays de tamaño 1
    call compute_normalized_roughness([2.0_real64], 1, Rtilde(1:1))
    call assert_true(abs(Rtilde(1) - 1.0) < 1e-12, 'compute_normalized_roughness: size 1')
  end subroutine

  !-------------------
  ! Test apply_elbow_rule
  !-------------------
subroutine test_apply_elbow_rule()
  real(real64) :: Rtilde(3)
  integer(int32) :: best_k_idx

  ! Caso normal
  Rtilde = (/1.0, 0.7, 0.6/)
  call apply_elbow_rule(Rtilde, 3, 0.05_real64, best_k_idx)
  call assert_true(best_k_idx == 3, 'apply_elbow_rule: normal case')

  ! Epsilon grande
  Rtilde = (/1.0, 0.9, 0.8/)
  call apply_elbow_rule(Rtilde, 3, 0.5_real64, best_k_idx)
  call assert_true(best_k_idx == 1, 'apply_elbow_rule: large epsilon')

  ! Todos iguales
  Rtilde = (/1.0, 1.0, 1.0/)
  call apply_elbow_rule(Rtilde, 3, 0.01_real64, best_k_idx)
  call assert_true(best_k_idx == 1, 'apply_elbow_rule: all equal')

  ! Arrays de tamaño 1
  call apply_elbow_rule([1.0_real64], 1, 0.1_real64, best_k_idx)
  call assert_true(best_k_idx == 1, 'apply_elbow_rule: size 1')
end subroutine



  !-------------------
  ! Test smooth_multiple_k (exhaustivo)
  !-------------------
  subroutine test_smooth_multiple_k()
    integer(int32), parameter :: n_vec=1, n_pts=5, n_k=3
    real(real64) :: vectors(n_vec, n_pts), smoothed(n_vec, n_pts)
    integer(int32) :: neighbors(n_pts), k_grid(n_k)
    real(real64) :: distances(n_pts), t_hat_k(n_vec, n_k), C(n_k)
    integer(int32) :: i
    vectors = reshape([1.0, 2.0, 3.0, 4.0, 5.0], [n_vec, n_pts])
    smoothed = 0.0
    neighbors = [1,2,3,4,5]
    distances = [0.0,1.0,2.0,3.0,4.0]
    k_grid = [2,3,5]
    ! Llenar smoothed con valores crecientes para simular historial
    do i=1,n_pts
      smoothed(:,i) = i
    end do
    call smooth_multiple_k(neighbors, distances, vectors, smoothed, n_vec, n_pts, 5, k_grid, n_k, 1.0_real64, t_hat_k, C)
    call assert_true(all(C >= 0.0), 'smooth_multiple_k: curvatures non-negative')
    call assert_true(all(t_hat_k(:,1) /= 0.0), 'smooth_multiple_k: output nonzero')
  end subroutine

  !-------------------
  ! Test compute_global_roughness (exhaustivo)
  !-------------------
  subroutine test_compute_global_roughness()
    integer(int32), parameter :: n_pts=5, n_k=3
    real(real64) :: C_all(n_pts, n_k), R(n_k)
    integer(int32) :: i, k
    C_all = 0.0
    ! Simula curvaturas crecientes
    do k=1,n_k
      do i=1,n_pts
        C_all(i,k) = real(i*k, real64)
      end do
    end do
    call compute_global_roughness(C_all, n_pts, n_k, R)
    call assert_true(all(R > 0.0), 'compute_global_roughness: all positive')
    ! Caso todo cero
    C_all = 0.0
    call compute_global_roughness(C_all, n_pts, n_k, R)
    call assert_true(all(R == 0.0), 'compute_global_roughness: all zero')
  end subroutine

  !-------------------
  ! Test smooth_vectors_gaussian_adaptive (básico)
  !-------------------
  subroutine test_smooth_vectors_gaussian_adaptive()
    integer(int32), parameter :: n_coord=1, n_vec=1, n_pts=5, n_k=2
    real(real64) :: coords(n_coord, n_pts), vectors(n_vec, n_pts), smoothed(n_vec, n_pts)
    integer(int32) :: kd_indices(n_pts), dimension_order(n_coord)
    integer(int32) :: workspace(n_pts), permutation(n_pts), left_stack(n_pts), right_stack(n_pts)
    real(real64) :: value_buffer(n_pts), distances(n_pts)
    integer(int32) :: neighbors(n_pts), k_grid(n_k), ierr
    real(real64) :: sigma_factor, epsilon
    integer(int32) :: i
    coords = reshape([1.0,2.0,3.0,4.0,5.0], [n_coord, n_pts])
    vectors = reshape([10.0,20.0,30.0,40.0,50.0], [n_vec, n_pts])
    k_grid = [1,3]
    sigma_factor = 1.0
    epsilon = 0.1
    call smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, &
      n_coord, n_vec, n_pts, 1, 3, 1, sigma_factor, epsilon, &
      kd_indices, dimension_order, neighbors, distances, workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    call assert_equal_int(ierr, 0, 'smooth_vectors_gaussian_adaptive: no error')
    call assert_true(all(smoothed(:,1:2) /= 0.0), 'smooth_vectors_gaussian_adaptive: output nonzero')
  end subroutine

  !-------------------
  ! Test smooth_multiple_k (más casos)
  !-------------------
  subroutine test_smooth_multiple_k_extra()
    integer(int32), parameter :: n_vec=2, n_pts=4, n_k=2
    real(real64) :: vectors(n_vec, n_pts), smoothed(n_vec, n_pts)
    integer(int32) :: neighbors(n_pts), k_grid(n_k)
    real(real64) :: distances(n_pts), t_hat_k(n_vec, n_k), C(n_k)
    ! Caso multidimensional
    vectors = reshape([1.0,2.0,3.0,4.0, 10.0,20.0,30.0,40.0], [n_vec, n_pts])
    smoothed = 0.0
    neighbors = [1,2,3,4]
    distances = [0.0,1.0,2.0,3.0]
    k_grid = [2,4]
    call smooth_multiple_k(neighbors, distances, vectors, smoothed, n_vec, n_pts, 4, k_grid, n_k, 1.0_real64, t_hat_k, C)
    call assert_true(all(C >= 0.0), 'smooth_multiple_k_extra: curvatures non-negative')
    call assert_true(all(t_hat_k(:,1) /= 0.0), 'smooth_multiple_k_extra: output nonzero')
    ! Caso con smoothed previo igual a t_hat_k
    vectors(:,:)  = 5.0_real64
  smoothed(:,:) = 5.0_real64
  neighbors(1:3) = [1,2,3]
  distances(1:3) = [0.0,1.0,2.0]
  k_grid = [1,3]

  call smooth_multiple_k(neighbors(1:3), distances(1:3), vectors, smoothed, &
                         n_vec, n_pts, 3, k_grid, n_k, 1.0_real64, t_hat_k, C)

  ! Como todo es 5.0, t_hat_k(:,kk) = 5.0, prev = 5.0, prev2 = 5.0 => C ≈ 0
  call assert_true(all(C < 1.0e-12_real64), 'smooth_multiple_k_extra: zero curvature for constant data')

  end subroutine

  !-------------------
  ! Test compute_global_roughness (más casos)
  !-------------------
  subroutine test_compute_global_roughness_extra()
    integer(int32), parameter :: n_pts=6, n_k=2
    real(real64) :: C_all(n_pts, n_k), R(n_k)
    ! Caso con valores negativos
    C_all = reshape([-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [n_pts, n_k])
    call compute_global_roughness(C_all, n_pts, n_k, R)
    call assert_true(R(1) < 0.0 .and. R(2) > 0.0, 'compute_global_roughness_extra: mixed sign')
    ! Caso n_pts pequeño
    C_all = 1.0
    call compute_global_roughness(C_all, 3, n_k, R)
    call assert_true(all(R == 2.0), 'compute_global_roughness_extra: n_pts=3')
  end subroutine

  !-------------------
  ! Test smooth_vectors_gaussian_adaptive (más casos)
  !-------------------
  subroutine test_smooth_vectors_gaussian_adaptive_extra()
    integer(int32), parameter :: n_coord=2, n_vec=2, n_pts=4, n_k=2
    real(real64) :: coords(n_coord, n_pts), vectors(n_vec, n_pts), smoothed(n_vec, n_pts)
    integer(int32) :: kd_indices(n_pts), dimension_order(n_coord)
    integer(int32) :: workspace(n_pts), permutation(n_pts), left_stack(n_pts), right_stack(n_pts)
    real(real64) :: value_buffer(n_pts), distances(n_pts)
    integer(int32) :: neighbors(n_pts), k_grid(n_k), ierr
    real(real64) :: sigma_factor, epsilon
    ! Caso multidimensional
    coords = reshape([1.0,2.0,3.0,4.0, 10.0,20.0,30.0,40.0], [n_coord, n_pts])
    vectors = reshape([1.0,2.0,3.0,4.0, 10.0,20.0,30.0,40.0], [n_vec, n_pts])
    k_grid = [1,2]
    sigma_factor = 1.0
    epsilon = 0.1
    call smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, n_coord, n_vec, n_pts, 1, 2, 1, sigma_factor, epsilon, &
      kd_indices, dimension_order, neighbors, distances, workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    call assert_equal_int(ierr, 0, 'smooth_vectors_gaussian_adaptive_extra: no error')
    call assert_true(all(smoothed /= 0.0), 'smooth_vectors_gaussian_adaptive_extra: output nonzero')
    ! Caso k_grid con un solo valor
    k_grid = [1,1]
    call smooth_vectors_gaussian_adaptive(coords, vectors, smoothed, n_coord, n_vec, n_pts, 1, 1, 1, sigma_factor, epsilon, &
      kd_indices, dimension_order, neighbors, distances, workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    call assert_equal_int(ierr, 0, 'smooth_vectors_gaussian_adaptive_extra: k_grid=1')
  end subroutine

end module mod_test_knn_smoothing
