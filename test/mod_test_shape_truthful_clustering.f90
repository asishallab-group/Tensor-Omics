! filepath: test/mod_test_shape_truthful_clustering.f90
!> Unit test suite for the shape_truthful_clustering module (STC, tangent-space variant).
module mod_test_shape_truthful_clustering
  use shape_truthful_clustering
  use tox_errors, only: is_ok, is_err
  use kd_tree, only: build_kd_index
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(43)

    all_tests(1)  = test_case("test_normal_error_basic", test_normal_error_basic)
    all_tests(2)  = test_case("test_normal_error_zero_tangent_dims", test_normal_error_zero_tangent_dims)
    all_tests(3)  = test_case("test_normal_error_all_tangent_dims", test_normal_error_all_tangent_dims)
    all_tests(4)  = test_case("test_normal_error_d_out_of_range", test_normal_error_d_out_of_range)
    all_tests(5)  = test_case("test_normal_error_negative_eigenvalue", test_normal_error_negative_eigenvalue)
    all_tests(6)  = test_case("test_normal_error_zero_dimensions", test_normal_error_zero_dimensions)
    all_tests(7)  = test_case("test_tangent_scales_basic", test_tangent_scales_basic)
    all_tests(8)  = test_case("test_tangent_scales_all_dims", test_tangent_scales_all_dims)
    all_tests(9)  = test_case("test_tangent_scales_zero_dims", test_tangent_scales_zero_dims)
    all_tests(10) = test_case("test_tangent_scales_d_out_of_range", test_tangent_scales_d_out_of_range)
    all_tests(11) = test_case("test_tangent_scales_negative_eigenvalue", test_tangent_scales_negative_eigenvalue)
    all_tests(12) = test_case("test_density_radius_alloc_default_percentile", test_density_radius_alloc_default_percentile)
    all_tests(13) = test_case("test_density_radius_alloc_custom_percentile", test_density_radius_alloc_custom_percentile)
    all_tests(14) = test_case("test_density_radius_preallocated", test_density_radius_preallocated)
    all_tests(15) = test_case("test_density_radius_invalid_percentile", test_density_radius_invalid_percentile)
    all_tests(16) = test_case("test_density_radius_zero_vectors", test_density_radius_zero_vectors)
    all_tests(17) = test_case("test_density_radius_zero_dimensions", test_density_radius_zero_dimensions)
    all_tests(18) = test_case("test_growth_radius_even_k", test_growth_radius_even_k)
    all_tests(19) = test_case("test_growth_radius_odd_k", test_growth_radius_odd_k)
    all_tests(20) = test_case("test_growth_radius_seed_index_out_of_range", test_growth_radius_seed_index_out_of_range)
    all_tests(21) = test_case("test_growth_radius_k_min_too_large", test_growth_radius_k_min_too_large)
    all_tests(22) = test_case("test_growth_radius_zero_dimensions", test_growth_radius_zero_dimensions)
    all_tests(23) = test_case("test_density_labels_basic", test_density_labels_basic)
    all_tests(24) = test_case("test_density_labels_zero_dimensions", test_density_labels_zero_dimensions)
    all_tests(25) = test_case("test_seeds_two_separated_clusters", test_seeds_two_separated_clusters)
    all_tests(26) = test_case("test_seeds_invalid_percentile", test_seeds_invalid_percentile)
    all_tests(27) = test_case("test_grow_ensemble_single_member", test_grow_ensemble_single_member)
    all_tests(28) = test_case("test_grow_ensemble_multi_member_union", test_grow_ensemble_multi_member_union)
    all_tests(29) = test_case("test_grow_ensemble_empty_ensemble", test_grow_ensemble_empty_ensemble)
    all_tests(30) = test_case("test_grow_ensemble_negative_radius", test_grow_ensemble_negative_radius)
    all_tests(31) = test_case("test_grow_ensemble_zero_dimensions", test_grow_ensemble_zero_dimensions)
    all_tests(32) = test_case("test_observable_full_rank_rectangle", test_observable_full_rank_rectangle)
    all_tests(33) = test_case("test_observable_low_rank_padding", test_observable_low_rank_padding)
    all_tests(34) = test_case("test_observable_too_few_members", test_observable_too_few_members)
    all_tests(35) = test_case("test_observable_dimension_too_small", test_observable_dimension_too_small)
    all_tests(36) = test_case("test_accept_ensemble_identical", test_accept_ensemble_identical)
    all_tests(37) = test_case("test_accept_ensemble_angle_exceeds_max", test_accept_ensemble_angle_exceeds_max)
    all_tests(38) = test_case("test_accept_ensemble_angle_within_max", test_accept_ensemble_angle_within_max)
    all_tests(39) = test_case("test_accept_ensemble_d_mismatch_within_dmax", test_accept_ensemble_d_mismatch_within_dmax)
    all_tests(40) = test_case("test_accept_ensemble_d_mismatch_exceeds_dmax", test_accept_ensemble_d_mismatch_exceeds_dmax)
    all_tests(41) = test_case("test_accept_ensemble_g_ratio_exceeds_max", test_accept_ensemble_g_ratio_exceeds_max)
    all_tests(42) = test_case("test_accept_ensemble_nonpositive_g", test_accept_ensemble_nonpositive_g)
    all_tests(43) = test_case("test_accept_ensemble_d_out_of_range", test_accept_ensemble_d_out_of_range)
  end function get_all_tests

  !> Run all shape_truthful_clustering tests.
  subroutine run_all_tests_shape_truthful_clustering()
    type(test_case) :: all_tests(43)
    integer(int32) :: i

    all_tests = get_all_tests()

    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All shape_truthful_clustering tests passed successfully."
  end subroutine run_all_tests_shape_truthful_clustering

  !> Run specific shape_truthful_clustering tests by name.
  subroutine run_named_tests_shape_truthful_clustering(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(43)
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
  end subroutine run_named_tests_shape_truthful_clustering

  ! --- normal_error -----------------------------------------------------

  !> Basic case: D=3, d=1 -- normal_error is the sum of the two smallest
  !| (normal-space) eigenvalues.
  subroutine test_normal_error_basic()
    integer(int32), parameter :: d_dim = 3, d = 1
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'normal_error failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(normal_error_value, 5.0d0, 1.0d-12, "normal_error basic case (d=1)")
  end subroutine test_normal_error_basic

  !> d=0: no tangent directions, everything is "normal" -- the full sum.
  subroutine test_normal_error_zero_tangent_dims()
    integer(int32), parameter :: d_dim = 3, d = 0
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'normal_error failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(normal_error_value, 14.0d0, 1.0d-12, "normal_error with d=0 (full sum)")
  end subroutine test_normal_error_zero_tangent_dims

  !> d=D: every direction is tangent, nothing left over -- sum over the
  !| empty range must be exactly zero.
  subroutine test_normal_error_all_tangent_dims()
    integer(int32), parameter :: d_dim = 3, d = 3
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'normal_error failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(normal_error_value, 0.0d0, 1.0d-12, "normal_error with d=D (empty sum)")
  end subroutine test_normal_error_all_tangent_dims

  !> d > D must be rejected by validation, not silently read out of bounds.
  subroutine test_normal_error_d_out_of_range()
    integer(int32), parameter :: d_dim = 3, d = 4
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    call assert_true(is_err(ierr), "normal_error should reject d > n_dimensions")
  end subroutine test_normal_error_d_out_of_range

  !> A negative eigenvalue is not a valid covariance eigenvalue and must be
  !| rejected by validation.
  subroutine test_normal_error_negative_eigenvalue()
    integer(int32), parameter :: d_dim = 3, d = 1
    real(real64) :: eigenvalues(d_dim) = [9.0d0, -1.0d0, 1.0d0]
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    call assert_true(is_err(ierr), "normal_error should reject a negative eigenvalue")
  end subroutine test_normal_error_negative_eigenvalue

  !> n_dimensions=0 must be rejected by validation.
  subroutine test_normal_error_zero_dimensions()
    integer(int32), parameter :: d_dim = 0, d = 0
    real(real64) :: eigenvalues(d_dim)
    real(real64) :: normal_error_value
    integer(int32) :: ierr

    call normal_error(d, eigenvalues, d_dim, normal_error_value, ierr)
    call assert_true(is_err(ierr), "normal_error should reject n_dimensions=0")
  end subroutine test_normal_error_zero_dimensions

  ! --- tangent_scales ----------------------------------------------------

  !> Basic case: D=3, d=2 -- tangent_scales is the square root of the two
  !| largest (tangent-space) eigenvalues.
  subroutine test_tangent_scales_basic()
    integer(int32), parameter :: d_dim = 3, d = 2
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    real(real64) :: expected(d) = [3.0d0, 2.0d0]
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'tangent_scales failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_array_real(tangent_scales_value, expected, d, 1.0d-12, "tangent_scales basic case (d=2)")
  end subroutine test_tangent_scales_basic

  !> d=D: every direction is tangent.
  subroutine test_tangent_scales_all_dims()
    integer(int32), parameter :: d_dim = 3, d = 3
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    real(real64) :: expected(d) = [3.0d0, 2.0d0, 1.0d0]
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'tangent_scales failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_array_real(tangent_scales_value, expected, d, 1.0d-12, "tangent_scales with d=D")
  end subroutine test_tangent_scales_all_dims

  !> d=0: no tangent directions -- must return a well-defined, empty (size
  !| zero) array rather than fail.
  subroutine test_tangent_scales_zero_dims()
    integer(int32), parameter :: d_dim = 3, d = 0
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'tangent_scales failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_int(size(tangent_scales_value), 0, "tangent_scales with d=0 should be an empty array")
  end subroutine test_tangent_scales_zero_dims

  !> d > D must be rejected by validation, not silently read out of bounds.
  subroutine test_tangent_scales_d_out_of_range()
    integer(int32), parameter :: d_dim = 3, d = 4
    real(real64) :: eigenvalues(d_dim) = [9.0d0, 4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "tangent_scales should reject d > n_dimensions")
  end subroutine test_tangent_scales_d_out_of_range

  !> A negative eigenvalue is not a valid covariance eigenvalue and must be
  !| rejected by validation.
  subroutine test_tangent_scales_negative_eigenvalue()
    integer(int32), parameter :: d_dim = 3, d = 2
    real(real64) :: eigenvalues(d_dim) = [9.0d0, -4.0d0, 1.0d0]
    real(real64) :: tangent_scales_value(d)
    integer(int32) :: ierr

    call tangent_scales(d, eigenvalues, d_dim, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "tangent_scales should reject a negative eigenvalue")
  end subroutine test_tangent_scales_negative_eigenvalue

  ! --- calculate_density_radius ------------------------------------------
  !
  ! Shared fixture for these tests: D=2, N=5 points on a line,
  ! (0,0),(1,0),(2,0),(3,0),(4,0). Mean = (2,0), so the mean-to-vector
  ! distances are [2,1,0,1,2], sorted ascending [0,1,1,2,2]. Expected
  ! percentile values below are hand-computed from calc_percentile_rank's
  ! documented linear-interpolation formula:
  !   rank = (percentile/100)*(n-1) + 1
  ! For the default 15th percentile: rank = 0.15*4+1 = 1.6 -> interpolate
  ! between sorted[1]=0 and sorted[2]=1 at fraction 0.6 -> 0.6.
  ! For the 50th percentile: rank = 0.5*4+1 = 3.0 exactly -> sorted[3] = 1.0.

  !> Default percentile (15%), via the allocating wrapper.
  subroutine test_density_radius_alloc_default_percentile()
    integer(int32), parameter :: d_dim = 2, n = 5
    real(real64) :: vectors(d_dim, n) = reshape( &
      [0.0d0,0.0d0, 1.0d0,0.0d0, 2.0d0,0.0d0, 3.0d0,0.0d0, 4.0d0,0.0d0], [d_dim, n])
    real(real64) :: radius
    integer(int32) :: ierr

    call calculate_density_radius_alloc(vectors, d_dim, n, radius, ierr=ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'calculate_density_radius_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(radius, 0.6d0, 1.0d-12, "density radius, default (15th) percentile")
  end subroutine test_density_radius_alloc_default_percentile

  !> Custom percentile (50%, the median), via the allocating wrapper.
  subroutine test_density_radius_alloc_custom_percentile()
    integer(int32), parameter :: d_dim = 2, n = 5
    real(real64) :: vectors(d_dim, n) = reshape( &
      [0.0d0,0.0d0, 1.0d0,0.0d0, 2.0d0,0.0d0, 3.0d0,0.0d0, 4.0d0,0.0d0], [d_dim, n])
    real(real64) :: radius
    integer(int32) :: ierr

    call calculate_density_radius_alloc(vectors, d_dim, n, radius, 0.5d0, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'calculate_density_radius_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(radius, 1.0d0, 1.0d-12, "density radius, custom (50th) percentile")
  end subroutine test_density_radius_alloc_custom_percentile

  !> Same 50th-percentile case, via the workspace-based validated entry
  !| point (the calling convention future SKGs like `seeds` will actually
  !| use, to avoid repeated allocation per call).
  subroutine test_density_radius_preallocated()
    integer(int32), parameter :: d_dim = 2, n = 5
    real(real64) :: vectors(d_dim, n) = reshape( &
      [0.0d0,0.0d0, 1.0d0,0.0d0, 2.0d0,0.0d0, 3.0d0,0.0d0, 4.0d0,0.0d0], [d_dim, n])
    real(real64) :: mean_vec(d_dim), distances(n), radius
    integer(int32) :: perm(n), ierr

    call calculate_density_radius(vectors, d_dim, n, mean_vec, distances, perm, radius, 0.5d0, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'calculate_density_radius failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(radius, 1.0d0, 1.0d-12, "density radius, preallocated workspace entry")
    call assert_equal_array_real(mean_vec, [2.0d0, 0.0d0], d_dim, 1.0d-12, "density radius, mean vector")
  end subroutine test_density_radius_preallocated

  !> A percentile outside [0,1] must be rejected by validation.
  subroutine test_density_radius_invalid_percentile()
    integer(int32), parameter :: d_dim = 2, n = 5
    real(real64) :: vectors(d_dim, n) = reshape( &
      [0.0d0,0.0d0, 1.0d0,0.0d0, 2.0d0,0.0d0, 3.0d0,0.0d0, 4.0d0,0.0d0], [d_dim, n])
    real(real64) :: radius
    integer(int32) :: ierr

    call calculate_density_radius_alloc(vectors, d_dim, n, radius, 1.5d0, ierr)
    call assert_true(is_err(ierr), "density radius should reject a percentile outside [0,1]")
  end subroutine test_density_radius_invalid_percentile

  !> n_vectors=0 must be rejected by validation.
  subroutine test_density_radius_zero_vectors()
    integer(int32), parameter :: d_dim = 2, n = 0
    real(real64) :: vectors(d_dim, n)
    real(real64) :: radius
    integer(int32) :: ierr

    call calculate_density_radius_alloc(vectors, d_dim, n, radius, ierr=ierr)
    call assert_true(is_err(ierr), "density radius should reject n_vectors=0")
  end subroutine test_density_radius_zero_vectors

  !> n_dimensions=0 must be rejected by validation.
  subroutine test_density_radius_zero_dimensions()
    integer(int32), parameter :: d_dim = 0, n = 5
    real(real64) :: vectors(d_dim, n)
    real(real64) :: radius
    integer(int32) :: ierr

    call calculate_density_radius_alloc(vectors, d_dim, n, radius, ierr=ierr)
    call assert_true(is_err(ierr), "density radius should reject n_dimensions=0")
  end subroutine test_density_radius_zero_dimensions

  ! --- calc_ensemble_growth_radius ---------------------------------------
  !
  ! Shared fixture: 11 points on a line, x = 0..10 (y=0), 1-indexed so point
  ! i has x = i-1. Seed = point 6 (x=5). Its nearest neighbors excluding
  ! itself, by absolute distance, are x=4 (dist 1), x=6 (dist 1), then a
  ! tie at distance 2 between x=3 and x=7. Which of the tied pair a k-NN
  ! query happens to return doesn't matter here: both have the same
  ! distance, so the resulting *distance* list is the same either way.
  !   k_min=4 (even): distances [1,1,2,2] -> median = 0.5*(1+2) = 1.5
  !   k_min=3 (odd):  distances [1,1,2]   -> median = 1.0

  !> Helper: build a k-d tree over the shared 11-point line fixture.
  subroutine build_line_fixture(vectors, kd_indices, dim_order)
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64),   intent(out) :: vectors(d_dim, n)
    integer(int32), intent(out) :: kd_indices(n)
    integer(int32), intent(out) :: dim_order(d_dim)
    integer(int32) :: work(n), perm(n), stack_left(n), stack_right(n), recursion_stack(3, n)
    integer(int32) :: i, ierr
    real(real64)   :: subarray(n)

    do i = 1, n
      vectors(1, i) = real(i - 1, real64)
      vectors(2, i) = 0.0d0
    end do
    dim_order = [1, 2]

    call build_kd_index(vectors, d_dim, n, kd_indices, dim_order, work, subarray, perm, &
                        stack_left, stack_right, recursion_stack, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'build_line_fixture: build_kd_index failed: ', ierr
      error stop
    end if
  end subroutine build_line_fixture

  !> Even k_min: median averages the two middle sorted distances.
  subroutine test_growth_radius_even_k()
    integer(int32), parameter :: d_dim = 2, n = 11, seed_index = 6
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: growth_radius
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    call calc_ensemble_growth_radius_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                                           seed_index, growth_radius, 4, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'calc_ensemble_growth_radius_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(growth_radius, 1.5d0, 1.0d-12, "growth radius, k_min=4 (even)")
  end subroutine test_growth_radius_even_k

  !> Odd k_min: median is the single middle sorted distance.
  subroutine test_growth_radius_odd_k()
    integer(int32), parameter :: d_dim = 2, n = 11, seed_index = 6
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: growth_radius
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    call calc_ensemble_growth_radius_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                                           seed_index, growth_radius, 3, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'calc_ensemble_growth_radius_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(growth_radius, 1.0d0, 1.0d-12, "growth radius, k_min=3 (odd)")
  end subroutine test_growth_radius_odd_k

  !> seed_index outside [1, n_vectors] must be rejected by validation.
  subroutine test_growth_radius_seed_index_out_of_range()
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: growth_radius
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    call calc_ensemble_growth_radius_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                                           n + 1, growth_radius, 4, ierr)
    call assert_true(is_err(ierr), "growth radius should reject seed_index > n_vectors")
  end subroutine test_growth_radius_seed_index_out_of_range

  !> k_min >= n_vectors leaves no room to exclude the seed itself and must
  !| be rejected by validation.
  subroutine test_growth_radius_k_min_too_large()
    integer(int32), parameter :: d_dim = 2, n = 11, seed_index = 6
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: growth_radius
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    call calc_ensemble_growth_radius_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                                           seed_index, growth_radius, n, ierr)
    call assert_true(is_err(ierr), "growth radius should reject k_min >= n_vectors")
  end subroutine test_growth_radius_k_min_too_large

  !> n_dimensions=0 must be rejected by validation.
  subroutine test_growth_radius_zero_dimensions()
    integer(int32), parameter :: d_dim = 0, n = 11, seed_index = 6
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: growth_radius
    integer(int32) :: ierr

    call calc_ensemble_growth_radius_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                                           seed_index, growth_radius, 4, ierr)
    call assert_true(is_err(ierr), "growth radius should reject n_dimensions=0")
  end subroutine test_growth_radius_zero_dimensions

  ! --- density_labels -----------------------------------------------------
  !
  ! Reuses the 11-point line fixture (x=0..10, y=0). With radius=1.5, an
  ! interior point's label counts itself plus its two immediate neighbors
  ! (3); an edge point counts itself plus its one immediate neighbor (2).

  !> Basic case: interior point vs. both edge points, radius=1.5.
  subroutine test_density_labels_basic()
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: labels(n)
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    call density_labels(vectors, d_dim, n, kd_indices, dim_order, 1.5d0, labels, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'density_labels failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_real(labels(1), 2.0d0, 1.0d-12, "density label at left edge (x=0)")
    call assert_equal_real(labels(6), 3.0d0, 1.0d-12, "density label at interior point (x=5)")
    call assert_equal_real(labels(11), 2.0d0, 1.0d-12, "density label at right edge (x=10)")
  end subroutine test_density_labels_basic

  !> n_dimensions=0 must be rejected by validation.
  subroutine test_density_labels_zero_dimensions()
    integer(int32), parameter :: d_dim = 0, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    real(real64)   :: labels(n)
    integer(int32) :: ierr

    call density_labels(vectors, d_dim, n, kd_indices, dim_order, 1.5d0, labels, ierr)
    call assert_true(is_err(ierr), "density_labels should reject n_dimensions=0")
  end subroutine test_density_labels_zero_dimensions

  ! --- seeds ---------------------------------------------------------------

  !> Two tight, well-separated clusters (spread ~0.2 apart internally,
  !| ~100 apart from each other) should yield exactly one seed per
  !| cluster: the default (15th percentile) density radius, computed from
  !| distances to the *global* mean, comes out far larger than the
  !| within-cluster spread but far smaller than the between-cluster gap,
  !| so the first pick in each cluster covers that whole cluster without
  !| reaching the other one. Which specific point within a cluster becomes
  !| the seed is not asserted (all three members of a cluster are
  !| density-tied at this radius, and heapsort's tie-breaking is not part
  !| of the SKG's contract) -- only that each cluster gets exactly one.
  subroutine test_seeds_two_separated_clusters()
    integer(int32), parameter :: d_dim = 2, n = 6
    real(real64) :: vectors(d_dim, n) = reshape( &
      [0.0d0,0.0d0, 0.1d0,0.0d0, 0.2d0,0.0d0, &
       100.0d0,0.0d0, 100.1d0,0.0d0, 100.2d0,0.0d0], [d_dim, n])
    logical :: is_seed_mask(n)
    integer(int32) :: ierr

    call seeds_alloc(vectors, d_dim, n, is_seed_mask, ierr=ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'seeds_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_int(count(is_seed_mask), 2, "seeds: expected exactly 2 seeds total")
    call assert_equal_int(count(is_seed_mask(1:3)), 1, "seeds: expected exactly 1 seed in cluster A")
    call assert_equal_int(count(is_seed_mask(4:6)), 1, "seeds: expected exactly 1 seed in cluster B")
  end subroutine test_seeds_two_separated_clusters

  !> A percentile outside [0,1] must be rejected by validation.
  subroutine test_seeds_invalid_percentile()
    integer(int32), parameter :: d_dim = 2, n = 6
    real(real64) :: vectors(d_dim, n) = reshape( &
      [0.0d0,0.0d0, 0.1d0,0.0d0, 0.2d0,0.0d0, &
       100.0d0,0.0d0, 100.1d0,0.0d0, 100.2d0,0.0d0], [d_dim, n])
    logical :: is_seed_mask(n)
    integer(int32) :: ierr

    call seeds_alloc(vectors, d_dim, n, is_seed_mask, 1.5d0, ierr)
    call assert_true(is_err(ierr), "seeds should reject a percentile outside [0,1]")
  end subroutine test_seeds_invalid_percentile

  ! --- grow_ensemble -------------------------------------------------------
  !
  ! Reuses the 11-point line fixture (x=0..10, y=0). With growth_radius=1.5:
  !   from a single member x=5 (index 6): covers x=4,5,6 (indices 5,6,7).
  !   from members x=4,5,6 (indices 5,6,7): the union of each member's own
  !   1.5-radius neighborhood covers x=3..7 (indices 4..8) -- e.g. x=4's own
  !   neighborhood alone already reaches x=3, which x=5's and x=6's
  !   neighborhoods (starting no earlier than x=3.5) do not.

  !> Growing from a single member.
  subroutine test_grow_ensemble_single_member()
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    logical        :: is_member_mask(n), is_member_mask_next(n), expected(n)
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    is_member_mask = .false.
    is_member_mask(6) = .true.

    call grow_ensemble_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                             is_member_mask, 1.5d0, is_member_mask_next, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'grow_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if

    expected = .false.
    expected(5:7) = .true.
    call assert_equal_array_logical(is_member_mask_next, expected, n, "grow_ensemble from a single member")
  end subroutine test_grow_ensemble_single_member

  !> Growing from several members: the result is the union of each
  !| member's own neighborhood, not just the neighborhood of one of them.
  subroutine test_grow_ensemble_multi_member_union()
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    logical        :: is_member_mask(n), is_member_mask_next(n), expected(n)
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)

    is_member_mask = .false.
    is_member_mask(5:7) = .true.

    call grow_ensemble_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                             is_member_mask, 1.5d0, is_member_mask_next, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'grow_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if

    expected = .false.
    expected(4:8) = .true.
    call assert_equal_array_logical(is_member_mask_next, expected, n, "grow_ensemble from three members (union)")
  end subroutine test_grow_ensemble_multi_member_union

  !> An empty ensemble (no members to grow from) must be rejected by
  !| validation.
  subroutine test_grow_ensemble_empty_ensemble()
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    logical        :: is_member_mask(n), is_member_mask_next(n)
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)
    is_member_mask = .false.

    call grow_ensemble_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                             is_member_mask, 1.5d0, is_member_mask_next, ierr)
    call assert_true(is_err(ierr), "grow_ensemble should reject an empty ensemble")
  end subroutine test_grow_ensemble_empty_ensemble

  !> A negative growth radius must be rejected by validation.
  subroutine test_grow_ensemble_negative_radius()
    integer(int32), parameter :: d_dim = 2, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    logical        :: is_member_mask(n), is_member_mask_next(n)
    integer(int32) :: ierr

    call build_line_fixture(vectors, kd_indices, dim_order)
    is_member_mask = .false.
    is_member_mask(6) = .true.

    call grow_ensemble_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                             is_member_mask, -1.5d0, is_member_mask_next, ierr)
    call assert_true(is_err(ierr), "grow_ensemble should reject a negative growth radius")
  end subroutine test_grow_ensemble_negative_radius

  !> n_dimensions=0 must be rejected by validation.
  subroutine test_grow_ensemble_zero_dimensions()
    integer(int32), parameter :: d_dim = 0, n = 11
    real(real64)   :: vectors(d_dim, n)
    integer(int32) :: kd_indices(n), dim_order(d_dim)
    logical        :: is_member_mask(n), is_member_mask_next(n)
    integer(int32) :: ierr

    is_member_mask = .false.
    is_member_mask(6) = .true.

    call grow_ensemble_alloc(vectors, d_dim, n, kd_indices, dim_order, &
                             is_member_mask, 1.5d0, is_member_mask_next, ierr)
    call assert_true(is_err(ierr), "grow_ensemble should reject n_dimensions=0")
  end subroutine test_grow_ensemble_zero_dimensions

  ! --- observable ------------------------------------------------------
  !
  ! observable does not need a k-d tree (only vectors + is_member_mask), so
  ! these tests build their own small fixtures directly rather than reusing
  ! build_line_fixture.

  !> A rectangle in the z=0 plane, embedded in 3D: full economy-mode rank
  !| (rank = min(D,k) = 3 = D, no zero-padding), with distinct, hand-
  !| computable eigenvalues -- the x/y cross term is exactly zero for a
  !| centered rectangle, so U is deterministically the standard basis
  !| (up to a per-column sign flip, which SVD never fixes).
  subroutine test_observable_full_rank_rectangle()
    integer(int32), parameter :: d_dim = 3, n = 4
    real(real64)   :: vectors(d_dim, n)
    logical        :: is_member_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr
    real(real64), parameter :: tol = 1.0d-6

    vectors(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    vectors(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    vectors(:, 3) = [0.0d0, 1.0d0, 0.0d0]
    vectors(:, 4) = [2.0d0, 1.0d0, 0.0d0]
    is_member_mask = .true.

    call observable_alloc(vectors, d_dim, n, is_member_mask, &
                          U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'observable_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_array_real(mu, [1.0d0, 0.5d0, 0.0d0], d_dim, tol, "observable: rectangle mean")
    call assert_equal_int(d, 2, "observable: rectangle intrinsic dimension")
    call assert_equal_real(eigenvalues(1), 4.0d0/3.0d0, tol, "observable: rectangle lambda_1")
    call assert_equal_real(eigenvalues(2), 1.0d0/3.0d0, tol, "observable: rectangle lambda_2")
    call assert_equal_real(eigenvalues(3), 0.0d0, tol, "observable: rectangle lambda_3")
    call assert_equal_real(normal_error_value, 0.0d0, tol, "observable: rectangle normal_error")
    call assert_equal_real(tangent_scales_value(1), sqrt(4.0d0/3.0d0), tol, "observable: rectangle tangent_scales(1)")
    call assert_equal_real(tangent_scales_value(2), sqrt(1.0d0/3.0d0), tol, "observable: rectangle tangent_scales(2)")
    call assert_equal_real(tangent_scales_value(3), 0.0d0, tol, "observable: rectangle tangent_scales(3)")
    call assert_true(G > 1.0d10, "observable: rectangle spectral gap should be huge at the true rank boundary")

    ! U columns are the standard basis up to sign (diagonal covariance with
    ! distinct eigenvalues 4, 1, 0 leaves no rotational freedom).
    call assert_equal_real(abs(U(1,1)), 1.0d0, tol, "observable: rectangle U column 1, x")
    call assert_equal_real(abs(U(2,1)), 0.0d0, tol, "observable: rectangle U column 1, y")
    call assert_equal_real(abs(U(3,1)), 0.0d0, tol, "observable: rectangle U column 1, z")
    call assert_equal_real(abs(U(1,2)), 0.0d0, tol, "observable: rectangle U column 2, x")
    call assert_equal_real(abs(U(2,2)), 1.0d0, tol, "observable: rectangle U column 2, y")
    call assert_equal_real(abs(U(3,2)), 0.0d0, tol, "observable: rectangle U column 2, z")
    call assert_equal_real(abs(U(1,3)), 0.0d0, tol, "observable: rectangle U column 3, x")
    call assert_equal_real(abs(U(2,3)), 0.0d0, tol, "observable: rectangle U column 3, y")
    call assert_equal_real(abs(U(3,3)), 1.0d0, tol, "observable: rectangle U column 3, z")
  end subroutine test_observable_full_rank_rectangle

  !> Three collinear points (intrinsic rank 1) embedded in a 5D ambient
  !| space: economy-mode rank = min(D,k) = 3 < D = 5, so U columns 4-5 and
  !| eigenvalues 4-5 must be observable's own zero-padding, not LAPACK
  !| output.
  subroutine test_observable_low_rank_padding()
    integer(int32), parameter :: d_dim = 5, n = 3
    real(real64)   :: vectors(d_dim, n)
    logical        :: is_member_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr
    real(real64), parameter :: tol = 1.0d-8

    vectors      = 0.0d0
    vectors(1,1) = 0.0d0
    vectors(1,2) = 1.0d0
    vectors(1,3) = 2.0d0
    is_member_mask = .true.

    call observable_alloc(vectors, d_dim, n, is_member_mask, &
                          U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'observable_alloc failed unexpectedly: ', ierr
      error stop
    end if

    call assert_equal_int(d, 1, "observable: collinear-in-5D intrinsic dimension")
    call assert_equal_real(eigenvalues(1), 1.0d0, tol, "observable: collinear-in-5D lambda_1")
    call assert_equal_array_real(eigenvalues(4:5), [0.0d0, 0.0d0], 2, tol, &
                                 "observable: collinear-in-5D padded eigenvalues (4:5) must be exactly zero")
    call assert_equal_array_real(U(:,4), [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0], d_dim, tol, &
                                 "observable: collinear-in-5D padded U column 4 must be exactly zero")
    call assert_equal_array_real(U(:,5), [0.0d0,0.0d0,0.0d0,0.0d0,0.0d0], d_dim, tol, &
                                 "observable: collinear-in-5D padded U column 5 must be exactly zero")
    call assert_equal_real(tangent_scales_value(1), 1.0d0, tol, "observable: collinear-in-5D tangent_scales(1)")
    call assert_equal_array_real(tangent_scales_value(2:5), [0.0d0,0.0d0,0.0d0,0.0d0], 4, tol, &
                                 "observable: collinear-in-5D padded tangent_scales(2:5) must be exactly zero")
    call assert_true(G > 1.0d10, "observable: collinear-in-5D spectral gap should be huge")
  end subroutine test_observable_low_rank_padding

  !> A single-member ensemble must be rejected: no meaningful covariance/SVD.
  subroutine test_observable_too_few_members()
    integer(int32), parameter :: d_dim = 3, n = 4
    real(real64)   :: vectors(d_dim, n)
    logical        :: is_member_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr

    vectors(:, 1) = [0.0d0, 0.0d0, 0.0d0]
    vectors(:, 2) = [2.0d0, 0.0d0, 0.0d0]
    vectors(:, 3) = [0.0d0, 1.0d0, 0.0d0]
    vectors(:, 4) = [2.0d0, 1.0d0, 0.0d0]
    is_member_mask = .false.
    is_member_mask(1) = .true.

    call observable_alloc(vectors, d_dim, n, is_member_mask, &
                          U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "observable should reject an ensemble with fewer than 2 members")
  end subroutine test_observable_too_few_members

  !> n_dimensions=1 leaves no room for even a single spectral-gap
  !| comparison (r ranges 1..D-1) and must be rejected.
  subroutine test_observable_dimension_too_small()
    integer(int32), parameter :: d_dim = 1, n = 3
    real(real64)   :: vectors(d_dim, n)
    logical        :: is_member_mask(n)
    real(real64)   :: U(d_dim, d_dim), eigenvalues(d_dim), mu(d_dim)
    real(real64)   :: normal_error_value, tangent_scales_value(d_dim), G
    integer(int32) :: d, ierr

    vectors(:, 1) = [0.0d0]
    vectors(:, 2) = [1.0d0]
    vectors(:, 3) = [2.0d0]
    is_member_mask = .true.

    call observable_alloc(vectors, d_dim, n, is_member_mask, &
                          U, eigenvalues, mu, d, G, normal_error_value, tangent_scales_value, ierr)
    call assert_true(is_err(ierr), "observable should reject n_dimensions=1")
  end subroutine test_observable_dimension_too_small

  ! --- accept_ensemble ---------------------------------------------------
  !
  ! accept_ensemble does not need a k-d tree or member vectors -- only the
  ! two observable tuples (U, d, G) at t and t+1 -- so these tests build
  ! small hand-picked bases/values directly.

  !> Identical basis, identical d, identical G: all three criteria trivially
  !| satisfied, even at zero tolerance for d and G.
  subroutine test_accept_ensemble_identical()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 5.0d0, U_tp1, 2_int32, 5.0d0, &
                               0.1d0, 0_int32, 0.0d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: identical observables must be accepted")
  end subroutine test_accept_ensemble_identical

  !> A 60-degree rotation of a 1D tangent basis, rejected by an alpha_max of 30 degrees.
  subroutine test_accept_ensemble_angle_exceeds_max()
    integer(int32), parameter :: d_dim = 2
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim), pi, theta, alpha_max
    logical        :: is_accepted
    integer(int32) :: ierr

    pi    = 4.0d0 * atan(1.0d0)
    theta = pi / 3.0d0   ! 60 degrees
    alpha_max = pi / 6.0d0  ! 30 degrees

    U_t = 0.0d0
    U_t(1,1) = 1.0d0
    U_t(2,2) = 1.0d0

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(theta)
    U_tp1(2,1) = sin(theta)
    U_tp1(2,2) = 1.0d0

    call accept_ensemble_alloc(d_dim, U_t, 1_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                               alpha_max, 0_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a 60-degree tangent rotation must be rejected at alpha_max=30deg")
  end subroutine test_accept_ensemble_angle_exceeds_max

  !> The same 60-degree rotation, accepted once alpha_max is raised to 70 degrees.
  subroutine test_accept_ensemble_angle_within_max()
    integer(int32), parameter :: d_dim = 2
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim), pi, theta, alpha_max
    logical        :: is_accepted
    integer(int32) :: ierr

    pi    = 4.0d0 * atan(1.0d0)
    theta = pi / 3.0d0        ! 60 degrees
    alpha_max = 7.0d0 * pi / 18.0d0  ! 70 degrees

    U_t = 0.0d0
    U_t(1,1) = 1.0d0
    U_t(2,2) = 1.0d0

    U_tp1 = 0.0d0
    U_tp1(1,1) = cos(theta)
    U_tp1(2,1) = sin(theta)
    U_tp1(2,2) = 1.0d0

    call accept_ensemble_alloc(d_dim, U_t, 1_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                               alpha_max, 0_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: a 60-degree tangent rotation must be accepted at alpha_max=70deg")
  end subroutine test_accept_ensemble_angle_within_max

  !> d changed by 1 (2 -> 1): the angle criterion is vacuously satisfied
  !| (no common dimension to compare), and d_max=1 tolerates the change.
  subroutine test_accept_ensemble_d_mismatch_within_dmax()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                               0.1d0, 1_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(is_accepted, "accept_ensemble: a change in d of 1 must be accepted when d_max=1")
  end subroutine test_accept_ensemble_d_mismatch_within_dmax

  !> The same change in d (2 -> 1), rejected once d_max=0.
  subroutine test_accept_ensemble_d_mismatch_exceeds_dmax()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 1.0d0, U_tp1, 1_int32, 1.0d0, &
                               0.1d0, 0_int32, 1.0d10, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a change in d of 1 must be rejected when d_max=0")
  end subroutine test_accept_ensemble_d_mismatch_exceeds_dmax

  !> G changes by a factor of 10 (ln(10) ~ 2.303): rejected at G_max=1.0.
  subroutine test_accept_ensemble_g_ratio_exceeds_max()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 1.0d0, U_tp1, 2_int32, 10.0d0, &
                               0.1d0, 0_int32, 1.0d0, is_accepted, ierr)
    if (.not. is_ok(ierr)) then
      write(*,*) 'accept_ensemble_alloc failed unexpectedly: ', ierr
      error stop
    end if
    call assert_true(.not. is_accepted, "accept_ensemble: a 10x change in G must be rejected at G_max=1.0")
  end subroutine test_accept_ensemble_g_ratio_exceeds_max

  !> Non-positive spectral gaps make log(G_tp1/G_t) undefined and must be rejected by validation.
  subroutine test_accept_ensemble_nonpositive_g()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, 2_int32, 0.0d0, U_tp1, 2_int32, 1.0d0, &
                               0.1d0, 0_int32, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject G_t <= 0")
  end subroutine test_accept_ensemble_nonpositive_g

  !> d_t out of [0, n_dimensions] must be rejected by validation.
  subroutine test_accept_ensemble_d_out_of_range()
    integer(int32), parameter :: d_dim = 3
    real(real64)   :: U_t(d_dim, d_dim), U_tp1(d_dim, d_dim)
    logical        :: is_accepted
    integer(int32) :: ierr, i

    U_t = 0.0d0
    do i = 1, d_dim
      U_t(i, i) = 1.0d0
    end do
    U_tp1 = U_t

    call accept_ensemble_alloc(d_dim, U_t, d_dim + 1, 1.0d0, U_tp1, 2_int32, 1.0d0, &
                               0.1d0, 0_int32, 1.0d10, is_accepted, ierr)
    call assert_true(is_err(ierr), "accept_ensemble should reject d_t > n_dimensions")
  end subroutine test_accept_ensemble_d_out_of_range

end module mod_test_shape_truthful_clustering
