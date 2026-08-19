!> Unit test suite for kd_tree module.
module mod_test_kd_tree
    use f42_kd_tree
    use f42_utils
    use tox_errors
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_kd_tree() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(27))

        all_tests(1) = test_case("test_kd_2d_cartesian", test_kd_2d_cartesian)
        all_tests(2) = test_case("test_kd_3d_spherical", test_kd_3d_spherical)
        all_tests(3) = test_case("test_kd_empty_array", test_kd_empty_array)
        all_tests(4) = test_case("test_kd_single_point", test_kd_single_point)
        all_tests(5) = test_case("test_kd_identical_points", test_kd_identical_points)
        all_tests(6) = test_case("test_kd_unit_vectors", test_kd_unit_vectors)
        all_tests(7) = test_case("test_kd_high_dim_low_points", test_kd_high_dim_low_points)
        all_tests(8) = test_case("test_kd_1d_sorted", test_kd_1d_sorted)
        all_tests(9) = test_case("test_kd_2d_minimal", test_kd_2d_minimal)
        all_tests(10) = test_case("test_kd_1d_minimal", test_kd_1d_minimal)
        all_tests(11) = test_case("test_kd_3d_large", test_kd_3d_large)
        all_tests(12) = test_case("test_kd_5d_medium", test_kd_5d_medium)
        all_tests(13) = test_case("test_vicinity_vectors_basic_2d", test_vicinity_vectors_basic_2d)
        all_tests(14) = test_case("test_vicinity_vectors_helper_matches_alloc", test_vicinity_vectors_helper_matches_alloc)
        all_tests(15) = test_case("test_vicinity_vectors_zero_radius_duplicates", test_vicinity_vectors_zero_radius_duplicates)
        all_tests(16) = test_case("test_vicinity_vectors_boundary_inclusive", test_vicinity_vectors_boundary_inclusive)
        all_tests(17) = test_case("test_vicinity_vectors_no_matches", test_vicinity_vectors_no_matches)
        all_tests(18) = test_case("test_vicinity_vectors_all_matches", test_vicinity_vectors_all_matches)
        all_tests(19) = test_case("test_vicinity_vectors_single_point", test_vicinity_vectors_single_point)
        all_tests(20) = test_case("test_vicinity_vectors_large_1d_bruteforce", test_vicinity_vectors_large_1d_bruteforce)
        all_tests(21) = test_case("test_vicinity_vectors_negative_radius", test_vicinity_vectors_negative_radius)
        all_tests(22) = test_case("test_vicinity_vectors_dimension_order_low", test_vicinity_vectors_dimension_order_low)
        all_tests(23) = test_case("test_vicinity_vectors_dimension_order_high", test_vicinity_vectors_dimension_order_high)
        all_tests(24) = test_case("test_vicinity_vectors_kd_index_low", test_vicinity_vectors_kd_index_low)
        all_tests(25) = test_case("test_vicinity_vectors_kd_index_high", test_vicinity_vectors_kd_index_high)
        all_tests(26) = test_case("test_vicinity_vectors_zero_points", test_vicinity_vectors_zero_points)
        all_tests(27) = test_case("test_vicinity_vectors_zero_dimensions", test_vicinity_vectors_zero_dimensions)
    end function get_all_tests_kd_tree

    !> Test 2D Cartesian KD-Tree.
    subroutine test_kd_2d_cartesian()
        integer(int32), parameter :: d = 2, n = 6
        real(real64) :: X(d, n) = reshape([1.0d0, 2.0d0, 2.0d0, 3.0d0, 3.0d0, 1.0d0, &
                                           4.0d0, 0.0d0, 0.0d0, 4.0d0, 5.0d0, 2.0d0], [d, n])
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(kd_ix, n, "2D Cartesian KD-Tree")
    end subroutine test_kd_2d_cartesian

    !> Test 3D Spherical KD-Tree.
    subroutine test_kd_3d_spherical()
        integer(int32), parameter :: d = 3, n = 8
        real(real64) :: V(d, n)
        integer(int32) :: sphere_ix(n), dim_order(d) = [1, 2, 3]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call random_unit_vectors(V, d, n)
        call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')
        call assert_permutation(sphere_ix, n, "3D Spherical KD-Tree")
    end subroutine test_kd_3d_spherical

    !> Test KD-Tree with empty array.
    subroutine test_kd_empty_array()
        integer(int32), parameter :: d = 2, n = 0
        real(real64) :: X(d, n)
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        ! This should return ERR_EMPTY_INPUT
        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, arg_pos=3_int32), 'Expected error, but did''nt get one')

        call assert_true(.true., "KD-Tree empty array handling")
    end subroutine test_kd_empty_array

    !> Test KD-Tree with single point.
    subroutine test_kd_single_point()
        integer(int32), parameter :: d = 2, n = 1
        real(real64) :: X(d, n) = reshape([1.0d0, 2.0d0], [d, n])
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_equal_int(kd_ix(1), 1, "KD-Tree single point index incorrect")
    end subroutine test_kd_single_point

    !> Test KD-Tree with identical points.
    subroutine test_kd_identical_points()
        integer(int32), parameter :: d = 3, n = 5
        real(real64) :: X(d, n) = 1.0d0
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2, 3]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(kd_ix, n, "KD-Tree identical points")
    end subroutine test_kd_identical_points

    !> Test KD-Tree with unit vectors.
    subroutine test_kd_unit_vectors()
        integer(int32), parameter :: d = 4, n = 4
        real(real64) :: V(d, n) = 0.0d0
        integer(int32) :: sphere_ix(n), dim_order(d) = [1, 2, 3, 4]
        integer(int32) :: work(n), perm(n)
        real(real64) :: subarray(n)
        integer(int32) :: i, ierr
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        do i = 1, d
            V(i, i) = 1.0d0
        end do

        call build_spherical_kd(V, d, n, sphere_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(sphere_ix, n, "KD-Tree unit vectors")
    end subroutine test_kd_unit_vectors

    !> Test KD-Tree with high dimension and few points.
    subroutine test_kd_high_dim_low_points()
        integer(int32), parameter :: d = 10, n = 3
        real(real64) :: X(d, n)
        integer :: i
        integer(int32) :: kd_ix(n), dim_order(d) = [(i, i=1, d)]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)
        call set_ok(ierr)

        call random_matrix(X, d, n)
        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')
        call assert_permutation(kd_ix, n, "KD-Tree high dimension low points")
    end subroutine test_kd_high_dim_low_points

    !> Test 1D sorted KD-Tree.
    subroutine test_kd_1d_sorted()
        integer(int32), parameter :: d = 1, n = 10
        real(real64) :: X(d, n)
        integer(int32) :: kd_ix(n), dim_order(d) = [1]
        integer(int32) :: work(n), perm(n)
        real(real64) :: subarray(n)
        integer(int32) :: i, ierr
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        do i = 1, n
            X(1, i) = i
        end do

        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(kd_ix, n, "1D sorted KD-Tree")
    end subroutine test_kd_1d_sorted

    !> Test 2D minimal KD-Tree.
    subroutine test_kd_2d_minimal()
        integer(int32), parameter :: d = 2, n = 2
        real(real64) :: X(d, n) = reshape([1.0d0, 2.0d0, 2.0d0, 1.0d0], [d, n])
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(kd_ix, n, "2D minimal KD-Tree")
    end subroutine test_kd_2d_minimal

    !> Test 1D minimal KD-Tree.
    subroutine test_kd_1d_minimal()
        integer(int32), parameter :: d = 1, n = 2
        real(real64) :: X(d, n) = reshape([1.0d0, 2.0d0], [d, n])
        integer(int32) :: kd_ix(n), dim_order(d) = [1]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(kd_ix, n, "1D minimal KD-Tree")
    end subroutine test_kd_1d_minimal

    !> Test 3D large KD-Tree.
    subroutine test_kd_3d_large()
        integer(int32), parameter :: d = 3, n = 100
        real(real64) :: X(d, n)
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2, 3]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        real(real64) :: val(d)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call random_matrix(X, d, n)
        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')
        call get_kd_point(X, kd_ix, 4, val, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')

        call assert_permutation(kd_ix, n, "3D large KD-Tree")
    end subroutine test_kd_3d_large

    !> Test 5D medium KD-Tree.
    subroutine test_kd_5d_medium()
        integer(int32), parameter :: d = 5, n = 10
        real(real64) :: X(d, n)
        integer(int32) :: kd_ix(n), dim_order(d) = [1, 2, 3, 4, 5]
        integer(int32) :: work(n), perm(n), ierr
        real(real64) :: subarray(n)
        integer(int32) :: recursion_stack(3, n)

        call set_ok(ierr)

        call random_matrix(X, d, n)
        call build_kd_index(X, d, n, kd_ix, dim_order, work, subarray, perm, recursion_stack, ierr)
        call assert_equal_int(ierr, ERR_OK, 'Unexpected error')
        call assert_permutation(kd_ix, n, "5D medium KD-Tree")
    end subroutine test_kd_5d_medium

    subroutine test_vicinity_vectors_basic_2d()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 6_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: radius, tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points), neighbor_count, ierr
        logical :: vicinity_mask(n_points), expected_mask(n_points)

        points(:, 1) = [0.0_real64, 0.0_real64]
        points(:, 2) = [1.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 1.0_real64]
        points(:, 4) = [2.0_real64, 0.0_real64]
        points(:, 5) = [0.0_real64, 2.0_real64]
        points(:, 6) = [3.0_real64, 3.0_real64]

        query_point = [0.0_real64, 0.0_real64]
        radius = 1.0_real64
        dimension_order = [1_int32, 2_int32]
        expected_mask = [.true., .true., .true., .false., .false., .false.]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_basic_2d: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, radius, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_basic_2d: vicinity_vectors_alloc")
        if (.not. is_ok(ierr)) return

        call assert_true(all(vicinity_mask .eqv. expected_mask), &
                         "test_vicinity_vectors_basic_2d: expected mask")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, radius, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_basic_2d: vicinity_vectors_count_alloc")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, 3_int32, &
                              "test_vicinity_vectors_basic_2d: expected count")

    end subroutine test_vicinity_vectors_basic_2d

    !> Test that helper and allocating routines produce identical results.
    subroutine test_vicinity_vectors_helper_matches_alloc()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 5_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: radius, tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points)
        integer(int32) :: tmp_mask_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        integer(int32) :: tmp_count_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        integer(int32) :: alloc_count, helper_count, ierr
        logical :: alloc_mask(n_points), helper_mask(n_points)

        points(:, 1) = [0.0_real64, 0.0_real64]
        points(:, 2) = [1.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 1.0_real64]
        points(:, 4) = [1.0_real64, 1.0_real64]
        points(:, 5) = [3.0_real64, 3.0_real64]

        query_point = [0.5_real64, 0.5_real64]
        radius = 0.75_real64
        dimension_order = [2_int32, 1_int32]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_helper_matches_alloc: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, radius, &
                                    dimension_order, kd_indices, alloc_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_helper_matches_alloc: mask alloc")
        if (.not. is_ok(ierr)) return

        tmp_mask_stack = -1_int32

        call vicinity_vectors_helper(query_point, points, n_dimensions, n_points, radius, &
                                     dimension_order, kd_indices, tmp_mask_stack, helper_mask)

        call assert_true(all(helper_mask .eqv. alloc_mask), &
                         "test_vicinity_vectors_helper_matches_alloc: masks agree")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, radius, &
                                          dimension_order, kd_indices, alloc_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_helper_matches_alloc: count alloc")
        if (.not. is_ok(ierr)) return

        tmp_count_stack = -1_int32

        call vicinity_vectors_count_helper(query_point, points, n_dimensions, n_points, radius, &
                                           dimension_order, kd_indices, tmp_count_stack, helper_count)

        call assert_equal_int(helper_count, alloc_count, &
                              "test_vicinity_vectors_helper_matches_alloc: counts agree")

        call assert_equal_int(helper_count, count(helper_mask, kind=int32), &
                              "test_vicinity_vectors_helper_matches_alloc: mask and count agree")

    end subroutine test_vicinity_vectors_helper_matches_alloc

    !> Test zero-radius search with duplicate points.
    subroutine test_vicinity_vectors_zero_radius_duplicates()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 5_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points), neighbor_count, ierr
        logical :: vicinity_mask(n_points), expected_mask(n_points)

        points(:, 1) = [0.0_real64, 0.0_real64]
        points(:, 2) = [0.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 0.0_real64]
        points(:, 4) = [1.0_real64, 0.0_real64]
        points(:, 5) = [0.0_real64, 1.0_real64]

        query_point = [0.0_real64, 0.0_real64]
        dimension_order = [1_int32, 2_int32]
        expected_mask = [.true., .true., .true., .false., .false.]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_zero_radius_duplicates: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 0.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_zero_radius_duplicates: mask call")
        if (.not. is_ok(ierr)) return

        call assert_true(all(vicinity_mask .eqv. expected_mask), &
                         "test_vicinity_vectors_zero_radius_duplicates: expected mask")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 0.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_zero_radius_duplicates: count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, 3_int32, &
                              "test_vicinity_vectors_zero_radius_duplicates: duplicate count")

    end subroutine test_vicinity_vectors_zero_radius_duplicates

    !> Test inclusion of points exactly on the radius boundary.
    subroutine test_vicinity_vectors_boundary_inclusive()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 5_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points), neighbor_count, ierr
        logical :: vicinity_mask(n_points), expected_mask(n_points)

        points(:, 1) = [1.0_real64, 0.0_real64]
        points(:, 2) = [-1.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 1.0_real64]
        points(:, 4) = [0.0_real64, -1.0_real64]
        points(:, 5) = [1.001_real64, 0.0_real64]

        query_point = [0.0_real64, 0.0_real64]
        dimension_order = [1_int32, 2_int32]
        expected_mask = [.true., .true., .true., .true., .false.]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_boundary_inclusive: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_boundary_inclusive: mask call")
        if (.not. is_ok(ierr)) return

        call assert_true(all(vicinity_mask .eqv. expected_mask), &
                         "test_vicinity_vectors_boundary_inclusive: expected mask")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_boundary_inclusive: count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, 4_int32, &
                              "test_vicinity_vectors_boundary_inclusive: expected count")

    end subroutine test_vicinity_vectors_boundary_inclusive

    !> Test a vicinity search with no matching points.
    subroutine test_vicinity_vectors_no_matches()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 4_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points), neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(:, 1) = [0.0_real64, 0.0_real64]
        points(:, 2) = [1.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 1.0_real64]
        points(:, 4) = [1.0_real64, 1.0_real64]

        query_point = [100.0_real64, 100.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_no_matches: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 0.5_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_no_matches: mask call")
        if (.not. is_ok(ierr)) return

        call assert_false(any(vicinity_mask), &
                          "test_vicinity_vectors_no_matches: mask is empty")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 0.5_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_no_matches: count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, 0_int32, &
                              "test_vicinity_vectors_no_matches: zero neighbors")

    end subroutine test_vicinity_vectors_no_matches

    !> Test a vicinity search containing all reference points.
    subroutine test_vicinity_vectors_all_matches()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 5_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points), neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(:, 1) = [-2.0_real64, -2.0_real64]
        points(:, 2) = [-1.0_real64, 1.0_real64]
        points(:, 3) = [0.0_real64, 0.0_real64]
        points(:, 4) = [1.0_real64, -1.0_real64]
        points(:, 5) = [2.0_real64, 2.0_real64]

        query_point = [0.0_real64, 0.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_all_matches: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 10.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_all_matches: mask call")
        if (.not. is_ok(ierr)) return

        call assert_true(all(vicinity_mask), &
                         "test_vicinity_vectors_all_matches: all points selected")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 10.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_all_matches: count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, n_points, &
                              "test_vicinity_vectors_all_matches: all points counted")

    end subroutine test_vicinity_vectors_all_matches

    !> Test vicinity searches using a single-point K-D tree.
    subroutine test_vicinity_vectors_single_point()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 1_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points), neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(:, 1) = [2.0_real64, -1.0_real64]
        dimension_order = [1_int32, 2_int32]

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_single_point: build K-D index")
        if (.not. is_ok(ierr)) return

        query_point = [2.0_real64, -1.0_real64]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 0.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_single_point: exact mask call")
        if (.not. is_ok(ierr)) return

        call assert_true(vicinity_mask(1), &
                         "test_vicinity_vectors_single_point: point selected")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 0.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_single_point: exact count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, 1_int32, &
                              "test_vicinity_vectors_single_point: one neighbor")

        query_point = [10.0_real64, 10.0_real64]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 0.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_single_point: distant mask call")
        if (.not. is_ok(ierr)) return

        call assert_false(vicinity_mask(1), &
                          "test_vicinity_vectors_single_point: distant point excluded")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 0.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_single_point: distant count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, 0_int32, &
                              "test_vicinity_vectors_single_point: distant count")

    end subroutine test_vicinity_vectors_single_point

    !> Test K-D vicinity search against direct brute-force evaluation.
    subroutine test_vicinity_vectors_large_1d_bruteforce()

        integer(int32), parameter :: n_dimensions = 1_int32, n_points = 201_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        real(real64) :: radius, tmp_value_buffer(n_points)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: tmp_workspace(n_points), tmp_permutation(n_points)
        integer(int32) :: tmp_recursion_stack(3, n_points)
        integer(int32) :: neighbor_count, expected_count, i_point, ierr
        logical :: vicinity_mask(n_points), expected_mask(n_points)

        do i_point = 1_int32, n_points
            points(1, i_point) = real(i_point - 1_int32, real64)
        end do

        query_point = [73.25_real64]
        radius = 10.25_real64
        dimension_order = [1_int32]

        do i_point = 1_int32, n_points
            expected_mask(i_point) = &
                abs(points(1, i_point) - query_point(1)) <= radius
        end do

        expected_count = count(expected_mask, kind=int32)

        call build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, &
                            tmp_recursion_stack, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_large_1d_bruteforce: build K-D index")
        if (.not. is_ok(ierr)) return

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, radius, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_large_1d_bruteforce: mask call")
        if (.not. is_ok(ierr)) return

        call assert_true(all(vicinity_mask .eqv. expected_mask), &
                         "test_vicinity_vectors_large_1d_bruteforce: brute-force mask")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, radius, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_true(is_ok(ierr), &
                         "test_vicinity_vectors_large_1d_bruteforce: count call")
        if (.not. is_ok(ierr)) return

        call assert_equal_int(neighbor_count, expected_count, &
                              "test_vicinity_vectors_large_1d_bruteforce: brute-force count")

    end subroutine test_vicinity_vectors_large_1d_bruteforce

    !> Test rejection of a negative radius.
    subroutine test_vicinity_vectors_negative_radius()

        integer(int32), parameter :: n_dimensions = 1_int32, n_points = 3_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(1, :) = [0.0_real64, 1.0_real64, 2.0_real64]
        query_point = [0.0_real64]
        dimension_order = [1_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, -1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_negative_radius: mask rejects radius")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, -1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_negative_radius: count rejects radius")

    end subroutine test_vicinity_vectors_negative_radius

    !> Test rejection of a dimension-order value below one.
    subroutine test_vicinity_vectors_dimension_order_low()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 3_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(:, 1) = [0.0_real64, 0.0_real64]
        points(:, 2) = [1.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 1.0_real64]

        query_point = [0.0_real64, 0.0_real64]
        dimension_order = [0_int32, 2_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_dimension_order_low: mask rejects order")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_dimension_order_low: count rejects order")

    end subroutine test_vicinity_vectors_dimension_order_low

    !> Test rejection of a dimension-order value above n_dimensions.
    subroutine test_vicinity_vectors_dimension_order_high()

        integer(int32), parameter :: n_dimensions = 2_int32, n_points = 3_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(:, 1) = [0.0_real64, 0.0_real64]
        points(:, 2) = [1.0_real64, 0.0_real64]
        points(:, 3) = [0.0_real64, 1.0_real64]

        query_point = [0.0_real64, 0.0_real64]
        dimension_order = [1_int32, 3_int32]
        kd_indices = [1_int32, 2_int32, 3_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_dimension_order_high: mask rejects order")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_dimension_order_high: count rejects order")

    end subroutine test_vicinity_vectors_dimension_order_high

    !> Test rejection of a K-D index below one.
    subroutine test_vicinity_vectors_kd_index_low()

        integer(int32), parameter :: n_dimensions = 1_int32, n_points = 3_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(1, :) = [0.0_real64, 1.0_real64, 2.0_real64]
        query_point = [0.0_real64]
        dimension_order = [1_int32]
        kd_indices = [1_int32, 0_int32, 3_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_kd_index_low: mask rejects K-D index")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_kd_index_low: count rejects K-D index")

    end subroutine test_vicinity_vectors_kd_index_low

    !> Test rejection of a K-D index above n_points.
    subroutine test_vicinity_vectors_kd_index_high()

        integer(int32), parameter :: n_dimensions = 1_int32, n_points = 3_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        points(1, :) = [0.0_real64, 1.0_real64, 2.0_real64]
        query_point = [0.0_real64]
        dimension_order = [1_int32]
        kd_indices = [1_int32, 4_int32, 3_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_kd_index_high: mask rejects K-D index")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_kd_index_high: count rejects K-D index")

    end subroutine test_vicinity_vectors_kd_index_high

    !> Test rejection of an empty reference-point collection.
    subroutine test_vicinity_vectors_zero_points()

        integer(int32), parameter :: n_dimensions = 1_int32, n_points = 0_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        query_point = [0.0_real64]
        dimension_order = [1_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_zero_points: mask rejects zero points")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_zero_points: count rejects zero points")

    end subroutine test_vicinity_vectors_zero_points

    !> Test rejection of a zero-dimensional coordinate space.
    subroutine test_vicinity_vectors_zero_dimensions()

        integer(int32), parameter :: n_dimensions = 0_int32, n_points = 1_int32
        real(real64) :: points(n_dimensions, n_points), query_point(n_dimensions)
        integer(int32) :: dimension_order(n_dimensions), kd_indices(n_points)
        integer(int32) :: neighbor_count, ierr
        logical :: vicinity_mask(n_points)

        kd_indices = [1_int32]

        call vicinity_vectors_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                    dimension_order, kd_indices, vicinity_mask, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_zero_dimensions: mask rejects zero dimensions")

        call vicinity_vectors_count_alloc(query_point, points, n_dimensions, n_points, 1.0_real64, &
                                          dimension_order, kd_indices, neighbor_count, ierr)

        call assert_false(is_ok(ierr), &
                          "test_vicinity_vectors_zero_dimensions: count rejects zero dimensions")

    end subroutine test_vicinity_vectors_zero_dimensions

    !> Helper: Generate random unit vectors.
    subroutine random_unit_vectors(V, d, n)
        integer(int32), intent(in) :: d, n
        real(real64), intent(out) :: V(d, n)
        integer(int32) :: i
        real(real64) :: norm
        call random_seed()
        do i = 1, n
            call random_number(V(:, i))
            V(:, i) = V(:, i) - 0.5d0
            norm = sqrt(sum(V(:, i)**2))
            if (norm > 0) V(:, i) = V(:, i)/norm
        end do
    end subroutine random_unit_vectors

    !> Helper: Generate random matrix.
    subroutine random_matrix(X, d, n)
        integer(int32), intent(in) :: d, n
        real(real64), intent(out) :: X(d, n)
        integer(int32) :: j
        call random_seed()
        do j = 1, n
            call random_number(X(:, j))
        end do
    end subroutine random_matrix

    !> Helper: Convert integer to string.
    function str(i) result(s)
        integer(int32), intent(in) :: i
        character(len=32) :: s
        write (s, *) i
        s = adjustl(s)
    end function str

end module mod_test_kd_tree
