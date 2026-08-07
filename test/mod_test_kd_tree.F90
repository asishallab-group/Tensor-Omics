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
        allocate (all_tests(21))

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
        all_tests(13) = test_case("test_kd_knn_query_basic", test_kd_knn_query_basic)
        all_tests(14) = test_case("test_kd_knn_query_k_equals_n", test_kd_knn_query_k_equals_n)
        all_tests(15) = test_case("test_kd_knn_query_invalid_k", test_kd_knn_query_invalid_k)
        all_tests(16) = test_case("test_kd_range_query_mask_basic", test_kd_range_query_mask_basic)
        all_tests(17) = test_case("test_kd_range_query_mask_negative_radius", test_kd_range_query_mask_negative_radius)
        all_tests(18) = test_case("test_kd_range_query_list_basic", test_kd_range_query_list_basic)
        all_tests(19) = test_case("test_kd_range_query_list_negative_radius", test_kd_range_query_list_negative_radius)
        all_tests(20) = test_case("test_kd_range_query_count_basic", test_kd_range_query_count_basic)
        all_tests(21) = test_case("test_kd_range_query_count_matches_list", test_kd_range_query_count_matches_list)
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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')
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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')
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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')
        call get_kd_point(X, kd_ix, 4, val, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')

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
        call assert_equal_int(get_err_code(ierr), ERR_OK, 'Unexpected error')
        call assert_permutation(kd_ix, n, "5D medium KD-Tree")
    end subroutine test_kd_5d_medium

    ! --- kd_knn_query / kd_range_query_mask / kd_range_query_list / kd_range_query_count ----
    !
    ! Shared fixture for these tests: D=2, N=11 points on a line, (0,0),(1,0),...,(10,0), so
    ! point i (1-indexed) has x=i-1. Distances and radius memberships are then trivial to
    ! hand-verify.

    !> Build the shared 11-point line fixture and its k-d tree.
    subroutine build_line_fixture(vectors, kd_indices, dim_order)
        real(real64), intent(out) :: vectors(2, 11)
        integer(int32), intent(out) :: kd_indices(11)
        integer(int32), intent(out) :: dim_order(2)
        integer(int32) :: i, ierr

        do i = 1, 11
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        dim_order = [1, 2]

        call build_kd_index_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'build_line_fixture: build_kd_index_alloc failed: ', ierr
            error stop
        end if
    end subroutine build_line_fixture

    !> Query at x=5 (point index 6), k=3: nearest are x=4,5,6 (indices 5,6,7), distances 1,0,1.
    subroutine test_kd_knn_query_basic()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbors(3), ierr
        real(real64)   :: distances(3)
        logical        :: found(11)
        integer(int32) :: i

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_knn_query_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                [5.0d0, 0.0d0], 3_int32, neighbors, distances, ierr)
        call assert_true(is_ok(ierr), "kd_knn_query_alloc should succeed")

        found = .false.
        do i = 1, 3
            found(neighbors(i)) = .true.
        end do
        call assert_true(found(5) .and. found(6) .and. found(7), &
                         "kd_knn_query: nearest 3 to x=5 should be indices {5,6,7}")
        call assert_equal_real(minval(distances), 0.0d0, 1.0d-9, "kd_knn_query: closest neighbor is exact match")
        call assert_equal_real(maxval(distances), 1.0d0, 1.0d-9, "kd_knn_query: farthest of the 3 is at distance 1")
    end subroutine test_kd_knn_query_basic

    !> k=n_points must return every point exactly once.
    subroutine test_kd_knn_query_k_equals_n()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbors(11), ierr
        real(real64)   :: distances(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_knn_query_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                [5.0d0, 0.0d0], 11_int32, neighbors, distances, ierr)
        call assert_true(is_ok(ierr), "kd_knn_query_alloc with k=n_points should succeed")
        call assert_permutation(neighbors, 11_int32, "kd_knn_query with k=n_points covers every point once")
    end subroutine test_kd_knn_query_k_equals_n

    !> k=0 and k>n_points must both be rejected by validation.
    subroutine test_kd_knn_query_invalid_k()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbors(11), ierr
        real(real64)   :: distances(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_knn_query_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                [5.0d0, 0.0d0], 0_int32, neighbors(1:0), distances(1:0), ierr)
        call assert_true(is_err(ierr), "kd_knn_query should reject k_neighbors=0")

        call kd_knn_query_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                [5.0d0, 0.0d0], 12_int32, neighbors, distances, ierr)
        call assert_true(is_err(ierr), "kd_knn_query should reject k_neighbors > n_points")
    end subroutine test_kd_knn_query_invalid_k

    !> Query at x=5, radius=1.5: matches x=4,5,6 (indices 5,6,7), nothing else.
    subroutine test_kd_range_query_mask_basic()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        logical        :: in_radius_mask(11), expected(11)
        integer(int32) :: ierr

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_range_query_mask_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                       [5.0d0, 0.0d0], 1.5d0, in_radius_mask, ierr)
        call assert_true(is_ok(ierr), "kd_range_query_mask_alloc should succeed")

        expected = .false.
        expected(5:7) = .true.
        call assert_equal_array_logical(in_radius_mask, expected, 11_int32, "kd_range_query_mask at x=5, radius=1.5")
    end subroutine test_kd_range_query_mask_basic

    !> A negative radius must be rejected by validation.
    subroutine test_kd_range_query_mask_negative_radius()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        logical        :: in_radius_mask(11)
        integer(int32) :: ierr

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_range_query_mask_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                       [5.0d0, 0.0d0], -1.0d0, in_radius_mask, ierr)
        call assert_true(is_err(ierr), "kd_range_query_mask should reject a negative radius")
    end subroutine test_kd_range_query_mask_negative_radius

    !> Same fixture/radius as the mask test: the compact list must name exactly {5,6,7}.
    subroutine test_kd_range_query_list_basic()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbors(11), n_found, ierr
        logical        :: found(11)
        integer(int32) :: i

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_range_query_list_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                       [5.0d0, 0.0d0], 1.5d0, neighbors, n_found, ierr)
        call assert_true(is_ok(ierr), "kd_range_query_list_alloc should succeed")
        call assert_equal_int(n_found, 3, "kd_range_query_list: 3 points within radius 1.5 of x=5")

        found = .false.
        do i = 1, n_found
            found(neighbors(i)) = .true.
        end do
        call assert_true(found(5) .and. found(6) .and. found(7), &
                         "kd_range_query_list: found points should be exactly indices {5,6,7}")
    end subroutine test_kd_range_query_list_basic

    !> A negative radius must be rejected by validation.
    subroutine test_kd_range_query_list_negative_radius()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbors(11), n_found, ierr

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_range_query_list_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                       [5.0d0, 0.0d0], -1.0d0, neighbors, n_found, ierr)
        call assert_true(is_err(ierr), "kd_range_query_list should reject a negative radius")
    end subroutine test_kd_range_query_list_negative_radius

    !> Same fixture/radius: the count-only form must agree with the mask/list forms (3).
    subroutine test_kd_range_query_count_basic()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbor_count, ierr

        call build_line_fixture(vectors, kd_indices, dim_order)

        call kd_range_query_count_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                        [5.0d0, 0.0d0], 1.5d0, neighbor_count, ierr)
        call assert_true(is_ok(ierr), "kd_range_query_count_alloc should succeed")
        call assert_equal_int(neighbor_count, 3, "kd_range_query_count: 3 points within radius 1.5 of x=5")
    end subroutine test_kd_range_query_count_basic

    !> Cross-check: for a range of query points/radii, count must always equal list's n_found
    !| and the number of .true. entries in the mask -- the three query forms must agree.
    subroutine test_kd_range_query_count_matches_list()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2)
        integer(int32) :: neighbors(11), n_found, neighbor_count, ierr
        logical        :: in_radius_mask(11)
        real(real64)   :: query_x
        integer(int32) :: qi

        call build_line_fixture(vectors, kd_indices, dim_order)

        do qi = 0, 10
            query_x = real(qi, real64)

            call kd_range_query_list_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                           [query_x, 0.0d0], 2.5d0, neighbors, n_found, ierr)
            call assert_true(is_ok(ierr), "kd_range_query_list_alloc should succeed")

            call kd_range_query_count_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                            [query_x, 0.0d0], 2.5d0, neighbor_count, ierr)
            call assert_true(is_ok(ierr), "kd_range_query_count_alloc should succeed")

            call kd_range_query_mask_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                           [query_x, 0.0d0], 2.5d0, in_radius_mask, ierr)
            call assert_true(is_ok(ierr), "kd_range_query_mask_alloc should succeed")

            call assert_equal_int(neighbor_count, n_found, "count must match list's n_found")
            call assert_equal_int(neighbor_count, count(in_radius_mask, kind=int32), "count must match mask's true count")
        end do
    end subroutine test_kd_range_query_count_matches_list

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
