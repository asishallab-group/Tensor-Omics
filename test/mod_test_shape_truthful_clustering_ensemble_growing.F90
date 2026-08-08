!> Unit test suite for tox_shape_truthful_clustering_ensemble_growing
!| (calc_ensemble_growth_radius, grow_ensemble), generated from
!| src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_ensemble_growing_kernel.F90.
module mod_test_shape_truthful_clustering_ensemble_growing
    use tox_shape_truthful_clustering_ensemble_growing, only: calc_ensemble_growth_radius_alloc, grow_ensemble_alloc
    use f42_kd_tree, only: build_kd_index_alloc
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_ensemble_growing() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(14))

        all_tests(1) = test_case("test_growth_radius_even_k", test_growth_radius_even_k)
        all_tests(2) = test_case("test_growth_radius_odd_k", test_growth_radius_odd_k)
        all_tests(3) = test_case("test_growth_radius_seed_index_out_of_range", test_growth_radius_seed_index_out_of_range)
        all_tests(4) = test_case("test_growth_radius_k_min_too_large", test_growth_radius_k_min_too_large)
        all_tests(5) = test_case("test_growth_radius_omitted_k_min_is_clamped", test_growth_radius_omitted_k_min_is_clamped)
        all_tests(6) = test_case("test_growth_radius_zero_dimensions", test_growth_radius_zero_dimensions)
        all_tests(7) = test_case("test_growth_radius_percentile_min", test_growth_radius_percentile_min)
        all_tests(8) = test_case("test_growth_radius_percentile_max", test_growth_radius_percentile_max)
        all_tests(9) = test_case("test_growth_radius_invalid_percentile", test_growth_radius_invalid_percentile)
        all_tests(10) = test_case("test_grow_ensemble_single_member", test_grow_ensemble_single_member)
        all_tests(11) = test_case("test_grow_ensemble_multi_member_union", test_grow_ensemble_multi_member_union)
        all_tests(12) = test_case("test_grow_ensemble_empty_ensemble_is_degenerate", test_grow_ensemble_empty_ensemble_is_degenerate)
        all_tests(13) = test_case("test_grow_ensemble_negative_radius", test_grow_ensemble_negative_radius)
        all_tests(14) = test_case("test_grow_ensemble_zero_dimensions", test_grow_ensemble_zero_dimensions)
    end function get_all_tests_shape_truthful_clustering_ensemble_growing

    ! --- calc_ensemble_growth_radius ---------------------------------------
    !
    ! Shared fixture: D=2, N=11 points on a line, (0,0),(1,0),...,(10,0). Seed at x=5 (index
    ! 6). Its 4 nearest neighbors are x=4,6 (distance 1) and x=3,7 (distance 2), sorted
    ! [1,1,2,2] -- even k_min=4, median = avg(2nd,3rd) = (1+2)/2 = 1.5. Its 3 nearest are
    ! x=4,6 (distance 1) and one of x=3/x=7 (distance 2, a tie), sorted [1,1,2] -- odd k_min=3,
    ! median = the middle (2nd) element = 1.

    subroutine test_growth_radius_even_k()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=4_int32, growth_radius=growth_radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calc_ensemble_growth_radius_alloc failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(growth_radius, 1.5d0, 1.0d-9, "calc_ensemble_growth_radius: even k_min=4")
    end subroutine test_growth_radius_even_k

    subroutine test_growth_radius_odd_k()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=3_int32, growth_radius=growth_radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calc_ensemble_growth_radius_alloc failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(growth_radius, 1.0d0, 1.0d-9, "calc_ensemble_growth_radius: odd k_min=3")
    end subroutine test_growth_radius_odd_k

    subroutine test_growth_radius_seed_index_out_of_range()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 12_int32, &
                                               growth_radius=growth_radius, ierr=ierr)
        call assert_true(is_err(ierr), "calc_ensemble_growth_radius should reject seed_index > n_vectors")
    end subroutine test_growth_radius_seed_index_out_of_range

    subroutine test_growth_radius_k_min_too_large()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=11_int32, growth_radius=growth_radius, ierr=ierr)
        call assert_true(is_err(ierr), "calc_ensemble_growth_radius should reject k_min > n_vectors-1")
    end subroutine test_growth_radius_k_min_too_large

    !> Regression test for a genuine crash: calling `calc_ensemble_growth_radius_alloc` with
    !| `k_min` *omitted* on this 11-point fixture (default 30, but DM_MAX(n_vectors - 1) is
    !| only 10) used to corrupt memory -- see `misc/code_gen_footgun.md`'s third entry for the
    !| generator-level root cause and `calc_ensemble_growth_radius_kernel`'s own
    !| `min(actual_k_min, n_vectors - 1)` clamp for the fix. Asserts the fix itself: omitting
    !| `k_min` here must resolve to exactly `n_vectors - 1 = 10`, matching an explicit call.
    subroutine test_growth_radius_omitted_k_min_is_clamped()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: radius_omitted, radius_explicit

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               growth_radius=radius_omitted, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calc_ensemble_growth_radius_alloc (k_min omitted) failed unexpectedly: ', ierr
            error stop
        end if

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=10_int32, growth_radius=radius_explicit, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calc_ensemble_growth_radius_alloc (k_min=10) failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_real(radius_omitted, radius_explicit, 1.0d-12, &
                               "calc_ensemble_growth_radius: an omitted k_min on N=11 clamps to exactly k_min=10")
    end subroutine test_growth_radius_omitted_k_min_is_clamped

    !> Same k_min=4 fixture as test_growth_radius_even_k, distances sorted [1,1,2,2].
    !| radius_percentile=0.0 is calc_percentile_helper's own "below the first index" edge
    !| case, which returns the smallest value in sorted order -- the nearest-neighbor
    !| distance, 1.0 -- irrespective of k_min's parity.
    subroutine test_growth_radius_percentile_min()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=4_int32, radius_percentile=0.0d0, &
                                               growth_radius=growth_radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calc_ensemble_growth_radius_alloc failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(growth_radius, 1.0d0, 1.0d-9, "calc_ensemble_growth_radius: radius_percentile=0.0 is the min")
    end subroutine test_growth_radius_percentile_min

    !> Same fixture, radius_percentile=100.0 -- the largest value in sorted order, the
    !| farthest of the k_min neighbors, 2.0.
    subroutine test_growth_radius_percentile_max()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=4_int32, radius_percentile=100.0d0, &
                                               growth_radius=growth_radius, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'calc_ensemble_growth_radius_alloc failed unexpectedly: ', ierr
            error stop
        end if
        call assert_equal_real(growth_radius, 2.0d0, 1.0d-9, "calc_ensemble_growth_radius: radius_percentile=100.0 is the max")
    end subroutine test_growth_radius_percentile_max

    subroutine test_growth_radius_invalid_percentile()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        real(real64)   :: growth_radius

        call build_line_fixture(vectors, kd_indices, dim_order)

        call calc_ensemble_growth_radius_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               k_min=4_int32, radius_percentile=101.0d0, &
                                               growth_radius=growth_radius, ierr=ierr)
        call assert_true(is_err(ierr), "calc_ensemble_growth_radius should reject radius_percentile > 100")
    end subroutine test_growth_radius_invalid_percentile

    subroutine test_growth_radius_zero_dimensions()
        real(real64)   :: vectors(0, 11)
        integer(int32) :: kd_indices(11), dim_order(0), ierr
        real(real64)   :: growth_radius

        call calc_ensemble_growth_radius_alloc(vectors, 0_int32, 11_int32, kd_indices, dim_order, 6_int32, &
                                               growth_radius=growth_radius, ierr=ierr)
        call assert_true(is_err(ierr), "calc_ensemble_growth_radius should reject n_dimensions=0")
    end subroutine test_growth_radius_zero_dimensions

    ! --- grow_ensemble -------------------------------------------------------
    !
    ! Same 11-point line fixture. With growth_radius=1.5: from a single member x=5 (index 6):
    ! covers x=4,5,6 (indices 5,6,7). From members x=4,5,6 (indices 5,6,7): the union of each
    ! member's own 1.5-radius neighborhood covers x=3..7 (indices 4..8).

    subroutine test_grow_ensemble_single_member()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        logical        :: is_member_mask(11), is_member_mask_next(11), expected(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        is_member_mask = .false.
        is_member_mask(6) = .true.

        call grow_ensemble_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                 is_member_mask, 1.5d0, is_member_mask_next, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'grow_ensemble_alloc failed unexpectedly: ', ierr
            error stop
        end if

        expected = .false.
        expected(5:7) = .true.
        call assert_equal_array_logical(is_member_mask_next, expected, 11_int32, "grow_ensemble from a single member")
    end subroutine test_grow_ensemble_single_member

    subroutine test_grow_ensemble_multi_member_union()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        logical        :: is_member_mask(11), is_member_mask_next(11), expected(11)

        call build_line_fixture(vectors, kd_indices, dim_order)

        is_member_mask = .false.
        is_member_mask(5:7) = .true.

        call grow_ensemble_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                 is_member_mask, 1.5d0, is_member_mask_next, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'grow_ensemble_alloc failed unexpectedly: ', ierr
            error stop
        end if

        expected = .false.
        expected(4:8) = .true.
        call assert_equal_array_logical(is_member_mask_next, expected, 11_int32, "grow_ensemble from three members (union)")
    end subroutine test_grow_ensemble_multi_member_union

    !> Under the new kernel-path validation, an all-.false. is_member_mask is a well-defined
    !| degenerate case (nothing to grow from), not a validation error -- unlike the old
    !| hand-written implementation, which rejected it explicitly. There is no clean DM_*
    !| annotation to express "at least one .true." for a plain logical array, and the
    !| computation is already well-defined for this input, so it is intentionally left
    !| unvalidated rather than forcing a bespoke runtime-error check for it.
    subroutine test_grow_ensemble_empty_ensemble_is_degenerate()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        logical        :: is_member_mask(11), is_member_mask_next(11)

        call build_line_fixture(vectors, kd_indices, dim_order)
        is_member_mask = .false.

        call grow_ensemble_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                 is_member_mask, 1.5d0, is_member_mask_next, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'grow_ensemble_alloc failed unexpectedly: ', ierr
            error stop
        end if
        call assert_true(.not. any(is_member_mask_next), "grow_ensemble: an empty ensemble stays empty")
    end subroutine test_grow_ensemble_empty_ensemble_is_degenerate

    subroutine test_grow_ensemble_negative_radius()
        real(real64)   :: vectors(2, 11)
        integer(int32) :: kd_indices(11), dim_order(2), ierr
        logical        :: is_member_mask(11), is_member_mask_next(11)

        call build_line_fixture(vectors, kd_indices, dim_order)
        is_member_mask = .false.
        is_member_mask(6) = .true.

        call grow_ensemble_alloc(vectors, 2_int32, 11_int32, kd_indices, dim_order, &
                                 is_member_mask, -1.5d0, is_member_mask_next, ierr)
        call assert_true(is_err(ierr), "grow_ensemble should reject a negative growth radius")
    end subroutine test_grow_ensemble_negative_radius

    subroutine test_grow_ensemble_zero_dimensions()
        real(real64)   :: vectors(0, 11)
        integer(int32) :: kd_indices(11), dim_order(0), ierr
        logical        :: is_member_mask(11), is_member_mask_next(11)

        is_member_mask = .false.
        is_member_mask(6) = .true.

        call grow_ensemble_alloc(vectors, 0_int32, 11_int32, kd_indices, dim_order, &
                                 is_member_mask, 1.5d0, is_member_mask_next, ierr)
        call assert_true(is_err(ierr), "grow_ensemble should reject n_dimensions=0")
    end subroutine test_grow_ensemble_zero_dimensions

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

end module mod_test_shape_truthful_clustering_ensemble_growing
