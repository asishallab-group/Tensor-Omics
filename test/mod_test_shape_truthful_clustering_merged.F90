!> Unit test suite for tox_shape_truthful_clustering (ensemble_identification_merged),
!| generated from
!| src/tox/shape_truthful_clustering/tox_shape_truthful_clustering_impl.F90.
module mod_test_shape_truthful_clustering_merged
    use tox_shape_truthful_clustering, only: ensemble_identification_merged
    use tox_shape_truthful_clustering_impl, only: STOP_REASON_FIXED_POINT, &
                                                    MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_SEED
    use f42_kd_tree, only: build_kd_index
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_merged() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(5))

        all_tests(1) = test_case("test_merged_single_seed_matches_per_seed_impl", &
                                 test_merged_single_seed_matches_per_seed_impl)
        all_tests(2) = test_case("test_merged_two_independent_seeds", test_merged_two_independent_seeds)
        all_tests(3) = test_case("test_merged_zero_seeds", test_merged_zero_seeds)
        all_tests(4) = test_case("test_merged_seed_count_mismatch", test_merged_seed_count_mismatch)
        all_tests(5) = test_case("test_merged_n_dimensions_too_small", test_merged_n_dimensions_too_small)
    end function get_all_tests_shape_truthful_clustering_merged

    !> D=2, N=7. A 5-point line (0,0)..(4,0), plus two far-away points a growth radius of
    !| 1.0 (k_min=1) never reaches -- matches fixture A in
    !| test/mod_test_shape_truthful_clustering.F90 exactly, so its expected trajectory
    !| (bootstrap->t=2->t=3->t=4->fixed point at t=5) is already worked out there.
    subroutine build_fixture_a(vectors, kd_indices, dim_order)
        real(real64), intent(out) :: vectors(2, 7)
        integer(int32), intent(out) :: kd_indices(7)
        integer(int32), intent(out) :: dim_order(2)
        integer(int32) :: i, ierr

        do i = 1, 5
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        vectors(:, 6) = [0.0d0, 1.5d0]
        vectors(:, 7) = [0.0d0, 3.0d0]
        dim_order = [1, 2]

        call build_kd_index(vectors, 2_int32, 7_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'build_fixture_a: build_kd_index failed: ', ierr
            error stop
        end if
    end subroutine build_fixture_a

    subroutine test_merged_single_seed_matches_per_seed_impl()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: seed_selection_mask(7)
        logical(c_bool)        :: ensemble_masks(7, 1), expected_mask(7)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(4, 1), ensemble_k_history(4, 1)
        integer(int32) :: ensemble_member_added_at_step(7, 1), expected_step(7)
        logical(c_bool)        :: ensemble_low_confidence_masks(7, 1)
        logical(c_bool)        :: ensemble_accepted_history(4, 1)
        real(real64)   :: ensemble_growth_radii(1), ensemble_G_history(4, 1), ensemble_mu_history(2, 4, 1)
        real(real64)   :: ensemble_S_history(2, 4, 1), ensemble_U_history(2, 2, 4, 1)
        real(real64)   :: ensemble_U_first(2, 2, 1)
        integer(int32) :: ensemble_d_first(1)

        call build_fixture_a(vectors, kd_indices, dim_order)
        seed_selection_mask = .false.
        seed_selection_mask(1) = .true.

        call ensemble_identification_merged(vectors, 2_int32, 7_int32, kd_indices, dim_order, seed_selection_mask, 1_int32, &
                                            k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=4_int32, &
                                            ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                            ensemble_growth_radii=ensemble_growth_radii, &
                                            ensemble_U_history=ensemble_U_history, ensemble_S_history=ensemble_S_history, &
                                            ensemble_d_history=ensemble_d_history, ensemble_G_history=ensemble_G_history, &
                                            ensemble_mu_history=ensemble_mu_history, ensemble_k_history=ensemble_k_history, &
                                            ensemble_accepted_history=ensemble_accepted_history, &
                                            ensemble_member_added_at_step=ensemble_member_added_at_step, &
                                            ensemble_low_confidence_masks=ensemble_low_confidence_masks, ensemble_U_first=ensemble_U_first, ensemble_d_first=ensemble_d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification_merged failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(ensemble_stop_reason(1), STOP_REASON_FIXED_POINT, "merged single seed: stop_reason")
        call assert_equal_real(ensemble_growth_radii(1), 1.0d0, 1.0d-9, "merged single seed: growth_radius")

        expected_mask = .false.
        expected_mask(1:5) = .true.
        call assert_equal_array_logical(ensemble_masks(:, 1), expected_mask, 7_int32, "merged single seed: mask")

        call assert_equal_array_int(ensemble_k_history(:, 1), [2, 3, 4, 5], 4_int32, "merged single seed: k_history")
        call assert_true(all(ensemble_accepted_history(:, 1)), "merged single seed: all accepted")

        expected_step = [MEMBER_ADDED_AT_STEP_SEED, 1, 2, 3, 4, MEMBER_ADDED_AT_STEP_NON_MEMBER, &
                         MEMBER_ADDED_AT_STEP_NON_MEMBER]
        call assert_equal_array_int(ensemble_member_added_at_step(:, 1), expected_step, 7_int32, &
                                    "merged single seed: member_added_at_step")

        ! Iteration 1's own bootstrap mask -- {seed=1, its one growth-radius neighbor=2}.
        expected_mask = .false.
        expected_mask(1:2) = .true.
        call assert_equal_array_logical(ensemble_low_confidence_masks(:, 1), expected_mask, 7_int32, &
                                        "merged single seed: low_confidence_mask is iteration 1's own mask")
    end subroutine test_merged_single_seed_matches_per_seed_impl

    !> Two independent copies of fixture A, translated far apart (+100,+100 in the second),
    !| seeded at both copies' own index-1 point. Each column must reproduce fixture A's
    !| trajectory independently -- neither seed's growth or history may leak into the
    !| other's column.
    subroutine test_merged_two_independent_seeds()
        real(real64)   :: vectors(2, 14)
        integer(int32) :: kd_indices(14), dim_order(2), ierr, i
        logical(c_bool)        :: seed_selection_mask(14)
        logical(c_bool)        :: ensemble_masks(14, 2), expected_mask_1(14), expected_mask_2(14)
        integer(int32) :: ensemble_stop_reason(2), ensemble_d_history(4, 2), ensemble_k_history(4, 2)
        integer(int32) :: ensemble_member_added_at_step(14, 2), expected_step_1(14), expected_step_2(14)
        logical(c_bool)        :: ensemble_low_confidence_masks(14, 2)
        logical(c_bool)        :: ensemble_accepted_history(4, 2)
        real(real64)   :: ensemble_growth_radii(2), ensemble_G_history(4, 2), ensemble_mu_history(2, 4, 2)
        real(real64)   :: ensemble_S_history(2, 4, 2), ensemble_U_history(2, 2, 4, 2)
        real(real64)   :: ensemble_U_first(2, 2, 2)
        integer(int32) :: ensemble_d_first(2)

        do i = 1, 5
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
        end do
        vectors(:, 6) = [0.0d0, 1.5d0]
        vectors(:, 7) = [0.0d0, 3.0d0]
        vectors(:, 8:14) = vectors(:, 1:7)
        vectors(1, 8:14) = vectors(1, 8:14) + 100.0d0
        vectors(2, 8:14) = vectors(2, 8:14) + 100.0d0
        dim_order = [1, 2]

        call build_kd_index(vectors, 2_int32, 14_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_merged_two_independent_seeds: build_kd_index failed: ', ierr
            error stop
        end if

        seed_selection_mask = .false.
        seed_selection_mask(1) = .true.
        seed_selection_mask(8) = .true.

        call ensemble_identification_merged(vectors, 2_int32, 14_int32, kd_indices, dim_order, seed_selection_mask, 2_int32, &
                                            k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=4_int32, &
                                            ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                            ensemble_growth_radii=ensemble_growth_radii, &
                                            ensemble_U_history=ensemble_U_history, ensemble_S_history=ensemble_S_history, &
                                            ensemble_d_history=ensemble_d_history, ensemble_G_history=ensemble_G_history, &
                                            ensemble_mu_history=ensemble_mu_history, ensemble_k_history=ensemble_k_history, &
                                            ensemble_accepted_history=ensemble_accepted_history, &
                                            ensemble_member_added_at_step=ensemble_member_added_at_step, &
                                            ensemble_low_confidence_masks=ensemble_low_confidence_masks, ensemble_U_first=ensemble_U_first, ensemble_d_first=ensemble_d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification_merged failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(ensemble_stop_reason(1), STOP_REASON_FIXED_POINT, "two seeds: column 1 stop_reason")
        call assert_equal_int(ensemble_stop_reason(2), STOP_REASON_FIXED_POINT, "two seeds: column 2 stop_reason")

        expected_mask_1 = .false.
        expected_mask_1(1:5) = .true.
        call assert_equal_array_logical(ensemble_masks(:, 1), expected_mask_1, 14_int32, "two seeds: column 1 mask")

        expected_mask_2 = .false.
        expected_mask_2(8:12) = .true.
        call assert_equal_array_logical(ensemble_masks(:, 2), expected_mask_2, 14_int32, "two seeds: column 2 mask")

        call assert_equal_array_int(ensemble_k_history(:, 1), [2, 3, 4, 5], 4_int32, "two seeds: column 1 k_history")
        call assert_equal_array_int(ensemble_k_history(:, 2), [2, 3, 4, 5], 4_int32, "two seeds: column 2 k_history")

        expected_step_1 = MEMBER_ADDED_AT_STEP_NON_MEMBER
        expected_step_1(1) = MEMBER_ADDED_AT_STEP_SEED
        expected_step_1(2:5) = [1, 2, 3, 4]
        call assert_equal_array_int(ensemble_member_added_at_step(:, 1), expected_step_1, 14_int32, &
                                    "two seeds: column 1 member_added_at_step")

        expected_step_2 = MEMBER_ADDED_AT_STEP_NON_MEMBER
        expected_step_2(8) = MEMBER_ADDED_AT_STEP_SEED
        expected_step_2(9:12) = [1, 2, 3, 4]
        call assert_equal_array_int(ensemble_member_added_at_step(:, 2), expected_step_2, 14_int32, &
                                    "two seeds: column 2 member_added_at_step")

        ! ensemble_U_first/ensemble_d_first: each column is its own seed's bootstrap basis,
        ! collinear along the x-axis in both copies -- must not leak across columns.
        call assert_equal_int(ensemble_d_first(1), 1_int32, "two seeds: column 1 d_first")
        call assert_equal_int(ensemble_d_first(2), 1_int32, "two seeds: column 2 d_first")
        call assert_equal_real(abs(ensemble_U_first(1,1,1)), 1.0d0, 1.0d-9, "two seeds: column 1 U_first, x")
        call assert_equal_real(abs(ensemble_U_first(2,1,1)), 0.0d0, 1.0d-9, "two seeds: column 1 U_first, y")
        call assert_equal_real(abs(ensemble_U_first(1,1,2)), 1.0d0, 1.0d-9, "two seeds: column 2 U_first, x")
        call assert_equal_real(abs(ensemble_U_first(2,1,2)), 0.0d0, 1.0d-9, "two seeds: column 2 U_first, y")
    end subroutine test_merged_two_independent_seeds

    subroutine test_merged_zero_seeds()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: seed_selection_mask(7)
        logical(c_bool)        :: ensemble_masks(7, 0)
        integer(int32) :: ensemble_stop_reason(0), ensemble_d_history(4, 0), ensemble_k_history(4, 0)
        integer(int32) :: ensemble_member_added_at_step(7, 0)
        logical(c_bool)        :: ensemble_low_confidence_masks(7, 0)
        logical(c_bool)        :: ensemble_accepted_history(4, 0)
        real(real64)   :: ensemble_growth_radii(0), ensemble_G_history(4, 0), ensemble_mu_history(2, 4, 0)
        real(real64)   :: ensemble_S_history(2, 4, 0), ensemble_U_history(2, 2, 4, 0)
        real(real64)   :: ensemble_U_first(2, 2, 0)
        integer(int32) :: ensemble_d_first(0)

        call build_fixture_a(vectors, kd_indices, dim_order)
        seed_selection_mask = .false.

        call ensemble_identification_merged(vectors, 2_int32, 7_int32, kd_indices, dim_order, seed_selection_mask, 0_int32, &
                                            k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=4_int32, &
                                            ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                            ensemble_growth_radii=ensemble_growth_radii, &
                                            ensemble_U_history=ensemble_U_history, ensemble_S_history=ensemble_S_history, &
                                            ensemble_d_history=ensemble_d_history, ensemble_G_history=ensemble_G_history, &
                                            ensemble_mu_history=ensemble_mu_history, ensemble_k_history=ensemble_k_history, &
                                            ensemble_accepted_history=ensemble_accepted_history, &
                                            ensemble_member_added_at_step=ensemble_member_added_at_step, &
                                            ensemble_low_confidence_masks=ensemble_low_confidence_masks, ensemble_U_first=ensemble_U_first, ensemble_d_first=ensemble_d_first, ierr=ierr)
        call assert_true(is_ok(ierr), "merged zero seeds: not an error")
    end subroutine test_merged_zero_seeds

    subroutine test_merged_seed_count_mismatch()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: seed_selection_mask(7)
        logical(c_bool)        :: ensemble_masks(7, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(4, 1), ensemble_k_history(4, 1)
        integer(int32) :: ensemble_member_added_at_step(7, 1)
        logical(c_bool)        :: ensemble_low_confidence_masks(7, 1)
        logical(c_bool)        :: ensemble_accepted_history(4, 1)
        real(real64)   :: ensemble_growth_radii(1), ensemble_G_history(4, 1), ensemble_mu_history(2, 4, 1)
        real(real64)   :: ensemble_S_history(2, 4, 1), ensemble_U_history(2, 2, 4, 1)
        real(real64)   :: ensemble_U_first(2, 2, 1)
        integer(int32) :: ensemble_d_first(1)

        call build_fixture_a(vectors, kd_indices, dim_order)
        seed_selection_mask = .false.
        seed_selection_mask(1) = .true.
        seed_selection_mask(6) = .true.

        call ensemble_identification_merged(vectors, 2_int32, 7_int32, kd_indices, dim_order, seed_selection_mask, 1_int32, &
                                            k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=4_int32, &
                                            ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                            ensemble_growth_radii=ensemble_growth_radii, &
                                            ensemble_U_history=ensemble_U_history, ensemble_S_history=ensemble_S_history, &
                                            ensemble_d_history=ensemble_d_history, ensemble_G_history=ensemble_G_history, &
                                            ensemble_mu_history=ensemble_mu_history, ensemble_k_history=ensemble_k_history, &
                                            ensemble_accepted_history=ensemble_accepted_history, &
                                            ensemble_member_added_at_step=ensemble_member_added_at_step, &
                                            ensemble_low_confidence_masks=ensemble_low_confidence_masks, ensemble_U_first=ensemble_U_first, ensemble_d_first=ensemble_d_first, ierr=ierr)
        call assert_true(is_err(ierr), "merged: n_selected_seed must match count(seed_selection_mask)")
    end subroutine test_merged_seed_count_mismatch

    subroutine test_merged_n_dimensions_too_small()
        real(real64)   :: vectors(1, 7)
        integer(int32) :: kd_indices(7), dim_order(1), ierr, i
        logical(c_bool)        :: seed_selection_mask(7)
        logical(c_bool)        :: ensemble_masks(7, 1)
        integer(int32) :: ensemble_stop_reason(1), ensemble_d_history(4, 1), ensemble_k_history(4, 1)
        integer(int32) :: ensemble_member_added_at_step(7, 1)
        logical(c_bool)        :: ensemble_low_confidence_masks(7, 1)
        logical(c_bool)        :: ensemble_accepted_history(4, 1)
        real(real64)   :: ensemble_growth_radii(1), ensemble_G_history(4, 1), ensemble_mu_history(1, 4, 1)
        real(real64)   :: ensemble_S_history(1, 4, 1), ensemble_U_history(1, 1, 4, 1)
        real(real64)   :: ensemble_U_first(1, 1, 1)
        integer(int32) :: ensemble_d_first(1)

        do i = 1, 7
            vectors(1, i) = real(i - 1, real64)
        end do
        dim_order = [1]
        call build_kd_index(vectors, 1_int32, 7_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_merged_n_dimensions_too_small: build_kd_index failed: ', ierr
            error stop
        end if

        seed_selection_mask = .false.
        seed_selection_mask(1) = .true.

        call ensemble_identification_merged(vectors, 1_int32, 7_int32, kd_indices, dim_order, seed_selection_mask, 1_int32, &
                                            k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=4_int32, &
                                            ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                            ensemble_growth_radii=ensemble_growth_radii, &
                                            ensemble_U_history=ensemble_U_history, ensemble_S_history=ensemble_S_history, &
                                            ensemble_d_history=ensemble_d_history, ensemble_G_history=ensemble_G_history, &
                                            ensemble_mu_history=ensemble_mu_history, ensemble_k_history=ensemble_k_history, &
                                            ensemble_accepted_history=ensemble_accepted_history, &
                                            ensemble_member_added_at_step=ensemble_member_added_at_step, &
                                            ensemble_low_confidence_masks=ensemble_low_confidence_masks, ensemble_U_first=ensemble_U_first, ensemble_d_first=ensemble_d_first, ierr=ierr)
        call assert_true(is_err(ierr), "merged: should reject n_dimensions=1")
    end subroutine test_merged_n_dimensions_too_small

end module mod_test_shape_truthful_clustering_merged
