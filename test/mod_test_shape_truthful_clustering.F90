!> Unit test suite for tox_shape_truthful_clustering (ensemble_identification), generated from
!| src/tox/shape_truthful_clustering/tox_shape_truthful_clustering_impl.F90.
module mod_test_shape_truthful_clustering
    use tox_shape_truthful_clustering, only: ensemble_identification
    use tox_shape_truthful_clustering_impl, only: STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE, &
                                                    STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT, &
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
    function get_all_tests_shape_truthful_clustering() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(9))

        all_tests(1) = test_case("test_ensemble_identification_natural_fixed_point", &
                                 test_ensemble_identification_natural_fixed_point)
        all_tests(2) = test_case("test_ensemble_identification_history_window_shifts", &
                                 test_ensemble_identification_history_window_shifts)
        all_tests(3) = test_case("test_ensemble_identification_max_size_at_bootstrap", &
                                 test_ensemble_identification_max_size_at_bootstrap)
        all_tests(4) = test_case("test_ensemble_identification_max_size_poisons_prior_accepts", &
                                 test_ensemble_identification_max_size_poisons_prior_accepts)
        all_tests(5) = test_case("test_ensemble_identification_rejected_immediately", &
                                 test_ensemble_identification_rejected_immediately)
        all_tests(6) = test_case("test_ensemble_identification_rejected_after_stable", &
                                 test_ensemble_identification_rejected_after_stable)
        all_tests(7) = test_case("test_ensemble_identification_seed_index_out_of_range", &
                                 test_ensemble_identification_seed_index_out_of_range)
        all_tests(8) = test_case("test_ensemble_identification_o_zero", &
                                 test_ensemble_identification_o_zero)
        all_tests(9) = test_case("test_ensemble_identification_n_dimensions_too_small", &
                                 test_ensemble_identification_n_dimensions_too_small)
    end function get_all_tests_shape_truthful_clustering

    ! --- Fixture A: D=2, N=7. A 5-point line (0,0)..(4,0), plus two far-away points that a
    ! growth radius of 1.0 (k_min=1 from the seed at (0,0), whose nearest neighbor is at
    ! distance 1.0) never reaches. Seed=1. Growth: {1}->{1,2}->{1,2,3}->{1,2,3,4}->
    ! {1,2,3,4,5}->{1,2,3,4,5} (unchanged, natural fixed point at the 5th grow_ensemble call,
    ! i.e. growth iteration t=5; 4 growth iterations were accepted: the bootstrap and t=2..4).
    ! Purely collinear throughout, so d=1 and the tangent basis never rotates: every real
    ! accept_ensemble check along the way succeeds.

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

    subroutine test_ensemble_identification_natural_fixed_point()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: final_ensemble_mask(7), expected_mask(7)
        integer(int32) :: stop_reason, d_history(4), k_history(4), expected_k(4)
        integer(int32) :: member_added_at_step(7), expected_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(4)
        real(real64)   :: growth_radius, G_history(4), mu_history(2, 4), S_history(2, 4), U_history(2, 2, 4)
        real(real64)   :: U_first(2, 2)
        integer(int32) :: d_first

        call build_fixture_a(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 2_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=4_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(stop_reason, STOP_REASON_FIXED_POINT, "natural fixed point: stop_reason")
        call assert_equal_real(growth_radius, 1.0d0, 1.0d-9, "natural fixed point: growth_radius")

        expected_mask = .false.
        expected_mask(1:5) = .true.
        call assert_equal_array_logical(final_ensemble_mask, expected_mask, 7_int32, "natural fixed point: final_ensemble_mask")

        expected_k = [2, 3, 4, 5]
        call assert_equal_array_int(k_history, expected_k, 4_int32, "natural fixed point: k_history")
        call assert_true(all(accepted_history), "natural fixed point: every retained iteration was accepted")
        call assert_true(all(d_history == 1), "natural fixed point: d stays 1 throughout (collinear)")

        expected_step = [MEMBER_ADDED_AT_STEP_SEED, 1, 2, 3, 4, MEMBER_ADDED_AT_STEP_NON_MEMBER, &
                         MEMBER_ADDED_AT_STEP_NON_MEMBER]
        call assert_equal_array_int(member_added_at_step, expected_step, 7_int32, "natural fixed point: member_added_at_step")

        ! Iteration 1's own bootstrap mask -- {seed=1, its one growth-radius neighbor=2} --
        ! regardless of how much further growth went on to reach.
        expected_mask = .false.
        expected_mask(1:2) = .true.
        call assert_equal_array_logical(low_confidence_mask, expected_mask, 7_int32, &
                                        "natural fixed point: low_confidence_mask is iteration 1's own mask")

        ! U_first/d_first: the bootstrap iteration's own tangent basis, {1,2} along the x-axis
        ! -- collinear, so d_first=1 and the first column is the x-axis unit vector (up to
        ! SVD's own sign ambiguity), never overwritten by the three further accepted iterations.
        call assert_equal_int(d_first, 1_int32, "natural fixed point: d_first is the bootstrap's own d")
        call assert_equal_real(abs(U_first(1,1)), 1.0d0, 1.0d-9, "natural fixed point: U_first column 1, x")
        call assert_equal_real(abs(U_first(2,1)), 0.0d0, 1.0d-9, "natural fixed point: U_first column 1, y")
    end subroutine test_ensemble_identification_natural_fixed_point

    !> Same fixture and trajectory, but o=2: only the last 2 of the 4 accepted iterations
    !| (k=4 and k=5) survive the shift-and-append ring buffer.
    subroutine test_ensemble_identification_history_window_shifts()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: final_ensemble_mask(7)
        integer(int32) :: stop_reason, d_history(2), k_history(2), expected_k(2)
        integer(int32) :: member_added_at_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(2)
        real(real64)   :: growth_radius, G_history(2), mu_history(2, 2), S_history(2, 2), U_history(2, 2, 2)
        real(real64)   :: U_first(2, 2)
        integer(int32) :: d_first

        call build_fixture_a(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 2_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=2_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(stop_reason, STOP_REASON_FIXED_POINT, "history window: stop_reason")
        expected_k = [4, 5]
        call assert_equal_array_int(k_history, expected_k, 2_int32, "history window: only the last 2 iterations survive")
        call assert_true(all(accepted_history), "history window: both retained iterations were accepted")
    end subroutine test_ensemble_identification_history_window_shifts

    !> f_max=0.2 (threshold 1.4 of N=7): the bootstrap's own 2-member candidate already meets
    !| it, so Stop Condition 1 fires before `observable`/`accept_ensemble` are ever called.
    subroutine test_ensemble_identification_max_size_at_bootstrap()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: final_ensemble_mask(7)
        integer(int32) :: stop_reason, d_history(3), k_history(3), expected_k(3)
        integer(int32) :: member_added_at_step(7), expected_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(3)
        real(real64)   :: growth_radius, G_history(3), mu_history(2, 3), S_history(2, 3), U_history(2, 2, 3)
        real(real64)   :: U_first(2, 2)
        integer(int32) :: d_first

        call build_fixture_a(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 2_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, f_max=0.2d0, &
                                     o=3_int32, final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(stop_reason, STOP_REASON_MAX_SIZE, "max size at bootstrap: stop_reason")
        call assert_true(.not. any(final_ensemble_mask), "max size at bootstrap: no ensemble is returned")
        expected_k = [0, 0, 0]
        call assert_equal_array_int(k_history, expected_k, 3_int32, "max size at bootstrap: k_history is empty")
        call assert_true(.not. any(accepted_history), "max size at bootstrap: accepted_history is empty")
        expected_step = MEMBER_ADDED_AT_STEP_NON_MEMBER
        call assert_equal_array_int(member_added_at_step, expected_step, 7_int32, &
                                    "max size at bootstrap: member_added_at_step is empty")

        ! Stop Condition 1 fires before observable is ever called for this seed at all (see
        ! the module comment above): no genuine SVD ever happened, so unlike the poisoned-
        ! after-acceptance case below, there is no real iteration-1 data to report either.
        call assert_true(.not. any(low_confidence_mask), &
                         "max size at bootstrap: low_confidence_mask stays empty, no SVD ever ran")
        call assert_equal_int(d_first, 0_int32, "max size at bootstrap: d_first stays zero, no SVD ever ran")
        call assert_true(.not. any(U_first /= 0.0d0), "max size at bootstrap: U_first stays zero, no SVD ever ran")
    end subroutine test_ensemble_identification_max_size_at_bootstrap

    !> f_max=0.35 (threshold 2.45 of N=7): the bootstrap's 2-member candidate is under it and
    !| gets accepted, but growth iteration t=2's 3-member candidate meets it -- Stop Condition
    !| 1 then discards even the already-accepted bootstrap: no ensemble at all is returned.
    subroutine test_ensemble_identification_max_size_poisons_prior_accepts()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: final_ensemble_mask(7), expected_mask(7)
        integer(int32) :: stop_reason, d_history(3), k_history(3), expected_k(3)
        integer(int32) :: member_added_at_step(7), expected_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(3)
        real(real64)   :: growth_radius, G_history(3), mu_history(2, 3), S_history(2, 3), U_history(2, 2, 3)
        real(real64)   :: U_first(2, 2)
        integer(int32) :: d_first

        call build_fixture_a(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 2_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, f_max=0.35d0, &
                                     o=3_int32, final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(stop_reason, STOP_REASON_MAX_SIZE, "max size poisons: stop_reason")
        call assert_true(.not. any(final_ensemble_mask), &
                         "max size poisons: the bootstrap's earlier acceptance is discarded too")
        expected_k = [0, 0, 0]
        call assert_equal_array_int(k_history, expected_k, 3_int32, "max size poisons: k_history is wiped")
        expected_step = MEMBER_ADDED_AT_STEP_NON_MEMBER
        call assert_equal_array_int(member_added_at_step, expected_step, 7_int32, &
                                    "max size poisons: member_added_at_step is wiped")

        ! The key behavior this output exists for: Stop Condition 1 wipes final_ensemble_mask
        ! (asserted above) and every history array, but iteration 1's own bootstrap mask --
        ! {seed=1, its one growth-radius neighbor=2} -- was genuinely SVD-backed and is
        ! deliberately never touched by that reset, so it must still be there.
        expected_mask = .false.
        expected_mask(1:2) = .true.
        call assert_equal_array_logical(low_confidence_mask, expected_mask, 7_int32, &
                                        "max size poisons: low_confidence_mask survives the reset that wipes everything else")

        ! U_first/d_first get the exact same treatment as low_confidence_mask, and for the
        ! same reason -- see the module comment on U_first in the kernel.
        call assert_equal_int(d_first, 1_int32, "max size poisons: d_first survives the reset")
        call assert_equal_real(abs(U_first(1,1)), 1.0d0, 1.0d-9, "max size poisons: U_first survives the reset, x")
        call assert_equal_real(abs(U_first(2,1)), 0.0d0, 1.0d-9, "max size poisons: U_first survives the reset, y")
    end subroutine test_ensemble_identification_max_size_poisons_prior_accepts

    ! --- Fixture B: D=3, N=7. A 5-point x-axis line (0,0,0)..(4,0,0), plus a 2-point branch
    ! at (1,1,0),(1,2,0) right next to the 2nd x-axis point. Growth radius 1.0 (k_min=1 from
    ! the seed's nearest neighbor at distance 1.0) already sweeps the branch's first point in
    ! at growth iteration t=2 (grown from {1,2}), the very first iteration `accept_ensemble`
    ! is ever invoked for: a 5-point, still-planar (z=0) set with one off-axis point forces
    ! d=2 (see `observable`'s own z=0-plane tests), a jump from the bootstrap's d=1 that a
    ! strict d_max=0 rejects immediately -- before the ensemble was ever "stably" accepted.

    subroutine build_fixture_b(vectors, kd_indices, dim_order)
        real(real64), intent(out) :: vectors(3, 7)
        integer(int32), intent(out) :: kd_indices(7)
        integer(int32), intent(out) :: dim_order(3)
        integer(int32) :: i, ierr

        do i = 1, 5
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
            vectors(3, i) = 0.0d0
        end do
        vectors(:, 6) = [1.0d0, 1.0d0, 0.0d0]
        vectors(:, 7) = [1.0d0, 2.0d0, 0.0d0]
        dim_order = [1, 2, 3]

        call build_kd_index(vectors, 3_int32, 7_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'build_fixture_b: build_kd_index failed: ', ierr
            error stop
        end if
    end subroutine build_fixture_b

    subroutine test_ensemble_identification_rejected_immediately()
        real(real64)   :: vectors(3, 7)
        integer(int32) :: kd_indices(7), dim_order(3), ierr
        logical(c_bool)        :: final_ensemble_mask(7), expected_mask(7)
        integer(int32) :: stop_reason, d_history(2), k_history(2), expected_k(2)
        integer(int32) :: member_added_at_step(7), expected_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(2), expected_accepted(2)
        real(real64)   :: growth_radius, G_history(2), mu_history(3, 2), S_history(3, 2), U_history(3, 3, 2)
        real(real64)   :: U_first(3, 3)
        integer(int32) :: d_first

        call build_fixture_b(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 3_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=2_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(stop_reason, STOP_REASON_REJECTED_IMMEDIATELY, "rejected immediately: stop_reason")

        expected_mask = .false.
        expected_mask(1:2) = .true.
        call assert_equal_array_logical(final_ensemble_mask, expected_mask, 7_int32, &
                                        "rejected immediately: final_ensemble_mask is the bootstrap only")

        expected_k = [2, 4]
        call assert_equal_array_int(k_history, expected_k, 2_int32, &
                                    "rejected immediately: history keeps the bootstrap and the rejected candidate")
        expected_accepted = [.true., .false.]
        call assert_equal_array_logical(accepted_history, expected_accepted, 2_int32, &
                                        "rejected immediately: only the last column is unaccepted")

        expected_step = [MEMBER_ADDED_AT_STEP_SEED, 1, MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER, &
                         MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER]
        call assert_equal_array_int(member_added_at_step, expected_step, 7_int32, &
                                    "rejected immediately: the rejected candidate's members never joined")
    end subroutine test_ensemble_identification_rejected_immediately

    ! --- Fixture C: D=3, N=7. Same idea as Fixture B, but the branch (2,1,0),(2,2,0) sits
    ! next to the 3rd x-axis point instead of the 2nd -- so growth iteration t=2 (grown from
    ! the bootstrap's {1,2}, reaching only as far as {1,2,3}) stays purely collinear and is
    ! accepted, and the branch is only swept in at t=3 (grown from {1,2,3}), by which point
    ! the ensemble has already been accepted `a`=2 times.

    subroutine build_fixture_c(vectors, kd_indices, dim_order)
        real(real64), intent(out) :: vectors(3, 7)
        integer(int32), intent(out) :: kd_indices(7)
        integer(int32), intent(out) :: dim_order(3)
        integer(int32) :: i, ierr

        do i = 1, 5
            vectors(1, i) = real(i - 1, real64)
            vectors(2, i) = 0.0d0
            vectors(3, i) = 0.0d0
        end do
        vectors(:, 6) = [2.0d0, 1.0d0, 0.0d0]
        vectors(:, 7) = [2.0d0, 2.0d0, 0.0d0]
        dim_order = [1, 2, 3]

        call build_kd_index(vectors, 3_int32, 7_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'build_fixture_c: build_kd_index failed: ', ierr
            error stop
        end if
    end subroutine build_fixture_c

    subroutine test_ensemble_identification_rejected_after_stable()
        real(real64)   :: vectors(3, 7)
        integer(int32) :: kd_indices(7), dim_order(3), ierr
        logical(c_bool)        :: final_ensemble_mask(7), expected_mask(7)
        integer(int32) :: stop_reason, d_history(3), k_history(3), expected_k(3)
        integer(int32) :: member_added_at_step(7), expected_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(3), expected_accepted(3)
        real(real64)   :: growth_radius, G_history(3), mu_history(3, 3), S_history(3, 3), U_history(3, 3, 3)
        real(real64)   :: U_first(3, 3)
        integer(int32) :: d_first

        call build_fixture_c(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 3_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=3_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_identification failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(stop_reason, STOP_REASON_REJECTED_AFTER_STABLE, "rejected after stable: stop_reason")

        expected_mask = .false.
        expected_mask(1:3) = .true.
        call assert_equal_array_logical(final_ensemble_mask, expected_mask, 7_int32, &
                                        "rejected after stable: final_ensemble_mask is the last accepted state")

        expected_k = [2, 3, 5]
        call assert_equal_array_int(k_history, expected_k, 3_int32, &
                                    "rejected after stable: history keeps both accepted steps and the rejection")
        expected_accepted = [.true., .true., .false.]
        call assert_equal_array_logical(accepted_history, expected_accepted, 3_int32, &
                                        "rejected after stable: only the last column is unaccepted")

        expected_step = [MEMBER_ADDED_AT_STEP_SEED, 1, 2, MEMBER_ADDED_AT_STEP_NON_MEMBER, &
                         MEMBER_ADDED_AT_STEP_NON_MEMBER, MEMBER_ADDED_AT_STEP_NON_MEMBER, &
                         MEMBER_ADDED_AT_STEP_NON_MEMBER]
        call assert_equal_array_int(member_added_at_step, expected_step, 7_int32, &
                                    "rejected after stable: the rejected candidate's new members never joined")
    end subroutine test_ensemble_identification_rejected_after_stable

    ! --- Validation (via the generated wrapper's ierr) ---

    subroutine test_ensemble_identification_seed_index_out_of_range()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: final_ensemble_mask(7)
        integer(int32) :: stop_reason, d_history(3), k_history(3)
        integer(int32) :: member_added_at_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(3)
        real(real64)   :: growth_radius, G_history(3), mu_history(2, 3), S_history(2, 3), U_history(2, 2, 3)
        real(real64)   :: U_first(2, 2)
        integer(int32) :: d_first

        call build_fixture_a(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 2_int32, 7_int32, kd_indices, dim_order, 8_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=3_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_identification should reject seed_index > n_vectors")
    end subroutine test_ensemble_identification_seed_index_out_of_range

    subroutine test_ensemble_identification_o_zero()
        real(real64)   :: vectors(2, 7)
        integer(int32) :: kd_indices(7), dim_order(2), ierr
        logical(c_bool)        :: final_ensemble_mask(7)
        integer(int32) :: stop_reason, d_history(0), k_history(0)
        integer(int32) :: member_added_at_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(0)
        real(real64)   :: growth_radius, G_history(0), mu_history(2, 0), S_history(2, 0), U_history(2, 2, 0)
        real(real64)   :: U_first(2, 2)
        integer(int32) :: d_first

        call build_fixture_a(vectors, kd_indices, dim_order)

        call ensemble_identification(vectors, 2_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=0_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_identification should reject o=0")
    end subroutine test_ensemble_identification_o_zero

    subroutine test_ensemble_identification_n_dimensions_too_small()
        real(real64)   :: vectors(1, 7)
        integer(int32) :: kd_indices(7), dim_order(1), ierr
        logical(c_bool)        :: final_ensemble_mask(7)
        integer(int32) :: stop_reason, d_history(3), k_history(3)
        integer(int32) :: member_added_at_step(7)
        logical(c_bool)        :: low_confidence_mask(7)
        logical(c_bool)        :: accepted_history(3)
        real(real64)   :: growth_radius, G_history(3), mu_history(1, 3), S_history(1, 3), U_history(1, 1, 3)
        real(real64)   :: U_first(1, 1)
        integer(int32) :: d_first
        integer(int32) :: i

        do i = 1, 7
            vectors(1, i) = real(i - 1, real64)
        end do
        dim_order = [1]
        call build_kd_index(vectors, 1_int32, 7_int32, kd_indices, dim_order, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'test_ensemble_identification_n_dimensions_too_small: build_kd_index failed: ', ierr
            error stop
        end if

        call ensemble_identification(vectors, 1_int32, 7_int32, kd_indices, dim_order, 1_int32, &
                                     k_min=1_int32, chordal_dist_max_as_prcnt_of_range=0.1d0, d_max=0_int32, G_max=1.0d10, RMSE_change_max=1.0d10, o=3_int32, &
                                     final_ensemble_mask=final_ensemble_mask, stop_reason=stop_reason, &
                                     growth_radius=growth_radius, U_history=U_history, S_history=S_history, &
                                     d_history=d_history, G_history=G_history, mu_history=mu_history, &
                                     k_history=k_history, accepted_history=accepted_history, &
                                     member_added_at_step=member_added_at_step, &
                                     low_confidence_mask=low_confidence_mask, U_first=U_first, d_first=d_first, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_identification should reject n_dimensions=1")
    end subroutine test_ensemble_identification_n_dimensions_too_small

end module mod_test_shape_truthful_clustering
