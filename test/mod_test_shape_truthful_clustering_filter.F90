!> Unit test suite for tox_shape_truthful_clustering_filter (filter_ensembles_by_stop_condition,
!| filter_ensembles_by_dimension, filter_ensembles_by_var_explained, filter_ensembles), generated
!| from src/tox/shape_truthful_clustering/tox_shape_truthful_clustering_filter_impl.F90.
module mod_test_shape_truthful_clustering_filter
    use tox_shape_truthful_clustering_filter, only: filter_ensembles_by_stop_condition, &
        filter_ensembles_by_dimension, filter_ensembles_by_var_explained, filter_ensembles
    use tox_shape_truthful_clustering_impl, only: STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE, &
        STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT
    use tox_errors, only: is_ok, is_err
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_filter() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(19))

        all_tests(1) = test_case("test_filter_stop_condition_exact_exclusion", &
                                 test_filter_stop_condition_exact_exclusion)
        all_tests(2) = test_case("test_filter_stop_condition_absent_is_noop", &
                                 test_filter_stop_condition_absent_is_noop)
        all_tests(3) = test_case("test_filter_stop_condition_all_four_values", &
                                 test_filter_stop_condition_all_four_values)
        all_tests(4) = test_case("test_filter_dimension_d_min_only", test_filter_dimension_d_min_only)
        all_tests(5) = test_case("test_filter_dimension_d_max_only", test_filter_dimension_d_max_only)
        all_tests(6) = test_case("test_filter_dimension_both_bounds", test_filter_dimension_both_bounds)
        all_tests(7) = test_case("test_filter_dimension_both_absent_is_noop", &
                                 test_filter_dimension_both_absent_is_noop)
        all_tests(8) = test_case("test_filter_dimension_no_final_excluded_once_bound_present", &
                                 test_filter_dimension_no_final_excluded_once_bound_present)
        all_tests(9) = test_case("test_filter_var_explained_clean_fixture", &
                                 test_filter_var_explained_clean_fixture)
        all_tests(10) = test_case("test_filter_var_explained_threshold_at_boundary", &
                                  test_filter_var_explained_threshold_at_boundary)
        all_tests(11) = test_case("test_filter_var_explained_k_le_one_guard", &
                                  test_filter_var_explained_k_le_one_guard)
        all_tests(12) = test_case("test_filter_var_explained_absent_is_noop", &
                                  test_filter_var_explained_absent_is_noop)
        all_tests(13) = test_case("test_filter_ensembles_combined_different_criteria", &
                                  test_filter_ensembles_combined_different_criteria)
        all_tests(14) = test_case("test_filter_ensembles_all_omitted_is_noop", &
                                  test_filter_ensembles_all_omitted_is_noop)
        all_tests(15) = test_case("test_filter_stop_condition_zero_ensembles", &
                                  test_filter_stop_condition_zero_ensembles)
        all_tests(16) = test_case("test_filter_dimension_invalid_n_dimensions", &
                                  test_filter_dimension_invalid_n_dimensions)
        all_tests(17) = test_case("test_filter_var_explained_invalid_threshold", &
                                  test_filter_var_explained_invalid_threshold)
        all_tests(18) = test_case("test_filter_stop_condition_invalid_value", &
                                  test_filter_stop_condition_invalid_value)
        all_tests(19) = test_case("test_filter_dimension_d_min_exceeds_d_max_still_computes", &
                                  test_filter_dimension_d_min_exceeds_d_max_still_computes)
    end function get_all_tests_shape_truthful_clustering_filter

    ! --- filter_ensembles_by_stop_condition ---------------------------------

    !> 4 ensembles, one of each Stop Condition. Excluding exactly
    !| STOP_REASON_REJECTED_IMMEDIATELY must leave that one ensemble (and only that one)
    !| ineligible.
    subroutine test_filter_stop_condition_exact_exclusion()
        integer(int32) :: ensemble_stop_reason(4)
        logical(c_bool)        :: allowed(4), eligible(4)
        integer(int32) :: ierr

        ensemble_stop_reason = [STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE, &
                                STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT]
        allowed = [.true., .true., .false., .true.]

        call filter_ensembles_by_stop_condition(4_int32, ensemble_stop_reason, allowed, eligible, ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_stop_condition failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(eligible(1), "stop condition: max_size ensemble still eligible")
        call assert_true(eligible(2), "stop condition: rejected_after_stable ensemble still eligible")
        call assert_true(.not. eligible(3), "stop condition: rejected_immediately ensemble excluded")
        call assert_true(eligible(4), "stop condition: fixed_point ensemble still eligible")
    end subroutine test_filter_stop_condition_exact_exclusion

    !> Omitting `allowed_stop_reasons` entirely must be a true no-op -- every ensemble eligible
    !| regardless of its actual Stop Condition.
    subroutine test_filter_stop_condition_absent_is_noop()
        integer(int32) :: ensemble_stop_reason(4)
        logical(c_bool)        :: eligible(4)
        integer(int32) :: ierr

        ensemble_stop_reason = [STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE, &
                                STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT]

        call filter_ensembles_by_stop_condition(4_int32, ensemble_stop_reason, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_stop_condition failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(all(eligible), "stop condition absent: every ensemble eligible")
    end subroutine test_filter_stop_condition_absent_is_noop

    !> All four STOP_REASON_* values exercised at once, each disallowed one at a time across 4
    !| ensembles sharing the same Stop Condition (STOP_REASON_FIXED_POINT) -- proves the mapping
    !| from `ensemble_stop_reason` index to `allowed_stop_reasons` slot is exactly 1:1, not
    !| off-by-one in either direction.
    subroutine test_filter_stop_condition_all_four_values()
        integer(int32) :: ensemble_stop_reason(1)
        logical(c_bool)        :: allowed(4), eligible(1)
        integer(int32) :: ierr, r

        do r = 1, 4
            ensemble_stop_reason(1) = r
            allowed = .true.
            allowed(r) = .false.

            call filter_ensembles_by_stop_condition(1_int32, ensemble_stop_reason, allowed, eligible, ierr)
            if (.not. is_ok(ierr)) then
                write (*, *) 'filter_ensembles_by_stop_condition failed unexpectedly: ', ierr
                error stop
            end if

            call assert_true(.not. eligible(1), "stop condition all values: disallowing its own value excludes it")
        end do
    end subroutine test_filter_stop_condition_all_four_values

    subroutine test_filter_stop_condition_zero_ensembles()
        integer(int32) :: ensemble_stop_reason(0)
        logical(c_bool)        :: eligible(0)
        integer(int32) :: ierr

        call filter_ensembles_by_stop_condition(0_int32, ensemble_stop_reason, eligible=eligible, ierr=ierr)
        call assert_true(is_ok(ierr), "stop condition: n_ensembles=0 must not be an error")
    end subroutine test_filter_stop_condition_zero_ensembles

    subroutine test_filter_stop_condition_invalid_value()
        integer(int32) :: ensemble_stop_reason(1)
        logical(c_bool)        :: eligible(1)
        integer(int32) :: ierr

        ensemble_stop_reason(1) = 5_int32 ! only 1..4 are valid Stop Condition indices

        call filter_ensembles_by_stop_condition(1_int32, ensemble_stop_reason, eligible=eligible, ierr=ierr)
        call assert_true(is_err(ierr), "stop condition: an out-of-range ensemble_stop_reason value must be rejected")
    end subroutine test_filter_stop_condition_invalid_value

    ! --- filter_ensembles_by_dimension --------------------------------------

    !> D=3. d_final = [0,1,2,3]. d_min=2 alone must exclude the first two (0,1), keep the last two.
    subroutine test_filter_dimension_d_min_only()
        integer(int32) :: ensemble_d_final(4)
        logical(c_bool)        :: ensemble_has_final(4), eligible(4)
        integer(int32) :: ierr

        ensemble_d_final = [0, 1, 2, 3]
        ensemble_has_final = .true.

        call filter_ensembles_by_dimension(3_int32, 4_int32, ensemble_d_final, ensemble_has_final, &
                                           d_min=2_int32, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_dimension failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(.not. eligible(1), "d_min only: d=0 excluded")
        call assert_true(.not. eligible(2), "d_min only: d=1 excluded")
        call assert_true(eligible(3), "d_min only: d=2 included (boundary, inclusive)")
        call assert_true(eligible(4), "d_min only: d=3 included")
    end subroutine test_filter_dimension_d_min_only

    !> D=3. d_final = [0,1,2,3]. d_max=1 alone must keep the first two (0,1), exclude the last two.
    subroutine test_filter_dimension_d_max_only()
        integer(int32) :: ensemble_d_final(4)
        logical(c_bool)        :: ensemble_has_final(4), eligible(4)
        integer(int32) :: ierr

        ensemble_d_final = [0, 1, 2, 3]
        ensemble_has_final = .true.

        call filter_ensembles_by_dimension(3_int32, 4_int32, ensemble_d_final, ensemble_has_final, &
                                           d_max=1_int32, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_dimension failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(eligible(1), "d_max only: d=0 included")
        call assert_true(eligible(2), "d_max only: d=1 included (boundary, inclusive)")
        call assert_true(.not. eligible(3), "d_max only: d=2 excluded")
        call assert_true(.not. eligible(4), "d_max only: d=3 excluded")
    end subroutine test_filter_dimension_d_max_only

    !> D=3. d_final = [0,1,2,3]. d_min=1, d_max=2 together must keep only the middle two.
    subroutine test_filter_dimension_both_bounds()
        integer(int32) :: ensemble_d_final(4)
        logical(c_bool)        :: ensemble_has_final(4), eligible(4)
        integer(int32) :: ierr

        ensemble_d_final = [0, 1, 2, 3]
        ensemble_has_final = .true.

        call filter_ensembles_by_dimension(3_int32, 4_int32, ensemble_d_final, ensemble_has_final, &
                                           d_min=1_int32, d_max=2_int32, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_dimension failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(.not. eligible(1), "both bounds: d=0 excluded (below d_min)")
        call assert_true(eligible(2), "both bounds: d=1 included")
        call assert_true(eligible(3), "both bounds: d=2 included")
        call assert_true(.not. eligible(4), "both bounds: d=3 excluded (above d_max)")
    end subroutine test_filter_dimension_both_bounds

    !> Both bounds absent must be a true no-op -- every ensemble eligible, even one with
    !| `ensemble_has_final=.false.` (the kernel's own doc comment: `ensemble_has_final` is "not
    !| even consulted" in this case).
    subroutine test_filter_dimension_both_absent_is_noop()
        integer(int32) :: ensemble_d_final(2)
        logical(c_bool)        :: ensemble_has_final(2), eligible(2)
        integer(int32) :: ierr

        ensemble_d_final = [0, 3]
        ensemble_has_final = [.true., .false.]

        call filter_ensembles_by_dimension(3_int32, 2_int32, ensemble_d_final, ensemble_has_final, eligible=eligible, &
                                           ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_dimension failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(all(eligible), "both bounds absent: every ensemble eligible, even has_final=.false.")
    end subroutine test_filter_dimension_both_absent_is_noop

    !> Once at least one bound is supplied, an ensemble with no final accepted state at all
    !| (`ensemble_has_final=.false.`) must be ineligible -- there is no `d` to judge, regardless
    !| of whatever `ensemble_d_final` happens to hold for it (here, deliberately a value that
    !| would otherwise pass).
    subroutine test_filter_dimension_no_final_excluded_once_bound_present()
        integer(int32) :: ensemble_d_final(1)
        logical(c_bool)        :: ensemble_has_final(1), eligible(1)
        integer(int32) :: ierr

        ensemble_d_final(1) = 1 ! would satisfy d_min=0/d_max=3 if has_final were true
        ensemble_has_final(1) = .false.

        call filter_ensembles_by_dimension(3_int32, 1_int32, ensemble_d_final, ensemble_has_final, &
                                           d_min=0_int32, d_max=3_int32, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_dimension failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(.not. eligible(1), "no final: ineligible once a bound is supplied, despite a passing d_final")
    end subroutine test_filter_dimension_no_final_excluded_once_bound_present

    subroutine test_filter_dimension_invalid_n_dimensions()
        integer(int32) :: ensemble_d_final(1)
        logical(c_bool)        :: ensemble_has_final(1), eligible(1)
        integer(int32) :: ierr

        ensemble_d_final(1) = 0
        ensemble_has_final(1) = .true.

        call filter_ensembles_by_dimension(1_int32, 1_int32, ensemble_d_final, ensemble_has_final, &
                                           d_min=0_int32, eligible=eligible, ierr=ierr)
        call assert_true(is_err(ierr), "dimension filter: n_dimensions=1 must be rejected (minimum is 2)")
    end subroutine test_filter_dimension_invalid_n_dimensions

    !> d_min > d_max is not itself validated (each bound is independently checked against
    !| [0, n_dimensions], not against each other) -- the combination is simply, and correctly,
    !| unsatisfiable by construction: every ensemble ends up ineligible.
    subroutine test_filter_dimension_d_min_exceeds_d_max_still_computes()
        integer(int32) :: ensemble_d_final(3)
        logical(c_bool)        :: ensemble_has_final(3), eligible(3)
        integer(int32) :: ierr

        ensemble_d_final = [0, 1, 2]
        ensemble_has_final = .true.

        call filter_ensembles_by_dimension(2_int32, 3_int32, ensemble_d_final, ensemble_has_final, &
                                           d_min=2_int32, d_max=1_int32, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_dimension failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(.not. any(eligible), "d_min > d_max: unsatisfiable, every ensemble ineligible")
    end subroutine test_filter_dimension_d_min_exceeds_d_max_still_computes

    ! --- filter_ensembles_by_var_explained ----------------------------------

    !> D=2, k=2 (k-1=1, so eigenvalues = S**2 exactly), d=1. S=[10,1] -> eigenvalues [100,1] ->
    !| variance explained = 100/(100+1) = 100/101, a clean, hand-computable fraction.
    !| var_explained_min=0.9 (< 100/101 ~ 0.9901) must pass.
    subroutine test_filter_var_explained_clean_fixture()
        real(real64)   :: ensemble_S_final(2, 1)
        integer(int32) :: ensemble_d_final(1), ensemble_k_final(1)
        logical(c_bool)        :: ensemble_has_final(1), eligible(1)
        integer(int32) :: ierr

        ensemble_S_final(:, 1) = [10.0d0, 1.0d0]
        ensemble_d_final(1) = 1
        ensemble_k_final(1) = 2
        ensemble_has_final(1) = .true.

        call filter_ensembles_by_var_explained(2_int32, 1_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, var_explained_min=0.9d0, &
                                               eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_var_explained failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(eligible(1), "var explained clean fixture: 100/101 >= 0.9")
    end subroutine test_filter_var_explained_clean_fixture

    !> Same clean 100/101 fixture as above, with three thresholds straddling it: strictly below
    !| (must pass), essentially at it (must pass, inclusive), and strictly above (must fail).
    subroutine test_filter_var_explained_threshold_at_boundary()
        real(real64)   :: ensemble_S_final(2, 1)
        integer(int32) :: ensemble_d_final(1), ensemble_k_final(1)
        logical(c_bool)        :: ensemble_has_final(1), eligible(1)
        integer(int32) :: ierr
        real(real64), parameter :: ve = 100.0d0/101.0d0

        ensemble_S_final(:, 1) = [10.0d0, 1.0d0]
        ensemble_d_final(1) = 1
        ensemble_k_final(1) = 2
        ensemble_has_final(1) = .true.

        call filter_ensembles_by_var_explained(2_int32, 1_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, var_explained_min=ve - 1.0d-9, &
                                               eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_var_explained failed unexpectedly: ', ierr
            error stop
        end if
        call assert_true(eligible(1), "var explained boundary: threshold just below ve must pass")

        call filter_ensembles_by_var_explained(2_int32, 1_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, var_explained_min=ve, &
                                               eligible=eligible, ierr=ierr)
        call assert_true(eligible(1), "var explained boundary: threshold exactly at ve must pass (inclusive)")

        call filter_ensembles_by_var_explained(2_int32, 1_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, var_explained_min=ve + 1.0d-9, &
                                               eligible=eligible, ierr=ierr)
        call assert_true(.not. eligible(1), "var explained boundary: threshold just above ve must fail")
    end subroutine test_filter_var_explained_threshold_at_boundary

    !> k_final<=1 leaves no meaningful k-1 denominator -- must be ineligible once
    !| var_explained_min is supplied, regardless of S_final/d_final/has_final.
    subroutine test_filter_var_explained_k_le_one_guard()
        real(real64)   :: ensemble_S_final(2, 2)
        integer(int32) :: ensemble_d_final(2), ensemble_k_final(2)
        logical(c_bool)        :: ensemble_has_final(2), eligible(2)
        integer(int32) :: ierr

        ensemble_S_final(:, 1) = [10.0d0, 1.0d0]
        ensemble_d_final(1) = 1
        ensemble_k_final(1) = 1 ! degenerate: k-1 = 0
        ensemble_has_final(1) = .true.

        ensemble_S_final(:, 2) = [10.0d0, 1.0d0]
        ensemble_d_final(2) = 1
        ensemble_k_final(2) = 0
        ensemble_has_final(2) = .true.

        call filter_ensembles_by_var_explained(2_int32, 2_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, var_explained_min=0.0d0, &
                                               eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_var_explained failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(.not. eligible(1), "k<=1 guard: k=1 ineligible even with var_explained_min=0.0")
        call assert_true(.not. eligible(2), "k<=1 guard: k=0 ineligible even with var_explained_min=0.0")
    end subroutine test_filter_var_explained_k_le_one_guard

    !> Omitting `var_explained_min` entirely must be a true no-op -- every ensemble eligible,
    !| even a degenerate k<=1 one that would otherwise be excluded.
    subroutine test_filter_var_explained_absent_is_noop()
        real(real64)   :: ensemble_S_final(2, 1)
        integer(int32) :: ensemble_d_final(1), ensemble_k_final(1)
        logical(c_bool)        :: ensemble_has_final(1), eligible(1)
        integer(int32) :: ierr

        ensemble_S_final(:, 1) = [10.0d0, 1.0d0]
        ensemble_d_final(1) = 1
        ensemble_k_final(1) = 1 ! degenerate k, irrelevant when the filter is absent
        ensemble_has_final(1) = .true.

        call filter_ensembles_by_var_explained(2_int32, 1_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, eligible=eligible, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles_by_var_explained failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(eligible(1), "var explained absent: no-op, eligible even for a degenerate k")
    end subroutine test_filter_var_explained_absent_is_noop

    subroutine test_filter_var_explained_invalid_threshold()
        real(real64)   :: ensemble_S_final(2, 1)
        integer(int32) :: ensemble_d_final(1), ensemble_k_final(1)
        logical(c_bool)        :: ensemble_has_final(1), eligible(1)
        integer(int32) :: ierr

        ensemble_S_final(:, 1) = [10.0d0, 1.0d0]
        ensemble_d_final(1) = 1
        ensemble_k_final(1) = 2
        ensemble_has_final(1) = .true.

        call filter_ensembles_by_var_explained(2_int32, 1_int32, ensemble_S_final, ensemble_d_final, &
                                               ensemble_k_final, ensemble_has_final, var_explained_min=1.5d0, &
                                               eligible=eligible, ierr=ierr)
        call assert_true(is_err(ierr), "var explained: var_explained_min > 1.0 must be rejected")
    end subroutine test_filter_var_explained_invalid_threshold

    ! --- filter_ensembles (combined orchestrator) ---------------------------

    !> D=2, o=1, 4 ensembles, each failing a *different* single criterion (or none):
    !| E1: fails stop condition (rejected_immediately, disallowed).
    !| E2: fails dimension (d=2, d_max=1).
    !| E3: fails variance explained (S=[1,10] -> ve=1/101 < 0.5).
    !| E4: passes every criterion.
    !| All four masks must be asserted independently, plus the combined `eligible`.
    subroutine test_filter_ensembles_combined_different_criteria()
        real(real64)   :: ensemble_U_history(2, 2, 1, 4), ensemble_S_history(2, 1, 4), ensemble_mu_history(2, 1, 4)
        real(real64)   :: ensemble_G_history(1, 4)
        integer(int32) :: ensemble_d_history(1, 4), ensemble_k_history(1, 4), ensemble_stop_reason(4)
        logical(c_bool)        :: ensemble_accepted_history(1, 4)
        logical(c_bool)        :: allowed(4)
        logical(c_bool)        :: eligible(4), eligible_by_stop_condition(4), eligible_by_dimension(4)
        logical(c_bool)        :: eligible_by_var_explained(4)
        integer(int32) :: ierr

        ensemble_U_history = 0.0d0
        ensemble_mu_history = 0.0d0
        ensemble_G_history = 0.0d0
        ensemble_k_history = 2 ! k-1=1 throughout
        ensemble_accepted_history = .true.
        ensemble_d_history(1, :) = 1
        ensemble_S_history(:, 1, 1) = [10.0d0, 1.0d0]
        ensemble_S_history(:, 1, 2) = [10.0d0, 1.0d0]
        ensemble_S_history(:, 1, 3) = [1.0d0, 10.0d0] ! ve = 1/101, fails var_explained_min=0.5
        ensemble_S_history(:, 1, 4) = [10.0d0, 1.0d0]

        ensemble_d_history(1, 2) = 2 ! fails d_max=1

        ensemble_stop_reason = STOP_REASON_FIXED_POINT
        ensemble_stop_reason(1) = STOP_REASON_REJECTED_IMMEDIATELY

        allowed = [.true., .true., .false., .true.]

        call filter_ensembles(2_int32, 1_int32, 4_int32, ensemble_U_history, ensemble_d_history, ensemble_S_history, &
                              ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, &
                              ensemble_stop_reason, allowed_stop_reasons=allowed, d_max=1_int32, &
                              var_explained_min=0.5d0, eligible=eligible, &
                              eligible_by_stop_condition=eligible_by_stop_condition, &
                              eligible_by_dimension=eligible_by_dimension, &
                              eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles failed unexpectedly: ', ierr
            error stop
        end if

        ! Ensemble 1: fails stop condition only.
        call assert_true(.not. eligible_by_stop_condition(1), "combined: ensemble 1 fails stop condition")
        call assert_true(eligible_by_dimension(1), "combined: ensemble 1 passes dimension")
        call assert_true(eligible_by_var_explained(1), "combined: ensemble 1 passes var explained")
        call assert_true(.not. eligible(1), "combined: ensemble 1 overall ineligible")

        ! Ensemble 2: fails dimension only.
        call assert_true(eligible_by_stop_condition(2), "combined: ensemble 2 passes stop condition")
        call assert_true(.not. eligible_by_dimension(2), "combined: ensemble 2 fails dimension")
        call assert_true(eligible_by_var_explained(2), "combined: ensemble 2 passes var explained")
        call assert_true(.not. eligible(2), "combined: ensemble 2 overall ineligible")

        ! Ensemble 3: fails variance explained only.
        call assert_true(eligible_by_stop_condition(3), "combined: ensemble 3 passes stop condition")
        call assert_true(eligible_by_dimension(3), "combined: ensemble 3 passes dimension")
        call assert_true(.not. eligible_by_var_explained(3), "combined: ensemble 3 fails var explained")
        call assert_true(.not. eligible(3), "combined: ensemble 3 overall ineligible")

        ! Ensemble 4: passes everything.
        call assert_true(eligible_by_stop_condition(4), "combined: ensemble 4 passes stop condition")
        call assert_true(eligible_by_dimension(4), "combined: ensemble 4 passes dimension")
        call assert_true(eligible_by_var_explained(4), "combined: ensemble 4 passes var explained")
        call assert_true(eligible(4), "combined: ensemble 4 overall eligible")
    end subroutine test_filter_ensembles_combined_different_criteria

    !> Omitting every optional filter (`allowed_stop_reasons`/`d_min`/`d_max`/
    !| `var_explained_min`) must leave every ensemble eligible under all four masks, matching
    !| each individual filter's own no-op convention.
    subroutine test_filter_ensembles_all_omitted_is_noop()
        real(real64)   :: ensemble_U_history(2, 2, 1, 2), ensemble_S_history(2, 1, 2), ensemble_mu_history(2, 1, 2)
        real(real64)   :: ensemble_G_history(1, 2)
        integer(int32) :: ensemble_d_history(1, 2), ensemble_k_history(1, 2), ensemble_stop_reason(2)
        logical(c_bool)        :: ensemble_accepted_history(1, 2)
        logical(c_bool)        :: eligible(2), eligible_by_stop_condition(2), eligible_by_dimension(2)
        logical(c_bool)        :: eligible_by_var_explained(2)
        integer(int32) :: ierr

        ensemble_U_history = 0.0d0
        ensemble_S_history = 0.0d0
        ensemble_mu_history = 0.0d0
        ensemble_G_history = 0.0d0
        ensemble_d_history = 0
        ensemble_k_history = 0 ! even a totally unpopulated history must not matter here
        ensemble_accepted_history = .false.
        ensemble_stop_reason = STOP_REASON_FIXED_POINT

        call filter_ensembles(2_int32, 1_int32, 2_int32, ensemble_U_history, ensemble_d_history, ensemble_S_history, &
                              ensemble_mu_history, ensemble_G_history, ensemble_k_history, ensemble_accepted_history, &
                              ensemble_stop_reason, eligible=eligible, &
                              eligible_by_stop_condition=eligible_by_stop_condition, &
                              eligible_by_dimension=eligible_by_dimension, &
                              eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'filter_ensembles failed unexpectedly: ', ierr
            error stop
        end if

        call assert_true(all(eligible), "all omitted: combined eligible all true")
        call assert_true(all(eligible_by_stop_condition), "all omitted: stop condition mask all true")
        call assert_true(all(eligible_by_dimension), "all omitted: dimension mask all true")
        call assert_true(all(eligible_by_var_explained), "all omitted: var explained mask all true")
    end subroutine test_filter_ensembles_all_omitted_is_noop

end module mod_test_shape_truthful_clustering_filter
