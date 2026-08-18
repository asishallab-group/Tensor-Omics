!> Unit test suite for tox_shape_truthful_clustering_reconciliation (ensemble_reconciliation,
!| merge_to_super_ensembles), generated from
!| src/tox/shape_truthful_clustering/tox_shape_truthful_clustering_reconciliation_impl.F90.
module mod_test_shape_truthful_clustering_reconciliation
    use tox_shape_truthful_clustering_reconciliation, only: ensemble_reconciliation, merge_to_super_ensembles
    use tox_shape_truthful_clustering_reconciliation_impl, only: MODE_REPORT, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_MERGE_ANY
    use tox_shape_truthful_clustering_impl, only: STOP_REASON_MAX_SIZE, STOP_REASON_REJECTED_AFTER_STABLE, &
        STOP_REASON_REJECTED_IMMEDIATELY, STOP_REASON_FIXED_POINT
    use tox_errors, only: is_ok, is_err, ERR_SIZE_MISMATCH
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_reconciliation() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(16))

        all_tests(1) = test_case("test_reconciliation_report_mode_no_transitive_merge", &
                                 test_reconciliation_report_mode_no_transitive_merge)
        all_tests(2) = test_case("test_reconciliation_merge_any_transitive", &
                                 test_reconciliation_merge_any_transitive)
        all_tests(3) = test_case("test_reconciliation_merge_overlap_coefficient_threshold_excludes_weak_chain", &
                                 test_reconciliation_merge_oc_excludes_weak_chain)
        all_tests(4) = test_case("test_reconciliation_merge_overlap_coefficient_threshold_includes_all", &
                                 test_reconciliation_merge_oc_includes_all)
        all_tests(5) = test_case("test_reconciliation_overlap_coefficient_not_computed_unless_requested", &
                                 test_reconciliation_oc_not_computed_unless_requested)
        all_tests(6) = test_case("test_reconciliation_group_exceeds_max_group_size", &
                                 test_reconciliation_group_exceeds_max_group_size)
        all_tests(7) = test_case("test_reconciliation_n_ensembles_too_small", &
                                 test_reconciliation_n_ensembles_too_small)
        all_tests(8) = test_case("test_reconciliation_invalid_mode", test_reconciliation_invalid_mode)
        all_tests(9) = test_case("test_reconciliation_stop_reason_filter_excludes_pair", &
                                 test_reconciliation_stop_reason_filter_excludes_pair)
        all_tests(10) = test_case("test_reconciliation_stop_reason_filter_absent_is_noop", &
                                  test_reconciliation_stop_reason_filter_absent_is_noop)
        all_tests(11) = test_case("test_reconciliation_stop_reason_filter_breaks_chain", &
                                  test_reconciliation_stop_reason_filter_breaks_chain)
        all_tests(12) = test_case("test_reconciliation_stop_reason_filter_unused_value_noop", &
                                  test_reconciliation_stop_reason_filter_unused_value_noop)
        all_tests(13) = test_case("test_reconciliation_dimension_filter_excludes_pair", &
                                  test_reconciliation_dimension_filter_excludes_pair)
        all_tests(14) = test_case("test_reconciliation_var_explained_filter_excludes_pair", &
                                  test_reconciliation_var_explained_filter_excludes_pair)
        all_tests(15) = test_case("test_merge_to_super_ensembles_all_eligible_matches_merge_any", &
                                  test_merge_to_super_ensembles_all_eligible_matches_merge_any)
        all_tests(16) = test_case("test_merge_to_super_ensembles_excludes_ineligible_ensemble", &
                                  test_merge_to_super_ensembles_excludes_ineligible_ensemble)
    end function get_all_tests_shape_truthful_clustering_reconciliation

    !> N=14 vectors, 6 ensembles.
    !| E1={1,2,3,4}, E2={3,4,5,6}, E3={5,6,7,8}: a chain, E1-E2 and E2-E3 each intersect at 2
    !| members, both ensembles size 4, so OC = 2/min(4,4) = 0.5. E1-E3 do not intersect at all.
    !| E4={9,10}: isolated, intersects nobody.
    !| E5={11,12,13}, E6={12,13,14}: a separate pair, intersecting at 2 members, both ensembles
    !| size 3, so OC = 2/min(3,3) = 2/3 -- deliberately the *stronger* of the two edge
    !| strengths (unlike this same fixture under the old Jaccard formula, where the 1-2-3
    !| chain's 1/3 and the 5-6 pair's 1/2 were already ordered the same way, this is not a
    !| coincidence carried over -- both formulas happen to preserve the same ranking here, but
    !| the actual values differ, see the individual test comments below).
    !|
    !| `ensemble_stop_reason` defaults every ensemble to STOP_REASON_FIXED_POINT, a neutral
    !| choice that (with no `allowed_stop_reasons` filter applied) changes nothing about any
    !| pre-existing test below; the `allowed_stop_reasons`-specific tests override individual
    !| entries via the optional `stop_reason_overrides` argument.
    !|
    !| Also returns a minimal, uniform D=2/o=1 history fixture (every ensemble "accepted" at a
    !| single column, k=2, d=0, everything else zero) -- `ensemble_reconciliation`'s own new
    !| required history arguments, needed only to drive `filter_ensembles_impl`'s dimension/
    !| variance-explained filters, neither of which any test below (other than
    !| `test_reconciliation_dimension_filter_excludes_pair`/
    !| `test_reconciliation_var_explained_filter_excludes_pair`, which build their own dedicated
    !| history fixtures) actually exercises.
    subroutine build_fixture(ensemble_masks, ensemble_stop_reason, &
                             ensemble_U_history, ensemble_d_history, ensemble_S_history, ensemble_mu_history, &
                             ensemble_G_history, ensemble_k_history, ensemble_accepted_history, &
                             stop_reason_overrides)
        logical(c_bool), intent(out) :: ensemble_masks(14, 6)
        integer(int32), intent(out) :: ensemble_stop_reason(6)
        real(real64), intent(out) :: ensemble_U_history(2, 2, 1, 6)
        integer(int32), intent(out) :: ensemble_d_history(1, 6)
        real(real64), intent(out) :: ensemble_S_history(2, 1, 6)
        real(real64), intent(out) :: ensemble_mu_history(2, 1, 6)
        real(real64), intent(out) :: ensemble_G_history(1, 6)
        integer(int32), intent(out) :: ensemble_k_history(1, 6)
        logical(c_bool), intent(out) :: ensemble_accepted_history(1, 6)
        integer(int32), intent(in), optional :: stop_reason_overrides(6)
            !! Non-zero entries override the corresponding ensemble's default
            !! STOP_REASON_FIXED_POINT; zero entries leave the default in place.

        integer(int32) :: k

        ensemble_masks = .false.
        ensemble_masks(1:4, 1) = .true.
        ensemble_masks(3:6, 2) = .true.
        ensemble_masks(5:8, 3) = .true.
        ensemble_masks(9:10, 4) = .true.
        ensemble_masks(11:13, 5) = .true.
        ensemble_masks(12:14, 6) = .true.

        ensemble_stop_reason = STOP_REASON_FIXED_POINT
        if (present(stop_reason_overrides)) then
            do k = 1, 6
                if (stop_reason_overrides(k) /= 0) ensemble_stop_reason(k) = stop_reason_overrides(k)
            end do
        end if

        ensemble_U_history = 0.0d0
        ensemble_d_history = 0
        ensemble_S_history = 0.0d0
        ensemble_mu_history = 0.0d0
        ensemble_G_history = 0.0d0
        ensemble_k_history = 2
        ensemble_accepted_history = .true.
    end subroutine build_fixture

    subroutine test_reconciliation_report_mode_no_transitive_merge()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_REPORT, report_overlap_coefficient=.true._c_bool, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 3_int32, "report mode: 3 intersecting pairs (1-2, 2-3, 5-6)")
        call assert_equal_array_int(super_ensembles(:, 1), [1, 2], 2_int32, "report mode: pair 1")
        call assert_equal_array_int(super_ensembles(:, 2), [2, 3], 2_int32, "report mode: pair 2")
        call assert_equal_array_int(super_ensembles(:, 3), [5, 6], 2_int32, "report mode: pair 3")
        call assert_equal_real(super_ensembles_overlap_coefficient(1, 1), 0.5d0, 1.0d-9, "report mode: OC(1,2)")
        call assert_equal_real(super_ensembles_overlap_coefficient(1, 2), 0.5d0, 1.0d-9, "report mode: OC(2,3)")
        call assert_equal_real(super_ensembles_overlap_coefficient(1, 3), 2.0d0/3.0d0, 1.0d-9, "report mode: OC(5,6)")
        call assert_true(all(eligible), "report mode: every ensemble eligible with no filter supplied")
    end subroutine test_reconciliation_report_mode_no_transitive_merge

    subroutine test_reconciliation_merge_any_transitive()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3), expected_group2(3)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_ANY, report_overlap_coefficient=.true._c_bool, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        ! 1-2-3 chain merges transitively into one 3-member group; 5-6 into a separate
        ! 2-member group (padded with 0); 4 is isolated and excluded entirely.
        call assert_equal_int(n_super_ensembles, 2_int32, "merge_any: two groups found")

        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, &
                                    "merge_any: group 1 is the transitively-merged 1-2-3 chain")
        call assert_equal_real(super_ensembles_overlap_coefficient(1, 1), 0.5d0, 1.0d-9, "merge_any: group 1 OC(1,2)")
        call assert_equal_real(super_ensembles_overlap_coefficient(2, 1), 0.5d0, 1.0d-9, "merge_any: group 1 OC(2,3)")

        expected_group2 = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 2), expected_group2, 3_int32, &
                                    "merge_any: group 2 is {5,6}, padded with 0")
        call assert_equal_real(super_ensembles_overlap_coefficient(1, 2), 2.0d0/3.0d0, 1.0d-9, "merge_any: group 2 OC(5,6)")
    end subroutine test_reconciliation_merge_any_transitive

    !> min_overlap_coefficient=0.6 exceeds the 1-2-3 chain's OC (0.5), excluding it entirely,
    !| but is still below the separate 5-6 pair's OC (2/3), which still qualifies.
    subroutine test_reconciliation_merge_oc_excludes_weak_chain()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group(3)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_OVERLAP_COEFFICIENT, &
                                     min_overlap_coefficient=0.6d0, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 1_int32, "merge_overlap_coefficient: only the weaker chain is excluded")
        expected_group = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group, 3_int32, &
                                    "merge_overlap_coefficient: the surviving group is {5,6}")
    end subroutine test_reconciliation_merge_oc_excludes_weak_chain

    !> min_overlap_coefficient=0.4 is below both the chain's OC (0.5) and the pair's OC (2/3),
    !| so all three edges qualify -- same result as merge_any.
    subroutine test_reconciliation_merge_oc_includes_all()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_OVERLAP_COEFFICIENT, &
                                     min_overlap_coefficient=0.4d0, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, &
                              "merge_overlap_coefficient: threshold below every edge includes all")
        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, &
                                    "merge_overlap_coefficient: group 1 unchanged from merge_any")
    end subroutine test_reconciliation_merge_oc_includes_all

    !> report_overlap_coefficient defaults to .false.: super_ensembles_overlap_coefficient must
    !| stay all-zero even though real intersections (with genuine nonzero Overlap Coefficient)
    !| exist.
    subroutine test_reconciliation_oc_not_computed_unless_requested()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_ANY, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, "overlap coefficient not requested: groups still found")
        call assert_true(all(super_ensembles_overlap_coefficient == 0.0d0), &
                         "overlap coefficient not requested: super_ensembles_overlap_coefficient stays all-zero")
    end subroutine test_reconciliation_oc_not_computed_unless_requested

    !> The transitively-merged {1,2,3} group has 3 members; max_group_size=2 cannot hold it.
    subroutine test_reconciliation_group_exceeds_max_group_size()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_ANY, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "a 3-member group must not fit in max_group_size=2")
    end subroutine test_reconciliation_group_exceeds_max_group_size

    subroutine test_reconciliation_n_ensembles_too_small()
        logical(c_bool)        :: ensemble_masks(14, 1)
        integer(int32) :: ensemble_stop_reason(1)
        real(real64)   :: ensemble_U_history(2, 2, 1, 1), ensemble_S_history(2, 1, 1), ensemble_mu_history(2, 1, 1)
        real(real64)   :: ensemble_G_history(1, 1)
        integer(int32) :: ensemble_d_history(1, 1), ensemble_k_history(1, 1)
        logical(c_bool)        :: ensemble_accepted_history(1, 1)
        integer(int32) :: super_ensembles(2, 0), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 0)
        logical(c_bool)        :: eligible(1), eligible_by_stop_condition(1), eligible_by_dimension(1)
        logical(c_bool)        :: eligible_by_var_explained(1)

        ensemble_masks = .false.
        ensemble_stop_reason = STOP_REASON_FIXED_POINT
        ensemble_U_history = 0.0d0
        ensemble_d_history = 0
        ensemble_S_history = 0.0d0
        ensemble_mu_history = 0.0d0
        ensemble_G_history = 0.0d0
        ensemble_k_history = 2
        ensemble_accepted_history = .true.

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=1_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_reconciliation should reject n_ensembles=1")
    end subroutine test_reconciliation_n_ensembles_too_small

    subroutine test_reconciliation_invalid_mode()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=99_int32, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_reconciliation should reject an unknown mode")
    end subroutine test_reconciliation_invalid_mode

    !> Ensemble 5 is STOP_REASON_REJECTED_IMMEDIATELY (the rest STOP_REASON_FIXED_POINT);
    !| excluding STOP_REASON_REJECTED_IMMEDIATELY from `allowed_stop_reasons` must drop the
    !| (5,6) pair from report mode's output entirely (5 is ineligible, so no edge touching it
    !| is ever considered, regardless of 6's own eligibility), leaving only the untouched 1-2-3
    !| chain's two pairs -- n_super_ensembles goes from 3 (see the no-filter report-mode test
    !| above) to 2.
    subroutine test_reconciliation_stop_reason_filter_excludes_pair()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6), overrides(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)
        logical(c_bool)        :: allowed(4)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        overrides = 0
        overrides(5) = STOP_REASON_REJECTED_IMMEDIATELY
        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history, overrides)

        allowed = [.true., .true., .false., .true.]
        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_REPORT, allowed_stop_reasons=allowed, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, &
                              "allowed_stop_reasons: excluding ensemble 5 drops the (5,6) pair, only 1-2 and 2-3 remain")
        call assert_equal_array_int(super_ensembles(:, 1), [1, 2], 2_int32, "allowed_stop_reasons: pair 1 unaffected")
        call assert_equal_array_int(super_ensembles(:, 2), [2, 3], 2_int32, "allowed_stop_reasons: pair 2 unaffected")
        call assert_true(.not. eligible_by_stop_condition(5), &
                         "allowed_stop_reasons: ensemble 5 flagged ineligible by stop condition")
        call assert_true(.not. eligible(5), "allowed_stop_reasons: ensemble 5's combined eligibility is false")
    end subroutine test_reconciliation_stop_reason_filter_excludes_pair

    !> Omitting `allowed_stop_reasons` entirely, even with a genuinely mixed set of Stop
    !| Conditions across ensembles (not all STOP_REASON_FIXED_POINT), must behave identically
    !| to the no-filter `test_reconciliation_merge_any_transitive` case above -- proving
    !| "absent" truly means "every Stop Condition allowed", not an accidental all-false default.
    subroutine test_reconciliation_stop_reason_filter_absent_is_noop()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6), overrides(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3), expected_group2(3)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        overrides = 0
        overrides(1) = STOP_REASON_MAX_SIZE
        overrides(4) = STOP_REASON_REJECTED_AFTER_STABLE
        overrides(6) = STOP_REASON_REJECTED_IMMEDIATELY
        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history, overrides)

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_ANY, report_overlap_coefficient=.true._c_bool, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, &
                              "allowed_stop_reasons absent: mixed stop reasons still merge exactly as merge_any")
        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, &
                                    "allowed_stop_reasons absent: group 1 unchanged")
        expected_group2 = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 2), expected_group2, 3_int32, &
                                    "allowed_stop_reasons absent: group 2 unchanged")
        call assert_true(all(eligible), "allowed_stop_reasons absent: every ensemble still eligible")
    end subroutine test_reconciliation_stop_reason_filter_absent_is_noop

    !> Ensemble 2 (the sole bridge of the 1-2-3 chain) is STOP_REASON_REJECTED_IMMEDIATELY;
    !| excluding it must fully break the chain into non-merged singletons -- 1 and 3 do not
    !| intersect directly (see build_fixture's own docstring), so nothing bridges them once 2
    !| is ineligible. Only the untouched {5,6} group should survive.
    subroutine test_reconciliation_stop_reason_filter_breaks_chain()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6), overrides(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group(3)
        logical(c_bool)        :: allowed(4)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        overrides = 0
        overrides(2) = STOP_REASON_REJECTED_IMMEDIATELY
        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history, overrides)

        allowed = [.true., .true., .false., .true.]
        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_ANY, allowed_stop_reasons=allowed, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 1_int32, &
                              "allowed_stop_reasons: excluding the chain's bridge (2) leaves only {5,6}")
        expected_group = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group, 3_int32, &
                                    "allowed_stop_reasons: surviving group is {5,6}, unmerged with the broken chain")
    end subroutine test_reconciliation_stop_reason_filter_breaks_chain

    !> Every ensemble in the fixture is STOP_REASON_FIXED_POINT; excluding a Stop Condition
    !| value that matches none of them (STOP_REASON_REJECTED_AFTER_STABLE) must be a true
    !| no-op, identical to the no-filter merge_any result.
    subroutine test_reconciliation_stop_reason_filter_unused_value_noop()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3), expected_group2(3)
        logical(c_bool)        :: allowed(4)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)

        allowed = [.true., .false., .true., .true.]
        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_MERGE_ANY, allowed_stop_reasons=allowed, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, &
                              "allowed_stop_reasons: excluding an unused Stop Condition is a no-op")
        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, &
                                    "allowed_stop_reasons no-op: group 1 unchanged")
        expected_group2 = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 2), expected_group2, 3_int32, &
                                    "allowed_stop_reasons no-op: group 2 unchanged")
    end subroutine test_reconciliation_stop_reason_filter_unused_value_noop

    !> Ensemble 5's final intrinsic dimension (2) exceeds filter_dim_max=1 (the rest are d=1);
    !| excluding it via `filter_dim_max` must drop the (5,6) pair from report mode's output
    !| entirely, leaving only the untouched 1-2-3 chain's two pairs -- mirrors
    !| `test_reconciliation_stop_reason_filter_excludes_pair`'s own shape, but through
    !| the dimension filter instead of the Stop Condition filter. Uses its own dedicated
    !| single-column (o=1) history fixture, since `build_fixture`'s own history is
    !| deliberately uniform/unused by the other tests above.
    subroutine test_reconciliation_dimension_filter_excludes_pair()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)
        ensemble_d_history(1, :) = 1
        ensemble_d_history(1, 5) = 2 ! only ensemble 5 exceeds filter_dim_max=1 below

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_REPORT, filter_dim_max=1_int32, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, &
                              "filter_dim_max filter: excluding ensemble 5 (d=2 > filter_dim_max=1) drops the (5,6) pair, 1-2/2-3 remain")
        call assert_equal_array_int(super_ensembles(:, 1), [1, 2], 2_int32, "filter_dim_max filter: pair 1 unaffected")
        call assert_equal_array_int(super_ensembles(:, 2), [2, 3], 2_int32, "filter_dim_max filter: pair 2 unaffected")
        call assert_true(.not. eligible_by_dimension(5), "filter_dim_max filter: ensemble 5 flagged ineligible by dimension")
        call assert_true(eligible_by_dimension(1), "filter_dim_max filter: ensemble 1 still eligible by dimension")
        call assert_true(.not. eligible(5), "filter_dim_max filter: ensemble 5's combined eligibility is false")
    end subroutine test_reconciliation_dimension_filter_excludes_pair

    !> Ensemble 5's final variance explained (eigenvalues [1,100], d=1 -> 1/101 ~ 0.0099) falls
    !| far short of var_explained_min=0.5; the rest ([100,1] -> 100/101 ~ 0.99) comfortably
    !| clear it. Excluding ensemble 5 must drop the (5,6) pair, mirroring the dimension-filter
    !| test above. k=2 throughout keeps the k-1=1 eigenvalue denominator simple.
    subroutine test_reconciliation_var_explained_filter_excludes_pair()
        logical(c_bool)        :: ensemble_masks(14, 6)
        integer(int32) :: ensemble_stop_reason(6)
        real(real64)   :: ensemble_U_history(2, 2, 1, 6), ensemble_S_history(2, 1, 6), ensemble_mu_history(2, 1, 6)
        real(real64)   :: ensemble_G_history(1, 6)
        integer(int32) :: ensemble_d_history(1, 6), ensemble_k_history(1, 6)
        logical(c_bool)        :: ensemble_accepted_history(1, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)
        logical(c_bool)        :: eligible(6), eligible_by_stop_condition(6), eligible_by_dimension(6)
        logical(c_bool)        :: eligible_by_var_explained(6)
        integer(int32) :: e

        call build_fixture(ensemble_masks, ensemble_stop_reason, ensemble_U_history, ensemble_d_history, &
                           ensemble_S_history, ensemble_mu_history, ensemble_G_history, ensemble_k_history, &
                           ensemble_accepted_history)
        ensemble_d_history(1, :) = 1
        do e = 1, 6
            ensemble_S_history(:, 1, e) = [10.0d0, 1.0d0] ! eigenvalues [100,1], ve = 100/101
        end do
        ensemble_S_history(:, 1, 5) = [1.0d0, 10.0d0] ! eigenvalues [1,100], ve = 1/101

        call ensemble_reconciliation(ensemble_masks=ensemble_masks, ensemble_stop_reason=ensemble_stop_reason, &
                                     n_dimensions=2_int32, n_vectors=14_int32, n_ensembles=6_int32, &
                                     ensemble_U_history=ensemble_U_history, ensemble_d_history=ensemble_d_history, &
                                     ensemble_S_history=ensemble_S_history, ensemble_mu_history=ensemble_mu_history, &
                                     ensemble_G_history=ensemble_G_history, ensemble_k_history=ensemble_k_history, &
                                     ensemble_accepted_history=ensemble_accepted_history, o=1_int32, &
                                     mode=MODE_REPORT, var_explained_min=0.5d0, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, &
                                     eligible=eligible, eligible_by_stop_condition=eligible_by_stop_condition, &
                                     eligible_by_dimension=eligible_by_dimension, &
                                     eligible_by_var_explained=eligible_by_var_explained, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, &
                              "var_explained filter: excluding ensemble 5 (ve~1/101 < 0.5) drops the (5,6) pair")
        call assert_equal_array_int(super_ensembles(:, 1), [1, 2], 2_int32, "var_explained filter: pair 1 unaffected")
        call assert_equal_array_int(super_ensembles(:, 2), [2, 3], 2_int32, "var_explained filter: pair 2 unaffected")
        call assert_true(.not. eligible_by_var_explained(5), &
                         "var_explained filter: ensemble 5 flagged ineligible by variance explained")
        call assert_true(eligible_by_var_explained(1), "var_explained filter: ensemble 1 still eligible")
    end subroutine test_reconciliation_var_explained_filter_excludes_pair

    !> Direct test of the newly-independent merge_to_super_ensembles: an all-`.true.` `eligible`
    !| mask over the same 6-ensemble fixture must reproduce
    !| `test_reconciliation_merge_any_transitive`'s own result exactly, with no history-array
    !| fixture involved at all (this kernel only ever consumes `ensemble_masks`/`eligible`).
    subroutine test_merge_to_super_ensembles_all_eligible_matches_merge_any()
        logical(c_bool)        :: ensemble_masks(14, 6), eligible(6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3), expected_group2(3)

        ensemble_masks = .false.
        ensemble_masks(1:4, 1) = .true.
        ensemble_masks(3:6, 2) = .true.
        ensemble_masks(5:8, 3) = .true.
        ensemble_masks(9:10, 4) = .true.
        ensemble_masks(11:13, 5) = .true.
        ensemble_masks(12:14, 6) = .true.
        eligible = .true.

        call merge_to_super_ensembles(ensemble_masks=ensemble_masks, eligible=eligible, n_vectors=14_int32, &
                                      n_ensembles=6_int32, mode=MODE_MERGE_ANY, report_overlap_coefficient=.true._c_bool, &
                                      max_group_size=3_int32, super_ensembles=super_ensembles, &
                                      n_super_ensembles=n_super_ensembles, &
                                      super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'merge_to_super_ensembles failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, "merge_to_super_ensembles: two groups found, all eligible")
        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, &
                                    "merge_to_super_ensembles: group 1 is the 1-2-3 chain")
        expected_group2 = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 2), expected_group2, 3_int32, &
                                    "merge_to_super_ensembles: group 2 is {5,6}")
    end subroutine test_merge_to_super_ensembles_all_eligible_matches_merge_any

    !> Ensemble 5 marked ineligible via a hand-constructed `eligible` mask (no filtering kernel
    !| involved at all): the (5,6) edge must never be considered, leaving only the untouched
    !| 1-2-3 chain -- proves `merge_to_super_ensembles` honors `eligible` as a plain input,
    !| independent of however a caller chose to compute it.
    subroutine test_merge_to_super_ensembles_excludes_ineligible_ensemble()
        logical(c_bool)        :: ensemble_masks(14, 6), eligible(6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3)

        ensemble_masks = .false.
        ensemble_masks(1:4, 1) = .true.
        ensemble_masks(3:6, 2) = .true.
        ensemble_masks(5:8, 3) = .true.
        ensemble_masks(9:10, 4) = .true.
        ensemble_masks(11:13, 5) = .true.
        ensemble_masks(12:14, 6) = .true.
        eligible = .true.
        eligible(5) = .false.

        call merge_to_super_ensembles(ensemble_masks=ensemble_masks, eligible=eligible, n_vectors=14_int32, &
                                      n_ensembles=6_int32, mode=MODE_MERGE_ANY, report_overlap_coefficient=.true._c_bool, &
                                      max_group_size=3_int32, super_ensembles=super_ensembles, &
                                      n_super_ensembles=n_super_ensembles, &
                                      super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'merge_to_super_ensembles failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 1_int32, &
                              "merge_to_super_ensembles: only the chain group survives ensemble 5's exclusion")
        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, &
                                    "merge_to_super_ensembles: group is the 1-2-3 chain, {5,6} never considered")
    end subroutine test_merge_to_super_ensembles_excludes_ineligible_ensemble

end module mod_test_shape_truthful_clustering_reconciliation
