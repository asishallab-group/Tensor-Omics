!> Unit test suite for tox_shape_truthful_clustering_reconciliation (ensemble_reconciliation),
!| generated from
!| src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_reconciliation_kernel.F90.
module mod_test_shape_truthful_clustering_reconciliation
    use tox_shape_truthful_clustering_reconciliation, only: ensemble_reconciliation
    use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_REPORT, MODE_MERGE_OVERLAP_COEFFICIENT, MODE_MERGE_ANY
    use tox_errors, only: is_ok, is_err, ERR_SIZE_MISMATCH
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_shape_truthful_clustering_reconciliation() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(8))

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
    subroutine build_fixture(ensemble_masks)
        logical, intent(out) :: ensemble_masks(14, 6)

        ensemble_masks = .false.
        ensemble_masks(1:4, 1) = .true.
        ensemble_masks(3:6, 2) = .true.
        ensemble_masks(5:8, 3) = .true.
        ensemble_masks(9:10, 4) = .true.
        ensemble_masks(11:13, 5) = .true.
        ensemble_masks(12:14, 6) = .true.
    end subroutine build_fixture

    subroutine test_reconciliation_report_mode_no_transitive_merge()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_REPORT, &
                                     report_overlap_coefficient=.true., max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
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
    end subroutine test_reconciliation_report_mode_no_transitive_merge

    subroutine test_reconciliation_merge_any_transitive()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3), expected_group2(3)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_ANY, &
                                     report_overlap_coefficient=.true., max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
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
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group(3)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_OVERLAP_COEFFICIENT, &
                                     min_overlap_coefficient=0.6d0, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
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
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)
        integer(int32) :: expected_group1(3)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_OVERLAP_COEFFICIENT, &
                                     min_overlap_coefficient=0.4d0, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
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
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_ANY, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
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
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_ANY, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "a 3-member group must not fit in max_group_size=2")
    end subroutine test_reconciliation_group_exceeds_max_group_size

    subroutine test_reconciliation_n_ensembles_too_small()
        logical        :: ensemble_masks(14, 1)
        integer(int32) :: super_ensembles(2, 0), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(1, 0)

        ensemble_masks = .false.

        call ensemble_reconciliation(ensemble_masks, 14_int32, 1_int32, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_reconciliation should reject n_ensembles=1")
    end subroutine test_reconciliation_n_ensembles_too_small

    subroutine test_reconciliation_invalid_mode()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_overlap_coefficient(2, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=99_int32, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_overlap_coefficient=super_ensembles_overlap_coefficient, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_reconciliation should reject an unknown mode")
    end subroutine test_reconciliation_invalid_mode

end module mod_test_shape_truthful_clustering_reconciliation
