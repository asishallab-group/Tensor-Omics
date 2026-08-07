!> Unit test suite for tox_shape_truthful_clustering_reconciliation (ensemble_reconciliation),
!| generated from
!| src/kernel/shape_truthful_clustering/tox_shape_truthful_clustering_reconciliation_kernel.F90.
module mod_test_shape_truthful_clustering_reconciliation
    use tox_shape_truthful_clustering_reconciliation, only: ensemble_reconciliation
    use tox_shape_truthful_clustering_reconciliation_kernel, only: MODE_REPORT, MODE_MERGE_JSI, MODE_MERGE_ANY
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
        all_tests(3) = test_case("test_reconciliation_merge_jsi_threshold_excludes_weak_chain", &
                                 test_reconciliation_merge_jsi_threshold_excludes_weak_chain)
        all_tests(4) = test_case("test_reconciliation_merge_jsi_threshold_includes_all", &
                                 test_reconciliation_merge_jsi_threshold_includes_all)
        all_tests(5) = test_case("test_reconciliation_jsi_not_computed_unless_requested", &
                                 test_reconciliation_jsi_not_computed_unless_requested)
        all_tests(6) = test_case("test_reconciliation_group_exceeds_max_group_size", &
                                 test_reconciliation_group_exceeds_max_group_size)
        all_tests(7) = test_case("test_reconciliation_n_ensembles_too_small", &
                                 test_reconciliation_n_ensembles_too_small)
        all_tests(8) = test_case("test_reconciliation_invalid_mode", test_reconciliation_invalid_mode)
    end function get_all_tests_shape_truthful_clustering_reconciliation

    !> N=14 vectors, 6 ensembles.
    !| E1={1,2,3,4}, E2={3,4,5,6}, E3={5,6,7,8}: a chain, E1-E2 and E2-E3 each intersect at 2
    !| members (JSI = 2/6 = 1/3), E1-E3 do not intersect at all.
    !| E4={9,10}: isolated, intersects nobody.
    !| E5={11,12,13}, E6={12,13,14}: a separate pair, intersecting at 2 members
    !| (JSI = 2/4 = 1/2).
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
        real(real64)   :: super_ensembles_JSI(1, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_REPORT, report_jsi=.true., max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 3_int32, "report mode: 3 intersecting pairs (1-2, 2-3, 5-6)")
        call assert_equal_array_int(super_ensembles(:, 1), [1, 2], 2_int32, "report mode: pair 1")
        call assert_equal_array_int(super_ensembles(:, 2), [2, 3], 2_int32, "report mode: pair 2")
        call assert_equal_array_int(super_ensembles(:, 3), [5, 6], 2_int32, "report mode: pair 3")
        call assert_equal_real(super_ensembles_JSI(1, 1), 1.0d0/3.0d0, 1.0d-9, "report mode: JSI(1,2)")
        call assert_equal_real(super_ensembles_JSI(1, 2), 1.0d0/3.0d0, 1.0d-9, "report mode: JSI(2,3)")
        call assert_equal_real(super_ensembles_JSI(1, 3), 0.5d0, 1.0d-9, "report mode: JSI(5,6)")
    end subroutine test_reconciliation_report_mode_no_transitive_merge

    subroutine test_reconciliation_merge_any_transitive()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(2, 30)
        integer(int32) :: expected_group1(3), expected_group2(3)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_ANY, report_jsi=.true., max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
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
        call assert_equal_real(super_ensembles_JSI(1, 1), 1.0d0/3.0d0, 1.0d-9, "merge_any: group 1 JSI(1,2)")
        call assert_equal_real(super_ensembles_JSI(2, 1), 1.0d0/3.0d0, 1.0d-9, "merge_any: group 1 JSI(2,3)")

        expected_group2 = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 2), expected_group2, 3_int32, &
                                    "merge_any: group 2 is {5,6}, padded with 0")
        call assert_equal_real(super_ensembles_JSI(1, 2), 0.5d0, 1.0d-9, "merge_any: group 2 JSI(5,6)")
    end subroutine test_reconciliation_merge_any_transitive

    !> min_jsi=0.4 exceeds both chain JSIs (1/3 ~= 0.333), excluding the 1-2-3 chain
    !| entirely, but is still below the separate 5-6 pair's JSI (1/2), which still qualifies.
    subroutine test_reconciliation_merge_jsi_threshold_excludes_weak_chain()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(2, 30)
        integer(int32) :: expected_group(3)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_JSI, min_jsi=0.4d0, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 1_int32, "merge_jsi: only the weaker chain is excluded")
        expected_group = [5, 6, 0]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group, 3_int32, &
                                    "merge_jsi: the surviving group is {5,6}")
    end subroutine test_reconciliation_merge_jsi_threshold_excludes_weak_chain

    !> min_jsi=0.3 is below both chain JSIs (1/3), so both qualify -- same result as merge_any.
    subroutine test_reconciliation_merge_jsi_threshold_includes_all()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(2, 30)
        integer(int32) :: expected_group1(3)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_JSI, min_jsi=0.3d0, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, "merge_jsi: threshold below every edge includes all")
        expected_group1 = [1, 2, 3]
        call assert_equal_array_int(super_ensembles(:, 1), expected_group1, 3_int32, "merge_jsi: group 1 unchanged from merge_any")
    end subroutine test_reconciliation_merge_jsi_threshold_includes_all

    !> report_jsi defaults to .false.: super_ensembles_JSI must stay all-zero even though real
    !| intersections (with genuine nonzero JSI) exist.
    subroutine test_reconciliation_jsi_not_computed_unless_requested()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(2, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_ANY, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        if (.not. is_ok(ierr)) then
            write (*, *) 'ensemble_reconciliation failed unexpectedly: ', ierr
            error stop
        end if

        call assert_equal_int(n_super_ensembles, 2_int32, "jsi not requested: groups still found")
        call assert_true(all(super_ensembles_JSI == 0.0d0), "jsi not requested: super_ensembles_JSI stays all-zero")
    end subroutine test_reconciliation_jsi_not_computed_unless_requested

    !> The transitively-merged {1,2,3} group has 3 members; max_group_size=2 cannot hold it.
    subroutine test_reconciliation_group_exceeds_max_group_size()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(2, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(1, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=MODE_MERGE_ANY, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "a 3-member group must not fit in max_group_size=2")
    end subroutine test_reconciliation_group_exceeds_max_group_size

    subroutine test_reconciliation_n_ensembles_too_small()
        logical        :: ensemble_masks(14, 1)
        integer(int32) :: super_ensembles(2, 0), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(1, 0)

        ensemble_masks = .false.

        call ensemble_reconciliation(ensemble_masks, 14_int32, 1_int32, max_group_size=2_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_reconciliation should reject n_ensembles=1")
    end subroutine test_reconciliation_n_ensembles_too_small

    subroutine test_reconciliation_invalid_mode()
        logical        :: ensemble_masks(14, 6)
        integer(int32) :: super_ensembles(3, 30), n_super_ensembles, ierr
        real(real64)   :: super_ensembles_JSI(2, 30)

        call build_fixture(ensemble_masks)

        call ensemble_reconciliation(ensemble_masks, 14_int32, 6_int32, mode=99_int32, max_group_size=3_int32, &
                                     super_ensembles=super_ensembles, n_super_ensembles=n_super_ensembles, &
                                     super_ensembles_JSI=super_ensembles_JSI, ierr=ierr)
        call assert_true(is_err(ierr), "ensemble_reconciliation should reject an unknown mode")
    end subroutine test_reconciliation_invalid_mode

end module mod_test_shape_truthful_clustering_reconciliation
