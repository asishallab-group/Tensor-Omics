!> Unit test suite for the calc_percentile family (f42_stats).
!|
!| These tests pin the contract introduced by issue #149: the `percentile`
!| argument is a fraction/probability in the statistical interval [0,1]
!| (e.g. 0.95 for the 95th percentile), NOT a percentage in [0,100].
!|
!| The three tiers are exercised separately, because each states the contract
!| in its own way: `calc_percentile` sorts for you, `calc_percentile_expert`
!| takes the permutation you supply, and `calc_percentile_impl` validates
!| nothing at all.
module mod_test_percentile
    use asserts
    use f42_utils_impl
    use f42_stats, only: calc_percentile, calc_percentile_expert
    use tox_errors
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 1d-12

    !> `percentile` is the 4th dummy of the expert tier and the 3rd of the allocating one,
    !| which is the position a rejected value is reported against.
    integer(int32), parameter :: PERCENTILE_ARG_EXPERT = 4_int32
    integer(int32), parameter :: PERCENTILE_ARG_PLAIN = 3_int32

contains

    !> Get array of all available percentile tests.
    function get_all_tests_percentile() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(8))
        all_tests(1) = test_case("test_percentile_min_max_median", test_percentile_min_max_median)
        all_tests(2) = test_case("test_percentile_interpolation", test_percentile_interpolation)
        all_tests(3) = test_case("test_percentile_single_element", test_percentile_single_element)
        all_tests(4) = test_case("test_percentile_rejects_percent_scale", test_percentile_rejects_percent_scale)
        all_tests(5) = test_case("test_percentile_boundaries_inclusive", test_percentile_boundaries_inclusive)
        all_tests(6) = test_case("test_percentile_sorts_unsorted_input", test_percentile_sorts_unsorted_input)
        all_tests(7) = test_case("test_percentile_plain_rejects_percent_scale", test_percentile_plain_rejects_percent_scale)
        all_tests(8) = test_case("test_percentile_impl_interpolation", test_percentile_impl_interpolation)
    end function get_all_tests_percentile

    !> A fraction in [0,1] must map min -> 0.0, max -> 1.0, median -> 0.5.
    subroutine test_percentile_min_max_median()
        real(real64) :: array(5) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
        integer(int32) :: perm(5)
        real(real64) :: value
        integer(int32) :: ierr

        ! array is already sorted -> identity permutation is the sorting order
        call init_perm(perm)

        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.0_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_min_max_median: p=0.0 ierr should be OK")
        call assert_equal_real(value, 10.0_real64, TOL, "test_percentile_min_max_median: p=0.0 should return the minimum")

        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.5_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_min_max_median: p=0.5 ierr should be OK")
        call assert_equal_real(value, 30.0_real64, TOL, "test_percentile_min_max_median: p=0.5 should return the median")

        call calc_percentile_expert(array, size(array, kind=int32), perm, 1.0_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_min_max_median: p=1.0 ierr should be OK")
        call assert_equal_real(value, 50.0_real64, TOL, "test_percentile_min_max_median: p=1.0 should return the maximum")
    end subroutine test_percentile_min_max_median

    !> Linear interpolation between adjacent sorted values for fractional ranks.
    subroutine test_percentile_interpolation()
        real(real64) :: array(5) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
        integer(int32) :: perm(5)
        real(real64) :: value
        integer(int32) :: ierr

        call init_perm(perm)

        ! rank = p*(n-1)+1 with n=5
        ! p=0.10 -> rank=1.4 -> 10 + 0.4*(20-10) = 14
        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.10_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_interpolation: p=0.10 ierr should be OK")
        call assert_equal_real(value, 14.0_real64, TOL, "test_percentile_interpolation: p=0.10 should interpolate to 14.0")

        ! p=0.25 -> rank=2.0 -> exactly 20
        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.25_real64, value, ierr=ierr)
        call assert_equal_real(value, 20.0_real64, TOL, "test_percentile_interpolation: p=0.25 should be 20.0")

        ! p=0.95 -> rank=4.8 -> 40 + 0.8*(50-40) = 48
        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.95_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_interpolation: p=0.95 ierr should be OK")
        call assert_equal_real(value, 48.0_real64, TOL, "test_percentile_interpolation: p=0.95 should interpolate to 48.0")
    end subroutine test_percentile_interpolation

    !> A single-element array returns that element for any percentile.
    subroutine test_percentile_single_element()
        real(real64) :: array(1) = [42.0_real64]
        integer(int32) :: perm(1)
        real(real64) :: value
        integer(int32) :: ierr

        call init_perm(perm)

        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.0_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_single_element: p=0.0 ierr should be OK")
        call assert_equal_real(value, 42.0_real64, TOL, "test_percentile_single_element: p=0.0 should return the only element")

        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.73_real64, value, ierr=ierr)
        call assert_equal_real(value, 42.0_real64, TOL, "test_percentile_single_element: p=0.73 should return the only element")
    end subroutine test_percentile_single_element

    !> Regression proof for issue #149: values on the old [0,100] percent
    !| scale (or otherwise outside [0,1]) must now be rejected as out of range.
    subroutine test_percentile_rejects_percent_scale()
        real(real64) :: array(5) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
        integer(int32) :: perm(5)
        real(real64) :: value
        integer(int32) :: ierr

        call init_perm(perm)

        ! 50.0 used to mean "the median"; on the [0,1] scale it is out of range.
        call calc_percentile_expert(array, size(array, kind=int32), perm, 50.0_real64, value, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_percentile_rejects_percent_scale: 50.0 must be rejected", &
                        arg_pos=PERCENTILE_ARG_EXPERT)

        ! Just above the upper bound.
        call calc_percentile_expert(array, size(array, kind=int32), perm, above(1.0_real64), value, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_percentile_rejects_percent_scale: >1 must be rejected", &
                        arg_pos=PERCENTILE_ARG_EXPERT)

        ! Just below the lower bound.
        call calc_percentile_expert(array, size(array, kind=int32), perm, below(0.0_real64), value, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_percentile_rejects_percent_scale: <0 must be rejected", &
                        arg_pos=PERCENTILE_ARG_EXPERT)
    end subroutine test_percentile_rejects_percent_scale

    !> The [0,1] interval is inclusive: exactly 0.0 and 1.0 are valid.
    subroutine test_percentile_boundaries_inclusive()
        real(real64) :: array(4) = [-1.0_real64, 0.0_real64, 2.0_real64, 5.0_real64]
        integer(int32) :: perm(4)
        real(real64) :: value
        integer(int32) :: ierr

        call init_perm(perm)

        call calc_percentile_expert(array, size(array, kind=int32), perm, 0.0_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_boundaries_inclusive: 0.0 must be accepted")
        call assert_equal_real(value, -1.0_real64, TOL, "test_percentile_boundaries_inclusive: p=0.0 returns minimum")

        call calc_percentile_expert(array, size(array, kind=int32), perm, 1.0_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_boundaries_inclusive: 1.0 must be accepted")
        call assert_equal_real(value, 5.0_real64, TOL, "test_percentile_boundaries_inclusive: p=1.0 returns maximum")
    end subroutine test_percentile_boundaries_inclusive

    !> The allocating tier sorts internally, so it accepts unsorted input and
    !| honours the same [0,1] contract.
    subroutine test_percentile_sorts_unsorted_input()
        real(real64) :: array(5) = [30.0_real64, 10.0_real64, 50.0_real64, 20.0_real64, 40.0_real64]
        real(real64) :: value
        integer(int32) :: ierr

        call calc_percentile(array, size(array, kind=int32), 0.5_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_sorts_unsorted_input: p=0.5 ierr should be OK")
        call assert_equal_real(value, 30.0_real64, TOL, "test_percentile_sorts_unsorted_input: p=0.5 should return the median")

        call calc_percentile(array, size(array, kind=int32), 0.95_real64, value, ierr=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_percentile_sorts_unsorted_input: p=0.95 ierr should be OK")
        call assert_equal_real(value, 48.0_real64, TOL, "test_percentile_sorts_unsorted_input: p=0.95 should interpolate to 48.0")
    end subroutine test_percentile_sorts_unsorted_input

    !> The allocating tier must also reject percent-scale input, at its own argument position.
    subroutine test_percentile_plain_rejects_percent_scale()
        real(real64) :: array(5) = [30.0_real64, 10.0_real64, 50.0_real64, 20.0_real64, 40.0_real64]
        real(real64) :: value
        integer(int32) :: ierr

        call calc_percentile(array, size(array, kind=int32), 95.0_real64, value, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_percentile_plain_rejects_percent_scale: 95.0 must be rejected", &
                        arg_pos=PERCENTILE_ARG_PLAIN)
    end subroutine test_percentile_plain_rejects_percent_scale

    !> The unvalidated implementation produces the same interpolated values for a
    !| caller-supplied sorted permutation.
    subroutine test_percentile_impl_interpolation()
        real(real64) :: array(5) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
        integer(int32) :: perm(5)
        real(real64) :: value

        call init_perm(perm)

        ! p=0.75 -> rank=4.0 -> exactly 40
        call calc_percentile_impl(array, size(array, kind=int32), perm, 0.75_real64, value)
        call assert_equal_real(value, 40.0_real64, TOL, "test_percentile_impl_interpolation: p=0.75 should be 40.0")

        ! p=0.5 -> rank=3.0 -> exactly 30 (median)
        call calc_percentile_impl(array, size(array, kind=int32), perm, 0.5_real64, value)
        call assert_equal_real(value, 30.0_real64, TOL, "test_percentile_impl_interpolation: p=0.5 should be 30.0")
    end subroutine test_percentile_impl_interpolation

end module mod_test_percentile
