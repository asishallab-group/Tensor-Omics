!> The test framework's own cases: an assertion must fail when it should.
!|
!| Every comparison with NaN is false, so an assertion written as "fail if the difference
!| exceeds the tolerance" lets NaN through. These cases call each real and complex assertion
!| with a NaN and check that exactly one failure was recorded, then take that failure back
!| with `forgive_assertion_failures`, so the case itself passes. The expected failures still
!| print their "ASSERTION FAILED" lines, each message beginning with "(expected)".
!| Controls check that equal values, equal infinities and "not equal" with NaN pass.
module mod_test_asserts
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use test_suite, only: test_case, assertion_failures_in_case, forgive_assertion_failures
    implicit none
    private
    public :: get_all_tests_asserts

contains

    !> Get array of all available tests.
    function get_all_tests_asserts() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(9))
        all_tests(1) = test_case("test_asserts_equal_real_fails_on_nan", test_asserts_equal_real_fails_on_nan)
        all_tests(2) = test_case("test_asserts_equal_array_real_fails_on_nan", test_asserts_equal_array_real_fails_on_nan)
        all_tests(3) = test_case("test_asserts_allclose_fails_on_nan", test_asserts_allclose_fails_on_nan)
        all_tests(4) = test_case("test_asserts_complex_fails_on_nan", test_asserts_complex_fails_on_nan)
        all_tests(5) = test_case("test_asserts_range_sorted_sum_fail_on_nan", test_asserts_range_sorted_sum_fail_on_nan)
        all_tests(6) = test_case("test_asserts_equal_infinities_pass", test_asserts_equal_infinities_pass)
        all_tests(7) = test_case("test_asserts_not_equal_passes_on_nan", test_asserts_not_equal_passes_on_nan)
        all_tests(8) = test_case("test_asserts_identical_accepts_nan_at_nan", test_asserts_identical_accepts_nan_at_nan)
        all_tests(9) = test_case("test_asserts_identical_rejects_nan_at_number", test_asserts_identical_rejects_nan_at_number)
    end function get_all_tests_asserts

    ! ------------------------------------------------------------------ helpers

    real(real64) function nan()
        nan = ieee_value(1.0_real64, ieee_quiet_nan)
    end function nan

    real(real64) function inf()
        inf = ieee_value(1.0_real64, ieee_positive_inf)
    end function inf

    !> The assertion called since `failures_before` must have recorded exactly one failure; it
    !| is taken back. Anything else fails the case, naming the assertion that let NaN through.
    subroutine expect_one_failure(failures_before, what)
        integer, intent(in) :: failures_before
        character(*), intent(in) :: what

        if (assertion_failures_in_case() == failures_before + 1) then
            call forgive_assertion_failures(1)
        else
            call assert_true(.false., what//" did not fail")
        end if
    end subroutine expect_one_failure

    ! ------------------------------------------------------------------ cases

    !> NaN as the value and as the expectation.
    subroutine test_asserts_equal_real_fails_on_nan()
        integer :: before

        before = assertion_failures_in_case()
        call assert_equal_real(nan(), 1.0_real64, 1.0e-12_real64, "(expected) NaN got")
        call expect_one_failure(before, "assert_equal_real(NaN, 1)")

        before = assertion_failures_in_case()
        call assert_equal_real(1.0_real64, nan(), 1.0e-12_real64, "(expected) NaN expected")
        call expect_one_failure(before, "assert_equal_real(1, NaN)")
    end subroutine test_asserts_equal_real_fails_on_nan

    !> A NaN among equal elements, on either side.
    subroutine test_asserts_equal_array_real_fails_on_nan()
        real(real64) :: good(3), bad(3)
        integer :: before

        good = [1.0_real64, 2.0_real64, 3.0_real64]
        bad = good
        bad(2) = nan()

        before = assertion_failures_in_case()
        call assert_equal_array_real(bad, good, 3, 1.0e-12_real64, "(expected) NaN got")
        call expect_one_failure(before, "assert_equal_array_real(NaN in a)")

        before = assertion_failures_in_case()
        call assert_equal_array_real(good, bad, 3, 1.0e-12_real64, "(expected) NaN expected")
        call expect_one_failure(before, "assert_equal_array_real(NaN in b)")
    end subroutine test_asserts_equal_array_real_fails_on_nan

    subroutine test_asserts_allclose_fails_on_nan()
        real(real64) :: good(3), bad(3)
        integer :: before

        good = [1.0_real64, 2.0_real64, 3.0_real64]
        bad = good
        bad(3) = nan()

        before = assertion_failures_in_case()
        call assert_allclose_array_real(bad, good, 3, 1.0e-12_real64, 1.0e-12_real64, "(expected) NaN got")
        call expect_one_failure(before, "assert_allclose_array_real(NaN in a)")
    end subroutine test_asserts_allclose_fails_on_nan

    !> A NaN real part, as a scalar and in an array.
    subroutine test_asserts_complex_fails_on_nan()
        complex(real64) :: good(2), bad(2)
        integer :: before

        good = [(1.0_real64, 1.0_real64), (2.0_real64, 0.0_real64)]
        bad = good
        bad(1) = cmplx(nan(), 1.0_real64, kind=real64)

        before = assertion_failures_in_case()
        call assert_equal_complex(bad(1), good(1), 1.0e-12_real64, "(expected) NaN got")
        call expect_one_failure(before, "assert_equal_complex(NaN, z)")

        before = assertion_failures_in_case()
        call assert_equal_array_complex(bad, good, 2, 1.0e-12_real64, "(expected) NaN got")
        call expect_one_failure(before, "assert_equal_array_complex(NaN in a)")
    end subroutine test_asserts_complex_fails_on_nan

    !> NaN is in no range, orders against nothing, and makes a sum NaN.
    subroutine test_asserts_range_sorted_sum_fail_on_nan()
        real(real64) :: values(3)
        integer :: before

        before = assertion_failures_in_case()
        call assert_in_range_real(nan(), 0.0_real64, 1.0_real64, "(expected) NaN in range")
        call expect_one_failure(before, "assert_in_range_real(NaN)")

        values = [1.0_real64, nan(), 3.0_real64]
        before = assertion_failures_in_case()
        call assert_sorted_real(values, 3, "(expected) NaN sorted")
        call expect_one_failure(before, "assert_sorted_real([1, NaN, 3])")

        before = assertion_failures_in_case()
        call assert_sum_equal(values, 3, 4.0_real64, "(expected) NaN sum")
        call expect_one_failure(before, "assert_sum_equal([1, NaN, 3])")
    end subroutine test_asserts_range_sorted_sum_fail_on_nan

    !> Control: two equal infinities are equal, although Inf - Inf is NaN.
    subroutine test_asserts_equal_infinities_pass()
        real(real64) :: values(2)

        call assert_equal_real(inf(), inf(), 0.0_real64, "equal infinities are equal")
        values = [1.0_real64, -inf()]
        call assert_equal_array_real(values, values, 2, 0.0_real64, "equal arrays with -Inf are equal")
    end subroutine test_asserts_equal_infinities_pass

    !> Control: NaN is not equal to anything, so assert_not_equal_real passes.
    subroutine test_asserts_not_equal_passes_on_nan()
        call assert_not_equal_real(nan(), 1.0_real64, 1.0e-12_real64, "NaN is not equal to 1")
    end subroutine test_asserts_not_equal_passes_on_nan

    !> Control: identical arrays, a NaN at the same position included, are identical.
    subroutine test_asserts_identical_accepts_nan_at_nan()
        real(real64) :: values(3)

        values = [1.0_real64, nan(), -inf()]
        call assert_identical_array_real(values, values, 3, "an array is identical to itself, NaN included")
    end subroutine test_asserts_identical_accepts_nan_at_nan

    !> A NaN where the other array holds a number is a difference.
    subroutine test_asserts_identical_rejects_nan_at_number()
        real(real64) :: values(3), changed(3)
        integer :: before

        values = [1.0_real64, 2.0_real64, 3.0_real64]
        changed = values
        changed(2) = nan()

        before = assertion_failures_in_case()
        call assert_identical_array_real(changed, values, 3, "(expected) NaN where a number was")
        call expect_one_failure(before, "assert_identical_array_real(NaN at a number)")
    end subroutine test_asserts_identical_rejects_nan_at_number

end module mod_test_asserts
