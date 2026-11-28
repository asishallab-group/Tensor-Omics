!> What the children of the `tox_paralog_analysis` suite share: the word whose only set bit is
!| bit 31, what an output holds before a call, the three non-finite values every procedure has
!| to refuse, and the check that compares the subsets a detection found with the expected ones.
!|
!| Bit masks are written out as words in f42's published layout (`f42_bit_masks_impl`): gene `i`
!| is bit `mod(i - 1, 32)` of word `(i - 1)/32 + 1`, least significant first, all 32 bits used.
!| So genes 1 and 3 are the word 5, gene 32 alone is the word `BIT_31`, and gene 33 is bit 0 of
!| the second word.
module mod_test_tox_paralog_analysis_fixtures
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use asserts, only: assert_equal_int, operator(//)
    use f42_math_impl, only: PI
    implicit none
    private
    public :: BIT_31, UNWRITTEN_COUNT, UNWRITTEN_WORD, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    public :: check_subsets
    public :: MIXED_N_DIMS, MIXED_N_GENES, MIXED_ANCESTOR, MIXED_GENES, MIXED_ANGLES

    !> The word whose only set bit is bit 31, 0x80000000: -2**31 as an int32, and `NA_integer_`
    !| in R. Not written `-huge(0_int32) - 1_int32`, which `-pedantic` flags as outside the
    !| symmetric range the standard implies.
    integer(int32), parameter :: BIT_31 = ibset(0_int32, 31)

    !> What a count holds before a call, so that a count the procedure never writes shows.
    integer(int32), parameter :: UNWRITTEN_COUNT = -777_int32
    !> What a bit mask word holds before a call: every bit set, which no expected word in this
    !| suite is, so a word the procedure never writes shows.
    integer(int32), parameter :: UNWRITTEN_WORD = -1_int32

    !> How many non-finite values `non_finite_value` hands out.
    integer(int32), parameter :: N_NON_FINITE = 3
    !> The names of the non-finite values, for the assertion messages.
    character(*), parameter :: NON_FINITE_NAMES(N_NON_FINITE) = ["NaN ", "+Inf", "-Inf"]

    !> A family holding both patterns, after the old suite's `test_detect_patterns_mixed_results`
    !| and issue #64's T6: the ancestor [1, 1, 0]; genes 1 and 2, [1, 0, 0] and [0, 1, 0], split it
    !| (their sum is the ancestor), and genes 3 and 4, [0.5, 0.5, 0] and [0.75, 0.75, 0], point
    !| along it and together are 1.25 times as long. Their angles to the ancestor, as a caller
    !| would hand them to the filters, are pi/4, pi/4, 0 and 0.
    integer(int32), parameter :: MIXED_N_DIMS = 3, MIXED_N_GENES = 4
    real(real64), parameter :: MIXED_ANCESTOR(MIXED_N_DIMS) = [1.0_real64, 1.0_real64, 0.0_real64]
    real(real64), parameter :: MIXED_GENES(MIXED_N_DIMS, MIXED_N_GENES) = reshape([ &
                               1.0_real64, 0.0_real64, 0.0_real64, &
                               0.0_real64, 1.0_real64, 0.0_real64, &
                               0.5_real64, 0.5_real64, 0.0_real64, &
                               0.75_real64, 0.75_real64, 0.0_real64], [MIXED_N_DIMS, MIXED_N_GENES])
    real(real64), parameter :: MIXED_ANGLES(MIXED_N_GENES) = [PI/4, PI/4, 0.0_real64, 0.0_real64]

contains

    !> The `i_value`-th non-finite value: NaN, +Inf, -Inf.
    function non_finite_value(i_value) result(value)
        integer(int32), intent(in) :: i_value
        real(real64) :: value

        select case (i_value)
        case (1)
            value = ieee_value(1.0_real64, ieee_quiet_nan)
        case (2)
            value = ieee_value(1.0_real64, ieee_positive_inf)
        case default
            value = ieee_value(1.0_real64, ieee_negative_inf)
        end select
    end function non_finite_value

    !> Checks that a detection found exactly the subsets in `expected`, in any order: `n_results`
    !| equals `n_expected`, and each expected subset is among the first `n_results` columns of
    !| `subsets` exactly once. The expected subsets are distinct, so together that is a one-to-one
    !| match. The order is left open because the detection stores a result when it finds it, and
    !| the order it finds them in is no documented property.
    subroutine check_subsets(n_mask_chunks, n_slots, subsets, n_results, n_expected, expected, msg)
        integer(int32), intent(in) :: n_mask_chunks, n_slots, n_results, n_expected
        integer(int32), intent(in) :: subsets(n_mask_chunks, n_slots), expected(n_mask_chunks, n_expected)
        character(*), intent(in) :: msg

        integer(int32) :: i_expected, i_result, n_found

        call assert_equal_int(n_results, n_expected, msg//": n_results")
        if (n_results /= n_expected .or. n_results > n_slots) return

        do i_expected = 1, n_expected
            n_found = 0
            do i_result = 1, n_results
                if (all(subsets(:, i_result) == expected(:, i_expected))) n_found = n_found + 1
            end do
            call assert_equal_int(n_found, 1_int32, msg//": expected subset "//i_expected//" (first word "// &
                                  expected(1, i_expected)//") found once")
        end do
    end subroutine check_subsets
end module mod_test_tox_paralog_analysis_fixtures
