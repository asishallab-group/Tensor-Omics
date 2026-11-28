!> The f42 array checks of `f42_array_checks_impl`: `all_equal`, `any_below` and `any_above`, each
!| on int32 and real64, rank 1 and rank 2.
!|
!| Every expected value is the one `all(array == value)`, `any(array < bound)` or
!| `any(array > bound)` gives for the same elements, read off the inputs by hand. Each check is
!| placed so that a wrong answer has one cause: an element that decides at the first position
!| shows the early exit, one only at the last position shows that the loop reaches the end, and
!| one equal to the bound shows that the comparison is strict. No input holds a NaN, which the
!| helpers' contract excludes.
module mod_test_f42_array_checks
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    use f42_array_checks_impl
    use test_suite, only: test_case
    implicit none
    public

    !> -2^31, the most negative int32. Not `-huge(0_int32) - 1_int32`, which `-pedantic` flags as
    !| outside the symmetric range the standard implies.
    integer(int32), parameter :: MOST_NEGATIVE_INT32 = ibset(0_int32, 31)

contains

    !> Get array of all available tests.
    function get_all_tests_f42_array_checks() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(8))
        all_tests(1) = test_case("test_all_equal_int", test_all_equal_int)
        all_tests(2) = test_case("test_all_equal_real", test_all_equal_real)
        all_tests(3) = test_case("test_any_below_int", test_any_below_int)
        all_tests(4) = test_case("test_any_above_int", test_any_above_int)
        all_tests(5) = test_case("test_any_below_real", test_any_below_real)
        all_tests(6) = test_case("test_any_above_real", test_any_above_real)
        all_tests(7) = test_case("test_array_checks_sections", test_array_checks_sections)
        all_tests(8) = test_case("test_array_checks_rank2", test_array_checks_rank2)
    end function get_all_tests_f42_array_checks

    !> `all_equal` on int32: `.true.` for no elements and for equal ones, `.false.` when the first
    !| or only the last element differs; the int32 extremes are ordinary values.
    subroutine test_all_equal_int()
        integer(int32), allocatable :: empty(:)
        integer(int32) :: values(4)

        allocate (empty(0))
        call assert_true(all_equal(empty, 7_int32), "test_all_equal_int: no elements")
        values = 7_int32
        call assert_true(all_equal(values, 7_int32), "test_all_equal_int: [7, 7, 7, 7] vs 7")
        values(1) = 8_int32
        call assert_false(all_equal(values, 7_int32), "test_all_equal_int: [8, 7, 7, 7] vs 7, the first differs")
        values(1) = 7_int32
        values(4) = 8_int32
        call assert_false(all_equal(values, 7_int32), "test_all_equal_int: [7, 7, 7, 8] vs 7, only the last differs")

        values = MOST_NEGATIVE_INT32
        call assert_true(all_equal(values, MOST_NEGATIVE_INT32), "test_all_equal_int: all -2^31 vs -2^31")
        values(4) = huge(0_int32)
        call assert_false(all_equal(values, MOST_NEGATIVE_INT32), "test_all_equal_int: [-2^31, ..., huge] vs -2^31")
        values = huge(0_int32)
        call assert_true(all_equal(values, huge(0_int32)), "test_all_equal_int: all huge vs huge")
        values(4) = MOST_NEGATIVE_INT32
        call assert_false(all_equal(values, huge(0_int32)), "test_all_equal_int: [huge, ..., -2^31] vs huge")
    end subroutine test_all_equal_int

    !> `all_equal` on real64 compares exactly with `==`: `-0.0` equals `0.0`, `tiny` (the smallest
    !| normal number, so no compiler flushes it) does not, and an infinity equals only itself.
    subroutine test_all_equal_real()
        real(real64), allocatable :: empty(:)
        real(real64) :: zeros(3), pair(2), values(4), inf

        inf = ieee_value(1.0_real64, ieee_positive_inf)

        allocate (empty(0))
        call assert_true(all_equal(empty, 1.0_real64), "test_all_equal_real: no elements")
        zeros = [0.0_real64, -0.0_real64, 0.0_real64]
        call assert_true(all_equal(zeros, 0.0_real64), "test_all_equal_real: [0, -0, 0] vs 0")
        pair = [0.0_real64, tiny(1.0_real64)]
        call assert_false(all_equal(pair, 0.0_real64), "test_all_equal_real: [0, tiny] vs 0")
        pair = inf
        call assert_true(all_equal(pair, inf), "test_all_equal_real: [Inf, Inf] vs Inf")
        pair(2) = -inf
        call assert_false(all_equal(pair, inf), "test_all_equal_real: [Inf, -Inf] vs Inf")

        values = [2.0_real64, 1.5_real64, 1.5_real64, 1.5_real64]
        call assert_false(all_equal(values, 1.5_real64), "test_all_equal_real: [2, 1.5, 1.5, 1.5] vs 1.5, the first differs")
        values = [1.5_real64, 1.5_real64, 1.5_real64, 2.0_real64]
        call assert_false(all_equal(values, 1.5_real64), &
                          "test_all_equal_real: [1.5, 1.5, 1.5, 2] vs 1.5, only the last differs")
        values = 1.5_real64
        call assert_true(all_equal(values, 1.5_real64), "test_all_equal_real: [1.5, 1.5, 1.5, 1.5] vs 1.5")
    end subroutine test_all_equal_real

    !> `any_below` on int32: `.false.` for no elements and for an element equal to the bound,
    !| `.true.` for one below it at the first or only at the last position; nothing is below -2^31,
    !| and -2^31 is below every other bound.
    subroutine test_any_below_int()
        integer(int32), allocatable :: empty(:)
        integer(int32) :: values(4)

        allocate (empty(0))
        call assert_false(any_below(empty, 0_int32), "test_any_below_int: no elements")
        values = [5_int32, 3_int32, 4_int32, 3_int32]
        call assert_false(any_below(values, 3_int32), "test_any_below_int: [5, 3, 4, 3] vs 3, equal is not below")
        values = [2_int32, 3_int32, 4_int32, 3_int32]
        call assert_true(any_below(values, 3_int32), "test_any_below_int: [2, 3, 4, 3] vs 3, the first is below")
        values = [5_int32, 3_int32, 4_int32, 2_int32]
        call assert_true(any_below(values, 3_int32), "test_any_below_int: [5, 3, 4, 2] vs 3, only the last is below")

        values = MOST_NEGATIVE_INT32
        call assert_false(any_below(values, MOST_NEGATIVE_INT32), "test_any_below_int: all -2^31 vs -2^31")
        values = huge(0_int32)
        values(4) = MOST_NEGATIVE_INT32
        call assert_true(any_below(values, huge(0_int32)), "test_any_below_int: [huge, ..., -2^31] vs huge")
        call assert_true(any_below(values, MOST_NEGATIVE_INT32 + 1_int32), &
                         "test_any_below_int: [huge, ..., -2^31] vs -2^31 + 1")
    end subroutine test_any_below_int

    !> `any_above` on int32: `.false.` for no elements and for an element equal to the bound,
    !| `.true.` for one above it at the first or only at the last position; nothing is above huge,
    !| and huge is above every other bound.
    subroutine test_any_above_int()
        integer(int32), allocatable :: empty(:)
        integer(int32) :: values(4)

        allocate (empty(0))
        call assert_false(any_above(empty, 0_int32), "test_any_above_int: no elements")
        values = [1_int32, 3_int32, 2_int32, 3_int32]
        call assert_false(any_above(values, 3_int32), "test_any_above_int: [1, 3, 2, 3] vs 3, equal is not above")
        values = [4_int32, 3_int32, 2_int32, 3_int32]
        call assert_true(any_above(values, 3_int32), "test_any_above_int: [4, 3, 2, 3] vs 3, the first is above")
        values = [1_int32, 3_int32, 2_int32, 4_int32]
        call assert_true(any_above(values, 3_int32), "test_any_above_int: [1, 3, 2, 4] vs 3, only the last is above")

        values = huge(0_int32)
        call assert_false(any_above(values, huge(0_int32)), "test_any_above_int: all huge vs huge")
        values = MOST_NEGATIVE_INT32
        values(4) = huge(0_int32)
        call assert_true(any_above(values, MOST_NEGATIVE_INT32), "test_any_above_int: [-2^31, ..., huge] vs -2^31")
        call assert_true(any_above(values, huge(0_int32) - 1_int32), "test_any_above_int: [-2^31, ..., huge] vs huge - 1")
    end subroutine test_any_above_int

    !> `any_below` on real64: `.false.` for no elements, for an element equal to the bound and for
    !| `-0.0` against `0.0` (equal under IEEE); `.true.` at the first or only the last position.
    !| `-Inf` is below every finite bound, and nothing is below `-Inf`.
    subroutine test_any_below_real()
        real(real64), allocatable :: empty(:)
        real(real64) :: values(4), inf

        inf = ieee_value(1.0_real64, ieee_positive_inf)

        allocate (empty(0))
        call assert_false(any_below(empty, 0.0_real64), "test_any_below_real: no elements")
        values = [5.0_real64, 3.0_real64, 4.0_real64, 3.0_real64]
        call assert_false(any_below(values, 3.0_real64), "test_any_below_real: [5, 3, 4, 3] vs 3, equal is not below")
        values = [1.0_real64, -0.0_real64, 2.0_real64, -0.0_real64]
        call assert_false(any_below(values, 0.0_real64), "test_any_below_real: [1, -0, 2, -0] vs 0, -0 is not below 0")
        values = [2.5_real64, 3.0_real64, 4.0_real64, 3.0_real64]
        call assert_true(any_below(values, 3.0_real64), "test_any_below_real: [2.5, 3, 4, 3] vs 3, the first is below")
        values = [5.0_real64, 3.0_real64, 4.0_real64, 2.5_real64]
        call assert_true(any_below(values, 3.0_real64), "test_any_below_real: [5, 3, 4, 2.5] vs 3, only the last is below")

        values = 1.0_real64
        values(4) = -inf
        call assert_true(any_below(values, -huge(1.0_real64)), "test_any_below_real: [1, 1, 1, -Inf] vs -huge")
        values = -inf
        values(1) = -huge(1.0_real64)
        values(4) = -huge(1.0_real64)
        call assert_false(any_below(values, -inf), "test_any_below_real: [-huge, -Inf, -Inf, -huge] vs -Inf")
    end subroutine test_any_below_real

    !> `any_above` on real64: `.false.` for no elements, for an element equal to the bound and for
    !| `0.0` against `-0.0` (equal under IEEE); `.true.` at the first or only the last position.
    !| `Inf` is above every finite bound, and nothing is above `Inf`.
    subroutine test_any_above_real()
        real(real64), allocatable :: empty(:)
        real(real64) :: values(4), inf

        inf = ieee_value(1.0_real64, ieee_positive_inf)

        allocate (empty(0))
        call assert_false(any_above(empty, 0.0_real64), "test_any_above_real: no elements")
        values = [1.0_real64, 3.0_real64, 2.0_real64, 3.0_real64]
        call assert_false(any_above(values, 3.0_real64), "test_any_above_real: [1, 3, 2, 3] vs 3, equal is not above")
        values = [-1.0_real64, 0.0_real64, -2.0_real64, 0.0_real64]
        call assert_false(any_above(values, -0.0_real64), "test_any_above_real: [-1, 0, -2, 0] vs -0, 0 is not above -0")
        values = [3.5_real64, 3.0_real64, 2.0_real64, 3.0_real64]
        call assert_true(any_above(values, 3.0_real64), "test_any_above_real: [3.5, 3, 2, 3] vs 3, the first is above")
        values = [1.0_real64, 3.0_real64, 2.0_real64, 3.5_real64]
        call assert_true(any_above(values, 3.0_real64), "test_any_above_real: [1, 3, 2, 3.5] vs 3, only the last is above")

        values = 1.0_real64
        values(4) = inf
        call assert_true(any_above(values, huge(1.0_real64)), "test_any_above_real: [1, 1, 1, Inf] vs huge")
        values = inf
        values(1) = huge(1.0_real64)
        values(4) = huge(1.0_real64)
        call assert_false(any_above(values, inf), "test_any_above_real: [huge, Inf, Inf, huge] vs Inf")
    end subroutine test_any_above_real

    !> A rank-1 section gives the section's answer, not the whole array's: a strided row `m(i, :)`,
    !| a column `m(:, j)`, and every second element `v(1::2)`. The one deciding element of `m` is
    !| `m(2, 3)`, which lies in row 2 and column 3 only.
    subroutine test_array_checks_sections()
        real(real64) :: m(3, 4)
        integer(int32) :: v(5)

        m = 1.0_real64
        m(2, 3) = -1.0_real64
        call assert_false(any_below(m(1, :), 0.0_real64), "test_array_checks_sections: row 1 has nothing below 0")
        call assert_true(any_below(m(2, :), 0.0_real64), "test_array_checks_sections: row 2 holds -1")
        call assert_false(any_below(m(:, 1), 0.0_real64), "test_array_checks_sections: column 1 has nothing below 0")
        call assert_true(any_below(m(:, 3), 0.0_real64), "test_array_checks_sections: column 3 holds -1")
        call assert_true(all_equal(m(3, :), 1.0_real64), "test_array_checks_sections: row 3 is all 1")
        call assert_false(all_equal(m(2, :), 1.0_real64), "test_array_checks_sections: row 2 is not all 1")
        call assert_false(any_above(m(2, :), 1.0_real64), "test_array_checks_sections: row 2 has nothing above 1")

        v = [1_int32, 9_int32, 1_int32, 9_int32, 1_int32]
        call assert_true(all_equal(v(1::2), 1_int32), "test_array_checks_sections: v(1::2) is all 1")
        call assert_false(any_above(v(1::2), 1_int32), "test_array_checks_sections: v(1::2) has nothing above 1")
        call assert_true(any_above(v, 1_int32), "test_array_checks_sections: v itself holds 9")
        call assert_false(all_equal(v(2::2), 1_int32), "test_array_checks_sections: v(2::2) is all 9")
    end subroutine test_array_checks_sections

    !> The six rank-2 specifics: the deciding element only in the last column (so every column is
    !| visited), a non-contiguous block `m(1:2, 2:)` that ignores a deciding element outside it,
    !| `0 x n` and `n x 0` giving the identity of `all` / `any`, and keyword calls, which resolve
    !| to the rank-2 specific because both ranks name their dummies alike.
    subroutine test_array_checks_rank2()
        integer(int32) :: ints(3, 4)
        real(real64) :: reals(3, 4)
        integer(int32), allocatable :: no_rows_int(:, :), no_columns_int(:, :)
        real(real64), allocatable :: no_rows_real(:, :), no_columns_real(:, :)

        ! the deciding element only in the last column
        ints = 5_int32
        call assert_true(all_equal(ints, 5_int32), "test_array_checks_rank2: int, all 5")
        call assert_false(any_below(ints, 5_int32), "test_array_checks_rank2: int, all 5, nothing below 5")
        call assert_false(any_above(ints, 5_int32), "test_array_checks_rank2: int, all 5, nothing above 5")
        ints(3, 4) = 4_int32
        call assert_false(all_equal(ints, 5_int32), "test_array_checks_rank2: int, 4 in the last column")
        call assert_true(any_below(ints, 5_int32), "test_array_checks_rank2: int, 4 in the last column is below 5")
        ints(3, 4) = 6_int32
        call assert_true(any_above(ints, 5_int32), "test_array_checks_rank2: int, 6 in the last column is above 5")

        reals = 5.0_real64
        call assert_true(all_equal(reals, 5.0_real64), "test_array_checks_rank2: real, all 5")
        call assert_false(any_below(reals, 5.0_real64), "test_array_checks_rank2: real, all 5, nothing below 5")
        call assert_false(any_above(reals, 5.0_real64), "test_array_checks_rank2: real, all 5, nothing above 5")
        reals(3, 4) = 4.0_real64
        call assert_false(all_equal(reals, 5.0_real64), "test_array_checks_rank2: real, 4 in the last column")
        call assert_true(any_below(reals, 5.0_real64), "test_array_checks_rank2: real, 4 in the last column is below 5")
        reals(3, 4) = 6.0_real64
        call assert_true(any_above(reals, 5.0_real64), "test_array_checks_rank2: real, 6 in the last column is above 5")

        ! a block that is not contiguous: the deciding element at (3, 1) lies outside it, the one
        ! at (2, 3) inside
        ints = 5_int32
        ints(3, 1) = 0_int32
        call assert_true(all_equal(ints(1:2, 2:), 5_int32), "test_array_checks_rank2: int block, all 5")
        call assert_false(any_below(ints(1:2, 2:), 5_int32), "test_array_checks_rank2: int block, nothing below 5")
        ints(3, 1) = 9_int32
        call assert_false(any_above(ints(1:2, 2:), 5_int32), "test_array_checks_rank2: int block, nothing above 5")
        ints(2, 3) = 0_int32
        call assert_false(all_equal(ints(1:2, 2:), 5_int32), "test_array_checks_rank2: int block holds 0")
        call assert_true(any_below(ints(1:2, 2:), 5_int32), "test_array_checks_rank2: int block, 0 is below 5")
        ints(2, 3) = 9_int32
        call assert_true(any_above(ints(1:2, 2:), 5_int32), "test_array_checks_rank2: int block, 9 is above 5")

        reals = 5.0_real64
        reals(3, 1) = 0.0_real64
        call assert_true(all_equal(reals(1:2, 2:), 5.0_real64), "test_array_checks_rank2: real block, all 5")
        call assert_false(any_below(reals(1:2, 2:), 5.0_real64), "test_array_checks_rank2: real block, nothing below 5")
        reals(3, 1) = 9.0_real64
        call assert_false(any_above(reals(1:2, 2:), 5.0_real64), "test_array_checks_rank2: real block, nothing above 5")
        reals(2, 3) = 0.0_real64
        call assert_false(all_equal(reals(1:2, 2:), 5.0_real64), "test_array_checks_rank2: real block holds 0")
        call assert_true(any_below(reals(1:2, 2:), 5.0_real64), "test_array_checks_rank2: real block, 0 is below 5")
        reals(2, 3) = 9.0_real64
        call assert_true(any_above(reals(1:2, 2:), 5.0_real64), "test_array_checks_rank2: real block, 9 is above 5")

        ! no elements: the identity of `all` and of `any`
        allocate (no_rows_int(0, 4), no_columns_int(3, 0), no_rows_real(0, 4), no_columns_real(3, 0))
        call assert_true(all_equal(no_rows_int, 1_int32), "test_array_checks_rank2: int, 0 x 4, all_equal")
        call assert_true(all_equal(no_columns_int, 1_int32), "test_array_checks_rank2: int, 3 x 0, all_equal")
        call assert_false(any_below(no_rows_int, 1_int32), "test_array_checks_rank2: int, 0 x 4, any_below")
        call assert_false(any_below(no_columns_int, 1_int32), "test_array_checks_rank2: int, 3 x 0, any_below")
        call assert_false(any_above(no_rows_int, 1_int32), "test_array_checks_rank2: int, 0 x 4, any_above")
        call assert_false(any_above(no_columns_int, 1_int32), "test_array_checks_rank2: int, 3 x 0, any_above")
        call assert_true(all_equal(no_rows_real, 1.0_real64), "test_array_checks_rank2: real, 0 x 4, all_equal")
        call assert_true(all_equal(no_columns_real, 1.0_real64), "test_array_checks_rank2: real, 3 x 0, all_equal")
        call assert_false(any_below(no_rows_real, 1.0_real64), "test_array_checks_rank2: real, 0 x 4, any_below")
        call assert_false(any_below(no_columns_real, 1.0_real64), "test_array_checks_rank2: real, 3 x 0, any_below")
        call assert_false(any_above(no_rows_real, 1.0_real64), "test_array_checks_rank2: real, 0 x 4, any_above")
        call assert_false(any_above(no_columns_real, 1.0_real64), "test_array_checks_rank2: real, 3 x 0, any_above")

        ! keyword calls
        reals = 2.0_real64
        call assert_true(all_equal(array=reals, value=2.0_real64), "test_array_checks_rank2: all_equal(array=, value=)")
        call assert_false(any_below(array=reals, bound=2.0_real64), "test_array_checks_rank2: any_below(array=, bound=)")
        ints = 2_int32
        call assert_false(any_above(array=ints, bound=2_int32), "test_array_checks_rank2: any_above(array=, bound=)")
    end subroutine test_array_checks_rank2

end module mod_test_f42_array_checks
