!> General assertion utilities for Fortran unit testing.
!! Provides a set of reusable assertion subroutines for verifying
!! expected behavior in tests of any kind (numeric, string, array, etc).
module asserts
    use, intrinsic :: iso_fortran_env, only: error_unit, real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use test_suite, only: COLOR_RED, COLOR_CREAM, COLOR_ERROR, COLOR_RESET, COLOR_GREEN, COLOR_YELLOW, COLOR_LIGHT_GRAY
    use tox_errors, only: get_err_code, get_err_arg_pos, ERR_OK
    implicit none
    private
    public :: assert_err
    public :: assert_true, assert_false, assert_equal_int, assert_not_equal_int, assert_array_int_contains
    public :: assert_equal_real, assert_not_equal_real, assert_equal_array_int
    public :: assert_equal_array_real, assert_no_nan_real, assert_no_inf_real
    public :: assert_in_range_real, assert_in_range_int, assert_contains_int, assert_sorted_int
    public :: assert_sorted_real, assert_same_shape, assert_string_equal
    public :: assert_string_contains, assert_allclose_array_real
    public :: assert_sum_equal, assert_unique_int, assert_permutation
    public :: assert_equal_array_char, assert_equal_array_logical
    public :: assert_equal_complex, assert_not_equal_complex, assert_equal_array_complex

    public :: operator(//)
    interface operator(//)
        module procedure str_concat_int, int_concat_str
        module procedure str_concat_real, real_concat_str
        module procedure str_concat_complex, complex_concat_str
    end interface operator(//)

contains

    subroutine assertion_error(msg, additional_msg, got, expected, tol, at)
        character(*), intent(in) :: msg
        character(*), intent(in), optional :: additional_msg, got, expected, tol, at

        logical :: comma

        write (error_unit, "(A)", advance="no") COLOR_RED // "ASSERTION FAILED" // COLOR_CREAM // ": " // COLOR_ERROR // trim(msg)
        if (present(additional_msg)) then
            write (error_unit, "(A)", advance="no") COLOR_RED // " --- " // COLOR_ERROR // trim(additional_msg)
        end if
        write (error_unit, "(A)") COLOR_RESET

        comma = .false.

        if (present(got) .or. present(expected) .or. present(tol) .or. present(at)) then
            write (error_unit, "(A)") COLOR_LIGHT_GRAY // "    ("

            if (present(got)) then
                write (error_unit, "(A)", advance="no") COLOR_RED // "      got" // COLOR_CREAM // ": " // got
                comma = .true.
            end if

            if (present(expected)) then
                if (comma) write (error_unit, "(A)") COLOR_LIGHT_GRAY // ","
                write (error_unit, "(A)", advance="no") COLOR_GREEN // " expected" // COLOR_CREAM // ": " // expected
                comma = .true.
            end if

            if (present(tol)) then
                if (comma) write (error_unit, "(A)") COLOR_LIGHT_GRAY // ","
                write (error_unit, "(A)", advance="no") COLOR_YELLOW // "      tol" // COLOR_CREAM // ": " // tol
                comma = .true.
            end if

            if (present(at)) then
                if (comma) write (error_unit, "(A)") COLOR_LIGHT_GRAY // ","
                write (error_unit, "(A)", advance="no") COLOR_CREAM // "       at: " // at
                comma = .true.
            end if

            write (error_unit, "()")
            write (error_unit, "(A)") COLOR_LIGHT_GRAY // "    )" // COLOR_RESET
        end if

        stop 1
    end subroutine assertion_error

    !> Assert that two complex numbers are equal within a tolerance.
    subroutine assert_equal_complex(a, b, tol, msg)
        complex(real64), intent(in) :: a, b
        real(real64), intent(in) :: tol
        character(*), intent(in) :: msg
        if (abs(a - b) > tol) then
            call assertion_error(msg, got=""//a, expected=""//b, tol=""//tol)
        end if
    end subroutine

    !> Assert that two complex numbers are not equal within a tolerance.
    subroutine assert_not_equal_complex(a, b, tol, msg)
        complex(real64), intent(in) :: a, b
        real(real64), intent(in) :: tol
        character(*), intent(in) :: msg
        if (abs(a - b) <= tol) then
            call assertion_error(msg, "(should not be equal)")
        end if
    end subroutine

    !> Assert that two complex arrays are equal within a tolerance.
    subroutine assert_equal_array_complex(a, b, n, tol, msg)
        complex(real64), intent(in) :: a(n), b(n)
        integer(int32), intent(in) :: n
        real(real64), intent(in) :: tol
        character(*), intent(in) :: msg
        integer(int32) :: i, n_diff
        n_diff = count(abs(a - b) > tol)
        if (n_diff > 0) then
            do i = 1, n
                if (abs(a(i) - b(i)) > tol) exit
            end do
            call assertion_error(msg, additional_msg=n_diff // " of " // n // " elements differ", &
                got=""//a(i), expected=""//b(i), tol=""//tol, at=""//i)
        end if
    end subroutine

    !> Assert that two logical arrays are equal within a tolerance.
    subroutine assert_equal_array_logical(a, b, n, msg)
        logical, intent(in) :: a(n), b(n)
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        integer(int32) :: i, n_diff
        n_diff = count(a .neqv. b)
        if (n_diff > 0) then
            do i = 1, n
                if (a(i) .neqv. b(i)) exit
            end do
            call assertion_error(msg, additional_msg=n_diff // " of " // n // " elements differ", &
                got=logical_to_str(a(i)), expected=logical_to_str(b(i)), at=""//i)
        end if
    end subroutine

    !> Assert that a logical condition is true.
    subroutine assert_true(cond, msg)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        if (.not. cond) then
            call assertion_error(msg, got=".false.", expected=".true.")
        end if
    end subroutine

    !> Assert that a logical condition is false.
    subroutine assert_false(cond, msg)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        if (cond) then
            call assertion_error(msg, got=".true.", expected=".false.")
        end if
    end subroutine

    !> Assert that two integers are equal.
    subroutine assert_equal_int(a, b, msg)
        integer(int32), intent(in) :: a, b
        character(*), intent(in) :: msg
        if (a /= b) then
            call assertion_error(trim(msg), got=""//a, expected=""//b)
        end if
    end subroutine

    !> Assert an error code, comparing the code and -- when asked -- the argument it blames.
    !|
    !| `ierr` packs both, so a mismatch reported as a single number ("expected 202, got 80201")
    !| costs the reader an unpacking step every time. This reports them apart. Omitting
    !| `arg_pos` checks the code alone, which is what a caller that does not care which dummy
    !| was rejected wants; `ERR_OK` needs no position either.
    subroutine assert_err(ierr, expected_code, msg, arg_pos)
        integer(int32), intent(in) :: ierr, expected_code
        character(*), intent(in) :: msg
        integer(int32), intent(in), optional :: arg_pos

        if (get_err_code(ierr) /= expected_code) then
            call assertion_error(trim(msg), additional_msg="wrong error code", &
                                 got=""//get_err_code(ierr)//" at argument "//get_err_arg_pos(ierr), &
                                 expected=""//expected_code)
            return
        end if

        if (present(arg_pos)) then
            if (get_err_arg_pos(ierr) /= arg_pos) then
                call assertion_error(trim(msg), additional_msg="right error code, wrong argument", &
                                     got=""//get_err_arg_pos(ierr), expected=""//arg_pos)
            end if
        end if
    end subroutine

    !> Assert that two integers are not equal.
    subroutine assert_not_equal_int(a, b, msg)
        integer(int32), intent(in) :: a, b
        character(*), intent(in) :: msg
        if (a == b) then
            call assertion_error(msg, additional_msg="should not be equal", got=a // " and " // b)
        end if
    end subroutine

    !> Assert that two real numbers are equal within a tolerance.
    subroutine assert_equal_real(a, b, tol, msg)
        real(real64), intent(in) :: a, b, tol
        character(*), intent(in) :: msg
        if (abs(a - b) > tol) then
            call assertion_error(msg, got=""//a, expected=""//b, tol=""//tol)
        end if
    end subroutine

    !> Assert that two real numbers are not equal within a tolerance.
    subroutine assert_not_equal_real(a, b, tol, msg)
        real(real64), intent(in) :: a, b, tol
        character(*), intent(in) :: msg
        if (abs(a - b) <= tol) then
            call assertion_error(msg, additional_msg="should not be equal", got=a // " and " // b)
        end if
    end subroutine

    !> Assert that two integer arrays are equal.
    subroutine assert_equal_array_int(a, b, n, msg)
        integer(int32), intent(in) :: a(n), b(n), n
        character(*), intent(in) :: msg
        integer(int32) :: i, n_diff
        n_diff = count(a /= b)
        if (n_diff > 0) then
            do i = 1, n
                if (a(i) /= b(i)) exit
            end do
            call assertion_error(msg, additional_msg=n_diff // " of " // n // " elements differ", &
                got=""//a(i), expected=""//b(i), at=""//i)
        end if
    end subroutine

    !> Assert that two real arrays are equal within a tolerance.
    subroutine assert_equal_array_real(a, b, n, tol, msg)
        real(real64), intent(in) :: a(n), b(n), tol
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        integer(int32) :: i, n_diff
        n_diff = count(abs(a - b) > tol)
        if (n_diff > 0) then
            do i = 1, n
                if (abs(a(i) - b(i)) > tol) exit
            end do
            call assertion_error(msg, additional_msg=n_diff // " of " // n // " elements differ", &
                got=""//a(i), expected=""//b(i), tol=""//tol, at=""//i)
        end if
    end subroutine

    !> Asserts that two character arrays are equal
    subroutine assert_equal_array_char(a, b, clen, n, msg)
        integer(int32), INTENT(IN) :: clen
        character(len=clen), intent(in) :: a(n), b(n)
        character(*), intent(in) :: msg
        integer, intent(in) :: n
        integer(int32) :: i, n_diff
        n_diff = count(a /= b)
        if (n_diff > 0) then
            do i = 1, n
                if (a(i) /= b(i)) exit
            end do
            call assertion_error(msg, additional_msg=n_diff // " of " // n // " elements differ", &
                got="'" // trim(a(i)) // "'", expected="'" // trim(b(i)) // "'", at=""//i)
        end if
    end subroutine

    !> Assert that a real array contains no NaN values.
    subroutine assert_no_nan_real(a, n, msg)
        real(real64), intent(in) :: a(n)
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        integer :: i

        do i = 1, n
            if (ieee_is_nan(a(i))) then
                call assertion_error(msg, additional_msg="NaN detected", got=""//a(i), at=""//i)
            end if
        end do
    end subroutine

    !> Assert that a real array contains no Inf values.
    subroutine assert_no_inf_real(a, n, msg)
        real(real64), intent(in) :: a(n)
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        integer :: i
        do i = 1, n
            if (abs(a(i)) > huge(1.0_real64)) then
                call assertion_error(msg, additional_msg="Inf detected", got=""//a(i), at=""//i)
            end if
        end do
    end subroutine

    !> Assert that a real value is within a given range [minval, maxval].
    subroutine assert_in_range_real(a, minval, maxval, msg)
        real(real64), intent(in) :: a, minval, maxval
        character(*), intent(in) :: msg
        if (a < minval .or. a > maxval) then
            call assertion_error(msg, expected="value in range [" // minval // "," // maxval // "]", got=""//a)
        end if
    end subroutine

    !> Assert that an integer value is within a given range [minval, maxval].
    subroutine assert_in_range_int(a, minval, maxval, msg)
        integer(int32), intent(in) :: a, minval, maxval
        character(*), intent(in) :: msg
        if (a < minval .or. a > maxval) then
            call assertion_error(msg, expected="value in range [" // minval // "," // maxval // "]", got=""//a)
        end if
    end subroutine

    !> Assert that an integer array contains a given value.
    subroutine assert_contains_int(arr, n, val, msg)
        integer(int32), intent(in) :: arr(n), n, val
        character(*), intent(in) :: msg
        if (.not. any(arr == val)) then
            call assertion_error(msg, additional_msg="value not found in array", &
                expected=""//val, got=array_preview_int(arr, n))
        end if
    end subroutine

    !> Assert that an integer array is sorted in non-decreasing order.
    subroutine assert_sorted_int(arr, n, msg)
        integer(int32), intent(in) :: arr(n), n
        character(*), intent(in) :: msg
        integer :: i
        do i = 2, n
            if (arr(i) < arr(i - 1)) then
                call assertion_error(msg, additional_msg="not sorted", &
                    got=arr(i - 1) // " > " // arr(i), at=""//i)
            end if
        end do
    end subroutine

    !> Assert that a real array is sorted in non-decreasing order.
    subroutine assert_sorted_real(arr, n, msg)
        real(real64), intent(in) :: arr(n)
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        integer :: i
        do i = 2, n
            if (arr(i) < arr(i - 1)) then
                call assertion_error(msg, additional_msg="not sorted", &
                    got=arr(i - 1) // " > " // arr(i), at=""//i)
            end if
        end do
    end subroutine

    !> Assert that two arrays have the same shape (1D only).
    subroutine assert_same_shape(n1, n2, msg)
        integer(int32), intent(in) :: n1, n2
        character(*), intent(in) :: msg
        if (n1 /= n2) then
            call assertion_error(msg, additional_msg="shapes differ", got=n1 // " and " // n2)
        end if
    end subroutine

    !> Assert that two strings are equal.
    subroutine assert_string_equal(a, b, msg)
        character(*), intent(in) :: a, b, msg
        if (trim(a) /= trim(b)) then
            call assertion_error(msg, got="'" // trim(a) // "'", expected="'" // trim(b) // "'")
        end if
    end subroutine

    !> Assert that string a contains string b.
    subroutine assert_string_contains(a, b, msg)
        character(*), intent(in) :: a, b, msg
        if (index(a, b) == 0) then
            call assertion_error(msg, additional_msg="substring not found", &
                got="string '" // trim(a) // "'", expected="substring '" // trim(b) // "'")
        end if
    end subroutine

    !> Assert that array a contains value b.
    subroutine assert_array_int_contains(a, b, n, msg)
        integer(int32), intent(in) :: n, a(n), b
        character(*), intent(in) :: msg
        if (findloc(a, b, 1) == 0) then
            call assertion_error(msg, additional_msg="value not found in array", &
                expected=""//b, got=array_preview_int(a, n))
        end if
    end subroutine

    !> Assert that two real arrays are close within relative and absolute tolerance.
    subroutine assert_allclose_array_real(a, b, n, rtol, atol, msg)
        real(real64), intent(in) :: a(n), b(n), rtol, atol
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        integer :: i
        real(real64) :: thresh
        do i = 1, n
            thresh = atol + rtol*abs(b(i))
            if (abs(a(i) - b(i)) > thresh) then
                call assertion_error(msg, additional_msg="|got - expected| = " // abs(a(i) - b(i)) // &
                    " exceeds atol + rtol*|expected| = " // thresh, got=""//a(i), expected=""//b(i), at=""//i)
            end if
        end do
    end subroutine

    !> Assert that the sum of an array equals an expected value.
    subroutine assert_sum_equal(arr, n, expected, msg)
        real(real64), intent(in) :: arr(n), expected
        integer(int32), intent(in) :: n
        character(*), intent(in) :: msg
        real(real64) :: s
        s = sum(arr)
        if (abs(s - expected) > 1e-12_real64) then
            call assertion_error(msg, got="sum=" // s, expected=""//expected)
        end if
    end subroutine

    !> Assert that all elements in an integer array are unique.
    subroutine assert_unique_int(arr, n, msg)
        integer(int32), intent(in) :: arr(n), n
        character(*), intent(in) :: msg
        integer :: i, j
        do i = 1, n - 1
            do j = i + 1, n
                if (arr(i) == arr(j)) then
                    call assertion_error(msg, additional_msg="duplicate value", got=""//arr(i), at=i // " and " // j)
                end if
            end do
        end do
    end subroutine

    !> Assert that an integer array is a permutation of 1..n.
    subroutine assert_permutation(arr, n, msg)
        integer(int32), intent(in) :: arr(n), n
        character(*), intent(in) :: msg
        integer :: i
        logical :: found(n)
        found = .false.
        do i = 1, n
            if (arr(i) < 1 .or. arr(i) > n) then
                call assertion_error(msg, additional_msg="value out of range for permutation", &
                    got=""//arr(i), expected="value in [1, " // n // "]", at=""//i)
            end if
            if (found(arr(i))) then
                call assertion_error(msg, additional_msg="duplicate value", got=""//arr(i), at=""//i)
            end if
            found(arr(i)) = .true.
        end do
    end subroutine

    !> Renders a logical as ".true."/".false." for use in assertion diagnostics.
    pure function logical_to_str(l) result(str_out)
        logical, intent(in) :: l
        character(len=:), allocatable :: str_out
        if (l) then
            str_out = ".true."
        else
            str_out = ".false."
        end if
    end function logical_to_str

    !> Renders an integer array as "[a, b, c, ...]" for assertion diagnostics, truncating
    !! long arrays so a failing "contains" assertion doesn't flood the terminal with a
    !! wall of numbers that isn't any more informative than the first few elements.
    function array_preview_int(arr, n) result(preview)
        integer(int32), intent(in) :: n, arr(n)
        character(len=:), allocatable :: preview
        integer(int32), parameter :: max_shown = 15
        integer(int32) :: i, n_shown

        n_shown = min(n, max_shown)
        preview = "["
        do i = 1, n_shown
            if (i > 1) preview = preview // ", "
            preview = preview // arr(i)
        end do
        if (n > max_shown) preview = preview // ", ... (" // (n - max_shown) // " more)"
        preview = preview // "]"
    end function array_preview_int

    !> Very efficient digit counting function for an integer. Having the absolute value, it needs only 4 cycles and one assignment in worst case
    pure integer(int32) function digit_count_int32(val) result(digit_count)
        integer(int32), intent(in) :: val
        integer(int32) :: abs_val
        integer(int32), parameter :: min_int = -huge(1_int32)-1, max_int=huge(1_int32)

        ! absolute value, safe for -2147483648
        if (val == min_int) then
            abs_val = max_int
        else
            abs_val = abs(val)
        end if

        ! --- perfectly balanced binary search over powers of 10 ---
        if (abs_val < 100000_int32) then
            if (abs_val < 1000_int32) then
                if (abs_val < 100_int32) then
                    if (abs_val < 10_int32) then
                        digit_count = 1
                    else
                        digit_count = 2
                    end if
                else
                    if (abs_val < 1000_int32) then
                        digit_count = 3
                    else
                        digit_count = 4
                    end if
                end if
            else
                if (abs_val < 10000_int32) then
                    digit_count = 4
                else
                    digit_count = 5
                end if
            end if
        else
            if (abs_val < 10000000_int32) then
                if (abs_val < 1000000_int32) then
                    digit_count = 6
                else
                    digit_count = 7
                end if
            else
                if (abs_val < 100000000_int32) then
                    digit_count = 8
                else
                    if (abs_val < 1000000000_int32) then
                        digit_count = 9
                    else
                        digit_count = 10
                    end if
                end if
            end if
        end if
    end function digit_count_int32

    function str_concat_complex(str, complex) result(str_out)
        character(len=*), intent(in) :: str
        complex(real64), intent(in) :: complex
        character(len=:), allocatable :: str_out

        str_out = str // "("//real(complex)//", "//aimag(complex)//")"
    end function str_concat_complex

    function complex_concat_str(complex, str) result(str_out)
        character(len=*), intent(in) :: str
        complex(real64), intent(in) :: complex
        character(len=:), allocatable :: str_out

        str_out = "" // complex // str
    end function complex_concat_str

    function str_concat_real(str, real) result(str_out)
        character(len=*), intent(in) :: str
        real(real64), intent(in) :: real
        character(len=:), allocatable :: str_out

        character(32) :: tmp

        write(tmp, '(ES25.17E3)') real
        str_out = str // trim(adjustl(tmp))
    end function str_concat_real

    function real_concat_str(real, str) result(str_out)
        character(len=*), intent(in) :: str
        real(real64), intent(in) :: real
        character(len=:), allocatable :: str_out

        str_out = "" // real // str
    end function real_concat_str

    function str_concat_int(str, int) result(str_out)
        character(len=*), intent(in) :: str
        integer(int32), intent(in) :: int
        character(len=:), allocatable :: str_out

        if (int < 0) then
            ! one more char for sign
            allocate(character(len=len(str) + digit_count_int32(int) + 1) :: str_out)
        else
            allocate(character(len=len(str) + digit_count_int32(int)) :: str_out)
        end if

        str_out(1:len(str)) = str
        write (str_out(len(str)+1:), "(I0)") int
    end function str_concat_int

    function int_concat_str(int, str) result(str_out)
        character(len=*), intent(in) :: str
        integer(int32), intent(in) :: int
        character(len=:), allocatable :: str_out

        str_out = "" // int // str
    end function int_concat_str
end module asserts
