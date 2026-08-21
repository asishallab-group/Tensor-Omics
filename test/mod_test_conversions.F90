!> Unit test for conversion routines.
module mod_test_conversions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use iso_c_binding
    use tox_conversions
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 0

contains

    !> Get array of all available tests.
    function get_all_tests_conversions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(5))

        all_tests(1) = test_case("test_tox_conversions_c_char_as_char", test_c_char_as_char)
        all_tests(2) = test_case("test_tox_conversions_c_char_as_view", test_c_char_as_view)
        all_tests(3) = test_case("test_tox_conversions_c_int_as_logical", test_c_int_as_logical)
        all_tests(4) = test_case("test_tox_conversions_char_as_c_char", test_char_as_c_char)
        all_tests(5) = test_case("test_tox_conversions_logical_as_c_int", test_logical_as_c_int)
    end function get_all_tests_conversions

    !> Test C int to logical conversion for scalar and array cases.
    subroutine test_c_int_as_logical
        integer(c_int) :: c_val
        logical :: casted_fortran
        integer(c_int) :: c_val_array(1)
        logical :: casted_array(1)

        c_val = 42_c_int
        call c_int_as_logical(c_val, casted_fortran)
        call assert_true(casted_fortran, "test_conversions_c_int_as_logical: value mismatch")

        c_val = -42_c_int
        call c_int_as_logical(c_val, casted_fortran)
        call assert_true(casted_fortran, "test_conversions_c_int_as_logical: value mismatch")

        c_val = 0_c_int
        call c_int_as_logical(c_val, casted_fortran)
        call assert_false(casted_fortran, "test_conversions_c_int_as_logical: value mismatch")

        c_val_array = [c_val]
        call c_int_as_logical(c_val_array, casted_array)
        call assert_false(casted_array(1), "test_conversions_c_int_as_logical: value mismatch")
    end subroutine test_c_int_as_logical

    !> Test C char to Fortran char conversion for scalar and array cases.
    subroutine test_c_char_as_char
        character(kind=c_char, len=1) :: c_val
        character(len=1) :: casted_fortran
        character(len=1) :: casted_array(1)
        character(kind=c_char, len=1) :: c_val_array(1)
        character(len=1) :: expected_fortran

        expected_fortran = "H"
        c_val = expected_fortran
        call c_char_as_char(c_val, casted_fortran)
        call assert_true(casted_fortran == expected_fortran, "test_conversions_c_char_as_char: value mismatch for fortran")

        c_val_array = [c_val]
        call c_char_as_char(c_val_array, casted_array)
        call assert_true(casted_array(1) == expected_fortran, "test_conversions_c_char_as_char: value mismatch for char array")

        expected_fortran = achar(0)
        call c_char_as_char(c_null_char, casted_fortran)
        call assert_true(casted_fortran == expected_fortran, "test_conversions_c_char_as_char: value mismatch for null char")
    end subroutine test_c_char_as_char

    !> A C buffer read as one Fortran string, without copying it.
    !| The two padding conventions are the point: a hand-written C caller NUL-terminates, the
    !| generated bindings blank-pad, and both have to reach the same answer. Only the second
    !| is exercised by the suites end to end, so the first is checked here or nowhere.
    subroutine test_c_char_as_view
        character(kind=c_char, len=1), target :: c_char_array(8)
        character(len=:), pointer :: view

        ! NUL-terminated, C's convention: the view stops at the NUL
        c_char_array = ["w", "a", "r", "d", c_null_char, c_null_char, c_null_char, c_null_char]
        view => c_char_as_view(c_char_array)
        call assert_equal_int(len(view), 4, "c_char_as_view: NUL-terminated length")
        call assert_true(view == "ward", "c_char_as_view: NUL-terminated value")

        ! blank-padded, what the generated bindings send: no NUL, so the whole buffer
        c_char_array = ["w", "a", "r", "d", " ", " ", " ", " "]
        view => c_char_as_view(c_char_array)
        call assert_equal_int(len(view), 8, "c_char_as_view: blank-padded length")
        ! and it still compares equal to the short literal, because Fortran blank-pads the
        ! shorter operand -- which is what lets `select case` accept either convention
        call assert_true(view == "ward", "c_char_as_view: blank-padded compares equal")

        ! no padding at all
        c_char_array = ["w", "e", "i", "g", "h", "t", "e", "d"]
        view => c_char_as_view(c_char_array)
        call assert_equal_int(len(view), 8, "c_char_as_view: full width length")
        call assert_true(view == "weighted", "c_char_as_view: full width value")

        ! a leading NUL is the empty string
        c_char_array = [c_null_char, "e", "l", "l", "o", " ", " ", " "]
        view => c_char_as_view(c_char_array)
        call assert_equal_int(len(view), 0, "c_char_as_view: leading NUL gives the empty string")

        ! it is a view, not a copy
        c_char_array = ["w", "a", "r", "d", c_null_char, c_null_char, c_null_char, c_null_char]
        view => c_char_as_view(c_char_array)
        c_char_array(1) = "h"
        call assert_true(view == "hard", "c_char_as_view: writing through the buffer changes the view")
    end subroutine test_c_char_as_view

    !> Test logical to C int conversion for scalar and array cases.
    subroutine test_logical_as_c_int
        logical :: f_val
        integer(c_int) :: casted_c
        logical :: f_val_array(1)
        integer(c_int) :: casted_array(1)

        f_val = .true.
        call logical_as_c_int(f_val, casted_c)
        call assert_true(casted_c /= 0_c_int, "test_tox_conversions_logical_as_c_int: expected non-zero for .true.")

        f_val = .false.
        call logical_as_c_int(f_val, casted_c)
        call assert_true(casted_c == 0_c_int, "test_tox_conversions_logical_as_c_int: expected zero for .false.")

        f_val_array = [.false.]
        call logical_as_c_int(f_val_array, casted_array)
        call assert_true(casted_array(1) == 0_c_int, "test_tox_conversions_logical_as_c_int: array value mismatch")
    end subroutine test_logical_as_c_int

    !> Test Fortran char to C char conversion for scalar and array cases.
    subroutine test_char_as_c_char
        character(len=1) :: f_val
        character(kind=c_char, len=1) :: casted_c
        character(len=1) :: f_val_array(1)
        character(kind=c_char, len=1) :: casted_array(1)

        f_val = "H"
        call char_as_c_char(f_val, casted_c)
        call assert_true(casted_c == "H", "test_tox_conversions_char_as_c_char: value mismatch")

        f_val_array = [f_val]
        call char_as_c_char(f_val_array, casted_array)
        call assert_true(casted_array(1) == "H", "test_tox_conversions_char_as_c_char: array value mismatch")

        f_val = achar(0)
        call char_as_c_char(f_val, casted_c)
        call assert_true(casted_c == c_null_char, "test_tox_conversions_char_as_c_char: null char mismatch")
    end subroutine test_char_as_c_char

end module mod_test_conversions
