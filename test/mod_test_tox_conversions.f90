!> Unit test suite for tox_conversions routine.
module mod_test_tox_conversions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use iso_c_binding
    use tox_conversions
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=64) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    integer(int32), parameter :: TEST_COUNT = 14
    real(real64), parameter :: TOL = 0

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(TEST_COUNT)

        all_tests(1) = test_case("test_tox_conversions_c_double_as_real64", test_tox_conversions_c_double_as_real64)
        all_tests(2) = test_case("test_tox_conversions_c_int_as_int32", test_tox_conversions_c_int_as_int32)
        all_tests(3) = test_case("test_tox_conversions_c_char_as_char", test_tox_conversions_c_char_as_char)
        all_tests(4) = test_case("test_tox_conversions_c_char_1d_as_string", test_tox_conversions_c_char_1d_as_string)
        all_tests(5) = test_case("test_tox_conversions_c_char_2d_as_string", test_tox_conversions_c_char_2d_as_string)
        all_tests(6) = test_case("test_tox_conversions_c_complex_as_complex", test_tox_conversions_c_complex_as_complex)
        all_tests(7) = test_case("test_tox_conversions_c_int_as_logical", test_tox_conversions_c_int_as_logical)
        all_tests(8) = test_case("test_tox_conversions_real64_as_c_double", test_tox_conversions_real64_as_c_double)
        all_tests(9) = test_case("test_tox_conversions_int32_as_c_int", test_tox_conversions_int32_as_c_int)
        all_tests(10) = test_case("test_tox_conversions_char_as_c_char", test_tox_conversions_char_as_c_char)
        all_tests(11) = test_case("test_tox_conversions_string_as_c_char_1d", test_tox_conversions_string_as_c_char_1d)
        all_tests(12) = test_case("test_tox_conversions_string_as_c_char_2d", test_tox_conversions_string_as_c_char_2d)
        all_tests(13) = test_case("test_tox_conversions_complex_as_c_complex", test_tox_conversions_complex_as_c_complex)
        all_tests(14) = test_case("test_tox_conversions_logical_as_c_int", test_tox_conversions_logical_as_c_int)
    end function get_all_tests

    subroutine test_tox_conversions_c_double_as_real64
        real(c_double) :: c_val
        real(real64) :: casted_fortran
        real(real64) :: expected_fortran
        real(c_double) :: c_val_array(1)
        real(real64) :: casted_array(1)

        c_val = 3.141592653589793_c_double
        expected_fortran = 3.141592653589793_real64
        call c_double_as_real64(c_val, casted_fortran)
        call assert_equal_real(casted_fortran, expected_fortran, TOL, "test_tox_conversions_c_double_as_real64: value mismatch")

        c_val_array = [c_val]
        call c_double_as_real64(c_val_array, casted_array)
        call assert_equal_array_real(casted_array, [expected_fortran], 1, TOL, "test_tox_conversions_c_double_as_real64: value mismatch")
    end subroutine test_tox_conversions_c_double_as_real64

    subroutine test_tox_conversions_c_complex_as_complex
        complex(c_double_complex) :: c_val
        complex(real64) :: casted_fortran
        complex(real64) :: expected_fortran
        complex(c_double_complex) :: c_val_array(1)
        complex(real64) :: casted_array(1)

        c_val = (3.141592653589793_c_double_complex, 3.141592653589793_c_double_complex)
        expected_fortran = (3.141592653589793_real64, 3.141592653589793_real64)
        call c_complex_as_complex(c_val, casted_fortran)
        call assert_equal_complex(casted_fortran, expected_fortran, TOL, "test_tox_conversions_c_complex_as_complex: value mismatch")

        c_val_array = [c_val]
        call c_complex_as_complex(c_val_array, casted_array)
        call assert_equal_array_complex(casted_array, [expected_fortran], 1, TOL, "test_tox_conversions_c_complex_as_complex: value mismatch")
    end subroutine test_tox_conversions_c_complex_as_complex

    subroutine test_tox_conversions_c_int_as_int32
        integer(c_int) :: c_val
        integer(int32) :: casted_fortran
        integer(c_int) :: c_val_array(1)
        integer(int32) :: casted_array(1)
        integer(int32) :: expected_fortran

        c_val = 42_c_int
        expected_fortran = 42_int32
        call c_int_as_int32(c_val, casted_fortran)
        call assert_equal_int(casted_fortran, expected_fortran, "test_tox_conversions_c_int_as_int32: value mismatch")

        c_val_array = [c_val]
        call c_int_as_int32(c_val_array, casted_array)
        call assert_equal_array_int(casted_array, [expected_fortran], 1, "test_tox_conversions_c_int_as_int32: value mismatch")
    end subroutine test_tox_conversions_c_int_as_int32

    subroutine test_tox_conversions_c_int_as_logical
        integer(c_int) :: c_val
        logical :: casted_fortran
        integer(c_int) :: c_val_array(1)
        logical :: casted_array(1)

        c_val = 42_c_int
        call c_int_as_logical(c_val, casted_fortran)
        call assert_true(casted_fortran, "test_tox_conversions_c_int_as_logical: value mismatch")

        c_val = -42_c_int
        call c_int_as_logical(c_val, casted_fortran)
        call assert_true(casted_fortran, "test_tox_conversions_c_int_as_logical: value mismatch")

        c_val = 0_c_int
        call c_int_as_logical(c_val, casted_fortran)
        call assert_false(casted_fortran, "test_tox_conversions_c_int_as_logical: value mismatch")

        c_val_array = [c_val]
        call c_int_as_logical(c_val_array, casted_array)
        call assert_false(casted_array(1), "test_tox_conversions_c_int_as_logical: value mismatch")
    end subroutine test_tox_conversions_c_int_as_logical

    subroutine test_tox_conversions_c_char_as_char
        character(c_char) :: c_val
        character(len=1) :: casted_fortran
        character(len=1) :: casted_array(1)
        character(c_char) :: c_val_array(1)
        character(len=1) :: expected_fortran

        expected_fortran = "H"
        c_val = expected_fortran
        call c_char_as_char(c_val, casted_fortran)
        call assert_true(casted_fortran == expected_fortran, "test_tox_conversions_c_char_as_char: value mismatch for fortran")

        c_val_array = [c_val]
        call c_char_as_char(c_val_array, casted_array)
        call assert_true(casted_array(1) == expected_fortran, "test_tox_conversions_c_char_as_char: value mismatch for char array")

        expected_fortran = achar(0)
        call c_char_as_char(c_null_char, casted_fortran)
        call assert_true(casted_fortran == expected_fortran, "test_tox_conversions_c_char_as_char: value mismatch for null char")
    end subroutine test_tox_conversions_c_char_as_char

    subroutine test_tox_conversions_c_char_1d_as_string
        character(c_char) :: c_char_array(5)
        character(len=:), allocatable :: f_char
        integer(int32) :: ierr

        call set_ok(ierr)

        c_char_array = ["H", "e", "l", "l", "o"]
        call c_char_1d_as_string(c_char_array, f_char, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_tox_conversions_c_char_1d_as_string: Unexpected Error Code")
        call assert_true(f_char == "Hello", "test_tox_conversions_c_char_1d_as_string: value mismatch")

        ! test empty string
        c_char_array = [c_null_char, "e", "l", "l", "o"]
        call c_char_1d_as_string(c_char_array, f_char, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_tox_conversions_c_char_1d_as_string: Unexpected Error Code")
        call assert_true(f_char == "", "test_tox_conversions_c_char_1d_as_string: value mismatch")
    end subroutine test_tox_conversions_c_char_1d_as_string

    subroutine test_tox_conversions_c_char_2d_as_string
        character(c_char) :: c_char_array(5, 2)
        character(len=:), allocatable :: f_char(:)
        integer(int32) :: ierr

        call set_ok(ierr)

        c_char_array(:, 1) = ["H", "e", "l", "l", "o"]
        c_char_array(:, 2) = [c_null_char, "e", "l", "l", "o"]

        call c_char_2d_as_string(c_char_array, f_char, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_tox_conversions_c_char_1d_as_string: Unexpected Error Code")
        call assert_equal_int(size(f_char, 1), 2, "test_tox_conversions_c_char_2d_as_string: Did not get two strings")
        call assert_true(f_char(1) == "Hello" .and. f_char(2) == "", "test_tox_conversions_c_char_2d_as_string: value mismatch")
    end subroutine test_tox_conversions_c_char_2d_as_string

    subroutine test_tox_conversions_real64_as_c_double
        real(real64) :: f_val
        real(c_double) :: casted_c
        real(c_double) :: expected_c
        real(real64) :: f_val_array(1)
        real(c_double) :: casted_array(1)

        f_val = 3.141592653589793_real64
        expected_c = 3.141592653589793_c_double
        call real64_as_c_double(f_val, casted_c)
        call assert_equal_real(casted_c, expected_c, TOL, "test_tox_conversions_real64_as_c_double: value mismatch")

        f_val_array = [f_val]
        call real64_as_c_double(f_val_array, casted_array)
        call assert_equal_array_real(casted_array, [expected_c], 1, TOL, "test_tox_conversions_real64_as_c_double: value mismatch")
    end subroutine test_tox_conversions_real64_as_c_double

    subroutine test_tox_conversions_complex_as_c_complex
        complex(real64) :: f_val
        complex(c_double_complex) :: casted_c
        complex(c_double_complex) :: expected_c
        complex(real64) :: f_val_array(1)
        complex(c_double_complex) :: casted_array(1)

        f_val = (3.141592653589793_real64, 3.141592653589793_real64)
        expected_c = (3.141592653589793_c_double_complex, 3.141592653589793_c_double_complex)
        call complex_as_c_complex(f_val, casted_c)
        call assert_equal_complex(casted_c, expected_c, TOL, "test_tox_conversions_complex_as_c_complex: value mismatch")

        f_val_array = [f_val]
        call complex_as_c_complex(f_val_array, casted_array)
        call assert_equal_array_complex(casted_array, [expected_c], 1, TOL, "test_tox_conversions_complex_as_c_complex: value mismatch")
    end subroutine test_tox_conversions_complex_as_c_complex

    subroutine test_tox_conversions_int32_as_c_int
        integer(int32) :: f_val
        integer(c_int) :: casted_c
        integer(c_int) :: expected_c
        integer(int32) :: f_val_array(1)
        integer(c_int) :: casted_array(1)

        f_val = 42_int32
        expected_c = 42_c_int
        call int32_as_c_int(f_val, casted_c)
        call assert_equal_int(casted_c, expected_c, "test_tox_conversions_int32_as_c_int: value mismatch")

        f_val_array = [f_val]
        call int32_as_c_int(f_val_array, casted_array)
        call assert_equal_array_int(casted_array, [expected_c], 1, "test_tox_conversions_int32_as_c_int: value mismatch")
    end subroutine test_tox_conversions_int32_as_c_int

    subroutine test_tox_conversions_logical_as_c_int
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
    end subroutine test_tox_conversions_logical_as_c_int

    subroutine test_tox_conversions_char_as_c_char
        character(len=1) :: f_val
        character(c_char) :: casted_c
        character(len=1) :: f_val_array(1)
        character(c_char) :: casted_array(1)

        f_val = "H"
        call char_as_c_char(f_val, casted_c)
        call assert_true(casted_c == "H", "test_tox_conversions_char_as_c_char: value mismatch")

        f_val_array = [f_val]
        call char_as_c_char(f_val_array, casted_array)
        call assert_true(casted_array(1) == "H", "test_tox_conversions_char_as_c_char: array value mismatch")

        f_val = achar(0)
        call char_as_c_char(f_val, casted_c)
        call assert_true(casted_c == c_null_char, "test_tox_conversions_char_as_c_char: null char mismatch")
    end subroutine test_tox_conversions_char_as_c_char

    subroutine test_tox_conversions_string_as_c_char_1d
        character(len=5) :: f_str
        character(c_char), allocatable :: c_array(:)

        f_str = "Hello"
        allocate(c_array(6))
        call string_as_c_char_1d(f_str, c_array)
        call assert_true(all(c_array == ["H", "e", "l", "l", "o", c_null_char]), "test_tox_conversions_string_as_c_char_1d: value mismatch")

        call string_as_c_char_1d(f_str, c_array(1:5))
        call assert_true(all(c_array(1:5) == ["H", "e", "l", "l", c_null_char]), "test_tox_conversions_string_as_c_char_1d: value mismatch")

        f_str = ""
        call string_as_c_char_1d(f_str, c_array)
        call assert_true(c_array(1) == c_null_char, "test_tox_conversions_string_as_c_char_1d: empty string mismatch")

        ! just check empty array output
        deallocate(c_array)
        allocate(c_array(0))
        call string_as_c_char_1d(f_str, c_array)
    end subroutine test_tox_conversions_string_as_c_char_1d

    subroutine test_tox_conversions_string_as_c_char_2d
        character(len=5) :: f_str(2)
        character(c_char) :: c_array(6, 2)

        f_str(1) = "Hello"
        f_str(2) = ""
        call string_as_c_char_2d(f_str, c_array)
        call assert_true(all(c_array(:, 1) == ["H", "e", "l", "l", "o", c_null_char]), "test_tox_conversions_string_as_c_char_2d: mismatch in first string")
        call assert_true(c_array(1, 2) == c_null_char, "test_tox_conversions_string_as_c_char_2d: mismatch in second string")
    end subroutine test_tox_conversions_string_as_c_char_2d

    !> Run all tox_conversions tests.
    subroutine run_all_tests_tox_conversions
        type(test_case) :: all_tests(TEST_COUNT)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print *, trim(all_tests(i)%name), " passed."
        end do
        print *, "All tox_conversions tests passed successfully."
    end subroutine run_all_tests_tox_conversions

    !> Run specific tox_conversions tests by name.
    subroutine run_named_tests_tox_conversions(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case) :: all_tests(TEST_COUNT)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print *, trim(test_names(i)), " passed."
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_conversions
end module mod_test_tox_conversions
