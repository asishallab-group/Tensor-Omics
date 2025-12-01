!> Unit test suite for tox_json routine.
module mod_test_json
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tox_json
    use tox_errors
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(1)

        all_tests(1) = test_case("test_tox_json_serialization", test_serialization)
    end function get_all_tests

    subroutine test_serialization
        ! Test idea:
        !   - an array of all possible types
        !   - for simplicity the array type is a self reference and tests recursion simultaneously
        !   - also, the self-containing array is reused as values of the object, and the object is also included in the array
        ! 
        type(json_value), dimension(6), target :: elements
        type(json_array), target :: array
        type(json_object), target :: object
        character(len=7), dimension(6), target :: keys
        integer(int32) :: unit, ierr
        integer(int32), target:: test_integer
        real(real64), target:: test_real
        complex(real64), target:: test_complex
        logical, target:: test_logical
        character(len=:), allocatable :: expected_fragment
        character(len=:), allocatable :: actual_fragment

        keys(1) = "integer"
        test_integer = -huge(1_int32)
        elements(1)%value => test_integer

        keys(2) = "real"
        test_real = huge(1.0_real64)
        elements(2)%value => test_real

        keys(3) = "logical"
        test_logical = .true.
        elements(3)%value => test_logical

        keys(4) = "complex"
        test_complex = cmplx(1.0_real64, -1.0_real64)
        elements(4)%value => test_complex

        keys(5) = "array"
        array%array => elements
        elements(5)%value => array

        keys(6) = "object"
        object%keys => keys
        object%values => elements
        elements(6)%value => object

        ! Case 1: basic object serialization test with recursion
        open(newunit=unit, file="test_json.json", form='formatted', access='stream', status='replace', iostat=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_serialization: could not open file")


        call serialize_json_object(object, unit)

        allocate(character(len=248) :: expected_fragment, actual_fragment)
        expected_fragment = '{"integer":-2147483647,"real": 1.7976931348623157E+308,"logical":true,"complex":[ 1.0000000000000000E+000,-1.0000000000000000E+000],"array":[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[-2147483647,'
        rewind(unit)
        read (unit, "(A)") actual_fragment
        call assert_string_equal(actual_fragment, expected_fragment, "test_serialization: case 1 Object fragments differ")
        deallocate(expected_fragment, actual_fragment)
        close(unit)

        ! Case 2: basic array serialization test with recursion
        open(newunit=unit, file="test_json.json", form='formatted', access='stream', status='replace', iostat=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_serialization: could not open file")

        call serialize_json_array(array, unit)

        allocate(character(len=108) :: expected_fragment, actual_fragment)
        expected_fragment = '[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[-2147483647,'
        rewind(unit)
        read (unit, "(A)") actual_fragment
        call assert_string_equal(actual_fragment, expected_fragment, "test_serialization: case 2 Array fragments differ")
        deallocate(expected_fragment, actual_fragment)
        close(unit)

        ! Case 3: NaN/Inf handling
        open(newunit=unit, file="test_json.json", form='formatted', access='stream', status='replace', iostat=ierr)
        call assert_equal_int(ierr, ERR_OK, "test_serialization: could not open file")

        test_complex = cmplx(ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf))
        test_real = ieee_value(1.0_real64, ieee_quiet_nan)

        call serialize_json_array(array, unit)

        allocate(character(len=48) :: expected_fragment, actual_fragment)
        expected_fragment = '[-2147483647,null,true,[null,null],[-2147483647,'
        rewind(unit)
        read (unit, "(A)") actual_fragment
        call assert_string_equal(actual_fragment, expected_fragment, "test_serialization: case 3 NaN/Inf fragments differ")
        deallocate(expected_fragment, actual_fragment)
        close(unit)
    end subroutine test_serialization

    !> Run all tox_json tests.
    subroutine run_all_tests_tox_json
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All tox_json tests passed successfully."
    end subroutine run_all_tests_tox_json

    !> Run specific tox_json tests by name.
    subroutine run_named_tests_tox_json(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_json
end module mod_test_json
