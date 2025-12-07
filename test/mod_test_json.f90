!> Unit test suite for f42_json routine.
module mod_test_json
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use f42_json
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

        all_tests(1) = test_case("test_f42_json_serialization", test_serialization)
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
        character(len=:), dimension(:), allocatable, target :: keys
        integer(int32) :: i_element
        integer(int32), target:: test_integer
        real(real64), target:: test_real
        complex(real64), target:: test_complex
        logical, target:: test_logical

        allocate(character(len=7) :: keys(6))

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
        array%elements => elements
        elements(5)%value => array

        keys(6) = "object"
        object%keys => keys
        object%values => elements
        elements(6)%value => object

        ! Case 1: basic object serialization test with recursion
        call helper_test_serialization(&
            object,&
            '{"integer":-2147483647,"real": 1.7976931348623157E+308,"logical":true,"complex":[ 1.0000000000000000E+000,-1.0000000000000000E+000],"array":[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[-2147483647,',&
            "test_serialization: case 1 Object"&
        )

        ! Case 2: basic array serialization test with recursion
        call helper_test_serialization(&
            array,&
            '[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[-2147483647,',&
            "test_serialization: case 2 Array"&
        )

        ! Case 3: basic array serialization test with custom recursion limit
        call helper_test_serialization(&
            array,&
            '[-2147483647, 1.7976931348623157E+308,true,[ 1.0000000000000000E+000,-1.0000000000000000E+000],[null,null,null,null,null,null],{"integer":null,"real":null,"logical":null,"complex":null,"array":null,"object":null}]',&
            "test_serialization: case 3 custom recursion limit",&
            2_int32&
        )

        ! Case 4: NaN/Inf handling
        test_complex = cmplx(ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf))
        test_real = ieee_value(1.0_real64, ieee_quiet_nan)

        call helper_test_serialization(&
            array,&
            '[-2147483647,null,true,[null,null],[-2147483647,',&
            "test_serialization: case 4 NaN/Inf handling"&
        )

        ! Case 5: Different lengths of json object's key-value arrays
        object%keys => keys(1:4)
        call helper_test_serialization(&
            object,&
            '{"integer":-2147483647,"real":null,"logical":true,"complex":[null,null]}',&
            "test_serialization: case 5 different size key-value arrays"&
        )

        ! Case 6: unassigned array values
        do i_element = 1, 4
            nullify(elements(i_element)%value)
        end do
        nullify(object%keys)
        nullify(elements(5)%value)

        call helper_test_serialization(&
            array,&
            '[null,null,null,null,null,{}]',&
            "test_serialization: case 6 unassigned array values"&
        )

        ! Case 7: unassigned object values
        elements(5)%value => array
        object%keys => keys
        nullify(array%elements)
        nullify(elements(6)%value)

        call helper_test_serialization(&
            object,&
            '{"integer":null,"real":null,"logical":null,"complex":null,"array":[],"object":null}',&
            "test_serialization: case 7 unassigned object values"&
        )

        ! Case 8: strings
        deallocate(keys)
        allocate(character(len=119) :: keys(1))
        keys(1) = 'we test "quoting", /slashing and \backslashing, 	tabbing, ' //&
                  achar(8) //  'backspacing, ' //&
                  achar(10) // 'new-lining, ' //&
                  achar(13) // 'carriage-returning, ' //&
                  achar(12) // 'form-feeding'
        elements(1)%value => keys(1)
        object%keys => keys(1:1)
        object%values => keys(1:1)


        call helper_test_serialization(&
            object,&
            '{"we test \"quoting\", \/slashing and \\backslashing, \ttabbing, \bbackspacing, \nnew-lining, \rcarriage-returning, \fform-feeding"' //&
            ':"we test \"quoting\", \/slashing and \\backslashing, \ttabbing, \bbackspacing, \nnew-lining, \rcarriage-returning, \fform-feeding"}',&
            "test_serialization: case 8 strings"&
        )
    end subroutine test_serialization

    subroutine helper_test_serialization(json_var, expected_fragment, test_case, max_depth)
        class(*), intent(in) :: json_var
        character(len=*), intent(in) :: expected_fragment
        character(len=*), intent(in) :: test_case
        integer(int32), intent(in), optional :: max_depth

        character(len=:), allocatable :: actual_fragment
        integer(int32) :: ierr, unit

        open(newunit=unit, file="test_json.json", form='formatted', access='stream', status='replace', iostat=ierr)
        call assert_equal_int(ierr, ERR_OK, trim(test_case) // ": could not open file")

        select type (json_var)
            type is (json_array)
                call serialize_json_array(json_var, unit, max_depth)
            type is (json_object)
                call serialize_json_object(json_var, unit, max_depth)
            class default
                error stop trim(test_case) // ": helper_test_serialization: Unsupported type"
        end select

        rewind(unit)

        allocate(actual_fragment, mold=expected_fragment)
        read (unit, "(A)") actual_fragment
        call assert_string_equal(actual_fragment, expected_fragment, trim(test_case) // ": fragments differ")
        close(unit, status="delete")
    end subroutine helper_test_serialization

    !> Run all f42_json tests.
    subroutine run_all_tests_f42_json
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All f42_json tests passed successfully."
    end subroutine run_all_tests_f42_json

    !> Run specific f42_json tests by name.
    subroutine run_named_tests_f42_json(test_names)
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
    end subroutine run_named_tests_f42_json
end module mod_test_json
