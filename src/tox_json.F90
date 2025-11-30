module tox_json
    use, intrinsic :: iso_fortran_env, only: int32

    private

    type :: json_value
        class(*), pointer :: value => null()
    end type json_value

    type :: json_array
        class(*), dimension(:), pointer :: array => null()
    end type json_array

    type :: json_object
        character(len=:), dimension(:), pointer :: keys => null()
        class(*), dimension(:), pointer :: values => null()
    end type json_object

    integer(int32) :: MAX_RECURSION_DEPTH = 20

    public :: json_array, json_object, json_value
    public :: serialize_json_array, serialize_json_object
contains

    subroutine serialize_scalar(scalar, unit)
        class(*), intent(in) :: scalar
        integer(int32), intent(in) :: unit

        select type(scalar)
            type is (integer)
                write (unit, "(G0)", advance="no") scalar
            type is (real)
                write (unit, "(G0)", advance="no") scalar
            type is (logical)
                write (unit, "(A)", advance="no") trim(merge("true ", "false", scalar))
            type is (character(*))
                write (unit, "('""', A, '""')", advance="no") trim(scalar)
            type is (complex)
                write (unit, "('[',G0,',',G0,']')", advance="no") scalar
            class default
                write (unit, "('null')", advance="no")
        end select
    end subroutine serialize_scalar

    recursive subroutine serialize_array(json_arr, unit, depth)
        type(json_array), intent(in) :: json_arr
        integer(int32), intent(in) :: unit
        integer(int32), intent(inout) :: depth

        depth = depth + 1
        if (depth < MAX_RECURSION_DEPTH) then
            write (unit, "('[')", advance="no")
            if (associated(json_arr%array)) then
                block
                    integer(int32) :: n_elements, i_element

                    n_elements = size(json_arr%array, dim=1, kind=int32)
                    do i_element = 1, n_elements
                        call serialize_json_value(json_arr%array(i_element), unit, depth)
                        if (i_element < n_elements) write (unit, "(',')")
                    end do
                end block
            end if
            write (unit, "(']')", advance="no")
        else
            write (unit, "('null')", advance="no")
        end if
        depth = depth - 1
    end subroutine serialize_array

    recursive subroutine serialize_json_value(element, unit, depth)
        class(*), intent(in) :: element
        integer(int32), intent(in) :: unit
        integer(int32), intent(inout) :: depth

        select type(element)
            type is (json_array)
                call serialize_array(element, unit, depth)
            type is (json_object)
                call serialize_object(element, unit, depth)
            type is (json_value)
                depth = depth + 1
                if (associated(element%value) .and. depth < MAX_RECURSION_DEPTH) then
                    call serialize_json_value(element%value, unit, depth)
                else
                    write (unit, "('null')", advance="no")
                end if
                depth = depth - 1
            class default
                call serialize_scalar(element, unit)
        end select
    end subroutine serialize_json_value

    recursive subroutine serialize_key_value_pair(key, value, unit, depth)
        character(len=*), intent(in) :: key
        class(*), intent(in) :: value
        integer(int32), intent(in) :: unit
        integer(int32), intent(inout) :: depth

        write (unit, "('""', A, '"":')", advance="no") trim(key)
        call serialize_json_value(value, unit, depth)
    end subroutine serialize_key_value_pair

    recursive subroutine serialize_object(json_obj, unit, depth)
        type(json_object), intent(in) :: json_obj
        integer(int32), intent(in) :: unit
        integer(int32), intent(inout) :: depth

        depth = depth + 1
        if (depth < MAX_RECURSION_DEPTH) then
            write (unit, "('{')", advance="no")
            if (associated(json_obj%keys) .and. associated(json_obj%values)) then
                block
                    integer(int32) :: n_entries, i_entry

                    n_entries = size(json_obj%keys, dim=1, kind=int32)
                    do i_entry = 1, n_entries
                        call serialize_key_value_pair(json_obj%keys(i_entry), json_obj%values(i_entry), unit, depth)
                        if (i_entry < n_entries) write (unit, "(',')")
                    end do
                end block
            end if
            write (unit, "('}')", advance="no")
        else
            write (unit, "('null')", advance="no")
        end if
        depth = depth - 1
    end subroutine serialize_object

    subroutine serialize_json_object(json_obj, unit)
        type(json_object), intent(in) :: json_obj
        integer(int32), intent(in) :: unit

        integer(int32) :: depth

        depth = 0_int32
        call serialize_object(json_obj, unit, depth)
    end subroutine serialize_json_object

    subroutine serialize_json_array(json_arr, unit)
        type(json_array), intent(in) :: json_arr
        integer(int32), intent(in) :: unit

        integer(int32) :: depth

        depth = 0_int32
        call serialize_array(json_arr, unit, depth)
    end subroutine serialize_json_array
end module tox_json

! Test idea:
!   - an array of all possible types
!   - for simplicity the array type is a self reference and tests recursion simultaneously
!   - also, the self-containing array is reused as values of the object, and the object is also included in the array
! 
! type(json_value), dimension(6), target :: elements
! type(json_array) :: array
! type(json_object) :: object
! character(len=7), dimension(4) :: keys
! keys(1) = "integer"
! elements(1)%value => -1_int32
! keys(2) = "real"
! elements(2)%value => 1.0_real64
! keys(3) = "logical"
! elements(3)%value => .true.
! keys(4) = "complex"
! elements(4)%value => cmplx(1.0_real64, -1.0_real64)
! keys(5) = "array"
! elements(5)%value => array
! array%array => elements
! keys(6) = "object"
! elements(6)%value => object
! object%keys => keys
! object%values => elements
! 
! ! Output for serializing `array` should be something like: [-1, 1.0, true, "(1.0, -1.0)", [...], {"integer": -1, ...}]
! ! Output for serializing `object` should be something like: {"integer": -1, "real": 1.0, "logical": true, "complex": "(1.0, -1.0)", "array": [...], "object": {"integer": -1, ...}}

! Test object first, then array. Just validate a fragment of the output
!   - for object: until complex value in array
!   - for array: until complex value
