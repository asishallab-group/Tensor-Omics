module tox_json
    use, intrinsic :: iso_fortran_env, only: int32

    type :: json_array
        class(*), dimension(:), pointer :: array => null()
    end type json_array

    type :: json_object
        character(len=:), dimension(:), pointer :: keys => null()
        class(*), dimension(:), pointer :: values => null()
    end type json_object
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

    recursive subroutine serialize_json_array(json_arr, unit)
        type(json_array), intent(in) :: json_arr
        integer(int32), intent(in) :: unit


        if (associated(json_arr%array)) then
            block
                integer(int32) :: n_elements, i_element

                n_elements = size(json_arr%array, dim=1, kind=int32)
                write (unit, "('[')", advance="no")
                do i_element = 1, n_elements
                    call serialize_element(json_arr%array(i_element), unit)
                    if (i_element < n_elements) write (unit, "(',')")
                end do
                write (unit, "(']')", advance="no")
            end block
        else
            write (unit, "('null')", advance="no")
        end if
    end subroutine serialize_json_array

    recursive subroutine serialize_element(element, unit)
        class(*), intent(in) :: element
        integer(int32), intent(in) :: unit

        select type(element)
            type is (json_array)
                call serialize_json_array(element, unit)
            type is (json_object)
                call serialize_json_object(element, unit)
            class default
                call serialize_scalar(element, unit)
        end select
    end subroutine serialize_element

    recursive subroutine serialize_key_value_pair(key, value, unit)
        character(len=*), intent(in) :: key
        class(*), intent(in) :: value
        integer(int32), intent(in) :: unit

        write (unit, "('""', A, '"":')", advance="no") trim(key)
        call serialize_element(value, unit)
    end subroutine serialize_key_value_pair

    recursive subroutine serialize_json_object(json_obj, unit)
        type(json_object), intent(in) :: json_obj
        integer(int32), intent(in) :: unit

        write (unit, "('{')", advance="no")
        if (associated(json_obj%keys) .and. associated(json_obj%values)) then
            block
                integer(int32) :: n_entries, i_entry

                n_entries = size(json_obj%keys, dim=1, kind=int32)
                do i_entry = 1, n_entries
                    call serialize_key_value_pair(json_obj%keys(i_entry), json_obj%values(i_entry), unit)
                    if (i_entry < n_entries) write (unit, "(',')")
                end do
            end block
        end if
        write (unit, "('}')", advance="no")
    end subroutine serialize_json_object
end module tox_json