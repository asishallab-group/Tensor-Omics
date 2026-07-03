#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[f42_serialize_char(module)]]
!| Module providing serialization and deserialization routines for character arrays
!| of up to 5 dimensions, arrays are serialized to a custom binary format with a magic number and type/dimension metadata.
module f42_serialize_char_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[f42_serialize_char(module):serialize_char_1d(subroutine)]]
    !| Serialize a 1D character array to a binary file.
    !| The file will contain a magic number, type code, dimension, shape, character length, and the array data.
    subroutine serialize_char_1d_c(&
            arr,&
            arr_strlen,&
            n_arr_elements,&
            filename,&
            filename_strlen,&
            ierr&
            ) bind(C, name="serialize_char_1d_c")
        use f42_serialize_char, only: serialize_char_1d
        use f42_serialize_char
        integer(c_int), intent(in), target :: arr_strlen
            !! String length of 'arr'
        integer(c_int), intent(in), target :: n_arr_elements
            !! Size of the 1. dimension/extent of `arr`
        integer(c_int), intent(in), target :: filename_strlen
            !! String length of 'filename'
        character(len=1, kind=c_char), intent(in), dimension(arr_strlen, n_arr_elements), target :: arr
            !! array to save
        character(len=1, kind=c_char), intent(in), dimension(filename_strlen), target :: filename
            !! output filename
        integer(c_int), intent(out), target :: ierr
            !! error code
        character(len=:), allocatable, dimension(:) :: arr_f
        character(len=:), allocatable :: filename_f

        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(arr)
        M_CHECK_NON_NULL(arr_strlen)
        M_CHECK_NON_NULL(n_arr_elements)
        M_CHECK_NON_NULL(filename)
        M_CHECK_NON_NULL(filename_strlen)

        call c_char_2d_as_string(arr, arr_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call serialize_char_1d(&
            arr = arr_f,&
            filename = filename_f,&
            ierr = ierr&
        )

    end subroutine serialize_char_1d_c

end module f42_serialize_char_c
#endif