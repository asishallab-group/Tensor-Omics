#include <src/macros.h>

!> Module for C-wrappers for [[f42_deserialize_char(module)]]
!| Module for deserializing character arrays from files
module f42_deserialize_char_c
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

    !> C-wrapper for [[f42_deserialize_char(module):deserialize_char_1d(subroutine)]]
    !| Directly deserialize a 1D character array from a file (array already allocated)
    subroutine deserialize_char_1d_c(filename_strlen, n_arr_elements, arr_strlen, arr, filename, ierr) bind(C, name="deserialize_char_1d_c")
        use f42_deserialize_char, only: deserialize_char_1d
        integer(c_int), intent(in), target :: filename_strlen
            !! String length of 'filename'
        integer(c_int), intent(in), target :: n_arr_elements
            !!  Size of the 1. dimension/extent of `arr`
        integer(c_int), intent(in), target :: arr_strlen
            !! String length of 'arr'
        character(len=1, kind=c_char), dimension(arr_strlen, n_arr_elements), intent(out), target :: arr
            !! Pre-allocated array to read the data into
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file to read from
        integer(c_int), intent(out), target :: ierr
            !! Error code
            !!
    
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(n_arr_elements)
        M_CHECK_NON_NULL(arr_strlen)
        M_CHECK_NON_NULL(arr)
        M_CHECK_NON_NULL(filename)
    
        call deserialize_char_1d(arr, filename, ierr)
    end subroutine deserialize_char_1d_c

end module f42_deserialize_char_c