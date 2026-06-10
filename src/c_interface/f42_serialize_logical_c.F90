#ifndef NO_C_INTERFACE
#include <src/macros.h>

!>  Module for C-wrappers for [[f42_serialize_logical(module)]]
!| Module for serializing logical arrays to binary files.
module f42_serialize_logical_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL
contains

    !> C-wrapper for [[f42_serialize_logical(module):serialize_logical_3d(subroutine)]]
    !| 
    !| Serialize a 3D logical array to a binary file.
    !| The file will contain a magic number, type code, dimension, shape, and the array data.
    subroutine serialize_logical_3d_c(arr, n_arr_elements_dim_1, n_arr_elements_dim_2, n_arr_elements_dim_3, filename, filename_strlen, ierr) bind(C, name="serialize_logical_3d_c")
        use f42_serialize_logical, only: serialize_logical_3d
        integer(c_int), intent(in), target :: n_arr_elements_dim_1
            !!  Size of the 1. dimension/extent of `arr`
        integer(c_int), intent(in), target :: n_arr_elements_dim_2
            !!  Size of the 2. dimension/extent of `arr`
        integer(c_int), intent(in), target :: n_arr_elements_dim_3
            !!  Size of the 3. dimension/extent of `arr`
        integer(c_int), intent(in), target :: filename_strlen
            !!  String length of 'filename'
        integer(c_int), intent(in), dimension(n_arr_elements_dim_1, n_arr_elements_dim_2, n_arr_elements_dim_3), target :: arr
            !! array to save
        character(len=1, kind=c_char), intent(in), dimension(filename_strlen), target :: filename
            !! output filename
        integer(c_int), intent(out), target :: ierr
            !! error code
        logical, allocatable, dimension(:, :, :) :: arr_f
        character(len=:), allocatable :: filename_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(arr)
        M_CHECK_NON_NULL(n_arr_elements_dim_1)
        M_CHECK_NON_NULL(n_arr_elements_dim_2)
        M_CHECK_NON_NULL(n_arr_elements_dim_3)
        M_CHECK_NON_NULL(filename)
        M_CHECK_NON_NULL(filename_strlen)
        call c_int_as_logical(arr, arr_f)
        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call serialize_logical_3d(arr_f, filename_f, ierr)
    end subroutine serialize_logical_3d_c

end module f42_serialize_logical_c
#endif