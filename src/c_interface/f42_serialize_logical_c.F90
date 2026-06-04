#include <src/macros.h>

!> Module for C-wrappers for [[f42_serialize_logical(module)]]
!| Module for serializing logical arrays to binary files.
module f42_serialize_logical_c
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err
contains

    !> C-wrapper for [[f42_serialize_logical(module):serialize_logical_3d(subroutine)]]
    !| Serialize a 3D logical array to a binary file.
    !| The file will contain a magic number, type code, dimension, shape, and the array data.
    subroutine serialize_logical_3d_c(filename_strlen, n_arr_elements_dim_1, n_arr_elements_dim_2, n_arr_elements_dim_3, arr, filename, ierr) bind(C, name="serialize_logical_3d_c")
        use f42_serialize_logical, only: serialize_logical_3d
        integer(c_int), intent(in), target :: filename_strlen
            !! String length of 'filename'
        integer(c_int), intent(in), target :: n_arr_elements_dim_1
            !!  Size of the 1. dimension/extent of `arr`
        integer(c_int), intent(in), target :: n_arr_elements_dim_2
            !!  Size of the 2. dimension/extent of `arr`
        integer(c_int), intent(in), target :: n_arr_elements_dim_3
            !!  Size of the 3. dimension/extent of `arr`
        integer(c_int), dimension(n_arr_elements_dim_1, n_arr_elements_dim_2, n_arr_elements_dim_3), intent(in), target :: arr
            !! array to save
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! output filename
        integer(c_int), intent(out), target :: ierr
            !! error code
    
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(n_arr_elements_dim_1)
        M_CHECK_NON_NULL(n_arr_elements_dim_2)
        M_CHECK_NON_NULL(n_arr_elements_dim_3)
        M_CHECK_NON_NULL(arr)
        M_CHECK_NON_NULL(filename)
    
        call serialize_logical_3d(arr, filename, ierr)
    end subroutine serialize_logical_3d_c

end module f42_serialize_logical_c