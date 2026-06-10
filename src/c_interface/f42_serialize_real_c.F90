#ifndef NO_C_INTERFACE
#include <src/macros.h>

!>  Module for C-wrappers for [[f42_serialize_real(module)]]
module f42_serialize_real_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL
contains

    !> C-wrapper for [[f42_serialize_real(module):serialize_real_nd(subroutine)]]
    !| 
    !| Writes serialized real array from R to file with metdata.
    subroutine serialize_real_nd_c(arr, n_arr_elements, dims, n_dims_elements, ndim, filename, filename_strlen, ierr) bind(C, name="serialize_real_nd_c")
        use f42_serialize_real, only: serialize_real_nd
        integer(c_int), intent(in), target :: n_arr_elements
            !!  Size of the 1. dimension/extent of `arr`
        integer(c_int), intent(in), target :: n_dims_elements
            !!  Size of the 1. dimension/extent of `dims`
        integer(c_int), intent(in), target :: filename_strlen
            !!  String length of 'filename'
        real(c_double), intent(in), dimension(n_arr_elements), target :: arr
            !! array to save
        integer(c_int), intent(in), dimension(n_dims_elements), target :: dims
            !! Dimensions of the array
        integer(c_int), intent(in), target :: ndim
            !! Number of dimensions
        character(len=1, kind=c_char), intent(in), dimension(filename_strlen), target :: filename
            !! filename
        integer(c_int), intent(out), target :: ierr
            !! error code
        character(len=:), allocatable :: filename_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(arr)
        M_CHECK_NON_NULL(n_arr_elements)
        M_CHECK_NON_NULL(dims)
        M_CHECK_NON_NULL(n_dims_elements)
        M_CHECK_NON_NULL(ndim)
        M_CHECK_NON_NULL(filename)
        M_CHECK_NON_NULL(filename_strlen)
        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call serialize_real_nd(arr, dims, ndim, filename_f, ierr)
    end subroutine serialize_real_nd_c

end module f42_serialize_real_c
#endif