#ifndef NO_C_INTERFACE
#include <src/macros.h>

!>  Module for C-wrappers for [[f42_serialize_char(module)]]
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

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL
contains

    !> C-wrapper for [[f42_serialize_char(module):serialize_char_nd(subroutine)]]
    !| 
    !| Serialize a character array of arbitrary dimensions to a binary file.
    !| The file will contain a magic number, type code, dimension, shape, character length, and the array data.
    !| @note This routine is only called by R and serializes only flat character arrays to the memory
    subroutine serialize_char_nd_c(flat, flat_strlen, flat_shape, n_flat_shape_elements, filename, filename_strlen, ierr) bind(C, name="serialize_char_nd_c")
        use f42_serialize_char, only: serialize_char_nd
        integer(c_int), intent(in), target :: flat_strlen
            !!  String length of 'flat'
        integer(c_int), intent(in), target :: n_flat_shape_elements
            !!  Size of the 1. dimension/extent of `flat_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !!  String length of 'filename'
        character(len=1, kind=c_char), intent(in), dimension(flat_strlen, *), target :: flat
            !! flat array to save
        integer(c_int), intent(in), dimension(n_flat_shape_elements), target :: flat_shape
            !! dimensions of the array
        character(len=1, kind=c_char), intent(in), dimension(filename_strlen), target :: filename
            !! output filename
        integer(c_int), intent(out), target :: ierr
            !! error code
        character(len=:), allocatable, dimension(:) :: flat_f
        character(len=:), allocatable :: filename_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(flat)
        M_CHECK_NON_NULL(flat_strlen)
        M_CHECK_NON_NULL(flat_shape)
        M_CHECK_NON_NULL(n_flat_shape_elements)
        M_CHECK_NON_NULL(filename)
        M_CHECK_NON_NULL(filename_strlen)
        call c_char_2d_as_string(flat(:, 1:size(flat_shape, kind=int32)), flat_f, ierr)
        if (is_err(ierr)) return
        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call serialize_char_nd(flat_f, flat_shape, filename_f, ierr)
    end subroutine serialize_char_nd_c

end module f42_serialize_char_c
#endif