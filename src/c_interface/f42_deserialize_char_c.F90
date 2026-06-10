#ifndef NO_C_INTERFACE
#include <src/macros.h>

!>  Module for C-wrappers for [[f42_deserialize_char(module)]]
!| Module for deserializing character arrays from files
module f42_deserialize_char_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL
contains

    !> C-wrapper for [[f42_deserialize_char(module):deserialize_char_nd(subroutine)]]
    !| 
    !| Subroutine to deserialize a flat character array from a file
    subroutine deserialize_char_nd_c(flat, flat_strlen, n_flat_elements, filename, filename_strlen, ierr) bind(C, name="deserialize_char_nd_c")
        use f42_deserialize_char, only: deserialize_char_nd
        integer(c_int), intent(in), target :: flat_strlen
            !!  String length of 'flat'
        integer(c_int), intent(in), target :: n_flat_elements
            !!  Size of the 1. dimension/extent of `flat`
        integer(c_int), intent(in), target :: filename_strlen
            !!  String length of 'filename'
        character(len=1, kind=c_char), intent(out), dimension(flat_strlen, n_flat_elements), target :: flat
            !! Output flat character array
        character(len=1, kind=c_char), intent(in), dimension(filename_strlen), target :: filename
            !! Name of the file to read
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable, dimension(:) :: flat_f
        character(len=:), allocatable :: filename_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(flat)
        M_CHECK_NON_NULL(flat_strlen)
        M_CHECK_NON_NULL(n_flat_elements)
        M_CHECK_NON_NULL(filename)
        M_CHECK_NON_NULL(filename_strlen)
        M_ALLOCATE(character(len=flat_strlen) :: flat_f(n_flat_elements))
        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return
        call deserialize_char_nd(flat_f, filename_f, ierr)
        call string_as_c_char_2d(flat_f, flat)
    end subroutine deserialize_char_nd_c

end module f42_deserialize_char_c
#endif