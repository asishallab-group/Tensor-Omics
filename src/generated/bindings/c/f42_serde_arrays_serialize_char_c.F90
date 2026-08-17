#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_serialize_char(module)]]
!| Module for serializing character arrays into files
module f42_serde_arrays_serialize_char_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: serialize_char_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_serialize_char(module):serialize_char_helper(subroutine)]]
    subroutine serialize_char_helper_c(&
            arr,&
            arr_strlen,&
            n_strings,&
            arr_shape,&
            n_arr_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="serialize_char_helper_c")
        use f42_serde_arrays_serialize_char, only: serialize_char_helper

        integer(c_int), intent(in), target :: arr_strlen
            !! length of the strings in `arr`
        integer(c_int), intent(in), target :: n_arr_shape_elements
            !! number of elements in `arr_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        character(len=1, kind=c_char), dimension(arr_strlen, *), intent(in), target :: arr
            !! Array to be serialized
        integer(c_int), intent(in), target :: n_strings
            !! Number of strings in `arr`
        integer(c_int), dimension(n_arr_shape_elements), intent(in), target :: arr_shape
            !! Extents of `arr`, one per dimension
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file to write to
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=arr_strlen), pointer, dimension(:) :: arr_f
        character(len=filename_strlen), pointer :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(arr_strlen)
        M_CHECK_NON_NULL(n_strings)
        M_CHECK_NON_NULL(n_arr_shape_elements)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(arr_shape, n_arr_shape_elements)
        M_CHECK_ARRAY_NON_NULL(arr, product(arr_shape))
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)

        call c_f_pointer(c_loc(arr), arr_f, [product(arr_shape)])
        call c_f_pointer(c_loc(filename), filename_f)

        call serialize_char_helper(&
            arr = arr_f,&
            n_strings = n_strings,&
            arr_shape = arr_shape,&
            filename = filename_f,&
            ierr = ierr&
        )
    end subroutine serialize_char_helper_c

end module f42_serde_arrays_serialize_char_c
#endif
