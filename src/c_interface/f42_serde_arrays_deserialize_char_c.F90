#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_deserialize_char(module)]]
!| Module for deserializing character arrays from files
module f42_serde_arrays_deserialize_char_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string, string_as_c_char_2d
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: deserialize_char_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_deserialize_char(module):deserialize_char_helper(subroutine)]]
    subroutine deserialize_char_helper_c(&
            arr,&
            arr_strlen,&
            n_strings,&
            orig_shape,&
            n_orig_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="deserialize_char_helper_c")
        use f42_serde_arrays_deserialize_char, only: deserialize_char_helper

        integer(c_int), intent(in), target :: arr_strlen
            !! length of the strings in `arr`
        integer(c_int), intent(in), target :: n_strings
            !! Number of strings in `arr`
        integer(c_int), intent(in), target :: n_orig_shape_elements
            !! number of elements in `orig_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        character(len=1, kind=c_char), dimension(arr_strlen, n_strings), intent(out), target :: arr
            !! Pre-allocated array to read the data into
        integer(c_int), dimension(n_orig_shape_elements), intent(in), target :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable, dimension(:) :: arr_f
        character(len=:), allocatable :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(arr_strlen)
        M_CHECK_NON_NULL(n_strings)
        M_CHECK_NON_NULL(n_orig_shape_elements)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(arr, arr_strlen * n_strings)
        M_CHECK_ARRAY_NON_NULL(orig_shape, n_orig_shape_elements)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)

        allocate(character(len=arr_strlen) :: arr_f(n_strings))
        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return

        call deserialize_char_helper(&
            arr = arr_f,&
            n_strings = n_strings,&
            orig_shape = orig_shape,&
            filename = filename_f,&
            ierr = ierr&
        )

        call string_as_c_char_2d(arr_f, arr)
    end subroutine deserialize_char_helper_c

end module f42_serde_arrays_deserialize_char_c
#endif
