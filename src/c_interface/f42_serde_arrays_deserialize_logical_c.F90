#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_deserialize_logical(module)]]
!| Module for deserializing logical arrays from files
module f42_serde_arrays_deserialize_logical_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: deserialize_logical_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_deserialize_logical(module):deserialize_logical_helper(subroutine)]]
    subroutine deserialize_logical_helper_c(&
            arr,&
            n_elements,&
            orig_shape,&
            n_orig_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="deserialize_logical_helper_c")
        use f42_serde_arrays_deserialize_logical, only: deserialize_logical_helper

        integer(c_int), intent(in), target :: n_elements
            !! Size of `arr`
        integer(c_int), intent(in), target :: n_orig_shape_elements
            !! number of elements in `orig_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        logical(c_bool), dimension(n_elements), intent(out), target :: arr
            !! Pre-allocated array to read the data into
        integer(c_int), dimension(n_orig_shape_elements), intent(in), target :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, dimension(n_elements) :: arr_f
        character(len=:), allocatable :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_elements)
        M_CHECK_NON_NULL(n_orig_shape_elements)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(arr, n_elements)
        M_CHECK_ARRAY_NON_NULL(orig_shape, n_orig_shape_elements)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return

        call deserialize_logical_helper(&
            arr = arr_f,&
            n_elements = n_elements,&
            orig_shape = orig_shape,&
            filename = filename_f,&
            ierr = ierr&
        )

        arr = arr_f
    end subroutine deserialize_logical_helper_c

end module f42_serde_arrays_deserialize_logical_c
#endif
