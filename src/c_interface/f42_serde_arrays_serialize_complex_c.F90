#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_serialize_complex(module)]]
!| Module for serializing complex arrays into files
module f42_serde_arrays_serialize_complex_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double_complex, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: serialize_complex_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_serialize_complex(module):serialize_complex_helper(subroutine)]]
    subroutine serialize_complex_helper_c(&
            arr,&
            n_elements,&
            orig_shape,&
            n_orig_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="serialize_complex_helper_c")
        use f42_serde_arrays_serialize_complex, only: serialize_complex_helper

        integer(c_int), intent(in), target :: n_elements
            !! Number of strings in `arr`
        integer(c_int), intent(in), target :: n_orig_shape_elements
            !! number of elements in `orig_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        complex(c_double_complex), dimension(n_elements), intent(in), target :: arr
            !! Array to be serialized
        integer(c_int), dimension(n_orig_shape_elements), intent(in), target :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file to write to
        integer(c_int), intent(out), target :: ierr
            !! Error code
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

        call serialize_complex_helper(&
            arr = arr,&
            n_elements = n_elements,&
            orig_shape = orig_shape,&
            filename = filename_f,&
            ierr = ierr&
        )
    end subroutine serialize_complex_helper_c

end module f42_serde_arrays_serialize_complex_c
#endif
