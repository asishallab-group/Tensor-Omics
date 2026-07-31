#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_deserialize_complex(module)]]
!| Module for deserializing complex arrays from files
module f42_serde_arrays_deserialize_complex_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double_complex, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: deserialize_complex_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_deserialize_complex(module):deserialize_complex_helper(subroutine)]]
    subroutine deserialize_complex_helper_c(&
            arr,&
            n_elements,&
            arr_shape,&
            n_arr_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="deserialize_complex_helper_c")
        use f42_serde_arrays_deserialize_complex, only: deserialize_complex_helper

        integer(c_int), intent(in), target :: n_arr_shape_elements
            !! number of elements in `arr_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        complex(c_double_complex), dimension(*), intent(out), target :: arr
            !! Pre-allocated array to read the data into
        integer(c_int), intent(in), target :: n_elements
            !! Size of `arr`
        integer(c_int), dimension(n_arr_shape_elements), intent(in), target :: arr_shape
            !! Extents of `arr`, one per dimension
            !! It is *VERY IMPORTANT* to compute this argument from the `dims_out` output produced by [[f42_serde_arrays_utils(module):get_array_metadata]].
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! | dims_out_capacity | 5_int32     |
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=:), allocatable :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_elements)
        M_CHECK_NON_NULL(n_arr_shape_elements)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(arr_shape, n_arr_shape_elements)
        M_CHECK_ARRAY_NON_NULL(arr, product(arr_shape))
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)

        call c_char_1d_as_string(filename, filename_f, ierr)
        if (is_err(ierr)) return

        call deserialize_complex_helper(&
            arr = arr(1:product(arr_shape)),&
            n_elements = n_elements,&
            arr_shape = arr_shape,&
            filename = filename_f,&
            ierr = ierr&
        )
    end subroutine deserialize_complex_helper_c

end module f42_serde_arrays_deserialize_complex_c
#endif
