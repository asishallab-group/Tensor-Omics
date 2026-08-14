#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_deserialize_int(module)]]
!| Module for deserializing integer arrays from files
module f42_serde_arrays_deserialize_int_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: deserialize_int_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_deserialize_int(module):deserialize_int_helper(subroutine)]]
    subroutine deserialize_int_helper_c(&
            arr,&
            n_elements,&
            arr_shape,&
            n_arr_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="deserialize_int_helper_c")
        use f42_serde_arrays_deserialize_int, only: deserialize_int_helper

        integer(c_int), intent(in), target :: n_arr_shape_elements
            !! number of elements in `arr_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), dimension(*), intent(out), target :: arr
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
        character(len=filename_strlen), pointer :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_elements)
        M_CHECK_NON_NULL(n_arr_shape_elements)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(arr_shape, n_arr_shape_elements)
        M_CHECK_ARRAY_NON_NULL(arr, product(arr_shape))
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)

        call c_f_pointer(c_loc(filename), filename_f)

        call deserialize_int_helper(&
            arr = arr(1:product(arr_shape)),&
            n_elements = n_elements,&
            arr_shape = arr_shape,&
            filename = filename_f,&
            ierr = ierr&
        )
    end subroutine deserialize_int_helper_c

end module f42_serde_arrays_deserialize_int_c
#endif
