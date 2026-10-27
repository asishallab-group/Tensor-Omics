#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_deserialize_char(module)]]
!| Module for deserializing character arrays from files
module f42_serde_arrays_deserialize_char_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, c_int, c_loc, c_ptr
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_EMPTY_INPUT, ERR_SIZE_MISMATCH
    M_IMPLICIT_NONE
    private

    character(kind=c_char, len=1), target :: empty_view_anchor

    public :: deserialize_char_helper_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_deserialize_char(module):deserialize_char_helper(subroutine)]]
    subroutine deserialize_char_helper_c(&
            arr,&
            n_strings,&
            strlen,&
            arr_shape,&
            n_arr_shape_elements,&
            filename,&
            filename_strlen,&
            ierr&
        ) bind(C, name="deserialize_char_helper_c")
        use f42_serde_arrays_deserialize_char, only: deserialize_char_helper

        integer(c_int), intent(in), target :: strlen
            !! Length of each string in `arr`.
            !! It is *VERY IMPORTANT* to compute this argument from the `type_code` output produced by [[f42_serde_arrays_utils(module):get_array_metadata]].
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! | dims_out_capacity | 5_int32     |
        integer(c_int), intent(in), target :: n_arr_shape_elements
            !! number of elements in `arr_shape`
        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        character(len=1, kind=c_char), dimension(strlen, *), intent(out), target :: arr
            !! Pre-allocated array to read the data into
        integer(c_int), intent(in), target :: n_strings
            !! Number of strings in `arr`
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
        character(len=strlen), pointer, dimension(:) :: arr_f
        type(c_ptr) :: arr_address
        character(len=filename_strlen), pointer :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_strings)
        M_CHECK_NON_NULL(strlen)
        M_CHECK_NON_NULL(n_arr_shape_elements)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_ARRAY_NON_NULL(arr_shape, n_arr_shape_elements)
        M_CHECK_ARRAY_NON_NULL(arr, product(arr_shape))
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        if (product(arr_shape) /= n_strings) then
            call set_err(ierr, ERR_SIZE_MISMATCH, 1_int32)
            return
        end if

        M_CHECK_CHARACTER_ARRAY_VIEW(arr, strlen, product(arr_shape))
        arr_address = c_loc(empty_view_anchor)
        if (product(arr_shape) > 0) arr_address = c_loc(arr)
        call c_f_pointer(arr_address, arr_f, [product(arr_shape)])
        M_CHECK_CHARACTER_VIEW(filename, filename_strlen)
        call c_f_pointer(c_loc(filename), filename_f)

        call deserialize_char_helper(&
            arr = arr_f,&
            n_strings = n_strings,&
            strlen = strlen,&
            arr_shape = arr_shape,&
            filename = filename_f,&
            ierr = ierr&
        )
    end subroutine deserialize_char_helper_c

end module f42_serde_arrays_deserialize_char_c
#endif
