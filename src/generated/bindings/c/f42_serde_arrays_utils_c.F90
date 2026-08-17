#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_serde_arrays_utils(module)]]
!| Module for array utilities.
!|
!| Defines the shared on-disk binary layout used by all typed array
!| serialize/deserialize modules (int/real/complex/logical/char) and the
!| header read/write/validate helpers that implement it. The file header is
!| a fixed sequence of unformatted stream records, written and read in this
!| order: magic number ([[f42_serde_arrays_utils(module):ARRAY_FILE_MAGIC(variable)]]),
!| type code, number of dimensions `ndim`, then `ndim` dimension sizes. The
!| raw array payload follows immediately after the header, written as one
!| contiguous block by the type-specific serializers.
!|
!| The header does NOT record the width of an element, so the payload is only
!| readable by a build that agrees with the writer on the storage size of the
!| type code. Logical arrays changed width once: they are now written as
!| `logical(c_bool)`, one byte per element, where earlier builds wrote the
!| default logical kind at four. A logical `.bin` file written before that
!| change therefore decodes to garbage here rather than failing, and one
!| written here does the same there. No other type code has moved.
module f42_serde_arrays_utils_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_f_pointer, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: get_array_metadata_c

contains

    !> summary: C-wrapper for [[f42_serde_arrays_utils(module):get_array_metadata(subroutine)]]
    subroutine get_array_metadata_c(&
            filename,&
            filename_strlen,&
            dims_out,&
            dims_out_capacity,&
            ndims,&
            type_code,&
            ierr&
        ) bind(C, name="get_array_metadata_c")
        use f42_serde_arrays_utils, only: get_array_metadata

        integer(c_int), intent(in), target :: filename_strlen
            !! length of the strings in `filename`
        integer(c_int), intent(in), target :: dims_out_capacity
            !! Capacity of the dims_out array
        character(len=1, kind=c_char), dimension(filename_strlen), intent(in), target :: filename
            !! Name of the file
        integer(c_int), dimension(dims_out_capacity), intent(out), target :: dims_out
            !! Array to store output dimensions
            !! The first `ndims` elements will hold the results.
        integer(c_int), intent(out), target :: ndims
            !! number of dimensions
        integer(c_int), intent(out), target :: type_code
            !! Type code of the serialized array
            !!
            !!
            !! | Type      | Code          |
            !! |-----------|---------------|
            !! | integer   | -1_int32      |
            !! | real      | -2_int32      |
            !! | complex   | -4_int32      |
            !! | logical   | -3_int32      |
            !! | character | string length |
        integer(c_int), intent(out), target :: ierr
            !! Error code
        character(len=filename_strlen), pointer :: filename_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(filename_strlen)
        M_CHECK_NON_NULL(dims_out_capacity)
        M_CHECK_NON_NULL(ndims)
        M_CHECK_NON_NULL(type_code)
        M_CHECK_ARRAY_NON_NULL(filename, filename_strlen)
        M_CHECK_ARRAY_NON_NULL(dims_out, dims_out_capacity)

        call c_f_pointer(c_loc(filename), filename_f)

        call get_array_metadata(&
            filename = filename_f,&
            dims_out = dims_out,&
            dims_out_capacity = dims_out_capacity,&
            ndims = ndims,&
            type_code = type_code,&
            ierr = ierr&
        )
    end subroutine get_array_metadata_c

end module f42_serde_arrays_utils_c
#endif
