#include <src/macros.h>

!> Module for serializing logical arrays into files
module f42_serde_arrays_serialize_logical
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_arrays_utils, only: write_file_header, LOGICAL_TYPE_CODE
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_WRITE_DATA, set_err
    implicit none

    private
    public :: serialize_logical_1d, serialize_logical_2d, serialize_logical_3d, &
              serialize_logical_4d, serialize_logical_5d, serialize_logical_helper

contains
    !> AUTHOR_AARON_SCHROEDER
    !| Subroutine to serialize a flat logical array into a file
    subroutine serialize_logical_helper(arr, n_elements, orig_shape, filename, ierr)
        integer(int32), intent(in) :: n_elements
            !! Number of strings in `arr`
        logical, dimension(n_elements), intent(in) :: arr
            !! Array to be serialized
        integer(int32), dimension(:), intent(in) :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=*), intent(in) :: filename
            !! Name of the file to write to
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_elements, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call write_file_header(filename, unit, LOGICAL_TYPE_CODE, size(orig_shape, kind=int32), orig_shape, ierr)
        if (is_err(ierr)) return

        ! Write the entire array as a contiguous block
        write (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_WRITE_DATA)

        close (unit)
    end subroutine serialize_logical_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 1D logical array into a file
    subroutine serialize_logical_1d(arr, filename, ierr)
        logical, dimension(:), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_logical_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 2D logical array into a file
    subroutine serialize_logical_2d(arr, filename, ierr)
        logical, dimension(:, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_logical_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 3D logical array into a file
    subroutine serialize_logical_3d(arr, filename, ierr)
        logical, dimension(:, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_logical_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 4D logical array into a file
    subroutine serialize_logical_4d(arr, filename, ierr)
        logical, dimension(:, :, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_logical_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 5D logical array into a file
    subroutine serialize_logical_5d(arr, filename, ierr)
        logical, dimension(:, :, :, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_logical_5d

end module f42_serde_arrays_serialize_logical

!> C binding for the subroutine to serialize a flat logical array into a file.
subroutine serialize_logical_nd_c(arr, orig_shape, n_dims, &
                                 filename, fn_len, ierr) bind(C, name="serialize_logical_nd_c")
    use, intrinsic :: iso_c_binding, only: c_char, c_int
    use f42_serde_arrays_serialize_logical, only: serialize_logical_helper
    use tox_conversions, only: c_char_1d_as_string, c_int_as_logical
    use tox_errors, only: is_err, map_err_arg_pos, ERR_ALLOC_FAIL, set_err
    M_USE_NULL_VALIDATION
    implicit none

    ! Arguments
    integer(c_int), intent(in), target :: n_dims
        !! Number of dimensions of the expected array (`size(orig_shape)`)
    integer(c_int), dimension(n_dims), intent(in), target :: orig_shape
        !! Original shape of the flattened array `arr` -> for a 2D array it would be `[n_elements, n_contained_arrays]`
    integer(c_int), dimension(product(orig_shape)), intent(in), target :: arr
        !! Input array to be serialized, (`.false.<=>0` and `.true.<=>1`)
    integer(c_int), intent(in), target :: fn_len
        !! Length of the filename
    character(kind=c_char, len=1), dimension(fn_len), intent(in), target  :: filename
        !! c_char array representing the filename
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: filename_f
    logical, dimension(:), allocatable :: arr_f
    integer(c_int) :: n_elements

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(arr)
    M_CHECK_NON_NULL(filename)
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(orig_shape)

    n_elements = size(arr, kind=c_int)

    ! Convert filename from c_char array to Fortran string
    call c_char_1d_as_string(filename, filename_f, ierr)
    if (is_err(ierr)) return

    n_elements = size(arr, kind=c_int)
    M_ALLOCATE(arr_f(n_elements))
    call c_int_as_logical(arr, arr_f)

    call serialize_logical_helper(arr_f, n_elements, orig_shape, filename_f, ierr)
    ! n_elements (helper arg 2) is derived from `arr` here rather than exposed as a C argument, so
    ! remap any error reported against it back onto `arr` (C-visible arg 1).
    call map_err_arg_pos(ierr, 2_c_int, 1_c_int)
end subroutine serialize_logical_nd_c
