#include <src/macros.h>

!> Module for deserializing character arrays from files
module f42_serde_arrays_deserialize_char
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_arrays_utils, only: check_file_header
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_READ_DATA, set_err
    M_IMPLICIT_NONE

    private
    public :: deserialize_char_1d, deserialize_char_2d, deserialize_char_3d, &
              deserialize_char_4d, deserialize_char_5d, deserialize_char_helper

contains
    !> M_EXPORT_C
    !| summary: Subroutine to deserialize a flat character array from a file
    !| AUTHOR_AARON_SCHROEDER
    subroutine deserialize_char_helper(arr, n_strings, strlen, arr_shape, filename, ierr)
        integer(int32), intent(in) :: n_strings
            !! Number of strings in `arr`
        integer(int32), intent(in) :: strlen
            !! Length of each string in `arr`.
            !! DM_OUTPUT_FROM(type_code, get_array_metadata, f42_serde_arrays_utils, AUTO)
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! | dims_out_capacity | 5_int32     |
        character(len=strlen), dimension(n_strings), intent(out) :: arr
            !! Pre-allocated array to read the data into
        integer(int32), dimension(:), intent(in) :: arr_shape
            !! Extents of `arr`, one per dimension
            !! DM_OUTPUT_FROM(dims_out, get_array_metadata, f42_serde_arrays_utils, AUTO)
            !!
            !! | Producer input    | Supplied by |
            !! |-------------------|-------------|
            !! | dims_out_capacity | 5_int32     |
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_strings, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call check_file_header(filename, len(arr, kind=int32), arr_shape, unit, ierr)
        if (is_err(ierr)) return

        ! Read the entire array as a contiguous block
        read (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_READ_DATA)

        close (unit)
    end subroutine deserialize_char_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 1D character array from a file (array already allocated)
    subroutine deserialize_char_1d(arr, filename, ierr)
        character(len=*), dimension(:), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in)  :: filename
            !! Name of the file to read from
        integer(int32), intent(out)   :: ierr
            !! Error code

        call deserialize_char_helper(arr, size(arr, kind=int32), len(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_char_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 2D character array from a file (array already allocated)
    subroutine deserialize_char_2d(arr, filename, ierr)
        character(len=*), dimension(:, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in)  :: filename
            !! Name of the file to read from
        integer(int32), intent(out)   :: ierr
            !! Error code

        call deserialize_char_helper(arr, size(arr, kind=int32), len(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_char_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 3D character array from a file (array already allocated)
    subroutine deserialize_char_3d(arr, filename, ierr)
        character(len=*), dimension(:, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in)  :: filename
            !! Name of the file to read from
        integer(int32), intent(out)   :: ierr
            !! Error code

        call deserialize_char_helper(arr, size(arr, kind=int32), len(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_char_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 4D character array from a file (array already allocated)
    subroutine deserialize_char_4d(arr, filename, ierr)
        character(len=*), dimension(:, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in)  :: filename
            !! Name of the file to read from
        integer(int32), intent(out)   :: ierr
            !! Error code

        call deserialize_char_helper(arr, size(arr, kind=int32), len(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_char_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 5D character array from a file (array already allocated)
    subroutine deserialize_char_5d(arr, filename, ierr)
        character(len=*), dimension(:, :, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in)  :: filename
            !! Name of the file to read from
        integer(int32), intent(out)   :: ierr
            !! Error code

        call deserialize_char_helper(arr, size(arr, kind=int32), len(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_char_5d

end module f42_serde_arrays_deserialize_char
