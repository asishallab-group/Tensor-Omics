#include <src/macros.h>

!> Module for deserializing logical arrays from files
module f42_serde_arrays_deserialize_logical
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_serde_arrays_utils, only: check_file_header, LOGICAL_TYPE_CODE
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_READ_DATA, set_err
    M_IMPLICIT_NONE

    private
    public :: deserialize_logical_1d, deserialize_logical_2d, &
              deserialize_logical_3d, deserialize_logical_4d, deserialize_logical_5d, deserialize_logical_helper

contains

    !> M_EXPORT_C
    !| summary: Deserialize a flat logical array from a file
    !| AUTHOR_AARON_SCHROEDER
    subroutine deserialize_logical_helper(arr, n_elements, arr_shape, filename, ierr)
        integer(int32), intent(in) :: n_elements
            !! Size of `arr`
        logical(c_bool), dimension(n_elements), intent(out) :: arr
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

        call validate_in_range_int(n_elements, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call check_file_header(filename, LOGICAL_TYPE_CODE, arr_shape, unit, ierr)
        if (is_err(ierr)) return

        ! Read the entire array as a contiguous block
        read (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_READ_DATA)

        close (unit)
    end subroutine deserialize_logical_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 1D logical array from a file
    subroutine deserialize_logical_1d(arr, filename, ierr)
        logical(c_bool), dimension(:), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_logical_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 2D logical array from a file
    subroutine deserialize_logical_2d(arr, filename, ierr)
        logical(c_bool), dimension(:, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_logical_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 3D logical array from a file
    subroutine deserialize_logical_3d(arr, filename, ierr)
        logical(c_bool), dimension(:, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_logical_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 4D logical array from a file
    subroutine deserialize_logical_4d(arr, filename, ierr)
        logical(c_bool), dimension(:, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_logical_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 5D logical array from a file
    subroutine deserialize_logical_5d(arr, filename, ierr)
        logical(c_bool), dimension(:, :, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_logical_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_logical_5d

end module f42_serde_arrays_deserialize_logical
