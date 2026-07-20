#include <src/macros.h>

!> Module for deserializing integer arrays from files
module f42_serde_arrays_deserialize_int
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_arrays_utils, only: check_file_header, INTEGER_TYPE_CODE
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_READ_DATA, set_err
    implicit none

    private
    public :: deserialize_int_1d, deserialize_int_2d, &
              deserialize_int_3d, deserialize_int_4d, deserialize_int_5d, deserialize_int_helper

contains

    !> M_EXPORT_C
    !| summary: Deserialize a flat integer array from a file
    !| AUTHOR_AARON_SCHROEDER
    subroutine deserialize_int_helper(arr, n_elements, orig_shape, filename, ierr)
        integer(int32), intent(in) :: n_elements
            !! Size of `arr`
        integer(int32), dimension(n_elements), intent(out) :: arr
            !! Pre-allocated array to read the data into
        integer(int32), dimension(:), intent(in) :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_elements, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call check_file_header(filename, INTEGER_TYPE_CODE, orig_shape, unit, ierr)
        if (is_err(ierr)) return

        ! Read the entire array as a contiguous block
        read (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_READ_DATA)

        close (unit)
    end subroutine deserialize_int_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 1D integer array from a file
    subroutine deserialize_int_1d(arr, filename, ierr)
        integer(int32), dimension(:), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_int_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_int_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 2D integer array from a file
    subroutine deserialize_int_2d(arr, filename, ierr)
        integer(int32), dimension(:, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_int_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_int_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 3D integer array from a file
    subroutine deserialize_int_3d(arr, filename, ierr)
        integer(int32), dimension(:, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_int_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_int_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 4D integer array from a file
    subroutine deserialize_int_4d(arr, filename, ierr)
        integer(int32), dimension(:, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_int_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_int_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 5D integer array from a file
    subroutine deserialize_int_5d(arr, filename, ierr)
        integer(int32), dimension(:, :, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_int_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_int_5d

end module f42_serde_arrays_deserialize_int
