#include <src/macros.h>

!> Module for serializing character arrays into files
module f42_serde_arrays_serialize_char
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_arrays_utils, only: write_file_header
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_WRITE_DATA, set_err
    M_IMPLICIT_NONE

    private
    public :: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
              serialize_char_4d, serialize_char_5d, serialize_char_helper

contains
    !> M_EXPORT_C
    !| summary: Subroutine to serialize a flat character array into a file
    !| AUTHOR_AARON_SCHROEDER
    subroutine serialize_char_helper(arr, n_strings, arr_shape, filename, ierr)
        integer(int32), intent(in) :: n_strings
            !! Number of strings in `arr`
        character(len=*), dimension(n_strings), intent(in) :: arr
            !! Array to be serialized
        integer(int32), dimension(:), intent(in) :: arr_shape
            !! Extents of `arr`, one per dimension
        character(len=*), intent(in) :: filename
            !! Name of the file to write to
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_strings, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call write_file_header(filename, unit, len(arr, kind=int32), size(arr_shape, kind=int32), arr_shape, ierr)
        if (is_err(ierr)) return

        ! Write the entire array as a contiguous block
        write (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_WRITE_DATA)

        close (unit)
    end subroutine serialize_char_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 1D character array into a file
    subroutine serialize_char_1d(arr, filename, ierr)
        character(len=*), dimension(:), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_char_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 2D character array into a file
    subroutine serialize_char_2d(arr, filename, ierr)
        character(len=*), dimension(:, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_char_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 3D character array into a file
    subroutine serialize_char_3d(arr, filename, ierr)
        character(len=*), dimension(:, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_char_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 4D character array into a file
    subroutine serialize_char_4d(arr, filename, ierr)
        character(len=*), dimension(:, :, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_char_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 5D character array into a file
    subroutine serialize_char_5d(arr, filename, ierr)
        character(len=*), dimension(:, :, :, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_char_5d

end module f42_serde_arrays_serialize_char
