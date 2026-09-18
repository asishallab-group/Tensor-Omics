#include <src/macros.h>

!> Module for serializing complex arrays into files
module f42_serde_arrays_serialize_complex
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_arrays_utils, only: write_file_header, COMPLEX_TYPE_CODE
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_WRITE_DATA, set_err
    M_IMPLICIT_NONE

    private
    public :: serialize_complex_1d, serialize_complex_2d, serialize_complex_3d, &
              serialize_complex_4d, serialize_complex_5d, serialize_complex_helper

contains
    !> M_EXPORT_C
    !| summary: Subroutine to serialize a flat complex array into a file
    !| AUTHOR_AARON_SCHROEDER
    subroutine serialize_complex_helper(arr, n_elements, arr_shape, filename, ierr)
        integer(int32), intent(in) :: n_elements
            !! Number of strings in `arr`
        complex(real64), dimension(n_elements), intent(in) :: arr
            !! Array to be serialized
        integer(int32), dimension(:), intent(in) :: arr_shape
            !! Extents of `arr`, one per dimension
        character(len=*), intent(in) :: filename
            !! Name of the file to write to
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_elements, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call write_file_header(filename, unit, COMPLEX_TYPE_CODE, size(arr_shape, kind=int32), arr_shape, ierr)
        if (is_err(ierr)) return

        ! Write the entire array as a contiguous block
        write (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_WRITE_DATA)

        close (unit)
    end subroutine serialize_complex_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 1D complex array into a file
    subroutine serialize_complex_1d(arr, filename, ierr)
        complex(real64), dimension(:), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_complex_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_complex_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 2D complex array into a file
    subroutine serialize_complex_2d(arr, filename, ierr)
        complex(real64), dimension(:, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_complex_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_complex_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 3D complex array into a file
    subroutine serialize_complex_3d(arr, filename, ierr)
        complex(real64), dimension(:, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_complex_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_complex_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 4D complex array into a file
    subroutine serialize_complex_4d(arr, filename, ierr)
        complex(real64), dimension(:, :, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_complex_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_complex_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly serialize a 5D complex array into a file
    subroutine serialize_complex_5d(arr, filename, ierr)
        complex(real64), dimension(:, :, :, :, :), contiguous, intent(in) :: arr
            !! Array to be serialized
        character(len=*), intent(in)  :: filename
            !! Name of the file to write to
        integer(int32), intent(out)   :: ierr
            !! Error code

        call serialize_complex_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine serialize_complex_5d

end module f42_serde_arrays_serialize_complex
