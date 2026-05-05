#include "../macros.h"

!> Module for deserializing character arrays from files
module f42_serialize_char
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_array_utils, only: write_file_header
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_WRITE_DATA, set_err
    implicit none

    private
    public :: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
              serialize_char_4d, serialize_char_5d, serialize_char_helper

contains
    !> AUTHOR_AARON_SCHROEDER
    !| Subroutine to serialize a flat character array into a file
    subroutine serialize_char_helper(arr, n_strings, orig_shape, filename, ierr)
        integer(int32), intent(in) :: n_strings
            !! Number of strings in `arr`
        character(len=*), dimension(n_strings), intent(in) :: arr
            !! Array to be serialized
        integer(int32), dimension(:), intent(in) :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=*), intent(in) :: filename
            !! Name of the file to write to
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_strings, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call write_file_header(filename, unit, len(arr, kind=int32), size(orig_shape, kind=int32), orig_shape, ierr)

        ! Read the entire array as a contiguous block
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

end module f42_serialize_char

!> C binding for the subroutine to serialize a flat character array into a file.
subroutine serialize_char_nd_c(raw_chars, clen, orig_shape, n_dims, &
                                 filename, fn_len, ierr) bind(C, name="serialize_char_nd_c")
    use, intrinsic :: iso_c_binding, only: c_char, c_int
    use f42_serialize_char, only: serialize_char_helper
    use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string
    use tox_errors, only: is_err, map_err_arg_pos
    M_USE_NULL_VALIDATION
    implicit none

    ! Arguments
    integer(c_int), intent(in), target :: clen
        !! Length of each character string
    integer(c_int), intent(in), target :: n_dims
        !! Number of dimensions of the expected array (`size(orig_shape)`)
    integer(c_int), dimension(n_dims), intent(in), target :: orig_shape
        !! Original shape of the flattened array `raw_chars` -> for a 1D array of 7 strings it would be `[7]` with `n_dims=1`
    character(kind=c_char, len=1), dimension(clen, product(orig_shape)), intent(in), target :: raw_chars
        !! Input array of c_chars (2D: clen x n_strings)
    integer(c_int), intent(in), target :: fn_len
        !! Length of the filename
    character(kind=c_char, len=1), dimension(fn_len), intent(in), target  :: filename
        !! c_char array representing the filename
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: filename_f
    character(len=:), dimension(:), allocatable :: raw_chars_f
    integer(c_int) :: n_strings


    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(clen)
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(raw_chars)
    M_CHECK_NON_NULL(filename)
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(orig_shape)

    n_strings = size(raw_chars, dim=2, kind=c_int)

    ! Convert filename from c_char array to Fortran string
    call c_char_1d_as_string(filename, filename_f, ierr)
    if (is_err(ierr)) return

    call c_char_2d_as_string(raw_chars, raw_chars_f, ierr)
    if (is_err(ierr)) return

    ! serialize
    call serialize_char_helper(raw_chars_f, n_strings, orig_shape, filename_f, ierr)
    call map_err_arg_pos(ierr, 2_c_int, 1_c_int)
    if (is_err(ierr)) return

end subroutine serialize_char_nd_c
