#include <src/macros.h>

!> Module for deserializing character arrays from files
module f42_serde_arrays_deserialize_char
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_arrays_utils, only: check_file_header
    use tox_errors, only: set_ok, is_err, validate_in_range_int, ERR_READ_DATA, set_err
    implicit none

    private
    public :: deserialize_char_1d, deserialize_char_2d, deserialize_char_3d, &
              deserialize_char_4d, deserialize_char_5d, deserialize_char_helper

contains
    !> AUTHOR_AARON_SCHROEDER
    !| Subroutine to deserialize a flat character array from a file
    subroutine deserialize_char_helper(arr, n_strings, orig_shape, filename, ierr)
        integer(int32), intent(in) :: n_strings
            !! Number of strings in `arr`
        character(len=*), dimension(n_strings), intent(out) :: arr
            !! Pre-allocated array to read the data into
        integer(int32), dimension(:), intent(in) :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_in_range_int(n_strings, ierr, min=0_int32, arg_pos=1_int32)

        if (is_err(ierr)) return

        call check_file_header(filename, len(arr, kind=int32), orig_shape, unit, ierr)

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

        call deserialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
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

        call deserialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
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

        call deserialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
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

        call deserialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
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

        call deserialize_char_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_char_5d

end module f42_serde_arrays_deserialize_char

!> C binding for the subroutine to deserialize a flat character array from a file.
subroutine deserialize_char_nd_c(strings, string_len, orig_shape, n_dims, &
                                 filename, fn_len, ierr) bind(C, name="deserialize_char_nd_c")
    use, intrinsic :: iso_c_binding, only: c_char, c_int
    use f42_serde_arrays_deserialize_char, only: deserialize_char_helper
    use tox_conversions, only: c_char_1d_as_string, string_as_c_char_2d
    use tox_errors, only: is_err, ERR_ALLOC_FAIL, set_err, map_err_arg_pos
    M_USE_NULL_VALIDATION
    implicit none

    ! Arguments
    integer(c_int), intent(in), target :: string_len
        !! Length of each character string
    integer(c_int), intent(in), target :: n_dims
        !! Number of dimensions of the expected array (`size(orig_shape)`)
    integer(c_int), dimension(n_dims), intent(in), target :: orig_shape
        !! Original shape of the flattened array `strings` -> for a 1D array of 7 strings it would be `[7]` with `n_dims=1`
    character(kind=c_char, len=1), dimension(string_len, product(orig_shape)), intent(out), target :: strings
        !! Output array of c_chars (2D: string_len x n_strings)
    integer(c_int), intent(in), target :: fn_len
        !! Length of the filename
    character(kind=c_char, len=1), dimension(fn_len), intent(in), target  :: filename
        !! Filename in raw bytes
    integer(c_int), intent(out), target :: ierr
        !! Error code

    character(len=:), allocatable :: filename_f
    character(len=:), dimension(:), allocatable :: strings_f
    integer(c_int) :: n_strings


    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(string_len)
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(strings)
    M_CHECK_NON_NULL(filename)
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(orig_shape)

    n_strings = size(strings, dim=2, kind=c_int)

    ! Convert filename from c_char array to Fortran string
    call c_char_1d_as_string(filename, filename_f, ierr)
    if (is_err(ierr)) return

    ! Deserialize
    M_ALLOCATE(character(len=string_len) :: strings_f(n_strings))
    call deserialize_char_helper(strings_f, n_strings, orig_shape, filename_f, ierr)
    call map_err_arg_pos(ierr, 2_c_int, 1_c_int)
    if (is_err(ierr)) return

    call string_as_c_char_2d(strings_f, strings)
end subroutine deserialize_char_nd_c
