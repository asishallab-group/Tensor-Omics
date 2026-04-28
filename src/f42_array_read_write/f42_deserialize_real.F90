#include "../macros.h"

!> Module for deserializing real arrays from files
module f42_deserialize_real
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_array_utils, only: check_file_header, REAL_TYPE_CODE
    use tox_errors, only: set_ok, is_err, validate_dimension_size, ERR_READ_DATA, set_err
    implicit none

    private
    public :: deserialize_real_1d, deserialize_real_2d, &
              deserialize_real_3d, deserialize_real_4d, deserialize_real_5d, deserialize_real_helper

contains

    !> AUTHOR_AARON_SCHROEDER
    !| Deserialize a flat real array from a file
    subroutine deserialize_real_helper(arr, n_elements, orig_shape, filename, ierr)
        integer(int32), intent(in) :: n_elements
            !! Size of `arr`
        real(real64), dimension(n_elements), intent(out) :: arr
            !! Pre-allocated array to read the data into
        integer(int32), dimension(:), intent(in) :: orig_shape
            !! Original shape of the flattened array `arr`
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit

        call set_ok(ierr)

        call validate_dimension_size(n_elements, ierr, arg_pos=2_int32)

        if (is_err(ierr)) return

        call check_file_header(filename, REAL_TYPE_CODE, orig_shape, unit, ierr)

        ! Read the entire array as a contiguous block
        read (unit, iostat=ierr) arr
        if (is_err(ierr)) call set_err(ierr, ERR_READ_DATA)

        close (unit)
    end subroutine deserialize_real_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 1D real array from a file
    subroutine deserialize_real_1d(arr, filename, ierr)
        real(real64), dimension(:), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_real_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_real_1d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 2D real array from a file
    subroutine deserialize_real_2d(arr, filename, ierr)
        real(real64), dimension(:, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_real_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_real_2d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 3D real array from a file
    subroutine deserialize_real_3d(arr, filename, ierr)
        real(real64), dimension(:, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_real_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_real_3d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 4D real array from a file
    subroutine deserialize_real_4d(arr, filename, ierr)
        real(real64), dimension(:, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_real_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_real_4d

    !> AUTHOR_AARON_SCHROEDER
    !| Directly deserialize a 5D real array from a file
    subroutine deserialize_real_5d(arr, filename, ierr)
        real(real64), dimension(:, :, :, :, :), contiguous, intent(out) :: arr
            !! Pre-allocated array to read the data into
        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ierr
            !! Error code

        call deserialize_real_helper(arr, size(arr, kind=int32), shape(arr, kind=int32), filename, ierr)
    end subroutine deserialize_real_5d

end module f42_deserialize_real

!> C binding for the subroutine to deserialize a real array from a file.
subroutine deserialize_real_nd_c(arr, orig_shape, n_dims, filename, fn_len, ierr) bind(C, name="deserialize_real_nd_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_char, c_double
    use f42_deserialize_real, only: deserialize_real_helper
    use tox_errors, only: is_err, map_err_arg_pos
    use tox_conversions, only: c_char_1d_as_string
    M_USE_NULL_VALIDATION
    implicit none

    ! Inputs / Outputs
    integer(c_int), intent(in), target :: n_dims
        !! Number of dimensions of the expected array (`size(orig_shape)`)
    integer(c_int), dimension(n_dims), intent(in), target :: orig_shape
        !! Original shape of the flattened array `arr` -> for a 2D array it would be `[n_elements, n_contained_arrays]`
    real(c_double), dimension(product(orig_shape)), intent(out), target :: arr
        !! Preallocated output array
    integer(c_int), intent(in), target :: fn_len
        !! Length of the filename
    character(kind=c_char, len=1), dimension(fn_len), intent(in), target :: filename
        !! Filename in raw bytes
    integer(c_int), intent(out), target :: ierr
        !! Error code

    ! Locals
    character(len=:), allocatable :: filename_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(arr)
    M_CHECK_NON_NULL(filename)
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(orig_shape)

    ! raw to String
    call c_char_1d_as_string(filename, filename_f, ierr)
    if (is_err(ierr)) return

    call deserialize_real_helper(arr, size(arr, kind=c_int), orig_shape, filename_f, ierr)
    call map_err_arg_pos(ierr, 2_c_int, 1_c_int)
end subroutine deserialize_real_nd_c
