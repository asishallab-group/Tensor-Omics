#include <src/macros.h>
#include <src/f42/serde/macros.h>

!> Module for array utilities.
module f42_serde_arrays_utils
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_utils
    use tox_errors
    implicit none

    public :: get_array_metadata, read_file_header_helper, write_file_header

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex
        !! Magic number for array files
contains

    !> AUTHOR_AARON_SCHROEDER
    !| Opens unit and writes fileheader with all metadata to the given filename
    subroutine write_file_header(filename, unit, type_code, ndim, dims, ierr)
        character(len=*), intent(in) :: filename
            !! filename to write to
        integer(int32), intent(in) :: type_code
            !! type code of the array
            !!
            !! |    Type   |           Code          |
            !! |-----------|-------------------------|
            !! |  integer  |   M_INTEGER_TYPE_CODE  |
            !! |    real   |   M_REAL_TYPE_CODE     |
            !! |  complex  |   M_COMPLEX_TYPE_CODE  |
            !! |  logical  |   M_LOGICAL_TYPE_CODE  |
            !! | character |      string length      |
            !!
        integer(int32), intent(in) :: ndim
            !! number of dimensions
        integer(int32), intent(in) :: dims(ndim)
            !! dimensions of the array
        integer(int32), intent(inout) :: ierr
            !! error code
        integer(int32), intent(out) :: unit
            !! Fortran unit number for the file

        call set_ok(ierr)

        open (newunit=unit, file=filename, form='unformatted', access='stream', status='replace', iostat=ierr)
        M_CHECK_IO_ERR(ERR_FILE_OPEN)

        write (unit, iostat=ierr) ARRAY_FILE_MAGIC
        M_CHECK_IO_ERR(ERR_WRITE_MAGIC)

        write (unit, iostat=ierr) type_code
        M_CHECK_IO_ERR(ERR_WRITE_TYPE)

        write (unit, iostat=ierr) ndim
        M_CHECK_IO_ERR(ERR_WRITE_NDIMS)

        write (unit, iostat=ierr) dims
        M_CHECK_IO_ERR(ERR_WRITE_DIMS)
    end subroutine write_file_header

    subroutine check_file_header(filename, expected_type_code, expected_shape, unit, ierr)
        character(len=*), intent(in) :: filename
            !! filename to read from
        integer(int32), intent(in) :: expected_type_code
            !! Expected type code in header
            !!
            !! |    Type   |           Code          |
            !! |-----------|-------------------------|
            !! |  integer  |   M_INTEGER_TYPE_CODE  |
            !! |    real   |   M_REAL_TYPE_CODE     |
            !! |  complex  |   M_COMPLEX_TYPE_CODE  |
            !! |  logical  |   M_LOGICAL_TYPE_CODE  |
            !! | character |  expected string length |
            !!
        integer(int32), dimension(:), intent(in) :: expected_shape
            !! Expected shape in header
        integer(int32), intent(out) :: unit
            !! Fortran unit number for the file
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: type_code, ndim
        integer(int32), dimension(:), allocatable :: dims

        call set_ok(ierr)
        ! open file and read header
        call read_file_header_helper(filename, unit, type_code, ndim, dims, ierr)
        if (is_err(ierr)) return

        call validate_type_code(type_code, expected_type_code, unit, ierr)
        if (ndim /= size(expected_shape, kind=int32)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
        else if (any(expected_shape /= dims)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
        end if
        if (is_err(ierr)) then
            close(unit)
        end if
    end subroutine check_file_header

    !> AUTHOR_AARON_SCHROEDER
    !| Opens unit and reads file header with all metadata from given file
    subroutine read_file_header_helper(filename, unit, type_code, ndims, dims, ierr)
        character(len=*), intent(in) :: filename
            !! filename to read from
        integer(int32), intent(out) :: unit
            !! Fortran unit number for the file
        integer(int32), intent(out) :: type_code
            !! type code of the array
            !!
            !! |    Type   |           Code          |
            !! |-----------|-------------------------|
            !! |  integer  |   M_INTEGER_TYPE_CODE  |
            !! |    real   |   M_REAL_TYPE_CODE     |
            !! |  complex  |   M_COMPLEX_TYPE_CODE  |
            !! |  logical  |   M_LOGICAL_TYPE_CODE  |
            !! | character |      string length      |
            !!
        integer(int32), intent(out) :: ndims
            !! number of dimensions
        integer(int32), intent(out) :: ierr
            !! error code
        integer(int32), allocatable, intent(out) :: dims(:)
            !! dimensions of the array

        integer(int32) :: magic

        call set_ok(ierr)

        open (newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ierr)
        M_CHECK_IO_ERR(ERR_FILE_OPEN)

        read (unit, iostat=ierr) magic
        M_CHECK_IO_ERR(ERR_READ_MAGIC)

        ! check magic (error if not same)
        ierr = magic - ARRAY_FILE_MAGIC
        M_CHECK_IO_ERR(ERR_INVALID_FORMAT)

        read (unit, iostat=ierr) type_code
        M_CHECK_IO_ERR(ERR_READ_TYPE)

        read (unit, iostat=ierr) ndims
        M_CHECK_IO_ERR(ERR_READ_NDIMS)

        if (ndims < 0 .or. ndims > 15) then
            call set_err(ierr, ERR_INVALID_FORMAT)
            close (unit)
            return
        end if

        M_ALLOCATE(dims(ndims))
        read (unit, iostat=ierr) dims
        M_CHECK_IO_ERR(ERR_READ_DIMS)
    end subroutine read_file_header_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Get the metadata of an array file
    subroutine get_array_metadata(filename, dims_out, dims_out_capacity, ndims, type_code, ierr)

        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ndims
            !! number of dimensions
        integer(int32), intent(in) :: dims_out_capacity
            !! Capacity of the dims_out array
        integer(int32), intent(out) :: dims_out(dims_out_capacity)
            !! Array to store output dimensions
        integer(int32), intent(out) :: type_code
            !! Type code of the serialized array
            !!
            !!
            !! |    Type   |           Code          |
            !! |-----------|-------------------------|
            !! |  integer  |   M_INTEGER_TYPE_CODE  |
            !! |    real   |   M_REAL_TYPE_CODE     |
            !! |  complex  |   M_COMPLEX_TYPE_CODE  |
            !! |  logical  |   M_LOGICAL_TYPE_CODE  |
            !! | character |      string length      |
            !!
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: unit
        integer(int32), allocatable :: dims(:)

        call read_file_header_helper(filename, unit, type_code, ndims, dims, ierr)
        close (unit)
        if (is_err(ierr)) return

        if (size(dims) > dims_out_capacity) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if

        dims_out = 1
        dims_out(1:ndims) = dims
    end subroutine get_array_metadata
end module f42_serde_arrays_utils

!> C binding for the subroutine to get the dimensions of an array file
subroutine get_array_metadata_c(filename, fn_len, dims_out, dims_out_capacity, ndims, type_code, ierr) bind(C, name="get_array_metadata_c")
    use iso_c_binding, only: c_int, c_char
    use f42_serde_arrays_utils, only: get_array_metadata
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION
    implicit none

    ! Input
    integer(c_int), intent(in), target :: fn_len
        !! Length of the filename array
    character(kind=c_char, len=1), intent(in), target :: filename(fn_len)
        !! Array of ASCII characters representing the filename
    integer(c_int), intent(in), target :: dims_out_capacity

    ! Output
    integer(c_int), intent(out), target :: dims_out(dims_out_capacity)
        !! Output array for dimensions
    integer(c_int), intent(out), target :: ndims
        !! Output variable for the number of dimensions
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), intent(out), target :: type_code
        !! Type code of the serialized array
        !!
        !!
        !! |    Type   |           Code          |
        !! |-----------|-------------------------|
        !! |  integer  |   M_INTEGER_TYPE_CODE  |
        !! |    real   |   M_REAL_TYPE_CODE     |
        !! |  complex  |   M_COMPLEX_TYPE_CODE  |
        !! |  logical  |   M_LOGICAL_TYPE_CODE  |
        !! | character |      string length      |
        !!

    ! Local variables
    character(len=:), allocatable :: filename_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(filename)
    M_CHECK_NON_NULL(dims_out)
    M_CHECK_NON_NULL(ndims)
    M_CHECK_NON_NULL(type_code)
    M_CHECK_NON_NULL(dims_out_capacity)

    call c_char_1d_as_string(filename, filename_f, ierr)
    if (is_err(ierr)) return

    call get_array_metadata(filename_f, dims_out, dims_out_capacity, ndims, type_code, ierr)
end subroutine get_array_metadata_c
