#include <src/macros.h>
#include <src/f42/serde/macros.h>

!> Module for array utilities.
!|
!| Defines the shared on-disk binary layout used by all typed array
!| serialize/deserialize modules (int/real/complex/logical/char) and the
!| header read/write/validate helpers that implement it. The file header is
!| a fixed sequence of unformatted stream records, written and read in this
!| order: magic number ([[f42_serde_arrays_utils(module):ARRAY_FILE_MAGIC(variable)]]),
!| type code, number of dimensions `ndim`, then `ndim` dimension sizes. The
!| raw array payload follows immediately after the header, written as one
!| contiguous block by the type-specific serializers.
module f42_serde_arrays_utils
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_serde_utils
    use tox_errors
    M_IMPLICIT_NONE

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

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Opens `filename`, reads its header via [[f42_serde_arrays_utils(module):read_file_header_helper(subroutine)]],
    !| and validates the stored type code and shape against the caller's expectations. On any
    !| mismatch the unit is closed and an error is set; on success the unit is left open (positioned
    !| right after the header, ready for the caller to read the payload) and it is the caller's
    !| responsibility to close it.
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

        ! Compare against the expected magic by storing the difference straight into ierr: this is
        ! zero (== ERR_OK) only when the header matches, so the M_CHECK_IO_ERR check below (which
        ! tests is_err(ierr)) doubles as the magic-number comparison without a separate if-block.
        ierr = magic - ARRAY_FILE_MAGIC
        M_CHECK_IO_ERR(ERR_INVALID_FORMAT)

        read (unit, iostat=ierr) type_code
        M_CHECK_IO_ERR(ERR_READ_TYPE)

        read (unit, iostat=ierr) ndims
        M_CHECK_IO_ERR(ERR_READ_NDIMS)

        ! Sanity-cap ndims before allocating dims(ndims): a corrupt or non-array file could yield an
        ! arbitrary/negative value here, so reject anything outside a generous but bounded range
        ! rather than risking a huge or invalid allocation request.
        if (ndims < 0 .or. ndims > 15) then
            call set_err(ierr, ERR_INVALID_FORMAT)
            close (unit)
            return
        end if

        M_ALLOCATE(dims(ndims))
        read (unit, iostat=ierr) dims
        M_CHECK_IO_ERR(ERR_READ_DIMS)
    end subroutine read_file_header_helper

    !> M_EXPORT_C
    !| summary: Get the metadata of an array file
    !| AUTHOR_AARON_SCHROEDER
    subroutine get_array_metadata(filename, dims_out, dims_out_capacity, ndims, type_code, ierr)

        character(len=*), intent(in) :: filename
            !! Name of the file
        integer(int32), intent(out) :: ndims
            !! number of dimensions
        integer(int32), intent(in) :: dims_out_capacity
            !! Capacity of the dims_out array
        integer(int32), intent(out) :: dims_out(dims_out_capacity)
            !! Array to store output dimensions
            !! DM_RESULT_SIZE_IS(ndims)
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
