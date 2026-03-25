#include "../macros.h"

!> Module for array utilities.

module f42_array_utils
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors
    implicit none

    public :: get_array_metadata, read_file_header, write_file_header, check_okay_ndims

    integer(int32), parameter :: ARRAY_FILE_MAGIC = int(z'46413230', int32) ! 'FA20' in hex
    !! Magic number for array files

   contains

  !> Check I/O error and set ierr accordingly
  subroutine check_okay_ioerror(ioerror, ierr, msg, unit)
    
    integer(int32), intent(in) :: ioerror
    !! IO error set by fortran
    integer(int32), intent(out) :: ierr
    !! Error code
    integer(int32), intent(in) :: msg
    !! Error code readable version used for setting
    integer(int32), intent(in) :: unit
    !! pass unit allowing it to be closed 

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, msg)
      close(unit)
      return
    end if
    
  end subroutine

  subroutine check_okay_ndims(ndims, expected, unit, ierr)
    integer(int32), intent(in):: ndims
    !! number of dimensions
    integer(int32), intent(in) :: expected
    !! expected number of dimensions
    integer(int32), intent(in) :: unit 
    !! unit representation
    integer(int32), intent(out) :: ierr 
    !! error code

    if(ndims /= expected) then
      close(unit)
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if

  end Subroutine


  !> Opens unit and writes fileheader with all metadata to the given filename
  subroutine write_file_header(filename, unit, type_code, ndim, dims, ierr, clen)
    character(len=*), intent(in) :: filename
    !! filename to write to
    integer(int32), intent(in) :: type_code
    !! type code of the array (1=int, 2=real, 3=char)
    integer(int32), intent(in) :: ndim
    !! number of dimensions
    integer(int32), intent(in) :: dims(ndim)
    !! dimensions of the array
    integer(int32), intent(in), optional :: clen
    !! character length (only for character arrays)
    integer(int32), intent(inout) :: ierr
    !! error code
    integer(int32), INTENT(OUT) :: unit
    !! Fortran unit number for the file
    integer(int32) :: ioerror
    !! internal I/O error code
    
    call set_ok(ierr)
    call set_ok(ioerror)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='replace', iostat=ioerror)

    call check_okay_ioerror(ioerror, ierr, ERR_FILE_OPEN, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) ARRAY_FILE_MAGIC
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_MAGIC, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) type_code
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_TYPE, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) ndim
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_NDIMS, unit)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) dims
    call check_okay_ioerror(ioerror, ierr, ERR_WRITE_DIMS, unit)
    if (.not. is_ok(ierr)) return

    if (type_code == 3 .and. present(clen)) then
      write(unit, iostat=ioerror) clen
      call check_okay_ioerror(ioerror, ierr, ERR_WRITE_CHARLEN, unit)
      if (.not. is_ok(ierr)) return
    end if
  end subroutine write_file_header

  !> Opens unit and reads file header with all metadata from given file
  subroutine read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    character(len=*), intent(in) :: filename
    !! filename to read from
    integer(int32), intent(out) :: unit
    !! Fortran unit number for the file
    integer(int32), intent(out) :: type_code
    !! type code of the array (1=int, 2=real, 3=char, 4=logical, 5=complex)
    integer(int32), intent(out) :: ndims
    !! number of dimensions
    integer(int32), intent(out) :: clen
    !! character length (only for character arrays)
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32), allocatable :: dims(:)
    !! dimensions of the array
    integer(int32) :: ioerror
    !! internal I/O error code

    integer(int32) :: magic
    !! magic number to verify file format


    call set_ok(ioerror)
    call set_ok(ierr)

    open(newunit=unit, file=filename, form='unformatted', access='stream', status='old', iostat=ioerror)
    call check_okay_ioerror(ioerror, ierr, ERR_FILE_OPEN, unit)
    if (.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) magic
    call check_okay_ioerror(ioerror, ierr, ERR_READ_MAGIC, unit)
    if (.not. is_ok(ierr)) return

    if (magic /= ARRAY_FILE_MAGIC) then
      call set_err_once(ierr, ERR_INVALID_FORMAT)
      close(unit)
      return
    end if

    read(unit, iostat=ioerror) type_code
    call check_okay_ioerror(ioerror, ierr, ERR_READ_TYPE, unit)
    if (.not. is_ok(ierr)) then
      close(unit)
      return
    end if

    read(unit, iostat=ioerror) ndims
    call check_okay_ioerror(ioerror, ierr, ERR_READ_NDIMS, unit)
    if (.not. is_ok(ierr)) then
      close(unit)
      return
    end if

    allocate(dims(ndims))
    read(unit, iostat=ioerror) dims
    call check_okay_ioerror(ioerror, ierr, ERR_READ_DIMS, unit)
    if (.not. is_ok(ierr)) then
      close(unit)
      return
    end if

    if(type_code==3) then
      read(unit, iostat=ioerror) clen
      call check_okay_ioerror(ioerror, ierr, ERR_READ_CHARLEN, unit)
      if (.not. is_ok(ierr)) then
        close(unit)
        return
      end if
    else
      clen = 0 ! Not applicable for non-character types
    end if

  end subroutine

  !> Get the metadata of an array file
  subroutine get_array_metadata(filename, dims_out, dims_out_capacity, ndims, ierr, clen)
    implicit none

    character(len=*), intent(in) :: filename
    !! Name of the file
    integer(int32), intent(out) :: ndims
    !! number of dimensions
    integer(int32), intent(in) :: dims_out_capacity
    !! Capacity of the dims_out array
    integer(int32), intent(out) :: dims_out(dims_out_capacity)
    !! Array to store output dimensions
    integer(int32), intent(out) :: ierr
    !! Error code
    integer(int32), intent(out), optional :: clen
    !! length of each string (needed for char arrays)

    integer(int32) :: unit
    integer(int32), allocatable :: dims(:)
    integer(int32) :: type_code
    integer(int32) :: i
    integer(int32) :: ioerror
    integer(int32) :: local_clen 

    call set_ok(ioerror)

    call read_file_header(filename, unit, type_code, ndims, dims, local_clen, ierr)
    close(unit)
    if (.not. is_ok(ierr)) return

    if(size(dims) > dims_out_capacity) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if

    do i = 1, ndims
      dims_out(i) = dims(i)
    end do

    if (present(clen)) then
      clen = local_clen
    end if
  end subroutine
end module f42_array_utils




!> C binding for the subroutine to get the dimensions of an array file
subroutine get_array_metadata_C(filename_raw, fn_len, dims_out, dims_out_capacity, ndims, ierr, clen) bind(C, name="get_array_metadata_C")
  use iso_c_binding, only: c_int, c_char
  use iso_fortran_env, only : int32
  use f42_array_utils, only : get_array_metadata
  use tox_conversions, only : c_char_1d_as_string
  use tox_errors, only : set_ok, is_ok
  M_USE_NULL_VALIDATION
  implicit none

  ! Input
  integer(c_int), intent(in), target :: fn_len
    !! Length of the filename array
  character(kind=c_char, len=1), intent(in), target :: filename_raw(fn_len)
    !! Array of ASCII characters representing the filename
  integer(c_int), intent(in), target :: dims_out_capacity
  
  ! Output
  integer(c_int), intent(out), target :: dims_out(dims_out_capacity)
    !! Output array for dimensions
  integer(c_int), intent(out), target :: ndims
    !! Output variable for the number of dimensions
  integer(c_int), intent(out), target :: ierr
    !! Error code
  integer(c_int), intent(out), target :: clen
    !! Character length (only for character arrays)

  ! Local variables
  character(len=:), allocatable :: filename
   !! Filename as a string

  M_CHECK_IERR_NON_NULL
  M_CHECK_NON_NULL(fn_len)
  M_CHECK_NON_NULL(filename_raw)
  M_CHECK_NON_NULL(dims_out)
  M_CHECK_NON_NULL(ndims)
  M_CHECK_NON_NULL(clen)
  

  call c_char_1d_as_string(filename_raw, filename, ierr)
  if( .not. is_ok(ierr)) return

  call get_array_metadata(filename, dims_out, dims_out_capacity, ndims, ierr, clen)
end subroutine get_array_metadata_C