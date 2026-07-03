#include "../macros.h"

!> Module providing serialization and deserialization routines for character arrays
!! of up to 5 dimensions, arrays are serialized to a custom binary format with a magic number and type/dimension metadata.

module f42_serialize_char
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use f42_array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
           serialize_char_4d, serialize_char_5d, serialize_char_nd

  integer(int32), parameter :: ARRAY_TYPE_CHAR = 3

contains

  !> category: C-interface
  !| Serialize a 1D character array to a binary file.
  !| The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_1d(arr, filename, ierr)
    character(len=*), contiguous, intent(in) :: arr(:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr 
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(1)
    integer(int32) :: ioerror
    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 1_int32, dims, ierr, clen)
    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if(.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_2d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(2)
    integer(int32) :: ioerror
    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 2_int32, dims, ierr, clen)

    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if(.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 3D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_3d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit, clen
    integer(int32) :: dims(3)
    dims = shape(arr)
    clen = len(arr)
    
    call set_ok(ierr)
    call set_ok(ioerror)
    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 3_int32, dims, ierr, clen)
    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 4D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_4d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(4)
    integer(int32) :: ioerror

    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 4_int32, dims, ierr, clen)
    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 5D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_5d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(5)
    integer(int32) :: ioerror
    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 5_int32, dims, ierr, clen)

    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a character array of arbitrary dimensions to a binary file.
  !| The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !| @note This routine is only called by R and serializes only flat character arrays to the memory
  subroutine serialize_char_nd(flat, flat_shape, filename, ierr)
    implicit none
    character(len=*), intent(in) :: flat(:)
    !! flat array to save
    integer(int32), intent(in) :: flat_shape(:)
    !! dimensions of the array
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: ioerror
    integer(int32) :: unit

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, size(flat_shape), flat_shape, ierr, len(flat))
    if(.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) flat
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if

    close(unit)
  end subroutine serialize_char_nd

end module f42_serialize_char
