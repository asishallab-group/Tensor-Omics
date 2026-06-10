#include "../macros.h"

!> Module for deserializing character arrays from files
module f42_deserialize_char
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use f42_array_utils, only : read_file_header, check_okay_ioerror, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_char_1d, deserialize_char_2d, deserialize_char_3d, deserialize_char_nd, &
          deserialize_char_4d, deserialize_char_5d

contains
  !> category: C-interface
  !| Subroutine to deserialize a flat character array from a file
  subroutine deserialize_char_nd(flat, filename, ierr)
    character(len=*), intent(out) :: flat(:)
      !! Output flat character array
    character(len=*), intent(in) :: filename
      !! Name of the file to read
    integer(int32), intent(out) :: ierr
      !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 1_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) flat
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_nd

  !> Directly deserialize a 1D character array from a file (array already allocated)
  subroutine deserialize_char_1d(arr, filename, ierr)
    character(len=*), contiguous, intent(out) :: arr(:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in)  :: filename
    !! Name of the file to read from
    integer(int32), intent(out)   :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 1_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_1d

  !> Directly deserialize a 2D character array from a file (array already allocated)
  subroutine deserialize_char_2d(arr, filename, ierr)
    character(len=*), intent(out) :: arr(:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in)  :: filename
    !! Name of the file
    integer(int32), intent(out)   :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 2_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_2d

  !> Directly deserialize a 3D character array from a file (array already allocated)
  subroutine deserialize_char_3d(arr, filename, ierr)
    character(len=*), intent(out) :: arr(:,:,:)
    !!Pre-allocated array to read the data into
    character(len=*), intent(in)  :: filename
    !! Name of the file
    integer(int32), intent(out)   :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 3_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_3d

  !> Directly deserialize a 4D character array from a file (array already allocated)
  subroutine deserialize_char_4d(arr, filename, ierr)
    character(len=*), intent(out) :: arr(:,:,:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in)  :: filename
    !! Name of the file
    integer(int32), intent(out)   :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 4_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_4d

  !> Directly deserialize a 5D character array from a file (array already allocated)
  subroutine deserialize_char_5d(arr, filename, ierr)
    character(len=*), intent(out) :: arr(:,:,:,:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in)  :: filename
    !! Name of the file
    integer(int32), intent(out)   :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 5_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_5d
  
end module f42_deserialize_char
