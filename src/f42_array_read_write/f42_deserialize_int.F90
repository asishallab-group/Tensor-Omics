#include "../macros.h"

!> Module for deserializing integer arrays from files
module f42_deserialize_int
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use f42_array_utils, only: read_file_header, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_int_1d, deserialize_int_2d, &
           deserialize_int_3d, deserialize_int_4d, deserialize_int_5d, deserialize_int_flat

contains

  !> Deserializes any array inot a flat array
  subroutine deserialize_int_flat(arr, filename, ierr)
    integer(int32), intent(out) :: arr(:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in) :: filename
    !! Name of the file
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_flat

  !> Deserialize a flat integer array from a file
  !> Directly deserialize a 1D integer array from a file
  subroutine deserialize_int_1d(arr, filename, ierr)
    integer(int32), intent(out) :: arr(:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in) :: filename
    !! Name of the file
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 1_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_1d

  !> Directly deserialize a 2D integer array from a file
  subroutine deserialize_int_2d(arr, filename, ierr)
    integer(int32), intent(out) :: arr(:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in) :: filename
    !! Name of the file
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 2_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_2d

  !> Directly deserialize a 3D integer array from a file
  subroutine deserialize_int_3d(arr, filename, ierr)
    integer(int32), intent(out) :: arr(:,:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in) :: filename
    !! Name of the file
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 3_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_3d

  !> Directly deserialize a 4D integer array from a file
  subroutine deserialize_int_4d(arr, filename, ierr)
    integer(int32), intent(out) :: arr(:,:,:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in) :: filename
    !! Name of the File
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 4_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_4d

  !> Directly deserialize a 5D integer array from a file
  subroutine deserialize_int_5d(arr, filename, ierr)
    integer(int32), intent(out) :: arr(:,:,:,:,:)
    !! Pre-allocated array to read the data into
    character(len=*), intent(in) :: filename
    !! Name of the file
    integer(int32), intent(out) :: ierr
    !! Error code

    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)

    call set_ok(ierr)
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 5_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_5d

end module f42_deserialize_int
