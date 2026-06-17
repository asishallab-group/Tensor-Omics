#include "../macros.h"

module f42_serialize_real
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  use f42_array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_real_1d, serialize_real_2d, serialize_real_3d, &
           serialize_real_4d, serialize_real_5d, serialize_real_nd

  integer(int32), parameter :: ARRAY_TYPE_REAL = 2
contains

  !> Serialize a 1D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_1d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! output filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: unit
    integer(int32) :: dims(1)
    integer(int32) :: ioerror
    dims = shape(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 1_int32, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_2d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: unit
    integer(int32) :: dims(2)
    integer(int32) :: ioerror

    call set_ok(ierr)
    call set_ok(ioerror)

    dims = shape(arr)
    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 2_int32, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 3D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_3d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit
    integer(int32) :: dims(3)

    call set_ok(ierr)
    call set_ok(ioerror)

    dims = shape(arr)
    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 3_int32, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 4D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_4d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), INTENT(OUT) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit
    integer(int32) :: dims(4)

    call set_ok(ierr)
    call set_ok(ioerror)
    dims = shape(arr)

    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 4_int32, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 5D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_5d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:,:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit
    integer(int32) :: dims(5)
    
    call set_ok(ierr)
    call set_ok(ioerror)

    dims = shape(arr)
    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 5_int32, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Writes serialized real array from R to file with metdata.
  subroutine serialize_real_nd(arr, dims, ndim, filename, ierr)
    real(real64), intent(in) :: arr(:)
      !! array to save
    integer(int32), intent(in) :: dims(:)
      !! Dimensions of the array
    integer(int32), intent(in) :: ndim
      !! Number of dimensions
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), INTENT(OUT) :: ierr
      !! error code

    integer(int32) :: ioerror, unit

    call set_ok(ierr)
    call set_ok(ioerror)

    if (size(dims) /= ndim) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if

    call write_file_header(filename, unit, ARRAY_TYPE_REAL, ndim, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine
 
end module f42_serialize_real
