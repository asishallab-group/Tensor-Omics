#include "../macros.h"

!> Module for deserializing real (double precision) arrays from binary files
module f42_deserialize_real
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use f42_array_utils, only: read_file_header, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_real_1d, deserialize_real_2d, &
           deserialize_real_3d, deserialize_real_4d, deserialize_real_5d, deserialize_real_flat

contains

  !> Deserializes any array into a flat real array
  subroutine deserialize_real_flat(arr, filename, ierr)
    real(real64), intent(out) :: arr(:)
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

    call validate_type_code(type_code, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_flat

  !> Directly deserialize a 1D real array from a file (array already allocated)
  subroutine deserialize_real_1d(arr, filename, ierr)
    real(real64), intent(out) :: arr(:)
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

    call validate_type_code(type_code, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 1_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_1d


  !> Directly deserialize a 2D real array from a file (array already allocated)
  subroutine deserialize_real_2d(arr, filename, ierr)
    real(real64), intent(out) :: arr(:,:)
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

    call validate_type_code(type_code, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 2_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_2d

  !> Directly deserialize a 3D real array from a file (array already allocated)
  subroutine deserialize_real_3d(arr, filename, ierr)
    real(real64), intent(out) :: arr(:,:,:)
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

    call validate_type_code(type_code, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 3_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_3d


  !> Directly deserialize a 4D real array from a file (array already allocated)
  subroutine deserialize_real_4d(arr, filename, ierr)
    real(real64), intent(out) :: arr(:,:,:,:)
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

    call validate_type_code(type_code, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 4_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_4d


  !> Directly deserialize a 5D real array from a file (array already allocated)
  subroutine deserialize_real_5d(arr, filename, ierr)
    real(real64), intent(out) :: arr(:,:,:,:,:)
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

    call validate_type_code(type_code, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 5_int32, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_5d

end module f42_deserialize_real




!> C binding for the subroutine to deserialize a real array from a file.
!> Deserializes an array of any dimension into a flat buffer.
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_nd_C(arr, arr_size, filename_raw, fn_len, ierr) bind(C, name="deserialize_real_nd_C")
    use iso_c_binding, only : c_int, c_double, c_char
    use iso_fortran_env, only: int32, real64
    use f42_deserialize_real, only : deserialize_real_flat
    use tox_errors, only : set_ok, is_ok
    use tox_conversions, only : c_char_1d_as_string
    M_USE_NULL_VALIDATION
    implicit none

    ! Inputs / Outputs
    integer(c_int), intent(in), target :: arr_size  
    !! size of the output array
    real(c_double), intent(out), target :: arr(arr_size)
    !! output array
    integer(c_int), intent(in), target :: fn_len
    !! length of the filename
    character(kind=c_char, len=1), intent(in), target :: filename_raw(fn_len)
    !! Filename in ascii
    integer(c_int), intent(out), target :: ierr
    !! error code

    ! Locals
    character(len=:), allocatable :: filename
    !! filename

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(arr_size)
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(arr)
    M_CHECK_NON_NULL(filename_raw)

    call set_ok(ierr)

    ! ASCII to String
    call c_char_1d_as_string(filename_raw, filename, ierr)
    if (.not. is_ok(ierr)) return

    call deserialize_real_flat(arr, filename, ierr)
end subroutine deserialize_real_nd_C