#include "../macros.h"

!> Module for deserializing complex arrays from files
module f42_deserialize_complex
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use f42_array_utils, only: read_file_header, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_complex_1d, deserialize_complex_2d, &
           deserialize_complex_3d, deserialize_complex_4d, deserialize_complex_5d, deserialize_complex_flat

contains

  !> Deserializes any array inot a flat array
  subroutine deserialize_complex_flat(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_flat
  !> Deserialize a flat complex array from a file
  !> Directly deserialize a 1D complex array from a file
  subroutine deserialize_complex_1d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_1d

  !> Directly deserialize a 2D complex array from a file
  subroutine deserialize_complex_2d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_2d

  !> Directly deserialize a 3D complex array from a file
  subroutine deserialize_complex_3d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_3d

  !> Directly deserialize a 4D complex array from a file
  subroutine deserialize_complex_4d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:,:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 4, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_4d

  !> Directly deserialize a 5D complex array from a file
  subroutine deserialize_complex_5d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:,:,:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_5d

end module f42_deserialize_complex



!> C binding for the subroutine to deserialize a complex array from a file.
!> The array is read into a flat buffer and reshaped in the calling language.
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_complex_nd_C(arr, arr_size, filename_raw, fn_len, ierr) bind(C, name="deserialize_complex_nd_C")
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32, real64
    use f42_deserialize_complex, only : deserialize_complex_flat
    use tox_errors, only : set_ok, is_ok
    use tox_conversions, only : c_char_1d_as_string
    M_USE_NULL_VALIDATION
    implicit none

    ! Inputs / Outputs
    integer(c_int), intent(in), target :: arr_size           ! Buffer length
    !! Size of the array
    complex(real64), intent(out), target :: arr(arr_size)      ! Preallocated buffer from C/Python
    !! preallocated array
    integer(c_int), intent(in), target :: fn_len
    !! length of the filename
    character(kind=c_char, len=1), intent(in), target :: filename_raw(fn_len)
    !! Filename in raw bytes
    integer(c_int), intent(out), target :: ierr
    !! Error code

    ! Locals
    character(len=:), allocatable :: filename

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(arr_size)
    M_CHECK_NON_NULL(fn_len)
    M_CHECK_NON_NULL(arr)
    M_CHECK_NON_NULL(filename_raw)

    call set_ok(ierr)

    ! raw to String
    call c_char_1d_as_string(filename_raw, filename, ierr)
    if (.not. is_ok(ierr)) return

    call deserialize_complex_flat(arr, filename, ierr)
end subroutine deserialize_complex_nd_C