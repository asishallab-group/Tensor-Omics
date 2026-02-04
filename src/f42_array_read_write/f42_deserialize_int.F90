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

!> R interface for deserializing an integer array from a file.
!> Deserializes an array of any dimension into a flat buffer.
!> @note The output array is handled and preallocated by R
subroutine deserialize_int_r(flat_arr, arr_size, filename_raw, fn_len, ierr)
  use iso_fortran_env, only: int32
  use iso_c_binding, only : c_char
  use f42_deserialize_int, only : deserialize_int_flat
  use tox_conversions, only : c_char_1d_as_string
  use tox_errors, only : set_ok, is_ok
  implicit none

  integer(int32), intent(in)  :: arr_size
  !! size of the array
  integer(int32), intent(out) :: flat_arr(arr_size)
  !! array passed by R
  integer(int32), intent(in)  :: fn_len
  !! length of the filename
  character(kind=c_char, len=1), intent(in)  :: filename_raw(fn_len)
  !! filename to read from
  integer(int32), intent(out) :: ierr
  !! error code
  
  character(len=:), allocatable :: filename
  !! filename in characters

  call set_ok(ierr)

  call c_char_1d_as_string(filename_raw, filename, ierr)
  if( .not. is_ok(ierr)) return

  call deserialize_int_flat(flat_arr, filename, ierr)
end subroutine

!> C binding for the subroutine to deserialize an integer array from a file
!> Deserializes an array of any dimension into a flat buffer.
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_int_nd_C(arr, arr_size, filename_raw, fn_len, ierr) bind(C, name="deserialize_int_nd_C")
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32
    use f42_deserialize_int, only : deserialize_int_flat
    use tox_errors, only : set_ok, is_ok
    use tox_conversions, only : c_char_1d_as_string
    implicit none

    ! Inputs / Outputs
    integer(c_int), value         :: arr_size           ! Buffer length
    !! Size of the array
    integer(c_int), intent(out)   :: arr(arr_size)      ! Preallocated buffer from C/Python
    !! preallocated array
    integer(c_int), value         :: fn_len
    !! length of the filename
    character(kind=c_char, len=1), intent(in)    :: filename_raw(fn_len)
    !! Filename in ascii
    integer(c_int), intent(out)   :: ierr
    !! Error code

    ! Locals
    character(len=:), allocatable :: filename

    call set_ok(ierr)

    ! raw → String
    call c_char_1d_as_string(filename_raw, filename, ierr)
    if (.not. is_ok(ierr)) return

    call deserialize_int_flat(arr, filename, ierr)
end subroutine deserialize_int_nd_C