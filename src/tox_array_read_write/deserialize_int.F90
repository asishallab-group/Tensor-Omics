!> Module for deserializing integer arrays from files
module int_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use array_utils, only: ascii_to_string, read_file_header, check_okay_dims, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_int_1d, deserialize_int_2d, &
           deserialize_int_3d, deserialize_int_4d, deserialize_int_5d

contains
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

    call check_okay_ndims(ndims, 1, unit, ierr)
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

    call check_okay_ndims(ndims, 2, unit, ierr)
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

    call check_okay_ndims(ndims, 3, unit, ierr)
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

    call check_okay_ndims(ndims, 4, unit, ierr)
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

    call check_okay_ndims(ndims, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_int_5d

end module int_deserialize_mod

!> R interface for deserializing an integer array from a file
!> @note The output array is handled and preallocated by R
subroutine deserialize_int_r(flat_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32
  use array_utils, only : ascii_to_string, read_file_header
  use tox_errors, only : set_err_once, set_ok, is_ok, ERR_SIZE_MISMATCH, ERR_READ_DATA
  implicit none

  integer(int32), intent(in)  :: arr_size
  !! size of the array
  integer(int32), intent(out) :: flat_arr(arr_size)
  !! array passed by R
  integer(int32), intent(in)  :: fn_len
  !! length of the filename
  integer(int32), intent(in)  :: filename_ascii(fn_len)
  !! filename to read from
  integer(int32), intent(out) :: ierr
  !! error code
  
  integer(int32) :: ioerror
  character(len=:), allocatable :: filename
  !! filename in characters
  integer(int32), allocatable   :: dims(:)
  !! dimensions of the array
  integer(int32)                :: unit, type_code, ndims, clen

  call set_ok(ierr)
  call set_ok(ioerror)

  call ascii_to_string(filename_ascii, fn_len, filename)

  call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
  if (.not. is_ok(ierr)) return

  if (product(dims) /= arr_size) then
    call set_err_once(ierr, ERR_SIZE_MISMATCH)
    return
  end if

  ! Read directly into R buffer
  read(unit, iostat=ioerror) flat_arr
  close(unit)
  if (.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_READ_DATA)
    RETURN
  end if
end subroutine


!> C binding for the subroutine to deserialize an integer array from a file
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_int_C(arr, arr_size, filename_ascii, fn_len, ierr) bind(C, name="deserialize_int_C")
    use iso_c_binding, only: c_int
    use iso_fortran_env, only: int32
    use array_utils, only: ascii_to_string, read_file_header
    use tox_errors, only : set_err_once, set_ok, is_ok, ERR_SIZE_MISMATCH, ERR_READ_DATA
    implicit none

    ! Inputs / Outputs
    integer(c_int), value         :: arr_size           ! Buffer length
    !! Size of the array
    integer(c_int), intent(out)   :: arr(arr_size)      ! Preallocated buffer from C/Python
    !! preallocated array
    integer(c_int), value         :: fn_len
    !! length of the filename
    integer(c_int), intent(in)    :: filename_ascii(fn_len)
    !! Filename in ascii
    integer(c_int), intent(out)   :: ierr
    !! Error code

    ! Locals
    character(len=:), allocatable :: filename
    integer(int32), allocatable   :: dims(:)
    integer(int32)                :: unit
    integer(int32)                :: ioerror
    integer(int32)                :: type_code, ndims, clen

    call set_ok(ierr)
    call set_ok(ioerror)

    ! ASCII → String
    call ascii_to_string(filename_ascii, fn_len, filename)

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    ! Safety check: ensure provided buffer matches size in file
    if (product(dims) /= arr_size) then
        call set_err_once(ierr, ERR_SIZE_MISMATCH)
        close(unit)
        return
    end if

    ! Read directly into C/Python-provided buffer → ZERO COPY
    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_READ_DATA)
        return
    end if
end subroutine deserialize_int_C