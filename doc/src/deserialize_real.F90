!> Module for deserializing real (double precision) arrays from binary files
module real_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use array_utils, only: ascii_to_string, read_file_header, check_okay_dims, check_okay_ndims
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors
  implicit none

  private
  public :: deserialize_real_1d, deserialize_real_2d, &
           deserialize_real_3d, deserialize_real_4d, deserialize_real_5d

contains

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

    call check_okay_ndims(ndims, 1, unit, ierr)
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

    call check_okay_ndims(ndims, 2, unit, ierr)
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

    call check_okay_ndims(ndims, 3, unit, ierr)
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

    call check_okay_ndims(ndims, 4, unit, ierr)
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

    call check_okay_ndims(ndims, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_real_5d

end module real_deserialize_mod

!> R binding for the subroutine to deserialize a flat real array from a file
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_flat_r(flat_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: real64, int32
  use array_utils, only : ascii_to_string, read_file_header
  use tox_errors, only : set_ok, set_err_once, is_ok, ERR_SIZE_MISMATCH, ERR_READ_DATA
  implicit none

  integer(int32), intent(in) :: fn_len
  !! length of the filename
  integer(int32), intent(in) :: arr_size
  !! size of the array
  real(real64), intent(out) :: flat_arr(arr_size)
  !! array provided by R
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! filename in ascii
  integer(int32), intent(out) :: ierr
  !! error code
  integer(int32) :: ioerror
  !! internal fortran error

  character(len=:), allocatable :: filename
  !! filename
  integer(int32), allocatable :: dims(:)
  !! dimensions
  integer(int32) :: unit, type_code, ndims, clen

  call set_ok(ierr)
  call set_ok(ioerror)
  call ascii_to_string(filename_ascii, fn_len, filename)

  call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
  if(.not. is_ok(ierr)) return

  if (product(dims) /= arr_size) then
    call set_err_once(ierr, ERR_SIZE_MISMATCH)
    close(unit)
    return
  end if

  read(unit, iostat=ioerror) flat_arr
  close(unit)
  if (.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_READ_DATA)
    return
  end if

end subroutine


!> C binding for the subroutine to deserialize a real array from a file
!> @note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_real_C(arr, arr_size, filename_ascii, fn_len, ierr) bind(C, name="deserialize_real_C")
    use iso_c_binding, only : c_int, c_double
    use iso_fortran_env, only: int32, real64
    use array_utils, only: ascii_to_string, read_file_header
    use tox_errors, only : set_ok, set_err_once, is_ok, ERR_SIZE_MISMATCH, ERR_READ_DATA
    implicit none

    ! Inputs / Outputs
    integer(c_int), value         :: arr_size  
    !! size of the output array
    real(c_double), intent(out)   :: arr(arr_size)
    !! output array
    integer(c_int), value         :: fn_len
    !! length of the filename
    integer(c_int), intent(in)    :: filename_ascii(fn_len)
    !! Filename in ascii
    integer(c_int), intent(out)   :: ierr
    !! error code

    integer(int32) :: ioerror
    !! internal fortran error

    ! Locals
    character(len=:), allocatable :: filename
    !! filename
    integer(int32), allocatable   :: dims(:)
    !! dimensions
    integer(int32)                :: unit
    integer(int32)                :: type_code, ndims, clen

    ierr = 0

    ! ASCII to String
    call ascii_to_string(filename_ascii, fn_len, filename)

    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if(.not. is_ok(ierr)) return

    ! Safety check: ensure provided buffer matches size in file
    if (product(dims) /= arr_size) then
        call set_err_once(ierr, ERR_SIZE_MISMATCH)
        close(unit)
        return
    end if

    ! Read directly into provided buffer
    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
        call set_err_once(ierr, ERR_READ_DATA)
        return
    end if
end subroutine deserialize_real_C