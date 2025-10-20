!> Module for deserializing character arrays from files
module char_deserialize_mod
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use array_utils, only : ascii_to_string, read_file_header, check_okay_ioerror, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_char_1d, deserialize_char_2d, deserialize_char_3d, deserialize_char_flat, &
          deserialize_char_4d, deserialize_char_5d

contains
  !> Subroutine to deserialize a flat character array from a file
  subroutine deserialize_char_flat(flat, dims, clen, filename, ierr)
    character(len=:), pointer, intent(out) :: flat(:)
      !! Output flat character array
    integer(int32), allocatable, intent(out) :: dims(:)
      !! Output dimensions of the array
    integer(int32), intent(out) :: clen
      !! Maximum length of character string
    character(len=*), intent(in) :: filename
      !! Name of the file to read
    integer(int32), intent(out) :: ierr
      !! Error code

    integer(int32) :: ioerror
      !! Internal I/O error code

    integer(int32) :: unit, magic, type_code, ndim

    call set_ok(ierr)
    call set_ok(ioerror)
    ! open file and read header
    call read_file_header(filename, unit, type_code, ndim, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    ! Allocate output array with stored character length
    allocate(character(len=clen) :: flat(product(dims)))
    
    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) flat
    call check_okay_ioerror(ioerror, ierr, ERR_READ_DATA, unit)
    if(.not. is_ok(ierr)) then
      deallocate(flat)
      return
    end if
    
    close(unit)
  end subroutine deserialize_char_flat

  !> Directly deserialize a 1D character array from a file (array already allocated)
  subroutine deserialize_char_1d(arr, filename, ierr)
    character(len=*), intent(out) :: arr(:)
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

    call check_okay_ndims(ndims, 1, unit, ierr)
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

    call check_okay_ndims(ndims, 2, unit, ierr)
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

    call check_okay_ndims(ndims, 3, unit, ierr)
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

    call check_okay_ndims(ndims, 4, unit, ierr)
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

    call check_okay_ndims(ndims, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    ! Read the entire array as a contiguous block
    read(unit, iostat=ioerror) arr
    close(unit)

    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_char_5d
  
end module char_deserialize_mod

!> Subroutine to deserialize a flat character array from a file and return it as an ASCII array callable by R
!> @note The array is returned flat and needs to be reshaped in R
subroutine deserialize_char_flat_r(ascii_arr, arr_size, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32
  use char_deserialize_mod, only: deserialize_char_flat
  use array_utils, only: ascii_to_string, string_to_ascii_arr
  use tox_errors, only : set_ok, is_ok
  implicit none

  ! Arrays are allocated by R
  integer(int32), intent(in) :: fn_len
  !! Length of the filename array
  integer(int32), intent(in) :: arr_size
  !! Size of the ASCII array
  integer(int32), intent(out) :: ascii_arr(arr_size)
  !! Output array of ASCII characters, preallocated by R
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(int32), intent(out) :: ierr
  !! Error code

  character(len=:), allocatable :: filename
  !! Filename as a string
  character(len=:), pointer :: flat(:)
  !! Flat character array
  integer(int32), allocatable :: dims(:)
  !! Output dimensions of the array
  integer(int32) :: clen
  !! Maximum length of character string
  integer(int32) :: total_array_size
  !! Total size of the ASCII array

  call set_ok(ierr)

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Deserialize flat character array (already trimmed in deserialize_char_flat)
  call deserialize_char_flat(flat, dims, clen, filename, ierr)
  if(.not. is_ok(ierr)) then
    if(associated(flat)) DEALLOCATE(flat)
    return
  end if
  total_array_size = product(dims)

  ! Write data to ASCII array
  call string_to_ascii_arr(flat, total_array_size, ascii_arr, clen)

  deallocate(flat)
end subroutine deserialize_char_flat_r

!> C binding for the subroutine to deserialize a flat character array from a file
!> @note The array is returned flat and needs to be reshaped in C/python
subroutine deserialize_char_flat_C(ascii_arr, clen, total_array_size, &
                                   filename_ascii, fn_len, ierr) bind(C, name="deserialize_char_flat_C")
  use iso_c_binding, only: c_int
  use iso_fortran_env, only: int32
  use char_deserialize_mod, only: deserialize_char_flat
  use array_utils, only : ascii_to_string, string_to_ascii_arr
  use tox_errors, only : is_ok, set_ok
  implicit none

  ! Arguments
  integer(c_int), value       :: clen
  !! Length of each character string
  integer(c_int), value       :: total_array_size
  !! Total size of the ASCII array
  integer(c_int), intent(out) :: ascii_arr(clen*total_array_size)
  !! Output array of ASCII characters, preallocated by C/Python (flat)
  integer(c_int), value       :: fn_len
  !! Length of the filename array
  integer(c_int), intent(in)  :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(c_int), intent(out) :: ierr
  !! error code

  character(len=:), allocatable :: filename
  character(len=:), pointer     :: flat(:)
  integer(c_int), allocatable   :: dims(:)
  integer(int32) :: i, j, actual_clen

  call set_ok(ierr)
  call ascii_to_string(filename_ascii, fn_len, filename)

  ! Deserialize
  call deserialize_char_flat(flat, dims, actual_clen, filename, ierr)
  if (.not. is_ok(ierr)) then
    if (associated(flat)) deallocate(flat)
    return
  end if

  ! Convert to ASCII (null-padding), flatten manually
  call string_to_ascii_arr(flat, total_array_size, ascii_arr, clen)

  deallocate(flat)
end subroutine