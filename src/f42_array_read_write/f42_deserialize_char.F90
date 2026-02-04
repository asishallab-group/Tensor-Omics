!> Module for deserializing character arrays from files
module f42_deserialize_char
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use f42_array_utils, only : read_file_header, check_okay_ioerror, check_okay_ndims
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

    integer(int32) :: unit, type_code, ndim

    call set_ok(ierr)
    call set_ok(ioerror)
    ! open file and read header
    call read_file_header(filename, unit, type_code, ndim, dims, clen, ierr)
    if (.not. is_ok(ierr)) return

    call validate_type_code(type_code, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

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

!> Subroutine to deserialize a flat character array from a file and return it as a raw array callable by R
!> Deserializes the array into a flat buffer.
!> @note The array is returned flat and needs to be reshaped in R
subroutine deserialize_char_flat_r(raw_arr, arr_size, filename_raw, fn_len, clen_in, ierr)
  use iso_fortran_env, only: int32
  use f42_deserialize_char, only: deserialize_char_flat
  use tox_errors, only : set_ok, is_ok
  use tox_conversions, only : string_as_c_char_1d_r, c_char_1d_as_string
  use iso_c_binding, only : c_char
  implicit none

  ! Arrays are allocated by R
  integer(int32), intent(in) :: fn_len
  !! Length of the filename array
  integer(int32), intent(in) :: arr_size
  !! Size of the ASCII array
  integer(int32), intent(in) :: clen_in
  !! Maximum length of character string
  character(kind=c_char, len=1), intent(out) :: raw_arr(clen_in, arr_size)
  !! Output array of ASCII characters, preallocated by R
  character(kind=c_char, len=1), intent(in) :: filename_raw(fn_len)
  !! Array of ASCII characters representing the filename
  integer(int32), intent(out) :: ierr
  !! Error code

  character(len=:), allocatable :: filename
  !! Filename as a string
  character(len=:), pointer :: flat(:)
  !! Flat character array
  integer(int32), allocatable :: dims(:)
  !! Output dimensions of the array
  integer(int32) :: clen, i
  integer(int32) :: total_array_size
  !! Total size of the ASCII array

  call set_ok(ierr)

  call c_char_1d_as_string(filename_raw, filename, ierr)
  if(.not. is_ok(ierr)) return

  ! Deserialize flat character array (already trimmed in deserialize_char_flat)
  call deserialize_char_flat(flat, dims, clen, filename, ierr)
  if(.not. is_ok(ierr)) then
    if(associated(flat)) DEALLOCATE(flat)
    return
  end if
  total_array_size = product(dims)

  ! Write data to raw array
  do i = 1, size(flat)
    call string_as_c_char_1d_r(flat(i), raw_arr(:, i))
  end do

  deallocate(flat)
end subroutine deserialize_char_flat_r

!> C binding for the subroutine to deserialize a flat character array from a file.
!> Deserializes the array into a flat buffer.
!> @note The array is returned flat and needs to be reshaped in C/python
subroutine deserialize_char_nd_C(raw_chars, clen, total_array_size, &
                                   filename_c, fn_len, ierr) bind(C, name="deserialize_char_nd_C")
  use iso_c_binding, only: c_char, c_int, c_null_char
  use iso_fortran_env, only: int32
  use f42_deserialize_char, only: deserialize_char_flat
  use tox_conversions, only: c_char_1d_as_string, string_as_c_char_1d
  use tox_errors, only : is_ok, set_ok, set_err_once, ERR_ALLOC_FAIL
  implicit none

  ! Arguments
  integer(c_int), value       :: clen
  !! Length of each character string
  integer(c_int), value       :: total_array_size
  !! Total size of the array
  character(kind=c_char, len=1), intent(out) :: raw_chars(clen, total_array_size)
  !! Output array of c_chars, preallocated by C/Python (2D: clen x total_array_size)
  integer(c_int), value       :: fn_len
  !! Length of the filename array
  character(kind=c_char, len=1), intent(in)  :: filename_c(fn_len)
  !! c_char array representing the filename
  integer(c_int), intent(out) :: ierr
  !! error code

  character(len=:), allocatable :: filename
  character(len=:), pointer     :: flat(:)
  integer(c_int), allocatable   :: dims(:)
  integer(int32) :: i, j, actual_clen

  call set_ok(ierr)
  
  ! Convert filename from c_char array to Fortran string
  call c_char_1d_as_string(filename_c, filename, ierr)
  if (.not. is_ok(ierr)) return

  ! Deserialize
  call deserialize_char_flat(flat, dims, actual_clen, filename, ierr)
  if (.not. is_ok(ierr)) then
    if (associated(flat)) deallocate(flat)
    return
  end if

  ! Check if allocated flat array matches expected size
  if (size(flat) /= total_array_size) then
    call set_err_once(ierr, ERR_ALLOC_FAIL)
    if (associated(flat)) deallocate(flat)
    return
  end if

  ! Convert Fortran strings to c_char array (2D: clen x total_array_size)
  do i = 1, total_array_size
    do j = 1, clen
      if (j <= len(flat(i)) .and. j <= actual_clen) then
        raw_chars(j, i) = flat(i)(j:j)
      else
        raw_chars(j, i) = c_null_char
      end if
    end do
  end do

  deallocate(flat)
end subroutine deserialize_char_nd_C