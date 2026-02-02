!> Module providing serialization and deserialization routines for character arrays
!! of up to 5 dimensions, arrays are serialized to a custom binary format with a magic number and type/dimension metadata.

module f42_serialize_char
  use safeguard
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use f42_array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
           serialize_char_4d, serialize_char_5d, serialize_char_nd

  integer(int32), parameter :: ARRAY_TYPE_CHAR = 3

contains

  !> Serialize a 1D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_1d(arr, filename, ierr)
    character(len=*), contiguous, intent(in) :: arr(:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr 
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(1)
    integer(int32) :: ioerror
    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 1_int32, dims, ierr, clen)
    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if(.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_2d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(2)
    integer(int32) :: ioerror
    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 2_int32, dims, ierr, clen)

    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if(.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 3D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_3d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit, clen
    integer(int32) :: dims(3)
    dims = shape(arr)
    clen = len(arr)
    
    call set_ok(ierr)
    call set_ok(ioerror)
    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 3_int32, dims, ierr, clen)
    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 4D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_4d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(4)
    integer(int32) :: ioerror

    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 4_int32, dims, ierr, clen)
    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 5D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_5d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:,:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit, clen
    integer(int32) :: dims(5)
    integer(int32) :: ioerror
    dims = shape(arr)
    clen = len(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 5_int32, dims, ierr, clen)

    if (.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) arr
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a character array of arbitrary dimensions to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  !! @note This routine is only called by R and serializes only flat character arrays to the memory
  subroutine serialize_char_nd(flat, dims, ndim, clen, filename, ierr)
    implicit none
    character(len=*), intent(in) :: flat(:)
    !! flat array to save
    integer(int32), intent(in) :: dims(:)
    !! dimensions of the array
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32), intent(in) :: ndim
    !! number of dimensions
    integer(int32), intent(in) :: clen
    !! Length of each string
    integer(int32) :: ioerror
    integer(int32) :: unit

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, ndim, dims, ierr, clen)
    if(.not. is_ok(ierr)) return

    ! Write the entire array as a contiguous block
    write(unit, iostat=ioerror) flat
    
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if

    close(unit)
  end subroutine serialize_char_nd

end module f42_serialize_char

!> serializes a flat character array to a binary file.
subroutine serialize_char_flat_r(raw_arr, array_size, dims, ndim, clen, filename_raw, fn_len, ierr)
  use iso_fortran_env, only: int32
  use iso_c_binding, only: c_char
  use f42_serialize_char, only: serialize_char_nd
  use tox_errors, only: set_err_once, set_ok, is_ok, ERR_ALLOC_FAIL
  use tox_conversions, only: c_char_1d_as_string, c_char_2d_as_string
  implicit none

  integer(int32), intent(in) :: ndim
  !! Number of dimensions
  integer(int32), intent(in) :: array_size
  !! size of the input array
  integer(int32), intent(in) :: clen
  !! character length
  character(kind=c_char, len=1), intent(in) :: raw_arr(clen, array_size)
  !! Flat character array in ASCII format
  integer(int32), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(int32), intent(in) :: fn_len
  !! length of the filename
  character(kind=c_char, len=1), intent(in) :: filename_raw(fn_len)
  !! filename in ascii
  integer(int32), intent(out) :: ierr
  !! error code
  integer(int32) :: ioerror

  character(len=:), allocatable :: filename
  character(len=:), allocatable :: flat(:)
  integer(int32) :: i, j, total

  call set_ok(ioerror)
  call set_ok(ierr)

  total = product(dims)

  if(.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_ALLOC_FAIL)
    RETURN
  end if  

  ! ASCII to character conversion
  call c_char_2d_as_string(raw_arr, flat, ierr)
  if (.not. is_ok(ierr)) return

  call c_char_1d_as_string(filename_raw, filename, ierr)
  if (.not. is_ok(ierr)) return

  call serialize_char_nd(flat, dims, ndim, clen, filename, ierr)

end subroutine serialize_char_flat_r

subroutine serialize_char_nd_C(raw_chars, dims, ndim, clen, filename_c, fn_len, ierr) bind(C, name="serialize_char_nd_C")
  use iso_c_binding, only: c_char, c_int, c_null_char
  use iso_fortran_env, only: int32
  use f42_serialize_char, only: serialize_char_nd
  use tox_conversions, only: c_char_1d_as_string, c_char_as_char
  use tox_errors
  implicit none

  integer(c_int), intent(in), value :: ndim
    !! Number of dimensions
  integer(c_int), intent(in) :: dims(ndim)
    !! Dimensions of the array
  integer(c_int), value :: clen
    !! Length of each string
  character(kind=c_char, len=1), intent(in) :: raw_chars(clen, product(dims))
    !! Raw chars array
  integer(c_int), value :: fn_len
    !! length of the filename array
  character(kind=c_char, len=1), intent(in) :: filename_c(fn_len)
    !! filename as c_char array
  integer(c_int), intent(out) :: ierr
    !! Error code

  character(len=:), allocatable :: filename
  character(len=clen), allocatable :: flat(:)
  integer(int32) :: i, j, total, ioerror
  character(len=1) :: fortran_char

  call set_ok(ierr)

  total = product(dims)
  allocate(flat(total), stat=ioerror)

  if(.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_ALLOC_FAIL)
    return
  end if

  do i = 1, total
    flat(i) = ""
    do j = 1, clen
      if (raw_chars(j, i) /= c_null_char) then
        call c_char_as_char(raw_chars(j, i), fortran_char)
        flat(i)(j:j) = fortran_char
      else
        exit
      end if
    end do
  end do

  ! Convert filename from c_char array to Fortran string
  call c_char_1d_as_string(filename_c, filename, ierr)
  if (is_err(ierr)) return

  call serialize_char_nd(flat, dims, ndim, clen, filename, ierr)
  
end subroutine serialize_char_nd_C