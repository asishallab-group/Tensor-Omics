!> Module providing serialization and deserialization routines for character arrays
!! of up to 5 dimensions, arrays are serialized to a custom binary format with a magic number and type/dimension metadata.

module serialize_char
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_char_1d, serialize_char_2d, serialize_char_3d, &
           serialize_char_4d, serialize_char_5d, serialize_char_nd

  integer(int32), parameter :: ARRAY_TYPE_CHAR = 3

contains

  !> Serialize a 1D character array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, character length, and the array data.
  subroutine serialize_char_1d(arr, filename, ierr)
    character(len=*), intent(in) :: arr(:)
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

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 1, dims, ierr, clen)
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
    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 2, dims, ierr, clen)

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
    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 3, dims, ierr, clen)
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

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 4, dims, ierr, clen)
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

    call write_file_header(filename, unit, ARRAY_TYPE_CHAR, 5, dims, ierr, clen)

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

end module serialize_char

!> serializes a flat character array to a binary file.
subroutine serialize_char_flat_r(ascii_arr, array_size, dims, ndim, clen, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32
  use serialize_char, only: serialize_char_nd
  use array_utils, only: ascii_to_string
  use tox_errors
  implicit none

  integer(int32), intent(in) :: ndim
  !! Number of dimensions
  integer(int32), intent(in) :: array_size
  !! size of the input array
  integer(int32), intent(in) :: clen
  !! character length
  integer(int32), intent(in) :: ascii_arr(clen, array_size)
  !! Flat character array in ASCII format
  integer(int32), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(int32), intent(in) :: fn_len
  !! length of the filename
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! filename in ascii
  integer(int32), intent(out) :: ierr
  !! error code
  integer(int32) :: ioerror

  character(len=:), allocatable :: filename
  character(len=clen), allocatable :: flat(:)
  integer(int32) :: i, j, total

  call set_ok(ioerror)
  call set_ok(ierr)

  total = product(dims)
  allocate(flat(total), stat=ioerror)

  if(.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_ALLOC_FAIL)
    RETURN
  end if  

  ! ASCII to character conversion
  do i = 1, total
    flat(i) = ""
    do j = 1, clen
      if (ascii_arr(j, i) > 0) then
        flat(i)(j:j) = char(ascii_arr(j, i))
      else
        exit
      end if
    end do
  end do

  call ascii_to_string(filename_ascii, fn_len, filename)

  call serialize_char_nd(flat, dims, ndim, clen, filename, ierr)

end subroutine serialize_char_flat_r

!> C binding for the subroutine to serialize a flat character array to a binary file.
subroutine serialize_char_flat_C(ascii_ptr, dims, ndim, clen, filename_ascii, fn_len, ierr) bind(C, name="serialize_char_flat_C")
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use iso_fortran_env, only: int32
  use serialize_char, only: serialize_char_nd
  use array_utils, only: ascii_to_string
  use tox_errors
  implicit none

  type(c_ptr), value :: ascii_ptr
  !! pointer to ascii array
  integer(c_int), value :: ndim
  !! Number of dimensions  
  integer(c_int), intent(in) :: dims(ndim)
  !! Dimensions of the array
  
  integer(c_int), value :: clen
  !! Character length
  integer(c_int), value :: fn_len
  !! Length of the filename array
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(c_int), intent(out) :: ierr
  !! error code

  integer(c_int), pointer :: ascii_arr(:)
  character(len=:), allocatable :: filename
  character(len=clen), allocatable :: flat(:)
  integer(int32) :: i, j, total, ioerror

  call set_ok(ierr)
  call set_ok(ioerror)

  total = product(dims)
  call c_f_pointer(ascii_ptr, ascii_arr, [clen * total])
  allocate(flat(total), stat=ioerror)

  if(.not. is_ok(ioerror)) then
    call set_err_once(ierr, ERR_ALLOC_FAIL)
    RETURN
  end if

  ! ASCII to Fortran character(len=clen)
  do i = 1, total
    flat(i) = ""
    do j = 1, clen
      if (ascii_arr((i - 1) * clen + j) > 0) then
        flat(i)(j:j) = char(ascii_arr((i - 1) * clen + j))
      else
        exit
      end if
    end do
  end do

  call ascii_to_string(filename_ascii, fn_len, filename)

  call serialize_char_nd(flat, dims, ndim, clen, filename ,ierr)
end subroutine serialize_char_flat_C