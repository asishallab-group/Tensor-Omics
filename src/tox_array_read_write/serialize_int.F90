!> Module for serializing integer arrays to binary files.
module serialize_int
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  use array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_int_1d, serialize_int_2d, serialize_int_3d, &
           serialize_int_4d, serialize_int_5d, serialize_int_nd

  integer(int32), parameter :: ARRAY_TYPE_INT = 1

contains

  !> Serialize a 1D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_1d(arr, filename, ierr)
    integer(int32), intent(in) :: arr(:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code
    integer(int32) :: unit
    integer(int32) :: ioerror
    integer(int32) :: dims(1)
    dims = shape(arr)
    
    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_INT, 1, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_2d(arr, filename, ierr)
    integer(int32), intent(in) :: arr(:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32) :: ierr
    !! error code
    integer(int32) :: unit
    integer(int32) :: ioerror
    integer(int32) :: dims(2)
    dims = shape(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_INT, 2, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 3D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_3d(arr, filename, ierr)
    integer(int32), intent(in) :: arr(:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32) :: ierr
    !! error code
    integer(int32) :: unit
    integer(int32) :: ioerror
    integer(int32) :: dims(3)
    dims = shape(arr)

    call set_ok(ierr)
    call set_ok(ioerror)
    call write_file_header(filename, unit, ARRAY_TYPE_INT, 3, dims, ierr)
    if (.not. is_ok(ierr)) return 

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 4D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_4d(arr, filename, ierr)
    integer(int32), intent(in) :: arr(:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: dims(4)
    integer(int32) :: unit
    dims = shape(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_INT, 4, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 5D integer(int32) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_int_5d(arr, filename, ierr)
    integer(int32), intent(in) :: arr(:,:,:,:,:)
    !! array to save
    character(len=*), intent(in) :: filename
    !! output filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: unit
    integer(int32) :: ioerror
    integer(int32) :: dims(5)
    
    dims = shape(arr)

    call set_ok(ierr)
    call set_ok(ioerror)
    call write_file_header(filename, unit, ARRAY_TYPE_INT, 5, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
  !> @note this is called by R
  subroutine serialize_int_nd(arr, dims, ndim, filename, ierr)
    integer(int32), intent(in) :: arr(:)
    !! Flat integer array to serialize
    integer(int32), intent(in) :: dims(:)
    !! Dimensions of the array
    integer(int32), intent(in) :: ndim
    !! Number of dimensions
    character(len=*), intent(in) :: filename
    !! filename
    integer(int32) :: unit
    integer(int32), INTENT(OUT) :: ierr
    !! error code
    integer(int32) :: ioerror

    call set_ok(ierr)
    call set_ok(ioerror)

    if (size(dims) /= ndim) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
    end if

    call write_file_header(filename, unit, ARRAY_TYPE_INT, ndim, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine
end module serialize_int


!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
subroutine serialize_int_flat_r(arr, array_size, dims, ndim, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32
  use array_utils
  use serialize_int, only: serialize_int_nd
  use tox_errors, only : set_ok
  implicit none

  integer(int32), intent(in) :: ndim 
  !! Number of dimensions
  integer(int32), intent(in) :: array_size
  !! Size of the flat array
  integer(int32), intent(in) :: arr(array_size)
  !! Flat integer array to serialize
  integer(int32), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(int32), intent(in) :: fn_len
  !! Length of the filename array
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !! Array of ASCII characters representing the filename
  integer(int32), intent(out) :: ierr
  !! Error code

  character(len=:), allocatable :: filename
  integer(int32) :: i, total_len

  call set_ok(ierr)

  call ascii_to_string(filename_ascii, fn_len, filename)

  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  call serialize_int_nd(arr(1:total_len), dims(1:ndim), ndim, filename, ierr)
end subroutine

!> C binding for the subroutine to serialize a flat integer array to a binary file.
subroutine serialize_int_nd_C(arr, dims, ndim, filename_ascii, fn_len, ierr) bind(C, name="serialize_int_nd_C")
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer
  use array_utils, only: ascii_to_string
  use serialize_int, only: serialize_int_nd
  use tox_errors, only : set_ok
  use iso_fortran_env, only : int32
  implicit none

  ! input
  type(c_ptr), value :: arr
    !! Pointer to the flat integer array
  integer(c_int), value :: ndim
    !! Number of dimensions
  integer(c_int), intent(in) :: dims(ndim)
    !! Dimensions of the array
  integer(c_int), value :: fn_len
    !! Length of the filename array
  integer(c_int), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(c_int), intent(out) :: ierr
    !! Error code

  ! Local
  character(len=:), allocatable :: filename
  integer(c_int), pointer :: arr_f(:)

  call set_ok(ierr)

  call ascii_to_string(filename_ascii, fn_len, filename)

  call c_f_pointer(arr, arr_f, [product(dims(1:ndim))])

  ! save
  call serialize_int_nd(arr_f, dims, ndim, filename, ierr)
end subroutine
