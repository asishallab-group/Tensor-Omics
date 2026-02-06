!> Module for serializing complex arrays to binary files.
module f42_serialize_complex
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  use f42_array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_complex_1d, serialize_complex_2d, serialize_complex_3d, &
           serialize_complex_4d, serialize_complex_5d, serialize_complex_nd

  integer(int32), parameter :: ARRAY_TYPE_COMPLEX = 5

contains

  !> Serialize a 1D complex(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_complex_1d(arr, filename, ierr)
    complex(real64), intent(in) :: arr(:)
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

    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, 1, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D complex(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_complex_2d(arr, filename, ierr)
    complex(real64), intent(in) :: arr(:,:)
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

    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, 2, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 3D complex(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_complex_3d(arr, filename, ierr)
    complex(real64), intent(in) :: arr(:,:,:)
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
    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, 3, dims, ierr)
    if (.not. is_ok(ierr)) return 

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 4D complex(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_complex_4d(arr, filename, ierr)
    complex(real64), intent(in) :: arr(:,:,:,:)
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

    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, 4, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 5D complex(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_complex_5d(arr, filename, ierr)
    complex(real64), intent(in) :: arr(:,:,:,:,:)
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
    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, 5, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a flat complex array with specified dimensions and number of dimensions to a binary file.
  !> @note this is called by R
  subroutine serialize_complex_nd(arr, dims, ndim, filename, ierr)
    complex(real64), intent(in) :: arr(:)
    !! Flat complex array to serialize
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

    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, ndim, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine
end module f42_serialize_complex


!> Serialize a flat complex array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_complex_*d to R is not needed.
subroutine serialize_complex_flat_r(arr, array_size, dims, ndim, filename_raw, fn_len, ierr)
  use iso_fortran_env, only: int32, real64
  use f42_serialize_complex, only: serialize_complex_nd
  use tox_conversions, only : c_char_1d_as_string
  use iso_c_binding, only: c_char
  use tox_errors, only : set_ok, is_ok
  implicit none

  integer(int32), intent(in) :: ndim 
  !! Number of dimensions
  integer(int32), intent(in) :: array_size
  !! Size of the flat array
  complex(real64), intent(in) :: arr(array_size)
  !! Flat complex array to serialize
  integer(int32), intent(in) :: dims(ndim)
  !! Dimensions of the array
  integer(int32), intent(in) :: fn_len
  !! Length of the filename array
  character(kind=c_char, len=1), intent(in) :: filename_raw(fn_len)
  !! Array of raw bytes characters representing the filename
  integer(int32), intent(out) :: ierr
  !! Error code

  character(len=:), allocatable :: filename
  integer(int32) :: i, total_len

  call set_ok(ierr)

  call c_char_1d_as_string(filename_raw, filename, ierr)
  if (.not. is_ok(ierr)) return

  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  call serialize_complex_nd(arr(1:total_len), dims(1:ndim), ndim, filename, ierr)
end subroutine serialize_complex_flat_r

!> C binding for the subroutine to serialize a flat complex array to a binary file.
subroutine serialize_complex_nd_C(arr, dims, ndim, filename_raw, fn_len, ierr) bind(C, name="serialize_complex_nd_C")
  use iso_c_binding, only: c_int, c_double_complex, c_char
  use f42_serialize_complex, only: serialize_complex_nd
  use tox_errors, only : set_ok, is_ok
  use iso_fortran_env, only : int32, real64
  use tox_conversions, only : c_char_1d_as_string
  implicit none

  ! input
  integer(c_int), value :: ndim
    !! Number of dimensions
  integer(c_int), intent(in) :: dims(ndim)
    !! Dimensions of the array
  complex(c_double_complex) :: arr(product(dims))  
    !! Pointer to the flat complex array
  integer(c_int), value :: fn_len
    !! Length of the filename array
  character(kind=c_char, len=1), intent(in) :: filename_raw(fn_len)
    !! Array of raw bytes characters representing the filename
  integer(c_int), intent(out) :: ierr
    !! Error code

  ! Local
  character(len=:), allocatable :: filename

  call set_ok(ierr)

  call c_char_1d_as_string(filename_raw, filename, ierr)
  if (.not. is_ok(ierr)) return

  ! save
  call serialize_complex_nd(arr, dims, ndim, filename, ierr)
end subroutine serialize_complex_nd_C