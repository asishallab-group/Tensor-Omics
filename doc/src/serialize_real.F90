module serialize_real
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only: c_loc
  use array_utils, only: write_file_header
  use tox_errors
  implicit none

  public:: serialize_real_1d, serialize_real_2d, serialize_real_3d, &
           serialize_real_4d, serialize_real_5d, serialize_real_nd

  integer(int32), parameter :: ARRAY_TYPE_REAL = 2
contains

  !> Serialize a 1D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_1d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! output filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: unit
    integer(int32) :: dims(1)
    integer(int32) :: ioerror
    dims = shape(arr)

    call set_ok(ierr)
    call set_ok(ioerror)

    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 1, dims, ierr)
    if (.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 2D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_2d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: unit
    integer(int32) :: dims(2)
    integer(int32) :: ioerror

    call set_ok(ierr)
    call set_ok(ioerror)

    dims = shape(arr)
    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 2, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 3D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_3d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit
    integer(int32) :: dims(3)

    call set_ok(ierr)
    call set_ok(ioerror)

    dims = shape(arr)
    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 3, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 4D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_4d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), INTENT(OUT) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit
    integer(int32) :: dims(4)

    call set_ok(ierr)
    call set_ok(ioerror)
    dims = shape(arr)

    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 4, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> Serialize a 5D real(real64) array to a binary file.
  !! The file will contain a magic number, type code, dimension, shape, and the array data.
  subroutine serialize_real_5d(arr, filename, ierr)
    real(real64), intent(in) :: arr(:,:,:,:,:)
      !! array to save
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), intent(out) :: ierr
    !! error code

    integer(int32) :: ioerror
    integer(int32) :: unit
    integer(int32) :: dims(5)
    
    call set_ok(ierr)
    call set_ok(ioerror)

    dims = shape(arr)
    call write_file_header(filename, unit, ARRAY_TYPE_REAL, 5, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine

  !> @brief Writes serialized real array from R to file with metdata.
  subroutine serialize_real_nd(arr, dims, ndim, filename, ierr)
    real(real64), intent(in) :: arr(:)
      !! array to save
    integer(int32), intent(in) :: dims(:)
      !! Dimensions of the array
    integer(int32), intent(in) :: ndim
      !! Number of dimensions
    character(len=*), intent(in) :: filename
      !! filename
    integer(int32), INTENT(OUT) :: ierr
      !! error code

    integer(int32) :: ioerror, unit

    call set_ok(ierr)
    call set_ok(ioerror)

    if (size(dims) /= ndim) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if

    call write_file_header(filename, unit, ARRAY_TYPE_REAL, ndim, dims, ierr)
    if(.not. is_ok(ierr)) return

    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_WRITE_DATA)
    end if
    close(unit)
  end subroutine
 
end module serialize_real

!> Serialize a flat integer array with specified dimensions and number of dimensions to a binary file.
!! R can not pass a multidimensional array directly, so we use a flat array and dimensions. Therefore, exposing serialize_int_*d to R is not needed.
subroutine serialize_real_flat_r(arr, array_size, dims, ndim, filename_ascii, fn_len, ierr)
  use iso_fortran_env, only: int32, real64
  use array_utils, only: ascii_to_string
  use serialize_real, only: serialize_real_nd
  use tox_errors, only : set_ok
  implicit none

  integer(int32), intent(in) :: ndim
    !! Number of dimensions
  integer(int32), intent(in) :: array_size
    !! Size of the flat array
  real(real64), intent(in) :: arr(array_size) 
    !! Flat real array to serialize
  integer(int32), intent(in) :: dims(ndim)
    !! Dimensions of the array
  integer(int32), intent(in) :: fn_len
    !! length of the filename array
  integer(int32), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  character(len=:), allocatable :: filename
  integer(int32), intent(out) :: ierr
    !! Error code

  integer(int32) :: i, total_len

  call set_ok(ierr)

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! calculate total size
  total_len = 1
  do i = 1, ndim
    total_len = total_len * dims(i)
  end do

  call serialize_real_nd(arr(1:total_len), dims(1:ndim), ndim, filename, ierr)
end subroutine

subroutine serialize_real_nd_C(arr, dims, ndim, filename_ascii, fn_len, ierr) bind(C, name="serialize_real_nd_C")
  use iso_c_binding, only: c_ptr, c_int, c_f_pointer, c_double
  use array_utils, only: ascii_to_string
  use serialize_real, only: serialize_real_nd
  use tox_errors, only : set_ok
  use iso_fortran_env, only : int32
  implicit none

  ! Input parameters
  type(c_ptr), value :: arr
    !! Pointer to the flat real array
  integer(c_int), value :: ndim
    !! Number of dimensions
  integer(c_int), intent(in) :: dims(ndim)
    !! Dimensions of the array
  integer(c_int), value :: fn_len
    !! Length of the filename array
  integer(c_int), intent(in) :: filename_ascii(fn_len)
    !! Array of ASCII characters representing the filename
  integer(c_int), intent(out) :: ierr

  ! Local
  character(len=:), allocatable :: filename
  real(c_double), pointer :: arr_f(:)
  integer(int32) :: i

  call ascii_to_string(filename_ascii, fn_len, filename)

  ! 1D-Array from C pointer
  call c_f_pointer(arr, arr_f, [product(dims(1:ndim))])

  ! save
  call serialize_real_nd(arr_f, dims, ndim, filename, ierr)
end subroutine