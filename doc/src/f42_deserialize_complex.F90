!> Module for deserializing complex arrays from files
module f42_deserialize_complex
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use iso_c_binding, only : c_loc, c_f_pointer
  use f42_array_utils, only: read_file_header, check_okay_ndims
  use tox_errors
  implicit none

  private
  public :: deserialize_complex_1d, deserialize_complex_2d, &
           deserialize_complex_3d, deserialize_complex_4d, deserialize_complex_5d, deserialize_complex_flat

contains

  !> Deserializes any array inot a flat array
  subroutine deserialize_complex_flat(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_flat
  !> Deserialize a flat complex array from a file
  !> Directly deserialize a 1D complex array from a file
  subroutine deserialize_complex_1d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 1, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_1d

  !> Directly deserialize a 2D complex array from a file
  subroutine deserialize_complex_2d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 2, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_2d

  !> Directly deserialize a 3D complex array from a file
  subroutine deserialize_complex_3d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 3, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_3d

  !> Directly deserialize a 4D complex array from a file
  subroutine deserialize_complex_4d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:,:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 4, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_4d

  !> Directly deserialize a 5D complex array from a file
  subroutine deserialize_complex_5d(arr, filename, ierr)
    complex(real64), intent(out) :: arr(:,:,:,:,:)
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

    call validate_type_code(type_code, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    call check_okay_ndims(ndims, 5, unit, ierr)
    if(.not. is_ok(ierr)) return

    read(unit, iostat=ioerror) arr
    close(unit)
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ERR_READ_DATA)
      return
    end if
  end subroutine deserialize_complex_5d

end module f42_deserialize_complex

!> R interface for deserializing a complex array from a file
!> Deserializes an array of any dimension into a flat buffer.
!> @note The output array is handled and preallocated by R
subroutine deserialize_complex_r(flat_arr, arr_size, filename_raw, fn_len, ierr)
  use iso_fortran_env, only: int32, real64
  use f42_deserialize_complex, only : deserialize_complex_flat
  use tox_errors, only : set_ok, is_ok
  use tox_conversions, only : c_char_1d_as_string
  use iso_c_binding, only : c_char
  implicit none

  integer(int32), intent(in)  :: arr_size
  !! size of the array
  complex(real64), intent(out) :: flat_arr(arr_size)
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
  if (.not. is_ok(ierr)) return

  call deserialize_complex_flat(flat_arr, filename, ierr)
end subroutine


!> C binding for the subroutine to deserialize a complex array from a file.
!> The array is read into a flat buffer and reshaped in the calling language.
!>@note It is assumed that the array is already allocated and passed together with its size
subroutine deserialize_complex_nd_C(arr, arr_size, filename_raw, fn_len, ierr) bind(C, name="deserialize_complex_nd_C")
    use iso_c_binding, only: c_int, c_char
    use iso_fortran_env, only: int32, real64
    use f42_deserialize_complex, only : deserialize_complex_flat
    use tox_errors, only : set_ok, is_ok
    use tox_conversions, only : c_char_1d_as_string
    implicit none

    ! Inputs / Outputs
    integer(c_int), value         :: arr_size           ! Buffer length
    !! Size of the array
    complex(real64), intent(out)  :: arr(arr_size)      ! Preallocated buffer from C/Python
    !! preallocated array
    integer(c_int), value         :: fn_len
    !! length of the filename
    character(kind=c_char, len=1), intent(in)    :: filename_raw(fn_len)
    !! Filename in raw bytes
    integer(c_int), intent(out)   :: ierr
    !! Error code

    ! Locals
    character(len=:), allocatable :: filename

    call set_ok(ierr)

    ! raw to String
    call c_char_1d_as_string(filename_raw, filename, ierr)
    if (.not. is_ok(ierr)) return

    call deserialize_complex_flat(arr, filename, ierr)
end subroutine deserialize_complex_nd_C