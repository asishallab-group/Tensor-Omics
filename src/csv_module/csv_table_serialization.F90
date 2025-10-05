!> Module for serializing CSV tables into binary files
module tox_csv_table_serialization
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_DIM_MISMATCH, set_ok, set_err_once, is_err
  use serialize_int, only: serialize_int_2d
  use serialize_real, only: serialize_real_2d
  use serialize_char, only: serialize_char_1d, serialize_char_2d
  use array_utils, only: write_file_header

contains

!// TODO Error handling for invalid or empty inputs

  !> Serialization routines using serial_array_module
  !> Serialize the type-banded arrays from a CSV table to binary files
  subroutine serialize_table(filename_prefix, int_cols, real_cols, char_cols, &
                           logical_cols, complex_cols, header, metadata, ierr)
    implicit none
    
    !| Base filename prefix (extensions will be added automatically)
    character(len=*), intent(in) :: filename_prefix
    !| 2D int array, integer columns to serialize
    integer(int32), intent(in) :: int_cols(:,:)
    !| 2D real array, real columns to serialize  
    real(real64), intent(in) :: real_cols(:,:)
    !| 2D char array, character columns to serialize
    character(len=*), intent(in) :: char_cols(:,:)
    !| 2D logical array, logical columns to serialize
    logical, intent(in) :: logical_cols(:,:)
    !| 2D complex array, complex columns to serialize
    complex(real64), intent(in) :: complex_cols(:,:)
    !| 1D char array, column names/headers
    character(len=*), intent(in) :: header(:)
    !| 2D int array, metadata (row 1: type, row 2: index in type array)
    integer(int32), intent(in) :: metadata(:,:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    character(len=1024) :: filename
    integer(int32) :: local_ierr
    
    call set_ok(ierr)
    call set_ok(local_ierr)
    
    ! Serialize integer columns if present
    if (size(int_cols, 1) > 0 .and. size(int_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_int_cols.dat'
      call serialize_int_2d(int_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize real columns if present
    if (size(real_cols, 1) > 0 .and. size(real_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_real_cols.dat'
      call serialize_real_2d(real_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize character columns if present
    if (size(char_cols, 1) > 0 .and. size(char_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_char_cols.dat'
      call serialize_char_2d(char_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize logical columns using direct file operations
    if (size(logical_cols, 1) > 0 .and. size(logical_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_logical_cols.dat'
      call serialize_logical_2d(logical_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize complex columns using direct file operations  
    if (size(complex_cols, 1) > 0 .and. size(complex_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_complex_cols.dat'
      call serialize_complex_2d(complex_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize header (column names)
    if (size(header) > 0) then
      filename = trim(filename_prefix) // '_header.dat'
      call serialize_char_1d(header, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize metadata
    filename = trim(filename_prefix) // '_metadata.dat'
    call serialize_int_2d(metadata, filename, local_ierr)
    if (is_err(local_ierr)) then
      call set_err_once(ierr, local_ierr)
      return
    end if
    
  end subroutine serialize_table

   ! //TODO Will be replaced by array_utils implementations in future releases
  !> Serialize a 2D logical array to a binary file using the array_utils pattern
  subroutine serialize_logical_2d(arr, filename, ierr)
    implicit none
    logical, intent(in) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: unit, ioerror
    integer(int32) :: dims(2)
    integer(int32), parameter :: ARRAY_TYPE_LOGICAL = 4
    
    dims = shape(arr)
    call set_ok(ierr)
    call set_ok(ioerror)
    
    ! Use array_utils to write header
    call write_file_header(filename, unit, ARRAY_TYPE_LOGICAL, 2, dims, ierr)
    if (is_err(ierr)) return
    
    ! Write the logical array data
    write(unit, iostat=ioerror) arr
    if (is_err(ioerror)) then
      call set_err_once(ierr, ioerror)
    end if
    close(unit)
  end subroutine serialize_logical_2d

  ! //TODO Will be replaced by array_utils implementations in future releases
  !> Serialize a 2D complex array to a binary file using the array_utils pattern
  subroutine serialize_complex_2d(arr, filename, ierr)
    implicit none
    complex(real64), intent(in) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: unit, ioerror
    integer(int32) :: dims(2)
    integer(int32), parameter :: ARRAY_TYPE_COMPLEX = 5
    
    dims = shape(arr)
    call set_ok(ierr)
    call set_ok(ioerror)
    
    ! Use array_utils to write header
    call write_file_header(filename, unit, ARRAY_TYPE_COMPLEX, 2, dims, ierr)
    if (is_err(ierr)) return
    
    ! Write the complex array data
    write(unit, iostat=ioerror) arr
    if (is_err(ioerror)) then
      call set_err_once(ierr, ioerror)
    end if
    close(unit)
  end subroutine serialize_complex_2d


end module tox_csv_table_serialization