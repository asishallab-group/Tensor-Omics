!> Module for deserializing CSV tables from binary files
module tox_csv_table_deserialization
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_FILE_OPEN, ERR_FILE_EMPTY, ERR_DIM_MISMATCH, set_ok, set_err_once, is_err
  use int_deserialize_mod, only: deserialize_int_2d
  use real_deserialize_mod, only: deserialize_real_2d
  use char_deserialize_mod, only: deserialize_char_1d, deserialize_char_2d
  use array_utils, only: read_file_header

contains

!// TODO: error handling for empty file

!> Deserialize type-banded arrays from binary files back into CSV table format
  subroutine deserialize_table(filename_prefix, int_cols, real_cols, char_cols, &
                              logical_cols, complex_cols, header, metadata, ierr)
    implicit none
    
    !| Base filename prefix (extensions will be added automatically)
    character(len=*), intent(in) :: filename_prefix
    !| 2D int array, output for deserialized integer columns
    integer(int32), intent(out) :: int_cols(:,:)
    !| 2D real array, output for deserialized real columns  
    real(real64), intent(out) :: real_cols(:,:)
    !| 2D char array, output for deserialized character columns
    character(len=*), intent(out) :: char_cols(:,:)
    !| 2D logical array, output for deserialized logical columns
    logical, intent(out) :: logical_cols(:,:)
    !| 2D complex array, output for deserialized complex columns
    complex(real64), intent(out) :: complex_cols(:,:)
    !| 1D char array, output for column names/headers
    character(len=*), intent(out) :: header(:)
    !| 2D int array, output for metadata (row 1: type, row 2: index in type array)
    integer(int32), intent(out) :: metadata(:,:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    character(len=1024) :: filename
    integer(int32) :: local_ierr
    logical :: file_exists
    
    call set_ok(ierr)
    call set_ok(local_ierr)
    
    ! First, always deserialize metadata to understand the table structure
    filename = trim(filename_prefix) // '_metadata.dat'
    inquire(file=filename, exist=file_exists)
    if (.not. file_exists) then
      call set_err_once(ierr, ERR_FILE_OPEN)
      return
    end if
    call deserialize_int_2d(metadata, filename, local_ierr)
    if (is_err(local_ierr)) then
      call set_err_once(ierr, local_ierr)
      return
    end if
    
    ! Deserialize header (column names)
    filename = trim(filename_prefix) // '_header.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      call deserialize_char_1d(header, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize integer columns if file exists
    filename = trim(filename_prefix) // '_int_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(int_cols, 1) > 0 .and. size(int_cols, 2) > 0) then
      call deserialize_int_2d(int_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize real columns if file exists
    filename = trim(filename_prefix) // '_real_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(real_cols, 1) > 0 .and. size(real_cols, 2) > 0) then
      call deserialize_real_2d(real_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize character columns if file exists
    filename = trim(filename_prefix) // '_char_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(char_cols, 1) > 0 .and. size(char_cols, 2) > 0) then
      call deserialize_char_2d(char_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize logical columns if file exists
    filename = trim(filename_prefix) // '_logical_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(logical_cols, 1) > 0 .and. size(logical_cols, 2) > 0) then
      call deserialize_logical_2d(logical_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize complex columns if file exists  
    filename = trim(filename_prefix) // '_complex_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(complex_cols, 1) > 0 .and. size(complex_cols, 2) > 0) then
      call deserialize_complex_2d(complex_cols, filename, local_ierr)
      if (is_err(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    end subroutine deserialize_table


  ! //TODO Will be replaced by array_utils implementations in future releases
  !> Deserialize a 2D logical array from a binary file using the array_utils pattern
  subroutine deserialize_logical_2d(arr, filename, ierr)
    implicit none
    logical, intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)
    
    call set_ok(ierr)
    
    ! Use array_utils to read header
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (is_err(ierr)) return
    
    ! Verify dimensions match
    if (ndims /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      close(unit)
      return
    end if
    
    ! Read the logical array data
    read(unit, iostat=ioerror) arr
    if (is_err(ioerror)) then
      call set_err_once(ierr, ioerror)
    end if
    close(unit)
  end subroutine deserialize_logical_2d

  ! //TODO Will be replaced by array_utils implementations in future releases
  !> Deserialize a 2D complex array from a binary file using the array_utils pattern
  subroutine deserialize_complex_2d(arr, filename, ierr)
    implicit none
    complex(real64), intent(out) :: arr(:,:)
    character(len=*), intent(in) :: filename
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: unit, type_code, ndims, clen, ioerror
    integer(int32), allocatable :: dims(:)
    
    call set_ok(ierr)
    
    ! Use array_utils to read header
    call read_file_header(filename, unit, type_code, ndims, dims, clen, ierr)
    if (is_err(ierr)) return
    
    ! Verify dimensions match
    if (ndims /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      close(unit)
      return
    end if
    
    ! Read the complex array data
    read(unit, iostat=ioerror) arr
    if (is_err(ioerror)) then
      call set_err_once(ierr, ioerror)
    end if
    close(unit)
  end subroutine deserialize_complex_2d

end module tox_csv_table_deserialization