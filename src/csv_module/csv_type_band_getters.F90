!> Module for CSV band type getters.
module tox_csv_type_band_getters
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: set_ok, set_err_once, ERR_INVALID_INPUT

  contains

  !> Get integer column by original table index
  pure subroutine get_int_column_by_index(int_cols, metadata, index, single_int_column, ierr)
    implicit none

    !| 2D int array, integer columns from read_table
    integer(int32), intent(in) :: int_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| Original table column index
    integer(int32), intent(in) :: index
    !| 1D int array, extracted column
    integer(int32), intent(out) :: single_int_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: type_band_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Check if index is valid
    if (index < 1 .or. index > size(metadata, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Check if the column at this index is actually an integer column
    if (metadata(1, index) /= 1) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column is not an integer type
      return
    end if
    
    ! Get the type-band index for this column
    type_band_index = metadata(2, index)
    
    ! Check if type-band index is valid
    if (type_band_index < 1 .or. type_band_index > size(int_cols, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Extract the column from int_cols
    single_int_column = int_cols(:, type_band_index)
    
  end subroutine get_int_column_by_index

  !> Get integer column by column name
  pure subroutine get_int_column_by_name(int_cols, metadata, header, name, single_int_column, ierr)
    implicit none

    !| 2D int array, integer columns from read_table
    integer(int32), intent(in) :: int_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| 1D char array, column names from read_table
    character(len=*), intent(in) :: header(:)
    !| Column name to search for
    character(len=*), intent(in) :: name
    !| 1D int array, extracted column
    integer(int32), intent(out) :: single_int_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: i, found_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Search for the column name in header
    found_index = 0
    do i = 1, size(header)
      if (trim(adjustl(header(i))) == trim(adjustl(name))) then
        found_index = i
        exit
      end if
    end do
    
    ! Check if column name was found
    if (found_index == 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column name not found
      return
    end if
    
    ! Call the by_index version with the found index
    call get_int_column_by_index(int_cols, metadata, found_index, single_int_column, ierr)
    
  end subroutine get_int_column_by_name

  !> Get real column by original table index
  pure subroutine get_real_column_by_index(real_cols, metadata, index, single_real_column, ierr)
    implicit none

    !| 2D real array, real columns from read_table
    real(real64), intent(in) :: real_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| Original table column index
    integer(int32), intent(in) :: index
    !| 1D real array, extracted column
    real(real64), intent(out) :: single_real_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: type_band_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Check if index is valid
    if (index < 1 .or. index > size(metadata, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Check if the column at this index is actually a real column
    if (metadata(1, index) /= 2) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column is not a real type
      return
    end if
    
    ! Get the type-band index for this column
    type_band_index = metadata(2, index)
    
    ! Check if type-band index is valid
    if (type_band_index < 1 .or. type_band_index > size(real_cols, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Extract the column from real_cols
    single_real_column = real_cols(:, type_band_index)
    
  end subroutine get_real_column_by_index

  !> Get real column by column name
  pure subroutine get_real_column_by_name(real_cols, metadata, header, name, single_real_column, ierr)
    implicit none

    !| 2D real array, real columns from read_table
    real(real64), intent(in) :: real_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| 1D char array, column names from read_table
    character(len=*), intent(in) :: header(:)
    !| Column name to search for
    character(len=*), intent(in) :: name
    !| 1D real array, extracted column
    real(real64), intent(out) :: single_real_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: i, found_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Search for the column name in header
    found_index = 0
    do i = 1, size(header)
      if (trim(adjustl(header(i))) == trim(adjustl(name))) then
        found_index = i
        exit
      end if
    end do
    
    ! Check if column name was found
    if (found_index == 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column name not found
      return
    end if
    
    ! Call the by_index version with the found index
    call get_real_column_by_index(real_cols, metadata, found_index, single_real_column, ierr)
    
  end subroutine get_real_column_by_name

  !> Get character column by original table index
  pure subroutine get_char_column_by_index(char_cols, metadata, index, single_char_column, ierr)
    implicit none

    !| 2D char array, character columns from read_table
    character(len=*), intent(in) :: char_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| Original table column index
    integer(int32), intent(in) :: index
    !| 1D char array, extracted column
    character(len=*), intent(out) :: single_char_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: type_band_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Check if index is valid
    if (index < 1 .or. index > size(metadata, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Check if the column at this index is actually a character column
    if (metadata(1, index) /= 3) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column is not a character type
      return
    end if
    
    ! Get the type-band index for this column
    type_band_index = metadata(2, index)
    
    ! Check if type-band index is valid
    if (type_band_index < 1 .or. type_band_index > size(char_cols, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Extract the column from char_cols
    single_char_column = char_cols(:, type_band_index)
    
  end subroutine get_char_column_by_index

  !> Get character column by column name
  pure subroutine get_char_column_by_name(char_cols, metadata, header, name, single_char_column, ierr)
    implicit none

    !| 2D char array, character columns from read_table
    character(len=*), intent(in) :: char_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| 1D char array, column names from read_table
    character(len=*), intent(in) :: header(:)
    !| Column name to search for
    character(len=*), intent(in) :: name
    !| 1D char array, extracted column
    character(len=*), intent(out) :: single_char_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: i, found_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Search for the column name in header
    found_index = 0
    do i = 1, size(header)
      if (trim(adjustl(header(i))) == trim(adjustl(name))) then
        found_index = i
        exit
      end if
    end do
    
    ! Check if column name was found
    if (found_index == 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column name not found
      return
    end if
    
    ! Call the by_index version with the found index
    call get_char_column_by_index(char_cols, metadata, found_index, single_char_column, ierr)
    
  end subroutine get_char_column_by_name

  !> Get logical column by original table index
  pure subroutine get_logical_column_by_index(logical_cols, metadata, index, single_logical_column, ierr)
    implicit none

    !| 2D logical array, logical columns from read_table
    logical, intent(in) :: logical_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| Original table column index
    integer(int32), intent(in) :: index
    !| 1D logical array, extracted column
    logical, intent(out) :: single_logical_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: type_band_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Check if index is valid
    if (index < 1 .or. index > size(metadata, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Check if the column at this index is actually a logical column
    if (metadata(1, index) /= 4) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column is not a logical type
      return
    end if
    
    ! Get the type-band index for this column
    type_band_index = metadata(2, index)
    
    ! Check if type-band index is valid
    if (type_band_index < 1 .or. type_band_index > size(logical_cols, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Extract the column from logical_cols
    single_logical_column = logical_cols(:, type_band_index)
    
  end subroutine get_logical_column_by_index

  !> Get logical column by column name
  pure subroutine get_logical_column_by_name(logical_cols, metadata, header, name, single_logical_column, ierr)
    implicit none

    !| 2D logical array, logical columns from read_table
    logical, intent(in) :: logical_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| 1D char array, column names from read_table
    character(len=*), intent(in) :: header(:)
    !| Column name to search for
    character(len=*), intent(in) :: name
    !| 1D logical array, extracted column
    logical, intent(out) :: single_logical_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: i, found_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Search for the column name in header
    found_index = 0
    do i = 1, size(header)
      if (trim(adjustl(header(i))) == trim(adjustl(name))) then
        found_index = i
        exit
      end if
    end do
    
    ! Check if column name was found
    if (found_index == 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column name not found
      return
    end if
    
    ! Call the by_index version with the found index
    call get_logical_column_by_index(logical_cols, metadata, found_index, single_logical_column, ierr)
    
  end subroutine get_logical_column_by_name

  !> Get complex column by original table index
  pure subroutine get_complex_column_by_index(complex_cols, metadata, index, single_complex_column, ierr)
    implicit none

    !| 2D complex array, complex columns from read_table
    complex(real64), intent(in) :: complex_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| Original table column index
    integer(int32), intent(in) :: index
    !| 1D complex array, extracted column
    complex(real64), intent(out) :: single_complex_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: type_band_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Check if index is valid
    if (index < 1 .or. index > size(metadata, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Check if the column at this index is actually a complex column
    if (metadata(1, index) /= 5) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column is not a complex type
      return
    end if
    
    ! Get the type-band index for this column
    type_band_index = metadata(2, index)
    
    ! Check if type-band index is valid
    if (type_band_index < 1 .or. type_band_index > size(complex_cols, 2)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Extract the column from complex_cols
    single_complex_column = complex_cols(:, type_band_index)
    
  end subroutine get_complex_column_by_index

  !> Get complex column by column name
  pure subroutine get_complex_column_by_name(complex_cols, metadata, header, name, single_complex_column, ierr)
    implicit none

    !| 2D complex array, complex columns from read_table
    complex(real64), intent(in) :: complex_cols(:,:)
    !| 2D int array, metadata from read_table
    integer(int32), intent(in) :: metadata(:,:)
    !| 1D char array, column names from read_table
    character(len=*), intent(in) :: header(:)
    !| Column name to search for
    character(len=*), intent(in) :: name
    !| 1D complex array, extracted column
    complex(real64), intent(out) :: single_complex_column(:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    
    ! Local variables
    integer(int32) :: i, found_index
    
    ! Initialize error code
    call set_ok(ierr)
    
    ! Search for the column name in header
    found_index = 0
    do i = 1, size(header)
      if (trim(adjustl(header(i))) == trim(adjustl(name))) then
        found_index = i
        exit
      end if
    end do
    
    ! Check if column name was found
    if (found_index == 0) then
      call set_err_once(ierr, ERR_INVALID_INPUT)  ! Column name not found
      return
    end if
    
    ! Call the by_index version with the found index
    call get_complex_column_by_index(complex_cols, metadata, found_index, single_complex_column, ierr)
    
  end subroutine get_complex_column_by_name

end module tox_csv_type_band_getters