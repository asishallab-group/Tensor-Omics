!> Module for reading heterogeneous CSV tables into type-banded arrays and providing fast accessors.
module tox_csv_file_reader
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_FILE_OPEN, ERR_FILE_EMPTY, ERR_DIM_MISMATCH, set_ok, set_err_once, is_ok
  use serialize_int, only: serialize_int_2d
  use serialize_real, only: serialize_real_2d
  use serialize_char, only: serialize_char_1d, serialize_char_2d
  use array_utils, only: write_file_header, read_file_header
  use int_deserialize_mod, only: deserialize_int_1d, deserialize_int_2d
  use real_deserialize_mod, only: deserialize_real_1d, deserialize_real_2d
  use char_deserialize_mod, only: deserialize_char_1d, deserialize_char_2d

contains

  !> Reads a table from a CSV file into type-banded arrays
  subroutine read_table(file_path, column_types, has_header, int_cols, real_cols, char_cols, &
                             logical_cols, complex_cols, header, metadata, ierr, sep_char, column_names)
    implicit none

    !| Path to the CSV file
    character(len=*), intent(in) :: file_path
    !| 1D int array, type for each column (1=int, 2=real, 3=char, 4=logical, 5=complex)
    integer(int32), intent(in) :: column_types(:)
    !| Logical, true if first line contains column names
    logical, intent(in) :: has_header
    !| 2D int array, output for integer columns
    integer(int32), intent(out) :: int_cols(:,:)
    !| 2D real array, output for real columns
    real(real64), intent(out) :: real_cols(:,:)
    !| 2D char array, output for character columns
    character(len=*), intent(out) :: char_cols(:,:)
    !| 2D logical array, output for logical columns
    logical, intent(out) :: logical_cols(:,:)
    !| 2D complex array, output for complex columns
    complex(real64), intent(out) :: complex_cols(:,:)
    !| 1D char array, output column names
    character(len=*), intent(out) :: header(:)
    !| 2D int array, output metadata (row 1: type, row 2: index in type array)
    integer(int32), intent(out) :: metadata(:,:)
    !| Error code: 0 - success, non-zero = error
    integer(int32), intent(out) :: ierr
    !| Optional character, column separator (default ',')
    character(len=*), intent(in), optional :: sep_char
    !| Optional 1D char array, overrides header line if provided
    character(len=*), intent(in), optional :: column_names(:)

    ! // TODO: Implement quotes handling for complex numbers
    ! // TODO: Check different reading approach: insted of large line buffer, use stream input of file for continous reading

    !| Local variables
    integer(int32) :: current_char, current_row, current_column, n_columns, col_type, type_index, file_unit, data_offset, io_err
    integer(int32) :: i_col_int, i_col_real, i_col_char, i_col_logical, i_col_complex, current_pos, sep_pos, close_paren_pos
    character(len=8192) :: line ! //TODO what size should this be? maybe as input argument?
    character(len=64) :: field ! //TODO what size should this be? maybe as input argument?
    character(len=1) :: sep

    ! Initialize error codes
    call set_ok(ierr)
    call set_ok(io_err)
    
    ! Error handling for dimension mismatches
    n_columns = size(column_types)
    
    ! Check header array dimension
    if (size(header) /= n_columns) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if
    
    ! Check metadata array dimensions
    if (size(metadata, 1) /= 2 .or. size(metadata, 2) /= n_columns) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if
    
    ! Check column_names dimension if present
    if (present(column_names)) then
      if (size(column_names) /= n_columns) then
        call set_err_once(ierr, ERR_DIM_MISMATCH)
        return
      end if
    end if

    ! Set separator (optional argument)
    if (present(sep_char)) then
      sep = sep_char(1:1)
    else
      sep = ','
    end if
    
    ! Open file
    file_unit = 10
    open(unit=file_unit, file=file_path, status='old', action='read', iostat=io_err)
    if (.not. is_ok(io_err)) then
      call set_err_once(ierr, ERR_FILE_OPEN)
      return
    end if

    ! Check for empty file.
    read (file_unit, '(A)', iostat=io_err) line
    if (.not. is_ok(io_err)) then
      call set_err_once(ierr, ERR_FILE_EMPTY)
      close (file_unit)
      return
    end if

    ! Go back to the beginning of the file.
    rewind(file_unit)

    ! Initialize data offset (number of lines to skip at the start of the file)
    data_offset = 0

    ! Skip initial comment lines starting with '#'.
    do
      read (file_unit, '(A)', iostat=io_err) line
      if (.not. is_ok(io_err)) exit
      if (len_trim(line) > 0) then
        if (line(1:1) /= '#') exit
      end if
      data_offset = data_offset + 1
    end do

    ! Read header or use provided column_names
    if (present(column_names)) then
      ! Use provided column names
      header = column_names
    end if
    if (has_header) then
      data_offset = data_offset + 1

      if(.not. present(column_names)) then
        ! Split header line by separator
        n_columns = size(column_types)
        current_char = 1
        do current_column = 1, n_columns
          ! Find next separator or end of line
          type_index = index(line(current_char:), sep)
          ! If no separator found, take rest of line
          if (type_index == 0) then
            header(current_column) = adjustl(trim(line(current_char:)))
            ! Check if we have fewer columns than expected
            if (current_column < n_columns) then
              call set_err_once(ierr, ERR_DIM_MISMATCH)
            end if
            exit
          else
            header(current_column) = adjustl(trim(line(current_char:current_char+type_index-2)))
            current_char = current_char + type_index
          end if
        end do
      end if
    end if

    ! Go back to the start of the file and skip first offset lines (comments and header if present)
    rewind(file_unit)
    do current_char = 1, data_offset
      read(file_unit, '(A)', iostat=io_err) line
      if (.not. is_ok(io_err)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        exit
      end if
    end do

    ! Read data rows and fill type-banded arrays
    metadata = 0
    n_columns = size(column_types)
    current_row = 0
    do while (.true.)
      ! Read the next line
      read(file_unit, '(A)', iostat=io_err) line
      if (.not. is_ok(io_err)) exit
      current_row = current_row + 1

      ! Reset column indices for each type
      i_col_int = 1
      i_col_real = 1
      i_col_char = 1
      i_col_logical = 1
      i_col_complex = 1
      ! Reset inline separator position
      current_pos = 1
      ! Go through each column in the line
      do current_column = 1, n_columns

        ! Find next separator or end of line
        ! First, extract the field normally to check for complex numbers
        sep_pos = index(line(current_pos:), sep)
        ! If no separator found, take rest of line
        if (sep_pos == 0) then
          field = adjustl(trim(line(current_pos:)))
          if (len_trim(field) == 0) then
            call set_err_once(ierr, ERR_INVALID_INPUT)  ! Empty field found
          end if
          if (current_column < n_columns) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)  ! Fewer columns than expected
          end if
        else
          field = adjustl(trim(line(current_pos:current_pos+sep_pos-2)))
        end if
        
        ! Special handling for complex numbers (field starts with '(' after trimming)
        if (len_trim(field) > 0 .and. field(1:1) == '(') then
          ! For complex numbers, find closing ')' instead of separator
          close_paren_pos = index(line(current_pos:), ')')
          if (close_paren_pos > 0) then
            ! Include the closing parenthesis in the field
            field = adjustl(trim(line(current_pos:current_pos+close_paren_pos-1)))
            ! Find the actual separator after the closing parenthesis
            sep_pos = index(line(current_pos+close_paren_pos:), sep)
            if (sep_pos > 0) then
              sep_pos = sep_pos + close_paren_pos
            else
              sep_pos = 0  ! No separator after complex number
            end if
          end if
        end if

        ! Read column type to store value in appropriate array
        col_type = column_types(current_column)

        select case (col_type)
        case (1) ! int
          read(field, *, iostat=io_err) int_cols(current_row, i_col_int)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 1
            metadata(2, current_column) = i_col_int
            ! Generate header name if none provided
            if (.not. has_header .and. .not. present(column_names)) then
              write(header(current_column), '(I0,A)') current_column, '_int'
            end if
          end if
          ! Increment column index for next int column
          i_col_int = i_col_int + 1
        
        case (2) ! real
          read(field, *, iostat=io_err) real_cols(current_row, i_col_real)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 2
            metadata(2, current_column) = i_col_real
            ! Generate header name if none provided
            if (.not. has_header .and. .not. present(column_names)) then
              write(header(current_column), '(I0,A)') current_column, '_real'
            end if
          end if
          ! Increment column index for next real column
          i_col_real = i_col_real + 1

        case (3) ! char
          read(field, *, iostat=io_err) char_cols(current_row, i_col_char)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 3
            metadata(2, current_column) = i_col_char
            ! Generate header name if none provided
            if (.not. has_header .and. .not. present(column_names)) then
              write(header(current_column), '(I0,A)') current_column, '_char'
            end if
          end if
          ! Increment column index for next char column
          i_col_char = i_col_char + 1

        case (4) ! logical
          ! Read logical column 
          ! (fortran interprets T, t, TRUE, true, .true. or 1 as .true.)
          ! (fortran interprets F, f, FALSE, false, .false. or 0 as .false.)
          read(field, *, iostat=io_err) logical_cols(current_row, i_col_logical)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 4
            metadata(2, current_column) = i_col_logical
            ! Generate header name if none provided
            if (.not. has_header .and. .not. present(column_names)) then
              write(header(current_column), '(I0,A)') current_column, '_logical'
            end if
          end if
          ! Increment column index for next logical column
          i_col_logical = i_col_logical + 1

        case (5) ! complex
          ! Complex numbers are expected in the form (a,b) for a + bi
          read(field, *, iostat=io_err) complex_cols(current_row, i_col_complex)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 5
            metadata(2, current_column) = i_col_complex
            ! Generate header name if none provided
            if (.not. has_header .and. .not. present(column_names)) then
              write(header(current_column), '(I0,A)') current_column, '_complex'
            end if
          end if
          ! Increment column index for next complex column
          i_col_complex = i_col_complex + 1

        end select

        ! Increment current position for next field if separator found
        if (sep_pos == 0) then
          exit
        else
          current_pos = current_pos + sep_pos
        end if

      end do
    end do

    close(file_unit)

  end subroutine read_table

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

  !> Serialization and deserialization routines using serial_array_module
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
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize real columns if present
    if (size(real_cols, 1) > 0 .and. size(real_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_real_cols.dat'
      call serialize_real_2d(real_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize character columns if present
    if (size(char_cols, 1) > 0 .and. size(char_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_char_cols.dat'
      call serialize_char_2d(char_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize logical columns using direct file operations
    if (size(logical_cols, 1) > 0 .and. size(logical_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_logical_cols.dat'
      call serialize_logical_2d(logical_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize complex columns using direct file operations  
    if (size(complex_cols, 1) > 0 .and. size(complex_cols, 2) > 0) then
      filename = trim(filename_prefix) // '_complex_cols.dat'
      call serialize_complex_2d(complex_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize header (column names)
    if (size(header) > 0) then
      filename = trim(filename_prefix) // '_header.dat'
      call serialize_char_1d(header, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Serialize metadata
    filename = trim(filename_prefix) // '_metadata.dat'
    call serialize_int_2d(metadata, filename, local_ierr)
    if (.not. is_ok(local_ierr)) then
      call set_err_once(ierr, local_ierr)
      return
    end if
    
  end subroutine serialize_table

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
    if (.not. is_ok(local_ierr)) then
      call set_err_once(ierr, local_ierr)
      return
    end if
    
    ! Deserialize header (column names)
    filename = trim(filename_prefix) // '_header.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists) then
      call deserialize_char_1d(header, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize integer columns if file exists
    filename = trim(filename_prefix) // '_int_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(int_cols, 1) > 0 .and. size(int_cols, 2) > 0) then
      call deserialize_int_2d(int_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize real columns if file exists
    filename = trim(filename_prefix) // '_real_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(real_cols, 1) > 0 .and. size(real_cols, 2) > 0) then
      call deserialize_real_2d(real_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize character columns if file exists
    filename = trim(filename_prefix) // '_char_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(char_cols, 1) > 0 .and. size(char_cols, 2) > 0) then
      call deserialize_char_2d(char_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize logical columns if file exists
    filename = trim(filename_prefix) // '_logical_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(logical_cols, 1) > 0 .and. size(logical_cols, 2) > 0) then
      call deserialize_logical_2d(logical_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    ! Deserialize complex columns if file exists  
    filename = trim(filename_prefix) // '_complex_cols.dat'
    inquire(file=filename, exist=file_exists)
    if (file_exists .and. size(complex_cols, 1) > 0 .and. size(complex_cols, 2) > 0) then
      call deserialize_complex_2d(complex_cols, filename, local_ierr)
      if (.not. is_ok(local_ierr)) then
        call set_err_once(ierr, local_ierr)
        return
      end if
    end if
    
    end subroutine deserialize_table

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
    if (.not. is_ok(ierr)) return
    
    ! Write the logical array data
    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
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
    if (.not. is_ok(ierr)) return
    
    ! Write the complex array data
    write(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ioerror)
    end if
    close(unit)
  end subroutine serialize_complex_2d

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
    if (.not. is_ok(ierr)) return
    
    ! Verify dimensions match
    if (ndims /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      close(unit)
      return
    end if
    
    ! Read the logical array data
    read(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
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
    if (.not. is_ok(ierr)) return
    
    ! Verify dimensions match
    if (ndims /= 2) then
      call set_err_once(ierr, ERR_DIM_MISMATCH)
      close(unit)
      return
    end if
    
    ! Read the complex array data
    read(unit, iostat=ioerror) arr
    if (.not. is_ok(ioerror)) then
      call set_err_once(ierr, ioerror)
    end if
    close(unit)
  end subroutine deserialize_complex_2d

end module tox_csv_file_reader