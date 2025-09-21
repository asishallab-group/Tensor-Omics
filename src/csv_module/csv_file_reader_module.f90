!> Module for reading heterogeneous CSV tables into type-banded arrays and providing fast accessors.
module tox_csv_file_reader
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_FILE_OPEN, ERR_FILE_EMPTY, ERR_DIM_MISMATCH, set_ok, set_err_once, is_ok

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

    !| Local variables
    integer(int32) :: current_row, current_column, n_columns, n_rows, col_type, type_index, file_unit, data_offset, io_err, i_col_int, i_col_real, i_col_char, i_col_logical, i_col_complex, current_pos, sep_pos, close_paren_pos, k, l
    character(len=1024) :: line, field
    character(len=1) :: sep

    ! Initialize error code
    call set_ok(ierr)
    
    ! Error handling for dimension mismatches
    n_columns = size(column_types)
    
    ! Check header array dimension
    if (size(header) /= n_columns) then
      set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if
    
    ! Check metadata array dimensions
    if (size(metadata, 1) /= 2 .or. size(metadata, 2) /= n_columns) then
      set_err_once(ierr, ERR_DIM_MISMATCH)
      return
    end if
    
    ! Check column_names dimension if present
    if (present(column_names)) then
      if (size(column_names) /= n_columns) then
        set_err_once(ierr, ERR_DIM_MISMATCH)
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
    open(unit=file_unit, file=file_path, status='old', action='read', iostat=io_err)
    if (.not. is_ok(io_err)) then
      set_err_once(ierr, ERR_FILE_OPEN)
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
    else if (has_header) then
      data_offset = data_offset + 1
      ! Split header line by separator
      n_columns = size(column_types)
      current_column = 1
      current_row = 1
      do while (current_column <= n_columns)
        ! Find next separator or end of line
        type_index = index(line(current_row:), sep)
        if (type_index == 0) then
          header(current_column) = adjustl(trim(line(current_row:)))
          exit
        else
          header(current_column) = adjustl(trim(line(current_row:current_row+type_index-2)))
          current_row = current_row + type_index
        end if
        current_column = current_column + 1
      end do
    end if

    ! Go back to the start of the file and skip first offset lines (comments and header if present)
    print *, "Data offset: ", data_offset
    rewind(file_unit)
    do current_row = 1, data_offset
      read(file_unit, '(A)', iostat=io_err) line
      if (.not. is_ok(io_err)) then
        set_err_once(ierr, ERR_INVALID_INPUT)
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
        if (sep_pos == 0) then
          field = adjustl(trim(line(current_pos:)))
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

  !> Get real column by original table index
  pure subroutine get_real_column(real_cols, metadata, index, single_real_column, ierr)
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
    
  end subroutine get_real_column

  !> Get integer column by original table index
  pure subroutine get_int_column(int_cols, metadata, index, single_int_column, ierr)
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
    set_ok(ierr)
    
    ! Check if index is valid
    if (index < 1 .or. index > size(metadata, 2)) then
      set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Check if the column at this index is actually an integer column
    if (metadata(1, index) /= 1) then
      set_err_once(ierr, ERR_INVALID_INPUT)  ! Column is not an integer type
      return
    end if
    
    ! Get the type-band index for this column
    type_band_index = metadata(2, index)
    
    ! Check if type-band index is valid
    if (type_band_index < 1 .or. type_band_index > size(int_cols, 2)) then
      set_err_once(ierr, ERR_INVALID_INPUT)
      return
    end if
    
    ! Extract the column from int_cols
    single_int_column = int_cols(:, type_band_index)
    
  end subroutine get_int_column

  !> Get character column by original table index
  pure subroutine get_char_column(char_cols, metadata, index, single_char_column, ierr)
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
    
  end subroutine get_char_column

  !> Get logical column by original table index
  pure subroutine get_logical_column(logical_cols, metadata, index, single_logical_column, ierr)
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
    
  end subroutine get_logical_column

  !> Get complex column by original table index
  pure subroutine get_complex_column(complex_cols, metadata, index, single_complex_column, ierr)
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
    
  end subroutine get_complex_column

  !> Serialization and deserialization routines using serial_array_module
  ! pure subroutine serialize_table()
  !   implicit none
  !   ! ...implementation to be added...
  ! end subroutine serialize_table

  ! pure subroutine deserialize_table()
  !   implicit none
  !   ! ...implementation to be added...
  ! end subroutine deserialize_table

  !> Serialization and deserialization routines using serial_array_module
  !pure subroutine serialize_table(...)
  !  implicit none
    ! ...existing code...
  !end subroutine serialize_table

  !pure subroutine deserialize_table(...)
  !  implicit none
    ! ...existing code...
  !end subroutine deserialize_table

end module tox_csv_file_reader
