!> Module for reading heterogeneous CSV tables into type-banded arrays and providing fast accessors.
module tox_csv_file_reader
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_FILE_OPEN, ERR_FILE_EMPTY, ERR_DIM_MISMATCH, set_ok, set_err_once, is_ok, is_err
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

    !| Local variables
    integer(int32) :: current_row, current_column, n_columns, col_type, type_index, file_unit, io_err
    integer(int32) :: i_col_int, i_col_real, i_col_char, i_col_logical, i_col_complex, stream_pos
    character(len=256) :: field
    character(len=1) :: sep, char_read
    logical :: end_of_file, end_of_line, in_complex

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
      if (len_trim(sep_char) == 0) then
        sep = ','
      else if (trim(sep_char) == '\t') then
        sep = char(9)  ! ASCII code 9 is TAB
      else 
        sep = sep_char(1:1)
      end if
    else
      sep = ','
    end if
    
    ! Open file as stream for character-by-character reading
    file_unit = 10
    open(unit=file_unit, file=file_path, status='old', action='read', access='stream', iostat=io_err)
    if (is_err(io_err)) then
      call set_err_once(ierr, ERR_FILE_OPEN)
      return
    end if

    stream_pos = 1  ! Position in stream for reading

    ! Check for empty file by trying to read first character
    read(file_unit, POS=stream_pos, iostat=io_err) char_read
    if (is_err(io_err)) then
      call set_err_once(ierr, ERR_FILE_EMPTY)
      close(file_unit)
      return
    end if

    ! Skip comment lines starting with '#'
    call skip_comment_lines(file_unit, stream_pos, io_err)
    if (is_err(io_err)) then
      call set_err_once(ierr, ERR_INVALID_INPUT)
      close(file_unit)
      return
    end if

    ! Read header if present
    if (present(column_names)) then
      ! Use provided column names
      header = column_names
      ! Skip header line if it exists in file
      if (has_header) then
        call skip_line(file_unit, stream_pos, io_err)
        if (is_err(io_err)) then
          call set_err_once(ierr, ERR_INVALID_INPUT)
          close(file_unit)
          return
        end if
      end if
    else if (has_header) then
      ! Read header from file
      call read_header_line(file_unit, sep, stream_pos, header, n_columns, io_err)
      if (is_err(io_err)) then
        call set_err_once(ierr, ERR_INVALID_INPUT)
        close(file_unit)
        return
      end if
    end if

    ! Read data rows and fill type-banded arrays
    metadata = 0
    n_columns = size(column_types)
    current_row = 1
    end_of_file = .false.
    
    do while (.not. end_of_file)
      
      ! Reset column indices for each type
      i_col_int = 1
      i_col_real = 1
      i_col_char = 1
      i_col_logical = 1
      i_col_complex = 1
      
      ! Read each field in the current row
      do current_column = 1, n_columns
        call read_field(file_unit, sep, stream_pos, field, end_of_line, end_of_file, io_err)
        if (is_err(io_err) .or. end_of_file) then
          ! End of file at start of line is normal
          if (current_column == 1) then 
            exit
          else
            ! End of file in middle of line is an error
            call set_err_once(ierr, ERR_INVALID_INPUT)
            close(file_unit)
            return
          end if
        end if

        ! Check for empty field
        if (len_trim(field) == 0) then
          call set_err_once(ierr, ERR_INVALID_INPUT)  ! Empty field found
          close(file_unit)
          return
        end if

        ! Read column type to store value in appropriate array
        col_type = column_types(current_column)

        ! Read field into appropriate type-banded array and increment column index
        select case (col_type)
        case (1) ! int
          read(field, *, iostat=io_err) int_cols(current_row, i_col_int)
          type_index = i_col_int
          i_col_int = i_col_int + 1
        case (2) ! real
          read(field, *, iostat=io_err) real_cols(current_row, i_col_real)
          type_index = i_col_real
          i_col_real = i_col_real + 1
        case (3) ! char
          read(field, *, iostat=io_err) char_cols(current_row, i_col_char)
          type_index = i_col_char
          i_col_char = i_col_char + 1
        case (4) ! logical
          ! Read logical column 
          ! (interprets T, t, TRUE, true, .true. or 1 as .true.)
          ! (interprets F, f, FALSE, false, .false. or 0 as .false.)
          read(field, *, iostat=io_err) logical_cols(current_row, i_col_logical)
          type_index = i_col_logical
          i_col_logical = i_col_logical + 1
        case (5) ! complex
          ! Complex numbers are expected in the form (a,b) for a + bi
          read(field, *, iostat=io_err) complex_cols(current_row, i_col_complex)
          type_index = i_col_complex
          i_col_complex = i_col_complex + 1
        end select

        ! Check for read error
        if (is_err(io_err)) then
          call set_err_once(ierr, ERR_INVALID_INPUT)
          close(file_unit)
          return
        end if

        ! Set metadata and generate header if first time this column type is encountered
        if (metadata(1, current_column) == 0) then
          metadata(1, current_column) = col_type
          metadata(2, current_column) = type_index
          ! Generate header name if none provided
          if (.not. has_header .and. .not. present(column_names)) then
            select case (col_type)
            case (1); write(header(current_column), '(I0,A)') current_column, '_int'
            case (2); write(header(current_column), '(I0,A)') current_column, '_real'
            case (3); write(header(current_column), '(I0,A)') current_column, '_char'
            case (4); write(header(current_column), '(I0,A)') current_column, '_logical'
            case (5); write(header(current_column), '(I0,A)') current_column, '_complex'
            end select
          end if
        end if

        ! Check if we've reached end of line
        if (end_of_line) then
          ! Verify we read all expected columns
          if (current_column < n_columns) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)  ! Fewer columns than expected
            close(file_unit)
            return
          end if
          exit
        end if
      end do
      
      ! If we read all columns but didn't hit end of line, there are too many columns
      if (.not. end_of_line .and. current_column == n_columns) then
        call set_err_once(ierr, ERR_DIM_MISMATCH)  ! More columns than expected
        close(file_unit)
        return
      end if

      current_row = current_row + 1
      stream_pos = stream_pos + 1
    end do

    close(file_unit)

  end subroutine read_table

  ! Helper subroutines for stream-based CSV reading

  !> Skip comment lines starting with '#' in stream
  subroutine skip_comment_lines(file_unit, stream_pos,io_err)
    implicit none
    integer(int32), intent(in) :: file_unit
    integer(int32), intent(inout) :: stream_pos
    integer(int32), intent(out) :: io_err
    
    character(len=1) :: char_read
    
    call set_ok(io_err)

    do
      read(file_unit, POS=stream_pos, iostat=io_err) char_read
      if (is_err(io_err)) return
      if (char_read == '#') then
        ! Skip rest of comment line
        call skip_line(file_unit, stream_pos, io_err)
        if (is_err(io_err)) return
      else
        return
      end if
    end do
  end subroutine skip_comment_lines

  !> Skip to end of current line
  subroutine skip_line(file_unit, stream_pos, io_err)
    implicit none
    integer(int32), intent(in) :: file_unit
    integer(int32), intent(out) :: io_err
    integer(int32), intent(inout) :: stream_pos

    character(len=1) :: char_read
    
    call set_ok(io_err)
    do
      read(file_unit, POS=stream_pos, iostat=io_err) char_read
      if (is_err(io_err)) return
      stream_pos = stream_pos + 1
      if (char_read == char(10)) then ! LF
        return
      else if (char_read == char(13)) then ! CR
        ! Check for CR-LF combination
        read(file_unit, POS=stream_pos, iostat=io_err) char_read
        if (is_ok(io_err) .and. char_read /= char(10)) then
          ! Not CR-LF, back up one character
          stream_pos = stream_pos - 1
        end if
        return
      end if
    end do
  end subroutine skip_line

  !> Read header line from stream
  subroutine read_header_line(file_unit, sep, stream_pos, header, n_columns, io_err)
    implicit none
    integer(int32), intent(in) :: file_unit
    character(len=1), intent(in) :: sep
    integer(int32), intent(inout) :: stream_pos
    character(len=*), intent(out) :: header(:)
    integer(int32), intent(in) :: n_columns
    integer(int32), intent(out) :: io_err
    
    character(len=256) :: field
    logical :: end_of_line, end_of_file
    integer(int32) :: col_idx
    
    call set_ok(io_err)
    
    do col_idx = 1, n_columns
      call read_field(file_unit, sep, stream_pos, field, end_of_line, end_of_file, io_err)
      if (is_err(io_err) .or. end_of_file) return
      header(col_idx) = adjustl(trim(field))
      
      if (end_of_line) then
        if (col_idx < n_columns) then
          call set_err_once(io_err, ERR_DIM_MISMATCH)
        end if
        stream_pos = stream_pos + 1
        return
      end if
    end do
  end subroutine read_header_line

  !> Read a single field from stream
  subroutine read_field(file_unit, sep, stream_pos, field, end_of_line, end_of_file, io_err)
    implicit none
    integer(int32), intent(in) :: file_unit
    character(len=1), intent(in) :: sep
    integer(int32), intent(inout) :: stream_pos
    character(len=*), intent(out) :: field
    logical, intent(out) :: end_of_line, end_of_file
    integer(int32), intent(out) :: io_err
    
    character(len=1) :: char_read
    integer(int32) :: field_pos, paren_count
    logical :: in_complex
    
    call set_ok(io_err)
    field = ''
    field_pos = 1
    end_of_line = .false.
    end_of_file = .false.
    in_complex = .false.
    paren_count = 0
    
    do
      read(file_unit, POS=stream_pos, iostat=io_err) char_read
      if (is_err(io_err)) then
        end_of_file = .true.
        return
      end if
      
      ! Handle line endings
      if (char_read == char(10)) then ! LF
        end_of_line = .true.
        return
      else if (char_read == char(13)) then ! CR
        ! Check for CR-LF combination
        stream_pos = stream_pos + 1
        read(file_unit, POS=stream_pos, iostat=io_err) char_read
        if (char_read == char(10)) then
          ! If it is CR-LF, just continue
          stream_pos = stream_pos + 1
        else if (is_err(io_err)) then
          call set_ok(io_err)
        end if
        stream_pos = stream_pos - 1
        end_of_line = .true.
        return
      end if
      
      ! Track parentheses for complex numbers
      if (char_read == '(') then
        in_complex = .true.
        paren_count = paren_count + 1
      else if (char_read == ')') then
        paren_count = paren_count - 1
        if (paren_count <= 0) then
          in_complex = .false.
          paren_count = 0
        end if
      end if
      
      ! Check for separator (only if not inside parentheses)
      if (.not. in_complex .and. char_read == sep) then
        stream_pos = stream_pos + 1
        return
      end if
      
      ! Add character to field
      if (field_pos <= len(field)) then
        field(field_pos:field_pos) = char_read
        field_pos = field_pos + 1
        stream_pos = stream_pos + 1
      else
        ! Field buffer overflow
        call set_err_once(io_err, ERR_INVALID_INPUT)
        return
      end if
    end do
  end subroutine read_field

end module tox_csv_file_reader