!> Module for reading heterogeneous CSV tables into type-banded arrays and providing fast accessors.
module tox_csv_file_reader
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use tox_errors, only: ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_FILE_EMPTY, set_ok, set_err_once, is_ok

contains

  !> Reads a table from a CSV file into type-banded arrays
  subroutine read_table(file_path, column_types, has_header, int_cols, real_cols, char_cols, &
                             logical_cols, complex_cols, header, metadata, ierr, sep_char, column_names)
    implicit none
    !-------------------------------------------------------------
    ! Reads a heterogeneous table from a CSV file and fills type-banded arrays.
    ! Args:
    !   file_path      : Path to the CSV file (input)
    !   column_types   : 1D int array, type for each column (1=int, 2=real, 3=char, 4=logical, 5=complex) (input)
    !   has_header     : logical, true if first line contains column names (input)
    !   int_cols       : 2D int array, output for integer columns
    !   real_cols      : 2D real array, output for real columns
    !   char_cols      : 2D char array, output for character columns
    !   logical_cols   : 2D logical array, output for logical columns
    !   complex_cols   : 2D complex array, output for complex columns
    !   header         : 1D char array, output column names
    !   metadata       : 2D int array, output metadata (row 1: type, row 2: index in type array)
    !   ierr           : int, error code (output)
    !   sep_char       : optional character, column separator (input, default ',')
    !   column_names   : optional 1D char array, overrides header line if provided (input)
    !-------------------------------------------------------------
    character(len=*), intent(in) :: file_path
    integer(int32), intent(in) :: column_types(:)
    logical, intent(in) :: has_header
    integer(int32), intent(out) :: int_cols(:,:)
    real(real64), intent(out) :: real_cols(:,:)
    character(len=*), intent(out) :: char_cols(:,:)
    logical, intent(out) :: logical_cols(:,:)
    complex(real64), intent(out) :: complex_cols(:,:)
    character(len=*), intent(out) :: header(:)
    integer(int32), intent(out) :: metadata(:,:)
    integer(int32), intent(out) :: ierr
    character(len=*), intent(in), optional :: sep_char
    character(len=*), intent(in), optional :: column_names(:) !//TODO clarify if both should be optional and in what order

    ! Local variables
    integer(int32) :: current_row, current_column, n_columns, n_rows, col_type, type_index, file_unit, data_offset, io_err, i_col_int, i_col_real, i_col_char, i_col_logical, i_col_complex, current_pos, sep_pos, k, l
    character(len=1024) :: line, field
    character(len=1) :: sep

    !//TODO Error handling for not equal lengths of column_types, column_names, header, metadata dimensions

    ! Step 1: Set separator (optional argument)
    if (present(sep_char)) then
      sep = sep_char(1:1)
    else
      sep = ','
    end if

    ! Step 2: Open file
    open(unit=file_unit, file=file_path, status='old', action='read', iostat=io_err)
    if (io_err /= 0) then
      print *, 'Error opening file: ', file_path !//TODO better error handling
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

    ! Step 3: Read header or use provided column_names
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
      if (io_err /= 0) exit !//TODO better error handling
    end do

    ! Step 4: Count number of rows
    n_rows = 0
    do
      read(file_unit, '(A)', iostat=io_err) line
      if (io_err /= 0) exit
      n_rows = n_rows + 1
    end do
    rewind(file_unit)
    if (has_header .and. .not. present(column_names)) then
      read(file_unit, '(A)', iostat=io_err) line ! skip header
    end if

    !//TODO clarify which structure complex numbers should have
    ! Step 5: Read data rows and fill type-banded arrays
    metadata = 0
    n_columns = size(column_types)
    do current_row = 1, n_rows
      ! Read the next line
      read(file_unit, '(A)', iostat=io_err) line
      if (io_err /= 0) exit !//TODO better error handling
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
        sep_pos = index(line(current_pos:), sep)

        ! If no separator found, take rest of line
        if (sep_pos == 0) then
          field = adjustl(trim(line(current_pos:)))
        else
          ! Extract field using current position and separator position
          field = adjustl(trim(line(current_pos:current_pos+sep_pos-2)))
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
          end if
          ! Increment column index for next int column
          i_col_int = i_col_int + 1
        
        case (2) ! real
          read(field, *, iostat=io_err) real_cols(current_row, i_col_real)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 2
            metadata(2, current_column) = i_col_real
          end if
          ! Increment column index for next real column
          i_col_real = i_col_real + 1

        case (3) ! char
          read(field, *, iostat=io_err) char_cols(current_row, i_col_char)
          ! If first time this column type is encountered, set metadata
          if (metadata(1, current_column) == 0) then
            metadata(1, current_column) = 3
            metadata(2, current_column) = i_col_char
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
          end if
          ! Increment column index for next logical column
          i_col_logical = i_col_logical + 1

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

  !> Getter functions for column access by name or index
  !pure function get_real_column(name_or_index, real_cols, header, metadata) result(col)
   ! implicit none
    ! ...existing code...
  !end function get_real_column

  !pure function get_int_column(name_or_index, int_cols, header, metadata) result(col)
   ! implicit none
    ! ...existing code...
  !end function get_int_column

  !pure function get_char_column(name_or_index, char_cols, header, metadata) result(col)
   ! implicit none
    ! ...existing code...
  !end function get_char_column

  !pure function get_logical_column(name_or_index, logical_cols, header, metadata) result(col)
   ! implicit none
    ! ...existing code...
!  end function get_logical_column

  !pure function get_complex_column(name_or_index, complex_cols, header, metadata) result(col)
  !  implicit none
    ! ...existing code...
!  end function get_complex_column

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
