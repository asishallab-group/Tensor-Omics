 !> Module to read a CSV file into a 2D array of strings.
module csv_file_reader_module
  use, intrinsic :: iso_fortran_env, only: int32
  use csv_parser_module, only: parse_line, MAX_FIELD_LEN
  use tox_errors, only: set_ok, set_err_once, is_ok, ERR_FILE_OPEN, ERR_EMPTY_INPUT
  implicit none

contains
  !> Reads a CSV file and returns its content as a 2D array of strings.
  !| Reads all non-empty lines in a first pass, allocates memory, and populates the output arrays in a second pass.
  subroutine read_csv_to_strings(filename, has_header, delimiter, header_out, data_out, ierr)
    !| Path to the input CSV file.
    character(len=*), intent(in) :: filename
    !| Flag indicating if the file has a header row.
    logical, intent(in) :: has_header
    !| Single character used to separate fields.
    character(len=1), intent(in) :: delimiter
    !| Output array for the header row.
    character(len=MAX_FIELD_LEN), allocatable, intent(out) :: header_out(:)
    !| Output 2D array for the data rows.
    character(len=MAX_FIELD_LEN), allocatable, intent(out) :: data_out(:, :)
    !| Output error code (0 for success).
    integer(int32), intent(out) :: ierr

    ! Local variables
    integer :: file_unit
    character(len=MAX_FIELD_LEN * 20) :: line
    integer(int32) :: i, data_offset
    integer(int32) :: num_rows, num_cols
    character(len=MAX_FIELD_LEN), allocatable :: fields(:)
    integer(int32) :: io_error

    ! Initialize error code.
    call set_ok(ierr)

    data_offset = 0

    ! Open file and check for error.
    open (newunit=file_unit, file=trim(filename), status='old', action='read', iostat=io_error)
    if (.not. is_ok(io_error)) then
      call set_err_once(ierr, ERR_FILE_OPEN)
      close (file_unit)
      return
    end if

    ! Check for empty file.
    read (file_unit, '(A)', iostat=io_error) line
    if (.not. is_ok(io_error)) then
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      close (file_unit)
      return
    end if

    ! Go back to the beginning of the file.
    rewind(file_unit)

    ! Skip initial comment lines starting with '#'.
    do
      read (file_unit, '(A)', iostat=io_error) line
      if (.not. is_ok(io_error)) exit
      data_offset = data_offset + 1
      if (len_trim(line) > 0) then
        if (line(1:1) /= '#') exit
      end if
    end do

    ! Check for empty lines after comments and before header / data.
    if (len_trim(line) == 0) then
      do
        read (file_unit, '(A)', iostat=io_error) line
        if ((.not. is_ok(io_error)) .or. len_trim(line) /= 0) exit
      end do
    end if

    ! Check for header and infer number of columns (either through header or first data row fields).
    if (has_header) then
      call parse_line(line, delimiter, header_out, ierr)
      data_offset = data_offset + 1
      num_cols = size(header_out)
    else
      call parse_line(line, delimiter, fields, ierr)
      num_cols = size(fields)
      backspace(file_unit)
    end if

    ! Count number of rows for allocation of data_out.
    num_rows = 0
    do
      read (file_unit, '(A)', iostat=io_error) line
      if (.not. is_ok(io_error)) exit
      if (len_trim(line) > 0) num_rows = num_rows + 1
    end do

    allocate(data_out(num_rows, num_cols))

    ! Go back in the file to the first line of data.
    rewind(file_unit)
    do i = 1, data_offset
      read (file_unit, '(A)', iostat=io_error)
    end do

    ! Read all non-empty lines into data_out.
    do i = 1, num_rows
      read (file_unit, '(A)', iostat=io_error) line
      if (len_trim(line) > 0) then
        call parse_line(line, delimiter, fields, ierr)
        if (size(fields) == num_cols) then
          data_out(i, :) = fields
        else
          data_out(i, :) = ""
        end if
      end if
      deallocate (fields)
    end do
  end subroutine read_csv_to_strings

end module csv_file_reader_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> C interface for getting CSV dimensions before reading.
subroutine get_csv_dims_c(filename_ascii, fn_len, has_header, delimiter_ascii, &
                          num_rows, num_cols, ierr) &
  bind(c, name='get_csv_dims_c')
  use, intrinsic :: iso_c_binding
  use csv_parser_module, only: parse_line, MAX_FIELD_LEN
  implicit none

  !| Input ASCII codes of the filename.
  integer(c_int), intent(in) :: filename_ascii(*)
  !| Length of the filename.
  integer(c_int), value, intent(in) :: fn_len
  !| Flag indicating if the file has a header.
  logical(c_bool), value, intent(in) :: has_header
  !| ASCII code of the delimiter character.
  integer(c_int), value, intent(in) :: delimiter_ascii
  !| Output for number of data rows.
  integer(c_int), intent(out) :: num_rows
  !| Output for number of columns.
  integer(c_int), intent(out) :: num_cols
  !| Output error code.
  integer(c_int), intent(out) :: ierr

  character(len=:), allocatable :: filename
  character(len=1) :: delimiter
  character(len=4096) :: line
  character(len=max_field_len), allocatable :: fields(:)
  integer :: i, unit, f_stat, line_count, header_offset
  logical :: has_header_f

  has_header_f = has_header
  allocate (character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do
  delimiter = char(delimiter_ascii)

  num_rows = 0; num_cols = 0; ierr = 0; line_count = 0; header_offset = 0

  open (newunit=unit, file=trim(filename), status='old', action='read', iostat=f_stat)
  if (f_stat /= 0) then; ierr = 10; return; end if

  do
    read (unit, '(A)', iostat=f_stat) line
    if (f_stat /= 0) exit
    if (len_trim(line) > 0) line_count = line_count + 1
  end do
  rewind (unit)

  if (line_count == 0) then; ierr = 4; close (unit); return; end if
  if (has_header_f) header_offset = 1
  num_rows = line_count - header_offset
  if (num_rows <= 0) then; ierr = 5; close (unit); return; end if

  if (has_header_f) read (unit, '(A)')
  read (unit, '(A)') line
  close (unit)

  call parse_line(line, delimiter, fields, f_stat)
  num_cols = size(fields)
  deallocate (fields)
end subroutine get_csv_dims_c

!> R interface for getting CSV dimensions before reading.
subroutine get_csv_dims_r(filename_ascii, fn_len, has_header, delimiter_ascii, &
                          num_rows, num_cols, ierr)
  use, intrinsic :: iso_fortran_env, only: int32
  use csv_parser_module, only: parse_line, MAX_FIELD_LEN
  implicit none

  !| Input ASCII codes of the filename.
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !| Length of the filename.
  integer(int32), intent(in) :: fn_len
  !| Flag indicating if the file has a header.
  logical, intent(in) :: has_header
  !| ASCII code of the delimiter character.
  integer(int32), intent(in) :: delimiter_ascii
  !| Output for number of data rows.
  integer(int32), intent(out) :: num_rows
  !| Output for number of columns.
  integer(int32), intent(out) :: num_cols
  !| Output error code.
  integer(int32), intent(out) :: ierr

  character(len=:), allocatable :: filename
  character(len=1) :: delimiter
  character(len=4096) :: line
  character(len=MAX_FIELD_LEN), allocatable :: fields(:)
  integer :: i, unit, f_stat, line_count, header_offset

  allocate (character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do
  delimiter = char(delimiter_ascii)

  num_rows = 0; num_cols = 0; ierr = 0; line_count = 0; header_offset = 0

  open (newunit=unit, file=trim(filename), status='old', action='read', iostat=f_stat)
  if (f_stat /= 0) then; ierr = 10; return; end if

  do
    read (unit, '(A)', iostat=f_stat) line
    if (f_stat /= 0) exit
    if (len_trim(line) > 0) line_count = line_count + 1
  end do
  rewind (unit)

  if (line_count == 0) then; ierr = 4; close (unit); return; end if
  if (has_header) header_offset = 1
  num_rows = line_count - header_offset
  if (num_rows <= 0) then; ierr = 5; close (unit); return; end if

  if (has_header) read (unit, '(A)')
  read (unit, '(A)') line
  close (unit)

  call parse_line(line, delimiter, fields, f_stat)
  num_cols = size(fields)
  deallocate (fields)
end subroutine get_csv_dims_r

!> C interface for reading a CSV into a flat character array.
subroutine read_csv_to_strings_c(filename_ascii, fn_len, has_header, delimiter_ascii, &
                                 header_out_ascii, data_out_ascii, ierr) &
  bind(c, name='read_csv_to_strings_c')
  use, intrinsic :: iso_c_binding
  use csv_file_reader_module, only: read_csv_to_strings
  use csv_parser_module, only: MAX_FIELD_LEN

  implicit none

  !| Input ASCII codes of the filename.
  integer(c_int), intent(in) :: filename_ascii(fn_len)
  !| Length of the filename.
  integer(c_int), value, intent(in) :: fn_len
  !| Flag indicating if the file has a header.
  logical(c_bool), value, intent(in) :: has_header
  !| ASCII code of the delimiter character.
  integer(c_int), value, intent(in) :: delimiter_ascii
  !| Output buffer for the header as a flat array of ASCII codes.
  integer(c_int), intent(out) :: header_out_ascii(*)
  !| Output buffer for the data as a flat array of ASCII codes.
  integer(c_int), intent(out) :: data_out_ascii(*)
  !| Output error code.
  integer(c_int), intent(out) :: ierr

  character(len=:), allocatable :: filename
  character(len=1) :: delimiter
  character(len=MAX_FIELD_LEN), allocatable :: header_out_f(:), data_out_f(:, :)
  integer(c_int) :: f_stat
  integer :: i, j, k, char_idx, num_rows, num_cols
  logical :: has_header_f

  has_header_f = has_header
  allocate (character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do
  delimiter = char(delimiter_ascii)

  call read_csv_to_strings(filename, has_header_f, delimiter, header_out_f, data_out_f, f_stat)
  ierr = f_stat
  if (ierr /= 0) return

  num_rows = size(data_out_f, dim=1)
  num_cols = size(data_out_f, dim=2)

  if (has_header_f) then
    do i = 1, num_cols
      do k = 1, len_trim(header_out_f(i))
        char_idx = (i - 1)*max_field_len + k
        header_out_ascii(char_idx) = ichar(header_out_f(i) (k:k))
      end do
    end do
  end if

  do j = 1, num_cols
    do i = 1, num_rows
      do k = 1, len_trim(data_out_f(i, j))
        char_idx = (i - 1)*num_cols*MAX_FIELD_LEN + (j - 1)*MAX_FIELD_LEN + k
        data_out_ascii(char_idx) = ichar(data_out_f(i, j) (k:k))
      end do
    end do
  end do
end subroutine read_csv_to_strings_c

!> R interface for reading a CSV into a 3D integer array of ASCII codes.
subroutine read_csv_to_strings_r(filename_ascii, fn_len, has_header, delimiter_ascii, &
                                 num_rows, num_cols, header_out_ascii, data_out_ascii, ierr)
  use, intrinsic :: iso_fortran_env, only: int32
  use csv_file_reader_module, only: read_csv_to_strings
  use csv_parser_module, only: MAX_FIELD_LEN
  implicit none

  !| Input ASCII codes of the filename.
  integer(int32), intent(in) :: filename_ascii(fn_len)
  !| Length of the filename.
  integer(int32), intent(in) :: fn_len
  !| Number of data rows in the file.
  integer(int32), intent(in) :: num_rows
  !| Number of columns in the file.
  integer(int32), intent(in) :: num_cols
  !| Flag indicating if the file has a header.
  logical, intent(in) :: has_header
  !| ASCII code of the delimiter character.
  integer(int32), intent(in) :: delimiter_ascii
  !| Output buffer for the header as a 2D array of ASCII codes.
  integer(int32), intent(out) :: header_out_ascii(MAX_FIELD_LEN, num_cols)
  !| Output buffer for the data as a 3D array of ASCII codes.
  integer(int32), intent(out) :: data_out_ascii(MAX_FIELD_LEN, num_rows, num_cols)
  !| Output error code.
  integer(int32), intent(out) :: ierr

  character(len=:), allocatable :: filename
  character(len=1) :: delimiter
  character(len=max_field_len), allocatable :: header_out_f(:), data_out_f(:, :)
  integer :: i, j, k

  allocate (character(len=fn_len) :: filename)
  do i = 1, fn_len
    filename(i:i) = char(filename_ascii(i))
  end do
  delimiter = char(delimiter_ascii)

  call read_csv_to_strings(filename, has_header, delimiter, header_out_f, data_out_f, ierr)
  if (ierr /= 0) return

  header_out_ascii = 0
  data_out_ascii = 0

  if (has_header) then
    do i = 1, num_cols
      do k = 1, len_trim(header_out_f(i))
        header_out_ascii(k, i) = ichar(header_out_f(i) (k:k))
      end do
    end do
  end if

  do j = 1, num_cols
    do i = 1, num_rows
      do k = 1, len_trim(data_out_f(i, j))
        data_out_ascii(k, i, j) = ichar(data_out_f(i, j) (k:k))
      end do
    end do
  end do
end subroutine read_csv_to_strings_r
