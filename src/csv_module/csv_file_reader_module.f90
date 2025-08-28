!> Module to read a CSV file into a 2D array of strings.
MODULE csv_file_reader_module
    USE, INTRINSIc :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_csv_to_strings

CONTAINS

    !> Reads a CSV file and returns its content as a 2D array of strings.
    !| Reads all non-empty lines in a first pass, allocates memory, and populates the output arrays in a second pass.
    SUBROUTINE read_csv_to_strings(filename, has_header, delimiter, header_out, data_out, status)
        !| Path to the input CSV file.
        CHARACTER(LEN=*), INTENT(IN) :: filename
        !| Flag indicating if the file has a header row.
        LOGICAL, INTENT(IN) :: has_header
        !| Single character used to separate fields.
        CHARACTER(LEN=1), INTENT(IN) :: delimiter
        !| Output array for the header row.
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, INTENT(OUT) :: header_out(:)
        !| Output 2D array for the data rows.
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, INTENT(OUT) :: data_out(:,:)
        !| Output status code (0 for success).
        INTEGER(INT32), INTENT(OUT) :: status

        INTEGER :: unit_num
        CHARACTER(LEN=MAX_FIELD_LEN*20) :: line
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: all_lines(:)
        INTEGER(INT32) :: line_count, i
        INTEGER(INT32) :: num_rows, num_cols, header_offset
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
        INTEGER(INT32) :: io_stat_check

        status = 0
        header_offset = 0

        OPEN(NEWUNIT=unit_num, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=io_stat_check)
        IF (io_stat_check /= 0) THEN
            status = 10
            RETURN
        END IF

        line_count = 0
        DO
            READ(unit_num, '(A)', IOSTAT=io_stat_check) line
            IF (io_stat_check /= 0) EXIT
            IF (LEN_TRIM(line) > 0) THEN
                line_count = line_count + 1
            END IF
        END DO

        IF (line_count == 0) THEN
            status = 4
            CLOSE(unit_num)
            RETURN
        END IF

        REWIND(unit_num)
        ALLOCATE(all_lines(line_count))

        i = 0
        DO
            READ(unit_num, '(A)', IOSTAT=io_stat_check) line
            IF (io_stat_check /= 0) EXIT
            IF (LEN_TRIM(line) > 0) THEN
                i = i + 1
                all_lines(i) = TRIM(line)
            END IF
        END DO
        CLOSE(unit_num)

        IF (has_header) header_offset = 1

        num_rows = line_count - header_offset
        IF (num_rows <= 0) THEN
            status = 5
            DEALLOCATE(all_lines)
            RETURN
        END IF

        CALL parse_line(all_lines(1 + header_offset), delimiter, fields, status)
        num_cols = SIZE(fields)
        DEALLOCATE(fields)

        ALLOCATE(data_out(num_rows, num_cols))
        IF (has_header) THEN
            CALL parse_line(all_lines(1), delimiter, header_out, status)
        END IF

        DO i = 1, num_rows
            CALL parse_line(all_lines(i + header_offset), delimiter, fields, status)
            IF (SIZE(fields) == num_cols) THEN
                data_out(i, :) = fields
            ELSE
                data_out(i, :) = ""
            END IF
            DEALLOCATE(fields)
        END DO

        DEALLOCATE(all_lines)
        status = 0
    END SUBROUTINE read_csv_to_strings

END MODULE csv_file_reader_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> C interface for getting CSV dimensions before reading.
SUBROUTINE get_csv_dims_c(filename_ascii, fn_len, has_header, delimiter_ascii, &
                          num_rows, num_cols, status) &
                          BIND(C, NAME='get_csv_dims_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input ASCII codes of the filename.
    INTEGER(C_INT), INTENT(IN) :: filename_ascii(*)
    !| Length of the filename.
    INTEGER(C_INT), VALUE, INTENT(IN) :: fn_len
    !| Flag indicating if the file has a header.
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: has_header
    !| ASCII code of the delimiter character.
    INTEGER(C_INT), VALUE, INTENT(IN) :: delimiter_ascii
    !| Output for number of data rows.
    INTEGER(C_INT), INTENT(OUT) :: num_rows
    !| Output for number of columns.
    INTEGER(C_INT), INTENT(OUT) :: num_cols
    !| Output status code.
    INTEGER(C_INT), INTENT(OUT) :: status

    CHARACTER(LEN=:), ALLOCATABLE :: filename
    CHARACTER(LEN=1) :: delimiter
    CHARACTER(LEN=4096) :: line
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
    INTEGER :: i, unit, f_stat, line_count, header_offset
    LOGICAL :: has_header_f

    has_header_f = has_header
    ALLOCATE(CHARACTER(LEN=fn_len) :: filename)
    DO i = 1, fn_len
        filename(i:i) = CHAR(filename_ascii(i))
    END DO
    delimiter = CHAR(delimiter_ascii)

    num_rows = 0; num_cols = 0; status = 0; line_count = 0; header_offset = 0

    OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=f_stat)
    IF (f_stat /= 0) THEN; status = 10; RETURN; END IF

    DO
        READ(unit, '(A)', IOSTAT=f_stat) line
        IF (f_stat /= 0) EXIT
        IF (LEN_TRIM(line) > 0) line_count = line_count + 1
    END DO
    REWIND(unit)

    IF (line_count == 0) THEN; status = 4; CLOSE(unit); RETURN; END IF
    IF (has_header_f) header_offset = 1
    num_rows = line_count - header_offset
    IF (num_rows <= 0) THEN; status = 5; CLOSE(unit); RETURN; END IF
    
    IF (has_header_f) READ(unit, '(A)')
    READ(unit, '(A)') line
    CLOSE(unit)

    CALL parse_line(line, delimiter, fields, f_stat)
    num_cols = SIZE(fields)
    DEALLOCATE(fields)
END SUBROUTINE get_csv_dims_c

!> R interface for getting CSV dimensions before reading.
SUBROUTINE get_csv_dims_r(filename_ascii, fn_len, has_header, delimiter_ascii, &
                          num_rows, num_cols, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input ASCII codes of the filename.
    INTEGER(INT32), INTENT(IN) :: filename_ascii(fn_len)
    !| Length of the filename.
    INTEGER(INT32), INTENT(IN) :: fn_len
    !| Flag indicating if the file has a header.
    LOGICAL, INTENT(IN) :: has_header
    !| ASCII code of the delimiter character.
    INTEGER(INT32), INTENT(IN) :: delimiter_ascii
    !| Output for number of data rows.
    INTEGER(INT32), INTENT(OUT) :: num_rows
    !| Output for number of columns.
    INTEGER(INT32), INTENT(OUT) :: num_cols
    !| Output status code.
    INTEGER(INT32), INTENT(OUT) :: status
    
    CHARACTER(LEN=:), ALLOCATABLE :: filename
    CHARACTER(LEN=1) :: delimiter
    CHARACTER(LEN=4096) :: line
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
    INTEGER :: i, unit, f_stat, line_count, header_offset

    ALLOCATE(CHARACTER(LEN=fn_len) :: filename)
    DO i = 1, fn_len
        filename(i:i) = CHAR(filename_ascii(i))
    END DO
    delimiter = CHAR(delimiter_ascii)

    num_rows = 0; num_cols = 0; status = 0; line_count = 0; header_offset = 0

    OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=f_stat)
    IF (f_stat /= 0) THEN; status = 10; RETURN; END IF

    DO
        READ(unit, '(A)', IOSTAT=f_stat) line
        IF (f_stat /= 0) EXIT
        IF (LEN_TRIM(line) > 0) line_count = line_count + 1
    END DO
    REWIND(unit)

    IF (line_count == 0) THEN; status = 4; CLOSE(unit); RETURN; END IF
    IF (has_header) header_offset = 1
    num_rows = line_count - header_offset
    IF (num_rows <= 0) THEN; status = 5; CLOSE(unit); RETURN; END IF
    
    IF (has_header) READ(unit, '(A)')
    READ(unit, '(A)') line
    CLOSE(unit)

    CALL parse_line(line, delimiter, fields, f_stat)
    num_cols = SIZE(fields)
    DEALLOCATE(fields)
END SUBROUTINE get_csv_dims_r

!> C interface for reading a CSV into a flat character array.
SUBROUTINE read_csv_to_strings_c(filename_ascii, fn_len, has_header, delimiter_ascii, &
                                 header_out_ascii, data_out_ascii, status) &
                                 BIND(C, NAME='read_csv_to_strings_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_file_reader_module, ONLY: read_csv_to_strings
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    
    IMPLICIT NONE
    
    !| Input ASCII codes of the filename.
    INTEGER(C_INT), INTENT(IN) :: filename_ascii(fn_len)
    !| Length of the filename.
    INTEGER(C_INT), VALUE, INTENT(IN) :: fn_len
    !| Flag indicating if the file has a header.
    LOGICAL(C_BOOL), VALUE, INTENT(IN) :: has_header
    !| ASCII code of the delimiter character.
    INTEGER(C_INT), VALUE, INTENT(IN) :: delimiter_ascii
    !| Output buffer for the header as a flat array of ASCII codes.
    INTEGER(C_INT), INTENT(OUT) :: header_out_ascii(*)
    !| Output buffer for the data as a flat array of ASCII codes.
    INTEGER(C_INT), INTENT(OUT) :: data_out_ascii(*)
    !| Output status code.
    INTEGER(C_INT), INTENT(OUT) :: status

    CHARACTER(LEN=:), ALLOCATABLE :: filename
    CHARACTER(LEN=1) :: delimiter
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: header_out_f(:), data_out_f(:,:)
    INTEGER(C_INT) :: f_stat
    INTEGER :: i, j, k, char_idx, num_rows, num_cols
    LOGICAL :: has_header_f

    has_header_f = has_header
    ALLOCATE(CHARACTER(LEN=fn_len) :: filename)
    DO i = 1, fn_len
        filename(i:i) = CHAR(filename_ascii(i))
    END DO
    delimiter = CHAR(delimiter_ascii)

    CALL read_csv_to_strings(filename, has_header_f, delimiter, header_out_f, data_out_f, f_stat)
    status = f_stat
    IF (status /= 0) RETURN

    num_rows = SIZE(data_out_f, DIM=1)
    num_cols = SIZE(data_out_f, DIM=2)

    IF (has_header_f) THEN
        DO i = 1, num_cols
            DO k = 1, LEN_TRIM(header_out_f(i))
                char_idx = (i-1) * MAX_FIELD_LEN + k
                header_out_ascii(char_idx) = ICHAR(header_out_f(i)(k:k))
            END DO
        END DO
    END IF

    DO j = 1, num_cols
        DO i = 1, num_rows
             DO k = 1, LEN_TRIM(data_out_f(i,j))
                char_idx = (i - 1) * num_cols * MAX_FIELD_LEN + (j - 1) * MAX_FIELD_LEN + k
                data_out_ascii(char_idx) = ICHAR(data_out_f(i,j)(k:k))
            END DO
        END DO
    END DO
END SUBROUTINE read_csv_to_strings_c

!> R interface for reading a CSV into a 3D integer array of ASCII codes.
SUBROUTINE read_csv_to_strings_r(filename_ascii, fn_len, has_header, delimiter_ascii, &
                                 num_rows, num_cols, header_out_ascii, data_out_ascii, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_file_reader_module, ONLY: read_csv_to_strings
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE
    
    !| Input ASCII codes of the filename.
    INTEGER(INT32), INTENT(IN) :: filename_ascii(fn_len)
    !| Length of the filename.
    INTEGER(INT32), INTENT(IN) :: fn_len
    !| Number of data rows in the file.
    INTEGER(INT32), INTENT(IN) :: num_rows
    !| Number of columns in the file.
    INTEGER(INT32), INTENT(IN) :: num_cols
    !| Flag indicating if the file has a header.
    LOGICAL, INTENT(IN) :: has_header
    !| ASCII code of the delimiter character.
    INTEGER(INT32), INTENT(IN) :: delimiter_ascii
    !| Output buffer for the header as a 2D array of ASCII codes.
    INTEGER(INT32), INTENT(OUT) :: header_out_ascii(MAX_FIELD_LEN, num_cols)
    !| Output buffer for the data as a 3D array of ASCII codes.
    INTEGER(INT32), INTENT(OUT) :: data_out_ascii(MAX_FIELD_LEN, num_rows, num_cols)
    !| Output status code.
    INTEGER(INT32), INTENT(OUT) :: status

    CHARACTER(LEN=:), ALLOCATABLE :: filename
    CHARACTER(LEN=1) :: delimiter
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: header_out_f(:), data_out_f(:,:)
    INTEGER :: i, j, k

    ALLOCATE(CHARACTER(LEN=fn_len) :: filename)
    DO i = 1, fn_len
         filename(i:i) = CHAR(filename_ascii(i))
    END DO
    delimiter = CHAR(delimiter_ascii)

    CALL read_csv_to_strings(filename, has_header, delimiter, header_out_f, data_out_f, status)
    IF (status /= 0) RETURN
    
    header_out_ascii = 0
    data_out_ascii = 0

    IF (has_header) THEN
        DO i = 1, num_cols
            DO k = 1, LEN_TRIM(header_out_f(i))
                header_out_ascii(k, i) = ICHAR(header_out_f(i)(k:k))
            END DO
        END DO
    END IF

    DO j = 1, num_cols
        DO i = 1, num_rows
            DO k = 1, LEN_TRIM(data_out_f(i,j))
                data_out_ascii(k, i, j) = ICHAR(data_out_f(i,j)(k:k))
            END DO
         END DO
    END DO
END SUBROUTINE read_csv_to_strings_r