!> Module for parsing lines from delimited text files.
MODULE csv_parser_module
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE tox_errors, ONLY: set_ok
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: parse_line, MAX_FIELD_LEN
    
    INTEGER, PARAMETER :: MAX_FIELD_LEN = 512

CONTAINS

    !> Parses a single line of text into an array of fields.
    !| Correctly handles delimiters inside parentheses, e.g., for complex numbers.
    SUBROUTINE parse_line(line, delim, fields, status)
        !| The input character string to be parsed.
        CHARACTER(LEN=*), INTENT(IN) :: line
        !| The single character delimiter to split the line by.
        CHARACTER(LEN=1), INTENT(IN) :: delim
        !| Output allocatable array of strings containing the parsed fields.
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, INTENT(OUT) :: fields(:)
        !| Output status code (0 for success).
        INTEGER(INT32), INTENT(OUT) :: status

        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: temp_fields(:)
        INTEGER(INT32) :: n_fields, i, start_pos
        INTEGER(INT32) :: paren_level
        CHARACTER(LEN=LEN(line)) :: trimmed_line

        CALL set_ok(status)
        trimmed_line = TRIM(line)
        
        n_fields = 1
        paren_level = 0
        DO i = 1, LEN_TRIM(trimmed_line)
            IF (trimmed_line(i:i) == '(') paren_level = paren_level + 1
            IF (trimmed_line(i:i) == ')') paren_level = paren_level - 1
            IF (trimmed_line(i:i) == delim .AND. paren_level == 0) THEN
                n_fields = n_fields + 1
            END IF
        END DO

        ALLOCATE(temp_fields(n_fields))
        
        start_pos = 1
        paren_level = 0
        n_fields = 1
        DO i = 1, LEN_TRIM(trimmed_line)
            IF (trimmed_line(i:i) == '(') paren_level = paren_level + 1
            IF (trimmed_line(i:i) == ')') paren_level = paren_level - 1
            IF (trimmed_line(i:i) == delim .AND. paren_level == 0) THEN
                temp_fields(n_fields) = TRIM(trimmed_line(start_pos:i-1))
                start_pos = i + 1
                n_fields = n_fields + 1
            END IF
        END DO
        temp_fields(n_fields) = TRIM(trimmed_line(start_pos:))

        IF (ALLOCATED(fields)) DEALLOCATE(fields)
        ALLOCATE(fields(n_fields))
        fields = temp_fields
        DEALLOCATE(temp_fields)
    END SUBROUTINE parse_line

END MODULE csv_parser_module

! =============================================================================
! C and R Wrapper Subroutines (External)
! =============================================================================

!> C interface for parsing a single line.
!| !! When using this C wrapper function, no copies of the arrays will be created. The Fortran routine will operate directly on the memory provided by the caller.
SUBROUTINE parse_line_c(line_ascii, line_len, delimiter_ascii, &
                        fields_out_ascii, num_fields_out, status) &
                        BIND(C, NAME='parse_line_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input ASCII codes of the line to parse.
    INTEGER(C_INT), INTENT(IN) :: line_ascii(*)
    !| Length of the input line.
    INTEGER(C_INT), VALUE, INTENT(IN) :: line_len
    !| ASCII code of the delimiter character.
    INTEGER(C_INT), VALUE, INTENT(IN) :: delimiter_ascii
    !| Output buffer for parsed fields as a flat array of ASCII codes.
    INTEGER(C_INT), INTENT(OUT) :: fields_out_ascii(*)
    !| Output number of fields found.
    INTEGER(C_INT), INTENT(OUT) :: num_fields_out
    !| Output status code.
    INTEGER(C_INT), INTENT(OUT) :: status

    CHARACTER(LEN=:), ALLOCATABLE :: line_f
    CHARACTER(LEN=1) :: delimiter_f
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields_f(:)
    INTEGER :: i, k, char_idx
    INTEGER(C_INT) :: total_out_size
    INTEGER(C_INT) :: fortran_status

    ALLOCATE(CHARACTER(LEN=line_len) :: line_f)
    DO i = 1, line_len
        line_f(i:i) = CHAR(line_ascii(i))
    END DO
    delimiter_f = CHAR(delimiter_ascii)

    CALL parse_line(line_f, delimiter_f, fields_f, fortran_status)
    status = fortran_status
    IF (status /= 0) THEN
        DEALLOCATE(line_f)
        RETURN
    END IF

    num_fields_out = SIZE(fields_f)
    total_out_size = num_fields_out * MAX_FIELD_LEN
    DO i = 1, total_out_size
        fields_out_ascii(i) = 0
    END DO
    
    DO i = 1, num_fields_out
        DO k = 1, LEN_TRIM(fields_f(i))
            char_idx = (i-1) * MAX_FIELD_LEN + k
            fields_out_ascii(char_idx) = ICHAR(fields_f(i)(k:k))
        END DO
    END DO

    DEALLOCATE(line_f, fields_f)
END SUBROUTINE parse_line_c

!> R interface for parsing a single line.
!| !! When using this R wrapper function, copies of the arrays will be created. No direct modification of the original R objects occurs.
SUBROUTINE parse_line_r(line_ascii, line_len, delimiter_ascii, &
                        fields_out_ascii, num_fields_out, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input ASCII codes of the line to parse.
    INTEGER(INT32), INTENT(IN) :: line_ascii(line_len)
    !| Length of the input line.
    INTEGER(INT32), INTENT(IN) :: line_len
    !| ASCII code of the delimiter character.
    INTEGER(INT32), INTENT(IN) :: delimiter_ascii
    !| Output buffer for parsed fields as a 2D array of ASCII codes.
    INTEGER(INT32), INTENT(OUT) :: fields_out_ascii(MAX_FIELD_LEN, *)
    !| Output number of fields found.
    INTEGER(INT32), INTENT(OUT) :: num_fields_out
    !| Output status code.
    INTEGER(INT32), INTENT(OUT) :: status

    CHARACTER(LEN=:), ALLOCATABLE :: line_f
    CHARACTER(LEN=1) :: delimiter_f
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields_f(:)
    INTEGER :: i, k
    INTEGER(INT32) :: fortran_status

    ALLOCATE(CHARACTER(LEN=line_len) :: line_f)
    DO i = 1, line_len
        line_f(i:i) = CHAR(line_ascii(i))
    END DO
    delimiter_f = CHAR(delimiter_ascii)

    CALL parse_line(line_f, delimiter_f, fields_f, fortran_status)
    status = fortran_status
    IF (status /= 0) THEN
        DEALLOCATE(line_f)
        RETURN
    END IF
    
    num_fields_out = SIZE(fields_f)
    DO i = 1, num_fields_out
        fields_out_ascii(:, i) = 0
    END DO
    DO i = 1, num_fields_out
         DO k = 1, LEN_TRIM(fields_f(i))
            fields_out_ascii(k, i) = ICHAR(fields_f(i)(k:k))
        END DO
    END DO

    DEALLOCATE(line_f, fields_f)
END SUBROUTINE parse_line_r