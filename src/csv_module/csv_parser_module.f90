!> @brief A module for parsing lines from a CSV file.
MODULE csv_parser_module
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: parse_line, MAX_FIELD_LEN
    
    INTEGER, PARAMETER :: MAX_FIELD_LEN = 512

CONTAINS

    !> @brief Parses a single line of text into an array of fields.
    SUBROUTINE parse_line(line, delim, fields, status)
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=1), INTENT(IN) :: delim
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, INTENT(OUT) :: fields(:)
        INTEGER(INT32), INTENT(OUT) :: status

        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: temp_fields(:)
        INTEGER(INT32) :: n_fields, i, start_pos
        INTEGER(INT32) :: paren_level
        CHARACTER(LEN=LEN(line)) :: trimmed_line

        status = 0
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

!> @brief C interface for parsing a single line.
SUBROUTINE parse_line_c(line_ascii, line_len, delimiter_ascii, &
                        fields_out_ascii, num_fields_out, status) &
                        BIND(C, NAME='parse_line_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    ! --- C Arguments ---
    INTEGER(C_INT), INTENT(IN) :: line_ascii(*)
    INTEGER(C_INT), VALUE, INTENT(IN) :: line_len, delimiter_ascii
    INTEGER(C_INT), INTENT(OUT) :: fields_out_ascii(*), num_fields_out, status

    ! --- Fortran Local Variables ---
    CHARACTER(LEN=:), ALLOCATABLE :: line_f
    CHARACTER(LEN=1) :: delimiter_f
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields_f(:)
    INTEGER :: i, k, char_idx
    INTEGER(C_INT) :: total_out_size
    INTEGER(C_INT) :: fortran_status

    ! --- 1. Convert C inputs to Fortran types ---
    ALLOCATE(CHARACTER(LEN=line_len) :: line_f)
    DO i = 1, line_len
        line_f(i:i) = CHAR(line_ascii(i))
    END DO
    delimiter_f = CHAR(delimiter_ascii)

    ! --- 2. Call the core Fortran subroutine ---
    CALL parse_line(line_f, delimiter_f, fields_f, fortran_status)
    status = fortran_status
    IF (status /= 0) THEN
        DEALLOCATE(line_f)
        RETURN
    END IF

    ! --- 3. Copy Fortran results back to C buffer ---
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

    ! --- 4. Cleanup ---
    DEALLOCATE(line_f, fields_f)
END SUBROUTINE parse_line_c

!> @brief R interface for parsing a single line.
SUBROUTINE parse_line_r(line_ascii, line_len, delimiter_ascii, &
                        fields_out_ascii, num_fields_out, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    IMPLICIT NONE

    ! --- R Arguments ---
    INTEGER(INT32), INTENT(IN) :: line_ascii(line_len)
    INTEGER(INT32), INTENT(IN) :: line_len, delimiter_ascii
    INTEGER(INT32), INTENT(OUT) :: fields_out_ascii(MAX_FIELD_LEN, *)
    INTEGER(INT32), INTENT(OUT) :: num_fields_out, status

    ! --- Fortran Local Variables ---
    CHARACTER(LEN=:), ALLOCATABLE :: line_f
    CHARACTER(LEN=1) :: delimiter_f
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields_f(:)
    INTEGER :: i, k
    INTEGER(INT32) :: fortran_status

    ! --- 1. Convert R inputs to Fortran types ---
    ALLOCATE(CHARACTER(LEN=line_len) :: line_f)
    DO i = 1, line_len
        line_f(i:i) = CHAR(line_ascii(i))
    END DO
    delimiter_f = CHAR(delimiter_ascii)

    ! --- 2. Call the core Fortran subroutine ---
    CALL parse_line(line_f, delimiter_f, fields_f, fortran_status)
    status = fortran_status
    IF (status /= 0) THEN
        DEALLOCATE(line_f)
        RETURN
    END IF
    
    ! --- 3. Copy Fortran results back to R buffer ---
    num_fields_out = SIZE(fields_f)
    DO i = 1, num_fields_out
        fields_out_ascii(:, i) = 0
    END DO
    DO i = 1, num_fields_out
         DO k = 1, LEN_TRIM(fields_f(i))
            fields_out_ascii(k, i) = ICHAR(fields_f(i)(k:k))
        END DO
    END DO

    ! --- 4. Cleanup ---
    DEALLOCATE(line_f, fields_f)
END SUBROUTINE parse_line_r