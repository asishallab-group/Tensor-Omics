MODULE csv_reader_interfaces
    IMPLICIT NONE
CONTAINS
    ! ===================================================================
    ! --- C Interface Functions (for R, Python, etc.)
    ! ===================================================================
    FUNCTION read_csv_file_C(filename_c, has_header_c, delimiter_c) RESULT(io_status_c) &
                               bind(C, name="read_csv_file_C")
        USE, INTRINSIC :: iso_c_binding
        USE csv_reader_module, ONLY: read_csv_file, MAX_LINE_LEN
        IMPLICIT NONE
        CHARACTER(KIND=c_char), DIMENSION(*), INTENT(IN) :: filename_c
        ! FIX: Added the VALUE attribute to scalar arguments
        LOGICAL(c_bool), VALUE, INTENT(IN) :: has_header_c
        CHARACTER(KIND=c_char), VALUE, INTENT(IN) :: delimiter_c
        INTEGER(c_int) :: io_status_c
        CHARACTER(LEN=MAX_LINE_LEN) :: filename_f
        LOGICAL :: has_header_f
        INTEGER :: i, j

        filename_f = ' '
        DO i = 1, MAX_LINE_LEN
            IF (filename_c(i) == C_NULL_CHAR) THEN
                DO j = 1, i - 1; filename_f(j:j) = filename_c(j); END DO
                EXIT
            END IF
        END DO
        IF (i > MAX_LINE_LEN) THEN
            DO j = 1, MAX_LINE_LEN; filename_f(j:j) = filename_c(j); END DO
        END IF

        has_header_f = has_header_c
        CALL read_csv_file(TRIM(filename_f), has_header_f, delimiter_c, io_status_c)
    END FUNCTION read_csv_file_C

    FUNCTION serialize_C(filename_c) RESULT(io_status_c) bind(C, name="serialize_C")
        USE, INTRINSIC :: iso_c_binding
        USE csv_reader_module, ONLY: serialize, MAX_LINE_LEN
        IMPLICIT NONE
        CHARACTER(KIND=c_char), DIMENSION(*), INTENT(IN) :: filename_c
        INTEGER(c_int) :: io_status_c
        CHARACTER(LEN=MAX_LINE_LEN) :: filename_f
        INTEGER :: i, j

        filename_f = ' '
        DO i = 1, MAX_LINE_LEN
            IF (filename_c(i) == C_NULL_CHAR) THEN
                DO j = 1, i - 1; filename_f(j:j) = filename_c(j); END DO
                EXIT
            END IF
        END DO
        IF (i > MAX_LINE_LEN) THEN
            DO j = 1, MAX_LINE_LEN; filename_f(j:j) = filename_c(j); END DO
        END IF
        CALL serialize(TRIM(filename_f), io_status_c)
    END FUNCTION serialize_C

    FUNCTION deserialize_C(filename_c) RESULT(io_status_c) bind(C, name="deserialize_C")
        USE, INTRINSIC :: iso_c_binding
        USE csv_reader_module, ONLY: deserialize, MAX_LINE_LEN
        IMPLICIT NONE
        CHARACTER(KIND=c_char), DIMENSION(*), INTENT(IN) :: filename_c
        INTEGER(c_int) :: io_status_c
        CHARACTER(LEN=MAX_LINE_LEN) :: filename_f
        INTEGER :: i, j

        filename_f = ' '
        DO i = 1, MAX_LINE_LEN
            IF (filename_c(i) == C_NULL_CHAR) THEN
                DO j = 1, i - 1; filename_f(j:j) = filename_c(j); END DO
                EXIT
            END IF
        END DO
        IF (i > MAX_LINE_LEN) THEN
            DO j = 1, MAX_LINE_LEN; filename_f(j:j) = filename_c(j); END DO
        END IF
        CALL deserialize(TRIM(filename_f), io_status_c)
    END FUNCTION deserialize_C

END MODULE csv_reader_interfaces