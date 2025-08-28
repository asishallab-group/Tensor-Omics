!> Module to convert specified columns from a 2D string array into a 2D real array.
MODULE csv_read_real_module
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    USE tox_errors, ONLY: ERR_OK, ERR_EMPTY_INPUT, ERR_INVALID_INPUT, ERR_INVALID_FORMAT, set_ok, set_err_once
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_real_columns

CONTAINS

    !> Reads specified columns from a string array and converts them to reals.
    !| Non-convertible fields are set to -9999.0 and an error status is flagged.
    SUBROUTINE read_real_columns(data_in, cols_to_read, real_data_out, status)
        !| Input 2D array of strings.
        CHARACTER(LEN=MAX_FIELD_LEN), INTENT(IN) :: data_in(:,:)
        !| 1D array of 1-based column indices to convert.
        INTEGER(INT32), INTENT(IN) :: cols_to_read(:)
        !| Output 2D array of double precision reals.
        REAL(REAL64), ALLOCATABLE, INTENT(OUT) :: real_data_out(:,:)
        !| Output status code.
        INTEGER(INT32), INTENT(OUT) :: status

        INTEGER(INT32) :: num_rows, num_cols_in, num_cols_out
        INTEGER(INT32) :: i, j, col_idx
        INTEGER(INT32) :: read_stat
        CHARACTER(LEN=MAX_FIELD_LEN) :: temp_field

        CALL set_ok(status)
        num_rows = SIZE(data_in, DIM=1)
        num_cols_in = SIZE(data_in, DIM=2)
        num_cols_out = SIZE(cols_to_read)

        IF (num_rows == 0 .OR. num_cols_out == 0) THEN
            status = ERR_EMPTY_INPUT
            RETURN
        END IF

        ALLOCATE(real_data_out(num_rows, num_cols_out))

        DO j = 1, num_cols_out
            col_idx = cols_to_read(j)
            IF (col_idx < 1 .OR. col_idx > num_cols_in) THEN
                CALL set_err_once(status, ERR_INVALID_INPUT)
                real_data_out(:, j) = -9999.0_REAL64
                CYCLE
            END IF

            DO i = 1, num_rows
                temp_field = TRIM(data_in(i, col_idx))
                READ(temp_field, *, IOSTAT=read_stat) real_data_out(i, j)
                IF (read_stat /= 0) THEN
                    real_data_out(i, j) = -9999.0_REAL64
                    CALL set_err_once(status, ERR_INVALID_FORMAT)
                END IF
            END DO
        END DO

    END SUBROUTINE read_real_columns

END MODULE csv_read_real_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> C interface for reading real columns from a flat character array.
!| !! When using this C wrapper function, no copies of the arrays will be created. The Fortran routine will operate directly on the memory provided by the caller.
SUBROUTINE read_real_columns_c(data_in_ascii, num_rows, num_cols_in, &
                               cols_to_read, num_cols_to_read, &
                               real_data_out, status) &
                               BIND(C, NAME='read_real_columns_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_read_real_module, ONLY: read_real_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input flat array of ASCII codes representing the 2D string data.
    INTEGER(C_INT), INTENT(IN) :: data_in_ascii(*)
    !| Number of rows and columns in the input data.
    INTEGER(C_INT), VALUE, INTENT(IN) :: num_rows, num_cols_in
    !| Number of columns to read.
    INTEGER(C_INT), VALUE, INTENT(IN) :: num_cols_to_read
    !| 1D array of 1-based column indices to convert.
    INTEGER(C_INT), INTENT(IN) :: cols_to_read(*)
    !| Output buffer for the converted reals.
    REAL(C_DOUBLE), INTENT(OUT) :: real_data_out(*)
    !| Output status code.
    INTEGER(C_INT), INTENT(OUT) :: status

    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in_fortran(:,:)
    INTEGER(C_INT), ALLOCATABLE :: cols_to_read_fortran(:)
    REAL(C_DOUBLE), ALLOCATABLE :: real_data_out_fortran(:,:)
    INTEGER(C_INT) :: fortran_status
    INTEGER :: i, j, k, char_idx

    ALLOCATE(data_in_fortran(num_rows, num_cols_in))
    data_in_fortran = ' '
    DO j = 1, num_cols_in
        DO i = 1, num_rows
            DO k = 1, MAX_FIELD_LEN
                char_idx = (i - 1) * num_cols_in * MAX_FIELD_LEN + (j - 1) * MAX_FIELD_LEN + k
                IF (data_in_ascii(char_idx) == 0) EXIT
                data_in_fortran(i, j)(k:k) = CHAR(data_in_ascii(char_idx))
            END DO
        END DO
    END DO

    ALLOCATE(cols_to_read_fortran(num_cols_to_read))
    cols_to_read_fortran = cols_to_read(1:num_cols_to_read)

    CALL read_real_columns(data_in_fortran, cols_to_read_fortran, real_data_out_fortran, fortran_status)
    status = fortran_status

    IF (ALLOCATED(real_data_out_fortran)) THEN
        DO j = 1, num_cols_to_read
            DO i = 1, num_rows
                ! Corrected formula for column-major flattening
                char_idx = (j - 1) * num_rows + i
                real_data_out(char_idx) = real_data_out_fortran(i, j)
            END DO
        END DO
    END IF

    DEALLOCATE(data_in_fortran, cols_to_read_fortran)
    IF (ALLOCATED(real_data_out_fortran)) DEALLOCATE(real_data_out_fortran)

END SUBROUTINE read_real_columns_c

!> R interface for reading real columns.
!| !! When using this R wrapper function, copies of the arrays will be created. No direct modification of the original R objects occurs.
SUBROUTINE read_real_columns_r(data_in_ascii, num_rows, num_cols_in, &
                               cols_to_read, num_cols_to_read, &
                               real_data_out, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    USE csv_read_real_module, ONLY: read_real_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    USE tox_errors, ONLY: set_ok
    IMPLICIT NONE

    !| Input 3D array of ASCII codes (MAX_FIELD_LEN, num_rows, num_cols_in).
    INTEGER(INT32), INTENT(IN) :: data_in_ascii(MAX_FIELD_LEN, num_rows, num_cols_in)
    !| Number of rows and columns in the input data.
    INTEGER(INT32), INTENT(IN) :: num_rows, num_cols_in
    !| Number of columns to read.
    INTEGER(INT32), INTENT(IN) :: num_cols_to_read
    !| 1D array of 1-based column indices to convert.
    INTEGER(INT32), INTENT(IN) :: cols_to_read(num_cols_to_read)
    !| Output buffer for the converted reals.
    REAL(REAL64), INTENT(OUT) :: real_data_out(num_rows, num_cols_to_read)
    !| Output status code.
    INTEGER(INT32), INTENT(OUT) :: status

    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in_fortran(:,:)
    REAL(REAL64), ALLOCATABLE :: real_data_out_fortran(:,:)
    INTEGER :: i, j, k

    CALL set_ok(status)

    ALLOCATE(data_in_fortran(num_rows, num_cols_in))
    DO j = 1, num_cols_in
        DO i = 1, num_rows
            data_in_fortran(i, j) = ' '
            DO k = 1, MAX_FIELD_LEN
                IF (data_in_ascii(k, i, j) == 0) EXIT
                data_in_fortran(i, j)(k:k) = CHAR(data_in_ascii(k, i, j))
            END DO
        END DO
    END DO

    CALL read_real_columns(data_in_fortran, cols_to_read, real_data_out_fortran, status)

    IF (ALLOCATED(real_data_out_fortran)) THEN
        real_data_out = real_data_out_fortran
    END IF

    DEALLOCATE(data_in_fortran)
    IF (ALLOCATED(real_data_out_fortran)) DEALLOCATE(real_data_out_fortran)

END SUBROUTINE read_real_columns_r