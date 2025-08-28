!> Module to extract specified character columns from a 2D string array.
MODULE csv_read_char_module
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    USE tox_errors, ONLY: ERR_OK, ERR_EMPTY_INPUT, ERR_INVALID_INPUT, set_ok
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_character_columns

CONTAINS

    !> Extracts specified columns from a 2D string array.
    SUBROUTINE read_character_columns(data_in, cols_to_read, char_data_out, status)
        !| Input 2D array of strings.
        CHARACTER(LEN=MAX_FIELD_LEN), INTENT(IN) :: data_in(:,:)
        !| 1D array of 1-based column indices to extract.
        INTEGER(INT32), INTENT(IN) :: cols_to_read(:)
        !| Output 2D array containing the extracted character columns.
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, INTENT(OUT) :: char_data_out(:,:)
        !| Output status code.
        INTEGER(INT32), INTENT(OUT) :: status

        INTEGER(INT32) :: num_rows, num_cols_in, num_cols_out
        INTEGER(INT32) :: j, col_idx

        CALL set_ok(status)
        num_rows = SIZE(data_in, DIM=1)
        num_cols_in = SIZE(data_in, DIM=2)
        num_cols_out = SIZE(cols_to_read)

        IF (num_rows == 0 .OR. num_cols_out == 0) THEN
            status = ERR_EMPTY_INPUT
            RETURN
        END IF

        ALLOCATE(char_data_out(num_rows, num_cols_out))

        DO j = 1, num_cols_out
            col_idx = cols_to_read(j)
            IF (col_idx < 1 .OR. col_idx > num_cols_in) THEN
                status = ERR_INVALID_INPUT ! Column index out of bounds
                char_data_out(:, j) = ""
                CYCLE
            END IF
            
            char_data_out(:, j) = data_in(:, col_idx)
        END DO

    END SUBROUTINE read_character_columns

END MODULE csv_read_char_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> C interface for reading character columns from a flat character array.
!| !! When using this C wrapper function, no copies of the arrays will be created. The Fortran routine will operate directly on the memory provided by the caller.
SUBROUTINE read_character_columns_c(data_in_ascii, num_rows, num_cols_in, &
                                    cols_to_read, num_cols_to_read, &
                                    char_data_out_ascii, status) &
                                    BIND(C, NAME='read_character_columns_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_read_char_module, ONLY: read_character_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input flat array of ASCII codes representing the 2D string data.
    INTEGER(C_INT), INTENT(IN), TARGET :: data_in_ascii(*)
    !| Number of rows in the input data.
    INTEGER(C_INT), VALUE, INTENT(IN) :: num_rows, num_cols_in
    !| Number of columns to read.
    INTEGER(C_INT), VALUE, INTENT(IN) :: num_cols_to_read
    !| 1D array of 1-based column indices to extract.
    INTEGER(C_INT), INTENT(IN) :: cols_to_read(*)
    !| Output buffer for the extracted columns as a flat array of ASCII codes.
    INTEGER(C_INT), INTENT(OUT) :: char_data_out_ascii(*)
    !| Output status code.
    INTEGER(C_INT), INTENT(OUT) :: status

    CHARACTER(LEN=MAX_FIELD_LEN), POINTER :: data_in_fortran(:,:)
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: char_data_out_fortran(:,:)
    INTEGER(C_INT), ALLOCATABLE :: cols_to_read_fortran(:)
    INTEGER(C_INT) :: fortran_status
    INTEGER :: i, j, k, char_idx, total_out_size

    ! Create a zero-copy Fortran pointer to the C data
    CALL c_f_pointer(c_loc(data_in_ascii(1)), data_in_fortran, &
                     SHAPE=[num_rows, num_cols_in])

    ALLOCATE(cols_to_read_fortran(num_cols_to_read))
    cols_to_read_fortran = cols_to_read(1:num_cols_to_read)

    CALL read_character_columns(data_in_fortran, cols_to_read_fortran, char_data_out_fortran, fortran_status)
    status = fortran_status

    IF (ALLOCATED(char_data_out_fortran)) THEN
        total_out_size = num_rows * num_cols_to_read * MAX_FIELD_LEN
        DO i = 1, total_out_size
            char_data_out_ascii(i) = 0
        END DO

        DO j = 1, num_cols_to_read
            DO i = 1, num_rows
                 DO k = 1, LEN_TRIM(char_data_out_fortran(i, j))
                    char_idx = (i - 1) * num_cols_to_read * MAX_FIELD_LEN + (j - 1) * MAX_FIELD_LEN + k
                    char_data_out_ascii(char_idx) = ICHAR(char_data_out_fortran(i, j)(k:k))
                END DO
            END DO
        END DO
    END IF

    DEALLOCATE(cols_to_read_fortran)
    IF (ALLOCATED(char_data_out_fortran)) DEALLOCATE(char_data_out_fortran)

END SUBROUTINE read_character_columns_c

!> R interface for reading character columns.
!| !! When using this R wrapper function, copies of the arrays will be created. No direct modification of the original R objects occurs.
SUBROUTINE read_character_columns_r(data_in_ascii, num_rows, num_cols_in, &
                                    cols_to_read, num_cols_to_read, &
                                    char_data_out_ascii, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_read_char_module, ONLY: read_character_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE

    !| Input 3D array of ASCII codes (MAX_FIELD_LEN, num_rows, num_cols_in).
    INTEGER(INT32), INTENT(IN) :: data_in_ascii(MAX_FIELD_LEN, num_rows, num_cols_in)
    !| Number of rows and columns in the input data.
    INTEGER(INT32), INTENT(IN) :: num_rows, num_cols_in
    !| Number of columns to read.
    INTEGER(INT32), INTENT(IN) :: num_cols_to_read
    !| 1D array of 1-based column indices to extract.
    INTEGER(INT32), INTENT(IN) :: cols_to_read(num_cols_to_read)
    !| Output buffer for extracted columns (MAX_FIELD_LEN, num_rows, num_cols_to_read).
    INTEGER(INT32), INTENT(OUT) :: char_data_out_ascii(MAX_FIELD_LEN, num_rows, num_cols_to_read)
    !| Output status code.
    INTEGER(INT32), INTENT(OUT) :: status

    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in_fortran(:,:), char_data_out_fortran(:,:)
    INTEGER :: i, j, k

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

    CALL read_character_columns(data_in_fortran, cols_to_read, char_data_out_fortran, status)

    IF (ALLOCATED(char_data_out_fortran)) THEN
        char_data_out_ascii = 0
        DO j = 1, num_cols_to_read
            DO i = 1, num_rows
                DO k = 1, LEN_TRIM(char_data_out_fortran(i, j))
                    char_data_out_ascii(k, i, j) = ICHAR(char_data_out_fortran(i, j)(k:k))
                END DO
            END DO
        END DO
    END IF

    DEALLOCATE(data_in_fortran)
    IF (ALLOCATED(char_data_out_fortran)) DEALLOCATE(char_data_out_fortran)

END SUBROUTINE read_character_columns_r