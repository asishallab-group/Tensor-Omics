!> @brief A module to convert specified columns from a 2D string array into a 2D logical array.
MODULE csv_read_logical_module
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: read_logical_columns

CONTAINS

    !> @brief Reads specified columns from a string array and converts them to logicals.
    SUBROUTINE read_logical_columns(data_in, cols_to_read, logical_data_out, status)
        CHARACTER(LEN=MAX_FIELD_LEN), INTENT(IN) :: data_in(:,:)
        INTEGER(INT32), INTENT(IN) :: cols_to_read(:)
        LOGICAL, ALLOCATABLE, INTENT(OUT) :: logical_data_out(:,:)
        INTEGER(INT32), INTENT(OUT) :: status

        INTEGER(INT32) :: num_rows, num_cols_in, num_cols_out
        INTEGER(INT32) :: i, j, col_idx
        INTEGER(INT32) :: read_stat
        CHARACTER(LEN=MAX_FIELD_LEN) :: temp_field

        status = 0
        num_rows = SIZE(data_in, DIM=1)
        num_cols_in = SIZE(data_in, DIM=2)
        num_cols_out = SIZE(cols_to_read)

        IF (num_rows == 0 .OR. num_cols_out == 0) THEN
            status = 1 ! No data to process
            RETURN
        END IF

        ALLOCATE(logical_data_out(num_rows, num_cols_out))

        DO j = 1, num_cols_out
            col_idx = cols_to_read(j)
            IF (col_idx < 1 .OR. col_idx > num_cols_in) THEN
                status = 2 ! Column index out of bounds
                logical_data_out(:, j) = .FALSE.
                CYCLE
            END IF

            DO i = 1, num_rows
                temp_field = TRIM(data_in(i, col_idx))
                READ(temp_field, *, IOSTAT=read_stat) logical_data_out(i, j)
                IF (read_stat /= 0) THEN
                    logical_data_out(i, j) = .FALSE.
                    status = 3 ! At least one conversion error occurred
                END IF
            END DO
        END DO

    END SUBROUTINE read_logical_columns

END MODULE csv_read_logical_module

! =============================================================================
! C and R Wrapper Subroutines
! =============================================================================

!> @brief C interface for reading logical columns from a flat character array.
SUBROUTINE read_logical_columns_c(data_in_ascii, num_rows, num_cols_in, &
                                  cols_to_read, num_cols_to_read, &
                                  logical_data_out, status) &
                                  BIND(C, NAME='read_logical_columns_c')
    USE, INTRINSIC :: iso_c_binding
    USE csv_read_logical_module, ONLY: read_logical_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE

    ! --- C Interoperable Arguments ---
    INTEGER(C_INT), INTENT(IN) :: data_in_ascii(*)
    INTEGER(C_INT), VALUE, INTENT(IN) :: num_rows, num_cols_in, num_cols_to_read
    INTEGER(C_INT), INTENT(IN) :: cols_to_read(*)
    LOGICAL(C_BOOL), INTENT(OUT) :: logical_data_out(*)
    INTEGER(C_INT), INTENT(OUT) :: status

    ! --- Fortran Local Variables ---
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in_fortran(:,:)
    INTEGER(C_INT), ALLOCATABLE :: cols_to_read_fortran(:)
    LOGICAL, ALLOCATABLE :: logical_data_out_fortran(:,:)
    INTEGER(C_INT) :: fortran_status
    INTEGER :: i, j, k, char_idx

    ! --- 1. Reconstruct Fortran arrays from flat C arrays ---
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

    ! --- 2. Call the core Fortran subroutine ---
    CALL read_logical_columns(data_in_fortran, cols_to_read_fortran, logical_data_out_fortran, fortran_status)
    status = fortran_status

    ! --- 3. Copy the Fortran result back to the C output buffer ---
    IF (ALLOCATED(logical_data_out_fortran)) THEN
        DO j = 1, num_cols_to_read
            DO i = 1, num_rows
                char_idx = (i - 1) * num_cols_to_read + j
                logical_data_out(char_idx) = logical_data_out_fortran(i, j)
            END DO
        END DO
    END IF

    ! --- 4. Cleanup ---
    DEALLOCATE(data_in_fortran, cols_to_read_fortran)
    IF (ALLOCATED(logical_data_out_fortran)) DEALLOCATE(logical_data_out_fortran)

END SUBROUTINE read_logical_columns_c

!> @brief R interface for reading logical columns.
SUBROUTINE read_logical_columns_r(data_in_ascii, num_rows, num_cols_in, &
                                  cols_to_read, num_cols_to_read, &
                                  logical_data_out, status)
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_read_logical_module, ONLY: read_logical_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    IMPLICIT NONE

    ! --- R Interoperable Arguments ---
    INTEGER(INT32), INTENT(IN) :: data_in_ascii(MAX_FIELD_LEN, num_rows, num_cols_in)
    INTEGER(INT32), INTENT(IN) :: num_rows, num_cols_in, num_cols_to_read
    INTEGER(INT32), INTENT(IN) :: cols_to_read(num_cols_to_read)
    LOGICAL, INTENT(OUT) :: logical_data_out(num_rows, num_cols_to_read)
    INTEGER(INT32), INTENT(OUT) :: status

    ! --- Fortran Local Variables ---
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in_fortran(:,:)
    LOGICAL, ALLOCATABLE :: logical_data_out_fortran(:,:)
    INTEGER :: i, j, k

    ! --- 1. Convert ASCII integer array to Fortran character array ---
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

    ! --- 2. Call the core Fortran subroutine ---
    CALL read_logical_columns(data_in_fortran, cols_to_read, logical_data_out_fortran, status)

    ! --- 3. Copy the result to the R output buffer ---
    IF (ALLOCATED(logical_data_out_fortran)) THEN
        logical_data_out = logical_data_out_fortran
    END IF

    ! --- 4. Cleanup ---
    DEALLOCATE(data_in_fortran)
    IF (ALLOCATED(logical_data_out_fortran)) DEALLOCATE(logical_data_out_fortran)

END SUBROUTINE read_logical_columns_r
