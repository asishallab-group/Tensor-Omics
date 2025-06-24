! =======================================================================
! INTERFACE WRAPPERS for csv_reader_module
! =======================================================================
! This file contains R and C language bindings for the procedures in
! the csv_reader_module. It should be compiled and linked with the
! original csv_reader_module.f90 file.
!
! The design follows the F42/Tensor Omics interfacing guidelines:
! 1. Core logic remains in the original Fortran module.
! 2. R wrappers are plain Fortran subroutines with an "_R" suffix.
! 3. C wrappers use iso_c_binding, have a "_C" suffix, and are kept
!    minimal to only handle data type and array convention mapping.
! =======================================================================

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! R INTERFACE WRAPPERS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! These subroutines are designed to be called from R via the
! .Fortran() function. They do not use bind(C).

SUBROUTINE read_csv_file_R(filename, has_header, delimiter, io_status)
    USE csv_reader_module
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL, INTENT(IN) :: has_header
    CHARACTER(LEN=1), INTENT(IN) :: delimiter
    INTEGER, INTENT(OUT) :: io_status
    CALL read_csv_file(TRIM(filename), has_header, delimiter, io_status)
END SUBROUTINE read_csv_file_R

SUBROUTINE serialize_R(filename, io_status)
    USE csv_reader_module
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: io_status
    CALL serialize(TRIM(filename), io_status)
END SUBROUTINE serialize_R

SUBROUTINE deserialize_R(filename, io_status)
    USE csv_reader_module
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: io_status
    CALL deserialize(TRIM(filename), io_status)
END SUBROUTINE deserialize_R

SUBROUTINE cleanup_csv_data_R()
    USE csv_reader_module
    IMPLICIT NONE
    CALL cleanup_csv_data()
END SUBROUTINE cleanup_csv_data_R

! R wrappers for functions are subroutines with an OUT parameter
SUBROUTINE get_num_rows_R(n)
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: n
    n = get_num_rows()
END SUBROUTINE get_num_rows_R

SUBROUTINE get_num_cols_R(n)
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER, INTENT(OUT) :: n
    n = get_num_cols()
END SUBROUTINE get_num_cols_R

SUBROUTINE get_header_R(j, header_name)
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: j
    CHARACTER(LEN=*), INTENT(OUT) :: header_name
    header_name = get_header(j)
END SUBROUTINE get_header_R

! For get_cell, the R wrapper will return all possible types.
! The R calling code will need to check the data_type output.
SUBROUTINE get_cell_R(i, j, data_type, i_val, r_val, c_val)
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i, j
    INTEGER, INTENT(OUT) :: data_type
    INTEGER(IK), INTENT(OUT) :: i_val
    REAL, INTENT(OUT) :: r_val
    CHARACTER(LEN=*), INTENT(OUT) :: c_val
    TYPE(generic_data_cell) :: cell

    cell = get_cell(i, j)
    data_type = cell%data_type
    i_val = 0_IK
    r_val = 0.0
    c_val = ''

    SELECT CASE (data_type)
    CASE (1)
        i_val = cell%i_val
    CASE (2)
        r_val = cell%r_val
    CASE (3)
        IF (ALLOCATED(cell%c_val)) THEN
            c_val = cell%c_val
            DEALLOCATE(cell%c_val)
        END IF
    END SELECT
END SUBROUTINE get_cell_R


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! C INTERFACE WRAPPERS
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! These subroutines use bind(C) and iso_c_binding to create a stable
! ABI for C, Python, and other languages. They handle the translation
! from C types and pointers to Fortran types and arrays.

SUBROUTINE read_csv_file_C(filename, n_chars, has_header, delimiter, io_status) &
    BIND(C, NAME='read_csv_file_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    CHARACTER(KIND=c_char), INTENT(IN) :: filename(*)
    INTEGER(c_int), VALUE :: n_chars
    LOGICAL(c_bool), VALUE :: has_header
    CHARACTER(KIND=c_char), VALUE :: delimiter
    INTEGER(c_int), INTENT(OUT) :: io_status

    CHARACTER(LEN=n_chars) :: f_filename
    LOGICAL :: f_has_header !<-- FIX: Declare a standard Fortran logical
    INTEGER :: f_io_status
    INTEGER :: i

    ! Convert C string to Fortran string
    DO i = 1, n_chars
        f_filename(i:i) = filename(i)
    END DO

    ! FIX: Assign C logical to Fortran logical to resolve kind mismatch
    f_has_header = has_header

    ! FIX: Call with the correctly-kinded logical
    CALL read_csv_file(TRIM(f_filename), f_has_header, delimiter, f_io_status)
    io_status = f_io_status
END SUBROUTINE read_csv_file_C

SUBROUTINE serialize_C(filename, n_chars, io_status) BIND(C, NAME='serialize_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    CHARACTER(KIND=c_char), INTENT(IN) :: filename(*)
    INTEGER(c_int), VALUE :: n_chars
    INTEGER(c_int), INTENT(OUT) :: io_status

    CHARACTER(LEN=n_chars) :: f_filename
    INTEGER :: f_io_status
    INTEGER :: i
    DO i = 1, n_chars
        f_filename(i:i) = filename(i)
    END DO
    CALL serialize(TRIM(f_filename), f_io_status)
    io_status = f_io_status
END SUBROUTINE serialize_C

SUBROUTINE deserialize_C(filename, n_chars, io_status) BIND(C, NAME='deserialize_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    CHARACTER(KIND=c_char), INTENT(IN) :: filename(*)
    INTEGER(c_int), VALUE :: n_chars
    INTEGER(c_int), INTENT(OUT) :: io_status

    CHARACTER(LEN=n_chars) :: f_filename
    INTEGER :: f_io_status
    INTEGER :: i
    DO i = 1, n_chars
        f_filename(i:i) = filename(i)
    END DO
    CALL deserialize(TRIM(f_filename), f_io_status)
    io_status = f_io_status
END SUBROUTINE deserialize_C

SUBROUTINE cleanup_csv_data_C() BIND(C, NAME='cleanup_csv_data_C')
    USE csv_reader_module
    IMPLICIT NONE
    CALL cleanup_csv_data()
END SUBROUTINE cleanup_csv_data_C

FUNCTION get_num_rows_C() RESULT(n) BIND(C, NAME='get_num_rows_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER(c_int) :: n
    n = get_num_rows()
END FUNCTION get_num_rows_C

FUNCTION get_num_cols_C() RESULT(n) BIND(C, NAME='get_num_cols_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER(c_int) :: n
    n = get_num_cols()
END FUNCTION get_num_cols_C

SUBROUTINE get_header_C(j, header_name_c, max_len) BIND(C, NAME='get_header_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER(c_int), VALUE :: j, max_len
    CHARACTER(KIND=c_char), INTENT(OUT) :: header_name_c(max_len)
    
    CHARACTER(LEN=MAX_FIELD_LEN) :: f_header
    INTEGER :: i, len_out

    f_header = get_header(j)
    len_out = MIN(LEN_TRIM(f_header), max_len - 1)
    
    ! Copy to C character array and null-terminate
    DO i = 1, len_out
        header_name_c(i) = f_header(i:i)
    END DO
    header_name_c(len_out + 1) = c_null_char
END SUBROUTINE get_header_C

SUBROUTINE get_cell_C(i, j, data_type, i_val, r_val, c_val_c, max_len) &
    BIND(C, NAME='get_cell_C')
    USE iso_c_binding
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER(c_int), VALUE :: i, j, max_len
    INTEGER(c_int), INTENT(OUT) :: data_type
    INTEGER(c_int64_t), INTENT(OUT) :: i_val  ! Match IK kind
    REAL(c_double), INTENT(OUT) :: r_val      ! Match REAL kind
    CHARACTER(KIND=c_char), INTENT(OUT) :: c_val_c(max_len)
    
    TYPE(generic_data_cell) :: cell
    INTEGER :: k, len_out

    ! Initialize output values
    data_type = 0
    i_val = 0
    r_val = 0.0D0
    c_val_c(1) = c_null_char

    cell = get_cell(i, j)
    data_type = cell%data_type

    SELECT CASE (data_type)
    CASE (1)
        i_val = cell%i_val
    CASE (2)
        r_val = cell%r_val
    CASE (3)
        IF (ALLOCATED(cell%c_val)) THEN
            len_out = MIN(LEN(cell%c_val), max_len - 1)
            DO k = 1, len_out
                c_val_c(k) = cell%c_val(k:k)
            END DO
            c_val_c(len_out + 1) = c_null_char
            DEALLOCATE(cell%c_val)
        END IF
    END SELECT
END SUBROUTINE get_cell_C
