!> Unit test suite for csv_read_char_module and its wrappers.
!| This file is already compatible with the new asserts.f90 module.
MODULE mod_test_csv_read_char
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE, INTRINSIC :: iso_c_binding
    USE csv_read_char_module, ONLY: read_character_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    USE asserts, ONLY: assert_equal_int, assert_string_equal

    IMPLICIT NONE
    PUBLIC

    INTERFACE
        SUBROUTINE read_character_columns_c(data_in_ascii, num_rows, num_cols_in, &
                                             cols_to_read, num_cols_to_read, &
                                             char_data_out_ascii, status) &
                                             BIND(C, NAME='read_character_columns_c')
            USE, INTRINSIC :: iso_c_binding
            INTEGER(C_INT), INTENT(IN), TARGET :: data_in_ascii(*)
            INTEGER(C_INT), VALUE, INTENT(IN) :: num_rows, num_cols_in, num_cols_to_read
            INTEGER(C_INT), INTENT(IN) :: cols_to_read(*)
            INTEGER(C_INT), INTENT(OUT) :: char_data_out_ascii(*)
            INTEGER(C_INT), INTENT(OUT) :: status
        END SUBROUTINE
    END INTERFACE

    ABSTRACT INTERFACE
        SUBROUTINE test_interface()
        END SUBROUTINE test_interface
    END INTERFACE

    TYPE :: test_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: test_proc => NULL()
    END TYPE test_case

CONTAINS

    FUNCTION get_all_tests_csv_read_char() RESULT(all_tests)
        TYPE(test_case) :: all_tests(2)
        all_tests(1) = test_case("test_read_char_cols_basic", test_read_char_cols_basic)
        all_tests(2) = test_case("test_read_char_cols_c_wrapper", test_read_char_cols_c_wrapper)
    END FUNCTION get_all_tests_csv_read_char

    SUBROUTINE run_all_tests_csv_read_char()
        TYPE(test_case) :: all_tests(2)
        INTEGER :: i
        all_tests = get_all_tests_csv_read_char()
        WRITE(*, '(A)') "--- Running Suite: csv_read_char ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: csv_read_char ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_csv_read_char

    SUBROUTINE run_named_tests_csv_read_char(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(2)
        INTEGER :: i, j
        LOGICAL :: found
        all_tests = get_all_tests_csv_read_char()
        WRITE(*, '(A)') "--- Running Suite: csv_read_char (named tests) ---"
        DO i = 1, SIZE(test_names)
            found = .FALSE.
            DO j = 1, SIZE(all_tests)
                IF (TRIM(test_names(i)) == TRIM(all_tests(j)%name)) THEN
                    WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(j)%name), "..."
                    CALL all_tests(j)%test_proc()
                    WRITE(*, '(A)') " PASSED"
                    found = .TRUE.
                    EXIT
                END IF
            END DO
            IF (.NOT. found) THEN
                WRITE(*,*) "WARNING: Unknown test '", TRIM(test_names(i)), "' in suite 'csv_read_char'"
            END IF
        END DO
        WRITE(*, '(A)') "--- Suite Finished: csv_read_char ---"
        WRITE(*,*)
    END SUBROUTINE run_named_tests_csv_read_char

    SUBROUTINE test_read_char_cols_basic()
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in(:,:), char_data_out(:,:)
        INTEGER(INT32) :: status
        INTEGER(INT32) :: cols_to_read(2)

        ALLOCATE(data_in(2, 3))
        data_in(1, 1) = "ID1"; data_in(1, 2) = "100"; data_in(1, 3) = "Val1"
        data_in(2, 1) = "ID2"; data_in(2, 2) = "200"; data_in(2, 3) = "Val2"
        cols_to_read = [1, 3]
        CALL read_character_columns(data_in, cols_to_read, char_data_out, status)

        CALL assert_equal_int(status, 0, "test_read_char_cols_basic: status")
        CALL assert_string_equal(char_data_out(1, 1), "ID1", "test_read_char_cols_basic: val(1,1)")
        CALL assert_string_equal(char_data_out(2, 2), "Val2", "test_read_char_cols_basic: val(2,2)")
        DEALLOCATE(data_in, char_data_out)
    END SUBROUTINE test_read_char_cols_basic

    SUBROUTINE test_read_char_cols_c_wrapper()
        CHARACTER(LEN=MAX_FIELD_LEN) :: str_arr(2,2), str_out(2,1)
        INTEGER(C_INT), ALLOCATABLE :: ascii_flat_in(:), ascii_flat_out(:)
        INTEGER(C_INT) :: status, i, j, k, char_idx
        INTEGER(C_INT) :: cols_to_read(1)
        INTEGER(C_INT), PARAMETER :: num_rows = 2, num_cols = 2, num_cols_out = 1

        ! Arrange
        str_arr(1,1) = "A1"; str_arr(1,2) = "B1"
        str_arr(2,1) = "A2"; str_arr(2,2) = "B2"
        ALLOCATE(ascii_flat_in(num_rows * num_cols * MAX_FIELD_LEN))
        ascii_flat_in = 0
        DO j = 1, num_cols
            DO i = 1, num_rows
                DO k = 1, LEN_TRIM(str_arr(i,j))
                    char_idx = (i - 1) * num_cols * MAX_FIELD_LEN + (j - 1) * MAX_FIELD_LEN + k
                    ascii_flat_in(char_idx) = ICHAR(str_arr(i,j)(k:k))
                END DO
            END DO
        END DO
        
        ALLOCATE(ascii_flat_out(num_rows * num_cols_out * MAX_FIELD_LEN))
        cols_to_read = [2]

        ! Act
        CALL read_character_columns_c(ascii_flat_in, num_rows, num_cols, cols_to_read, num_cols_out, ascii_flat_out, status)

        ! Assert: Reconstruct Fortran string from C output for comparison
        CALL assert_equal_int(status, 0, "c_wrapper: status")
        str_out = ' '
        DO i = 1, num_rows
            DO k = 1, MAX_FIELD_LEN
                char_idx = (i - 1) * num_cols_out * MAX_FIELD_LEN + k
                IF (ascii_flat_out(char_idx) == 0) EXIT
                str_out(i,1)(k:k) = CHAR(ascii_flat_out(char_idx))
            END DO
        END DO
        CALL assert_string_equal(str_out(1,1), "B1", "c_wrapper: val(1,1)")
        CALL assert_string_equal(str_out(2,1), "B2", "c_wrapper: val(2,1)")

        DEALLOCATE(ascii_flat_in, ascii_flat_out)
    END SUBROUTINE test_read_char_cols_c_wrapper

END MODULE mod_test_csv_read_char