!> @brief Unit test suite for csv_read_complex_module and its wrappers.
MODULE mod_test_csv_read_complex
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    USE, INTRINSIC :: iso_c_binding
    USE csv_read_complex_module, ONLY: read_complex_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    USE asserts, ONLY: assert_equal_int, assert_equal_size, assert_equal_complex

    IMPLICIT NONE
    PUBLIC

    INTERFACE
        SUBROUTINE read_complex_columns_c(data_in_ascii, num_rows, num_cols_in, &
                                          cols_to_read, num_cols_to_read, &
                                          complex_data_out, status) &
                                          BIND(C, NAME='read_complex_columns_c')
            USE, INTRINSIC :: iso_c_binding
            INTEGER(C_INT), INTENT(IN) :: data_in_ascii(*)
            INTEGER(C_INT), VALUE, INTENT(IN) :: num_rows, num_cols_in, num_cols_to_read
            INTEGER(C_INT), INTENT(IN) :: cols_to_read(*)
            COMPLEX(C_DOUBLE), INTENT(OUT) :: complex_data_out(*)
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

    FUNCTION get_all_tests_csv_read_complex() RESULT(all_tests)
        TYPE(test_case) :: all_tests(2)
        all_tests(1) = test_case("test_read_complex_cols_basic", test_read_complex_cols_basic)
        all_tests(2) = test_case("test_read_complex_cols_c_wrapper", test_read_complex_cols_c_wrapper)
    END FUNCTION get_all_tests_csv_read_complex

    SUBROUTINE run_all_tests_csv_read_complex()
        TYPE(test_case) :: all_tests(2)
        INTEGER :: i
        all_tests = get_all_tests_csv_read_complex()
        WRITE(*, '(A)') "--- Running Suite: csv_read_complex ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: csv_read_complex ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_csv_read_complex

    SUBROUTINE run_named_tests_csv_read_complex(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(2)
        INTEGER :: i, j
        LOGICAL :: found
        all_tests = get_all_tests_csv_read_complex()
        WRITE(*, '(A)') "--- Running Suite: csv_read_complex (named tests) ---"
        DO i = 1, SIZE(test_names)
            found = .false.
            DO j = 1, SIZE(all_tests)
                IF (TRIM(test_names(i)) == TRIM(all_tests(j)%name)) THEN
                    WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(j)%name), "..."
                    CALL all_tests(j)%test_proc()
                    WRITE(*, '(A)') " PASSED"
                    found = .true.
                    EXIT
                END IF
            END DO
            IF (.NOT. found) THEN
                WRITE(*,*) "WARNING: Unknown test '", TRIM(test_names(i)), "' in suite 'csv_read_complex'"
            END IF
        END DO
        WRITE(*, '(A)') "--- Suite Finished: csv_read_complex ---"
        WRITE(*,*)
    END SUBROUTINE run_named_tests_csv_read_complex

    SUBROUTINE test_read_complex_cols_basic()
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in(:,:)
        COMPLEX(REAL64), ALLOCATABLE :: complex_data_out(:,:)
        INTEGER(INT32) :: status
        INTEGER(INT32) :: cols_to_read(1)

        ALLOCATE(data_in(2, 2))
        data_in(1, 1) = "1"; data_in(1, 2) = "(1.5, -2.5)"
        data_in(2, 1) = "2"; data_in(2, 2) = "( 3.0, 4.0 )"
        cols_to_read = [2]
        CALL read_complex_columns(data_in, cols_to_read, complex_data_out, status)

        CALL assert_equal_int(status, 0, "test_read_complex_cols_basic: status")
        CALL assert_equal_complex(complex_data_out(1, 1), (1.5_REAL64, -2.5_REAL64), 1e-9_REAL64, &
                                  "test_read_complex_cols_basic: val(1,1)")
        CALL assert_equal_complex(complex_data_out(2, 1), (3.0_REAL64, 4.0_REAL64), 1e-9_REAL64, &
                                  "test_read_complex_cols_basic: val(2,1)")
        DEALLOCATE(data_in, complex_data_out)
    END SUBROUTINE test_read_complex_cols_basic

    SUBROUTINE test_read_complex_cols_c_wrapper()
        CHARACTER(LEN=MAX_FIELD_LEN) :: str_arr(2,1)
        INTEGER(C_INT), ALLOCATABLE :: ascii_flat_in(:)
        COMPLEX(C_DOUBLE), ALLOCATABLE :: complex_data_out(:)
        INTEGER(C_INT) :: status, i, k, char_idx
        INTEGER(C_INT) :: cols_to_read(1)

        ! Arrange
        str_arr(1,1) = "(1.1, 2.2)"
        str_arr(2,1) = "(3.3, -4.4)"
        ALLOCATE(ascii_flat_in(2 * 1 * MAX_FIELD_LEN))
        ascii_flat_in = 0
        DO i = 1, 2
            DO k = 1, LEN_TRIM(str_arr(i,1))
                char_idx = (i - 1) * 1 * MAX_FIELD_LEN + (1 - 1) * MAX_FIELD_LEN + k
                ascii_flat_in(char_idx) = ICHAR(str_arr(i,1)(k:k))
            END DO
        END DO
        
        ALLOCATE(complex_data_out(2 * 1))
        cols_to_read = [1]

        ! Act
        CALL read_complex_columns_c(ascii_flat_in, 2, 1, cols_to_read, 1, complex_data_out, status)

        ! Assert
        CALL assert_equal_int(status, 0, "c_wrapper: status")
        CALL assert_equal_complex(complex_data_out(1), (1.1_REAL64, 2.2_REAL64), 1e-9_REAL64, "c_wrapper: val(1,1)")
        CALL assert_equal_complex(complex_data_out(2), (3.3_REAL64, -4.4_REAL64), 1e-9_REAL64, "c_wrapper: val(2,1)")

        DEALLOCATE(ascii_flat_in, complex_data_out)
    END SUBROUTINE test_read_complex_cols_c_wrapper

END MODULE mod_test_csv_read_complex
