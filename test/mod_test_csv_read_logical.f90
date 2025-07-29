!> @brief Unit test suite for csv_read_logical_module and its wrappers.
MODULE mod_test_csv_read_logical
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE, INTRINSIC :: iso_c_binding
    USE csv_read_logical_module, ONLY: read_logical_columns
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
    USE asserts, ONLY: assert_equal_int, assert_equal_size, assert_true, assert_false

    IMPLICIT NONE
    PUBLIC

    INTERFACE
        SUBROUTINE read_logical_columns_c(data_in_ascii, num_rows, num_cols_in, &
                                          cols_to_read, num_cols_to_read, &
                                          logical_data_out, status) &
                                          BIND(C, NAME='read_logical_columns_c')
            USE, INTRINSIC :: iso_c_binding
            INTEGER(C_INT), INTENT(IN) :: data_in_ascii(*)
            INTEGER(C_INT), VALUE, INTENT(IN) :: num_rows, num_cols_in, num_cols_to_read
            INTEGER(C_INT), INTENT(IN) :: cols_to_read(*)
            LOGICAL(C_BOOL), INTENT(OUT) :: logical_data_out(*)
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

    FUNCTION get_all_tests_csv_read_logical() RESULT(all_tests)
        TYPE(test_case) :: all_tests(2)
        all_tests(1) = test_case("test_read_logical_cols_basic", test_read_logical_cols_basic)
        all_tests(2) = test_case("test_read_logical_cols_c_wrapper", test_read_logical_cols_c_wrapper)
    END FUNCTION get_all_tests_csv_read_logical

    SUBROUTINE run_all_tests_csv_read_logical()
        TYPE(test_case) :: all_tests(2)
        INTEGER :: i
        all_tests = get_all_tests_csv_read_logical()
        WRITE(*, '(A)') "--- Running Suite: csv_read_logical ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: csv_read_logical ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_csv_read_logical

    SUBROUTINE run_named_tests_csv_read_logical(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(2)
        INTEGER :: i, j
        LOGICAL :: found
        all_tests = get_all_tests_csv_read_logical()
        WRITE(*, '(A)') "--- Running Suite: csv_read_logical (named tests) ---"
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
                WRITE(*,*) "WARNING: Unknown test '", TRIM(test_names(i)), "' in suite 'csv_read_logical'"
            END IF
        END DO
        WRITE(*, '(A)') "--- Suite Finished: csv_read_logical ---"
        WRITE(*,*)
    END SUBROUTINE run_named_tests_csv_read_logical

    SUBROUTINE test_read_logical_cols_basic()
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: data_in(:,:)
        LOGICAL, ALLOCATABLE :: logical_data_out(:,:)
        INTEGER(INT32) :: status
        INTEGER(INT32) :: cols_to_read(1)

        ALLOCATE(data_in(4, 2))
        data_in(1, 1) = "10"; data_in(1, 2) = ".TRUE."
        data_in(2, 1) = "20"; data_in(2, 2) = "T"
        data_in(3, 1) = "30"; data_in(3, 2) = ".FALSE."
        data_in(4, 1) = "40"; data_in(4, 2) = "F"
        cols_to_read = [2]
        CALL read_logical_columns(data_in, cols_to_read, logical_data_out, status)

        CALL assert_equal_int(status, 0, "test_read_logical_cols_basic: status")
        CALL assert_true(logical_data_out(1, 1), "test_read_logical_cols_basic: val(1,1)")
        CALL assert_true(logical_data_out(2, 1), "test_read_logical_cols_basic: val(2,1)")
        CALL assert_false(logical_data_out(3, 1), "test_read_logical_cols_basic: val(3,1)")
        CALL assert_false(logical_data_out(4, 1), "test_read_logical_cols_basic: val(4,1)")
        DEALLOCATE(data_in, logical_data_out)
    END SUBROUTINE test_read_logical_cols_basic

    SUBROUTINE test_read_logical_cols_c_wrapper()
        CHARACTER(LEN=MAX_FIELD_LEN) :: str_arr(2,2)
        INTEGER(C_INT), ALLOCATABLE :: ascii_flat_in(:)
        LOGICAL(C_BOOL), ALLOCATABLE :: logical_data_out(:)
        INTEGER(C_INT) :: status, i, j, k, char_idx
        INTEGER(C_INT) :: cols_to_read(2)

        ! Arrange
        str_arr(1,1) = "T"; str_arr(1,2) = "F"
        str_arr(2,1) = ".TRUE."; str_arr(2,2) = ".FALSE."
        ALLOCATE(ascii_flat_in(2 * 2 * MAX_FIELD_LEN))
        ascii_flat_in = 0
        DO j = 1, 2
            DO i = 1, 2
                DO k = 1, LEN_TRIM(str_arr(i,j))
                    char_idx = (i - 1) * 2 * MAX_FIELD_LEN + (j - 1) * MAX_FIELD_LEN + k
                    ascii_flat_in(char_idx) = ICHAR(str_arr(i,j)(k:k))
                END DO
            END DO
        END DO
        
        ALLOCATE(logical_data_out(2 * 2))
        cols_to_read = [1, 2]

        ! Act
        CALL read_logical_columns_c(ascii_flat_in, 2, 2, cols_to_read, 2, logical_data_out, status)

        ! Assert
        CALL assert_equal_int(status, 0, "c_wrapper: status")
        ! CORRECTED: Explicitly convert LOGICAL(C_BOOL) to default LOGICAL for assertion.
        CALL assert_true(LOGICAL(logical_data_out(1)), "c_wrapper: val(1,1)")  ! T
        CALL assert_true(LOGICAL(logical_data_out(3)), "c_wrapper: val(2,1)")  ! .TRUE.
        CALL assert_false(LOGICAL(logical_data_out(2)), "c_wrapper: val(1,2)") ! F
        CALL assert_false(LOGICAL(logical_data_out(4)), "c_wrapper: val(2,2)") ! .FALSE.

        DEALLOCATE(ascii_flat_in, logical_data_out)
    END SUBROUTINE test_read_logical_cols_c_wrapper

END MODULE mod_test_csv_read_logical
