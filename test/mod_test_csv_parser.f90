!> Unit test suite for csv_parser_module.
MODULE mod_test_csv_parser
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_parser_module, ONLY: parse_line, MAX_FIELD_LEN
    USE asserts, ONLY: assert_true, assert_equal_int, assert_string_equal

    IMPLICIT NONE
    PUBLIC

    ABSTRACT INTERFACE
        SUBROUTINE test_interface()
        END SUBROUTINE test_interface
    END INTERFACE

    TYPE :: test_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: test_proc => NULL()
    END TYPE test_case

CONTAINS

    !> Get array of all available tests for the csv_parser suite.
    FUNCTION get_all_tests_csv_parser() RESULT(all_tests)
        TYPE(test_case) :: all_tests(3)
        all_tests(1) = test_case("test_simple_line", test_simple_line)
        all_tests(2) = test_case("test_empty_field", test_empty_field)
        all_tests(3) = test_case("test_complex_field", test_complex_field)
    END FUNCTION get_all_tests_csv_parser

    !> Run all tests in this module.
    SUBROUTINE run_all_tests_csv_parser()
        TYPE(test_case) :: all_tests(3)
        INTEGER :: i
        all_tests = get_all_tests_csv_parser()
        WRITE(*, '(A)') "--- Running Suite: csv_parser ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: csv_parser ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_csv_parser

    !> Run specific tests by name.
    SUBROUTINE run_named_tests_csv_parser(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(3)
        INTEGER :: i, j
        LOGICAL :: found
        all_tests = get_all_tests_csv_parser()
        WRITE(*, '(A)') "--- Running Suite: csv_parser (named tests) ---"
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
                WRITE(*,*) "WARNING: Unknown test '", TRIM(test_names(i)), "' in suite 'csv_parser'"
            END IF
        END DO
        WRITE(*, '(A)') "--- Suite Finished: csv_parser ---"
        WRITE(*,*)
    END SUBROUTINE run_named_tests_csv_parser

    !> Test parsing a simple line.
    SUBROUTINE test_simple_line()
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
        INTEGER(INT32) :: status
        CALL parse_line("a,b,c", ",", fields, status)
        CALL assert_equal_int(status, 0, "test_simple_line: status was not 0")
        CALL assert_equal_int(SIZE(fields), 3, "test_simple_line: incorrect number of fields")
        CALL assert_string_equal(fields(1), "a", "test_simple_line: field 1 incorrect")
        CALL assert_string_equal(fields(2), "b", "test_simple_line: field 2 incorrect")
        CALL assert_string_equal(fields(3), "c", "test_simple_line: field 3 incorrect")
    END SUBROUTINE test_simple_line

    !> Test parsing a line with an empty field.
    SUBROUTINE test_empty_field()
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
        INTEGER(INT32) :: status
        CALL parse_line("a,,c", ",", fields, status)
        CALL assert_equal_int(status, 0, "test_empty_field: status was not 0")
        CALL assert_equal_int(SIZE(fields), 3, "test_empty_field: incorrect number of fields")
        CALL assert_string_equal(fields(2), "", "test_empty_field: empty field was not empty")
    END SUBROUTINE test_empty_field

    !> Test parsing a complex number correctly.
    SUBROUTINE test_complex_field()
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
        INTEGER(INT32) :: status
        CALL parse_line("1,(1.0, 2.0),3", ",", fields, status)
        CALL assert_equal_int(status, 0, "test_complex_field: status was not 0")
        CALL assert_equal_int(SIZE(fields), 3, "test_complex_field: incorrect number of fields")
        CALL assert_string_equal(fields(2), "(1.0, 2.0)", "test_complex_field: complex field parsed incorrectly")
    END SUBROUTINE test_complex_field

END MODULE mod_test_csv_parser