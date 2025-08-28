!> Unit test suite for csv_file_reader_module.
MODULE mod_test_csv_file_reader
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    USE csv_file_reader_module, ONLY: read_csv_to_strings
    USE csv_parser_module, ONLY: MAX_FIELD_LEN
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

    !> Helper to create a dummy file for testing.
    SUBROUTINE create_dummy_file(filename, content)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=*), INTENT(IN) :: content
        INTEGER :: unit
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
        WRITE(unit, '(A)') content
        CLOSE(unit)
    END SUBROUTINE create_dummy_file

    !> Helper to delete a dummy file after testing.
    SUBROUTINE delete_dummy_file(filename)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        LOGICAL :: ex
        INTEGER :: ierr
        INQUIRE(FILE=TRIM(filename), EXIST=ex)
        IF (ex) THEN
            OPEN(UNIT=10, FILE=TRIM(filename), STATUS='OLD')
            CLOSE(UNIT=10, STATUS='DELETE', IOSTAT=ierr)
        END IF
    END SUBROUTINE delete_dummy_file

    !> Get array of all available tests for this suite.
    FUNCTION get_all_tests_csv_file_reader() RESULT(all_tests)
        TYPE(test_case) :: all_tests(1)
        all_tests(1) = test_case("test_read_basic_csv", test_read_basic_csv)
    END FUNCTION get_all_tests_csv_file_reader

    !> Run all tests in this module.
    SUBROUTINE run_all_tests_csv_file_reader()
        TYPE(test_case) :: all_tests(1)
        INTEGER :: i
        all_tests = get_all_tests_csv_file_reader()
        WRITE(*, '(A)') "--- Running Suite: csv_file_reader ---"
        DO i = 1, SIZE(all_tests)
            WRITE(*, '(A, A, A)', ADVANCE='NO') "  Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*, '(A)') " PASSED"
        END DO
        WRITE(*, '(A)') "--- Suite PASSED: csv_file_reader ---"
        WRITE(*,*)
    END SUBROUTINE run_all_tests_csv_file_reader

    !> Run specific tests by name.
    SUBROUTINE run_named_tests_csv_file_reader(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(1)
        INTEGER :: i, j
        LOGICAL :: found
        all_tests = get_all_tests_csv_file_reader()
        WRITE(*, '(A)') "--- Running Suite: csv_file_reader (named tests) ---"
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
                WRITE(*,*) "WARNING: Unknown test '", TRIM(test_names(i)), "' in suite 'csv_file_reader'"
            END IF
        END DO
        WRITE(*, '(A)') "--- Suite Finished: csv_file_reader ---"
        WRITE(*,*)
    END SUBROUTINE run_named_tests_csv_file_reader

    !> Test basic CSV file reading with a header.
    SUBROUTINE test_read_basic_csv()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=1024) :: content
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: header(:), data(:,:)
        INTEGER(INT32) :: status

        filename = "test.csv"
        content = "ID,Name,Value" // ACHAR(10) // "1,Alpha,100" // ACHAR(10) // "2,Beta,200"
        CALL create_dummy_file(filename, content)

        CALL read_csv_to_strings(filename, .TRUE., ",", header, data, status)

        CALL assert_equal_int(status, 0, "test_read_basic_csv: status was not 0")
        CALL assert_equal_int(SIZE(header), 3, "test_read_basic_csv: incorrect header size")
        CALL assert_equal_int(SIZE(data, DIM=1), 2, "test_read_basic_csv: incorrect number of rows")
        CALL assert_equal_int(SIZE(data, DIM=2), 3, "test_read_basic_csv: incorrect number of columns")
        CALL assert_string_equal(header(2), "Name", "test_read_basic_csv: incorrect header value")
        CALL assert_string_equal(data(2, 3), "200", "test_read_basic_csv: incorrect data value")

        CALL delete_dummy_file(filename)
    END SUBROUTINE test_read_basic_csv

END MODULE mod_test_csv_file_reader