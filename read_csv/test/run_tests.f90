!> @brief Main program to run the entire test suite.
!> @details This program acts as a test runner. It can execute all tests,
!>          a specific suite of tests, or named tests within a suite based on
!>          command-line arguments.
PROGRAM run_tests
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32
    ! Import the test suites that should be included in the runner.
    USE mod_test_csv_parser
    USE mod_test_csv_file_reader
    USE mod_test_csv_read_int
    USE mod_test_csv_read_real
    USE mod_test_csv_read_logical
    USE mod_test_csv_read_char
    USE mod_test_csv_read_complex

    IMPLICIT NONE

    ! Abstract interface for the procedure pointers in the test_suite type.
    ABSTRACT INTERFACE
        SUBROUTINE run_all_tests_interface()
        END SUBROUTINE
        SUBROUTINE run_named_tests_interface(test_names)
            CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        END SUBROUTINE
    END INTERFACE

    ! Holds the name and procedure pointers for a test suite.
    TYPE :: test_suite
        CHARACTER(LEN=64) :: name
        PROCEDURE(run_all_tests_interface), POINTER, NOPASS :: run_all => NULL()
        PROCEDURE(run_named_tests_interface), POINTER, NOPASS :: run_named => NULL()
    END TYPE test_suite

    ! A registry to hold all available test suites.
    TYPE(test_suite), ALLOCATABLE :: available_suites(:)

    ! --- Variables for the Main Program Body ---
    INTEGER :: i, n_args
    CHARACTER(LEN=256) :: arg
    LOGICAL :: suite_found
    ! CORRECTED: Declaration moved to the top of the program block.
    CHARACTER(LEN=64), ALLOCATABLE :: test_names_to_run(:)

    ! --- Main Program Executable Logic ---
    CALL initialize_suites()
    n_args = command_argument_count()

    ! Case 1: No arguments, run all tests from all suites.
    IF (n_args == 0) THEN
        WRITE(*, '(A)') "Running all tests from all suites..."
        WRITE(*,*)
        DO i = 1, SIZE(available_suites)
            CALL available_suites(i)%run_all()
        END DO
        WRITE(*, '(A)') "All tests completed."
        STOP 0
    END IF

    ! Case 2: One or two arguments (suite name and optional test names).
    IF (n_args >= 1) THEN
        CALL get_command_argument(1, arg)
        suite_found = .FALSE.
        DO i = 1, SIZE(available_suites)
            IF (TRIM(arg) == TRIM(available_suites(i)%name)) THEN
                suite_found = .TRUE.
                ! Case 2a: Just a suite name, run all tests in that suite.
                IF (n_args == 1) THEN
                    CALL available_suites(i)%run_all()
                ! Case 2b: Suite name and test names, run specific tests.
                ELSE IF (n_args == 2) THEN
                    CALL get_command_argument(2, arg)
                    CALL parse_test_names(arg, test_names_to_run)
                    CALL available_suites(i)%run_named(test_names_to_run)
                END IF
                EXIT
            END IF
        END DO
        IF (.NOT. suite_found) THEN
            WRITE(*,*) "ERROR: Test suite '", TRIM(arg), "' not found."
            STOP 1
        END IF
    END IF

CONTAINS

    !> @brief Adds a new test suite to the global registry.
    SUBROUTINE add_suite(name, run_all_proc, run_named_proc)
        CHARACTER(LEN=*), INTENT(IN) :: name
        PROCEDURE(run_all_tests_interface) :: run_all_proc
        PROCEDURE(run_named_tests_interface) :: run_named_proc
        
        TYPE(test_suite), ALLOCATABLE :: new_suites(:)
        INTEGER :: old_size
        
        old_size = 0
        IF (ALLOCATED(available_suites)) THEN
            old_size = SIZE(available_suites)
        END IF
        
        ALLOCATE(new_suites(old_size + 1))
        IF (old_size > 0) THEN
            new_suites(1:old_size) = available_suites
        END IF
        
        new_suites(old_size + 1)%name = name
        new_suites(old_size + 1)%run_all => run_all_proc
        new_suites(old_size + 1)%run_named => run_named_proc
        
        CALL move_alloc(new_suites, available_suites)
    END SUBROUTINE add_suite

    !> @brief Initializes the registry with all test suites.
    !> @note To add a new suite to the test runner, register it here.
    SUBROUTINE initialize_suites()
        ! This is where you register each test suite module.
    CALL add_suite("csv_parser", run_all_tests_csv_parser, run_named_tests_csv_parser)
    CALL add_suite("csv_file_reader", run_all_tests_csv_file_reader, run_named_tests_csv_file_reader)
    CALL add_suite("csv_read_int", run_all_tests_csv_read_int, run_named_tests_csv_read_int)
    CALL add_suite("csv_read_real",run_all_tests_csv_read_real, run_named_tests_csv_read_real)
    CALL add_suite("csv_read_logical",run_all_tests_csv_read_logical, run_named_tests_csv_read_logical)
    CALL add_suite("csv_read_char",run_all_tests_csv_read_char, run_named_tests_csv_read_char)
    CALL add_suite("csv_read_complex",run_all_tests_csv_read_complex, run_named_tests_csv_read_complex)

        ! Add future suites here, e.g.:
        ! CALL add_suite("another_module", run_all_tests_another_module, run_named_tests_another_module)
    END SUBROUTINE initialize_suites

    !> @brief Parses a comma-separated string of test names into an array.
    SUBROUTINE parse_test_names(names_str, names_arr)
        CHARACTER(LEN=*), INTENT(IN) :: names_str
        CHARACTER(LEN=64), ALLOCATABLE, INTENT(OUT) :: names_arr(:)
        
        CHARACTER(LEN=LEN(names_str)) :: temp_str
        INTEGER :: count, i, start
        
        temp_str = TRIM(names_str)
        count = 1
        DO i = 1, LEN_TRIM(temp_str)
            IF (temp_str(i:i) == ',') THEN
                count = count + 1
            END IF
        END DO
        
        ALLOCATE(names_arr(count))
        start = 1
        count = 1
        DO i = 1, LEN_TRIM(temp_str)
            IF (temp_str(i:i) == ',') THEN
                names_arr(count) = temp_str(start:i-1)
                start = i + 1
                count = count + 1
            END IF
        END DO
        names_arr(count) = temp_str(start:)
    END SUBROUTINE parse_test_names

END PROGRAM run_tests
