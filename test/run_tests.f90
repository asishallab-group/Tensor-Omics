! filepath: test/run_tests.f90
PROGRAM main
    USE asserts
    USE mod_test_csv_reader ! Your new test suite
    ! Add other test modules here if you create more, e.g.:
    ! USE mod_test_normalize_by_std_dev
    ! USE mod_test_sorting

    IMPLICIT NONE

    ! Type to hold a suite's name and its run functions
    TYPE :: test_suite_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: run_all_proc => NULL()
        PROCEDURE(test_named_interface), POINTER, NOPASS :: run_named_proc => NULL()
    END TYPE test_suite_case

    ! Abstract interface for running named tests
    ABSTRACT INTERFACE
        SUBROUTINE test_named_interface(test_names)
            CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        END SUBROUTINE test_named_interface
    END INTERFACE

    TYPE(test_suite_case), ALLOCATABLE :: available_suites(:)
    
    CHARACTER(LEN=256), ALLOCATABLE :: args(:)
    CHARACTER(LEN=64) :: suite_name_arg
    CHARACTER(LEN=64), ALLOCATABLE :: test_names_arg(:)
    INTEGER :: num_args
    INTEGER :: i, j

    ! Private procedures for this program
    PRIVATE :: initialize_suites, add_suite, parse_command_line_args

CONTAINS

    !> @brief Initializes the registry of available test suites.
    SUBROUTINE initialize_suites()
        ! Start with empty registry
        IF (ALLOCATED(available_suites)) DEALLOCATE(available_suites)
        ALLOCATE(available_suites(0))

        ! Add each suite
        CALL add_suite("csv_reader", run_all_tests_csv_reader, run_named_tests_csv_reader)
        ! Add other suites here as you create them:
        ! CALL add_suite("normalization", run_all_tests_normalize_by_std_dev, run_named_tests_normalize_by_std_dev)
        ! CALL add_suite("sorting", run_all_tests_sorting, run_named_tests_sorting)
    END SUBROUTINE initialize_suites

    !> @brief Adds a new test suite to the registry.
    SUBROUTINE add_suite(name, run_all_proc, run_named_proc)
        CHARACTER(LEN=*), INTENT(IN) :: name
        PROCEDURE(test_interface), POINTER :: run_all_proc
        PROCEDURE(test_named_interface), POINTER :: run_named_proc

        INTEGER :: old_size
        TYPE(test_suite_case), ALLOCATABLE :: temp_suites(:)

        old_size = SIZE(available_suites)
        ALLOCATE(temp_suites(old_size + 1))
        IF (old_size > 0) temp_suites(1:old_size) = available_suites
        
        temp_suites(old_size + 1)%name = name
        temp_suites(old_size + 1)%run_all_proc => run_all_proc
        temp_suites(old_size + 1)%run_named_proc => run_named_proc

        IF (ALLOCATED(available_suites)) DEALLOCATE(available_suites)
        available_suites = temp_suites
        IF (ALLOCATED(temp_suites)) DEALLOCATE(temp_suites)
    END SUBROUTINE add_suite

    !> @brief Parses command line arguments to determine tests to run.
    SUBROUTINE parse_command_line_args(suite_name, test_names, num_found_args)
        CHARACTER(LEN=*), INTENT(OUT) :: suite_name
        CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: test_names(:)
        INTEGER, INTENT(OUT) :: num_found_args

        CHARACTER(LEN=256) :: arg_str
        INTEGER :: i_arg, start_pos, end_pos, num_commas

        num_found_args = COMMAND_ARGUMENT_COUNT()
        suite_name = ""
        IF (ALLOCATED(test_names)) DEALLOCATE(test_names)
        ALLOCATE(test_names(0))

        IF (num_found_args == 0) THEN
            ! Run all tests, no specific suite or tests
            RETURN
        END IF

        ! First argument is suite name
        CALL GET_COMMAND_ARGUMENT(1, VALUE=arg_str)
        suite_name = TRIM(arg_str)

        IF (num_found_args == 1) THEN
            ! Run all tests in specified suite
            RETURN
        END IF

        ! Remaining arguments are specific test names (comma-separated if only one arg)
        ! This assumes "test1,test2,test3" is passed as a single argument string
        CALL GET_COMMAND_ARGUMENT(2, VALUE=arg_str)
        
        num_commas = COUNT(arg_str(1:LEN_TRIM(arg_str)) == ',')
        IF (num_commas == 0) THEN
            ! Single test name or first part of comma-separated string
            ALLOCATE(test_names(1))
            test_names(1) = TRIM(arg_str)
        ELSE
            ! Comma-separated list of test names
            ALLOCATE(test_names(num_commas + 1))
            start_pos = 1
            DO i_arg = 1, num_commas + 1
                end_pos = INDEX(arg_str(start_pos:), ',')
                IF (end_pos == 0) THEN ! Last part
                    test_names(i_arg) = TRIM(arg_str(start_pos:))
                    EXIT
                ELSE
                    test_names(i_arg) = TRIM(arg_str(start_pos : start_pos + end_pos - 2))
                    start_pos = start_pos + end_pos
                END IF
            END DO
        END IF

    END SUBROUTINE parse_command_line_args

    ! --- Main Test Execution Logic ---
    BLOCK MAIN_TEST_EXECUTION
        INTEGER :: i_suite
        LOGICAL :: suite_found
        CHARACTER(LEN=64) :: parsed_suite_name
        CHARACTER(LEN=64), ALLOCATABLE :: parsed_test_names(:)
        INTEGER :: num_parsed_args

        CALL initialize_suites()
        CALL parse_command_line_args(parsed_suite_name, parsed_test_names, num_parsed_args)

        IF (num_parsed_args == 0) THEN
            ! Run all tests from all suites
            WRITE(*,*) "Running all tests from all suites..."
            DO i_suite = 1, SIZE(available_suites)
                CALL available_suites(i_suite)%run_all_proc()
            END DO
        ELSE
            ! Run specific suite or specific tests within a suite
            suite_found = .FALSE.
            DO i_suite = 1, SIZE(available_suites)
                IF (TRIM(parsed_suite_name) == TRIM(available_suites(i_suite)%name)) THEN
                    suite_found = .TRUE.
                    IF (num_parsed_args == 1) THEN
                        ! Run all tests in the specified suite
                        WRITE(*,*) "Running all tests in suite: ", TRIM(parsed_suite_name), "..."
                        CALL available_suites(i_suite)%run_all_proc()
                    ELSE
                        ! Run specific tests within the specified suite
                        WRITE(*,*) "Running tests in suite '", TRIM(parsed_suite_name), &
                                   "': ", TRIM(parsed_test_names(1)), &
                                   " (and possibly others) ..."
                        CALL available_suites(i_suite)%run_named_proc(parsed_test_names)
                    END IF
                    EXIT
                END IF
            END DO
            IF (.NOT. suite_found) THEN
                WRITE(*,*) "Error: Unknown test suite: ", TRIM(parsed_suite_name)
                ERROR STOP "Unknown suite."
            END IF
        END IF

        WRITE(*,*) "All requested tests completed successfully."
    END BLOCK MAIN_TEST_EXECUTION

END PROGRAM main