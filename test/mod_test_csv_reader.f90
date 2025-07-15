! filepath: test/mod_test_csv_reader.f90
MODULE mod_test_csv_reader
    USE asserts
    USE csv_reader_module, ONLY: read_csv_file, cleanup_csv_data, get_num_rows, get_num_cols, &
                                 get_header, get_column_data_type_by_index, get_column_data_type_by_name, &
                                 get_int_column, get_real_column, get_char_column, get_logical_column, &
                                 get_complex_column, get_cell, serialize, deserialize, &
                                 IK, MAX_FIELD_LEN, DEFAULT_DELIMITER

    IMPLICIT NONE
    
    ! Abstract interface for all test procedures
    ABSTRACT INTERFACE
        SUBROUTINE test_interface()
        END SUBROUTINE test_interface
    END INTERFACE

    ! Type to hold test name and procedure pointer
    TYPE :: test_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: test_proc => NULL()
    END TYPE test_case

    PRIVATE

    PUBLIC :: get_all_tests_csv_reader, run_all_tests_csv_reader, run_named_tests_csv_reader

CONTAINS

    ! Helper subroutine to create a dummy CSV file
    SUBROUTINE create_dummy_csv(filename, content)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=*), INTENT(IN) :: content
        INTEGER :: unit
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE')
        WRITE(unit, '(A)') TRIM(content)
        CLOSE(unit)
    END SUBROVEINE create_dummy_csv

    ! Helper subroutine to clean up dummy CSV files
    SUBROUTINE delete_dummy_csv(filename)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER :: stat
        INQUIRE(FILE=TRIM(filename), EXIST=stat)
        IF (stat == 0) THEN
            CLOSE(FILE=TRIM(filename)) ! Ensure closed before deleting
            CALL SYSTEM('rm ' // TRIM(filename), STATUS=stat) ! For Linux/macOS
            ! For Windows: CALL SYSTEM('del ' // TRIM(filename), STATUS=stat)
        END IF
    END SUBROUTINE delete_dummy_csv


    !> @brief Test basic CSV reading with header.
    SUBROUTINE test_basic_read_with_header()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, num_rows, num_cols
        CHARACTER(LEN=MAX_FIELD_LEN) :: header_name
        REAL :: r_val
        LOGICAL :: l_val
        COMPLEX :: co_val
        INTEGER(IK) :: i_val

        filename = "test_basic_header.csv"
        content = "ID,Name,Value,Active,Coord" // NEW_LINE // &
                  "1,Alpha,10.5,.TRUE.,(1.0,2.0)" // NEW_LINE // &
                  "2,Beta,20.0,.FALSE.,(3.0,-4.0)"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_basic_read_with_header: read status")
        
        num_rows = get_num_rows()
        num_cols = get_num_cols()
        CALL assert_equal_int(num_rows, 2, "test_basic_read_with_header: num_rows")
        CALL assert_equal_int(num_cols, 5, "test_basic_read_with_header: num_cols")

        header_name = get_header(2)
        CALL assert_true(TRIM(header_name) == "Name", "test_basic_read_with_header: header name")

        i_val = get_cell(1, 1)%i_val
        CALL assert_equal_int(i_val, 1, "test_basic_read_with_header: cell(1,1) int")
        
        r_val = get_cell(1, 3)%r_val
        CALL assert_equal_real(r_val, 10.5, 1.0E-6, "test_basic_read_with_header: cell(1,3) real")

        l_val = get_cell(1,4)%l_val
        CALL assert_true(l_val, "test_basic_read_with_header: cell(1,4) logical")

        co_val = get_cell(1,5)%co_val
        CALL assert_equal_real(REAL(co_val), 1.0, 1.0E-6, "test_basic_read_with_header: cell(1,5) complex real part")
        CALL assert_equal_real(AIMAG(co_val), 2.0, 1.0E-6, "test_basic_read_with_header: cell(1,5) complex imag part")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_basic_read_with_header

    !> @brief Test CSV reading without header.
    SUBROUTINE test_read_no_header()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, num_rows, num_cols
        INTEGER(IK) :: i_val
        REAL :: r_val

        filename = "test_no_header.csv"
        content = "1,Alpha,10.5" // NEW_LINE // &
                  "2,Beta,20.0"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .FALSE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_read_no_header: read status")
        
        num_rows = get_num_rows()
        num_cols = get_num_cols()
        CALL assert_equal_int(num_rows, 2, "test_read_no_header: num_rows")
        CALL assert_equal_int(num_cols, 3, "test_read_no_header: num_cols")

        i_val = get_cell(1, 1)%i_val
        CALL assert_equal_int(i_val, 1, "test_read_no_header: cell(1,1) int")
        
        r_val = get_cell(1, 3)%r_val
        CALL assert_equal_real(r_val, 10.5, 1.0E-6, "test_read_no_header: cell(1,3) real")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_read_no_header

    !> @brief Test CSV reading with TAB delimiter.
    SUBROUTINE test_read_tab_delimiter()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, num_rows, num_cols
        CHARACTER(LEN=MAX_FIELD_LEN) :: char_val
        
        filename = "test_tab_delim.csv"
        content = "ID" // DEFAULT_DELIMITER // "Name" // NEW_LINE // &
                  "1" // DEFAULT_DELIMITER // "Alpha"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=DEFAULT_DELIMITER, io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_read_tab_delimiter: read status")
        
        num_rows = get_num_rows()
        num_cols = get_num_cols()
        CALL assert_equal_int(num_rows, 1, "test_read_tab_delimiter: num_rows")
        CALL assert_equal_int(num_cols, 2, "test_read_tab_delimiter: num_cols")

        char_val = get_cell(1, 2)%c_val
        CALL assert_true(TRIM(char_val) == "Alpha", "test_read_tab_delimiter: cell(1,2) char")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_read_tab_delimiter

    !> @brief Test error handling for non-existent file.
    SUBROUTINE test_read_file_not_found()
        CHARACTER(LEN=64) :: filename
        INTEGER :: io_status

        filename = "non_existent.csv"
        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 1, "test_read_file_not_found: io_status for missing file")
    END SUBROUTINE test_read_file_not_found

    !> @brief Test parsing with empty fields.
    SUBROUTINE test_empty_fields()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        CHARACTER(LEN=MAX_FIELD_LEN) :: char_val

        filename = "test_empty_fields.csv"
        content = "Col1,Col2,Col3" // NEW_LINE // &
                  "1,,Text"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_empty_fields: read status")
        CALL assert_equal_int(get_num_rows(), 1, "test_empty_fields: num_rows")
        CALL assert_equal_int(get_num_cols(), 3, "test_empty_fields: num_cols")

        char_val = get_cell(1, 2)%c_val
        CALL assert_true(TRIM(char_val) == "", "test_empty_fields: empty cell should be empty string")
        
        char_val = get_cell(1, 3)%c_val
        CALL assert_true(TRIM(char_val) == "Text", "test_empty_fields: cell(1,3) char")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_empty_fields

    !> @brief Test get_column_data_type_by_index and get_column_data_type_by_name.
    SUBROUTINE test_get_column_types()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, col_type

        filename = "test_col_types.csv"
        content = "IntCol,RealCol,CharCol,LogicalCol,ComplexCol" // NEW_LINE // &
                  "1,1.0,A,.TRUE.,(1.0,1.0)"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_get_column_types: read status")
        
        col_type = get_column_data_type_by_index(1)
        CALL assert_equal_int(col_type, 1, "test_get_column_types: IntCol type by index")
        
        col_type = get_column_data_type_by_name("RealCol")
        CALL assert_equal_int(col_type, 2, "test_get_column_types: RealCol type by name")

        col_type = get_column_data_type_by_name("LogicalCol")
        CALL assert_equal_int(col_type, 4, "test_get_column_types: LogicalCol type by name")
        
        col_type = get_column_data_type_by_name("ComplexCol")
        CALL assert_equal_int(col_type, 5, "test_get_column_types: ComplexCol type by name")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_get_column_types

    !> @brief Test column getters.
    SUBROUTINE test_column_getters()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        INTEGER(IK), ALLOCATABLE :: int_col(:)
        REAL, ALLOCATABLE :: real_col(:)
        LOGICAL, ALLOCATABLE :: logical_col(:)
        COMPLEX, ALLOCATABLE :: complex_col(:)
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: char_col(:)

        filename = "test_column_getters.csv"
        content = "Ints,Reals,Chars,Logicals,Complexes" // NEW_LINE // &
                  "1,1.1,A,.TRUE.,(1.0,1.0)" // NEW_LINE // &
                  "2,2.2,B,.FALSE.,(2.0,2.0)"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_column_getters: read status")

        int_col = get_int_column(1)
        CALL assert_equal_int(int_col(1), 1, "test_column_getters: int_col(1)")
        CALL assert_equal_int(int_col(2), 2, "test_column_getters: int_col(2)")
        IF (ALLOCATED(int_col)) DEALLOCATE(int_col)

        real_col = get_real_column(2)
        CALL assert_equal_real(real_col(1), 1.1, 1E-6, "test_column_getters: real_col(1)")
        CALL assert_equal_real(real_col(2), 2.2, 1E-6, "test_column_getters: real_col(2)")
        IF (ALLOCATED(real_col)) DEALLOCATE(real_col)

        char_col = get_char_column(3)
        CALL assert_true(TRIM(char_col(1)) == "A", "test_column_getters: char_col(1)")
        CALL assert_true(TRIM(char_col(2)) == "B", "test_column_getters: char_col(2)")
        IF (ALLOCATED(char_col)) DEALLOCATE(char_col)

        logical_col = get_logical_column(4)
        CALL assert_true(logical_col(1), "test_column_getters: logical_col(1)")
        CALL assert_false(logical_col(2), "test_column_getters: logical_col(2)")
        IF (ALLOCATED(logical_col)) DEALLOCATE(logical_col)

        complex_col = get_complex_column(5)
        CALL assert_equal_real(REAL(complex_col(1)), 1.0, 1E-6, "test_column_getters: complex_col(1) real")
        CALL assert_equal_real(AIMAG(complex_col(1)), 1.0, 1E-6, "test_column_getters: complex_col(1) imag")
        CALL assert_equal_real(REAL(complex_col(2)), 2.0, 1E-6, "test_column_getters: complex_col(2) real")
        CALL assert_equal_real(AIMAG(complex_col(2)), 2.0, 1E-6, "test_column_getters: complex_col(2) imag")
        IF (ALLOCATED(complex_col)) DEALLOCATE(complex_col)

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_column_getters

    !> @brief Test complex number parsing with whitespace.
    SUBROUTINE test_complex_whitespace_parsing()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        COMPLEX :: co_val

        filename = "test_complex_ws.csv"
        content = "Coord" // NEW_LINE // &
                  "( 1.0 , 2.0 )" // NEW_LINE // &
                  " ( 3.0 , -4.0 ) "
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_complex_whitespace_parsing: read status")
        CALL assert_equal_int(get_num_rows(), 2, "test_complex_whitespace_parsing: num_rows")
        
        co_val = get_cell(1, 1)%co_val
        CALL assert_equal_real(REAL(co_val), 1.0, 1.0E-6, "test_complex_whitespace_parsing: cell(1,1) real part")
        CALL assert_equal_real(AIMAG(co_val), 2.0, 1.0E-6, "test_complex_whitespace_parsing: cell(1,1) imag part")

        co_val = get_cell(2, 1)%co_val
        CALL assert_equal_real(REAL(co_val), 3.0, 1.0E-6, "test_complex_whitespace_parsing: cell(2,1) real part")
        CALL assert_equal_real(AIMAG(co_val), -4.0, 1.0E-6, "test_complex_whitespace_parsing: cell(2,1) imag part")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_complex_whitespace_parsing

    !> @brief Test serialization and deserialization.
    SUBROUTINE test_serialize_deserialize()
        CHARACTER(LEN=64) :: csv_filename, bin_filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        INTEGER(IK) :: i_val

        csv_filename = "test_serial.csv"
        bin_filename = "test_serial.bin"
        content = "A,B" // NEW_LINE // &
                  "1,10.5" // NEW_LINE // &
                  "2,20.0"
        CALL create_dummy_csv(csv_filename, content)

        ! Read data
        CALL read_csv_file(csv_filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_serialize_deserialize: read status")
        CALL assert_equal_int(get_num_rows(), 2, "test_serialize_deserialize: num_rows before serial")
        
        ! Serialize
        CALL serialize(bin_filename, io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_serialize_deserialize: serialize status")
        
        ! Cleanup and deserialize
        CALL cleanup_csv_data()
        CALL assert_equal_int(get_num_rows(), 0, "test_serialize_deserialize: num_rows after cleanup")

        CALL deserialize(bin_filename, io_status=io_status)
        CALL assert_equal_int(io_status, 0, "test_serialize_deserialize: deserialize status")
        
        CALL assert_equal_int(get_num_rows(), 2, "test_serialize_deserialize: num_rows after deserial")
        CALL assert_equal_int(get_num_cols(), 2, "test_serialize_deserialize: num_cols after deserial")
        
        i_val = get_cell(1, 1)%i_val
        CALL assert_equal_int(i_val, 1, "test_serialize_deserialize: cell(1,1) after deserial")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(csv_filename)
        CALL delete_dummy_csv(bin_filename)
    END SUBROUTINE test_serialize_deserialize

    !> @brief Get array of all available tests for csv_reader.
    FUNCTION get_all_tests_csv_reader() RESULT(all_tests)
        TYPE(test_case) :: all_tests(9) ! Number of tests in this module

        all_tests(1) = test_case("test_basic_read_with_header", test_basic_read_with_header)
        all_tests(2) = test_case("test_read_no_header", test_read_no_header)
        all_tests(3) = test_case("test_read_tab_delimiter", test_read_tab_delimiter)
        all_tests(4) = test_case("test_read_file_not_found", test_read_file_not_found)
        all_tests(5) = test_case("test_empty_fields", test_empty_fields)
        all_tests(6) = test_case("test_get_column_types", test_get_column_types)
        all_tests(7) = test_case("test_column_getters", test_column_getters)
        all_tests(8) = test_case("test_complex_whitespace_parsing", test_complex_whitespace_parsing)
        all_tests(9) = test_case("test_serialize_deserialize", test_serialize_deserialize)
    END FUNCTION get_all_tests_csv_reader

    !> @brief Run all tests in this module.
    SUBROUTINE run_all_tests_csv_reader()
        TYPE(test_case) :: all_tests(9)
        INTEGER :: i

        all_tests = get_all_tests_csv_reader()
        CALL reset_test_failed_flag() ! Reset flag for current run

        DO i = 1, SIZE(all_tests)
            WRITE(*,*) "Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*,*) TRIM(all_tests(i)%name), " passed."
        END DO
        IF (.NOT. test_failed_flag) THEN
            WRITE(*,*) "All csv_reader_module tests passed successfully."
        ELSE
            WRITE(*,*) "Some csv_reader_module tests FAILED."
            CALL error_stop_test() ! Re-throw error if any test failed
        END IF
    END SUBROUTINE run_all_tests_csv_reader

    !> @brief Run specific tests by name.
    SUBROUTINE run_named_tests_csv_reader(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(9)
        INTEGER :: i, j
        LOGICAL :: found

        all_tests = get_all_tests_csv_reader()
        CALL reset_test_failed_flag() ! Reset flag for current run

        DO i = 1, SIZE(test_names)
            found = .false.
            DO j = 1, SIZE(all_tests)
                IF (TRIM(test_names(i)) == TRIM(all_tests(j)%name)) THEN
                    WRITE(*,*) "Running test: ", TRIM(test_names(i)), "..."
                    CALL all_tests(j)%test_proc()
                    WRITE(*,*) TRIM(test_names(i)), " passed."
                    found = .true.
                    EXIT
                END IF
            END DO
            IF (.NOT. found) THEN
                WRITE(*,*) "WARNING: Unknown test: ", TRIM(test_names(i))
            END IF
        END DO
        IF (.NOT. test_failed_flag) THEN
            WRITE(*,*) "Specified csv_reader_module tests passed successfully."
        ELSE
            WRITE(*,*) "Some specified csv_reader_module tests FAILED."
            CALL error_stop_test() ! Re-throw error if any test failed
        END IF
    END SUBROUTINE run_named_tests_csv_reader

END MODULE mod_test_csv_reader