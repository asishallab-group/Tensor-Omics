! filepath: test/mod_test_csv_reader.f90
MODULE mod_test_csv_reader
    USE asserts, ONLY: assert_true, assert_false, assert_equal_int, assert_equal_real, &
                        assert_string_equal, assert_equal_array_real, assert_no_nan_real, &
                        assert_no_inf_real, assert_in_range_real, assert_contains_int, &
                        assert_sorted_int, assert_sorted_real, assert_same_shape, &
                        assert_string_contains, assert_allclose_array_real, assert_sum_equal, &
                        assert_unique_int, assert_permutation, &
                        test_failed_flag, reset_test_failed_flag, error_stop_test_wrapper

    ! Import only necessary public interfaces from your core module
    USE csv_reader_module, ONLY: read_csv_file, cleanup_csv_data, serialize, deserialize, &
                                 get_num_rows, get_num_cols, get_header, &
                                 get_column_data_type_by_index, get_column_data_type_by_name, &
                                 get_int_column, get_real_column, get_char_column, &
                                 get_logical_column, get_complex_column, get_cell, &
                                 IK, MAX_FIELD_LEN, DEFAULT_DELIMITER, generic_data_cell

    USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: int64, real64

    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = real64
    
    ABSTRACT INTERFACE
        SUBROUTINE test_interface()
        END SUBROUTINE test_interface
    END INTERFACE

    TYPE :: test_case
        CHARACTER(LEN=64) :: name
        PROCEDURE(test_interface), POINTER, NOPASS :: test_proc => NULL()
    END TYPE test_case

    PRIVATE

    PUBLIC :: get_all_tests_csv_reader, run_all_tests_csv_reader, run_named_tests_csv_reader

CONTAINS

    SUBROUTINE create_dummy_csv(filename, content)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=*), INTENT(IN) :: content
        INTEGER :: unit
        OPEN(NEWUNIT=unit, FILE=TRIM(filename), STATUS='REPLACE', ACTION='WRITE', FORM='FORMATTED')
        WRITE(unit, '(A)') TRIM(content)
        CLOSE(unit)
    END SUBROUTINE create_dummy_csv

    SUBROUTINE delete_dummy_csv(filename)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        LOGICAL :: file_exists
        INTEGER :: system_status
        INQUIRE(FILE=TRIM(filename), EXIST=file_exists)
        IF (file_exists) THEN
            CALL SYSTEM('rm ' // TRIM(filename), STATUS=system_status)
        END IF
    END SUBROUTINE delete_dummy_csv

    !> @brief Test basic CSV reading with header, including all data types.
    SUBROUTINE test_basic_read_with_header()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, num_rows, num_cols
        CHARACTER(LEN=MAX_FIELD_LEN) :: header_name_temp
        REAL(dp) :: r_val
        LOGICAL :: l_val
        COMPLEX(dp) :: co_val
        INTEGER(IK) :: i_val
        CHARACTER(LEN=MAX_FIELD_LEN) :: c_val_temp
        TYPE(generic_data_cell) :: temp_cell

        filename = "test_basic_header.csv"
        content = "ID,Name,Value,Active,Coord" // NEW_LINE('') // &
                  "1,Alpha,10.5,.TRUE.,(1.0,2.0)" // NEW_LINE('') // &
                  "2,Beta,20.0,.FALSE.,(3.0,-4.0)" // NEW_LINE('') // &
                  "3,Gamma,-5.2,T,(-0.5,0.5)" // NEW_LINE('') // &
                  "4,Delta,1.0,F,(0.0,0.0)"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, &
                              "test_basic_read_with_header: read status")
        
        num_rows = get_num_rows()
        num_cols = get_num_cols()
        CALL assert_equal_int(INT(num_rows, int64), INT(4, int64), &
                              "test_basic_read_with_header: num_rows")
        CALL assert_equal_int(INT(num_cols, int64), INT(5, int64), &
                              "test_basic_read_with_header: num_cols")

        header_name_temp = get_header(2)
        CALL assert_string_equal(TRIM(header_name_temp), "Name", "test_basic_read_with_header: header name")

        temp_cell = get_cell(1, 1)
        i_val = temp_cell%i_val
        CALL assert_equal_int(i_val, 1_int64, "test_basic_read_with_header: cell(1,1) int")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(1, 3)
        r_val = temp_cell%r_val
        CALL assert_equal_real(r_val, 10.5_dp, 1.0E-6_dp, "test_basic_read_with_header: cell(1,3) real")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(1,4)
        l_val = temp_cell%l_val
        CALL assert_true(l_val, "test_basic_read_with_header: cell(1,4) logical .TRUE.")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(1,5)
        co_val = temp_cell%co_val
        CALL assert_equal_real(REAL(co_val, kind=dp), 1.0_dp, 1.0E-6_dp, &
                               "test_basic_read_with_header: cell(1,5) complex real part")
        CALL assert_equal_real(AIMAG(co_val), 2.0_dp, 1.0E-6_dp, &
                               "test_basic_read_with_header: cell(1,5) complex imag part")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(2,2)
        c_val_temp = temp_cell%c_val
        CALL assert_string_equal(TRIM(c_val_temp), "Beta", "test_basic_read_with_header: cell(2,2) char")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(3,4)
        l_val = temp_cell%l_val
        CALL assert_true(l_val, "test_basic_read_with_header: cell(3,4) logical T")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(4,4)
        l_val = temp_cell%l_val
        CALL assert_false(l_val, "test_basic_read_with_header: cell(4,4) logical F")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_basic_read_with_header

    !> @brief Test CSV reading without header.
    SUBROUTINE test_read_no_header()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, num_rows, num_cols
        INTEGER(IK) :: i_val
        REAL(dp) :: r_val
        TYPE(generic_data_cell) :: temp_cell

        filename = "test_no_header.csv"
        content = "1,Alpha,10.5" // NEW_LINE('') // &
                  "2,Beta,20.0"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .FALSE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, &
                              "test_read_no_header: read status")
        
        num_rows = get_num_rows()
        num_cols = get_num_cols()
        CALL assert_equal_int(INT(num_rows, int64), INT(2, int64), "test_read_no_header: num_rows")
        CALL assert_equal_int(INT(num_cols, int64), INT(3, int64), "test_read_no_header: num_cols")

        temp_cell = get_cell(1, 1)
        i_val = temp_cell%i_val
        CALL assert_equal_int(i_val, 1_int64, "test_read_no_header: cell(1,1) int")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)
        
        temp_cell = get_cell(1, 3)
        r_val = temp_cell%r_val
        CALL assert_equal_real(r_val, 10.5_dp, 1.0E-6_dp, "test_read_no_header: cell(1,3) real")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_read_no_header

    !> @brief Test CSV reading with TAB delimiter (default).
    SUBROUTINE test_read_tab_delimiter()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        CHARACTER(LEN=MAX_FIELD_LEN) :: char_val
        TYPE(generic_data_cell) :: temp_cell
        
        filename = "test_tab_delim.csv"
        content = "ID" // DEFAULT_DELIMITER // "Name" // NEW_LINE('') // &
                  "1" // DEFAULT_DELIMITER // "Alpha"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=DEFAULT_DELIMITER, io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_read_tab_delimiter: read status with explicit TAB")
        CALL assert_equal_int(INT(get_num_rows(), int64), INT(1, int64), "test_read_tab_delimiter: num_rows explicit TAB")
        CALL assert_equal_int(INT(get_num_cols(), int64), INT(2, int64), "test_read_tab_delimiter: num_cols explicit TAB")
        temp_cell = get_cell(1, 2)
        char_val = temp_cell%c_val
        CALL assert_string_equal(TRIM(char_val), "Alpha", "test_read_tab_delimiter: cell(1,2) char explicit TAB")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)
        CALL cleanup_csv_data()

        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_read_tab_delimiter: read status with default TAB")
        CALL assert_equal_int(INT(get_num_rows(), int64), INT(1, int64), "test_read_tab_delimiter: num_rows default TAB")
        CALL assert_equal_int(INT(get_num_cols(), int64), INT(2, int64), "test_read_tab_delimiter: num_cols default TAB")
        temp_cell = get_cell(1, 2)
        char_val = temp_cell%c_val
        CALL assert_string_equal(TRIM(char_val), "Alpha", "test_read_tab_delimiter: cell(1,2) char default TAB")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)
        CALL cleanup_csv_data()

        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_read_tab_delimiter

    !> @brief Test error handling for non-existent file.
    SUBROUTINE test_read_file_not_found()
        CHARACTER(LEN=64) :: filename
        INTEGER :: io_status

        filename = "non_existent_file_xyz.csv"
        CALL read_csv_file(filename, .TRUE., io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), INT(1, int64), "test_read_file_not_found: io_status for missing file")
        CALL cleanup_csv_data()
    END SUBROUTINE test_read_file_not_found

    !> @brief Test parsing with empty fields.
    SUBROUTINE test_empty_fields()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        CHARACTER(LEN=MAX_FIELD_LEN) :: char_val
        COMPLEX(dp) :: comp_val
        TYPE(generic_data_cell) :: temp_cell

        filename = "test_empty_fields.csv"
        content = "Col1,Col2,Col3,Col4,Col5" // NEW_LINE('') // &
                  "1,,Text,,( , )" ! Empty number and complex part
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_empty_fields: read status")
        CALL assert_equal_int(INT(get_num_rows(), int64), INT(1, int64), "test_empty_fields: num_rows")
        CALL assert_equal_int(INT(get_num_cols(), int64), INT(5, int64), "test_empty_fields: num_cols")

        temp_cell = get_cell(1, 2)
        char_val = temp_cell%c_val
        CALL assert_string_equal(TRIM(char_val), "", "test_empty_fields: Col2 empty string")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)
        ! CORRECTED: Line broken up
        CALL assert_equal_int(INT(get_column_data_type_by_index(2), int64), INT(3, int64), &
                              "test_empty_fields: Col2 type should be Character")
        
        temp_cell = get_cell(1, 3)
        char_val = temp_cell%c_val
        CALL assert_string_equal(TRIM(char_val), "Text", "test_empty_fields: Col3 char")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(1, 4)
        char_val = temp_cell%c_val
        CALL assert_string_equal(TRIM(char_val), "", "test_empty_fields: Col4 empty string")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)
        ! CORRECTED: Line broken up
        CALL assert_equal_int(INT(get_column_data_type_by_index(4), int64), INT(3, int64), &
                              "test_empty_fields: Col4 type should be Character")

        temp_cell = get_cell(1, 5)
        comp_val = temp_cell%co_val
        CALL assert_equal_real(REAL(comp_val, kind=dp), 0.0_dp, 1E-6_dp, "test_empty_fields: Col5 complex real default")
        CALL assert_equal_real(AIMAG(comp_val), 0.0_dp, 1E-6_dp, "test_empty_fields: Col5 complex imag default")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)
        ! CORRECTED: Line broken up
        CALL assert_equal_int(INT(get_column_data_type_by_index(5), int64), INT(5, int64), &
                              "test_empty_fields: Col5 type should be Complex")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_empty_fields

    !> @brief Test get_column_data_type_by_index and get_column_data_type_by_name.
    SUBROUTINE test_get_column_types()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status, col_type

        filename = "test_col_types.csv"
        content = "IntCol,RealCol,CharCol,LogicalCol,ComplexCol" // NEW_LINE('') // &
                  "1,1.0,A,.TRUE.,(1.0,1.0)"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_get_column_types: read status")
        
        col_type = get_column_data_type_by_index(1)
        CALL assert_equal_int(INT(col_type, int64), INT(1, int64), "test_get_column_types: IntCol type by index")
        
        col_type = get_column_data_type_by_name("RealCol")
        CALL assert_equal_int(INT(col_type, int64), INT(2, int64), "test_get_column_types: RealCol type by name")

        col_type = get_column_data_type_by_name("CharCol")
        CALL assert_equal_int(INT(col_type, int64), INT(3, int64), "test_get_column_types: CharCol type by name")

        col_type = get_column_data_type_by_name("LogicalCol")
        CALL assert_equal_int(INT(col_type, int64), INT(4, int64), "test_get_column_types: LogicalCol type by name")
        
        col_type = get_column_data_type_by_name("ComplexCol")
        CALL assert_equal_int(INT(col_type, int64), INT(5, int64), "test_get_column_types: ComplexCol type by name")

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_get_column_types

    !> @brief Test column getters.
    SUBROUTINE test_column_getters()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        INTEGER(IK), ALLOCATABLE :: int_col(:)
        REAL(dp), ALLOCATABLE :: real_col(:)
        LOGICAL, ALLOCATABLE :: logical_col(:)
        COMPLEX(dp), ALLOCATABLE :: complex_col(:)
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: char_col(:)

        filename = "test_column_getters.csv"
        content = "Ints,Reals,Chars,Logicals,Complexes" // NEW_LINE('') // &
                  "1,1.1,A,.TRUE.,(1.0,1.0)" // NEW_LINE('') // &
                  "2,2.2,B,.FALSE.,(2.0,2.0)"
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_column_getters: read status")

        int_col = get_int_column(1)
        CALL assert_equal_int(INT(SIZE(int_col), int64), INT(2, int64), "test_column_getters: int_col size")
        CALL assert_equal_int(int_col(1), 1_int64, "test_column_getters: int_col(1)")
        CALL assert_equal_int(int_col(2), 2_int64, "test_column_getters: int_col(2)")
        IF (ALLOCATED(int_col)) DEALLOCATE(int_col)

        real_col = get_real_column(2)
        CALL assert_equal_int(INT(SIZE(real_col), int64), INT(2, int64), "test_column_getters: real_col size")
        CALL assert_equal_real(real_col(1), 1.1_dp, 1E-6_dp, "test_column_getters: real_col(1)")
        CALL assert_equal_real(real_col(2), 2.2_dp, 1E-6_dp, "test_column_getters: real_col(2)")
        IF (ALLOCATED(real_col)) DEALLOCATE(real_col)

        char_col = get_char_column(3)
        CALL assert_equal_int(INT(SIZE(char_col), int64), INT(2, int64), "test_column_getters: char_col size")
        CALL assert_string_equal(TRIM(char_col(1)), "A", "test_column_getters: char_col(1)")
        CALL assert_string_equal(TRIM(char_col(2)), "B", "test_column_getters: char_col(2)")
        IF (ALLOCATED(char_col)) DEALLOCATE(char_col)

        logical_col = get_logical_column(4)
        CALL assert_equal_int(INT(SIZE(logical_col), int64), INT(2, int64), "test_column_getters: logical_col size")
        CALL assert_true(logical_col(1), "test_column_getters: logical_col(1)")
        CALL assert_false(logical_col(2), "test_column_getters: logical_col(2)")
        IF (ALLOCATED(logical_col)) DEALLOCATE(logical_col)

        complex_col = get_complex_column(5)
        CALL assert_equal_int(INT(SIZE(complex_col), int64), INT(2, int64), "test_column_getters: complex_col size")
        CALL assert_equal_real(REAL(complex_col(1), kind=dp), 1.0_dp, 1E-6_dp, &
                               "test_column_getters: complex_col(1) real")
        CALL assert_equal_real(AIMAG(complex_col(1)), 1.0_dp, 1E-6_dp, &
                               "test_column_getters: complex_col(1) imag")
        CALL assert_equal_real(REAL(complex_col(2), kind=dp), 2.0_dp, 1E-6_dp, &
                               "test_column_getters: complex_col(2) real")
        CALL assert_equal_real(AIMAG(complex_col(2)), 2.0_dp, 1E-6_dp, &
                               "test_column_getters: complex_col(2) imag")
        IF (ALLOCATED(complex_col)) DEALLOCATE(complex_col)

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_column_getters

    !> @brief Test complex number parsing with whitespace.
    SUBROUTINE test_complex_whitespace_parsing()
        CHARACTER(LEN=64) :: filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        COMPLEX(dp) :: co_val
        TYPE(generic_data_cell) :: temp_cell

        filename = "test_complex_ws.csv"
        content = "Coord" // NEW_LINE('') // &
                  "( 1.0 , 2.0 )" // NEW_LINE('') // &
                  " ( 3.0 , -4.0 ) " // NEW_LINE('') // &
                  "(+5.0,+6.0)" ! Test explicit plus signs
        CALL create_dummy_csv(filename, content)

        CALL read_csv_file(filename, .TRUE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_complex_whitespace_parsing: read status")
        CALL assert_equal_int(INT(get_num_rows(), int64), INT(3, int64), "test_complex_whitespace_parsing: num_rows")
        
        temp_cell = get_cell(1, 1)
        co_val = temp_cell%co_val
        CALL assert_equal_real(REAL(co_val, kind=dp), 1.0_dp, 1.0E-6_dp, &
                               "test_complex_whitespace_parsing: cell(1,1) real part")
        CALL assert_equal_real(AIMAG(co_val), 2.0_dp, 1.0E-6_dp, &
                               "test_complex_whitespace_parsing: cell(1,1) imag part")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(2, 1)
        co_val = temp_cell%co_val
        CALL assert_equal_real(REAL(co_val, kind=dp), 3.0_dp, 1.0E-6_dp, &
                               "test_complex_whitespace_parsing: cell(2,1) real part")
        CALL assert_equal_real(AIMAG(co_val), -4.0_dp, 1.0E-6_dp, &
                               "test_complex_whitespace_parsing: cell(2,1) imag part")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(3, 1)
        co_val = temp_cell%co_val
        CALL assert_equal_real(REAL(co_val, kind=dp), 5.0_dp, 1.0E-6_dp, &
                               "test_complex_whitespace_parsing: cell(3,1) real part")
        CALL assert_equal_real(AIMAG(co_val), 6.0_dp, 1.0E-6_dp, &
                               "test_complex_whitespace_parsing: cell(3,1) imag part")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(filename)
    END SUBROUTINE test_complex_whitespace_parsing

    !> @brief Test serialization and deserialization.
    SUBROUTINE test_serialize_deserialize()
        CHARACTER(LEN=64) :: csv_filename, bin_filename
        CHARACTER(LEN=200) :: content
        INTEGER :: io_status
        INTEGER(IK) :: i_val_check
        REAL(dp) :: r_val_check
        LOGICAL :: l_val_check
        COMPLEX(dp) :: co_val_check
        CHARACTER(LEN=MAX_FIELD_LEN) :: c_val_check
        TYPE(generic_data_cell) :: temp_cell

        csv_filename = "test_serial.csv"
        bin_filename = "test_serial.bin"
        content = "Ints,Reals,Chars,Logicals,Complexes" // NEW_LINE('') // &
                  "1,1.1,A,.TRUE.,(1.0,1.0)" // NEW_LINE('') // &
                  "2,2.2,B,.FALSE.,(2.0,2.0)"
        CALL create_dummy_csv(filename=csv_filename, content=content)

        CALL read_csv_file(filename=csv_filename, has_header=.TRUE., delimiter=',', io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_serialize_deserialize: read status")
        CALL assert_equal_int(INT(get_num_rows(), int64), INT(2, int64), "test_serialize_deserialize: num_rows before serial")
        
        CALL serialize(filename=bin_filename, io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_serialize_deserialize: serialize status")
        
        CALL cleanup_csv_data()
        CALL assert_equal_int(INT(get_num_rows(), int64), 0_int64, "test_serialize_deserialize: num_rows after cleanup")

        CALL deserialize(filename=bin_filename, io_status=io_status)
        CALL assert_equal_int(INT(io_status, int64), 0_int64, "test_serialize_deserialize: deserialize status")
        
        CALL assert_equal_int(INT(get_num_rows(), int64), INT(2, int64), "test_serialize_deserialize: num_rows after deserial")
        CALL assert_equal_int(INT(get_num_cols(), int64), INT(5, int64), "test_serialize_deserialize: num_cols after deserial")
        
        temp_cell = get_cell(1, 1)
        i_val_check = temp_cell%i_val
        CALL assert_equal_int(i_val_check, 1_int64, "test_serialize_deserialize: cell(1,1) after deserial")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(2, 2)
        r_val_check = temp_cell%r_val
        CALL assert_equal_real(r_val_check, 2.2_dp, 1E-6_dp, "test_serialize_deserialize: cell(2,2) after deserial")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(1, 3)
        c_val_check = temp_cell%c_val
        CALL assert_string_equal(TRIM(c_val_check), "A", "test_serialize_deserialize: cell(1,3) after deserial")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(2, 4)
        l_val_check = temp_cell%l_val
        CALL assert_false(l_val_check, "test_serialize_deserialize: cell(2,4) after deserial")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        temp_cell = get_cell(1, 5)
        co_val_check = temp_cell%co_val
        CALL assert_equal_real(REAL(co_val_check, kind=dp), 1.0_dp, 1E-6_dp, &
                               "test_serialize_deserialize: cell(1,5) real after deserial")
        CALL assert_equal_real(AIMAG(co_val_check), 1.0_dp, 1E-6_dp, &
                               "test_serialize_deserialize: cell(1,5) imag after deserial")
        IF (ALLOCATED(temp_cell%c_val)) DEALLOCATE(temp_cell%c_val)

        CALL cleanup_csv_data()
        CALL delete_dummy_csv(csv_filename)
        CALL delete_dummy_csv(bin_filename)
    END SUBROUTINE test_serialize_deserialize

    FUNCTION get_all_tests_csv_reader() RESULT(all_tests)
        TYPE(test_case) :: all_tests(9)

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

    SUBROUTINE run_all_tests_csv_reader()
        TYPE(test_case) :: all_tests(9)
        INTEGER :: i

        all_tests = get_all_tests_csv_reader()
        CALL reset_test_failed_flag()

        DO i = 1, SIZE(all_tests)
            WRITE(*,*) "Running test: ", TRIM(all_tests(i)%name), "..."
            CALL all_tests(i)%test_proc()
            WRITE(*,*) TRIM(all_tests(i)%name), " passed."
        END DO
        IF (.NOT. test_failed_flag) THEN
            WRITE(*,*) "All csv_reader_module tests passed successfully."
        ELSE
            WRITE(*,*) "Some csv_reader_module tests FAILED."
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE run_all_tests_csv_reader

    SUBROUTINE run_named_tests_csv_reader(test_names)
        CHARACTER(LEN=*), INTENT(IN) :: test_names(:)
        TYPE(test_case) :: all_tests(9)
        INTEGER :: i, j
        LOGICAL :: found

        all_tests = get_all_tests_csv_reader()
        CALL reset_test_failed_flag()

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
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE run_named_tests_csv_reader

END MODULE mod_test_csv_reader
