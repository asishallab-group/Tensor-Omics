! F42-Compliant CSV Reader and DataTable Module
! Implements binary serialization and deserialization.
! All memory cleanup is handled by a dedicated public subroutine.
MODULE csv_reader_module
    IMPLICIT NONE

    ! --- Public Constants ---
    INTEGER, PARAMETER, PUBLIC :: MAX_LINE_LEN = 4096
    INTEGER, PARAMETER, PUBLIC :: MAX_FIELD_LEN = 512
    INTEGER, PARAMETER, PUBLIC :: IK = SELECTED_INT_KIND(18)
    CHARACTER(LEN=1), PARAMETER, PUBLIC :: DEFAULT_DELIMITER = CHAR(9) ! TAB character

    ! --- Derived Type for public access, but not for internal storage ---
    TYPE, PUBLIC :: generic_data_cell
        INTEGER :: data_type = 0 ! 0=Empty, 1=Integer, 2=Real, 3=Character, 4=Logical, 5=Complex
        INTEGER(IK) :: i_val = 0
        REAL :: r_val = 0.0
        LOGICAL :: l_val = .FALSE.
        COMPLEX :: co_val = (0.0, 0.0)
        CHARACTER(LEN=:), ALLOCATABLE :: c_val
    END TYPE generic_data_cell

    ! --- NEW: Header type for binary serialization ---
    TYPE :: data_table_header
        CHARACTER(LEN=8) :: magic_number = 'F42DTBL'
        INTEGER :: version = 1
        INTEGER :: num_rows = 0
        INTEGER :: num_cols = 0
        INTEGER :: n_int_cols = 0
        INTEGER :: n_real_cols = 0
        INTEGER :: n_char_cols = 0
        INTEGER :: n_logical_cols = 0
        INTEGER :: n_complex_cols = 0
        INTEGER :: max_char_len = 0
    END TYPE data_table_header

    ! --- PRIVATE: Refactored internal data storage ---
    PRIVATE
    INTEGER :: num_rows = 0
    INTEGER :: num_cols = 0
    INTEGER :: n_int_cols = 0
    INTEGER :: n_real_cols = 0
    INTEGER :: n_char_cols = 0
    INTEGER :: n_logical_cols = 0
    INTEGER :: n_complex_cols = 0

    ! Typed data containers
    INTEGER(IK), ALLOCATABLE :: int_array(:,:)
    REAL, ALLOCATABLE :: real_array(:,:)
    LOGICAL, ALLOCATABLE :: logical_array(:,:)
    COMPLEX, ALLOCATABLE :: complex_array(:,:)
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: char_array(:,:)

    ! Metadata and mapping arrays
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: header_fields(:)
    INTEGER, ALLOCATABLE :: column_data_type(:)  ! 1=Int, 2=Real, 3=Char, 4=Logical, 5=Complex
    INTEGER, ALLOCATABLE :: column_typed_index(:) ! Index within the typed array

    ! --- Public Interface ---
    PUBLIC :: read_csv_file, cleanup_csv_data, serialize, deserialize
    PUBLIC :: get_num_rows, get_num_cols, get_cell, get_header
    PUBLIC :: get_column_data_type_by_index, get_column_data_type_by_name
    PUBLIC :: get_int_column, get_real_column, get_char_column, get_logical_column, get_complex_column

    ! --- Private procedures ---
    PRIVATE :: generate_header, parse_and_store_line, infer_type
    PRIVATE :: replace_char_all ! Declaring local REPLACE function

CONTAINS

    !=======================================================================
    ! Public Subroutine: read_csv_file (Refactored)
    ! Reads a CSV, determines types, and populates the typed arrays.
    !=======================================================================
    SUBROUTINE read_csv_file(filename, has_header, delimiter, io_status)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        LOGICAL, INTENT(IN) :: has_header
        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: delimiter
        INTEGER, INTENT(OUT) :: io_status
        INTEGER :: unit_num = 10, stat, line_count, i, j
        CHARACTER(LEN=MAX_LINE_LEN) :: line, first_data_line
        CHARACTER(LEN=1) :: current_delimiter
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
        INTEGER :: inferred_type
        CHARACTER(LEN=MAX_FIELD_LEN) :: temp_field_str ! Declare a temporary variable

        io_status = 0
        CALL cleanup_csv_data()

        IF (PRESENT(delimiter)) THEN
            current_delimiter = delimiter
        ELSE
            current_delimiter = DEFAULT_DELIMITER ! Set default to TAB
        END IF

        ! --- Pass 1: Determine dimensions and infer types ---
        OPEN(UNIT=unit_num, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=stat)
        IF (stat /= 0) THEN
            io_status = 1; RETURN
        END IF

        ! Read header if it exists
        IF (has_header) THEN
            READ(unit_num, '(A)', IOSTAT=stat) line
            IF (stat == 0) THEN
               CALL parse_and_store_line(line, current_delimiter, fields, stat)
               num_cols = SIZE(fields)
               ALLOCATE(header_fields(num_cols)); header_fields = fields
            END IF
        END IF

        ! Read first data line to determine types and column count
        first_data_line = ''
        DO
            READ(unit_num, '(A)', IOSTAT=stat) line
            IF (stat < 0) EXIT
            IF (LEN_TRIM(line) > 0) THEN
                first_data_line = line
                num_rows = num_rows + 1
                EXIT
            END IF
        END DO

        IF (first_data_line == '') THEN
            io_status = 4; CLOSE(unit_num); RETURN ! Empty file
        END IF

        IF (.NOT. has_header) THEN
            CALL parse_and_store_line(first_data_line, current_delimiter, fields, stat)
            num_cols = SIZE(fields)
        END IF

        ! Count remaining rows
        DO
            READ(unit_num, '(A)', IOSTAT=stat) line
            IF (stat < 0) EXIT
            IF (LEN_TRIM(line) > 0) num_rows = num_rows + 1
        END DO
        REWIND(unit_num)

        ! --- Allocate and set up mapping arrays ---
        ALLOCATE(column_data_type(num_cols), column_typed_index(num_cols))
        CALL parse_and_store_line(first_data_line, current_delimiter, fields, stat)
        n_int_cols = 0; n_real_cols = 0; n_char_cols = 0; n_logical_cols = 0; n_complex_cols = 0
        DO i = 1, num_cols
            CALL infer_type(TRIM(fields(i)), inferred_type)
            column_data_type(i) = inferred_type
            SELECT CASE(inferred_type)
            CASE(1)
                n_int_cols = n_int_cols + 1; column_typed_index(i) = n_int_cols
            CASE(2)
                n_real_cols = n_real_cols + 1; column_typed_index(i) = n_real_cols
            CASE(3)
                n_char_cols = n_char_cols + 1; column_typed_index(i) = n_char_cols
            CASE(4)
                n_logical_cols = n_logical_cols + 1; column_typed_index(i) = n_logical_cols
            CASE(5)
                n_complex_cols = n_complex_cols + 1; column_typed_index(i) = n_complex_cols
            END SELECT
        END DO
        DEALLOCATE(fields)

        ! --- Allocate typed data arrays ---
        IF (n_int_cols > 0) ALLOCATE(int_array(num_rows, n_int_cols))
        IF (n_real_cols > 0) ALLOCATE(real_array(num_rows, n_real_cols))
        IF (n_logical_cols > 0) ALLOCATE(logical_array(num_rows, n_logical_cols))
        IF (n_complex_cols > 0) ALLOCATE(complex_array(num_rows, n_complex_cols))
        IF (n_char_cols > 0) ALLOCATE(char_array(num_rows, n_char_cols))

        ! --- Pass 2: Read and store data ---
        IF (has_header) READ(unit_num, '(A)', IOSTAT=stat) ! Skip header
        DO i = 1, num_rows
            READ(unit_num, '(A)', IOSTAT=stat) line
            IF (stat < 0 .OR. LEN_TRIM(line) == 0) CYCLE
            CALL parse_and_store_line(line, current_delimiter, fields, stat)
            IF (SIZE(fields) /= num_cols) THEN; DEALLOCATE(fields); CYCLE; END IF ! Skip mismatched rows
            DO j = 1, num_cols
                temp_field_str = TRIM(fields(j)) ! Assign the trimmed field to the temporary variable
                SELECT CASE(column_data_type(j))
                CASE(1)
                    READ(temp_field_str, *, IOSTAT=stat) int_array(i, column_typed_index(j))
                    IF (stat /= 0) int_array(i, column_typed_index(j)) = 0
                CASE(2)
                    READ(temp_field_str, *, IOSTAT=stat) real_array(i, column_typed_index(j))
                    IF (stat /= 0) real_array(i, column_typed_index(j)) = 0.0
                CASE(3)
                    char_array(i, column_typed_index(j)) = temp_field_str
                CASE(4)
                    READ(temp_field_str, *, IOSTAT=stat) logical_array(i, column_typed_index(j))
                    IF (stat /= 0) logical_array(i, column_typed_index(j)) = .FALSE.
                CASE(5)
                    READ(temp_field_str, *, IOSTAT=stat) complex_array(i, column_typed_index(j))
                    IF (stat /= 0) complex_array(i, column_typed_index(j)) = (0.0, 0.0)
                END SELECT
            END DO
            DEALLOCATE(fields)
        END DO
        CLOSE(unit_num)
    END SUBROUTINE read_csv_file

    !=======================================================================
    ! Public Subroutine: serialize
    !=======================================================================
    SUBROUTINE serialize(filename, io_status)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(OUT) :: io_status
        INTEGER :: unit_num = 20
        TYPE(data_table_header) :: header

        io_status = 0
        header = generate_header()

        OPEN(UNIT=unit_num, FILE=TRIM(filename), FORM='UNFORMATTED', &
             ACCESS='STREAM', STATUS='REPLACE', IOSTAT=io_status)
        IF (io_status /= 0) RETURN

        WRITE(unit_num, IOSTAT=io_status) header
        IF (io_status /= 0) THEN; CLOSE(unit_num); RETURN; END IF

        IF (ALLOCATED(header_fields)) WRITE(unit_num, IOSTAT=io_status) header_fields
        WRITE(unit_num, IOSTAT=io_status) column_data_type
        WRITE(unit_num, IOSTAT=io_status) column_typed_index
        IF (header%n_int_cols > 0) WRITE(unit_num, IOSTAT=io_status) int_array
        IF (header%n_real_cols > 0) WRITE(unit_num, IOSTAT=io_status) real_array
        IF (header%n_logical_cols > 0) WRITE(unit_num, IOSTAT=io_status) logical_array
        IF (header%n_complex_cols > 0) WRITE(unit_num, IOSTAT=io_status) complex_array
        IF (header%n_char_cols > 0) WRITE(unit_num, IOSTAT=io_status) char_array

        CLOSE(unit_num)
    END SUBROUTINE serialize

    !=======================================================================
    ! Public Subroutine: deserialize
    !=======================================================================
    SUBROUTINE deserialize(filename, io_status)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(OUT) :: io_status
        INTEGER :: unit_num = 20
        TYPE(data_table_header) :: header

        io_status = 0
        CALL cleanup_csv_data()

        OPEN(UNIT=unit_num, FILE=TRIM(filename), FORM='UNFORMATTED', &
             ACCESS='STREAM', STATUS='OLD', IOSTAT=io_status)
        IF (io_status /= 0) RETURN

        READ(unit_num, IOSTAT=io_status) header
        IF (io_status /= 0) THEN; CLOSE(unit_num); RETURN; END IF

        ! Verify header
        IF (TRIM(header%magic_number) /= 'F42DTBL') THEN
            io_status = -1 ! Bad magic number
            CLOSE(unit_num)
            RETURN
        END IF

        ! Populate module variables from header
        num_rows = header%num_rows; num_cols = header%num_cols
        n_int_cols = header%n_int_cols; n_real_cols = header%n_real_cols
        n_char_cols = header%n_char_cols
        n_logical_cols = header%n_logical_cols
        n_complex_cols = header%n_complex_cols

        ! Allocate arrays based on header info
        ALLOCATE(header_fields(num_cols))
        ALLOCATE(column_data_type(num_cols))
        ALLOCATE(column_typed_index(num_cols))
        IF (n_int_cols > 0) ALLOCATE(int_array(num_rows, n_int_cols))
        IF (n_real_cols > 0) ALLOCATE(real_array(num_rows, n_real_cols))
        IF (n_logical_cols > 0) ALLOCATE(logical_array(num_rows, n_logical_cols))
        IF (n_complex_cols > 0) ALLOCATE(complex_array(num_rows, n_complex_cols))
        IF (n_char_cols > 0) ALLOCATE(char_array(num_rows, n_char_cols))

        ! Read raw data blocks
        IF (num_cols > 0) READ(unit_num, IOSTAT=io_status) header_fields
        READ(unit_num, IOSTAT=io_status) column_data_type
        READ(unit_num, IOSTAT=io_status) column_typed_index
        IF (n_int_cols > 0) READ(unit_num, IOSTAT=io_status) int_array
        IF (n_real_cols > 0) READ(unit_num, IOSTAT=io_status) real_array
        IF (n_logical_cols > 0) READ(unit_num, IOSTAT=io_status) logical_array
        IF (n_complex_cols > 0) READ(unit_num, IOSTAT=io_status) complex_array
        IF (n_char_cols > 0) READ(unit_num, IOSTAT=io_status) char_array

        CLOSE(unit_num)
    END SUBROUTINE deserialize

    !=======================================================================
    ! Public "Destructor" Subroutine for manual memory cleanup
    !=======================================================================
    SUBROUTINE cleanup_csv_data()
        IF (ALLOCATED(int_array)) DEALLOCATE(int_array)
        IF (ALLOCATED(real_array)) DEALLOCATE(real_array)
        IF (ALLOCATED(logical_array)) DEALLOCATE(logical_array)
        IF (ALLOCATED(complex_array)) DEALLOCATE(complex_array)
        IF (ALLOCATED(char_array)) DEALLOCATE(char_array)
        IF (ALLOCATED(header_fields)) DEALLOCATE(header_fields)
        IF (ALLOCATED(column_data_type)) DEALLOCATE(column_data_type)
        IF (ALLOCATED(column_typed_index)) DEALLOCATE(column_typed_index)
        num_rows=0; num_cols=0; n_int_cols=0; n_real_cols=0; n_char_cols=0
        n_logical_cols=0; n_complex_cols=0
    END SUBROUTINE cleanup_csv_data

    !=======================================================================
    ! Accessor Functions (Getters)
    !=======================================================================
    FUNCTION get_num_rows() RESULT(n); INTEGER :: n; n = num_rows; END FUNCTION
    FUNCTION get_num_cols() RESULT(n); INTEGER :: n; n = num_cols; END FUNCTION

    FUNCTION get_header(j) RESULT(header_name)
        INTEGER, INTENT(IN) :: j
        CHARACTER(LEN=:), ALLOCATABLE :: header_name ! Make return allocatable
        CHARACTER(LEN=MAX_FIELD_LEN) :: temp_header_name
        
        IF (ALLOCATED(header_fields) .AND. j > 0 .AND. j <= num_cols) THEN
            temp_header_name = header_fields(j)
            ALLOCATE(CHARACTER(LEN=LEN_TRIM(temp_header_name)) :: header_name)
            header_name = TRIM(temp_header_name)
        ELSE
            ALLOCATE(CHARACTER(LEN=0) :: header_name)
            header_name = ''
        END IF
    END FUNCTION get_header

    FUNCTION get_cell(i, j) RESULT(cell_data)
        INTEGER, INTENT(IN) :: i, j
        TYPE(generic_data_cell) :: cell_data
        INTEGER :: typed_idx
        ! Initialize cell_data to empty state
        cell_data%data_type = 0
        cell_data%i_val = 0
        cell_data%r_val = 0.0
        cell_data%l_val = .FALSE.
        cell_data%co_val = (0.0, 0.0)
        IF (ALLOCATED(cell_data%c_val)) DEALLOCATE(cell_data%c_val)

        IF (i > 0 .AND. i <= num_rows .AND. j > 0 .AND. j <= num_cols) THEN
            cell_data%data_type = column_data_type(j)
            typed_idx = column_typed_index(j)
            SELECT CASE(cell_data%data_type)
            CASE(1)
                cell_data%i_val = int_array(i, typed_idx)
            CASE(2)
                cell_data%r_val = real_array(i, typed_idx)
            CASE(3)
                ALLOCATE(CHARACTER(LEN=LEN_TRIM(char_array(i, typed_idx))) :: cell_data%c_val)
                cell_data%c_val = TRIM(char_array(i, typed_idx))
            CASE(4)
                cell_data%l_val = logical_array(i, typed_idx)
            CASE(5)
                cell_data%co_val = complex_array(i, typed_idx)
            END SELECT
        END IF
    END FUNCTION get_cell

    ! --- NEW: Getters for column data type ---
    FUNCTION get_column_data_type_by_index(j) RESULT(col_type)
        INTEGER, INTENT(IN) :: j
        INTEGER :: col_type
        col_type = 0 ! Default to empty/unknown
        IF (j > 0 .AND. j <= num_cols .AND. ALLOCATED(column_data_type)) THEN
            col_type = column_data_type(j)
        END IF
    END FUNCTION get_column_data_type_by_index

    FUNCTION get_column_data_type_by_name(col_name) RESULT(col_type)
        CHARACTER(LEN=*), INTENT(IN) :: col_name
        INTEGER :: col_type, i
        col_type = 0 ! Default to empty/unknown
        IF (ALLOCATED(header_fields)) THEN
            DO i = 1, num_cols
                IF (TRIM(header_fields(i)) == TRIM(col_name)) THEN
                    col_type = column_data_type(i)
                    RETURN
                END IF
            END DO
        END IF
    END FUNCTION get_column_data_type_by_name

    ! --- NEW: Getters for scalar column types ---
    FUNCTION get_int_column(j) RESULT(column_data)
        INTEGER, INTENT(IN) :: j
        INTEGER(IK), ALLOCATABLE :: column_data(:)
        INTEGER :: typed_idx

        IF (j > 0 .AND. j <= num_cols .AND. &
            column_data_type(j) == 1 .AND. ALLOCATED(int_array)) THEN
            typed_idx = column_typed_index(j)
            ALLOCATE(column_data(num_rows))
            column_data = int_array(:, typed_idx)
        ELSE
            ! Return unallocated array if conditions not met
            IF (ALLOCATED(column_data)) DEALLOCATE(column_data)
        END IF
    END FUNCTION get_int_column

    FUNCTION get_real_column(j) RESULT(column_data)
        INTEGER, INTENT(IN) :: j
        REAL, ALLOCATABLE :: column_data(:)
        INTEGER :: typed_idx

        IF (j > 0 .AND. j <= num_cols .AND. &
            column_data_type(j) == 2 .AND. ALLOCATED(real_array)) THEN
            typed_idx = column_typed_index(j)
            ALLOCATE(column_data(num_rows))
            column_data = real_array(:, typed_idx)
        ELSE
            IF (ALLOCATED(column_data)) DEALLOCATE(column_data)
        END IF
    END FUNCTION get_real_column

    FUNCTION get_char_column(j) RESULT(column_data)
        INTEGER, INTENT(IN) :: j
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: column_data(:)
        INTEGER :: typed_idx

        IF (j > 0 .AND. j <= num_cols .AND. &
            column_data_type(j) == 3 .AND. ALLOCATED(char_array)) THEN
            typed_idx = column_typed_index(j)
            ALLOCATE(column_data(num_rows))
            column_data = char_array(:, typed_idx)
        ELSE
            IF (ALLOCATED(column_data)) DEALLOCATE(column_data)
        END IF
    END FUNCTION get_char_column

    FUNCTION get_logical_column(j) RESULT(column_data)
        INTEGER, INTENT(IN) :: j
        LOGICAL, ALLOCATABLE :: column_data(:)
        INTEGER :: typed_idx

        IF (j > 0 .AND. j <= num_cols .AND. &
            column_data_type(j) == 4 .AND. ALLOCATED(logical_array)) THEN
            typed_idx = column_typed_index(j)
            ALLOCATE(column_data(num_rows))
            column_data = logical_array(:, typed_idx)
        ELSE
            IF (ALLOCATED(column_data)) DEALLOCATE(column_data)
        END IF
    END FUNCTION get_logical_column

    FUNCTION get_complex_column(j) RESULT(column_data)
        INTEGER, INTENT(IN) :: j
        COMPLEX, ALLOCATABLE :: column_data(:)
        INTEGER :: typed_idx

        IF (j > 0 .AND. j <= num_cols .AND. &
            column_data_type(j) == 5 .AND. ALLOCATED(complex_array)) THEN
            typed_idx = column_typed_index(j)
            ALLOCATE(column_data(num_rows))
            column_data = complex_array(:, typed_idx)
        ELSE
            IF (ALLOCATED(column_data)) DEALLOCATE(column_data)
        END IF
    END FUNCTION get_complex_column

    !=======================================================================
    ! Private Helper Subroutines
    !=======================================================================
    FUNCTION generate_header() RESULT(header)
        TYPE(data_table_header) :: header
        header%num_rows = num_rows
        header%num_cols = num_cols
        header%n_int_cols = n_int_cols
        header%n_real_cols = n_real_cols
        header%n_char_cols = n_char_cols
        header%n_logical_cols = n_logical_cols
        header%n_complex_cols = n_complex_cols
        header%max_char_len = MAX_FIELD_LEN
    END FUNCTION generate_header

    SUBROUTINE parse_and_store_line(line, delim, fields, status)
        CHARACTER(LEN=*), INTENT(IN) :: line
        CHARACTER(LEN=1), INTENT(IN) :: delim
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, INTENT(OUT) :: fields(:)
        INTEGER, INTENT(OUT) :: status
        INTEGER :: n_parts, i
        CHARACTER(LEN=LEN(line)) :: temp_line
        status = 0
        temp_line = line
        n_parts = 0
        IF(LEN_TRIM(line) > 0) THEN
            n_parts = 1
            DO i = 1, LEN_TRIM(line)
                IF (line(i:i) == delim) n_parts = n_parts + 1
            END DO
        END IF
        ALLOCATE(fields(n_parts))
        DO i = 1, n_parts -1
            fields(i) = temp_line(: SCAN(temp_line, delim) - 1)
            temp_line = temp_line(SCAN(temp_line, delim) + 1 :)
        END DO
        IF (n_parts > 0) fields(n_parts) = temp_line
    END SUBROUTINE parse_and_store_line

    ! New private helper to infer type more robustly
    SUBROUTINE infer_type(field, inferred_type)
        CHARACTER(LEN=*), INTENT(IN) :: field
        INTEGER, INTENT(OUT) :: inferred_type
        INTEGER :: stat
        REAL :: temp_r
        INTEGER(IK) :: temp_i
        LOGICAL :: temp_l
        COMPLEX :: temp_co
        CHARACTER(LEN=LEN(field)) :: stripped_field_tmp

        inferred_type = 3 ! Default to Character

        ! Remove all whitespace for robust parsing of complex and logical types
        stripped_field_tmp = ADJUSTL(ADJUSTR(field))
        stripped_field_tmp = replace_char_all(stripped_field_tmp, ' ', '') ! Use our local replace

        IF (LEN_TRIM(stripped_field_tmp) == 0) THEN
            RETURN ! Empty field, remains character
        END IF

        ! Try Logical first (e.g., .TRUE., T, .FALSE., F)
        READ(stripped_field_tmp, *, IOSTAT=stat) temp_l
        IF (stat == 0) THEN
            inferred_type = 4
            RETURN
        END IF

        ! Try Complex next (e.g., (1.0, 2.0))
        READ(stripped_field_tmp, *, IOSTAT=stat) temp_co
        IF (stat == 0) THEN
            inferred_type = 5
            RETURN
        END IF

        ! Try Integer (if no decimal point)
        IF (INDEX(stripped_field_tmp, '.') == 0) THEN
            READ(stripped_field_tmp, *, IOSTAT=stat) temp_i
            IF (stat == 0) THEN
                inferred_type = 1
                RETURN
            END IF
        END IF

        ! Try Real (if integer parsing failed or decimal point present)
        READ(stripped_field_tmp, *, IOSTAT=stat) temp_r
        IF (stat == 0) THEN
            inferred_type = 2
            RETURN
        END IF

        ! If all fail, it's a character (default)
    END SUBROUTINE infer_type

    ! --- NEW: Local REPLACE function ---
    ! Replaces all occurrences of 'old_char' with 'new_char' in 'string'
    FUNCTION replace_char_all(string, old_char, new_char) RESULT(res_string)
        CHARACTER(LEN=*), INTENT(IN) :: string
        CHARACTER(LEN=1), INTENT(IN) :: old_char, new_char
        CHARACTER(LEN=LEN(string)) :: res_string
        INTEGER :: i

        res_string = string
        DO i = 1, LEN(res_string)
            IF (res_string(i:i) == old_char) THEN
                res_string(i:i) = new_char
            END IF
        END DO
    END FUNCTION replace_char_all

END MODULE csv_reader_module