! F42-Compliant CSV Reader and DataTable Module
! Implements binary serialization and deserialization.
! All memory cleanup is handled by a dedicated public subroutine.
MODULE csv_reader_module
    IMPLICIT NONE

    ! --- Public Constants ---
    INTEGER, PARAMETER, PUBLIC :: MAX_LINE_LEN = 4096
    INTEGER, PARAMETER, PUBLIC :: MAX_FIELD_LEN = 512
    INTEGER, PARAMETER, PUBLIC :: IK = SELECTED_INT_KIND(18)

    ! --- Derived Type for public access, but not for internal storage ---
    TYPE, PUBLIC :: generic_data_cell
        INTEGER :: data_type = 0 ! 0=Empty, 1=Integer, 2=Real, 3=Character
        INTEGER(IK) :: i_val = 0
        REAL :: r_val = 0.0
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
        INTEGER :: max_char_len = 0
    END TYPE data_table_header

    ! --- PRIVATE: Refactored internal data storage ---
    PRIVATE
    INTEGER :: num_rows = 0
    INTEGER :: num_cols = 0
    INTEGER :: n_int_cols = 0
    INTEGER :: n_real_cols = 0
    INTEGER :: n_char_cols = 0

    ! Typed data containers
    INTEGER(IK), ALLOCATABLE :: int_array(:,:)
    REAL, ALLOCATABLE :: real_array(:,:)
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: char_array(:,:)

    ! Metadata and mapping arrays
    CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: header_fields(:)
    INTEGER, ALLOCATABLE :: column_data_type(:)  ! 1=Int, 2=Real, 3=Char
    INTEGER, ALLOCATABLE :: column_typed_index(:) ! Index within the typed array

    ! --- Public Interface ---
    PUBLIC :: read_csv_file, cleanup_csv_data, serialize, deserialize
    PUBLIC :: get_num_rows, get_num_cols, get_cell, get_header

    ! --- Private procedures ---
    PRIVATE :: generate_header, parse_and_store_line

CONTAINS

    !=======================================================================
    ! Public Subroutine: read_csv_file (Refactored)
    ! Reads a CSV, determines types, and populates the typed arrays.
    !=======================================================================
    SUBROUTINE read_csv_file(filename, has_header, delimiter, io_status)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        LOGICAL, INTENT(IN) :: has_header
        CHARACTER(LEN=1), INTENT(IN) :: delimiter
        INTEGER, INTENT(OUT) :: io_status
        INTEGER :: unit_num = 10, stat, line_count, i, j
        CHARACTER(LEN=MAX_LINE_LEN) :: line, first_data_line
        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE :: fields(:)
        REAL :: temp_r
        INTEGER(IK) :: temp_i
        CHARACTER(LEN=MAX_FIELD_LEN) :: temp_field

        io_status = 0
        CALL cleanup_csv_data()

        ! --- Pass 1: Determine dimensions and infer types ---
        OPEN(UNIT=unit_num, FILE=TRIM(filename), STATUS='OLD', ACTION='READ', IOSTAT=stat)
        IF (stat /= 0) THEN
            io_status = 1; RETURN
        END IF

        ! Read header if it exists
        IF (has_header) THEN
            READ(unit_num, '(A)', IOSTAT=stat) line
            IF (stat == 0) THEN
               CALL parse_and_store_line(line, delimiter, fields, stat)
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
             CALL parse_and_store_line(first_data_line, delimiter, fields, stat)
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
        CALL parse_and_store_line(first_data_line, delimiter, fields, stat)
        n_int_cols = 0; n_real_cols = 0; n_char_cols = 0
        DO i = 1, num_cols
            temp_field = TRIM(fields(i))
            IF (LEN_TRIM(temp_field) == 0) THEN ! Handle empty fields
                column_data_type(i) = 3; n_char_cols = n_char_cols + 1; column_typed_index(i) = n_char_cols
            ELSE IF (INDEX(temp_field, '.') > 0) THEN 
                 READ(temp_field, *, IOSTAT=stat) temp_r
                 IF (stat == 0) THEN
                    column_data_type(i) = 2; n_real_cols = n_real_cols + 1; column_typed_index(i) = n_real_cols
                 ELSE
                    column_data_type(i) = 3; n_char_cols = n_char_cols + 1; column_typed_index(i) = n_char_cols
                 END IF
            ELSE
                 READ(temp_field, *, IOSTAT=stat) temp_i
                 IF (stat == 0) THEN
                    column_data_type(i) = 1; n_int_cols = n_int_cols + 1; column_typed_index(i) = n_int_cols
                 ELSE
                    column_data_type(i) = 3; n_char_cols = n_char_cols + 1; column_typed_index(i) = n_char_cols
                 END IF
            END IF
        END DO
        DEALLOCATE(fields)

        ! --- Allocate typed data arrays ---
        IF (n_int_cols > 0) ALLOCATE(int_array(num_rows, n_int_cols))
        IF (n_real_cols > 0) ALLOCATE(real_array(num_rows, n_real_cols))
        IF (n_char_cols > 0) ALLOCATE(char_array(num_rows, n_char_cols))

        ! --- Pass 2: Read and store data ---
        IF (has_header) READ(unit_num, '(A)', IOSTAT=stat) ! Skip header
        DO i = 1, num_rows
            READ(unit_num, '(A)', IOSTAT=stat) line
            IF (stat < 0 .OR. LEN_TRIM(line) == 0) CYCLE
            CALL parse_and_store_line(line, delimiter, fields, stat)
            IF (SIZE(fields) /= num_cols) THEN; DEALLOCATE(fields); CYCLE; END IF ! Skip mismatched rows
            DO j = 1, num_cols
                temp_field = TRIM(fields(j))
                SELECT CASE(column_data_type(j))
                CASE(1)
                    READ(temp_field,*, IOSTAT=stat) int_array(i, column_typed_index(j))
                    IF (stat /= 0) int_array(i, column_typed_index(j)) = 0
                CASE(2)
                    READ(temp_field,*, IOSTAT=stat) real_array(i, column_typed_index(j))
                    IF (stat /= 0) real_array(i, column_typed_index(j)) = 0.0
                CASE(3)
                    char_array(i, column_typed_index(j)) = temp_field
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

        ! Allocate arrays based on header info
        ALLOCATE(header_fields(num_cols))
        ALLOCATE(column_data_type(num_cols))
        ALLOCATE(column_typed_index(num_cols))
        IF (n_int_cols > 0) ALLOCATE(int_array(num_rows, n_int_cols))
        IF (n_real_cols > 0) ALLOCATE(real_array(num_rows, n_real_cols))
        IF (n_char_cols > 0) ALLOCATE(char_array(num_rows, n_char_cols))

        ! Read raw data blocks
        IF (num_cols > 0) READ(unit_num, IOSTAT=io_status) header_fields
        READ(unit_num, IOSTAT=io_status) column_data_type
        READ(unit_num, IOSTAT=io_status) column_typed_index
        IF (n_int_cols > 0) READ(unit_num, IOSTAT=io_status) int_array
        IF (n_real_cols > 0) READ(unit_num, IOSTAT=io_status) real_array
        IF (n_char_cols > 0) READ(unit_num, IOSTAT=io_status) char_array

        CLOSE(unit_num)
    END SUBROUTINE deserialize

    !=======================================================================
    ! Public "Destructor" Subroutine for manual memory cleanup
    !=======================================================================
    SUBROUTINE cleanup_csv_data()
        IF (ALLOCATED(int_array)) DEALLOCATE(int_array)
        IF (ALLOCATED(real_array)) DEALLOCATE(real_array)
        IF (ALLOCATED(char_array)) DEALLOCATE(char_array)
        IF (ALLOCATED(header_fields)) DEALLOCATE(header_fields)
        IF (ALLOCATED(column_data_type)) DEALLOCATE(column_data_type)
        IF (ALLOCATED(column_typed_index)) DEALLOCATE(column_typed_index)
        num_rows=0; num_cols=0; n_int_cols=0; n_real_cols=0; n_char_cols=0
    END SUBROUTINE cleanup_csv_data

    !=======================================================================
    ! Accessor Functions (Getters)
    !=======================================================================
    FUNCTION get_num_rows() RESULT(n); INTEGER :: n; n = num_rows; END FUNCTION
    FUNCTION get_num_cols() RESULT(n); INTEGER :: n; n = num_cols; END FUNCTION

    FUNCTION get_header(j) RESULT(header_name)
        INTEGER, INTENT(IN) :: j
        CHARACTER(LEN=MAX_FIELD_LEN) :: header_name
        header_name = ''
        IF (ALLOCATED(header_fields) .AND. j > 0 .AND. j <= num_cols) THEN
            header_name = header_fields(j)
        END IF
    END FUNCTION get_header

    FUNCTION get_cell(i, j) RESULT(cell_data)
        INTEGER, INTENT(IN) :: i, j
        TYPE(generic_data_cell) :: cell_data
        INTEGER :: typed_idx
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
            END SELECT
        END IF
    END FUNCTION get_cell

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

END MODULE csv_reader_module