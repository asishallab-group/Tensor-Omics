! =====================================================================
! Module: tox_datatable_mod
!compile with:
!   'gfortran tox_demo.f90 -o tox_demo'
! =====================================================================
MODULE tox_datatable_mod
    IMPLICIT NONE
    PRIVATE 

    ! public
    INTEGER, PARAMETER, PUBLIC :: TYPE_INT = 1
    INTEGER, PARAMETER, PUBLIC :: TYPE_REAL = 2
    INTEGER, PARAMETER, PUBLIC :: TYPE_CHAR = 3
    INTEGER, PARAMETER, PUBLIC :: COL_NAME_LEN = 256 ! max length for column names
    INTEGER, PARAMETER, PUBLIC :: IK = SELECTED_INT_KIND(18) ! memory calculation value

    ! parameters for serialisation
    CHARACTER(LEN=8), PARAMETER :: FILE_MAGIC_NUMBER = "F42DTBL " ! 42 bring a towel, why didnt the eagles just bring the ring to mordor?
    INTEGER, PARAMETER :: FILE_VERSION = 1

    ! derived type for file header
    TYPE :: tox_datatable_header
        CHARACTER(LEN=8) :: magic = FILE_MAGIC_NUMBER
        INTEGER :: version = FILE_VERSION
        INTEGER :: header_size = 0 ! store header if skipped
        INTEGER :: n_rows = 0
        INTEGER :: k_cols = 0
        INTEGER :: i_cols = 0
        INTEGER :: r_cols = 0
        INTEGER :: c_cols = 0
        INTEGER :: max_c_len = 0
        INTEGER(KIND=IK) :: offset_type_map = 0 ! byte offset from start
        INTEGER(KIND=IK) :: offset_int_data = 0
        INTEGER(KIND=IK) :: offset_real_data = 0
        INTEGER(KIND=IK) :: offset_char_data = 0
        INTEGER(KIND=IK) :: offset_col_names = 0
        INTEGER(KIND=IK) :: total_file_size = 0
    END TYPE tox_datatable_header

    ! --- Public Derived Type ---
    TYPE, PUBLIC :: tox_datatable
        INTEGER, ALLOCATABLE :: type_mapping_array(:,:) ! 2 x k (k=num_cols)
        INTEGER, ALLOCATABLE :: int_data(:,:)           ! n_rows x i (i=num_int_cols)
        REAL, ALLOCATABLE    :: real_data(:,:)          ! n_rows x r (r=num_real_cols)
        CHARACTER(LEN=:), ALLOCATABLE :: char_data(:,:) ! n_rows x c (c=num_char_cols), length=max_char_len
        CHARACTER(LEN=COL_NAME_LEN), ALLOCATABLE :: column_names(:) ! k

        ! --- Metadata (stored after reading/allocation) ---
        INTEGER :: num_rows = 0
        INTEGER :: num_cols = 0
        INTEGER :: num_int_cols = 0
        INTEGER :: num_real_cols = 0
        INTEGER :: num_char_cols = 0
        INTEGER :: max_char_len_used = 0
    CONTAINS
        PROCEDURE, PUBLIC :: read_from_file => datatable_read_from_file
        PROCEDURE, PUBLIC :: calculate_memory => datatable_calculate_memory
        PROCEDURE, PUBLIC :: destroy => datatable_destroy
        PROCEDURE, PUBLIC :: print_summary => datatable_print_summary
        PROCEDURE, PUBLIC :: print_data => datatable_print_data
        PROCEDURE, PUBLIC :: serialize => datatable_serialize      ! CAREFUL NOT TESTED IN DEPTH
        PROCEDURE, PUBLIC :: deserialize => datatable_deserialize  ! CAREFUL NOT TESTED IN DEPTH
    END TYPE tox_datatable

CONTAINS

    ! ==
    ! subroutine: datatable_read_from_file
    ! ==
    SUBROUTINE datatable_read_from_file(self, filename, sep_char, column_types, &
                                        max_char_len, n_rows, opt_column_names)
        CLASS(tox_datatable), INTENT(INOUT) :: self
        CHARACTER(LEN=*), INTENT(IN)       :: filename
        CHARACTER(LEN=1), INTENT(IN)       :: sep_char
        INTEGER, INTENT(IN)                :: column_types(:) ! k: 1=Int, 2=Real, 3=Char
        INTEGER, INTENT(IN)                :: max_char_len    ! max length for char columns
        INTEGER, INTENT(IN)                :: n_rows          ! expected number of data rows
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: opt_column_names(:) ! k: Optional provided names

        ! --- lcal variable declarations ---
        INTEGER :: k                  ! total number of columns
        INTEGER :: i_count, r_count, c_count ! Counts of each type
        INTEGER :: i, j, current_row
        INTEGER :: unit_num = 11      ! file unit number
        INTEGER :: io_stat, alloc_stat, parse_stat
        CHARACTER(LEN=4096) :: line   ! reading line buffer
        CHARACTER(LEN=:), ALLOCATABLE :: fields(:) ! For split line parts
        LOGICAL :: use_provided_names
        LOGICAL :: read_header_from_file
        CHARACTER(LEN=COL_NAME_LEN) :: default_name_fmt = '(A,I0)' ! Format for default names like "col1"
        INTEGER :: target_idx
        INTEGER :: col_type
        CHARACTER(LEN=:), ALLOCATABLE :: trimmed_field ! allocatable for flexibility
        INTEGER :: alloc_char_len ! char length

        ! --- executable Statments -

        ! validation stuff
        k = SIZE(column_types)
        IF (k <= 0) THEN
            WRITE(*,*) 'Error (read_from_file): column_types array must not be empty.'
            STOP 101
        END IF
        IF (n_rows <= 0) THEN
            WRITE(*,*) 'Error (read_from_file): n_rows must be positive.'
            STOP 102
        END IF
        ! Allow max_char_len to be 0 if there are no character columns
        IF (max_char_len <= 0 .AND. COUNT(column_types == TYPE_CHAR) > 0) THEN
            WRITE(*,*) 'Error (read_from_file): max_char_len must be positive if character columns exist.'
            STOP 103
        ELSE IF (max_char_len < 0) THEN ! no negative max_char_length, fix -> i have no idea how this even happened
             WRITE(*,*) 'Error (read_from_file): max_char_len cannot be negative.'
             STOP 103
        END IF
        DO i = 1, k
            IF (column_types(i) < TYPE_INT .OR. column_types(i) > TYPE_CHAR) THEN
                WRITE(*,'(A,I0,A,I0)') 'Error (read_from_file): Invalid column type ', &
                                       column_types(i), ' at index ', i
                STOP 104
            END IF
        END DO

        use_provided_names = PRESENT(opt_column_names)
        read_header_from_file = .NOT. use_provided_names

        IF (use_provided_names) THEN
            IF (SIZE(opt_column_names) /= k) THEN
                WRITE(*,*) 'Error (read_from_file): Size of opt_column_names must match size of column_types.'
                STOP 105
            END IF
        END IF

        ! --- Deallocate existing data if any ---
        CALL self%destroy()

        ! --- Calculate type counts and store metadata ---
        self%num_rows = n_rows
        self%num_cols = k
        ! Handle case where max_char_len is 0 because no char columns exist
        IF (COUNT(column_types == TYPE_CHAR) > 0) THEN
           self%max_char_len_used = max_char_len
        ELSE
           self%max_char_len_used = 0
        END IF
        i_count = COUNT(column_types == TYPE_INT)
        r_count = COUNT(column_types == TYPE_REAL)
        c_count = COUNT(column_types == TYPE_CHAR)
        self%num_int_cols = i_count
        self%num_real_cols = r_count
        self%num_char_cols = c_count

        WRITE(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') 'Info: k=', k, ', i=', i_count, &
              ', r=', r_count, ', c=', c_count, ', max_c_len=', self%max_char_len_used

        ! --- Allocate Internal Arrays ---
        ALLOCATE(self%type_mapping_array(2, k), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 901
        IF (i_count > 0) THEN
            ALLOCATE(self%int_data(n_rows, i_count), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 901
        END IF
        IF (r_count > 0) THEN
            ALLOCATE(self%real_data(n_rows, r_count), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 901
        END IF
        IF (c_count > 0) THEN
             ! Ensure max_char_len is at least 1 if c_count > 0 for allocation syntax
            alloc_char_len = MAX(1, self%max_char_len_used) ! Use 1 if actual len is 0 but c_count > 0 though this shouldnt happen with validation)
            ALLOCATE(CHARACTER(LEN=alloc_char_len) :: self%char_data(n_rows, c_count), STAT=alloc_stat)
            IF (alloc_stat /= 0) GOTO 901
        END IF
        ALLOCATE(self%column_names(k), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 901

        ! --- Populate Type Mapping Array ---
        i_count = 0; r_count = 0; c_count = 0 ! counter reset for indexing
        DO j = 1, k
            self%type_mapping_array(1, j) = column_types(j)
            SELECT CASE (column_types(j))
            CASE (TYPE_INT)
                i_count = i_count + 1
                self%type_mapping_array(2, j) = i_count
            CASE (TYPE_REAL)
                r_count = r_count + 1
                self%type_mapping_array(2, j) = r_count
            CASE (TYPE_CHAR)
                c_count = c_count + 1
                self%type_mapping_array(2, j) = c_count
            END SELECT
        END DO

        ! --- Open File ---
        OPEN(UNIT=unit_num, FILE=filename, STATUS='OLD', ACTION='READ', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'Error (read_from_file): Cannot open file "', TRIM(filename), '". IOSTAT=', io_stat
            CALL self%destroy() ! clean up allocated memory
            STOP 106
        END IF

        ! --- Handle Column Names ---
        IF (use_provided_names) THEN
            WRITE(*,*) 'Info: Using provided column names.'
            DO j = 1, k
                 self%column_names(j) = opt_column_names(j) ! this relies on intrinsic assignment truncation/padding - careful
            END DO
            ! Decide if file also has a header to skip (so far assumes no skip if names provided)
            ! 'if (skip_header_flag) read(unit_num, '(A)') line'
        ELSE IF (read_header_from_file) THEN
            WRITE(*,*) 'Info: Reading column names from file header.'
            READ(unit_num, '(A)', IOSTAT=io_stat) line
            IF (io_stat /= 0) THEN
                WRITE(*,'(A,A,I0)') 'Error (read_from_file): Cannot read header line from file "', &
                                    TRIM(filename), '". IOSTAT=', io_stat
                CLOSE(unit_num); CALL self%destroy(); STOP 107
            END IF
            CALL split_line_basic(TRIM(line), sep_char, fields, io_stat)
            IF (io_stat /= 0 .OR. SIZE(fields) /= k) THEN
                 WRITE(*,'(A,I0,A,I0)') 'Error: Header line has ', SIZE(fields), &
                       ' fields, expected ', k, '. Cannot parse header.'
                 IF (ALLOCATED(fields)) DEALLOCATE(fields)
                 CLOSE(unit_num); CALL self%destroy(); STOP 108
            END IF
            DO j = 1, k
                self%column_names(j) = TRIM(fields(j)) ! assign header field, rely on padding/truncation
            END DO
            DEALLOCATE(fields)
        ELSE ! generate default names (no header in file, names not provided)
            WRITE(*,*) 'Info: Generating default column names (col1, col2, ...).'
             DO j = 1, k
                 WRITE(self%column_names(j), fmt=default_name_fmt) 'col', j
             END DO
        END IF

        ! --- Read Data Rows ---
        WRITE(*,*) 'Info: Reading data rows...'
        DO current_row = 1, n_rows
            READ(unit_num, '(A)', IOSTAT=io_stat) line
            IF (io_stat < 0) THEN ! EOF
                WRITE(*,'(A,I0,A,I0)') 'Error: Unexpected end of file reached at row ', &
                                       current_row, '. Expected ', n_rows, ' data rows.'
                 CLOSE(unit_num); CALL self%destroy(); STOP 109
            ELSE IF (io_stat > 0) THEN
                 WRITE(*,'(A,I0,A,I0)') 'Error: Failed reading data line ', current_row, &
                                        '. IOSTAT=', io_stat
                 CLOSE(unit_num); CALL self%destroy(); STOP 110
            END IF

            ! handle blank lines - Fill row with defaults and continue
            IF (LEN_TRIM(line) == 0) THEN
                WRITE(*,'(A,I0)') 'Warning: Skipping blank line encountered at effective data row ', current_row
                ! Fill this row with defaults
                DO j = 1, k
                    col_type = self%type_mapping_array(1, j)
                    target_idx = self%type_mapping_array(2, j)
                    SELECT CASE(col_type)
                    CASE (TYPE_INT);  IF (ALLOCATED(self%int_data))  self%int_data(current_row, target_idx) = 0
                    CASE (TYPE_REAL); IF (ALLOCATED(self%real_data)) self%real_data(current_row, target_idx) = 0.0
                    CASE (TYPE_CHAR); IF (ALLOCATED(self%char_data)) self%char_data(current_row, target_idx) = ' '
                    END SELECT
                END DO
                CYCLE ! and next row
            END IF

            CALL split_line_basic(TRIM(line), sep_char, fields, io_stat)
            ! Hl andle rows with incorrect number of fields - Filrow with defaults and continue
            IF (io_stat /= 0 .OR. SIZE(fields) /= k) THEN
                 WRITE(*,'(A,I0,A,I0,A,I0)') 'Warning: Data row ', current_row, ' has ', SIZE(fields), &
                                             ' fields, expected ', k, '. Filling row with defaults.'
                 ! Fill this row with defaults
                 DO j = 1, k
                    col_type = self%type_mapping_array(1, j)
                    target_idx = self%type_mapping_array(2, j)
                    SELECT CASE(col_type)
                    CASE (TYPE_INT);  IF (ALLOCATED(self%int_data))  self%int_data(current_row, target_idx) = 0
                    CASE (TYPE_REAL); IF (ALLOCATED(self%real_data)) self%real_data(current_row, target_idx) = 0.0
                    CASE (TYPE_CHAR); IF (ALLOCATED(self%char_data)) self%char_data(current_row, target_idx) = ' '
                    END SELECT
                 END DO
                 IF (ALLOCATED(fields)) DEALLOCATE(fields)
                 CYCLE ! Move to the next row
            END IF

            ! parse fields based on predefined types
            DO j = 1, k
                ! allocate trimmed_field with the actual length of the current field
                ALLOCATE(CHARACTER(LEN=MAX(1,LEN_TRIM(fields(j)))) :: trimmed_field, STAT=alloc_stat) ! Len >= 1
                IF (alloc_stat /= 0) THEN
                     WRITE(*,'(A,I0,A,I0)') 'Error: Failed allocating trimmed_field at row', current_row, ', col', j
                     trimmed_field = ' ' ! assign blank -> avoids crash somehow even if there is no trimmed field dont remove
                ELSE
                     trimmed_field = TRIM(fields(j))
                     IF (LEN_TRIM(trimmed_field) == 0) trimmed_field = ' ' ! Ensure len>=1 for reads
                END IF

                col_type = self%type_mapping_array(1, j)
                target_idx = self%type_mapping_array(2, j)

                SELECT CASE(col_type)
                CASE (TYPE_INT)
                    IF (LEN_TRIM(trimmed_field) > 0) THEN ! avoid reading empty string as int
                        READ(trimmed_field, *, IOSTAT=parse_stat) self%int_data(current_row, target_idx)
                        IF (parse_stat /= 0) THEN
                            WRITE(*,'(A,I0,A,I0,A,A,A)') 'Warning: Failed parsing INT at row ', current_row, &
                                 ', col ', j, '. Field="', TRIM(fields(j)), '". Setting to 0.'
                            self%int_data(current_row, target_idx) = 0 ! default
                        END IF
                    ELSE
                         self%int_data(current_row, target_idx) = 0 ! default
                    END IF
                CASE (TYPE_REAL)
                     IF (LEN_TRIM(trimmed_field) > 0) THEN ! dont read empty string as real
                         READ(trimmed_field, *, IOSTAT=parse_stat) self%real_data(current_row, target_idx)
                         IF (parse_stat /= 0) THEN
                             WRITE(*,'(A,I0,A,I0,A,A,A)') 'Warning: Failed parsing REAL at row ', current_row, &
                                 ', col ', j, '. Field="', TRIM(fields(j)), '". Setting to 0.0.'
                             self%real_data(current_row, target_idx) = 0.0 ! default value or NaN
                         END IF
                     ELSE
                          self%real_data(current_row, target_idx) = 0.0 ! default for empty field
                     END IF
                CASE (TYPE_CHAR)
                     ! assign trimed original field -> Fortran handles padding/truncation(realy cool actually :D)
                     IF (ALLOCATED(self%char_data)) THEN
                         self%char_data(current_row, target_idx) = TRIM(fields(j))
                     END IF
                END SELECT

                ! deallocate trimmed_field for this iteration
                IF (ALLOCATED(trimmed_field)) DEALLOCATE(trimmed_field)
            END DO
            DEALLOCATE(fields) ! Deallocate fields array for this row
        END DO

        WRITE(*,*) 'Info: Finished reading data rows.'
        CLOSE(unit_num)
        RETURN

        ! --- error handling ---
 901    CONTINUE
        WRITE(*,*) 'Fatal Error (read_from_file): Memory allocation failed. ALLOC_STAT=', alloc_stat
        CALL self%destroy() ! Attempt cleanup
        STOP 901

    END SUBROUTINE datatable_read_from_file

    ! =========
    ! Function: datatable_calculate_memory
    ! =========
     FUNCTION datatable_calculate_memory(self) RESULT(total_bytes)
        CLASS(tox_datatable), INTENT(IN) :: self
        INTEGER(KIND=IK) :: total_bytes ! Use large integer kind

        ! --- local variable declarations ---
        INTEGER :: k, n
        INTEGER :: i_count, r_count, c_count
        INTEGER :: char_len
        INTEGER(KIND=IK) :: size_type_map, size_int_data, size_real_data
        INTEGER(KIND=IK) :: size_char_data, size_col_names
        INTEGER(KIND=IK) :: sizeof_int, sizeof_real, sizeof_char
        INTEGER :: calc_char_len ! Moved declaration here

        ! --- Executable Statements ---

        ! Get dimensions from the object (assuming they are set)
        n = self%num_rows
        k = self%num_cols
        i_count = self%num_int_cols
        r_count = self%num_real_cols
        c_count = self%num_char_cols
        char_len = self%max_char_len_used

        ! Basic type sizes (in bytes)
        sizeof_int = STORAGE_SIZE(0) / 8_IK       ! assumes 8 bits per storage unit 
        sizeof_real = STORAGE_SIZE(0.0) / 8_IK
        sizeof_char = STORAGE_SIZE('A') / 8_IK    ! typically 1 byte per char

        ! calculate size of each component 
        size_type_map = 2_IK * INT(k, KIND=IK) * sizeof_int

        size_int_data = INT(n, KIND=IK) * INT(i_count, KIND=IK) * sizeof_int
        size_real_data = INT(n, KIND=IK) * INT(r_count, KIND=IK) * sizeof_real

        ! Use 1 as min length for calculation if char columns exist but length is 0
        calc_char_len = char_len
        IF (c_count > 0 .AND. calc_char_len <= 0) THEN
             WRITE(*,*) 'Warning (calculate_memory): max_char_len_used is zero ', &
                        'or negative, using 1 for char data size calculation.'
             calc_char_len = 1
        END IF
        size_char_data = INT(n, KIND=IK) * INT(c_count, KIND=IK) * INT(calc_char_len, KIND=IK) * sizeof_char

        size_col_names = INT(k, KIND=IK) * INT(COL_NAME_LEN, KIND=IK) * sizeof_char

        ! Sum the components
        total_bytes = size_type_map + size_int_data + size_real_data + size_char_data + size_col_names

        ! Optional: Print breakdown
        WRITE(*,'(A, G0)') '* MemCalc: Type Map Bytes: ', size_type_map
        WRITE(*,'(A, G0)') '* MemCalc: Int Data Bytes: ', size_int_data
        WRITE(*,'(A, G0)') '* MemCalc: Real Data Bytes: ', size_real_data
        WRITE(*,'(A, G0)') '* MemCalc: Char Data Bytes: ', size_char_data
        WRITE(*,'(A, G0)') '* MemCalc: Col Names Bytes: ', size_col_names
        WRITE(*,'(A, G0)') '* MemCalc: Total Bytes:      ', total_bytes

    END FUNCTION datatable_calculate_memory

    ! =====================================================================
    ! Subroutine: datatable_destroy
    ! =====================================================================
    SUBROUTINE datatable_destroy(self)
        CLASS(tox_datatable), INTENT(INOUT) :: self

        ! --- Executable Statements Start Here ---
        IF (ALLOCATED(self%type_mapping_array)) DEALLOCATE(self%type_mapping_array)
        IF (ALLOCATED(self%int_data)) DEALLOCATE(self%int_data)
        IF (ALLOCATED(self%real_data)) DEALLOCATE(self%real_data)
        IF (ALLOCATED(self%char_data)) DEALLOCATE(self%char_data)
        IF (ALLOCATED(self%column_names)) DEALLOCATE(self%column_names)

        ! Reset metadata
        self%num_rows = 0
        self%num_cols = 0
        self%num_int_cols = 0
        self%num_real_cols = 0
        self%num_char_cols = 0
        self%max_char_len_used = 0

    END SUBROUTINE datatable_destroy

    ! =====================================================================
    ! Subroutine: datatable_print_summary
    ! =====================================================================
     SUBROUTINE datatable_print_summary(self)
        CLASS(tox_datatable), INTENT(IN) :: self
        ! --- Local Variable Declarations ---
        INTEGER :: j

        ! --- Executable Statements Start Here ---
        WRITE(*,'(A)')    '--- DataTable Summary ---'
        WRITE(*,'(A,I0)') 'Num Rows: ', self%num_rows
        WRITE(*,'(A,I0)') 'Num Cols: ', self%num_cols
        WRITE(*,'(A,I0)') 'Num Int Cols:  ', self%num_int_cols
        WRITE(*,'(A,I0)') 'Num Real Cols: ', self%num_real_cols
        WRITE(*,'(A,I0)') 'Num Char Cols: ', self%num_char_cols
        WRITE(*,'(A,I0)') 'Max Char Len Used: ', self%max_char_len_used

        IF (self%num_cols > 0 .AND. ALLOCATED(self%column_names) .AND. ALLOCATED(self%type_mapping_array)) THEN
             WRITE(*,'(A)') 'Column Names & Types (Type Map Index):'
             DO j = 1, self%num_cols
                 WRITE(*,'(2X, I3, A, A, A, I0, A, I0, A)') j, ': "', TRIM(self%column_names(j)), &
                                             '" (Type=', self%type_mapping_array(1,j), &
                                             ', Idx=', self%type_mapping_array(2,j), ')'
             END DO
        ELSE
             WRITE(*,*) "Column info not available (arrays not allocated)."
        END IF
        WRITE(*,'(A)')    '-------------------------'

    END SUBROUTINE datatable_print_summary

    ! =====================================================================
    ! Subroutine: datatable_print_data
    ! =====================================================================
    SUBROUTINE datatable_print_data(self, max_rows_to_print)
        CLASS(tox_datatable), INTENT(IN) :: self
        INTEGER, INTENT(IN), OPTIONAL    :: max_rows_to_print

        ! --- Local Variable Declarations ---
        INTEGER :: i, j, rows_to_print
        INTEGER :: target_idx, col_type
        INTEGER, PARAMETER :: default_width = 12 ! Moved declaration here
        INTEGER :: name_len, col_width          ! Moved declaration here
        CHARACTER(LEN=256) :: formatted_value    ! Buffer for formatted value length check

        ! --- executables ---
        IF (self%num_rows == 0 .OR. self%num_cols == 0 .OR. .NOT. ALLOCATED(self%type_mapping_array)) THEN
            WRITE(*,*) 'Info (print_data): Table is empty or not fully initialized.'
            RETURN
        END IF

        IF (PRESENT(max_rows_to_print)) THEN
            rows_to_print = MIN(self%num_rows, max_rows_to_print)
        ELSE
            rows_to_print = MIN(self%num_rows, 10) ! Default to printing max 10 rows
        END IF
        IF (rows_to_print <= 0) RETURN ! Nothing to print

        WRITE(*,'(A, I0, A)') '--- DataTable Data (First ', rows_to_print, ' Rows) ---'

        ! print header
        IF (ALLOCATED(self%column_names)) THEN
            WRITE(*,'(A)', ADVANCE='NO') 'Row | '
            DO j = 1, self%num_cols
                ! need to figure out a reasonable width, perhaps based on max(len(name), default_width)
                name_len = LEN_TRIM(self%column_names(j))
                col_width = MAX(name_len, default_width)
                WRITE(*,'(A,A)', ADVANCE='NO') RPAD(TRIM(ADJUSTL(self%column_names(j))), col_width), ' | '
            END DO
            WRITE(*,*)
        ELSE
            WRITE(*,*) '(Header names not allocated)'
        END IF

        ! Print Data Rows
        DO i = 1, rows_to_print
            WRITE(*,'(I3, A)', ADVANCE='NO') i, ' | '
            DO j = 1, self%num_cols
                col_type = self%type_mapping_array(1, j)
                target_idx = self%type_mapping_array(2, j)
                ! Calculate column width based on header (if available) or default
                IF (ALLOCATED(self%column_names)) THEN
                   name_len = LEN_TRIM(self%column_names(j))
                ELSE
                   name_len = 0
                END IF
                col_width = MAX(name_len, default_width)

                ! Get formatted value to determine padding needed
                formatted_value = get_formatted_value(self, i, j)

                SELECT CASE(col_type)
                CASE (TYPE_INT)
                    IF (ALLOCATED(self%int_data)) THEN
                        WRITE(*,'(I0)', ADVANCE='NO') self%int_data(i, target_idx) ! Use I0 for auto width
                    ELSE
                        WRITE(*,'(A)', ADVANCE='NO') 'N/A_I'
                    END IF
                CASE (TYPE_REAL)
                     IF (ALLOCATED(self%real_data)) THEN
                         WRITE(*,'(G0.5)', ADVANCE='NO') self%real_data(i, target_idx) ! G0.5 for auto width, 5 decimal places
                     ELSE
                         WRITE(*,'(A)', ADVANCE='NO') 'N/A_R'
                     END IF
                CASE (TYPE_CHAR)
                     IF (ALLOCATED(self%char_data)) THEN
                         WRITE(*,'(A)', ADVANCE='NO') TRIM(self%char_data(i, target_idx))
                     ELSE
                         WRITE(*,'(A)', ADVANCE='NO') 'N/A_C'
                     END IF
                END SELECT
                ! Pad after value (crude alignment based on formatted length)
                WRITE(*,'(A)', ADVANCE='NO') REPEAT(' ', MAX(0, col_width - LEN_TRIM(formatted_value) ) )
                WRITE(*,'(A)', ADVANCE='NO') ' | '
            END DO
            WRITE(*,*) ! Newline
        END DO
        WRITE(*,'(A)')    '----------------------------------'

    END SUBROUTINE datatable_print_data


    ! ===========================================
    ! Subroutine: datatable_serialize
    ! =====================
    SUBROUTINE datatable_serialize(self, filename, iostat_out)
        CLASS(tox_datatable), INTENT(IN) :: self
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(OUT) :: iostat_out ! 0 on success, non-zero on error

        ! --- local variables ---
        INTEGER :: unit_num = 12
        INTEGER :: io_stat
        TYPE(tox_datatable_header) :: header
        INTEGER(KIND=IK) :: current_pos, header_bytes

        ! --- Executable Statements Start Here ---
        iostat_out = 0 ! obviously assuming success

        ! --- validtae tbl state ---
        IF (self%num_rows <= 0 .OR. self%num_cols <= 0 .OR. &
            .NOT. ALLOCATED(self%type_mapping_array) .OR. &
            .NOT. ALLOCATED(self%column_names) ) THEN
            WRITE(*,*) 'Error (serialize): DataTable is not fully initialized or is empty.'
            iostat_out = -1 ! Custom error code for invalid state
            RETURN
        END IF
        ! check consistency of allocated arrays vs counts
        IF (self%num_int_cols > 0 .AND. .NOT. ALLOCATED(self%int_data)) GOTO 991
        IF (self%num_real_cols > 0 .AND. .NOT. ALLOCATED(self%real_data)) GOTO 991
        IF (self%num_char_cols > 0 .AND. .NOT. ALLOCATED(self%char_data)) GOTO 991

        ! --- header prep ---
        header%magic = FILE_MAGIC_NUMBER
        header%version = FILE_VERSION
        header%n_rows = self%num_rows
        header%k_cols = self%num_cols
        header%i_cols = self%num_int_cols
        header%r_cols = self%num_real_cols
        header%c_cols = self%num_char_cols
        header%max_c_len = self%max_char_len_used

        ! --- Open Output File ---
        OPEN(NEWUNIT=unit_num, FILE=filename, STATUS='REPLACE', ACTION='WRITE', &
             ACCESS='STREAM', FORM='UNFORMATTED', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'Error (serialize): Cannot open file "', TRIM(filename), '" for writing. IOSTAT=', io_stat
            iostat_out = io_stat
            RETURN
        END IF

        ! --- writes header (initially with offsets=0) and calculate its size ---
        current_pos = 1_IK ! start position
        WRITE(unit_num, POS=current_pos, IOSTAT=io_stat) header ! dummy header
        IF(io_stat /= 0) GOTO 993
        INQUIRE(UNIT=unit_num, POS=current_pos) ! position after header write
        header_bytes = current_pos - 1_IK
        header%header_size = INT(header_bytes, KIND=KIND(header%header_size)) ! store calculated size

        ! --- write data blocks and calculating offsets - THIS IS THE IMPORTANT PART ---
        header%offset_type_map = current_pos
        ! record the starting position of the upcoming data block in the corresponding header field

        WRITE(unit_num, POS=current_pos, IOSTAT=io_stat) self%type_mapping_array! 
        ! writes the entire array to the file starting at the current position(see line above)

        IF(io_stat /= 0) GOTO 993 ! error handling

        INQUIRE(UNIT=unit_num, POS=current_pos) ! updates current position to the point after the just written block 


        !now this is repeated for each datatype and the corresponding arrays are filled with their respective values
        IF (ALLOCATED(self%int_data)) THEN
            header%offset_int_data = current_pos
            WRITE(unit_num, POS=current_pos, IOSTAT=io_stat) self%int_data
             IF(io_stat /= 0) GOTO 993
            INQUIRE(UNIT=unit_num, POS=current_pos)
        ELSE
            header%offset_int_data = 0_IK ! if there is nothing of this value
        END IF

        IF (ALLOCATED(self%real_data)) THEN
            header%offset_real_data = current_pos
            WRITE(unit_num, POS=current_pos, IOSTAT=io_stat) self%real_data
             IF(io_stat /= 0) GOTO 993
            INQUIRE(UNIT=unit_num, POS=current_pos)
        ELSE
             header%offset_real_data = 0_IK
        END IF

        IF (ALLOCATED(self%char_data)) THEN
             header%offset_char_data = current_pos
             WRITE(unit_num, POS=current_pos, IOSTAT=io_stat) self%char_data
              IF(io_stat /= 0) GOTO 993
             INQUIRE(UNIT=unit_num, POS=current_pos)
        ELSE
             header%offset_char_data = 0_IK
        END IF

        header%offset_col_names = current_pos
        WRITE(unit_num, POS=current_pos, IOSTAT=io_stat) self%column_names
        IF(io_stat /= 0) GOTO 993
        INQUIRE(UNIT=unit_num, POS=current_pos)

        header%total_file_size = current_pos - 1_IK ! Total size is final position - 1

        ! --- Rewrite Header with Correct Offsets and Size ---
        WRITE(unit_num, POS=1_IK, IOSTAT=io_stat) header
        IF(io_stat /= 0) GOTO 993

        ! --- Finish ---
        CLOSE(unit_num, IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'Error (serialize): Failed closing file "', TRIM(filename), '". IOSTAT=', io_stat
            iostat_out = io_stat
            RETURN
        END IF

        WRITE(*,'(A,A)') 'Info (serialize): DataTable successfully serialized to "', TRIM(filename), '"'
        RETURN

        ! --- Error Handling ---
 991    CONTINUE
        WRITE(*,*) 'Error (serialize): Inconsistent DataTable state (missing allocated array for non-zero count).'
        iostat_out = -2
        RETURN

 993    CONTINUE
        WRITE(*,'(A,A,A,I0)') 'Error (serialize): Failed writing to file "', TRIM(filename), '". IOSTAT=', io_stat
        CLOSE(unit_num, STATUS='DELETE') ! clean up its mistakes
        iostat_out = io_stat
        RETURN

    END SUBROUTINE datatable_serialize

    ! =====================================================================
    ! Subroutine: datatable_deserialize
    ! =====================================================================
    SUBROUTINE datatable_deserialize(self, filename, iostat_out)
        CLASS(tox_datatable), INTENT(INOUT) :: self
        CHARACTER(LEN=*), INTENT(IN) :: filename
        INTEGER, INTENT(OUT) :: iostat_out ! 0 on success, non-zero on error

        ! --- Local Variable Declarations ---
        INTEGER :: unit_num = 13
        INTEGER :: io_stat, alloc_stat
        TYPE(tox_datatable_header) :: header
        INTEGER :: n, k, i_count, r_count, c_count, max_c_len
        INTEGER :: alloc_char_len ! Moved declaration here

        ! --- Executablas ---
        iostat_out = 0 ! Assume success(what else could it be)

        ! --- Destroy existing data ---
        CALL self%destroy()

        ! --- Open Input File ---
        OPEN(NEWUNIT=unit_num, FILE=filename, STATUS='OLD', ACTION='READ', &
             ACCESS='STREAM', FORM='UNFORMATTED', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'Error (deserialize): Cannot open file "', TRIM(filename), '" for reading. IOSTAT=', io_stat
            iostat_out = io_stat
            RETURN
        END IF

        ! --- Read and Verify Header ---
        READ(unit_num, POS=1_IK, IOSTAT=io_stat) header
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'Error (deserialize): Failed reading header from "', TRIM(filename), '". IOSTAT=', io_stat
            CLOSE(unit_num); iostat_out = io_stat; RETURN
        END IF

        IF (TRIM(header%magic) /= TRIM(FILE_MAGIC_NUMBER)) THEN
            WRITE(*,'(A,A,A,A)') 'Error (deserialize): Invalid magic number in file "', TRIM(filename), &
                                 '". Expected: "', TRIM(FILE_MAGIC_NUMBER), '", Got: "', TRIM(header%magic), '"'
            CLOSE(unit_num); iostat_out = -10; RETURN ! Custom error code
        END IF

        IF (header%version > FILE_VERSION) THEN
            WRITE(*,'(A,I0,A,I0)') 'Warning (deserialize): File version (', header%version, &
                                  ') is newer than supported version (', FILE_VERSION, '). Attempting to read anyway.'
            ! Potentially add more robust version checking/handling here
        ELSE IF (header%version < FILE_VERSION) THEN
             WRITE(*,'(A,I0,A,I0)') 'Warning (deserialize): File version (', header%version, &
                                  ') is older than supported version (', FILE_VERSION, '). Attempting to read anyway.'
        END IF

        ! --- Extract Dimensions from Header ---
        n = header%n_rows
        k = header%k_cols
        i_count = header%i_cols
        r_count = header%r_cols
        c_count = header%c_cols
        max_c_len = header%max_c_len
        IF (n <= 0 .OR. k <= 0 .OR. i_count < 0 .OR. r_count < 0 .OR. c_count < 0 .OR. &
           (i_count + r_count + c_count /= k) .OR. (c_count > 0 .AND. max_c_len <= 0) ) THEN
            WRITE(*,*) 'Error (deserialize): Invalid dimensions or counts read from header.'
            ! Break long write statement
            WRITE(*,'(A,I0,A,I0)') 'n=', n, ' k=', k
            WRITE(*,'(A,I0,A,I0)') 'i=', i_count, ' r=', r_count
            WRITE(*,'(A,I0,A,I0)') 'c=', c_count, ' max_c_len=', max_c_len
            CLOSE(unit_num); iostat_out = -11; RETURN
        END IF

        ! --- Allocate Arrays ---
        self%num_rows = n
        self%num_cols = k
        self%num_int_cols = i_count
        self%num_real_cols = r_count
        self%num_char_cols = c_count
        self%max_char_len_used = max_c_len

        ALLOCATE(self%type_mapping_array(2, k), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 911
        IF (i_count > 0) THEN
            ALLOCATE(self%int_data(n, i_count), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 911
        END IF
        IF (r_count > 0) THEN
            ALLOCATE(self%real_data(n, r_count), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 911
        END IF
        IF (c_count > 0) THEN
             alloc_char_len = MAX(1, max_c_len) ! Assign value here, after max_c_len is known
             ALLOCATE(CHARACTER(LEN=alloc_char_len) :: self%char_data(n, c_count), STAT=alloc_stat)
             IF (alloc_stat /= 0) GOTO 911
        END IF
        ALLOCATE(self%column_names(k), STAT=alloc_stat); IF (alloc_stat /= 0) GOTO 911

        ! --- Read Data Blocks ---
        READ(unit_num, POS=header%offset_type_map, IOSTAT=io_stat) self%type_mapping_array
        IF(io_stat /= 0) GOTO 912

        IF (header%offset_int_data > 0_IK .AND. ALLOCATED(self%int_data)) THEN
            READ(unit_num, POS=header%offset_int_data, IOSTAT=io_stat) self%int_data
            IF(io_stat /= 0) GOTO 912
        ELSE IF (header%offset_int_data > 0_IK .AND. .NOT. ALLOCATED(self%int_data)) THEN
            WRITE(*,*) 'Error (deserialize): Header indicates int_data but array not allocated.'
            GOTO 913
        END IF

        IF (header%offset_real_data > 0_IK .AND. ALLOCATED(self%real_data)) THEN
             READ(unit_num, POS=header%offset_real_data, IOSTAT=io_stat) self%real_data
             IF(io_stat /= 0) GOTO 912
        ELSE IF (header%offset_real_data > 0_IK .AND. .NOT. ALLOCATED(self%real_data)) THEN
             WRITE(*,*) 'Error (deserialize): Header indicates real_data but array not allocated.'
             GOTO 913
        END IF

        IF (header%offset_char_data > 0_IK .AND. ALLOCATED(self%char_data)) THEN
              READ(unit_num, POS=header%offset_char_data, IOSTAT=io_stat) self%char_data
              IF(io_stat /= 0) GOTO 912
        ELSE IF (header%offset_char_data > 0_IK .AND. .NOT. ALLOCATED(self%char_data)) THEN
              WRITE(*,*) 'Error (deserialize): Header indicates char_data but array not allocated.'
              GOTO 913
        END IF

        READ(unit_num, POS=header%offset_col_names, IOSTAT=io_stat) self%column_names
        IF(io_stat /= 0) GOTO 912

        ! --- Finish ---
        CLOSE(unit_num, IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A,A,I0)') 'Error (deserialize): Failed closing file "', TRIM(filename), '". IOSTAT=', io_stat
            iostat_out = io_stat
            RETURN ! Data might be partially loaded but inconsistent
        END IF

        WRITE(*,'(A,A)') 'Info (deserialize): DataTable successfully deserialized from "', TRIM(filename), '"'
        RETURN

        ! --- Error Handling ---
 911    CONTINUE
        WRITE(*,*) 'Fatal Error (deserialize): Memory allocation failed. ALLOC_STAT=', alloc_stat
        CLOSE(unit_num); CALL self%destroy(); iostat_out = alloc_stat; RETURN
 912    CONTINUE
        WRITE(*,'(A,A,A,I0)') 'Error (deserialize): Failed reading data block from "', TRIM(filename), '". IOSTAT=', io_stat
        CLOSE(unit_num); CALL self%destroy(); iostat_out = io_stat; RETURN
 913    CONTINUE
         WRITE(*,*) 'Error (deserialize): Mismatch between header offsets and allocated arrays.'
         CLOSE(unit_num); CALL self%destroy(); iostat_out = -12; RETURN

    END SUBROUTINE datatable_deserialize

    ! =====================================================================
    ! Helper Subroutine: split_line_basic
    ! =====================================================================
    SUBROUTINE split_line_basic(input_line, delim, parts, status)
        CHARACTER(LEN=*), INTENT(IN)  :: input_line
        CHARACTER(LEN=1), INTENT(IN)  :: delim
        CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: parts(:)
        INTEGER, INTENT(OUT)          :: status

        ! --- Local Variable Declarations ---
        CHARACTER(LEN=LEN(input_line)) :: temp_line
        INTEGER :: i, start_pos
        INTEGER, PARAMETER :: MAX_FIELDS_ESTIMATE = 2048 ! Adjust if needed
        CHARACTER(LEN=MAX_FIELDS_ESTIMATE) :: temp_parts(MAX_FIELDS_ESTIMATE) ! Fixed intermediate buffer
        INTEGER :: num_parts_found

        ! --- Executable Statements Start Here ---
        status = 0
        temp_line = TRIM(input_line) ! Work with trimmed line
        num_parts_found = 0
        start_pos = 1

        IF (LEN(temp_line) == 0) THEN ! Handle empty line
             IF (ALLOCATED(parts)) DEALLOCATE(parts)
             ALLOCATE(CHARACTER(LEN=1) :: parts(0), STAT=status) ! Allocate zero-size, length 1 is arbitrary
             RETURN
        END IF

        DO i = 1, LEN(temp_line)
            IF (temp_line(i:i) == delim) THEN
                num_parts_found = num_parts_found + 1
                IF (num_parts_found > MAX_FIELDS_ESTIMATE) THEN
                    status = -1; WRITE(*,*) 'Error(split): Exceeded MAX_FIELDS_ESTIMATE'; RETURN
                END IF
                temp_parts(num_parts_found) = temp_line(start_pos : i-1)
                start_pos = i + 1
            END IF
        END DO

        ! Add the last part
        num_parts_found = num_parts_found + 1
        IF (num_parts_found > MAX_FIELDS_ESTIMATE) THEN
            status = -1; WRITE(*,*) 'Error(split): Exceeded MAX_FIELDS_ESTIMATE'; RETURN
        END IF
        temp_parts(num_parts_found) = temp_line(start_pos : LEN(temp_line))

        ! Allocate 'parts' to the exact size found and copy
        IF (ALLOCATED(parts)) DEALLOCATE(parts)
        ! Allocate with sufficient length (max possible field length is input_line length)
        ALLOCATE(CHARACTER(LEN=MAX(1,LEN(input_line))) :: parts(num_parts_found), STAT=status) ! Ensure len >= 1
        IF (status /= 0) THEN
            WRITE(*,*) 'Error(split): Allocation failed for parts array.'; RETURN
        END IF

        DO i = 1, num_parts_found
             parts(i) = TRIM(temp_parts(i)) ! Copy trimmed part
        END DO

    END SUBROUTINE split_line_basic

    ! =====================================================================
    ! Helper Function: RPAD - Right pad a string (Internal)
    ! =====================================================================
    PURE FUNCTION RPAD(str_in, total_len) RESULT(str_out)
        CHARACTER(LEN=*), INTENT(IN) :: str_in
        INTEGER, INTENT(IN) :: total_len
        CHARACTER(LEN=total_len) :: str_out
        ! --- local variable declarations ---
        INTEGER :: len_in
        ! --- executable statements start here ---
        len_in = LEN_TRIM(str_in)
        IF (total_len <= 0) THEN
            str_out = ''
        ELSE IF (len_in >= total_len) THEN
            str_out = str_in(1:total_len)
        ELSE
            str_out = TRIM(str_in) // REPEAT(' ', total_len - len_in)
        END IF
    END FUNCTION RPAD

    ! =====================================================================
    ! Helper Function: get formatted value(internal)
    ! =====================================================================
     FUNCTION get_formatted_value(self, row, col_idx) RESULT(val_str)
        CLASS(tox_datatable), INTENT(IN) :: self
        INTEGER, INTENT(IN) :: row, col_idx
        CHARACTER(LEN=256) :: val_str ! buffer

        ! --- loval variables ---
        INTEGER :: target_idx, col_type

        ! --- executables---
        val_str = 'ERROR' ! Default
        IF (row > self%num_rows .OR. col_idx > self%num_cols .OR. &
            row < 1 .OR. col_idx < 1 .OR. &
           .NOT. ALLOCATED(self%type_mapping_array)) RETURN

        col_type = self%type_mapping_array(1, col_idx)
        target_idx = self%type_mapping_array(2, col_idx)

        SELECT CASE(col_type)
            CASE (TYPE_INT)
                 IF (ALLOCATED(self%int_data)) THEN
                     WRITE(val_str,'(I0)') self%int_data(row, target_idx)
                 ELSE
                     val_str = 'N/A_I'
                 END IF
            CASE (TYPE_REAL)
                  IF (ALLOCATED(self%real_data)) THEN
                      WRITE(val_str,'(G0.5)') self%real_data(row, target_idx)
                  ELSE
                      val_str = 'N/A_R'
                  END IF
            CASE (TYPE_CHAR)
                  IF (ALLOCATED(self%char_data)) THEN
                       val_str = TRIM(self%char_data(row, target_idx))
                  ELSE
                       val_str = 'N/A_C'
                  END IF
        END SELECT
        val_str = TRIM(val_str)

     END FUNCTION get_formatted_value


END MODULE tox_datatable_mod


! =====================================================================
! --- tox mod demo
! =====================================================================
PROGRAM main_tox_demo
    USE tox_datatable_mod 
    IMPLICIT NONE

    ! --- config ---
    CHARACTER(LEN=256) :: csv_filename = 'data_tox.csv'
    CHARACTER(LEN=256) :: bin_filename = 'data_tox.txdata' ! Binary file name
    CHARACTER(LEN=1)   :: separator = ','
    INTEGER, PARAMETER :: num_data_rows = 6 ! Expected number of data rows
    INTEGER, PARAMETER :: max_str_len = 50  ! Max length for string columns

    ! Define column types (k=4 columns: Int, Real, Char, Int)
    INTEGER, PARAMETER :: column_def_types(4) = [TYPE_INT, TYPE_REAL, TYPE_CHAR, TYPE_INT]
    CHARACTER(LEN=COL_NAME_LEN) :: column_def_names(4)

    ! --- Variables ---
    TYPE(tox_datatable) :: my_table
    TYPE(tox_datatable) :: loaded_table ! For deserialization test
    INTEGER(KIND=IK) :: mem_req
    INTEGER :: k_def, i_def, r_def, c_def
    INTEGER :: serialize_status, deserialize_status ! For checking IO results

    ! === Execution ===

    WRITE(*,*) '--- TOX DataTable Demo ---'

    ! Assign column names individually
    column_def_names(1) = 'SampleID'
    column_def_names(2) = 'Measurement'
    column_def_names(3) = 'Category'
    column_def_names(4) = 'QualityScore'

    ! --- Create a sample CSV file (data_tox.csv) ---
    CALL create_sample_csv(csv_filename, separator)

    ! --- Calculate expected memory requirements ---
    k_def = SIZE(column_def_types)
    i_def = COUNT(column_def_types == TYPE_INT)
    r_def = COUNT(column_def_types == TYPE_REAL)
    c_def = COUNT(column_def_types == TYPE_CHAR)
    mem_req = calculate_memory_static(num_data_rows, k_def, i_def, r_def, c_def, max_str_len)
    WRITE(*,'(A, G0, A)') 'Static Memory Estimate: ', mem_req, ' bytes'

    ! --- Initialize and Read Data from CSV ---
    CALL my_table%read_from_file( filename=csv_filename, &
                                  sep_char=separator, &
                                  column_types=column_def_types, &
                                  max_char_len=max_str_len, &
                                  n_rows=num_data_rows, &
                                  opt_column_names=column_def_names ) ! currently need to provide names

    ! --- Print Summary of Original Table ---
    WRITE(*,*) ; WRITE(*,*) '--- Original Table (from CSV) ---'
    CALL my_table%print_summary()

    ! --- Calculate Memory (using object state) ---
    mem_req = my_table%calculate_memory()
    WRITE(*,'(A, G0, A)') 'Memory from Object State: ', mem_req, ' bytes'

    ! --- Print Data of Original Table ---
    CALL my_table%print_data()

    ! --- Serialize the Table ---
    WRITE(*,*) ; WRITE(*,*) '--- Serializing Table ---'
    CALL my_table%serialize(bin_filename, serialize_status)

    ! --- Destroy the original table to ensure clean load ---
    CALL my_table%destroy()
    WRITE(*,*) 'Info: Original table destroyed.'

    ! --- Deserialize the Table into a new object ---
    WRITE(*,*) ; WRITE(*,*) '--- Deserializing Table ---'
    CALL loaded_table%deserialize(bin_filename, deserialize_status)

    ! --- Print Summary of Loaded Table ---
    WRITE(*,*) ; WRITE(*,*) '--- Loaded Table (from Binary) ---'
    CALL loaded_table%print_summary()

    ! --- Calculate Memory (using loaded object state) ---
    mem_req = loaded_table%calculate_memory()
    WRITE(*,'(A, G0, A)') 'Memory from Loaded Object State: ', mem_req, ' bytes'

    ! --- Print Data of Loaded Table ---
    CALL loaded_table%print_data()

    ! --- Cleanup ---
    CALL loaded_table%destroy() ! Destroy the loaded table too
    WRITE(*,*)
    WRITE(*,*) '--- Demo Finished ---'


CONTAINS

    ! =====================================================================
    ! Helper: Create a sample CSV for testing
    ! =====================================================================
    SUBROUTINE create_sample_csv(filename, sep)
        CHARACTER(LEN=*), INTENT(IN) :: filename
        CHARACTER(LEN=1), INTENT(IN) :: sep
        ! --- Local Variable Declarations ---
        INTEGER :: unit_num = 10
        INTEGER :: io_stat
        ! --- Executable Statements Start Here ---
        OPEN(UNIT=unit_num, FILE=filename, STATUS='REPLACE', ACTION='WRITE', IOSTAT=io_stat)
        IF (io_stat /= 0) THEN
            WRITE(*,'(A,A)') 'Error creating sample file: ', TRIM(filename)
            RETURN
        END IF

        WRITE(unit_num,'(I0,A,G0,A,A,A,I0)') 102, sep, 67.8, sep, 'Beta', sep, 88
        WRITE(unit_num,'(I0,A,G0,A,A,A,I0)') 101, sep, 123.45, sep, 'Alpha', sep, 95
        WRITE(unit_num,'(I0,A,G0,A,A,A,I0)') 102, sep, 67.8, sep, 'Beta', sep, 88
        WRITE(unit_num,'(I0,A,G0,A,A,A,I0)') 103, sep, 90.12, sep, 'Gamma', sep, 70
        WRITE(unit_num,'(A,A,A,A,A,A,A)') '104',sep,'-5.5',sep,'Alpha',sep,'91' ! testing string parsing for numeric
        WRITE(unit_num,'(I0,A,G0,A,A,A,I0)') 105, sep, 100.0, sep, 'Delta Long Name Test', sep, 85 ! testing long string

        CLOSE(unit_num)
        WRITE(*,'(A,A)') 'Created sample file: ', TRIM(filename)

    END SUBROUTINE create_sample_csv

    ! =====================================================================
    ! calculate memory based on i/r/c counts
    ! =====================================================================
     FUNCTION calculate_memory_static(n_rows, k_cols, i_count, r_count, c_count, max_char_len) RESULT(total_bytes)
        INTEGER, INTENT(IN) :: n_rows, k_cols, i_count, r_count, c_count, max_char_len
        INTEGER(KIND=IK) :: total_bytes ! large integer kind

        ! --- Local Variable Declarations ---
        INTEGER(KIND=IK) :: size_type_map, size_int_data, size_real_data
        INTEGER(KIND=IK) :: size_char_data, size_col_names
        INTEGER(KIND=IK) :: sizeof_int, sizeof_real, sizeof_char
        INTEGER :: local_max_char_len

        ! --- execs ---
        sizeof_int = STORAGE_SIZE(0) / 8_IK
        sizeof_real = STORAGE_SIZE(0.0) / 8_IK
        sizeof_char = STORAGE_SIZE('A') / 8_IK

        size_type_map = 2_IK * INT(k_cols, KIND=IK) * sizeof_int
        size_int_data = INT(n_rows, KIND=IK) * INT(i_count, KIND=IK) * sizeof_int
        size_real_data = INT(n_rows, KIND=IK) * INT(r_count, KIND=IK) * sizeof_real

        local_max_char_len = max_char_len
        IF (c_count > 0 .AND. local_max_char_len <= 0) THEN
             local_max_char_len = 1 ! 1 as minimal input
        END IF
        size_char_data = INT(n_rows, KIND=IK) * INT(c_count, KIND=IK) * INT(local_max_char_len, KIND=IK) * sizeof_char
        size_col_names = INT(k_cols, KIND=IK) * INT(COL_NAME_LEN, KIND=IK) * sizeof_char

        total_bytes = size_type_map + size_int_data + size_real_data + size_char_data + size_col_names
    END FUNCTION calculate_memory_static

END PROGRAM main_tox_demo
