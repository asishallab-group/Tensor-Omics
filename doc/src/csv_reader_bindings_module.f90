! csv_reader_bindings_module.f90
MODULE csv_reader_bindings_module
    USE iso_c_binding      ! Essential for C interoperability
    USE csv_reader_module  ! Import your main CSV reader module
                           ! This gives access to MAX_FIELD_LEN, IK, etc.

CONTAINS

    ! --- Wrapper for read_csv_file ---
    ! Original: SUBROUTINE read_csv_file(filename, has_header, delimiter, io_status)
    SUBROUTINE read_csv_file_C(filename_ptr, filename_len, has_header_c, &
                               delimiter_char_ptr, io_status_ptr) BIND(C, &
                               NAME='read_csv_file_C') ! Line continuation
        TYPE(C_PTR) :: filename_ptr          ! Pointer to C string
        INTEGER(C_INT) :: filename_len       ! Length of C string
        LOGICAL(C_BOOL) :: has_header_c       ! Boolean flag (named differently to avoid conflict)
        TYPE(C_PTR) :: delimiter_char_ptr    ! Pointer to optional delimiter char
        INTEGER(C_INT) :: io_status_ptr      ! Pointer for OUT status

        CHARACTER(KIND=C_CHAR), POINTER :: filename_c_ptr(:) ! Pointer to array of C_CHAR
        CHARACTER(LEN=MAX_LINE_LEN) :: filename_fortran_local ! Local copy for processing
        CHARACTER(LEN=1) :: current_delimiter
        INTEGER :: iostat_val
        CHARACTER(LEN=1, KIND=C_CHAR), POINTER :: p_delimiter_char
        LOGICAL :: has_header_fortran
        INTEGER :: i_char

        ! Associate C string pointer with Fortran POINTER array of C_CHARs
        CALL C_F_POINTER(filename_ptr, filename_c_ptr, [filename_len])
        
        ! Convert array of C_CHARs to Fortran CHARACTER string (fixed)
        filename_fortran_local = ' ' ! Initialize with spaces
        IF (filename_len > 0) THEN
            DO i_char = 1, MIN(filename_len, LEN(filename_fortran_local))
                filename_fortran_local(i_char:i_char) = filename_c_ptr(i_char)
            END DO
        END IF

        ! Convert C_BOOL to Fortran LOGICAL
        has_header_fortran = has_header_c

        ! Handle optional delimiter
        IF (C_ASSOCIATED(delimiter_char_ptr)) THEN
            CALL C_F_POINTER(delimiter_char_ptr, p_delimiter_char)
            current_delimiter = p_delimiter_char
            CALL read_csv_file(TRIM(filename_fortran_local), has_header_fortran, &
                               current_delimiter, iostat_val)
        ELSE
            ! Call with default delimiter (which is now handled by your Fortran module's optional argument)
            CALL read_csv_file(TRIM(filename_fortran_local), has_header_fortran, &
                               io_status=iostat_val)
        END IF

        io_status_ptr = iostat_val ! Assign the status back to the C pointer

    END SUBROUTINE read_csv_file_C

    ! --- Wrapper for cleanup_csv_data ---
    ! Original: SUBROUTINE cleanup_csv_data()
    SUBROUTINE cleanup_csv_data_C() BIND(C, NAME='cleanup_csv_data_C')
        CALL cleanup_csv_data()
    END SUBROUTINE cleanup_csv_data_C

    ! --- Wrapper for serialize ---
    ! Original: SUBROUTINE serialize(filename, io_status)
    SUBROUTINE serialize_C(filename_ptr, filename_len, io_status_ptr) BIND(C, &
                               NAME='serialize_C') ! Line continuation
        TYPE(C_PTR) :: filename_ptr
        INTEGER(C_INT) :: filename_len
        INTEGER(C_INT) :: io_status_ptr

        CHARACTER(KIND=C_CHAR), POINTER :: filename_c_ptr(:)
        CHARACTER(LEN=MAX_LINE_LEN) :: filename_fortran_local
        INTEGER :: iostat_val
        INTEGER :: i_char

        CALL C_F_POINTER(filename_ptr, filename_c_ptr, [filename_len])
        filename_fortran_local = ' '
        IF (filename_len > 0) THEN
            DO i_char = 1, MIN(filename_len, LEN(filename_fortran_local))
                filename_fortran_local(i_char:i_char) = filename_c_ptr(i_char)
            END DO
        END IF

        CALL serialize(TRIM(filename_fortran_local), iostat_val)
        io_status_ptr = iostat_val
    END SUBROUTINE serialize_C

    ! --- Wrapper for deserialize ---
    ! Original: SUBROUTINE deserialize(filename, io_status)
    SUBROUTINE deserialize_C(filename_ptr, filename_len, io_status_ptr) BIND(C, &
                               NAME='deserialize_C') ! Line continuation
        TYPE(C_PTR) :: filename_ptr
        INTEGER(C_INT) :: filename_len
        INTEGER(C_INT) :: io_status_ptr

        CHARACTER(KIND=C_CHAR), POINTER :: filename_c_ptr(:)
        CHARACTER(LEN=MAX_LINE_LEN) :: filename_fortran_local
        INTEGER :: iostat_val
        INTEGER :: i_char

        CALL C_F_POINTER(filename_ptr, filename_c_ptr, [filename_len])
        filename_fortran_local = ' '
        IF (filename_len > 0) THEN
            DO i_char = 1, MIN(filename_len, LEN(filename_fortran_local))
                filename_fortran_local(i_char:i_char) = filename_c_ptr(i_char)
            END DO
        END IF

        CALL deserialize(TRIM(filename_fortran_local), iostat_val)
        io_status_ptr = iostat_val
    END SUBROUTINE deserialize_C

    ! --- Wrapper for get_num_rows ---
    ! Original: FUNCTION get_num_rows() RESULT(n)
    FUNCTION get_num_rows_C() RESULT(rows) BIND(C, NAME='get_num_rows_C')
        INTEGER(C_INT) :: rows
        rows = get_num_rows()
    END FUNCTION get_num_rows_C

    ! --- Wrapper for get_num_cols ---
    ! Original: FUNCTION get_num_cols() RESULT(n)
    FUNCTION get_num_cols_C() RESULT(cols) BIND(C, NAME='get_num_cols_C')
        INTEGER(C_INT) :: cols
        cols = get_num_cols()
    END FUNCTION get_num_cols_C

    ! --- Wrapper for get_header ---
    ! Original: FUNCTION get_header(j) RESULT(header_name)
    SUBROUTINE get_header_C(j, header_name_ptr, header_name_len_ptr) BIND(C, &
                               NAME='get_header_C') ! Line continuation
        INTEGER(C_INT) :: j
        TYPE(C_PTR) :: header_name_ptr    ! Pointer to C string (Fortran-allocated)
        INTEGER(C_INT) :: header_name_len_ptr ! Pointer to return length

        ! Temporary allocatable to hold the result of get_header, must be TARGET
        CHARACTER(LEN=:), ALLOCATABLE, TARGET :: header_name_fortran_alloc

        ! Call your Fortran getter, which returns an allocated character string
        header_name_fortran_alloc = get_header(j)

        IF (ALLOCATED(header_name_fortran_alloc)) THEN
            ! Get a C pointer to the start of the Fortran character string's data
            header_name_ptr = C_LOC(header_name_fortran_alloc)
            ! Return the trimmed length of the string
            header_name_len_ptr = LEN_TRIM(header_name_fortran_alloc)
        ELSE
            header_name_ptr = C_NULL_PTR
            header_name_len_ptr = 0
        END IF
        ! IMPORTANT: 'header_name_fortran_alloc' is allocated here. Its memory
        ! must be explicitly DEALLOCATED later by a Python call to free_fortran_char_scalar_C.
    END SUBROUTINE get_header_C

    ! Deallocation wrapper for CHARACTER(LEN=:) scalars (like header names or single cell strings)
    SUBROUTINE free_fortran_char_scalar_C(data_ptr) BIND(C, &
                               NAME='free_fortran_char_scalar_C') ! Line continuation
        TYPE(C_PTR) :: data_ptr
        CHARACTER(LEN=:), POINTER :: p_char_scalar_dealloc

        IF (C_ASSOCIATED(data_ptr)) THEN
            ! C_F_POINTER can handle deferred-length pointers, assuming the length info is embedded
            CALL C_F_POINTER(data_ptr, p_char_scalar_dealloc)
            IF (ASSOCIATED(p_char_scalar_dealloc)) THEN
                DEALLOCATE(p_char_scalar_dealloc)
            END IF
        END IF
    END SUBROUTINE free_fortran_char_scalar_C

    ! --- Wrapper for get_cell's functionality: Per-type access ---
    ! These are similar to get_X_column, but for a single element.
    ! Python should call get_column_data_type_by_index_C first, then the appropriate getter.

    FUNCTION get_int_cell_by_index_C(i, j) RESULT(val) BIND(C, &
                               NAME='get_int_cell_by_index_C') ! Line continuation
        INTEGER(C_INT) :: i, j
        INTEGER(IK) :: val
        TYPE(generic_data_cell) :: cell_data_temp
        cell_data_temp = get_cell(i, j)
        IF (cell_data_temp%data_type == 1) THEN
            val = cell_data_temp%i_val
        ELSE
            val = 0_IK ! Default for incorrect type or error
        END IF
        ! Deallocate temporary character component if it was allocated by get_cell
        IF (ALLOCATED(cell_data_temp%c_val)) DEALLOCATE(cell_data_temp%c_val)
    END FUNCTION get_int_cell_by_index_C

    FUNCTION get_real_cell_by_index_C(i, j) RESULT(val) BIND(C, &
                               NAME='get_real_cell_by_index_C') ! Line continuation
        INTEGER(C_INT) :: i, j
        REAL :: val
        TYPE(generic_data_cell) :: cell_data_temp
        cell_data_temp = get_cell(i, j)
        IF (cell_data_temp%data_type == 2) THEN
            val = cell_data_temp%r_val
        ELSE
            val = 0.0 ! Default for incorrect type or error
        END IF
        IF (ALLOCATED(cell_data_temp%c_val)) DEALLOCATE(cell_data_temp%c_val)
    END FUNCTION get_real_cell_by_index_C

    FUNCTION get_logical_cell_by_index_C(i, j) RESULT(val) BIND(C, &
                               NAME='get_logical_cell_by_index_C') ! Line continuation
        INTEGER(C_INT) :: i, j
        LOGICAL(C_BOOL) :: val
        TYPE(generic_data_cell) :: cell_data_temp
        cell_data_temp = get_cell(i, j)
        IF (cell_data_temp%data_type == 4) THEN
            val = cell_data_temp%l_val
        ELSE
            val = .FALSE._C_BOOL ! Use _C_BOOL for logical as it's part of iso_c_binding
        END IF
        IF (ALLOCATED(cell_data_temp%c_val)) DEALLOCATE(cell_data_temp%c_val)
    END FUNCTION get_logical_cell_by_index_C

    FUNCTION get_complex_cell_by_index_C(i, j) RESULT(val) BIND(C, &
                               NAME='get_complex_cell_by_index_C') ! Line continuation
        INTEGER(C_INT) :: i, j
        COMPLEX :: val
        TYPE(generic_data_cell) :: cell_data_temp
        cell_data_temp = get_cell(i, j)
        IF (cell_data_temp%data_type == 5) THEN
            val = cell_data_temp%co_val
        ELSE
            val = (0.0, 0.0) ! Default for incorrect type or error
        END IF
        IF (ALLOCATED(cell_data_temp%c_val)) DEALLOCATE(cell_data_temp%c_val)
    END FUNCTION get_complex_cell_by_index_C

    ! For character cells, we will use a dedicated function that allocates and returns a pointer
    SUBROUTINE get_char_cell_by_index_C(i, j, char_ptr, char_len_ptr) BIND(C, &
                               NAME='get_char_cell_by_index_C') ! Line continuation
        INTEGER(C_INT) :: i, j
        TYPE(C_PTR) :: char_ptr
        INTEGER(C_INT) :: char_len_ptr

        CHARACTER(LEN=:), ALLOCATABLE, TARGET :: allocated_char_val ! Allocate new memory
        TYPE(generic_data_cell) :: cell_data_temp

        char_ptr = C_NULL_PTR
        char_len_ptr = 0

        cell_data_temp = get_cell(i, j)
        IF (cell_data_temp%data_type == 3 .AND. ALLOCATED(cell_data_temp%c_val)) THEN
            ALLOCATE(CHARACTER(LEN=LEN_TRIM(cell_data_temp%c_val)) :: allocated_char_val)
            allocated_char_val = cell_data_temp%c_val
            char_ptr = C_LOC(allocated_char_val)
            char_len_ptr = LEN_TRIM(allocated_char_val)
        END IF
        IF (ALLOCATED(cell_data_temp%c_val)) DEALLOCATE(cell_data_temp%c_val)
        ! IMPORTANT: 'allocated_char_val' is allocated here. Its memory
        ! must be explicitly DEALLOCATED later by a Python call to free_fortran_char_scalar_C.
    END SUBROUTINE get_char_cell_by_index_C

    ! --- Wrapper for get_column_data_type_by_name ---
    ! Original: FUNCTION get_column_data_type_by_name(col_name) RESULT(col_type)
    SUBROUTINE get_column_data_type_by_name_C(col_name_ptr, col_name_len, &
                                             col_type_ptr) BIND(C, &
                                             NAME='get_column_data_type_by_name_C') ! Line continuation
        TYPE(C_PTR) :: col_name_ptr
        INTEGER(C_INT) :: col_name_len
        INTEGER(C_INT) :: col_type_ptr ! Pointer for OUT argument

        CHARACTER(KIND=C_CHAR), POINTER :: col_name_c_ptr(:) ! Pointer to array of C_CHAR
        CHARACTER(LEN=MAX_FIELD_LEN) :: col_name_fortran_local ! Local copy for processing
        INTEGER :: f_col_type_val
        INTEGER :: i_char

        CALL C_F_POINTER(col_name_ptr, col_name_c_ptr, [col_name_len])
        col_name_fortran_local = ' '
        IF (col_name_len > 0) THEN
            DO i_char = 1, MIN(col_name_len, LEN(col_name_fortran_local))
                col_name_fortran_local(i_char:i_char) = col_name_c_ptr(i_char)
            END DO
        END IF

        f_col_type_val = get_column_data_type_by_name(TRIM(col_name_fortran_local))
        col_type_ptr = f_col_type_val
    END SUBROUTINE get_column_data_type_by_name_C

    ! --- Wrappers for column getters (returning Fortran-allocated memory) ---
    ! (These were already implemented in the previous step, included here for completeness)

    ! Wrapper for get_int_column
    ! Original: FUNCTION get_int_column(j) RESULT(column_data)
    SUBROUTINE get_int_column_C(j, data_ptr, num_elements_ptr) BIND(C, &
                               NAME='get_int_column_C') ! Line continuation
        INTEGER(C_INT) :: j                          ! Column index
        TYPE(C_PTR) :: data_ptr                      ! Pointer to allocated data
        INTEGER(C_INT) :: num_elements_ptr           ! Pointer to return array size

        INTEGER(IK), ALLOCATABLE, TARGET :: column_data_fortran(:)
        INTEGER :: actual_col_type

        actual_col_type = get_column_data_type_by_index(j)

        IF (actual_col_type == 1) THEN ! Check if it's an integer column (1)
            column_data_fortran = get_int_column(j) ! Call your Fortran getter
            IF (ALLOCATED(column_data_fortran)) THEN
                data_ptr = C_LOC(column_data_fortran(1)) ! Get C pointer to first element
                num_elements_ptr = SIZE(column_data_fortran) ! Get size
            ELSE
                data_ptr = C_NULL_PTR
                num_elements_ptr = 0
            END IF
        ELSE
            data_ptr = C_NULL_PTR
            num_elements_ptr = 0
        END IF
    END SUBROUTINE get_int_column_C

    ! Deallocation wrapper for INTEGER(IK) arrays returned by get_int_column_C
    SUBROUTINE free_fortran_int_array_C(data_ptr, num_elements) BIND(C, &
                               NAME='free_fortran_int_array_C') ! Line continuation
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements
        INTEGER(IK), POINTER :: p_array(:)

        IF (C_ASSOCIATED(data_ptr) .AND. num_elements > 0) THEN
            CALL C_F_POINTER(data_ptr, p_array, [num_elements])
            IF (ASSOCIATED(p_array)) THEN
                DEALLOCATE(p_array)
            END IF
        END IF
    END SUBROUTINE free_fortran_int_array_C

    ! Wrapper for get_real_column
    ! Original: FUNCTION get_real_column(j) RESULT(column_data)
    SUBROUTINE get_real_column_C(j, data_ptr, num_elements_ptr) BIND(C, &
                               NAME='get_real_column_C') ! Line continuation
        INTEGER(C_INT) :: j
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements_ptr

        REAL, ALLOCATABLE, TARGET :: column_data_fortran(:)
        INTEGER :: actual_col_type

        actual_col_type = get_column_data_type_by_index(j)
        IF (actual_col_type == 2) THEN ! 2 is Real type
            column_data_fortran = get_real_column(j)
            IF (ALLOCATED(column_data_fortran)) THEN
                data_ptr = C_LOC(column_data_fortran(1))
                num_elements_ptr = SIZE(column_data_fortran)
            ELSE
                data_ptr = C_NULL_PTR
                num_elements_ptr = 0
            END IF
        ELSE
            data_ptr = C_NULL_PTR
            num_elements_ptr = 0
        END IF
    END SUBROUTINE get_real_column_C

    ! Deallocation wrapper for REAL arrays
    SUBROUTINE free_fortran_real_array_C(data_ptr, num_elements) BIND(C, &
                               NAME='free_fortran_real_array_C') ! Line continuation
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements
        REAL, POINTER :: p_array(:)
        IF (C_ASSOCIATED(data_ptr) .AND. num_elements > 0) THEN
            CALL C_F_POINTER(data_ptr, p_array, [num_elements])
            IF (ASSOCIATED(p_array)) THEN
                DEALLOCATE(p_array)
            END IF
        END IF
    END SUBROUTINE free_fortran_real_array_C

    ! Wrapper for get_logical_column
    ! Original: FUNCTION get_logical_column(j) RESULT(column_data)
    SUBROUTINE get_logical_column_C(j, data_ptr, num_elements_ptr) BIND(C, &
                               NAME='get_logical_column_C') ! Line continuation
        INTEGER(C_INT) :: j
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements_ptr

        LOGICAL, ALLOCATABLE, TARGET :: column_data_fortran(:)
        INTEGER :: actual_col_type

        actual_col_type = get_column_data_type_by_index(j)
        IF (actual_col_type == 4) THEN ! 4 is Logical type
            column_data_fortran = get_logical_column(j)
            IF (ALLOCATED(column_data_fortran)) THEN
                data_ptr = C_LOC(column_data_fortran(1))
                num_elements_ptr = SIZE(column_data_fortran)
            ELSE
                data_ptr = C_NULL_PTR
                num_elements_ptr = 0
            END IF
        ELSE
            data_ptr = C_NULL_PTR
            num_elements_ptr = 0
        END IF
    END SUBROUTINE get_logical_column_C

    ! Deallocation wrapper for LOGICAL arrays
    SUBROUTINE free_fortran_logical_array_C(data_ptr, num_elements) BIND(C, &
                               NAME='free_fortran_logical_array_C') ! Line continuation
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements
        LOGICAL, POINTER :: p_array(:)
        IF (C_ASSOCIATED(data_ptr) .AND. num_elements > 0) THEN
            CALL C_F_POINTER(data_ptr, p_array, [num_elements])
            IF (ASSOCIATED(p_array)) THEN
                DEALLOCATE(p_array)
            END IF
        END IF
    END SUBROUTINE free_fortran_logical_array_C

    ! Wrapper for get_complex_column
    ! Original: FUNCTION get_complex_column(j) RESULT(column_data)
    SUBROUTINE get_complex_column_C(j, data_ptr, num_elements_ptr) BIND(C, &
                               NAME='get_complex_column_C') ! Line continuation
        INTEGER(C_INT) :: j
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements_ptr

        COMPLEX, ALLOCATABLE, TARGET :: column_data_fortran(:)
        INTEGER :: actual_col_type

        actual_col_type = get_column_data_type_by_index(j)
        IF (actual_col_type == 5) THEN ! 5 is Complex type
            column_data_fortran = get_complex_column(j)
            IF (ALLOCATED(column_data_fortran)) THEN
                data_ptr = C_LOC(column_data_fortran(1))
                num_elements_ptr = SIZE(column_data_fortran)
            ELSE
                data_ptr = C_NULL_PTR
                num_elements_ptr = 0
            END IF
        ELSE
            data_ptr = C_NULL_PTR
            num_elements_ptr = 0
        END IF
    END SUBROUTINE get_complex_column_C

    ! Deallocation wrapper for COMPLEX arrays
    SUBROUTINE free_fortran_complex_array_C(data_ptr, num_elements) BIND(C, &
                               NAME='free_fortran_complex_array_C') ! Line continuation
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements
        COMPLEX, POINTER :: p_array(:)
        IF (C_ASSOCIATED(data_ptr) .AND. num_elements > 0) THEN
            CALL C_F_POINTER(data_ptr, p_array, [num_elements])
            IF (ASSOCIATED(p_array)) THEN
                DEALLOCATE(p_array)
            END IF
        END IF
    END SUBROUTINE free_fortran_complex_array_C

    ! Wrapper for get_char_column
    ! Original: FUNCTION get_char_column(j) RESULT(column_data)
    SUBROUTINE get_char_column_C(j, data_ptr, num_elements_ptr, max_char_len_ptr) BIND(C, &
                               NAME='get_char_column_C') ! Line continuation
        INTEGER(C_INT) :: j
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements_ptr
        INTEGER(C_INT) :: max_char_len_ptr

        CHARACTER(LEN=MAX_FIELD_LEN), ALLOCATABLE, TARGET :: column_data_fortran(:)
        INTEGER :: actual_col_type

        actual_col_type = get_column_data_type_by_index(j)
        IF (actual_col_type == 3) THEN ! 3 is Character type
            column_data_fortran = get_char_column(j)
            IF (ALLOCATED(column_data_fortran)) THEN
                data_ptr = C_LOC(column_data_fortran(1))
                num_elements_ptr = SIZE(column_data_fortran)
                max_char_len_ptr = MAX_FIELD_LEN ! Pass the fixed length of each string
            ELSE
                data_ptr = C_NULL_PTR
                num_elements_ptr = 0
                max_char_len_ptr = 0
            END IF
        ELSE
            data_ptr = C_NULL_PTR
            num_elements_ptr = 0
            max_char_len_ptr = 0
        END IF
    END SUBROUTINE get_char_column_C

    ! Deallocation wrapper for CHARACTER arrays (fixed-length, contiguous block)
    SUBROUTINE free_fortran_char_array_C(data_ptr, num_elements, max_char_len) BIND(C, &
                               NAME='free_fortran_char_array_C') ! Line continuation
        TYPE(C_PTR) :: data_ptr
        INTEGER(C_INT) :: num_elements
        INTEGER(C_INT) :: max_char_len
        CHARACTER(LEN=MAX_FIELD_LEN), POINTER :: p_array_dealloc(:)
        
        IF (C_ASSOCIATED(data_ptr) .AND. num_elements > 0 .AND. max_char_len > 0) THEN
            ! Provide the shape explicitly.
            CALL C_F_POINTER(data_ptr, p_array_dealloc, [num_elements])
            IF (ASSOCIATED(p_array_dealloc)) THEN
                DEALLOCATE(p_array_dealloc)
            END IF
        END IF
    END SUBROUTINE free_fortran_char_array_C

END MODULE csv_reader_bindings_module