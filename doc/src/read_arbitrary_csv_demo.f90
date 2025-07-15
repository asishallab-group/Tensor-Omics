PROGRAM data_serialization_demo
    USE csv_reader_module
    IMPLICIT NONE
    INTEGER :: io_stat, i, j
    TYPE(generic_data_cell) :: cell

    ! --- NEW: Variables for command-line argument handling ---
    CHARACTER(LEN=1024) :: input_filename, bin_filename
    INTEGER :: arg_count

    ! 1. Get filename from command-line arguments
    arg_count = COMMAND_ARGUMENT_COUNT()
    IF (arg_count /= 1) THEN
        WRITE(*,'(A)') 'Error: Incorrect number of arguments.'
        WRITE(*,'(A)') 'Usage: ./your_program_name <input_csv_file>'
        STOP 1 ! Stop with an error code
    END IF
    CALL GET_COMMAND_ARGUMENT(1, input_filename)
    bin_filename = TRIM(input_filename) // '.bin'

    ! 2. Read the specified CSV file into memory
    WRITE(*,'(A,A,A)') "--- 1. Reading data from '", TRIM(input_filename), "' ---"
    CALL read_csv_file(TRIM(input_filename), .TRUE., ',', io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,'(A,A,A,I0)') "ERROR: Failed to read '", TRIM(input_filename), "'. Status:", io_stat
        STOP
    END IF
    WRITE(*,*) "Read complete."
    WRITE(*,*)

    ! 3. Serialize the in-memory data to a binary file
    WRITE(*,'(A,A,A)') "--- 2. Serializing data to '", TRIM(bin_filename), "' ---"
    CALL serialize(TRIM(bin_filename), io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,'(A,I0)') "ERROR: Failed to serialize data. Status:", io_stat
        STOP
    END IF
    WRITE(*,*) "Serialization complete."
    WRITE(*,*)

    ! 4. Clean up the current in-memory data
    WRITE(*,*) "--- 3. Cleaning up in-memory data ---"
    CALL cleanup_csv_data()
    WRITE(*,*) "Cleanup complete. Rows in memory:", get_num_rows()
    WRITE(*,*)

    ! 5. Deserialize the data from the binary file
    WRITE(*,'(A,A,A)') "--- 4. Deserializing data from '", TRIM(bin_filename), "' ---"
    CALL deserialize(TRIM(bin_filename), io_stat)
    IF (io_stat /= 0) THEN
        WRITE(*,'(A,I0)') "ERROR: Failed to deserialize data. Status:", io_stat
        STOP
    END IF
    WRITE(*,*) "Deserialization complete."
    WRITE(*,*)

    ! 6. Print the deserialized data to verify it matches the original
    WRITE(*,*) "--- 5. Verifying deserialized data ---"
    WRITE(*,'(A,I0,A,I0,A)') "Table has ", get_num_rows(), " rows and ", get_num_cols(), " columns."
    ! Print header
    IF (get_num_cols() > 0) THEN
        WRITE(*,'(A)', ADVANCE='NO') 'RowH: '
        DO j = 1, get_num_cols()
            WRITE(*,'(A,A,A)', ADVANCE='NO') '[', TRIM(ADJUSTL(get_header(j))), '] '
        ENDDO
        WRITE(*,*)
    END IF
    ! Print data rows
    DO i = 1, get_num_rows()
        WRITE(*,'(A,I0,A)', ADVANCE='NO') 'Row', i, ': '
        DO j = 1, get_num_cols()
            cell = get_cell(i, j)
            SELECT CASE (cell%data_type)
            CASE (1) ! Integer
                WRITE(*,'(A,I0,A)', ADVANCE='NO') '[INT:', cell%i_val, '] '
            CASE (2) ! Real
                WRITE(*,'(A,G12.5,A)', ADVANCE='NO') '[REAL:', cell%r_val, '] '
            CASE (3) ! Character
                WRITE(*,'(A,A,A)', ADVANCE='NO') '[STR:"', TRIM(cell%c_val), '"] '
                IF(ALLOCATED(cell%c_val)) DEALLOCATE(cell%c_val)
            CASE DEFAULT ! Empty/Error
                WRITE(*,'(A)', ADVANCE='NO') '[ERR] '
            END SELECT
        END DO
        WRITE(*,*) ! Newline after each row
    END DO

    ! Final cleanup
    CALL cleanup_csv_data()
    WRITE(*,*)
    WRITE(*,*) "--- Program finished successfully ---"

END PROGRAM data_serialization_demo