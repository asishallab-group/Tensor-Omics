!=======================================================================
! F42-Compliant Test Program for csv_reader_module
!
! This program serves as the main test executable for fpm. It will:
! 1. Create a dummy CSV file for testing.
! 2. Call read_csv_file and verify the dimensions and data.
! 3. Call serialize to write a binary version.
! 4. Call cleanup_csv_data to clear memory.
! 5. Call deserialize to read the binary file back.
! 6. Verify the data again to ensure the cycle was successful.
! 7. Clean up generated files.
!
! It will halt with ERROR STOP if any verification fails.
!=======================================================================
PROGRAM test_runner
    USE csv_reader_module
    IMPLICIT NONE

    INTEGER :: io_status, unit_num = 15
    TYPE(generic_data_cell) :: cell
    CHARACTER(LEN=MAX_FIELD_LEN) :: header
    REAL, PARAMETER :: TOL = 1e-6

    PRINT *, "F42 TEST: Starting csv_reader_module verification..."

    ! --- Setup: Create a temporary CSV file for the test ---
    OPEN(UNIT=unit_num, FILE="test_data.csv", STATUS="REPLACE", ACTION="WRITE")
    WRITE(unit_num, '(A)') "ID,Name,Value,Status"
    WRITE(unit_num, '(A)') "101,First,1.23,Active"
    WRITE(unit_num, '(A)') "102,Second,4.56,Inactive"
    WRITE(unit_num, '(A)') "103,Third,7.89,Active"
    CLOSE(unit_num)

    ! --- Test 1: Read CSV File ---
    CALL read_csv_file("test_data.csv", has_header=.TRUE., delimiter=',', io_status=io_status)
    IF (io_status /= 0) ERROR STOP "TEST FAILED: read_csv_file returned a non-zero status."
    IF (get_num_rows() /= 3) ERROR STOP "TEST FAILED: Incorrect number of rows found."
    IF (get_num_cols() /= 4) ERROR STOP "TEST FAILED: Incorrect number of columns found."
    PRINT *, "--> PASSED: read_csv_file"

    ! --- Test 2: Verify Header and Cell Data ---
    header = get_header(3)
    IF (TRIM(header) /= "Value") ERROR STOP "TEST FAILED: get_header(3) did not return 'Value'."

    cell = get_cell(2, 3) ! Row 2, Column 3 should be 4.56
    IF (cell%data_type /= 2) ERROR STOP "TEST FAILED: Data type for cell (2,3) is not REAL."
    IF (ABS(cell%r_val - 4.56) > TOL) ERROR STOP "TEST FAILED: Incorrect real value for cell (2,3)."

    cell = get_cell(3, 4) ! Row 3, Column 4 should be "Active"
    IF (cell%data_type /= 3) ERROR STOP "TEST FAILED: Data type for cell (3,4) is not CHARACTER."
    IF (TRIM(cell%c_val) /= "Active") ERROR STOP "TEST FAILED: Incorrect string value for cell (3,4)."
    IF(ALLOCATED(cell%c_val)) DEALLOCATE(cell%c_val)
    PRINT *, "--> PASSED: Data verification"

    ! --- Test 3: Serialize Data ---
    CALL serialize("test_data.bin", io_status)
    IF (io_status /= 0) ERROR STOP "TEST FAILED: serialize returned a non-zero status."
    PRINT *, "--> PASSED: serialize"

    ! --- Test 4: Cleanup and Deserialize ---
    CALL cleanup_csv_data()
    IF (get_num_rows() /= 0) ERROR STOP "TEST FAILED: cleanup_csv_data did not reset row count."

    CALL deserialize("test_data.bin", io_status)
    IF (io_status /= 0) ERROR STOP "TEST FAILED: deserialize returned a non-zero status."
    PRINT *, "--> PASSED: cleanup and deserialize"

    ! --- Test 5: Verify Deserialized Data ---
    IF (get_num_rows() /= 3) ERROR STOP "TEST FAILED: Row count incorrect after deserialize."
    IF (get_num_cols() /= 4) ERROR STOP "TEST FAILED: Column count incorrect after deserialize."
    cell = get_cell(1, 1) ! Row 1, Column 1 should be 101
    IF (cell%i_val /= 101) ERROR STOP "TEST FAILED: Incorrect integer value for cell (1,1) after deserialize."
    PRINT *, "--> PASSED: Deserialized data verification"

    ! --- Final Cleanup ---
    CALL cleanup_csv_data()
    ! CORRECTED WAY TO DELETE FILES:
    ! Open the file and immediately close it with STATUS='DELETE'
    OPEN(UNIT=unit_num, FILE="test_data.csv", STATUS="OLD", IOSTAT=io_status)
    IF (io_status == 0) CLOSE(unit_num, STATUS="DELETE")

    OPEN(UNIT=unit_num, FILE="test_data.bin", STATUS="OLD", IOSTAT=io_status)
    IF (io_status == 0) CLOSE(unit_num, STATUS="DELETE")

    PRINT *, "==============================================="
    PRINT *, "SUCCESS: All tests for csv_reader_module passed."
    PRINT *, "==============================================="

END PROGRAM test_runner