! filepath: test/asserts.f90
MODULE asserts
    IMPLICIT NONE

    SAVE ! Ensure state persists across calls if needed, though not strictly for these assertions

    ! Global variable to track test status
    LOGICAL :: test_failed_flag = .FALSE.

CONTAINS

    !> @brief Verifies that a condition is true.
    SUBROUTINE assert_true(condition, message)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (.NOT. condition) THEN
            test_failed_flag = .TRUE.
            WRITE(*,*) "FAIL: ", TRIM(message)
            CALL error_stop_test() ! Stop on first failure
        END IF
    END SUBROUTINE assert_true

    !> @brief Verifies that a condition is false.
    SUBROUTINE assert_false(condition, message)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (condition) THEN
            test_failed_flag = .TRUE.
            WRITE(*,*) "FAIL: ", TRIM(message)
            CALL error_stop_test()
        END IF
    END SUBROUTINE assert_false

    !> @brief Compares real numbers.
    SUBROUTINE assert_equal_real(actual, expected, tolerance, message)
        REAL, INTENT(IN) :: actual, expected, tolerance
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (ABS(actual - expected) > tolerance) THEN
            test_failed_flag = .TRUE.
            WRITE(*,*) "FAIL: ", TRIM(message), &
                      " Actual: ", actual, " Expected: ", expected, " Tolerance: ", tolerance
            CALL error_stop_test()
        END IF
    END SUBROUTINE assert_equal_real

    !> @brief Compares integers.
    SUBROUTINE assert_equal_int(actual, expected, message)
        INTEGER, INTENT(IN) :: actual, expected
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (actual /= expected) THEN
            test_failed_flag = .TRUE.
            WRITE(*,*) "FAIL: ", TRIM(message), &
                      " Actual: ", actual, " Expected: ", expected
            CALL error_stop_test()
        END IF
    END SUBROUTINE assert_equal_int

    !> @brief Verifies a value is within a range.
    SUBROUTINE assert_in_range_real(value, min_val, max_val, message)
        REAL, INTENT(IN) :: value, min_val, max_val
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (value < min_val .OR. value > max_val) THEN
            test_failed_flag = .TRUE.
            WRITE(*,*) "FAIL: ", TRIM(message), &
                      " Value: ", value, " Min: ", min_val, " Max: ", max_val
            CALL error_stop_test()
        END IF
    END SUBROUTINE assert_in_range_real

    !> @brief Compares real arrays.
    SUBROUTINE assert_equal_array_real(actual, expected, num_elements, tolerance, message)
        REAL, DIMENSION(:), INTENT(IN) :: actual, expected
        INTEGER, INTENT(IN) :: num_elements
        REAL, INTENT(IN) :: tolerance
        CHARACTER(LEN=*), INTENT(IN) :: message
        INTEGER :: i
        IF (SIZE(actual) /= num_elements .OR. SIZE(expected) /= num_elements) THEN
            test_failed_flag = .TRUE.
            WRITE(*,*) "FAIL: ", TRIM(message), " Array sizes mismatch."
            CALL error_stop_test()
        END IF
        DO i = 1, num_elements
            IF (ABS(actual(i) - expected(i)) > tolerance) THEN
                test_failed_flag = .TRUE.
                WRITE(*,*) "FAIL: ", TRIM(message), " Mismatch at index ", i, &
                          " Actual: ", actual(i), " Expected: ", expected(i), " Tolerance: ", tolerance
                CALL error_stop_test()
            END IF
        END DO
    END SUBROUTINE assert_equal_array_real

    !> @brief Verifies absence of NaN.
    SUBROUTINE assert_no_nan_real(array, num_elements, message)
        REAL, DIMENSION(:), INTENT(IN) :: array
        INTEGER, INTENT(IN) :: num_elements
        CHARACTER(LEN=*), INTENT(IN) :: message
        INTEGER :: i
        DO i = 1, num_elements
            IF (IS_NAN(array(i))) THEN
                test_failed_flag = .TRUE.
                WRITE(*,*) "FAIL: ", TRIM(message), " NaN found at index ", i
                CALL error_stop_test()
            END IF
        END DO
    END SUBROUTINE assert_no_nan_real

    !> @brief Utility to stop execution on test failure
    SUBROUTINE error_stop_test()
        ! In a full framework, this might set a global flag and return,
        ! but for simple stopping, error stop is effective.
        ERROR STOP "Test failure detected."
    END SUBROUTINE error_stop_test

    !> @brief Function to reset the test_failed_flag
    SUBROUTINE reset_test_failed_flag()
        test_failed_flag = .FALSE.
    END SUBROUTINE reset_test_failed_flag

END MODULE asserts