!> @brief A library of assertion subroutines for testing.
MODULE asserts
    USE, INTRINSIC :: iso_fortran_env, ONLY: INT32, REAL64
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: assert_true, assert_false, assert_equal_int, assert_string_equal, &
              assert_equal_size, assert_equal_real, assert_equal_complex

CONTAINS

    !> @brief Fails the test if the condition is not true.
    SUBROUTINE assert_true(condition, message)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (.NOT. condition) THEN
            WRITE(*, '(A)') 'ASSERTION FAILED: ' // TRIM(message)
            ERROR STOP 1
        END IF
    END SUBROUTINE assert_true

    !> @brief Fails the test if the condition is not false.
    SUBROUTINE assert_false(condition, message)
        LOGICAL, INTENT(IN) :: condition
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (condition) THEN
            WRITE(*, '(A)') 'ASSERTION FAILED: ' // TRIM(message)
            ERROR STOP 1
        END IF
    END SUBROUTINE assert_false

    !> @brief Fails the test if two integers are not equal.
    SUBROUTINE assert_equal_int(actual, expected, message)
        INTEGER(INT32), INTENT(IN) :: actual, expected
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (actual /= expected) THEN
            WRITE(*, '(A)') 'ASSERTION FAILED: ' // TRIM(message)
            WRITE(*, '("  Got: ", I0, ", Expected: ", I0)') actual, expected
            ERROR STOP 1
        END IF
    END SUBROUTINE assert_equal_int

    !> @brief Fails the test if two real numbers are not equal within a tolerance.
    SUBROUTINE assert_equal_real(actual, expected, tolerance, message)
        REAL(REAL64), INTENT(IN) :: actual, expected, tolerance
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (ABS(actual - expected) > tolerance) THEN
            WRITE(*, '(A)') 'ASSERTION FAILED: ' // TRIM(message)
            WRITE(*, '("  Got: ", F0.10, ", Expected: ", F0.10)') actual, expected
            ERROR STOP 1
        END IF
    END SUBROUTINE assert_equal_real
    
    !> @brief Fails the test if two complex numbers are not equal within a tolerance.
    SUBROUTINE assert_equal_complex(actual, expected, tolerance, message)
        COMPLEX(REAL64), INTENT(IN) :: actual, expected
        REAL(REAL64), INTENT(IN) :: tolerance
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (ABS(REAL(actual) - REAL(expected)) > tolerance .OR. &
            ABS(AIMAG(actual) - AIMAG(expected)) > tolerance) THEN
            WRITE(*, '(A)') 'ASSERTION FAILED: ' // TRIM(message)
            WRITE(*, '("  Got: (", F0.10, ",", F0.10, "), Expected: (", F0.10, ",", F0.10, ")")') &
                REAL(actual), AIMAG(actual), REAL(expected), AIMAG(expected)
            ERROR STOP 1
        END IF
    END SUBROUTINE assert_equal_complex

    !> @brief Fails the test if two strings are not equal.
    SUBROUTINE assert_string_equal(actual, expected, message)
        CHARACTER(LEN=*), INTENT(IN) :: actual, expected
        CHARACTER(LEN=*), INTENT(IN) :: message
        IF (TRIM(actual) /= TRIM(expected)) THEN
            WRITE(*, '(A)') 'ASSERTION FAILED: ' // TRIM(message)
            WRITE(*, '("  Got: ''", A, "'', Expected: ''", A, "''")') TRIM(actual), TRIM(expected)
            ERROR STOP 1
        END IF
    END SUBROUTINE assert_string_equal

    !> @brief Fails the test if an array size is not equal to the expected size.
    SUBROUTINE assert_equal_size(actual_size, expected_size, message)
        INTEGER(INT32), INTENT(IN) :: actual_size, expected_size
        CHARACTER(LEN=*), INTENT(IN) :: message
        CALL assert_equal_int(actual_size, expected_size, message)
    END SUBROUTINE assert_equal_size

END MODULE asserts
