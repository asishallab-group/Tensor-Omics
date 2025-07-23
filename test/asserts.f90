! filepath: test/asserts.f90
MODULE asserts
    USE, INTRINSIC :: iso_fortran_env, ONLY: error_unit, real64, int64

    IMPLICIT NONE
    PRIVATE

    ! --- NEW: Define a parameter for the integer kind to match IK from csv_reader_module ---
    INTEGER, PARAMETER :: IK_ASSERT = SELECTED_INT_KIND(18) ! Should match IK in csv_reader_module

    PUBLIC :: assert_true, assert_false, assert_equal_int, assert_not_equal_int
    PUBLIC :: assert_equal_real, assert_not_equal_real, assert_equal_array_int
    PUBLIC :: assert_equal_array_real, assert_no_nan_real, assert_no_inf_real
    PUBLIC :: assert_in_range_real, assert_contains_int, assert_sorted_int
    PUBLIC :: assert_sorted_real, assert_same_shape, assert_string_equal
    PUBLIC :: assert_string_contains, assert_allclose_array_real
    PUBLIC :: assert_sum_equal, assert_unique_int, assert_permutation
    PUBLIC :: test_failed_flag, reset_test_failed_flag, error_stop_test_wrapper

    ! Add global flag for test runner to track overall suite status
    SAVE
    LOGICAL :: test_failed_flag = .FALSE.

CONTAINS

    !> Assert that a logical condition is true.
    SUBROUTINE assert_true(cond, msg)
        LOGICAL, INTENT(IN) :: cond
        CHARACTER(*), INTENT(IN) :: msg
        IF (.NOT. cond) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "FAIL: ", TRIM(msg)
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_true

    !> Assert that a logical condition is false.
    SUBROUTINE assert_false(cond, msg)
        LOGICAL, INTENT(IN) :: cond
        CHARACTER(*), INTENT(IN) :: msg
        IF (cond) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED (expected false): ", TRIM(msg)
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_false

    !> Assert that two integers are equal.
    SUBROUTINE assert_equal_int(a, b, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: a, b ! FIX: Use IK_ASSERT
        CHARACTER(*), INTENT(IN) :: msg
        IF (a /= b) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (got ", a, ", expected ", b, ")"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_equal_int

    !> Assert that two integers are not equal.
    SUBROUTINE assert_not_equal_int(a, b, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: a, b ! FIX: Use IK_ASSERT
        CHARACTER(*), INTENT(IN) :: msg
        IF (a == b) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED (should not be equal): ", TRIM(msg)
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_not_equal_int

    !> Assert that two real numbers are equal within a tolerance.
    SUBROUTINE assert_equal_real(a, b, tol, msg)
        REAL(real64), INTENT(IN) :: a, b, tol
        CHARACTER(*), INTENT(IN) :: msg
        IF (ABS(a-b) > tol) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (got ", a, ", expected ", b, ", tol=", tol, ")"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_equal_real

    !> Assert that two real numbers are not equal within a tolerance.
    SUBROUTINE assert_not_equal_real(a, b, tol, msg)
        REAL(real64), INTENT(IN) :: a, b, tol
        CHARACTER(*), INTENT(IN) :: msg
        IF (ABS(a-b) <= tol) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED (should not be equal): ", TRIM(msg)
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_not_equal_real

    !> Assert that two integer arrays are equal.
    SUBROUTINE assert_equal_array_int(a, b, n, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: a(n), b(n) ! FIX: Use IK_ASSERT
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        IF (ANY(a /= b)) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (integer arrays differ)"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_equal_array_int

    !> Assert that two real arrays are equal within a tolerance.
    SUBROUTINE assert_equal_array_real(a, b, n, tol, msg)
        REAL(real64), INTENT(IN) :: a(n), b(n), tol
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        IF (ANY(ABS(a-b) > tol)) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (real arrays differ, tol=", tol, ")"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_equal_array_real

    !> Assert that a real array contains no NaN values.
    SUBROUTINE assert_no_nan_real(a, n, msg)
        REAL(real64), INTENT(IN) :: a(n)
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i
        DO i = 1, n
            IF (a(i) /= a(i)) THEN ! This is the portable NaN check
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: NaN detected - ", TRIM(msg)
                CALL error_stop_test_wrapper()
            END IF
        END DO
    END SUBROUTINE assert_no_nan_real

    !> Assert that a real array contains no Inf values.
    SUBROUTINE assert_no_inf_real(a, n, msg)
        REAL(real64), INTENT(IN) :: a(n)
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i
        DO i = 1, n
            IF (ABS(a(i)) > HUGE(1.0_real64)) THEN
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: Inf detected - ", TRIM(msg)
                CALL error_stop_test_wrapper()
            END IF
        END DO
    END SUBROUTINE assert_no_inf_real

    !> Assert that a real value is within a given range [minval, maxval].
    SUBROUTINE assert_in_range_real(a, minval, maxval, msg)
        REAL(real64), INTENT(IN) :: a, minval, maxval
        CHARACTER(*), INTENT(IN) :: msg
        IF (a < minval .OR. a > maxval) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (value ", a, " not in [", minval, ",", maxval, "])"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_in_range_real

    !> Assert that an integer array contains a given value.
    SUBROUTINE assert_contains_int(arr, n, val, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: arr(n), val ! FIX: Use IK_ASSERT
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        IF (.NOT. ANY(arr == val)) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (value ", val, " not found)"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_contains_int

    !> Assert that an integer array is sorted in non-decreasing order.
    SUBROUTINE assert_sorted_int(arr, n, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: arr(n) ! FIX: Use IK_ASSERT
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i
        DO i = 2, n
            IF (arr(i) < arr(i-1)) THEN
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (not sorted at position ", i, ")"
                CALL error_stop_test_wrapper()
            END IF
        END DO
    END SUBROUTINE assert_sorted_int

    !> Assert that a real array is sorted in non-decreasing order.
    SUBROUTINE assert_sorted_real(arr, n, msg)
        REAL(real64), INTENT(IN) :: arr(n)
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i
        DO i = 2, n
            IF (arr(i) < arr(i-1)) THEN
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (not sorted at position ", i, ")"
                CALL error_stop_test_wrapper()
            END IF
        END DO
    END SUBROUTINE assert_sorted_real

    !> Assert that two arrays have the same shape (1D only).
    SUBROUTINE assert_same_shape(n1, n2, msg)
        INTEGER, INTENT(IN) :: n1, n2
        CHARACTER(*), INTENT(IN) :: msg
        IF (n1 /= n2) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (shapes differ: ", n1, " vs ", n2, ")"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_same_shape

    !> Assert that two strings are equal.
    SUBROUTINE assert_string_equal(a, b, msg)
        CHARACTER(*), INTENT(IN) :: a, b, msg
        IF (TRIM(a) /= TRIM(b)) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (got '"//TRIM(a)//"', expected '"//TRIM(b)//"')"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_string_equal

    !> Assert that string a contains string b.
    SUBROUTINE assert_string_contains(a, b, msg)
        CHARACTER(*), INTENT(IN) :: a, b, msg
        IF (INDEX(a, b) == 0) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (substring '"//TRIM(b)//"' not found in '"//TRIM(a)//"')"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_string_contains

    !> Assert that two real arrays are close within relative and absolute tolerance.
    SUBROUTINE assert_allclose_array_real(a, b, n, rtol, atol, msg)
        REAL(real64), INTENT(IN) :: a(n), b(n), rtol, atol
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i
        DO i = 1, n
            IF (ABS(a(i) - b(i)) > atol + rtol * ABS(b(i))) THEN
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (arrays differ at ", i, ")"
                CALL error_stop_test_wrapper()
            END IF
        END DO
    END SUBROUTINE assert_allclose_array_real

    !> Assert that the sum of an array equals an expected value.
    SUBROUTINE assert_sum_equal(arr, n, expected, msg)
        REAL(real64), INTENT(IN) :: arr(n), expected
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        REAL(real64) :: s
        s = SUM(arr)
        IF (ABS(s - expected) > 1e-12_real64) THEN
            test_failed_flag = .TRUE.
            WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (sum=", s, ", expected=", expected, ")"
            CALL error_stop_test_wrapper()
        END IF
    END SUBROUTINE assert_sum_equal

    !> Assert that all elements in an integer array are unique.
    SUBROUTINE assert_unique_int(arr, n, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: arr(n) ! FIX: Use IK_ASSERT
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i, j
        DO i = 1, n-1
            DO j = i+1, n
                IF (arr(i) == arr(j)) THEN
                    test_failed_flag = .TRUE.
                    WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), &
                                        " (duplicate value ", arr(i), &
                                        " at positions ", i, " and ", j, ")"
                    CALL error_stop_test_wrapper()
                END IF
            END DO
        END DO
    END SUBROUTINE assert_unique_int

    !> Assert that an integer array is a permutation of 1..n.
    SUBROUTINE assert_permutation(arr, n, msg)
        INTEGER(IK_ASSERT), INTENT(IN) :: arr(n) ! FIX: Use IK_ASSERT
        INTEGER, INTENT(IN) :: n
        CHARACTER(*), INTENT(IN) :: msg
        INTEGER :: i
        LOGICAL :: found(n)
        found = .FALSE.
        DO i = 1, n
            IF (arr(i) < 1 .OR. arr(i) > n) THEN
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (value out of range: ", arr(i), ")"
                CALL error_stop_test_wrapper()
            END IF
            IF (found(arr(i))) THEN
                test_failed_flag = .TRUE.
                WRITE(error_unit,*) "ASSERTION FAILED: ", TRIM(msg), " (duplicate value: ", arr(i), ")"
                CALL error_stop_test_wrapper()
            END IF
            found(arr(i)) = .TRUE.
        END DO
    END SUBROUTINE assert_permutation

    !> @brief Wrapper to stop execution on test failure.
    !> This allows `ERROR STOP` to be controlled globally by the test runner.
    SUBROUTINE error_stop_test_wrapper()
        ERROR STOP "Test failure detected by assertion."
    END SUBROUTINE error_stop_test_wrapper

    !> @brief Resets the global test_failed_flag.
    SUBROUTINE reset_test_failed_flag()
        test_failed_flag = .FALSE.
    END SUBROUTINE reset_test_failed_flag

END MODULE asserts