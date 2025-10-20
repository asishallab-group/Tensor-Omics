!> General assertion utilities for Fortran unit testing.
!! Provides a set of reusable assertion subroutines for verifying
!! expected behavior in tests of any kind (numeric, string, array, etc).
module asserts
  use, intrinsic :: iso_fortran_env, only: error_unit, real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
  implicit none
  private
  public :: assert_true, assert_false, assert_equal_int, assert_not_equal_int
  public :: assert_equal_real, assert_not_equal_real, assert_equal_array_int
  public :: assert_equal_array_real, assert_no_nan_real, assert_no_inf_real
  public :: assert_in_range_real, assert_contains_int, assert_sorted_int
  public :: assert_sorted_real, assert_same_shape, assert_string_equal
  public :: assert_string_contains, assert_allclose_array_real
  public :: assert_sum_equal, assert_unique_int, assert_permutation
  public :: assert_equal_array_char
  public :: assert_equal_complex, assert_not_equal_complex, assert_equal_array_complex


contains
  
  !> Asserts that two character arrays are equal
  subroutine assert_equal_array_char(a, b, clen, n, msg)
    integer(int32), INTENT(IN) :: clen
    character(len=clen), intent(in) :: a(n), b(n)
    character(*), intent(in) :: msg
    integer, intent(in) :: n
    if (any(a /= b)) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (character arrays differ)"
      stop 1
    end if
  end subroutine

  !> Assert that two complex numbers are equal within a tolerance.
  subroutine assert_equal_complex(a, b, tol, msg)
    complex(real64), intent(in) :: a, b
    real(real64), intent(in) :: tol
    character(*), intent(in) :: msg
    if (abs(a - b) > tol) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), &
           " (got ", a, ", expected ", b, ", tol=", tol, ")"
      stop 1
    end if
  end subroutine

  !> Assert that two complex numbers are not equal within a tolerance.
  subroutine assert_not_equal_complex(a, b, tol, msg)
    complex(real64), intent(in) :: a, b
    real(real64), intent(in) :: tol
    character(*), intent(in) :: msg
    if (abs(a - b) <= tol) then
      write(error_unit,*) "ASSERTION FAILED (should not be equal): ", trim(msg)
      stop 1
    end if
  end subroutine

  !> Assert that two complex arrays are equal within a tolerance.
  subroutine assert_equal_array_complex(a, b, n, tol, msg)
    complex(real64), intent(in) :: a(n), b(n)
    integer, intent(in) :: n
    real(real64), intent(in) :: tol
    character(*), intent(in) :: msg
    if (any(abs(a - b) > tol)) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), &
           " (complex arrays differ, tol=", tol, ")"
      stop 1
    end if
  end subroutine

  !> Assert that a logical condition is true.
  subroutine assert_true(cond, msg)
    logical, intent(in) :: cond
    character(*), intent(in) :: msg
    if (.not. cond) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg)
      stop 1
    end if
  end subroutine

  !> Assert that a logical condition is false.
  subroutine assert_false(cond, msg)
    logical, intent(in) :: cond
    character(*), intent(in) :: msg
    if (cond) then
      write(error_unit,*) "ASSERTION FAILED (expected false): ", trim(msg)
      stop 1
    end if
  end subroutine

  !> Assert that two integers are equal.
  subroutine assert_equal_int(a, b, msg)
    integer, intent(in) :: a, b
    character(*), intent(in) :: msg
    if (a /= b) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (got ", a, ", expected ", b, ")"
      stop 1
    end if
  end subroutine

  !> Assert that two integers are not equal.
  subroutine assert_not_equal_int(a, b, msg)
    integer, intent(in) :: a, b
    character(*), intent(in) :: msg
    if (a == b) then
      write(error_unit,*) "ASSERTION FAILED (should not be equal): ", trim(msg)
      stop 1
    end if
  end subroutine

  !> Assert that two real numbers are equal within a tolerance.
  subroutine assert_equal_real(a, b, tol, msg)
    real(real64), intent(in) :: a, b, tol
    character(*), intent(in) :: msg
    if (abs(a-b) > tol) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (got ", a, ", expected ", b, ", tol=", tol, ")"
      stop 1
    end if
  end subroutine

  !> Assert that two real numbers are not equal within a tolerance.
  subroutine assert_not_equal_real(a, b, tol, msg)
    real(real64), intent(in) :: a, b, tol
    character(*), intent(in) :: msg
    if (abs(a-b) <= tol) then
      write(error_unit,*) "ASSERTION FAILED (should not be equal): ", trim(msg)
      stop 1
    end if
  end subroutine

  !> Assert that two integer arrays are equal.
  subroutine assert_equal_array_int(a, b, n, msg)
    integer, intent(in) :: a(n), b(n), n
    character(*), intent(in) :: msg
    if (any(a /= b)) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (integer arrays differ)"
      stop 1
    end if
  end subroutine

  !> Assert that two real arrays are equal within a tolerance.
  subroutine assert_equal_array_real(a, b, n, tol, msg)
    real(real64), intent(in) :: a(n), b(n), tol
    integer, intent(in) :: n
    character(*), intent(in) :: msg
    if (any(abs(a-b) > tol)) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (real arrays differ, tol=", tol, ")"
      stop 1
    end if
  end subroutine

  !> Assert that a real array contains no NaN values.
  subroutine assert_no_nan_real(a, n, msg)
    real(real64), intent(in) :: a(n)
    integer, intent(in) :: n
    character(*), intent(in) :: msg
    integer :: i
    
    do i = 1, n
      if (ieee_is_nan(a(i))) then
        write(error_unit,*) "ASSERTION FAILED: NaN detected - ", trim(msg)
        stop 1
      end if
    end do
  end subroutine

  !> Assert that a real array contains no Inf values.
  subroutine assert_no_inf_real(a, n, msg)
    real(real64), intent(in) :: a(n)
    integer, intent(in) :: n
    character(*), intent(in) :: msg
    integer :: i
    do i = 1, n
      if (abs(a(i)) > huge(1.0_real64)) then
        write(error_unit,*) "ASSERTION FAILED: Inf detected - ", trim(msg)
        stop 1
      end if
    end do
  end subroutine

  !> Assert that a real value is within a given range [minval, maxval].
  subroutine assert_in_range_real(a, minval, maxval, msg)
    real(real64), intent(in) :: a, minval, maxval
    character(*), intent(in) :: msg
    if (a < minval .or. a > maxval) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (value ", a, " not in [", minval, ",", maxval, "])"
      stop 1
    end if
  end subroutine

  !> Assert that an integer array contains a given value.
  subroutine assert_contains_int(arr, n, val, msg)
    integer, intent(in) :: arr(n), n, val
    character(*), intent(in) :: msg
    if (.not. any(arr == val)) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (value ", val, " not found)"
      stop 1
    end if
  end subroutine

  !> Assert that an integer array is sorted in non-decreasing order.
  subroutine assert_sorted_int(arr, n, msg)
    integer, intent(in) :: arr(n), n
    character(*), intent(in) :: msg
    integer :: i
    do i = 2, n
      if (arr(i) < arr(i-1)) then
        write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (not sorted at position ", i, ")"
        stop 1
      end if
    end do
  end subroutine

  !> Assert that a real array is sorted in non-decreasing order.
  subroutine assert_sorted_real(arr, n, msg)
    real(real64), intent(in) :: arr(n)
    integer, intent(in) :: n
    character(*), intent(in) :: msg
    integer :: i
    do i = 2, n
      if (arr(i) < arr(i-1)) then
        write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (not sorted at position ", i, ")"
        stop 1
      end if
    end do
  end subroutine

  !> Assert that two arrays have the same shape (1D only).
  subroutine assert_same_shape(n1, n2, msg)
    integer, intent(in) :: n1, n2
    character(*), intent(in) :: msg
    if (n1 /= n2) then
        write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (shapes differ: ", n1, " vs ", n2, ")"
        stop 1
    end if
  end subroutine

  !> Assert that two strings are equal.
  subroutine assert_string_equal(a, b, msg)
    character(*), intent(in) :: a, b, msg
    if (trim(a) /= trim(b)) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (got '"//trim(a)//"', expected '"//trim(b)//"')"
      stop 1
    end if
  end subroutine

  !> Assert that string a contains string b.
  subroutine assert_string_contains(a, b, msg)
    character(*), intent(in) :: a, b, msg
    if (index(a, b) == 0) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (substring '"//trim(b)//"' not found in '"//trim(a)//"')"
      stop 1
    end if
  end subroutine

  !> Assert that two real arrays are close within relative and absolute tolerance.
  subroutine assert_allclose_array_real(a, b, n, rtol, atol, msg)
    real(real64), intent(in) :: a(n), b(n), rtol, atol
    integer, intent(in) :: n
    character(*), intent(in) :: msg
    integer :: i
    do i = 1, n
      if (abs(a(i) - b(i)) > atol + rtol * abs(b(i))) then
        write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (arrays differ at ", i, ")"
        stop 1
      end if
    end do
  end subroutine

  !> Assert that the sum of an array equals an expected value.
  subroutine assert_sum_equal(arr, n, expected, msg)
    real(real64), intent(in) :: arr(n), expected
    integer, intent(in) :: n
    character(*), intent(in) :: msg
    real(real64) :: s
    s = sum(arr)
    if (abs(s - expected) > 1e-12_real64) then
      write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (sum=", s, ", expected=", expected, ")"
      stop 1
    end if
  end subroutine

  !> Assert that all elements in an integer array are unique.
  subroutine assert_unique_int(arr, n, msg)
    integer, intent(in) :: arr(n), n
    character(*), intent(in) :: msg
    integer :: i, j
    do i = 1, n-1
      do j = i+1, n
        if (arr(i) == arr(j)) then
          write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (duplicate value ", arr(i), " at positions ", i, " and ", j, ")"
          stop 1
        end if
      end do
    end do
  end subroutine

  !> Assert that an integer array is a permutation of 1..n.
  subroutine assert_permutation(arr, n, msg)
    integer, intent(in) :: arr(n), n
    character(*), intent(in) :: msg
    integer :: i
    logical :: found(n)
    found = .false.
    do i = 1, n
      if (arr(i) < 1 .or. arr(i) > n) then
        write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (value out of range: ", arr(i), ")"
        stop 1
      end if
      if (found(arr(i))) then
        write(error_unit,*) "ASSERTION FAILED: ", trim(msg), " (duplicate value: ", arr(i), ")"
        stop 1
      end if
      found(arr(i)) = .true.
    end do
  end subroutine

end module asserts