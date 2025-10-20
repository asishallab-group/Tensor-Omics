!> Unit tests for compute_relative_axis_contributions and wrappers
module mod_test_relative_axis_contributions
  use relative_axis_plane_tools
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic
  implicit none
  public 


  ! Abstract interface for all test procedures
  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  ! Type to hold test name and procedure pointer
  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case


contains

  !> Get array of all available tests for relative axis contributions
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(11)
    all_tests(1) = test_case("test_positive_vector", test_positive_vector)
    all_tests(2) = test_case("test_negative_vector", test_negative_vector)
    all_tests(3) = test_case("test_mixed_vector", test_mixed_vector)
    all_tests(4) = test_case("test_zero_vector", test_zero_vector)
    all_tests(5) = test_case("test_one_nonzero_axis", test_one_nonzero_axis)
    all_tests(6) = test_case("test_all_equal", test_all_equal)
    all_tests(7) = test_case("test_large_vector", test_large_vector)
    all_tests(8) = test_case("test_wrappers", test_wrappers)
    all_tests(9) = test_case("test_nan_vector", test_nan_vector)
    all_tests(10) = test_case("test_inf_vector", test_inf_vector)
    all_tests(11) = test_case("test_empty_vector", test_empty_vector)
  end function get_all_tests
  !> Run all relative axis contribution tests
  subroutine run_all_tests_relative_axis()
    type(test_case) :: all_tests(11)
    integer :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All relative axis contribution tests passed successfully."
  end subroutine run_all_tests_relative_axis

  !> Run specific relative axis tests by name
  subroutine run_named_tests_relative_axis(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(11)
    integer :: i, j
    logical :: found
    all_tests = get_all_tests()
    do i = 1, size(test_names)
      found = .false.
      do j = 1, size(all_tests)
        if (trim(test_names(i)) == trim(all_tests(j)%name)) then
          call all_tests(j)%test_proc()
          print *, trim(test_names(i)), " passed."
          found = .true.
          exit
        end if
      end do
      if (.not. found) then
        print *, "Unknown test: ", trim(test_names(i))
      end if
    end do
  end subroutine run_named_tests_relative_axis

  !> Test: vector with all positive values
  subroutine test_positive_vector()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n) = [1.0_real64, 2.0_real64, 3.0_real64]
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: positive vector')
    call assert_sum_equal(contrib, n, 1.0_real64, 'positive vector: sum')
    call assert_no_nan_real(contrib, n, 'positive vector: nan')
    call assert_no_inf_real(contrib, n, 'positive vector: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'positive vector: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'positive vector: max')
  end subroutine

  !> Test: vector with all negative values
  subroutine test_negative_vector()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n) = [-1.0_real64, -2.0_real64, -3.0_real64]
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: negative vector')
    call assert_sum_equal(contrib, n, 1.0_real64, 'negative vector: sum')
    call assert_no_nan_real(contrib, n, 'negative vector: nan')
    call assert_no_inf_real(contrib, n, 'negative vector: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'negative vector: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'negative vector: max')
  end subroutine

  !> Test: vector with mixed positive and negative values
  subroutine test_mixed_vector()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n) = [2.0_real64, -2.0_real64, 4.0_real64]
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: mixed vector')
    call assert_sum_equal(contrib, n, 1.0_real64, 'mixed vector: sum')
    call assert_no_nan_real(contrib, n, 'mixed vector: nan')
    call assert_no_inf_real(contrib, n, 'mixed vector: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'mixed vector: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'mixed vector: max')
  end subroutine

  !> Test: vector with all zeros (edge case)
  subroutine test_zero_vector()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n) = [0.0_real64, 0.0_real64, 0.0_real64]
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_true(ierr /= 0, 'ierr should be nonzero for zero vector')
    call assert_sum_equal(contrib, n, 0.0_real64, 'zero vector: sum')
    call assert_no_nan_real(contrib, n, 'zero vector: nan')
    call assert_no_inf_real(contrib, n, 'zero vector: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'zero vector: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'zero vector: max')
  end subroutine

  !> Test: vector with only one nonzero axis
  subroutine test_one_nonzero_axis()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n) = [0.0_real64, 5.0_real64, 0.0_real64]
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: one nonzero axis')
    call assert_sum_equal(contrib, n, 1.0_real64, 'one nonzero axis: sum')
    call assert_no_nan_real(contrib, n, 'one nonzero axis: nan')
    call assert_no_inf_real(contrib, n, 'one nonzero axis: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'one nonzero axis: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'one nonzero axis: max')
  end subroutine

  !> Test: vector with all axes equal
  subroutine test_all_equal()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n) = [2.0_real64, 2.0_real64, 2.0_real64]
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: all equal')
    call assert_sum_equal(contrib, n, 1.0_real64, 'all equal: sum')
    call assert_no_nan_real(contrib, n, 'all equal: nan')
    call assert_no_inf_real(contrib, n, 'all equal: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'all equal: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'all equal: max')
  end subroutine

  !> Test: large vector with all values equal (performance and correctness)
  subroutine test_large_vector()
    integer(int32), parameter :: n = 100
    real(real64) :: vec(n)
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    vec = 1.0_real64
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: large vector')
    call assert_sum_equal(contrib, n, 1.0_real64, 'large vector: sum')
    call assert_no_nan_real(contrib, n, 'large vector: nan')
    call assert_no_inf_real(contrib, n, 'large vector: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'large vector: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'large vector: max')
  end subroutine

  !> Test: wrapper subroutines for shift and expression vectors
  subroutine test_wrappers()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n), contrib(n)
    integer(int32) :: ierr
    vec = [1.0_real64, -2.0_real64, 3.0_real64]
    call relative_axes_changes_from_shift_vector(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: wrapper shift')
    call assert_sum_equal(contrib, n, 1.0_real64, 'wrapper shift: sum')
    call assert_no_nan_real(contrib, n, 'wrapper shift: nan')
    call assert_no_inf_real(contrib, n, 'wrapper shift: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'wrapper shift: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'wrapper shift: max')

    vec = [4.0_real64, 0.0_real64, 2.0_real64]
    call relative_axes_expression_from_expression_vector(vec, n, contrib, ierr)
    call assert_equal_int(ierr, 0, 'ierr should be 0 for valid input: wrapper expression')
    call assert_sum_equal(contrib, n, 1.0_real64, 'wrapper expression: sum')
    call assert_no_nan_real(contrib, n, 'wrapper expression: nan')
    call assert_no_inf_real(contrib, n, 'wrapper expression: inf')
    call assert_in_range_real(minval(contrib), 0.0_real64, 1.0_real64, 'wrapper expression: min')
    call assert_in_range_real(maxval(contrib), 0.0_real64, 1.0_real64, 'wrapper expression: max')
  end subroutine

  !> Test: vector containing NaN
  subroutine test_nan_vector()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n)
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    vec = [1.0_real64, 0.0_real64, ieee_value(1.0_real64, ieee_quiet_nan)]
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_true(ierr /= 0, 'ierr should be nonzero for NaN vector')
    call assert_no_nan_real(contrib, n, 'NaN vector: nan')
    call assert_no_inf_real(contrib, n, 'NaN vector: inf')
  end subroutine

  !> Test: vector containing Inf
  subroutine test_inf_vector()
    integer(int32), parameter :: n = 3
    real(real64) :: vec(n)
    real(real64) :: contrib(n)
    integer(int32) :: ierr
    vec = [1.0_real64, 0.0_real64, ieee_value(1.0_real64, ieee_positive_inf)]
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_true(ierr /= 0, 'ierr should be nonzero for Inf vector')
    call assert_no_nan_real(contrib, n, 'Inf vector: nan')
    call assert_no_inf_real(contrib, n, 'Inf vector: inf')
  end subroutine

  !> Test: n = 0 (empty vector)
  subroutine test_empty_vector()
    integer(int32), parameter :: n = 0
    real(real64), allocatable :: vec(:)
    real(real64), allocatable :: contrib(:)
    integer(int32) :: ierr
    allocate(vec(n), contrib(n))
    call compute_relative_axis_contributions(vec, n, contrib, ierr)
    call assert_true(ierr /= 0, 'ierr should be nonzero for n=0 (empty vector)')
  end subroutine

end module mod_test_relative_axis_contributions
