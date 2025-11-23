! filepath: test/mod_test_quantile_normalization.f90
!> Unit test suite for quantile_normalization routine.
module mod_test_quantile_normalization
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
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

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(13)
    
    all_tests(1) = test_case("test_qn_preserves_dimensions", test_qn_preserves_dimensions)
    all_tests(2) = test_case("test_qn_identical_rows", test_qn_identical_rows)
    all_tests(3) = test_case("test_qn_no_nans_and_standardizes", test_qn_no_nans_and_standardizes)
    all_tests(4) = test_case("test_qn_single_row", test_qn_single_row)
    all_tests(5) = test_case("test_qn_single_column", test_qn_single_column)
    all_tests(6) = test_case("test_qn_all_equal", test_qn_all_equal)
    all_tests(7) = test_case("test_qn_large_random", test_qn_large_random)
    all_tests(8) = test_case("test_qn_negative_values", test_qn_negative_values)
    all_tests(9) = test_case("test_qn_zero_matrix", test_qn_zero_matrix)
    all_tests(10) = test_case("test_qn_sorted_input", test_qn_sorted_input)
    all_tests(11) = test_case("test_qn_reverse_sorted", test_qn_reverse_sorted)
    all_tests(12) = test_case("test_qn_edge_cases", test_qn_edge_cases)
    all_tests(13) = test_case("test_qn_empty_matrix", test_qn_empty_matrix)
  end function get_all_tests

  !> Run all quantile_normalization tests.
  subroutine run_all_tests_quantile_normalization()
    type(test_case) :: all_tests(13)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All quantile_normalization tests passed successfully."
  end subroutine run_all_tests_quantile_normalization

  !> Run specific quantile_normalization tests by name.
  subroutine run_named_tests_quantile_normalization(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(13)
    integer(int32) :: i, j
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
  end subroutine run_named_tests_quantile_normalization

  !> Test that quantile normalization preserves matrix dimensions.
  subroutine test_qn_preserves_dimensions()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(2,3) :: mat, result
    real(real64) :: temp_col(2), rank_means(2)
    integer(int32) :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(:,1) = [0.1d0, 0.2d0]
    mat(:,2) = [0.3d0, 0.4d0]
    mat(:,3) = [0.5d0, 0.6d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_int(size(result,1), nrow, "qn_preserves_dimensions: row count not preserved")
    call assert_equal_int(size(result,2), ncol, "qn_preserves_dimensions: col count not preserved")
  end subroutine test_qn_preserves_dimensions

  !> Test that quantile normalization handles identical rows correctly.
  subroutine test_qn_identical_rows()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(2,3) :: mat, result, expected
    real(real64) :: temp_col(2), rank_means(2)
    integer(int32) :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(1,:) = [5.0d0, 5.0d0, 5.0d0]
    mat(2,:) = [5.0d0, 5.0d0, 5.0d0]
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, expected, 6, 1d-12, "qn_identical_rows: identical rows not preserved")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), "qn_identical_rows: result has non-finite values")
  end subroutine test_qn_identical_rows

  !> Test that quantile normalization produces finite values and standardizes distributions.
  subroutine test_qn_no_nans_and_standardizes()
    integer(int32) :: nrow, ncol, i, max_stack, ierr
    real(real64), dimension(2,3) :: mat, result
    real(real64), allocatable :: sorted_cols(:,:)
    real(real64) :: temp_col(2), rank_means(2)
    integer(int32) :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(:,1) = [2.0d0, 0.0d0]
    mat(:,2) = [5.0d0, 3.0d0]
    mat(:,3) = [7.0d0, 1.0d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_true(.not. any(reshape(isnan_mat(result), [size(result)])), "qn_no_nans_and_standardizes: result has NaN")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), &
                        "qn_no_nans_and_standardizes: result has non-finite values")
    call assert_equal_int(size(result,1), nrow, "qn_no_nans_and_standardizes: row count not preserved")
    call assert_equal_int(size(result,2), ncol, "qn_no_nans_and_standardizes: col count not preserved")
    allocate(sorted_cols(nrow, ncol))
    do i = 1, ncol
        sorted_cols(:,i) = sort_vec(result(:,i))
    end do
    do i = 2, ncol
        call assert_equal_array_real(sorted_cols(:,i), sorted_cols(:,1), nrow, 1d-6, &
                                    "qn_no_nans_and_standardizes: column distributions differ")
    end do
    deallocate(sorted_cols)
  end subroutine test_qn_no_nans_and_standardizes

  !> Test quantile normalization with single row matrix.
  subroutine test_qn_single_row()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(1,4) :: mat, result, expected
    real(real64) :: temp_col(1), rank_means(1)
    integer(int32) :: perm(1), stack_left(10), stack_right(10)
    real(real64) :: row_mean
    nrow = 1; ncol = 4; max_stack = 10
    mat(1,:) = [1.0d0, 2.0d0, 3.0d0, 4.0d0]
    row_mean = sum(mat(1,:)) / ncol
    expected(1,:) = row_mean
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, expected, 4, 1d-12, "qn_single_row: values should all equal row mean")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), "qn_single_row: result has non-finite values")
  end subroutine test_qn_single_row

  !> Test quantile normalization with single column matrix.
  subroutine test_qn_single_column()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(3,1) :: mat, result, expected
    real(real64) :: temp_col(3), rank_means(3)
    integer(int32) :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 1; max_stack = 10
    mat(:,1) = [7.0d0, 2.0d0, 5.0d0]
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, expected, 3, 1d-12, "qn_single_column: single column should remain unchanged")
  end subroutine test_qn_single_column

  !> Test quantile normalization with all equal values.
  subroutine test_qn_all_equal()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(4,4) :: mat, result, expected
    real(real64) :: temp_col(4), rank_means(4)
    integer(int32) :: perm(4), stack_left(10), stack_right(10)
    nrow = 4; ncol = 4; max_stack = 10
    mat = 3.14d0
    expected = mat
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, expected, 16, 1d-12, "qn_all_equal: all equal values not preserved")
  end subroutine test_qn_all_equal

  !> Test quantile normalization with large random matrix.
  subroutine test_qn_large_random()
    integer(int32), parameter :: nrow=10, ncol=10
    integer(int32) :: i, max_stack, ierr
    real(real64), dimension(nrow,ncol) :: mat, result
    real(real64), allocatable :: sorted_cols(:,:)
    real(real64) :: temp_col(nrow), rank_means(nrow)
    integer(int32) :: perm(nrow), stack_left(100), stack_right(100)
    integer(int32) :: n_seed
    integer(int32), allocatable :: seed_array(:)
    max_stack = 100
    call random_seed(size=n_seed)
    allocate(seed_array(n_seed))
    seed_array = 42
    call random_seed(put=seed_array)
    deallocate(seed_array)
    call random_number(mat)
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    allocate(sorted_cols(nrow, ncol))
    do i = 1, ncol
        sorted_cols(:,i) = sort_vec(result(:,i))
    end do
    do i = 2, ncol
        call assert_equal_array_real(sorted_cols(:,i), sorted_cols(:,1), nrow, 1d-6, &
                                "qn_large_random: random matrix column distributions differ")
    end do
    call assert_no_nan_real(result, nrow*ncol, "qn_large_random: NaN in result")
    deallocate(sorted_cols)
  end subroutine test_qn_large_random

  !> Test quantile normalization with negative values.
  subroutine test_qn_negative_values()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(2,3) :: mat, result
    real(real64) :: temp_col(2), rank_means(2)
    integer(int32) :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 3; max_stack = 10
    mat(:,1) = [-1.0d0, -2.0d0]
    mat(:,2) = [-3.0d0, -4.0d0]
    mat(:,3) = [-5.0d0, -6.0d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_no_nan_real(result, 6, "qn_negative_values: NaN in result")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), "qn_negative_values: non-finite values in result")
  end subroutine test_qn_negative_values

  !> Test quantile normalization with zero matrix.
  subroutine test_qn_zero_matrix()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(3,3) :: mat, result, expected
    real(real64) :: temp_col(3), rank_means(3)
    integer(int32) :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 3; max_stack = 10
    mat = 0.0d0
    expected = 0.0d0
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, expected, 9, 1d-12, "qn_zero_matrix: zero matrix not preserved")
  end subroutine test_qn_zero_matrix

  !> Test quantile normalization with pre-sorted input.
  subroutine test_qn_sorted_input()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(3,3) :: mat, result
    real(real64) :: temp_col(3), rank_means(3)
    integer(int32) :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 3; max_stack = 10
    mat(:,1) = [1.0d0, 2.0d0, 3.0d0]
    mat(:,2) = [1.0d0, 2.0d0, 3.0d0]
    mat(:,3) = [1.0d0, 2.0d0, 3.0d0]
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, mat, 9, 1d-12, "qn_sorted_input: sorted input not preserved")
  end subroutine test_qn_sorted_input

  !> Test quantile normalization with reverse sorted input.
  subroutine test_qn_reverse_sorted()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(3,3) :: mat, result, expected
    real(real64) :: temp_col(3), rank_means(3)
    integer(int32) :: perm(3), stack_left(10), stack_right(10)
    nrow = 3; ncol = 3; max_stack = 10
    mat(:,1) = [3.0d0, 3.0d0, 3.0d0]
    mat(:,2) = [2.0d0, 2.0d0, 2.0d0]
    mat(:,3) = [1.0d0, 1.0d0, 1.0d0]
    expected = 2.0d0
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_equal_array_real(result, expected, 9, 1d-12, "qn_reverse_sorted: reverse sorted input not normalized to mean")
  end subroutine test_qn_reverse_sorted

  !> Test edge cases for quantile normalization.
  subroutine test_qn_edge_cases()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(2,2) :: mat, result
    real(real64) :: temp_col(2), rank_means(2)
    integer(int32) :: perm(2), stack_left(10), stack_right(10)
    nrow = 2; ncol = 2; max_stack = 10
    mat = reshape([1e-10, 1e10, 1e-5, 1e5], [2,2])
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 0, "quantile_normalization_r returned error")
    call assert_no_nan_real(result, 4, "qn_edge_cases: NaN in result")
    call assert_true(all(reshape(isfinite_mat(result), [size(result)])), "qn_edge_cases: non-finite values in result")
  end subroutine test_qn_edge_cases

  !> Test with empty input matrix.
  subroutine test_qn_empty_matrix()
    integer(int32) :: nrow, ncol, max_stack, ierr
    real(real64), dimension(0,0) :: mat, result
    real(real64), dimension(0) :: temp_col, rank_means
    integer(int32), dimension(0) :: perm, stack_left, stack_right
    nrow = 0; ncol = 0; max_stack = 0
    call quantile_normalization_r(nrow, ncol, mat, result, temp_col, rank_means, perm, stack_left, stack_right, max_stack, ierr)
    call assert_equal_int(ierr, 202, "quantile_normalization_r should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_qn_empty_matrix

  !> Helper function to check if all values are finite (matrix version).
  function isfinite_mat(arr) result(mask)
    real(real64), intent(in) :: arr(:,:)
    logical :: mask(size(arr,1), size(arr,2))
    integer(int32) :: i, j
    do i = 1, size(arr,1)
      do j = 1, size(arr,2)
        mask(i,j) = abs(arr(i,j)) < huge(arr(i,j)) .and. arr(i,j) == arr(i,j)
      end do
    end do
  end function isfinite_mat

  !> Helper function to check for NaN values (matrix version).
  function isnan_mat(arr) result(mask)
    real(real64), intent(in) :: arr(:,:)
    logical :: mask(size(arr,1), size(arr,2))
    integer(int32) :: i, j
    do i = 1, size(arr,1)
      do j = 1, size(arr,2)
        mask(i,j) = ieee_is_nan(arr(i,j))
      end do
    end do
  end function isnan_mat

  !> Helper function to sort a vector (simple bubble sort).
  function sort_vec(vec) result(sorted)
    real(real64), intent(in) :: vec(:)
    real(real64) :: sorted(size(vec))
    integer(int32) :: i, j
    sorted = vec
    do i = 1, size(vec)-1
      do j = i+1, size(vec)
        if (sorted(j) < sorted(i)) then
          call swap(sorted(i), sorted(j))
        end if
      end do
    end do
  end function sort_vec

  !> Helper subroutine to swap two real values.
  subroutine swap(a, b)
    real(real64), intent(inout) :: a, b
    real(real64) :: tmp
    tmp = a
    a = b
    b = tmp
  end subroutine swap

end module mod_test_quantile_normalization