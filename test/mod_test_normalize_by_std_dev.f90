! filepath: test/mod_test_normalize_by_std_dev.f90
!> Unit test suite for normalize_by_std_dev routine.
module mod_test_normalize_by_std_dev
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
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
    
    all_tests(1) = test_case("test_normalize_by_std_dev_basic", test_normalize_by_std_dev_basic)
    all_tests(2) = test_case("test_normalize_by_std_dev_constant_rows", test_normalize_by_std_dev_constant_rows)
    all_tests(3) = test_case("test_normalize_by_std_dev_large_numbers", test_normalize_by_std_dev_large_numbers)
    all_tests(4) = test_case("test_identity_matrix", test_identity_matrix)
    all_tests(5) = test_case("test_zero_rows", test_zero_rows)
    all_tests(6) = test_case("test_negative_rows", test_negative_rows)
    all_tests(7) = test_case("test_large_random_matrix", test_large_random_matrix)
    all_tests(8) = test_case("test_single_nonzero", test_single_nonzero)
    all_tests(9) = test_case("test_small_large_values", test_small_large_values)
    all_tests(10) = test_case("test_nan_inf_input", test_nan_inf_input)
    all_tests(11) = test_case("test_single_row_col", test_single_row_col)
    all_tests(12) = test_case("test_empty_matrix", test_empty_matrix)
    all_tests(13) = test_case("test_symmetric_rows", test_symmetric_rows)
  end function get_all_tests

  !> Run all normalize_by_std_dev tests.
  subroutine run_all_tests_normalize_by_std_dev()
    type(test_case) :: all_tests(13)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All normalize_by_std_dev tests passed successfully."
  end subroutine run_all_tests_normalize_by_std_dev

  !> Run specific normalize_by_std_dev tests by name.
  subroutine run_named_tests_normalize_by_std_dev(test_names)
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
  end subroutine run_named_tests_normalize_by_std_dev

  !> Test that normalize_by_std_dev normalizes values correctly.
  subroutine test_normalize_by_std_dev_basic()
    real(real64), dimension(2,2) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i, j, ierr

    mat = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")

    do i = 1, 2
      std_dev(i) = sqrt((mat(i,1)**2 + mat(i,2)**2) / 2.0d0)
      do j = 1, 2
        expected(i,j) = mat(i,j) / std_dev(i)
      end do
    end do

    call assert_equal_array_real(result, expected, 4, 1d-12, "normalize_by_std_dev: basic normalization failed")
  end subroutine test_normalize_by_std_dev_basic

  !> Test that normalize_by_std_dev handles constant rows (should normalize to 1).
  subroutine test_normalize_by_std_dev_constant_rows()
    real(real64), dimension(2,2) :: mat, result, expected
    integer(int32) :: i, j, ierr

    mat = reshape([5.0d0, 5.0d0, 5.0d0, 5.0d0], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")

    expected = 1.0d0

    call assert_true(all(result == 1.0d0), "normalize_by_std_dev: not all values are 1 for constant rows")
    call assert_no_nan_real(result, 4, "normalize_by_std_dev: NaN in result for constant rows")
  end subroutine test_normalize_by_std_dev_constant_rows

  !> Test that normalize_by_std_dev normalizes large numbers properly.
  subroutine test_normalize_by_std_dev_large_numbers()
    real(real64), dimension(2,2) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i, j, ierr

    mat = reshape([1e6, 2e6, 1e6, 2e6], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")

    do i = 1, 2
      std_dev(i) = sqrt((mat(i,1)**2 + mat(i,2)**2) / 2.0d0)
      do j = 1, 2
        expected(i,j) = mat(i,j) / std_dev(i)
      end do
    end do

    call assert_equal_array_real(result, expected, 4, 1d-12, "normalize_by_std_dev: large numbers normalization failed")
    call assert_no_nan_real(result, 4, "normalize_by_std_dev: NaN in result for large numbers")
    call assert_true(all(isfinite_mat(result)), "normalize_by_std_dev: Inf in result for large numbers")
  end subroutine test_normalize_by_std_dev_large_numbers

  !> Test normalization of the identity matrix.
  subroutine test_identity_matrix()
    real(real64), dimension(3,3) :: mat, result
    integer(int32) :: i, j, ierr
    mat = 0.0d0
    do i = 1, 3
      mat(i,i) = 1.0d0
    end do
    call normalize_by_std_dev_r(3, 3, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    do i = 1, 3
      call assert_in_range_real(sum(result(i,:)**2)/3.0d0, 1d0-1d-12, 1d0+1d-12, "identity: RMS not 1")
      do j = 1, 3
        if (i /= j) call assert_equal_real(result(i,j), 0.0d0, 1d-12, "identity: off-diagonal not zero")
      end do
    end do
  end subroutine test_identity_matrix

  !> Test normalization of rows with all zeros.
  subroutine test_zero_rows()
    real(real64), dimension(2,3) :: mat, result
    integer(int32) :: ierr
    mat = 0.0d0
    call normalize_by_std_dev_r(2, 3, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    call assert_true(all(result == 0.0d0), "zero rows: not all zeros")
    call assert_no_nan_real(result, 6, "zero rows: NaN in result")
  end subroutine test_zero_rows

  !> Test normalization of rows with negative values.
  subroutine test_negative_rows()
    real(real64), dimension(2,3) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i, j, ierr
    mat = reshape([-2.0d0, -4.0d0, -6.0d0, -8.0d0, -10.0d0, -12.0d0], [2,3])
    call normalize_by_std_dev_r(2, 3, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    do i = 1, 2
      std_dev(i) = sqrt(sum(mat(i,:)**2)/3.0d0)
      do j = 1, 3
        expected(i,j) = mat(i,j)/std_dev(i)
      end do
    end do
    call assert_equal_array_real(result, expected, 6, 1d-12, "negative rows: normalization failed")
  end subroutine test_negative_rows

  !> Test normalization of a large random matrix.
  subroutine test_large_random_matrix()
    integer(int32), parameter :: nrow=20, ncol=30
    real(real64), dimension(nrow,ncol) :: mat, result
    integer(int32) :: i, ierr
    integer(int32) :: n_seed
    integer(int32), allocatable :: seed_array(:)
    call random_seed(size=n_seed)
    allocate(seed_array(n_seed))
    seed_array = 42
    call random_seed(put=seed_array)
    deallocate(seed_array)
    call random_number(mat)
    call normalize_by_std_dev_r(nrow, ncol, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    do i = 1, nrow
      call assert_in_range_real(sqrt(sum(result(i,:)**2)/ncol), 1d0-1d-10, 1d0+1d-10, "large random: RMS not 1")
    end do
    call assert_no_nan_real(result, nrow*ncol, "large random: NaN in result")
  end subroutine test_large_random_matrix

  !> Test normalization of rows with a single nonzero value.
  subroutine test_single_nonzero()
    real(real64), dimension(2,4) :: mat, result, expected
    integer(int32) :: i, ierr
    mat = 0.0d0
    mat(1,3) = 5.0d0
    mat(2,2) = -7.0d0
    call normalize_by_std_dev_r(2, 4, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    do i = 1, 2
      expected(i,:) = mat(i,:) / sqrt(sum(mat(i,:)**2)/4.0d0)
    end do
    call assert_equal_array_real(result, expected, 8, 1d-12, "single nonzero: normalization failed")
  end subroutine test_single_nonzero

  !> Test normalization with very small and very large values.
  subroutine test_small_large_values()
    real(real64), dimension(2,2) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i, j, ierr
    mat = reshape([1e-10, 1e10, 1e-10, 1e10], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    do i = 1, 2
      std_dev(i) = sqrt(sum(mat(i,:)**2)/2.0d0)
      do j = 1, 2
        expected(i,j) = mat(i,j)/std_dev(i)
      end do
    end do
    call assert_equal_array_real(result, expected, 4, 1d-10, "small/large values: normalization failed")
    call assert_no_nan_real(result, 4, "small/large values: NaN in result")
  end subroutine test_small_large_values

  !> Test normalization when input contains NaN or Inf.
  subroutine test_nan_inf_input()
    real(real64), dimension(2,2) :: mat, result
    integer(int32) :: ierr
    mat = reshape([1.0d0, 2.0d0, huge(1.0d0), 4.0d0], [2,2])
    call normalize_by_std_dev_r(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    call assert_true(all(isfinite_mat(result)), "normalize_by_std_dev: output contains NaN/Inf unexpectedly")
  end subroutine test_nan_inf_input

  !> Test normalization of a single row and a single column matrix.
  subroutine test_single_row_col()
    real(real64), dimension(1,4) :: mat1, result1, expected1
    real(real64), dimension(4,1) :: mat2, result2
    real(real64) :: std_dev
    integer(int32) :: j, ierr
    mat1 = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [1,4])
    call normalize_by_std_dev_r(1, 4, mat1, result1, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    std_dev = sqrt(sum(mat1(1,:)**2)/4.0d0)
    do j = 1, 4
      expected1(1,j) = mat1(1,j)/std_dev
    end do
    call assert_equal_array_real(result1, expected1, 4, 1d-12, "single row: normalization failed")
    mat2 = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [4,1])
    call normalize_by_std_dev_r(4, 1, mat2, result2, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    call assert_true(all(abs(result2) == 1.0d0), "single col: normalization failed")
  end subroutine test_single_row_col

  !> Test normalization of an empty matrix.
  subroutine test_empty_matrix()
    real(real64), allocatable :: mat(:,:), result(:,:)
    integer(int32) :: ierr
    allocate(mat(0,0), result(0,0))
    call normalize_by_std_dev_r(0, 0, mat, result, ierr)
    call assert_equal_int(ierr, 202, "normalize_by_std_dev_r returned error")
    ! No assertion needed: just check no crash
  end subroutine test_empty_matrix

  !> Test normalization of symmetric rows.
  subroutine test_symmetric_rows()
    real(real64), dimension(2,3) :: mat, result
    integer(int32) :: j, ierr
    mat(1,:) = [1.0d0, 2.0d0, 3.0d0]
    mat(2,:) = [2.0d0, 4.0d0, 6.0d0]
    call normalize_by_std_dev_r(2, 3, mat, result, ierr)
    call assert_equal_int(ierr, 0, "normalize_by_std_dev_r returned error")
    do j = 1, 3
      call assert_equal_real(result(2,j), result(1,j), 1d-12, "symmetric rows: not equal after normalization")
    end do
  end subroutine test_symmetric_rows

  !> Helper function to check if all values are finite (matrix version).
  function isfinite_mat(arr) result(mask)
    real(real64), intent(in) :: arr(:,:)
    logical :: mask(size(arr,1), size(arr,2))
    integer(int32) :: i, j
    do i = 1, size(arr,1)
      do j = 1, size(arr,2)
        mask(i,j) = abs(arr(i,j)) < huge(1.0d0)
      end do
    end do
  end function isfinite_mat

end module mod_test_normalize_by_std_dev