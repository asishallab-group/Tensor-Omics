!> Unit test suite for root_mean_sq_normalization routine.
module mod_test_root_mean_sq_normalization
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use tox_normalization
  use test_suite
  use tox_errors
  implicit none
  public

contains

  !> Get array of all available tests.
  function get_all_tests_root_mean_sq_normalization() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(13))
    
    all_tests(1) = test_case("test_root_mean_sq_normalization_basic", test_root_mean_sq_normalization_basic)
    all_tests(2) = test_case("test_root_mean_sq_normalization_constant_rows", test_root_mean_sq_normalization_constant_rows)
    all_tests(3) = test_case("test_root_mean_sq_normalization_large_numbers", test_root_mean_sq_normalization_large_numbers)
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
  end function get_all_tests_root_mean_sq_normalization

  !> Test that root_mean_sq_normalization normalizes values correctly.
  subroutine test_root_mean_sq_normalization_basic()
    real(real64), dimension(2,2) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i_tissue, i_gene, ierr

    mat = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [2,2])
    call root_mean_sq_normalization(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")

    do i_gene = 1, 2
      std_dev(i_gene) = sqrt((mat(1,i_gene)**2 + mat(2,i_gene)**2) / 2.0d0)
      do i_tissue = 1, 2
        expected(i_tissue, i_gene) = mat(i_tissue, i_gene) / std_dev(i_gene)
      end do
    end do
    call assert_equal_array_real(result, expected, 4, 1d-12, "root_mean_sq_normalization: basic normalization failed")
  end subroutine test_root_mean_sq_normalization_basic

  !> Test that root_mean_sq_normalization handles constant rows (should normalize to 1).
  subroutine test_root_mean_sq_normalization_constant_rows()
    real(real64), dimension(2,2) :: mat, result, expected
    integer(int32) :: ierr

    mat = reshape([5.0d0, 5.0d0, 5.0d0, 5.0d0], [2,2])
    call root_mean_sq_normalization(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")

    expected = 1.0d0

    call assert_true(all(result == 1.0d0), "root_mean_sq_normalization: not all values are 1 for constant rows")
    call assert_no_nan_real(result, 4, "root_mean_sq_normalization: NaN in result for constant rows")
  end subroutine test_root_mean_sq_normalization_constant_rows

  !> Test that root_mean_sq_normalization normalizes large numbers properly.
  subroutine test_root_mean_sq_normalization_large_numbers()
    real(real64), dimension(2,2) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i_tissue, i_gene, ierr

    mat = reshape([1e6, 2e6, 1e6, 2e6], [2,2])
    call root_mean_sq_normalization(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")

    do i_gene = 1, 2
      std_dev(i_gene) = sqrt((mat(1, i_gene)**2 + mat(2, i_gene)**2) / 2.0d0)
      do i_tissue = 1, 2
        expected(i_tissue, i_gene) = mat(i_tissue, i_gene) / std_dev(i_gene)
      end do
    end do

    call assert_equal_array_real(result, expected, 4, 1d-12, "root_mean_sq_normalization: large numbers normalization failed")
    call assert_no_nan_real(result, 4, "root_mean_sq_normalization: NaN in result for large numbers")
    call assert_true(all(ieee_is_finite(result)), "root_mean_sq_normalization: Inf in result for large numbers")
  end subroutine test_root_mean_sq_normalization_large_numbers

  !> Test normalization of the identity matrix.
  subroutine test_identity_matrix()
    real(real64), dimension(3,3) :: mat, result
    integer(int32) :: i_tissue, i_gene, ierr
    mat = 0.0d0
    do i_gene = 1, 3
      mat(i_gene,i_gene) = 1.0d0
    end do
    call root_mean_sq_normalization(3, 3, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    do i_gene = 1, 3
      call assert_in_range_real(sum(result(:, i_gene)**2)/3.0d0, 1d0-1d-12, 1d0+1d-12, "identity: RMS not 1")
      do i_tissue = 1, 3
        if (i_tissue /= i_gene) then 
          call assert_equal_real(result(i_tissue, i_gene), 0.0d0, 1d-12, "identity: off-diagonal not zero")
        else
          call assert_equal_real(result(i_tissue, i_gene), 1d0 / sqrt(1d0 / 3.0d0), 1d-12, "identity: diagonal value not number of tissues/genes")
        end if
      end do
    end do
  end subroutine test_identity_matrix

  !> Test normalization of rows with all zeros.
  subroutine test_zero_rows()
    real(real64), dimension(2,3) :: mat, result
    integer(int32) :: ierr
    mat = 0.0d0
    call root_mean_sq_normalization(2, 3, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    call assert_true(all(result == 0.0d0), "zero rows: not all zeros")
    call assert_no_nan_real(result, 6, "zero rows: NaN in result")
  end subroutine test_zero_rows

  !> Test normalization of rows with negative values.
  subroutine test_negative_rows()
    real(real64), dimension(2,3) :: mat, result, expected
    real(real64), dimension(3) :: std_dev
    integer(int32) :: i_tissue, i_gene, ierr
    mat = reshape([-2.0d0, -4.0d0, -6.0d0, -8.0d0, -10.0d0, -12.0d0], [2,3])
    call root_mean_sq_normalization(3, 2, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    do i_gene = 1, 3
      std_dev(i_gene) = sqrt(sum(mat(:,i_gene)**2)/2.0d0)
      do i_tissue = 1, 2
        expected(i_tissue, i_gene) = mat(i_tissue, i_gene)/std_dev(i_gene)
      end do
    end do
    call assert_equal_array_real(result, expected, 6, 1d-12, "negative rows: normalization failed")
  end subroutine test_negative_rows

  !> Test normalization of a large random matrix.
  subroutine test_large_random_matrix()
    integer(int32), parameter :: n_genes=20, n_tissues=30
    real(real64), dimension(n_tissues, n_genes) :: mat, result
    integer(int32) :: i_gene, ierr
    integer(int32) :: n_seed
    integer(int32), allocatable :: seed_array(:)
    call random_seed(size=n_seed)
    allocate(seed_array(n_seed))
    seed_array = 42
    call random_seed(put=seed_array)
    deallocate(seed_array)
    call random_number(mat)
    call root_mean_sq_normalization(n_genes, n_tissues, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    do i_gene = 1, n_genes
      call assert_in_range_real(sqrt(sum(result(:, i_gene)**2)/n_tissues), 1d0-1d-10, 1d0+1d-10, "large random: RMS not 1")
    end do
    call assert_no_nan_real(result, n_genes*n_tissues, "large random: NaN in result")
  end subroutine test_large_random_matrix

  !> Test normalization of rows with a single nonzero value.
  subroutine test_single_nonzero()
    real(real64), dimension(2,4) :: mat, result, expected
    integer(int32) :: i, ierr
    mat = 0.0d0
    mat(1,3) = 5.0d0
    mat(2,2) = -7.0d0
    call root_mean_sq_normalization(2, 4, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    do i = 1, 2
      expected(i,:) = mat(i,:) / sqrt(sum(mat(i,:)**2)/4.0d0)
    end do
    call assert_equal_array_real(result, expected, 8, 1d-12, "single nonzero: normalization failed")
  end subroutine test_single_nonzero

  !> Test normalization with very small and very large values.
  subroutine test_small_large_values()
    real(real64), dimension(2,2) :: mat, result, expected
    real(real64), dimension(2) :: std_dev
    integer(int32) :: i_tissue, i_gene, ierr
    mat = reshape([1e-10, 1e10, 1e-10, 1e10], [2,2])
    call root_mean_sq_normalization(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    do i_gene = 1, 2
      std_dev(i_gene) = sqrt(sum(mat(:, i_gene)**2)/2.0d0)
      do i_tissue = 1, 2
        expected(i_tissue, i_gene) = mat(i_tissue, i_gene)/std_dev(i_gene)
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
    call root_mean_sq_normalization(2, 2, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    call assert_true(all(ieee_is_finite(result)), "root_mean_sq_normalization: output contains NaN/Inf unexpectedly")
  end subroutine test_nan_inf_input

  !> Test normalization of a single row and a single column matrix.
  subroutine test_single_row_col()
    real(real64), dimension(4, 1) :: mat1, result1, expected1
    real(real64), dimension(1, 4) :: mat2, result2
    real(real64) :: std_dev
    integer(int32) :: i_tissue, ierr
    mat1 = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [4, 1])
    call root_mean_sq_normalization(1, 4, mat1, result1, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    std_dev = sqrt(sum(mat1(:, 1)**2)/4.0d0)
    do i_tissue = 1, 4
      expected1(i_tissue, 1) = mat1(i_tissue, 1)/std_dev
    end do
    call assert_equal_array_real(result1, expected1, 4, 1d-12, "single row: normalization failed")
    mat2 = reshape([2.0d0, 4.0d0, 6.0d0, 8.0d0], [1, 4])
    call root_mean_sq_normalization(4, 1, mat2, result2, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    call assert_true(all(abs(result2) == 1.0d0), "single col: normalization failed")
  end subroutine test_single_row_col

  !> Test normalization of an empty matrix.
  subroutine test_empty_matrix()
    real(real64), allocatable :: mat(:,:), result(:,:)
    integer(int32) :: ierr
    allocate(mat(1,1), result(1,1))
    call root_mean_sq_normalization(0, 0, mat, result, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "root_mean_sq_normalization returned error")
    ! No assertion needed: just check no crash
  end subroutine test_empty_matrix

  !> Test normalization of symmetric rows.
  subroutine test_symmetric_rows()
    real(real64), dimension(3,2) :: mat, result
    integer(int32) :: i_tissue, ierr
    mat(:, 1) = [1.0d0, 2.0d0, 3.0d0]
    mat(:, 2) = [2.0d0, 4.0d0, 6.0d0]
    call root_mean_sq_normalization(2, 3, mat, result, ierr)
    call assert_equal_int(ierr, ERR_OK, "root_mean_sq_normalization returned error")
    do i_tissue = 1, 3
      call assert_equal_real(result(i_tissue, 2), result(i_tissue, 1), 1d-12, "symmetric rows: not equal after normalization")
    end do
  end subroutine test_symmetric_rows

end module mod_test_root_mean_sq_normalization