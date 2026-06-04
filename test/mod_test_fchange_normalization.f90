! filepath: test/mod_test_calc_fchange.f90
!> Unit test suite for calc_fchange routine.
module mod_test_calc_fchange
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public


contains

  !> Get array of all available tests.
  function get_all_tests_calc_fchange() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(9))
    
    all_tests(1) = test_case("test_calc_fchange_basic_calculation", test_calc_fchange_basic_calculation)
    all_tests(2) = test_case("test_calc_fchange_single_pair", test_calc_fchange_single_pair)
    all_tests(3) = test_case("test_calc_fchange_multiple_pairs", test_calc_fchange_multiple_pairs)
    all_tests(4) = test_case("test_calc_fchange_negative_values", test_calc_fchange_negative_values)
    all_tests(5) = test_case("test_calc_fchange_zero_values", test_calc_fchange_zero_values)
    all_tests(6) = test_case("test_calc_fchange_large_values", test_calc_fchange_large_values)
    all_tests(7) = test_case("test_calc_fchange_identical_values", test_calc_fchange_identical_values)
    all_tests(8) = test_case("test_calc_fchange_mixed_values", test_calc_fchange_mixed_values)
    all_tests(9) = test_case("test_calc_fchange_empty_matrix", test_calc_fchange_empty_matrix)
  end function get_all_tests_calc_fchange

  !> Test basic fold change calculation.
  subroutine test_calc_fchange_basic_calculation()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(2, 2) :: i_matrix  ! 2 genes × 2 samples
    real(real64), dimension(1, 2) :: o_matrix, expected_matrix  ! 2 genes × 1 pair
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    ! Input matrix (column-major): [1,2; 4,8] → geneA=[1,4], geneB=[2,8]
    i_matrix(:, 1) = [1.0d0, 4.0d0]
    i_matrix(:, 2) = [2.0d0, 8.0d0]
    
    ! Control column = 1 (muscle_dietM), Condition column = 2 (muscle_dietP)
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_basic_calculation

  !> Test with single gene, single pair.
  subroutine test_calc_fchange_single_pair()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(2, 1) :: i_matrix  ! 1 gene × 2 samples
    real(real64), dimension(1, 1) :: o_matrix, expected_matrix  ! 1 gene × 1 pair
    
    n_genes = 1; n_cols = 2; n_pairs = 1
    i_matrix(:, 1) = [5.0d0, 15.0d0]  ! gene1: control=5, condition=15
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_single_pair

  !> Test with multiple condition-control pairs.
  subroutine test_calc_fchange_multiple_pairs()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(3) :: control_cols, cond_cols
    real(real64), dimension(6, 2) :: i_matrix
    real(real64), dimension(3, 2) :: o_matrix, expected_matrix
    
    n_genes = 2; n_cols = 6; n_pairs = 3
    ! 2 genes, 6 samples (3 control-condition pairs)
    i_matrix(:, 1) = [1.0d0, 3.0d0, 5.0d0, 7.0d0, 9.0d0, 11.0d0]
    i_matrix(:, 2) = [2.0d0, 4.0d0, 6.0d0, 8.0d0, 10.0d0, 12.0d0]
    
    control_cols = [1, 3, 5]  ! Control samples
    cond_cols = [2, 4, 6]     ! Condition samples
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_multiple_pairs

  !> Test with negative expression values.
  subroutine test_calc_fchange_negative_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(2) :: control_cols, cond_cols
    real(real64), dimension(4, 2) :: i_matrix  ! 2 genes × 4 samples
    real(real64), dimension(2, 2) :: o_matrix, expected_matrix  ! 2 genes × 2 pairs
    
    n_genes = 2; n_cols = 4; n_pairs = 2
    i_matrix(:, 1) = [-1.0d0, 3.0d0, -5.0d0, 7.0d0]
    i_matrix(:, 2) = [-2.0d0, 4.0d0, -6.0d0, 8.0d0]
    
    control_cols = [1, 3]  ! Columns with negative values
    cond_cols = [2, 4]     ! Columns with positive values
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_negative_values

  !> Test with zero expression values.
  subroutine test_calc_fchange_zero_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(2, 2) :: i_matrix  ! 2 genes × 2 samples
    real(real64), dimension(1, 2) :: o_matrix, expected_matrix  ! 2 genes × 1 pair
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    i_matrix(:, 1) = [0.0d0, 5.0d0]  ! Control=0, Condition>0
    i_matrix(:, 2) = [0.0d0, 10.0d0]  ! Control=0, Condition>0
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_zero_values

  !> Test with large expression values for numerical stability.
  subroutine test_calc_fchange_large_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(2, 2) :: i_matrix  ! 2 genes × 2 samples
    real(real64), dimension(1, 2) :: o_matrix, expected_matrix  ! 2 genes × 1 pair
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    i_matrix(:, 1) = [1d6, 1d9]  ! Very large values
    i_matrix(:, 2) = [2d6, 2d9]  ! Very large values
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_large_values

  !> Test with identical control and condition values (zero fold change).
  subroutine test_calc_fchange_identical_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(2) :: control_cols, cond_cols
    real(real64), dimension(4, 2) :: i_matrix  ! 2 genes × 4 samples
    real(real64), dimension(2, 2) :: o_matrix, expected_matrix  ! 2 genes × 2 pairs
    
    n_genes = 2; n_cols = 4; n_pairs = 2
    i_matrix(:, 1) = [5.0d0, 5.0d0, 7.0d0, 7.0d0]
    i_matrix(:, 2) = [10.0d0, 10.0d0, 14.0d0, 14.0d0]
    
    control_cols = [1, 3]  ! Identical to condition columns
    cond_cols = [2, 4]
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_identical_values

  !> Test with mixed positive, negative, and zero values.
  subroutine test_calc_fchange_mixed_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(3) :: control_cols, cond_cols
    real(real64), dimension(6, 3) :: i_matrix  ! 3 genes × 6 samples
    real(real64), dimension(3, 3) :: o_matrix, expected_matrix   ! 3 genes × 3 pairs
    
    n_genes = 3; n_cols = 6; n_pairs = 3
    i_matrix = reshape([-1.0d0, 0.0d0, 5.0d0,   & ! Control samples (cols 1,3,5)
                2.0d0, -3.0d0, 0.0d0,   & ! Condition samples (cols 2,4,6)
                10.0d0, -5.0d0, 2.0d0,  &
                -8.0d0, 7.0d0, -1.0d0,  &
                0.0d0, 0.0d0, 0.0d0,    &
                15.0d0, -10.0d0, 3.0d0], shape(i_matrix))
    
    control_cols = [1, 3, 5]
    cond_cols = [2, 4, 6]
    
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_fchange returned error")

    expected_matrix = i_matrix(cond_cols, :) - i_matrix(control_cols, :)

    call assert_equal_array_real(o_matrix, expected_matrix, size(o_matrix, kind=int32), 1d-12, &
                            "test_calc_fchange_multiple_pairs: output doesn't match expected")
  end subroutine test_calc_fchange_mixed_values

  !> Test with empty input matrix.
  subroutine test_calc_fchange_empty_matrix()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(0) :: control_cols, cond_cols
    real(real64), dimension(0) :: i_matrix, o_matrix
    n_genes = 0; n_cols = 0; n_pairs = 0
    call calc_fchange(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "calc_fchange should return error for empty input")
  end subroutine test_calc_fchange_empty_matrix

end module mod_test_calc_fchange