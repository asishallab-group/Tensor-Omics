! filepath: test/mod_test_calc_fchange.f90
!> Unit test suite for calc_fchange routine.
module mod_test_calc_fchange
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
    type(test_case) :: all_tests(11)
    
    all_tests(1) = test_case("test_calc_fchange_basic_calculation", test_calc_fchange_basic_calculation)
    all_tests(2) = test_case("test_calc_fchange_result_correctness", test_calc_fchange_result_correctness)
    all_tests(3) = test_case("test_calc_fchange_preserves_dimensions", test_calc_fchange_preserves_dimensions)
    all_tests(4) = test_case("test_calc_fchange_single_pair", test_calc_fchange_single_pair)
    all_tests(5) = test_case("test_calc_fchange_multiple_pairs", test_calc_fchange_multiple_pairs)
    all_tests(6) = test_case("test_calc_fchange_negative_values", test_calc_fchange_negative_values)
    all_tests(7) = test_case("test_calc_fchange_zero_values", test_calc_fchange_zero_values)
    all_tests(8) = test_case("test_calc_fchange_large_values", test_calc_fchange_large_values)
    all_tests(9) = test_case("test_calc_fchange_identical_values", test_calc_fchange_identical_values)
    all_tests(10) = test_case("test_calc_fchange_mixed_values", test_calc_fchange_mixed_values)
    all_tests(11) = test_case("test_calc_fchange_empty_matrix", test_calc_fchange_empty_matrix)
  end function get_all_tests

  !> Run all calc_fchange tests.
  subroutine run_all_tests_calc_fchange()
    type(test_case) :: all_tests(11)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All calc_fchange tests passed successfully."
  end subroutine run_all_tests_calc_fchange

  !> Run specific calc_fchange tests by name.
  subroutine run_named_tests_calc_fchange(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(11)
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
  end subroutine run_named_tests_calc_fchange

  !> Test basic fold change calculation.
  subroutine test_calc_fchange_basic_calculation()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(4) :: i_matrix  ! 2 genes × 2 samples
    real(real64), dimension(2) :: o_matrix  ! 2 genes × 1 pair
    real(real64), dimension(2) :: expected_matrix
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    ! Input matrix (column-major): [1,2; 4,8] → geneA=[1,4], geneB=[2,8]
    i_matrix = [1.0d0, 2.0d0, 4.0d0, 8.0d0]
    
    ! Control column = 1 (muscle_dietM), Condition column = 2 (muscle_dietP)
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    ! Expected fold changes: geneA = 4-1 = 3, geneB = 8-2 = 6
    expected_matrix = [3.0d0, 6.0d0]
    
    call assert_equal_array_real(o_matrix, expected_matrix, 2, 1d-12, &
                            "test_calc_fchange_basic_calculation: fold change calculation incorrect")
    call assert_no_nan_real(o_matrix, 2, "test_calc_fchange_basic_calculation: NaN in result")
  end subroutine test_calc_fchange_basic_calculation

  !> Test result correctness with specific expected values (from R test).
  subroutine test_calc_fchange_result_correctness()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(4) :: i_matrix
    real(real64), dimension(2) :: o_matrix
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    ! Same data as basic test but focusing on correctness
    i_matrix = [1.0d0, 2.0d0, 4.0d0, 8.0d0]
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    ! Check specific values as in R test
    call assert_equal_real(o_matrix(1), 3.0d0, 1d-12, &
                      "test_calc_fchange_result_correctness: geneA fold change incorrect")
    call assert_equal_real(o_matrix(2), 6.0d0, 1d-12, &
                      "test_calc_fchange_result_correctness: geneB fold change incorrect")
  end subroutine test_calc_fchange_result_correctness

  !> Test that output dimensions are correct.
  subroutine test_calc_fchange_preserves_dimensions()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(2) :: control_cols, cond_cols
    real(real64), dimension(12) :: i_matrix  ! 3 genes × 4 samples
    real(real64), dimension(6) :: o_matrix   ! 3 genes × 2 pairs
    
    n_genes = 3; n_cols = 4; n_pairs = 2
    i_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, &
                7.0d0, 8.0d0, 9.0d0, 10.0d0, 11.0d0, 12.0d0]
    
    control_cols = [1, 3]  ! Control columns
    cond_cols = [2, 4]     ! Condition columns
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    ! Output should have n_genes * n_pairs elements
    call assert_equal_int(size(o_matrix), n_genes * n_pairs, &
                     "test_calc_fchange_preserves_dimensions: output size incorrect")
    call assert_no_nan_real(o_matrix, n_genes * n_pairs, &
                       "test_calc_fchange_preserves_dimensions: NaN in result")
  end subroutine test_calc_fchange_preserves_dimensions

  !> Test with single gene, single pair.
  subroutine test_calc_fchange_single_pair()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(2) :: i_matrix  ! 1 gene × 2 samples
    real(real64), dimension(1) :: o_matrix  ! 1 gene × 1 pair
    
    n_genes = 1; n_cols = 2; n_pairs = 1
    i_matrix = [5.0d0, 15.0d0]  ! gene1: control=5, condition=15
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    ! Expected: 15 - 5 = 10
    call assert_equal_real(o_matrix(1), 10.0d0, 1d-12, &
                      "test_calc_fchange_single_pair: single pair calculation incorrect")
  end subroutine test_calc_fchange_single_pair

  !> Test with multiple condition-control pairs.
  subroutine test_calc_fchange_multiple_pairs()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(3) :: control_cols, cond_cols
    real(real64), dimension(12) :: i_matrix  ! 2 genes × 6 samples
    real(real64), dimension(6) :: o_matrix   ! 2 genes × 3 pairs
    
    n_genes = 2; n_cols = 6; n_pairs = 3
    ! 2 genes, 6 samples (3 control-condition pairs)
    i_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, &  ! samples 1,2,3
                7.0d0, 8.0d0, 9.0d0, 10.0d0, 11.0d0, 12.0d0]  ! samples 4,5,6
    
    control_cols = [1, 3, 5]  ! Control samples
    cond_cols = [2, 4, 6]     ! Condition samples
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    call assert_no_nan_real(o_matrix, 6, "test_calc_fchange_multiple_pairs: NaN in result")
    
    ! Check first pair: gene1: 3-1=2, gene2: 4-2=2
    call assert_equal_real(o_matrix(1), 2.0d0, 1d-12, &
                      "test_calc_fchange_multiple_pairs: first pair gene1 incorrect")
    call assert_equal_real(o_matrix(2), 2.0d0, 1d-12, &
                      "test_calc_fchange_multiple_pairs: first pair gene2 incorrect")
  end subroutine test_calc_fchange_multiple_pairs

  !> Test with negative expression values.
  subroutine test_calc_fchange_negative_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(2) :: control_cols, cond_cols
    real(real64), dimension(8) :: i_matrix  ! 2 genes × 4 samples
    real(real64), dimension(4) :: o_matrix  ! 2 genes × 2 pairs
    
    n_genes = 2; n_cols = 4; n_pairs = 2
    i_matrix = [-1.0d0, -2.0d0, 3.0d0, 4.0d0, -5.0d0, -6.0d0, 7.0d0, 8.0d0]
    
    control_cols = [1, 3]  ! Columns with negative values
    cond_cols = [2, 4]     ! Columns with positive values
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    call assert_no_nan_real(o_matrix, 4, "test_calc_fchange_negative_values: NaN in result")
    
    ! First pair: gene1: 3-(-1)=4, gene2: 4-(-2)=6
    call assert_equal_real(o_matrix(1), 4.0d0, 1d-12, &
                      "test_calc_fchange_negative_values: negative control handling incorrect")
    call assert_equal_real(o_matrix(2), 6.0d0, 1d-12, &
                      "test_calc_fchange_negative_values: negative control handling incorrect")
  end subroutine test_calc_fchange_negative_values

  !> Test with zero expression values.
  subroutine test_calc_fchange_zero_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(4) :: i_matrix  ! 2 genes × 2 samples
    real(real64), dimension(2) :: o_matrix  ! 2 genes × 1 pair
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    i_matrix = [0.0d0, 0.0d0, 5.0d0, 10.0d0]  ! Control=0, Condition>0
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    ! Expected: gene1: 5-0=5, gene2: 10-0=10
    call assert_equal_real(o_matrix(1), 5.0d0, 1d-12, &
                      "test_calc_fchange_zero_values: zero control handling incorrect")
    call assert_equal_real(o_matrix(2), 10.0d0, 1d-12, &
                      "test_calc_fchange_zero_values: zero control handling incorrect")
  end subroutine test_calc_fchange_zero_values

  !> Test with large expression values for numerical stability.
  subroutine test_calc_fchange_large_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(1) :: control_cols, cond_cols
    real(real64), dimension(4) :: i_matrix  ! 2 genes × 2 samples
    real(real64), dimension(2) :: o_matrix  ! 2 genes × 1 pair
    
    n_genes = 2; n_cols = 2; n_pairs = 1
    i_matrix = [1d6, 2d6, 1d9, 2d9]  ! Very large values
    control_cols = [1]
    cond_cols = [2]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    call assert_no_nan_real(o_matrix, 2, "test_calc_fchange_large_values: NaN in result")
    
    ! Expected: gene1: 1e9-1e6 ≈ 1e9, gene2: 2e9-2e6 ≈ 2e9
    call assert_in_range_real(o_matrix(1), 9d8, 1.1d9, &
                         "test_calc_fchange_large_values: large value calculation incorrect")
    call assert_in_range_real(o_matrix(2), 1.9d9, 2.1d9, &
                         "test_calc_fchange_large_values: large value calculation incorrect")
  end subroutine test_calc_fchange_large_values

  !> Test with identical control and condition values (zero fold change).
  subroutine test_calc_fchange_identical_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(2) :: control_cols, cond_cols
    real(real64), dimension(8) :: i_matrix  ! 2 genes × 4 samples
    real(real64), dimension(4) :: o_matrix  ! 2 genes × 2 pairs
    
    n_genes = 2; n_cols = 4; n_pairs = 2
    i_matrix = [5.0d0, 10.0d0, 5.0d0, 10.0d0, 7.0d0, 14.0d0, 7.0d0, 14.0d0]
    
    control_cols = [1, 3]  ! Identical to condition columns
    cond_cols = [2, 4]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    ! All fold changes should be zero (identical values)
    call assert_equal_real(o_matrix(1), 0.0d0, 1d-12, &
                      "test_calc_fchange_identical_values: identical values should give zero fold change")
    call assert_equal_real(o_matrix(2), 0.0d0, 1d-12, &
                      "test_calc_fchange_identical_values: identical values should give zero fold change")
    call assert_equal_real(o_matrix(3), 0.0d0, 1d-12, &
                      "test_calc_fchange_identical_values: identical values should give zero fold change")
    call assert_equal_real(o_matrix(4), 0.0d0, 1d-12, &
                      "test_calc_fchange_identical_values: identical values should give zero fold change")
  end subroutine test_calc_fchange_identical_values

  !> Test with mixed positive, negative, and zero values.
  subroutine test_calc_fchange_mixed_values()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(3) :: control_cols, cond_cols
    real(real64), dimension(18) :: i_matrix  ! 3 genes × 6 samples
    real(real64), dimension(9) :: o_matrix   ! 3 genes × 3 pairs
    
    n_genes = 3; n_cols = 6; n_pairs = 3
    i_matrix = [-1.0d0, 0.0d0, 5.0d0,   & ! Control samples (cols 1,3,5)
                2.0d0, -3.0d0, 0.0d0,   & ! Condition samples (cols 2,4,6)
                10.0d0, -5.0d0, 2.0d0,  &
                -8.0d0, 7.0d0, -1.0d0,  &
                0.0d0, 0.0d0, 0.0d0,    &
                15.0d0, -10.0d0, 3.0d0]
    
    control_cols = [1, 3, 5]
    cond_cols = [2, 4, 6]
    
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_fchange_r returned error")
    
    call assert_no_nan_real(o_matrix, 9, "test_calc_fchange_mixed_values: NaN in result")
    
    ! Check some specific mixed calculations
    ! Pair 1: gene1: 2-(-1)=3, gene2: -3-0=-3, gene3: 0-5=-5
    call assert_equal_real(o_matrix(1), 3.0d0, 1d-12, &
                      "test_calc_fchange_mixed_values: mixed calculation incorrect")
    call assert_equal_real(o_matrix(2), -3.0d0, 1d-12, &
                      "test_calc_fchange_mixed_values: mixed calculation incorrect")
    call assert_equal_real(o_matrix(3), -5.0d0, 1d-12, &
                      "test_calc_fchange_mixed_values: mixed calculation incorrect")
  end subroutine test_calc_fchange_mixed_values

  !> Test with empty input matrix.
  subroutine test_calc_fchange_empty_matrix()
    integer(int32) :: n_genes, n_cols, n_pairs, ierr
    integer(int32), dimension(0) :: control_cols, cond_cols
    real(real64), dimension(0) :: i_matrix, o_matrix
    n_genes = 0; n_cols = 0; n_pairs = 0
    call calc_fchange_r(n_genes, n_cols, n_pairs, control_cols, cond_cols, i_matrix, o_matrix, ierr)
    call assert_equal_int(ierr, 202, "calc_fchange_r should return error for empty input")
  end subroutine test_calc_fchange_empty_matrix

end module mod_test_calc_fchange