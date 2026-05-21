!> Unit test suite for calc_tiss_avg routine.
module mod_test_calc_tiss_avg
  use asserts
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_normalization
  use test_suite, only: test_case
  use tox_errors
  implicit none
  public

 
contains

  !> Get array of all available tests.
  function get_all_tests_tiss_avg() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(8))
    
    all_tests(1) = test_case("test_calc_tiss_avg_three_tissues", test_calc_tiss_avg_three_tissues)
    all_tests(2) = test_case("test_calc_tiss_avg_single_group", test_calc_tiss_avg_single_group)
    all_tests(3) = test_case("test_calc_tiss_avg_unequal_replicates", test_calc_tiss_avg_unequal_replicates)
    all_tests(4) = test_case("test_calc_tiss_avg_single_replicate", test_calc_tiss_avg_single_replicate)
    all_tests(5) = test_case("test_calc_tiss_avg_large_values", test_calc_tiss_avg_large_values)
    all_tests(6) = test_case("test_calc_tiss_avg_negative_values", test_calc_tiss_avg_negative_values)
    all_tests(7) = test_case("test_calc_tiss_avg_zero_values", test_calc_tiss_avg_zero_values)
    all_tests(8) = test_case("test_calc_tiss_avg_empty_matrix", test_calc_tiss_avg_empty_matrix)
  end function get_all_tests_tiss_avg

  !> Test tissue averaging with 3 tissues and 2 replicates each (from R test).
  subroutine test_calc_tiss_avg_three_tissues()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_c
    real(real64), dimension(3, 2) :: expected_matrix
    real(real64), dimension(6, 2) :: input_matrix
    real(real64), dimension(3, 2) :: output_matrix
    
    n_gene = 2; n_grps = 3
    input_matrix(:, 1) = [1d0, 3d0, 5d0, 2d0, 4d0, 6d0]
    input_matrix(:, 2) = [7d0, 9d0, 11d0, 8d0, 10d0, 12d0]
    group_c = [2, 2, 2]
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    ! Expected results (column-major):
    ! Tissue1: Gene1=mean(1,3)=2.0, Gene2=mean(7,9)=8.0
    ! Tissue2: Gene1=mean(5,2)=3.5, Gene2=mean(11,8)=9.5  
    ! Tissue3: Gene1=mean(4,6)=5.0, Gene2=mean(10,12)=11.0
    expected_matrix(:, 1) = [2d0, 3.5d0, 5d0]
    expected_matrix(:, 2) = [8d0, 9.5d0, 11d0]
    
    call assert_equal_array_real(output_matrix, expected_matrix, n_grps * n_gene, 1d-12, &
                            "test_calc_tiss_avg_three_tissues: averaging calculation incorrect")
  end subroutine test_calc_tiss_avg_three_tissues

  !> Test with single tissue group.
  subroutine test_calc_tiss_avg_single_group()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(1) :: group_c
    real(real64), dimension(1, 3) :: expected_matrix
    real(real64), dimension(2, 3) :: input_matrix
    real(real64), dimension(1, 3) :: output_matrix
    
    n_gene = 3; n_grps = 1
    input_matrix(:, 1) = [1d0, 4d0]
    input_matrix(:, 2) = [2d0, 5d0]
    input_matrix(:, 3) = [3d0, 6d0]
    group_c = [2]  ! 2 replicates for single tissue
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    ! Expected: average of samples 1 and 2 (columns 1 and 2)
    ! Gene1: mean(1,4)=2.5, Gene2: mean(2,5)=3.5, Gene3: mean(3,6)=4.5
    expected_matrix(1, :) = [2.5d0, 3.5d0, 4.5d0]
    
    call assert_equal_array_real(output_matrix, expected_matrix, n_grps * n_gene, 1d-12, &
                            "test_calc_tiss_avg_single_group: single tissue averaging incorrect")
  end subroutine test_calc_tiss_avg_single_group

  !> Test with unequal number of replicates per tissue.
  subroutine test_calc_tiss_avg_unequal_replicates()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_c
    real(real64), dimension(3, 2) :: expected_matrix
    real(real64), dimension(7, 2) :: input_matrix
    real(real64), dimension(3, 2) :: output_matrix
    
    n_gene = 2; n_grps = 3
    input_matrix(:, 1) = [1d0, 3d0, 5d0, 7d0, 9d0, 11d0, 13d0]
    input_matrix(:, 2) = [2d0, 4d0, 6d0, 8d0, 10d0, 12d0, 14d0]
    
    ! Tissue1: 2 replicates (cols 1,2), Tissue2: 3 replicates (cols 3,4,5), Tissue3: 2 replicates (cols 6,7)
    group_c = [2, 3, 2]
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    expected_matrix(:, 1) = [2d0, 7d0, 12d0]
    expected_matrix(:, 2) = [3d0, 8d0, 13d0]

    call assert_equal_array_real(output_matrix, expected_matrix, n_grps * n_gene, 1d-12, &
                            "test_calc_tiss_avg_unequal_replicates: averaging incorrect")
  end subroutine test_calc_tiss_avg_unequal_replicates

  !> Test with single replicate per tissue (no averaging needed).
  subroutine test_calc_tiss_avg_single_replicate()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_c
    real(real64), dimension(3, 2) :: input_matrix, output_matrix, expected_matrix
    
    n_gene = 2; n_grps = 3
    input_matrix(:, 1) = [1.0d0, 2.0d0, 3.0d0]
    input_matrix(:, 2) = [4.0d0, 5.0d0, 6.0d0]
    group_c = [1, 1, 1]  ! Single replicate per tissue
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    ! With single replicates, output should equal input
    expected_matrix = input_matrix
    
    call assert_equal_array_real(output_matrix, expected_matrix, n_gene * n_grps, 1d-12, &
                            "test_calc_tiss_avg_single_replicate: single replicate should preserve values")
  end subroutine test_calc_tiss_avg_single_replicate

  !> Test with large values for numerical stability.
  subroutine test_calc_tiss_avg_large_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_c
    real(real64), dimension(4, 2) :: input_matrix
    real(real64), dimension(2, 2) :: output_matrix, expected_matrix
    
    n_gene = 2; n_grps = 2
    input_matrix(:, 1) = [1d6, 1d9, 1d12, 1d15]
    input_matrix(:, 2) = [2d3, 2d9, 2d9, 2d15]
    group_c = [2, 2]
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    expected_matrix(:, 1) = [5.005d8, 5.005d14]
    expected_matrix(:, 2) = [1.000001d9, 1.000001d15]

    call assert_equal_array_real(output_matrix, expected_matrix, n_grps * n_gene, 1d-12, &
                            "test_calc_tiss_avg_large_values: averaging incorrect")
  end subroutine test_calc_tiss_avg_large_values

  !> Test with negative values.
  subroutine test_calc_tiss_avg_negative_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_c
    real(real64), dimension(4, 2) :: input_matrix
    real(real64), dimension(2, 2) :: output_matrix, expected_matrix
    
    n_gene = 2; n_grps = 2
    input_matrix(:, 1) = [-1d0, -3d0, -7d0, 5d0]
    input_matrix(:, 2) = [0d0, 0d0, -1d0, 1d0]
    group_c = [2, 2]
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    expected_matrix(:, 1) = [-2d0, -1d0]
    expected_matrix(:, 2) = [0d0, 0d0]

    call assert_equal_array_real(output_matrix, expected_matrix, n_grps * n_gene, 1d-12, &
                            "test_calc_tiss_avg_negative_values: averaging incorrect")
  end subroutine test_calc_tiss_avg_negative_values

  !> Test with zero values.
  subroutine test_calc_tiss_avg_zero_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_c
    real(real64), dimension(4, 2) :: input_matrix
    real(real64), dimension(2, 2) :: output_matrix, expected_matrix
    
    n_gene = 2; n_grps = 2
    input_matrix = 0d0
    group_c = [2, 2]
    expected_matrix = 0d0
    
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_OK, "calc_tiss_avg returned error")
    
    call assert_equal_array_real(output_matrix, expected_matrix, n_grps * n_gene, 1d-12, &
                            "test_calc_tiss_avg_zero_values: zero averaging incorrect")
  end subroutine test_calc_tiss_avg_zero_values

  !> Test with empty input matrix.
  subroutine test_calc_tiss_avg_empty_matrix()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(1) :: group_c
    real(real64), dimension(1) :: input_matrix, output_matrix
    n_gene = 0; n_grps = 0
    call calc_tiss_avg(n_gene, n_grps, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "calc_tiss_avg should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_calc_tiss_avg_empty_matrix

end module mod_test_calc_tiss_avg