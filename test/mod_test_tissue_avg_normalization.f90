! filepath: test/mod_test_calc_tiss_avg.f90
!> Unit test suite for calc_tiss_avg routine.
module mod_test_calc_tiss_avg
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
    
    all_tests(1) = test_case("test_calc_tiss_avg_three_tissues", test_calc_tiss_avg_three_tissues)
    all_tests(2) = test_case("test_calc_tiss_avg_generic_tissues", test_calc_tiss_avg_generic_tissues)
    all_tests(3) = test_case("test_calc_tiss_avg_preserves_dimensions", test_calc_tiss_avg_preserves_dimensions)
    all_tests(4) = test_case("test_calc_tiss_avg_single_tissue", test_calc_tiss_avg_single_tissue)
    all_tests(5) = test_case("test_calc_tiss_avg_unequal_replicates", test_calc_tiss_avg_unequal_replicates)
    all_tests(6) = test_case("test_calc_tiss_avg_single_replicate", test_calc_tiss_avg_single_replicate)
    all_tests(7) = test_case("test_calc_tiss_avg_large_values", test_calc_tiss_avg_large_values)
    all_tests(8) = test_case("test_calc_tiss_avg_negative_values", test_calc_tiss_avg_negative_values)
    all_tests(9) = test_case("test_calc_tiss_avg_zero_values", test_calc_tiss_avg_zero_values)
    all_tests(10) = test_case("test_calc_tiss_avg_mixed_values", test_calc_tiss_avg_mixed_values)
    all_tests(11) = test_case("test_calc_tiss_avg_empty_matrix", test_calc_tiss_avg_empty_matrix)
  end function get_all_tests

  !> Run all calc_tiss_avg tests.
  subroutine run_all_tests_calc_tiss_avg()
    type(test_case) :: all_tests(11)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All calc_tiss_avg tests passed successfully."
  end subroutine run_all_tests_calc_tiss_avg

  !> Run specific calc_tiss_avg tests by name.
  subroutine run_named_tests_calc_tiss_avg(test_names)
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
  end subroutine run_named_tests_calc_tiss_avg

! filepath: test/mod_test_calc_tiss_avg.f90
  !> Test tissue averaging with 3 tissues and 2 replicates each (from R test).
  subroutine test_calc_tiss_avg_three_tissues()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_s, group_c
    real(real64), dimension(12) :: input_matrix, output_matrix
    real(real64), dimension(6) :: expected_matrix
    
    n_gene = 2; n_grps = 3
    ! Input matrix (2 genes × 6 samples): column-major
    input_matrix = [1.0d0, 7.0d0, 3.0d0, 9.0d0, 5.0d0, 11.0d0, &  ! samples 1,2,3
                    2.0d0, 8.0d0, 4.0d0, 10.0d0, 6.0d0, 12.0d0]   ! samples 4,5,6
    
    ! Group starts: tissue1 starts at col 1, tissue2 at col 3, tissue3 at col 5
    group_s = [1, 3, 5]
    ! Group counts: each tissue has 2 replicates
    group_c = [2, 2, 2]

    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    ! Expected results (column-major):
    ! Tissue1: Gene1=mean(1,3)=2.0, Gene2=mean(7,9)=8.0
    ! Tissue2: Gene1=mean(5,2)=3.5, Gene2=mean(11,8)=9.5  
    ! Tissue3: Gene1=mean(4,6)=5.0, Gene2=mean(10,12)=11.0
    expected_matrix = [2.0d0, 8.0d0, 3.5d0, 9.5d0, 5.0d0, 11.0d0]
    
    call assert_equal_array_real(output_matrix(1:6), expected_matrix, 6, 1d-12, &
                            "test_calc_tiss_avg_three_tissues: averaging calculation incorrect")
    call assert_no_nan_real(output_matrix(1:6), 6, "test_calc_tiss_avg_three_tissues: NaN in result")
  end subroutine test_calc_tiss_avg_three_tissues

  !> Test tissue averaging with generic tissue names (from R test).
  subroutine test_calc_tiss_avg_generic_tissues()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_s, group_c
    real(real64), dimension(12) :: input_matrix, output_matrix
    real(real64), dimension(6) :: expected_matrix
    
    n_gene = 2; n_grps = 3
    ! Same data as previous test but representing generic tissue naming
    input_matrix = [1.0d0, 7.0d0, 3.0d0, 9.0d0, 5.0d0, 11.0d0, &
                    2.0d0, 8.0d0, 4.0d0, 10.0d0, 6.0d0, 12.0d0]
    
    group_s = [1, 3, 5]
    group_c = [2, 2, 2]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    ! Same expected results as previous test
    expected_matrix = [2.0d0, 8.0d0, 3.5d0, 9.5d0, 5.0d0, 11.0d0]
    
    call assert_equal_array_real(output_matrix(1:6), expected_matrix, 6, 1d-12, &
                            "test_calc_tiss_avg_generic_tissues: averaging calculation incorrect")
  end subroutine test_calc_tiss_avg_generic_tissues

  !> Test with single tissue group.
  subroutine test_calc_tiss_avg_single_tissue()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(1) :: group_s, group_c
    real(real64), dimension(6) :: input_matrix, output_matrix
    real(real64), dimension(3) :: expected_matrix
    
    n_gene = 3; n_grps = 1
    input_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0]
    group_s = [1]
    group_c = [2]  ! 2 replicates for single tissue
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    ! Expected: average of samples 1 and 2 (columns 1 and 2)
    ! Gene1: mean(1,4)=2.5, Gene2: mean(2,5)=3.5, Gene3: mean(3,6)=4.5
    expected_matrix = [2.5d0, 3.5d0, 4.5d0]
    
    call assert_equal_array_real(output_matrix(1:3), expected_matrix, 3, 1d-12, &
                            "test_calc_tiss_avg_single_tissue: single tissue averaging incorrect")
  end subroutine test_calc_tiss_avg_single_tissue

  !> Test that dimensions are properly calculated.
  subroutine test_calc_tiss_avg_preserves_dimensions()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_s, group_c
    real(real64), dimension(8) :: input_matrix
    real(real64), dimension(8) :: output_matrix  ! n_gene * n_grps = 4 * 2 = 8, same as input
    
    n_gene = 4; n_grps = 2
    input_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, 8.0d0]
    group_s = [1, 3]
    group_c = [2, 2]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    ! Output should have n_gene * n_grps elements
    call assert_equal_int(size(output_matrix), n_gene * n_grps, &
                     "test_calc_tiss_avg_preserves_dimensions: output size incorrect")
    call assert_no_nan_real(output_matrix, n_gene * n_grps, &
                       "test_calc_tiss_avg_preserves_dimensions: NaN in result")
  end subroutine test_calc_tiss_avg_preserves_dimensions

  !> Test with unequal number of replicates per tissue.
  subroutine test_calc_tiss_avg_unequal_replicates()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_s, group_c
    real(real64), dimension(14) :: input_matrix
    real(real64), dimension(6) :: output_matrix  ! n_gene * n_grps = 2 * 3 = 6, NOT 14
    
    n_gene = 2; n_grps = 3
    input_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0, 7.0d0, &
                    8.0d0, 9.0d0, 10.0d0, 11.0d0, 12.0d0, 13.0d0, 14.0d0]
    
    ! Tissue1: 2 replicates (cols 1,2), Tissue2: 3 replicates (cols 3,4,5), Tissue3: 2 replicates (cols 6,7)
    group_s = [1, 3, 6]
    group_c = [2, 3, 2]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    call assert_no_nan_real(output_matrix, 6, "test_calc_tiss_avg_unequal_replicates: NaN in result")
    call assert_equal_int(size(output_matrix), n_gene * n_grps, &
                     "test_calc_tiss_avg_unequal_replicates: output size incorrect")
    
    ! Check that tissue 2 average accounts for 3 replicates
    ! Gene 1, tissue 2: mean(5, 7, 11) = 7.67
    call assert_in_range_real(output_matrix(3), 7.0d0, 8.0d0, &
                         "test_calc_tiss_avg_unequal_replicates: unequal replicate averaging incorrect")
  end subroutine test_calc_tiss_avg_unequal_replicates

  !> Test with single replicate per tissue (no averaging needed).
  subroutine test_calc_tiss_avg_single_replicate()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_s, group_c
    real(real64), dimension(6) :: input_matrix, output_matrix, expected_matrix
    
    n_gene = 2; n_grps = 3
    input_matrix = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0]
    group_s = [1, 2, 3]
    group_c = [1, 1, 1]  ! Single replicate per tissue
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    ! With single replicates, output should equal input
    expected_matrix = input_matrix
    
    call assert_equal_array_real(output_matrix, expected_matrix, 6, 1d-12, &
                            "test_calc_tiss_avg_single_replicate: single replicate should preserve values")
  end subroutine test_calc_tiss_avg_single_replicate

  !> Test with large values for numerical stability.
  subroutine test_calc_tiss_avg_large_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_s, group_c
    real(real64), dimension(8) :: input_matrix
    real(real64), dimension(4) :: output_matrix  ! n_gene * n_grps = 2 * 2 = 4, NOT 8
    
    n_gene = 2; n_grps = 2
    input_matrix = [1d6, 2d6, 1d9, 2d9, 1d12, 2d12, 1d15, 2d15]
    group_s = [1, 3]
    group_c = [2, 2]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    call assert_no_nan_real(output_matrix, 4, "test_calc_tiss_avg_large_values: NaN in result")
    call assert_true(all(output_matrix > 0.0d0), "test_calc_tiss_avg_large_values: all results should be positive")
    
    ! Check that averaging works with large numbers
    ! Gene1, Tissue1: mean(1d6, 1d9) = 5.005d8
    call assert_in_range_real(output_matrix(1), 5d8, 5.1d8, &
                         "test_calc_tiss_avg_large_values: large number averaging incorrect")
  end subroutine test_calc_tiss_avg_large_values

  !> Test with negative values.
  subroutine test_calc_tiss_avg_negative_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_s, group_c
    real(real64), dimension(8) :: input_matrix
    real(real64), dimension(4) :: output_matrix  ! n_gene * n_grps = 2 * 2 = 4, NOT 8
    
    n_gene = 2; n_grps = 2
    input_matrix = [-1.0d0, -2.0d0, -3.0d0, -4.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0]
    group_s = [1, 3]
    group_c = [2, 2]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    call assert_no_nan_real(output_matrix, 4, "test_calc_tiss_avg_negative_values: NaN in result")
    
    ! Check specific negative averages
    ! Gene1, Tissue1: mean(-1, -3) = -2.0
    call assert_equal_real(output_matrix(1), -2.0d0, 1d-12, &
                      "test_calc_tiss_avg_negative_values: negative averaging incorrect")
  end subroutine test_calc_tiss_avg_negative_values

  !> Test with zero values.
  subroutine test_calc_tiss_avg_zero_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(2) :: group_s, group_c
    real(real64), dimension(8) :: input_matrix
    real(real64), dimension(4) :: output_matrix, expected_matrix
    
    n_gene = 2; n_grps = 2
    input_matrix = [0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0]
    group_s = [1, 3]
    group_c = [2, 2]
    expected_matrix = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    
    call assert_equal_array_real(output_matrix, expected_matrix, 4, 1d-12, &
                            "test_calc_tiss_avg_zero_values: zero averaging incorrect")
  end subroutine test_calc_tiss_avg_zero_values

  !> Test with mixed positive, negative, and zero values.
  subroutine test_calc_tiss_avg_mixed_values()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(3) :: group_s, group_c
    real(real64), dimension(12) :: input_matrix
    real(real64), dimension(6) :: output_matrix
    
    n_gene = 2; n_grps = 3
    input_matrix = [-5.0d0, 5.0d0, 0.0d0, 0.0d0, 10.0d0, -10.0d0, &
                    2.0d0, -2.0d0, 8.0d0, -8.0d0, 0.0d0, 0.0d0]
    group_s = [1, 3, 5]
    group_c = [2, 2, 2]
    
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 0, "calc_tiss_avg_r returned error")
    call assert_no_nan_real(output_matrix, 6, "test_calc_tiss_avg_mixed_values: NaN in result")
    
    ! Expected calculations:
    ! Tissue1: Gene1=mean(-5,0)=-2.5, Gene2=mean(5,0)=2.5
    ! Tissue2: Gene1=mean(10,2)=6.0, Gene2=mean(-10,-2)=-6.0
    ! Tissue3: Gene1=mean(0,8)=4.0, Gene2=mean(0,-8)=-4.0
    call assert_equal_real(output_matrix(1), -2.5d0, 1d-12, &
                      "test_calc_tiss_avg_mixed_values: tissue1 gene1 average incorrect")
    call assert_equal_real(output_matrix(2), 2.5d0, 1d-12, &
                      "test_calc_tiss_avg_mixed_values: tissue1 gene2 average incorrect")
  end subroutine test_calc_tiss_avg_mixed_values

  !> Test with empty input matrix.
  subroutine test_calc_tiss_avg_empty_matrix()
    integer(int32) :: n_gene, n_grps, ierr
    integer(int32), dimension(0) :: group_s, group_c
    real(real64), dimension(0) :: input_matrix, output_matrix
    n_gene = 0; n_grps = 0
    call calc_tiss_avg_r(n_gene, n_grps, group_s, group_c, input_matrix, output_matrix, ierr)
    call assert_equal_int(ierr, 202, "calc_tiss_avg_r should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_calc_tiss_avg_empty_matrix

end module mod_test_calc_tiss_avg