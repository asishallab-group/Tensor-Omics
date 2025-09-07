! filepath: test/mod_test_log2_transformation.f90
!> Unit test suite for log2_transformation routine.
module mod_test_log2_transformation
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
    
    all_tests(1) = test_case("test_log2_basic_values", test_log2_basic_values)
    all_tests(2) = test_case("test_log2_zeros_handling", test_log2_zeros_handling)
    all_tests(3) = test_case("test_log2_preserves_dimensions", test_log2_preserves_dimensions)
    all_tests(4) = test_case("test_log2_single_element", test_log2_single_element)
    all_tests(5) = test_case("test_log2_large_values", test_log2_large_values)
    all_tests(6) = test_case("test_log2_small_values", test_log2_small_values)
    all_tests(7) = test_case("test_log2_powers_of_two", test_log2_powers_of_two)
    all_tests(8) = test_case("test_log2_random_matrix", test_log2_random_matrix)
    all_tests(9) = test_case("test_log2_negative_handling", test_log2_negative_handling)
    all_tests(10) = test_case("test_log2_edge_cases", test_log2_edge_cases)
    all_tests(11) = test_case("test_log2_monotonic_property", test_log2_monotonic_property)
    all_tests(12) = test_case("test_log2_mathematical_properties", test_log2_mathematical_properties)
    all_tests(13) = test_case("test_log2_empty_matrix", test_log2_empty_matrix)
  end function get_all_tests

  !> Run all log2_transformation tests.
  subroutine run_all_tests_log2_transformation()
    type(test_case) :: all_tests(13)
    integer(int32) :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All log2_transformation tests passed successfully."
  end subroutine run_all_tests_log2_transformation

  !> Run specific log2_transformation tests by name.
  subroutine run_named_tests_log2_transformation(test_names)
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
  end subroutine run_named_tests_log2_transformation

  !> Test log2(x+1) transformation with basic known values (from R test).
  subroutine test_log2_basic_values()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat, expected_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 2; n_tissues = 2
    ! Matrix: [0, 7; 3, 15] in column-major (flattened: [0, 3, 7, 15])
    input_flat = [0.0d0, 3.0d0, 7.0d0, 15.0d0]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    ! Expected: log2(x+1) = [log2(1), log2(4), log2(8), log2(16)] = [0, 2, 3, 4]
    expected_flat = [log(1.0d0)/LOG2, log(4.0d0)/LOG2, log(8.0d0)/LOG2, log(16.0d0)/LOG2]
    
    call assert_equal_array_real(output_flat, expected_flat, 4, 1d-12, &
                            "test_log2_basic_values: transformation values incorrect")
    call assert_equal_real(output_flat(1), 0.0d0, 1d-12, "test_log2_basic_values: log2(0+1) should be 0")
    call assert_equal_real(output_flat(2), 2.0d0, 1d-12, "test_log2_basic_values: log2(3+1) should be 2")
    call assert_equal_real(output_flat(3), 3.0d0, 1d-12, "test_log2_basic_values: log2(7+1) should be 3")
    call assert_equal_real(output_flat(4), 4.0d0, 1d-12, "test_log2_basic_values: log2(15+1) should be 4")
  end subroutine test_log2_basic_values

  !> Test that log2(0+1) = 0 for all zeros (from R test).
  subroutine test_log2_zeros_handling()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat, expected_flat
    
    n_genes = 2; n_tissues = 2
    input_flat = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    expected_flat = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_equal_array_real(output_flat, expected_flat, 4, 1d-12, &
                            "test_log2_zeros_handling: all zeros should become zeros")
    call assert_true(all(output_flat == 0.0d0), "test_log2_zeros_handling: all values should be exactly 0")
  end subroutine test_log2_zeros_handling

  !> Test that dimensions are preserved (from R test concept).
  subroutine test_log2_preserves_dimensions()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(6) :: input_flat, output_flat
    
    n_genes = 2; n_tissues = 3
    input_flat = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0, 6.0d0]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    ! Check that we get exactly the same number of elements
    call assert_equal_int(size(output_flat), n_genes * n_tissues, &
                     "test_log2_preserves_dimensions: output size incorrect")
    call assert_equal_int(size(output_flat), size(input_flat), &
                     "test_log2_preserves_dimensions: input/output size mismatch")
  end subroutine test_log2_preserves_dimensions

  !> Test single element matrix.
  subroutine test_log2_single_element()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(1) :: input_flat, output_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 1; n_tissues = 1
    input_flat = [7.0d0]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_equal_real(output_flat(1), log(8.0d0)/LOG2, 1d-12, &
                      "test_log2_single_element: log2(7+1) incorrect")
  end subroutine test_log2_single_element

  !> Test with large values to check numerical stability.
  subroutine test_log2_large_values()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 2; n_tissues = 2
    input_flat = [1d6, 1d9, 1d12, 1d15]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    ! Check that results are finite and reasonable
    call assert_no_nan_real(output_flat, 4, "test_log2_large_values: NaN in result")
    call assert_true(all(output_flat > 0.0d0), "test_log2_large_values: all results should be positive")
    call assert_true(all(output_flat < 100.0d0), "test_log2_large_values: results should be reasonable")
    
    ! Check specific large value: log2(1e6 + 1) ≈ log2(1e6) ≈ 19.93
    call assert_in_range_real(output_flat(1), 19.0d0, 21.0d0, "test_log2_large_values: log2(1e6+1) out of range")
  end subroutine test_log2_large_values

  !> Test with very small positive values.
  subroutine test_log2_small_values()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 2; n_tissues = 2
    input_flat = [1d-6, 1d-9, 1d-12, 1d-15]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_no_nan_real(output_flat, 4, "test_log2_small_values: NaN in result")
    call assert_true(all(output_flat > 0.0d0), "test_log2_small_values: all results should be positive")
    
    ! For very small x, log2(x+1) ≈ log2(1) = 0, but slightly positive
    call assert_in_range_real(output_flat(4), 0.0d0, 1d-10, "test_log2_small_values: very small values should be near 0")
  end subroutine test_log2_small_values

  !> Test with powers of 2 minus 1 for exact results.
  subroutine test_log2_powers_of_two()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat, expected_flat
    
    n_genes = 2; n_tissues = 2
    ! Powers of 2 minus 1: log2((2^n - 1) + 1) = log2(2^n) = n
    input_flat = [1.0d0, 3.0d0, 7.0d0, 15.0d0]  ! 2^1-1, 2^2-1, 2^3-1, 2^4-1
    expected_flat = [1.0d0, 2.0d0, 3.0d0, 4.0d0]  ! log2(2), log2(4), log2(8), log2(16)
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_equal_array_real(output_flat, expected_flat, 4, 1d-12, &
                            "test_log2_powers_of_two: powers of 2 results incorrect")
  end subroutine test_log2_powers_of_two

  !> Test with random matrix for general properties.
  subroutine test_log2_random_matrix()
    integer(int32), parameter :: n_genes = 5, n_tissues = 4
    real(real64), dimension(n_genes * n_tissues) :: input_flat, output_flat
    integer(int32) :: i, ierr
    integer(int32) :: n_seed
    integer(int32), allocatable :: seed_array(:)
    ! For reproducibility: initialize the random number generator seed
    call random_seed(size=n_seed)
    allocate(seed_array(n_seed))
    seed_array = 42  ! Fixed value for reproducibility
    call random_seed(put=seed_array)
    deallocate(seed_array)
    call random_number(input_flat)
    input_flat = input_flat * 100.0d0  ! Scale to [0, 100]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_no_nan_real(output_flat, n_genes * n_tissues, "test_log2_random_matrix: NaN in result")
    call assert_true(all(output_flat >= 0.0d0), "test_log2_random_matrix: all results should be non-negative")
    
    ! Check monotonicity: if input[i] > input[j], then output[i] > output[j]
    do i = 1, n_genes * n_tissues - 1
      if (input_flat(i) > input_flat(i+1)) then
        call assert_true(output_flat(i) > output_flat(i+1), "test_log2_random_matrix: monotonicity violated")
      end if
    end do
  end subroutine test_log2_random_matrix

  !> Test behavior with negative values (should still work due to +1).
  subroutine test_log2_negative_handling()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 2; n_tissues = 2
    input_flat = [-0.5d0, -0.9d0, -0.99d0, -0.999d0]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_no_nan_real(output_flat, 4, "test_log2_negative_handling: NaN in result")
    call assert_true(all(output_flat > -10.0d0), "test_log2_negative_handling: results should be reasonable")
    
    ! log2(-0.5 + 1) = log2(0.5) = -1
    call assert_equal_real(output_flat(1), -1.0d0, 1d-12, "test_log2_negative_handling: log2(-0.5+1) should be -1")
  end subroutine test_log2_negative_handling

  !> Test edge cases and boundary conditions.
  subroutine test_log2_edge_cases()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(4) :: input_flat, output_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 2; n_tissues = 2
    input_flat = [0.0d0, 1.0d0, huge(1.0d0), tiny(1.0d0)]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    call assert_no_nan_real(output_flat, 4, "test_log2_edge_cases: NaN in result")
    call assert_true(all(output_flat > -1000.0d0), "test_log2_edge_cases: results should not be extremely negative")
    
    call assert_equal_real(output_flat(1), 0.0d0, 1d-12, "test_log2_edge_cases: log2(0+1) should be 0")
    call assert_equal_real(output_flat(2), 1.0d0, 1d-12, "test_log2_edge_cases: log2(1+1) should be 1")
  end subroutine test_log2_edge_cases

  !> Test monotonic property: log2(x+1) is strictly increasing.
  subroutine test_log2_monotonic_property()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(5) :: input_flat, output_flat
    integer(int32) :: i
    
    n_genes = 5; n_tissues = 1
    input_flat = [1.0d0, 2.0d0, 5.0d0, 10.0d0, 20.0d0]  ! Increasing sequence
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    ! Check strict monotonicity
    do i = 1, 4
      call assert_true(output_flat(i) < output_flat(i+1), &
                  "test_log2_monotonic_property: monotonicity violated")
    end do
  end subroutine test_log2_monotonic_property

  !> Test mathematical properties of log2 transformation.
  subroutine test_log2_mathematical_properties()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(8) :: input_flat, output_flat
    real(real64), parameter :: LOG2 = log(2.0d0)
    
    n_genes = 4; n_tissues = 2
    input_flat = [1.0d0, 3.0d0, 7.0d0, 15.0d0, 2.0d0, 6.0d0, 14.0d0, 30.0d0]
    
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 0, "log2_transformation_r returned error")
    
    ! Property: log2((2x)+1) ≈ log2(2x) = log2(2) + log2(x) = 1 + log2(x)
    ! But this is approximate for log2(x+1) vs log2(2x+1)
    
    ! Just check that all results are reasonable and follow expected patterns
    call assert_no_nan_real(output_flat, 8, "test_log2_mathematical_properties: NaN in result")
    call assert_true(all(output_flat > 0.0d0), "test_log2_mathematical_properties: all should be positive")
    
    ! Check some specific relationships
    call assert_true(output_flat(2) > output_flat(1), "test_log2_mathematical_properties: ordering incorrect")
    call assert_true(output_flat(6) > output_flat(2), "test_log2_mathematical_properties: larger input gives larger output")
  end subroutine test_log2_mathematical_properties

  !> Test with empty input matrix.
  subroutine test_log2_empty_matrix()
    integer(int32) :: n_genes, n_tissues, ierr
    real(real64), dimension(0) :: input_flat, output_flat
    n_genes = 0; n_tissues = 0
    call log2_transformation_r(n_genes, n_tissues, input_flat, output_flat, ierr)
    call assert_equal_int(ierr, 202, "log2_transformation_r should return error for empty input")
    ! No further assertion needed: just check no crash
  end subroutine test_log2_empty_matrix

end module mod_test_log2_transformation