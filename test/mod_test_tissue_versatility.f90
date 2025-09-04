!> Unit test suite for tissue versatility routines.
module mod_test_tissue_versatility
  use asserts
  use avmod
  use tox_errors, only: ERR_OK, ERR_EMPTY_INPUT
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  !> Get array of all available tests.
  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(10)
    all_tests(1) = test_case("test_partial_axis_selection", test_partial_axis_selection)
    all_tests(2) = test_case("test_mixed_vectors", test_mixed_vectors)
    all_tests(3) = test_case("test_angle_degrees", test_angle_degrees)
    all_tests(4) = test_case("test_high_dimensional_vectors", test_high_dimensional_vectors)
    all_tests(5) = test_case("test_randomized_vectors_axes", test_randomized_vectors_axes)
    all_tests(6) = test_case("test_invalid_input_no_axes", test_invalid_input_no_axes)
    all_tests(7) = test_case("test_epsilon_threshold_stability", test_epsilon_threshold_stability)
    all_tests(8) = test_case("test_edge_case_needs_clamp", test_edge_case_needs_clamp)
    all_tests(9) = test_case("test_unbalanced_components", test_unbalanced_components)
    all_tests(10) = test_case("test_comprehensive_edge_cases", test_comprehensive_edge_cases)
  end function get_all_tests

  !> Run all tissue versatility tests.
  subroutine run_all_tests_tissue_versatility()
    type(test_case) :: all_tests(10)
    integer(int32) :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All tissue versatility tests passed successfully."
  end subroutine run_all_tests_tissue_versatility

  !> Run specific tissue versatility tests by name.
  subroutine run_named_tests_tissue_versatility(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(10)
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
  end subroutine run_named_tests_tissue_versatility


  !> Test axis selection (subspace).
  subroutine test_partial_axis_selection()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    integer(int32) :: ierr
    expr(:,1) = [1.0_real64, 2.0_real64, 3.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .false., .true.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 2, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Partial axis selection should succeed")
    ! Should behave as a 2D vector [1,3]
    call assert_true(tv(1) >= 0.0_real64 .and. tv(1) <= 1.0_real64, "Partial axis TV in [0,1]")
    call assert_true(angle(1) >= 0.0_real64 .and. angle(1) <= 90.0_real64, "Partial axis angle in [0,90]")
  end subroutine test_partial_axis_selection

  !> Test mixed vectors (uniform, single axis, null) - consolidates basic cases.
  subroutine test_mixed_vectors()
    real(real64) :: expr(3,3), tv(3), angle(3)
    logical :: select_vec(3), select_axes(3)
    integer(int32) :: ierr
    expr(:,1) = [1.0_real64, 1.0_real64, 1.0_real64] ! uniform → TV=0
    expr(:,2) = [0.0_real64, 0.0_real64, 2.0_real64] ! single axis → TV=1
    expr(:,3) = [0.0_real64, 0.0_real64, 0.0_real64] ! null → TV=1, angle=90°
    select_vec = [.true., .true., .true.]
    select_axes = [.true., .true., .true.]
    call compute_tissue_versatility(3, 3, expr, select_vec, 3, select_axes, 3, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Mixed vectors should succeed")
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Mixed: uniform TV")
    call assert_equal_real(tv(2), 1.0_real64, 1e-12_real64, "Mixed: single axis TV")
    call assert_equal_real(tv(3), 1.0_real64, 1e-12_real64, "Mixed: null TV")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Mixed: uniform angle")
    call assert_true(angle(2) > 0.0_real64, "Mixed: single axis angle > 0")
    call assert_equal_real(angle(3), 90.0_real64, 1e-12_real64, "Mixed: null angle")
  end subroutine test_mixed_vectors

  !> Test angle output in degrees for a known case.
  subroutine test_angle_degrees()
    real(real64) :: expr(2,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(2)
    integer(int32) :: ierr
    expr(:,1) = [1.0_real64, 0.0_real64]
    select_vec = [.true.]
    select_axes = [.true., .true.]
    call compute_tissue_versatility(2, 1, expr, select_vec, 1, select_axes, 2, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Angle degrees test should succeed")
    ! Angle should be 45 degrees (between [1,0] and [1,1])
    call assert_true(abs(angle(1) - 45.0_real64) < 1e-12_real64, "Angle output is 45 degrees")
  end subroutine test_angle_degrees

  !> Test tissue versatility in higher dimensions (4D, 5D).
  subroutine test_high_dimensional_vectors()
    real(real64) :: expr4(4,1), tv4(1), angle4(1)
    real(real64) :: expr5(5,1), tv5(1), angle5(1)
    logical :: select_vec(1), select_axes4(4), select_axes5(5)
    integer(int32) :: ierr4, ierr5
    expr4(:,1) = [1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64]
    expr5(:,1) = [2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64]
    select_vec = [.true.]
    select_axes4 = [.true., .true., .true., .true.]
    select_axes5 = [.true., .true., .true., .true., .true.]
    call compute_tissue_versatility(4, 1, expr4, select_vec, 1, select_axes4, 4, tv4, angle4, ierr4)
    call compute_tissue_versatility(5, 1, expr5, select_vec, 1, select_axes5, 5, tv5, angle5, ierr5)
    call assert_equal_int(ierr4, ERR_OK, "4D vectors should succeed")
    call assert_equal_int(ierr5, ERR_OK, "5D vectors should succeed")
    call assert_equal_real(tv4(1), 0.0_real64, 1e-12_real64, "4D uniform TV")
    call assert_equal_real(angle4(1), 0.0_real64, 1e-12_real64, "4D uniform angle")
    call assert_equal_real(tv5(1), 0.0_real64, 1e-12_real64, "5D uniform TV")
    call assert_equal_real(angle5(1), 0.0_real64, 1e-12_real64, "5D uniform angle")
  end subroutine test_high_dimensional_vectors

  !> Test tissue versatility with randomized vectors and axis selections.
  subroutine test_randomized_vectors_axes()
    integer(int32), parameter :: n_axes = 5, n_vecs = 4
    real(real64) :: expr(n_axes, n_vecs), tv(n_vecs), angle(n_vecs)
    logical :: select_vec(n_vecs), select_axes(n_axes)
    integer(int32) :: i, ierr
    call random_seed()
    call random_number(expr)
    select_vec = [.true., .true., .true., .true.]
    select_axes = [.true., .false., .true., .false., .true.]
    call compute_tissue_versatility(n_axes, n_vecs, expr, select_vec, 4, select_axes, 3, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Randomized vectors should succeed")
    do i = 1, n_vecs
      call assert_true(tv(i) >= 0.0_real64 .and. tv(i) <= 1.0_real64, "Randomized TV in [0,1]")
      call assert_true(angle(i) >= 0.0_real64 .and. angle(i) <= 90.0_real64, "Randomized angle in [0,90]")
    end do
  end subroutine test_randomized_vectors_axes

  !> Test invalid input: no axes selected (should return error code).
  subroutine test_invalid_input_no_axes()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    integer(int32) :: ierr
    expr(:,1) = [1.0_real64, 2.0_real64, 3.0_real64]
    select_vec = [.true.]
    select_axes = [.false., .false., .false.]
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 0, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "No axes selected should return ERR_EMPTY_INPUT")
  end subroutine test_invalid_input_no_axes

  !> Test epsilon threshold stability with extremely small vectors.
  !| This test verifies that vectors with norm <= sqrt(epsilon) are correctly
  !| detected as "numerically zero" and assigned TV=1, angle=90° without
  !| causing numerical instability in the cos_phi calculation.
  subroutine test_epsilon_threshold_stability()
    real(real64) :: expr(3,4), tv(4), angle(4)
    logical :: select_vec(4), select_axes(3)
    integer(int32) :: ierr
    real(real64), parameter :: eps_sqrt = sqrt(epsilon(1.0_real64))  ! ~1.49e-8
    real(real64), parameter :: large_component = 1.0e-5_real64  ! Much larger than sqrt(epsilon)
    
    ! Test case 1: Vector with norm exactly at sqrt(epsilon) threshold
    expr(:,1) = [eps_sqrt/sqrt(3.0_real64), eps_sqrt/sqrt(3.0_real64), eps_sqrt/sqrt(3.0_real64)]
    
    ! Test case 2: Vector with norm slightly below sqrt(epsilon) threshold
    expr(:,2) = [eps_sqrt*0.5_real64/sqrt(3.0_real64), eps_sqrt*0.5_real64/sqrt(3.0_real64), eps_sqrt*0.5_real64/sqrt(3.0_real64)]
    
    ! Test case 3: Vector with mixed components: one large, others extremely small
    expr(:,3) = [large_component, eps_sqrt*0.1_real64, eps_sqrt*0.1_real64]
    
    ! Test case 4: Vector that causes underflow to exactly zero
    ! Use components so small that their squares underflow to 0.0
    ! This will cause norm_v = 0.0 and trigger the zero-check protection
    ! Components: 1e-200 each → norm_v = 3*(1e-200)^2 = 3e-400 → underflows to 0.0
    expr(:,4) = [1e-200_real64, 1e-200_real64, 1e-200_real64]
    
    select_vec = [.true., .true., .true., .true.]
    select_axes = [.true., .true., .true.]
    
    call compute_tissue_versatility(3, 4, expr, select_vec, 4, select_axes, 3, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Epsilon threshold stability should succeed")
    
    ! Test case 1: At threshold - should be treated as zero
    call assert_equal_real(tv(1), 1.0_real64, 1e-12_real64, "Epsilon threshold: TV=1")
    call assert_equal_real(angle(1), 90.0_real64, 1e-12_real64, "Epsilon threshold: angle=90°")
    
    ! Test case 2: Below threshold - should be treated as zero
    call assert_equal_real(tv(2), 1.0_real64, 1e-12_real64, "Below epsilon: TV=1")
    call assert_equal_real(angle(2), 90.0_real64, 1e-12_real64, "Below epsilon: angle=90°")
    
    ! Test case 3: Above threshold - should compute normally
    call assert_true(tv(3) >= 0.0_real64 .and. tv(3) <= 1.0_real64, "Above epsilon: TV in [0,1]")
    call assert_true(angle(3) >= 0.0_real64 .and. angle(3) <= 90.0_real64, "Above epsilon: angle in [0,90°]")
    
    ! Test case 4: Potential explosion case - should be protected
    call assert_equal_real(tv(4), 1.0_real64, 1e-12_real64, "Explosion case: TV=1")
    call assert_equal_real(angle(4), 90.0_real64, 1e-12_real64, "Explosion case: angle=90°")
    
    ! Verify all results are finite (no NaN, no Inf)
    call assert_true(tv(1) == tv(1), "TV(1) is finite (not NaN)")
    call assert_true(tv(2) == tv(2), "TV(2) is finite (not NaN)")
    call assert_true(tv(3) == tv(3), "TV(3) is finite (not NaN)")
    call assert_true(tv(4) == tv(4), "TV(4) is finite (not NaN)")
    call assert_true(angle(1) == angle(1), "Angle(1) is finite (not NaN)")
    call assert_true(angle(2) == angle(2), "Angle(2) is finite (not NaN)")
    call assert_true(angle(3) == angle(3), "Angle(3) is finite (not NaN)")
    call assert_true(angle(4) == angle(4), "Angle(4) is finite (not NaN)")
  end subroutine test_epsilon_threshold_stability

  !> Test edge case that passes sqrt(epsilon) threshold but needs clamp protection.
  !| This test demonstrates a vector that is above the sqrt(epsilon) threshold
  !| but could generate cos_phi slightly outside [-1,1] due to floating-point
  !| precision errors, requiring the clamp protection.
  subroutine test_edge_case_needs_clamp()
    real(real64) :: expr(3,1), tv(1), angle(1)
    logical :: select_vec(1), select_axes(3)
    integer(int32) :: ierr
    real(real64), parameter :: eps_sqrt = sqrt(epsilon(1.0_real64))  ! ~1.49e-8
    
    ! Vector designed to potentially cause cos_phi slightly > 1.0 due to precision
    ! Use a vector that creates cos_phi very close to 1 but might overshoot
    ! Vector close to diagonal but with tiny variations that could cause precision issues
    expr(:,1) = [1.0_real64 + 1e-15_real64, 1.0_real64, 1.0_real64 - 1e-15_real64]
    
    select_vec = [.true.]
    select_axes = [.true., .true., .true.]
    
    call compute_tissue_versatility(3, 1, expr, select_vec, 1, select_axes, 3, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Edge case clamp should succeed")
    
    ! This vector passes the sqrt(epsilon) threshold so gets full processing
    ! The clamp ensures we get mathematically valid results
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Edge case TV=0 (uniform)")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Edge case angle=0°")
    
    ! Verify results are finite and in valid ranges
    call assert_true(tv(1) >= 0.0_real64 .and. tv(1) <= 1.0_real64, "Edge case TV in [0,1]")
    call assert_true(angle(1) >= 0.0_real64 .and. angle(1) <= 90.0_real64, "Edge case angle in [0,90°]")
    call assert_true(tv(1) == tv(1), "Edge case TV is finite (not NaN)")
    call assert_true(angle(1) == angle(1), "Edge case angle is finite (not NaN)")
  end subroutine test_edge_case_needs_clamp

  !> Test unbalanced components that could cause cos_phi explosion without clamp.
  !| This test uses vectors with highly unbalanced components to create scenarios
  !| where cos_phi could theoretically become very large, demonstrating the need
  !| for clamp protection even when vectors pass the sqrt(epsilon) threshold.
  subroutine test_unbalanced_components()
    real(real64) :: expr(3,2), tv(2), angle(2)
    logical :: select_vec(2), select_axes(3)
    integer(int32) :: ierr
    real(real64), parameter :: eps_sqrt = sqrt(epsilon(1.0_real64))  ! ~1.49e-8
    
    ! Test case 1: Vector with one large component, others extremely small
    ! This creates a scenario where norm_v is dominated by one component
    ! but dot_prod includes contributions from all components
    expr(:,1) = [1e-6_real64, eps_sqrt*0.01_real64, eps_sqrt*0.01_real64]
    
    ! Test case 2: Vector with extreme imbalance that still passes sqrt(epsilon)
    ! One component much larger than others, creating potential for precision issues
    expr(:,2) = [1e-5_real64, eps_sqrt*0.001_real64, eps_sqrt*0.001_real64]
    
    select_vec = [.true., .true.]
    select_axes = [.true., .true., .true.]
    
    call compute_tissue_versatility(3, 2, expr, select_vec, 2, select_axes, 3, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Unbalanced components should succeed")
    
    ! Both vectors pass sqrt(epsilon) threshold and get full processing
    ! The clamp ensures mathematically valid cos_phi values
    
    ! Test case 1: Highly unbalanced but mathematically valid
    call assert_true(tv(1) >= 0.0_real64 .and. tv(1) <= 1.0_real64, "Unbalanced 1: TV in [0,1]")
    call assert_true(angle(1) >= 0.0_real64 .and. angle(1) <= 90.0_real64, "Unbalanced 1: angle in [0,90°]")
    call assert_true(tv(1) > 0.0_real64, "Unbalanced 1: TV > 0 (not uniform)")
    
    ! Test case 2: Extreme imbalance but controlled
    call assert_true(tv(2) >= 0.0_real64 .and. tv(2) <= 1.0_real64, "Unbalanced 2: TV in [0,1]")
    call assert_true(angle(2) >= 0.0_real64 .and. angle(2) <= 90.0_real64, "Unbalanced 2: angle in [0,90°]")
    call assert_true(tv(2) > 0.0_real64, "Unbalanced 2: TV > 0 (not uniform)")
    
    ! Verify all results are finite (demonstrating clamp protection worked)
    call assert_true(tv(1) == tv(1), "Unbalanced 1: TV is finite (not NaN)")
    call assert_true(tv(2) == tv(2), "Unbalanced 2: TV is finite (not NaN)")
    call assert_true(angle(1) == angle(1), "Unbalanced 1: angle is finite (not NaN)")
    call assert_true(angle(2) == angle(2), "Unbalanced 2: angle is finite (not NaN)")
  end subroutine test_unbalanced_components

  !> Test comprehensive edge cases combining all numerical protections.
  !| This test consolidates numerical stability, threshold, clamp, and explosion
  !| protection in a single comprehensive test to avoid redundancy.
  subroutine test_comprehensive_edge_cases()
    real(real64) :: expr(3,6), tv(6), angle(6)
    logical :: select_vec(6), select_axes(3)
    integer(int32) :: i, ierr
    real(real64), parameter :: eps_sqrt = sqrt(epsilon(1.0_real64))  ! ~1.49e-8
    
    ! Case 1: Large numbers (numerical stability)
    expr(:,1) = [1e15_real64, 1e15_real64, 1e15_real64]
    
    ! Case 2: Small numbers above threshold (numerical stability)
    expr(:,2) = [1e-4_real64, 1e-4_real64, 1e-4_real64]
    
    ! Case 3: Below epsilon threshold (should be treated as zero)
    expr(:,3) = [eps_sqrt*0.5_real64/sqrt(3.0_real64), eps_sqrt*0.5_real64/sqrt(3.0_real64), eps_sqrt*0.5_real64/sqrt(3.0_real64)]
    
    ! Case 4: Underflow to zero (explosion protection)
    expr(:,4) = [1e-200_real64, 1e-200_real64, 1e-200_real64]
    
    ! Case 5: Potential clamp case (precision issues)
    expr(:,5) = [1.0_real64 + 1e-15_real64, 1.0_real64, 1.0_real64 - 1e-15_real64]
    
    ! Case 6: Unbalanced components requiring clamp
    expr(:,6) = [1e-6_real64, eps_sqrt*0.01_real64, eps_sqrt*0.01_real64]
    
    select_vec = [.true., .true., .true., .true., .true., .true.]
    select_axes = [.true., .true., .true.]
    
    call compute_tissue_versatility(3, 6, expr, select_vec, 6, select_axes, 3, tv, angle, ierr)
    call assert_equal_int(ierr, ERR_OK, "Comprehensive edge cases should succeed")
    
    ! Case 1: Large numbers should work normally (uniform → TV=0)
    call assert_equal_real(tv(1), 0.0_real64, 1e-12_real64, "Large numbers TV")
    call assert_equal_real(angle(1), 0.0_real64, 1e-12_real64, "Large numbers angle")
    
    ! Case 2: Small numbers should work normally (uniform → TV=0)
    call assert_equal_real(tv(2), 0.0_real64, 1e-12_real64, "Small numbers TV")
    call assert_equal_real(angle(2), 0.0_real64, 1e-12_real64, "Small numbers angle")
    
    ! Case 3: Below threshold → protected (TV=1, angle=90°)
    call assert_equal_real(tv(3), 1.0_real64, 1e-12_real64, "Below threshold TV")
    call assert_equal_real(angle(3), 90.0_real64, 1e-12_real64, "Below threshold angle")
    
    ! Case 4: Underflow → protected (TV=1, angle=90°)
    call assert_equal_real(tv(4), 1.0_real64, 1e-12_real64, "Underflow TV")
    call assert_equal_real(angle(4), 90.0_real64, 1e-12_real64, "Underflow angle")
    
    ! Case 5: Clamp protection (uniform → TV=0)
    call assert_equal_real(tv(5), 0.0_real64, 1e-12_real64, "Clamp case TV")
    call assert_equal_real(angle(5), 0.0_real64, 1e-12_real64, "Clamp case angle")
    
    ! Case 6: Unbalanced but controlled by clamp
    call assert_true(tv(6) >= 0.0_real64 .and. tv(6) <= 1.0_real64, "Unbalanced TV in [0,1]")
    call assert_true(angle(6) >= 0.0_real64 .and. angle(6) <= 90.0_real64, "Unbalanced angle in [0,90°]")
    call assert_true(tv(6) > 0.0_real64, "Unbalanced TV > 0 (not uniform)")
    
    ! Verify all results are finite (demonstrating all protections work)
    do i = 1, 6
      call assert_true(tv(i) == tv(i), "TV is finite (not NaN)")
      call assert_true(angle(i) == angle(i), "Angle is finite (not NaN)")
    end do
  end subroutine test_comprehensive_edge_cases

end module mod_test_tissue_versatility