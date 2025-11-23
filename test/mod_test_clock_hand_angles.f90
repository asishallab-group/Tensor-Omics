module mod_test_clock_hand_angles
  use asserts
  use relative_axis_plane_tools
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  public

  ! Mathematical constants
  real(real64), parameter :: PI = 3.141592653589793_real64
  real(real64), parameter :: TOL = 1e-12_real64

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
    type(test_case) :: all_tests(25)
    
    all_tests(1) = test_case("test_identical_vectors_2d", test_identical_vectors_2d)
    all_tests(2) = test_case("test_opposite_vectors_2d", test_opposite_vectors_2d)
    all_tests(3) = test_case("test_perpendicular_vectors_2d", test_perpendicular_vectors_2d)
    all_tests(4) = test_case("test_45_degree_rotation_2d", test_45_degree_rotation_2d)
    all_tests(5) = test_case("test_clockwise_vs_counterclockwise_2d", test_clockwise_vs_counterclockwise_2d)
    all_tests(6) = test_case("test_identical_vectors_3d", test_identical_vectors_3d)
    all_tests(7) = test_case("test_perpendicular_vectors_3d", test_perpendicular_vectors_3d)
    all_tests(8) = test_case("test_arbitrary_3d_rotation", test_arbitrary_3d_rotation)
    all_tests(9) = test_case("test_high_dimensional_basic", test_high_dimensional_basic)
    all_tests(10) = test_case("test_high_dimensional_selected_axes", test_high_dimensional_selected_axes)
    all_tests(11) = test_case("test_zero_vectors", test_zero_vectors)
    all_tests(12) = test_case("test_denormalized_vectors", test_denormalized_vectors)
    all_tests(13) = test_case("test_tiny_vectors_precision", test_tiny_vectors_precision)
    all_tests(14) = test_case("test_huge_vectors_precision", test_huge_vectors_precision)
    all_tests(15) = test_case("test_nearly_identical_vectors", test_nearly_identical_vectors)
    all_tests(16) = test_case("test_nearly_opposite_vectors", test_nearly_opposite_vectors)
    all_tests(17) = test_case("test_mixed_positive_negative", test_mixed_positive_negative)
    all_tests(18) = test_case("test_single_pair_shift_vectors", test_single_pair_shift_vectors)
    all_tests(19) = test_case("test_multiple_pairs_shift_vectors", test_multiple_pairs_shift_vectors)
    all_tests(20) = test_case("test_shift_vectors_with_selection_mask", test_shift_vectors_with_selection_mask)
    all_tests(21) = test_case("test_invalid_selected_axes", test_invalid_selected_axes)
    all_tests(22) = test_case("test_out_of_bounds_selected_axes", test_out_of_bounds_selected_axes)
    all_tests(23) = test_case("test_performance_large_scale", test_performance_large_scale)
    all_tests(24) = test_case("test_consistency_between_functions", test_consistency_between_functions)
    all_tests(25) = test_case("test_mathematical_properties", test_mathematical_properties)
  end function get_all_tests

  !> Run all clock hand angle tests.
  subroutine run_all_tests_clock_hand_angles()
    type(test_case) :: all_tests(25)
    integer :: i
    
    all_tests = get_all_tests()
    
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All clock hand angle tests passed successfully."
  end subroutine run_all_tests_clock_hand_angles

  !> Run specific clock hand angle tests by name.
  subroutine run_named_tests_clock_hand_angles(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(25)
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
  end subroutine run_named_tests_clock_hand_angles

  ! ==================== 2D TESTS ====================

  !> Test identical vectors in 2D (should give 0 angle).
  subroutine test_identical_vectors_2d()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64]
    v2 = [1.0_real64, 0.0_real64]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: identical 2D vectors")
    call assert_equal_real(signed_angle, 0.0_real64, TOL, "Identical 2D vectors")
  end subroutine test_identical_vectors_2d

  !> Test opposite vectors in 2D (should give ±π).
  subroutine test_opposite_vectors_2d()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64]
    v2 = [-1.0_real64, 0.0_real64]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: opposite 2D vectors")
    call assert_equal_real(abs(signed_angle), PI, TOL, "Opposite 2D vectors magnitude")
  end subroutine test_opposite_vectors_2d

  !> Test perpendicular vectors in 2D (should give ±π/2).
  subroutine test_perpendicular_vectors_2d()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64]
    v2 = [0.0_real64, 1.0_real64]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: perpendicular 2D vectors")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, TOL, "Perpendicular 2D vectors magnitude")
    call assert_true(signed_angle > 0.0_real64, "Counterclockwise rotation should be positive")
  end subroutine test_perpendicular_vectors_2d

  !> Test 45-degree rotation in 2D.
  subroutine test_45_degree_rotation_2d()
    real(real64) :: v1(2), v2(2), signed_angle, expected
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64]
    v2 = [sqrt(2.0_real64)/2.0_real64, sqrt(2.0_real64)/2.0_real64]
    expected = PI/4.0_real64
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: 45-degree rotation")
    call assert_equal_real(signed_angle, expected, TOL, "45-degree counterclockwise rotation")
  end subroutine test_45_degree_rotation_2d

  !> Test clockwise vs counterclockwise rotations in 2D.
  subroutine test_clockwise_vs_counterclockwise_2d()
    real(real64) :: v1(2), v2_ccw(2), v2_cw(2), angle_ccw, angle_cw
    integer :: selected_axes(3), ierr_ccw, ierr_cw
    v1 = [1.0_real64, 0.0_real64]
    v2_ccw = [0.0_real64, 1.0_real64]
    v2_cw = [0.0_real64, -1.0_real64]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2_ccw, 2, angle_ccw, selected_axes, ierr_ccw)
    call clock_hand_angle_between_vectors(v1, v2_cw, 2, angle_cw, selected_axes, ierr_cw)
    call assert_equal_int(ierr_ccw, 0, "ierr should be 0 for valid input: ccw")
    call assert_equal_int(ierr_cw, 0, "ierr should be 0 for valid input: cw")
    call assert_true(angle_ccw > 0.0_real64, "Counterclockwise should be positive")
    call assert_true(angle_cw < 0.0_real64, "Clockwise should be negative")
    call assert_equal_real(abs(angle_ccw), abs(angle_cw), TOL, "Magnitudes should be equal")
  end subroutine test_clockwise_vs_counterclockwise_2d

  ! ==================== 3D TESTS ====================

  !> Test identical vectors in 3D.
  subroutine test_identical_vectors_3d()
    real(real64) :: v1(3), v2(3), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 1.0_real64, 1.0_real64]
    v2 = [1.0_real64, 1.0_real64, 1.0_real64]
    selected_axes = [1, 2, 3]
    call clock_hand_angle_between_vectors(v1, v2, 3, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: identical 3D vectors")
    call assert_equal_real(signed_angle, 0.0_real64, TOL, "Identical 3D vectors")
  end subroutine test_identical_vectors_3d

  !> Test perpendicular vectors in 3D.
  subroutine test_perpendicular_vectors_3d()
    real(real64) :: v1(3), v2(3), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64, 0.0_real64]
    v2 = [0.0_real64, 1.0_real64, 0.0_real64]
    selected_axes = [1, 2, 3]
    call clock_hand_angle_between_vectors(v1, v2, 3, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: perpendicular 3D vectors")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, TOL, "Perpendicular 3D vectors")
  end subroutine test_perpendicular_vectors_3d

  !> Test arbitrary 3D rotation.
  subroutine test_arbitrary_3d_rotation()
    real(real64) :: v1(3), v2(3), signed_angle, dot_product, expected_magnitude
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 2.0_real64, 3.0_real64]
    v2 = [2.0_real64, 1.0_real64, 3.0_real64]
    selected_axes = [1, 2, 3]
    v1 = v1 / sqrt(sum(v1**2))
    v2 = v2 / sqrt(sum(v2**2))
    call clock_hand_angle_between_vectors(v1, v2, 3, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: arbitrary 3D rotation")
    dot_product = sum(v1 * v2)
    expected_magnitude = acos(max(-1.0_real64, min(1.0_real64, dot_product)))
    call assert_equal_real(abs(signed_angle), expected_magnitude, TOL, "3D arbitrary rotation magnitude")
  end subroutine test_arbitrary_3d_rotation

  ! ==================== HIGH DIMENSIONAL TESTS ====================

  !> Test basic high-dimensional vectors.
  subroutine test_high_dimensional_basic()
    real(real64) :: v1(5), v2(5), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    v2 = [0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    selected_axes = [1, 2, 3]
    call clock_hand_angle_between_vectors(v1, v2, 5, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: high-dimensional basic")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, TOL, "High-dimensional perpendicular vectors")
  end subroutine test_high_dimensional_basic

  !> Test high-dimensional vectors with specific selected axes.
  subroutine test_high_dimensional_selected_axes()
    real(real64) :: v1(7), v2(7), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    v2 = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64]
    selected_axes = [3, 5, 1]
    call clock_hand_angle_between_vectors(v1, v2, 7, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: high-dimensional selected axes")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, TOL, "High-dimensional with selected axes")
  end subroutine test_high_dimensional_selected_axes

  ! ==================== EDGE CASES ====================

  !> Test zero vectors (potential division by zero).
  subroutine test_zero_vectors()
    real(real64) :: v1(3), v2(3), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [0.0_real64, 0.0_real64, 0.0_real64]
    v2 = [1.0_real64, 0.0_real64, 0.0_real64]
    selected_axes = [1, 2, 3]
    call clock_hand_angle_between_vectors(v1, v2, 3, signed_angle, selected_axes, ierr)
    call assert_true(ierr == 0, "ierr should be zero for zero vector input")
    call assert_true(.not. (signed_angle /= signed_angle), "Zero vector should not produce NaN")
  end subroutine test_zero_vectors

  !> Test denormalized vectors (should still work).
  subroutine test_denormalized_vectors()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [100.0_real64, 0.0_real64]
    v2 = [0.0_real64, 50.0_real64]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: denormalized vectors")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, TOL, "Denormalized vectors")
  end subroutine test_denormalized_vectors

  !> Test tiny vectors near machine precision.
  subroutine test_tiny_vectors_precision()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    real(real64), parameter :: tiny = 1e-14_real64
    v1 = [tiny, 0.0_real64]
    v2 = [0.0_real64, tiny]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: tiny vectors")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, 1e-10_real64, "Tiny vectors precision")
  end subroutine test_tiny_vectors_precision

  !> Test huge vectors near overflow.
  subroutine test_huge_vectors_precision()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    real(real64), parameter :: huge_val = 1e14_real64
    v1 = [huge_val, 0.0_real64]
    v2 = [0.0_real64, huge_val]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: huge vectors")
    call assert_equal_real(abs(signed_angle), PI/2.0_real64, TOL, "Huge vectors precision")
  end subroutine test_huge_vectors_precision

  !> Test nearly identical vectors (precision boundary).
  subroutine test_nearly_identical_vectors()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    real(real64), parameter :: epsilon = 1e-15_real64
    v1 = [1.0_real64, 0.0_real64]
    v2 = [1.0_real64, epsilon]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: nearly identical vectors")
    call assert_true(abs(signed_angle) < 1e-10_real64, "Nearly identical vectors should have tiny angle")
  end subroutine test_nearly_identical_vectors

  !> Test nearly opposite vectors (precision boundary).
  subroutine test_nearly_opposite_vectors()
    real(real64) :: v1(2), v2(2), signed_angle
    integer :: selected_axes(3), ierr
    real(real64), parameter :: epsilon = 1e-15_real64
    v1 = [1.0_real64, 0.0_real64]
    v2 = [-1.0_real64, epsilon]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: nearly opposite vectors")
    call assert_true(abs(abs(signed_angle) - PI) < 1e-10_real64, "Nearly opposite vectors should be close to π")
  end subroutine test_nearly_opposite_vectors

  !> Test vectors with mixed positive/negative components.
  subroutine test_mixed_positive_negative()
    real(real64) :: v1(3), v2(3), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, -2.0_real64, 3.0_real64]
    v2 = [-2.0_real64, 1.0_real64, -3.0_real64]
    selected_axes = [1, 2, 3]
    v1 = v1 / sqrt(sum(v1**2))
    v2 = v2 / sqrt(sum(v2**2))
    call clock_hand_angle_between_vectors(v1, v2, 3, signed_angle, selected_axes, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: mixed positive/negative vectors")
    call assert_true(abs(signed_angle) >= 0.0_real64 .and. abs(signed_angle) <= PI, "Mixed sign vectors in valid range")
  end subroutine test_mixed_positive_negative

  ! ==================== SHIFT VECTORS TESTS ====================

  !> Test single pair of shift vectors.
  subroutine test_single_pair_shift_vectors()
    real(real64) :: origins(2, 1), targets(2, 1), signed_angles(1)
    logical :: vecs_selection_mask(1)
    integer :: selected_axes(3), ierr
    origins(:, 1) = [1.0_real64, 0.0_real64]
    targets(:, 1) = [0.0_real64, 1.0_real64]
    vecs_selection_mask = [.true.]
    selected_axes = [1, 2, 1]
    call clock_hand_angles_for_shift_vectors(origins, targets, 2, 1, vecs_selection_mask, 1, &
                                           selected_axes, signed_angles, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: single pair shift vectors")
    call assert_equal_real(abs(signed_angles(1)), PI/2.0_real64, TOL, "Single pair shift vectors")
  end subroutine test_single_pair_shift_vectors

  !> Test multiple pairs of shift vectors.
  subroutine test_multiple_pairs_shift_vectors()
    real(real64) :: origins(2, 3), targets(2, 3), signed_angles(3)
    logical :: vecs_selection_mask(3)
    integer :: selected_axes(3), ierr
    origins(:, 1) = [1.0_real64, 0.0_real64]
    targets(:, 1) = [0.0_real64, 1.0_real64]
    origins(:, 2) = [1.0_real64, 0.0_real64]
    targets(:, 2) = [-1.0_real64, 0.0_real64]
    origins(:, 3) = [1.0_real64, 0.0_real64]
    targets(:, 3) = [0.0_real64, -1.0_real64]
    vecs_selection_mask = [.true., .true., .true.]
    selected_axes = [1, 2, 1]
    call clock_hand_angles_for_shift_vectors(origins, targets, 2, 3, vecs_selection_mask, 3, selected_axes, signed_angles, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: multiple pairs shift vectors")
    call assert_equal_real(signed_angles(1), PI/2.0_real64, TOL, "First rotation (90° CCW)")
    call assert_equal_real(abs(signed_angles(2)), PI, TOL, "Second rotation (180°)")
    call assert_equal_real(signed_angles(3), -PI/2.0_real64, TOL, "Third rotation (90° CW)")
  end subroutine test_multiple_pairs_shift_vectors

  !> Test shift vectors with selection mask.
  subroutine test_shift_vectors_with_selection_mask()
    real(real64) :: origins(2, 4), targets(2, 4), signed_angles(2)
    logical :: vecs_selection_mask(4)
    integer :: selected_axes(3), ierr
    origins(:, 1) = [1.0_real64, 0.0_real64]
    targets(:, 1) = [0.0_real64, 1.0_real64]
    origins(:, 2) = [1.0_real64, 0.0_real64]
    targets(:, 2) = [-1.0_real64, 0.0_real64]
    origins(:, 3) = [1.0_real64, 0.0_real64]
    targets(:, 3) = [0.0_real64, -1.0_real64]
    origins(:, 4) = [1.0_real64, 0.0_real64]
    targets(:, 4) = [sqrt(2.0_real64)/2.0_real64, sqrt(2.0_real64)/2.0_real64]
    vecs_selection_mask = [.false., .true., .false., .true.]
    selected_axes = [1, 2, 1]
    call clock_hand_angles_for_shift_vectors(origins, targets, 2, 4, vecs_selection_mask, 2, selected_axes, signed_angles, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: shift vectors with selection mask")
    call assert_equal_real(abs(signed_angles(1)), PI, TOL, "Second vector (180°)")
    call assert_equal_real(signed_angles(2), PI/4.0_real64, TOL, "Fourth vector (45°)")
  end subroutine test_shift_vectors_with_selection_mask

  ! ==================== ERROR CONDITIONS ====================

  !> Test invalid selected axes (same indices).
  subroutine test_invalid_selected_axes()
    real(real64) :: v1(5), v2(5), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    v2 = [0.0_real64, 1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
    selected_axes = [1, 1, 1]
    call clock_hand_angle_between_vectors(v1, v2, 5, signed_angle, selected_axes, ierr)
    call assert_true(ierr /= 0, "ierr should be nonzero for invalid selected axes")
    call assert_true(.not. (signed_angle /= signed_angle), "Invalid axes should not produce NaN")
  end subroutine test_invalid_selected_axes

  !> Test out-of-bounds selected axes.
  subroutine test_out_of_bounds_selected_axes()
    real(real64) :: v1(3), v2(3), signed_angle
    integer :: selected_axes(3), ierr
    v1 = [1.0_real64, 0.0_real64, 0.0_real64]
    v2 = [0.0_real64, 1.0_real64, 0.0_real64]
    selected_axes = [1, 2, 5]
    call clock_hand_angle_between_vectors(v1, v2, 3, signed_angle, selected_axes, ierr)
    call assert_true(ierr /= 0, "ierr should be nonzero for out-of-bounds selected axes")
    call assert_true(.not. (signed_angle /= signed_angle), "Out-of-bounds axes should not produce NaN")
  end subroutine test_out_of_bounds_selected_axes

  ! ==================== PERFORMANCE TESTS ====================

  !> Test performance with large-scale data.
  subroutine test_performance_large_scale()
    integer, parameter :: n_dims = 100, n_vecs = 1000
    real(real64) :: origins(n_dims, n_vecs), targets(n_dims, n_vecs), signed_angles(n_vecs)
    logical :: vecs_selection_mask(n_vecs)
    integer :: selected_axes(3), i, ierr
    do i = 1, n_vecs
      origins(:, i) = real(i, real64) / real(n_vecs, real64)
      targets(:, i) = real(n_vecs - i + 1, real64) / real(n_vecs, real64)
    end do
    vecs_selection_mask = .true.
    selected_axes = [1, 2, 3]
    call clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_vecs, selected_axes, signed_angles, ierr)
    call assert_equal_int(ierr, 0, "ierr should be 0 for valid input: performance large scale")
    do i = 1, n_vecs
      call assert_true(abs(signed_angles(i)) <= PI, "Large-scale angles in valid range")
    end do
  end subroutine test_performance_large_scale

  ! ==================== CONSISTENCY TESTS ====================

  !> Test consistency between single and batch functions.
  subroutine test_consistency_between_functions()
    real(real64) :: v1(3), v2(3), single_angle
    real(real64) :: origins(3, 1), targets(3, 1), batch_angles(1)
    logical :: vecs_selection_mask(1)
    integer :: selected_axes(3), ierr_single, ierr_batch
    v1 = [1.0_real64, 2.0_real64, 3.0_real64]
    v2 = [3.0_real64, 2.0_real64, 1.0_real64]
    selected_axes = [1, 2, 3]
    v1 = v1 / sqrt(sum(v1**2))
    v2 = v2 / sqrt(sum(v2**2))
    call clock_hand_angle_between_vectors(v1, v2, 3, single_angle, selected_axes, ierr_single)
    origins(:, 1) = v1
    targets(:, 1) = v2
    vecs_selection_mask = [.true.]
    call clock_hand_angles_for_shift_vectors(origins, targets, 3, 1, vecs_selection_mask, 1, selected_axes, batch_angles, ierr_batch)
    call assert_equal_int(ierr_single, 0, "ierr should be 0 for valid input: single function")
    call assert_equal_int(ierr_batch, 0, "ierr should be 0 for valid input: batch function")
    call assert_equal_real(single_angle, batch_angles(1), TOL, "Single vs batch consistency")
  end subroutine test_consistency_between_functions

  !> Test mathematical properties (commutativity, etc.).
  subroutine test_mathematical_properties()
    real(real64) :: v1(2), v2(2), angle_12, angle_21
    integer :: selected_axes(3), ierr_12, ierr_21
    v1 = [1.0_real64, 0.0_real64]
    v2 = [sqrt(2.0_real64)/2.0_real64, sqrt(2.0_real64)/2.0_real64]
    selected_axes = [1, 2, 1]
    call clock_hand_angle_between_vectors(v1, v2, 2, angle_12, selected_axes, ierr_12)
    call clock_hand_angle_between_vectors(v2, v1, 2, angle_21, selected_axes, ierr_21)
    call assert_equal_int(ierr_12, 0, "ierr should be 0 for valid input: commutativity 1->2")
    call assert_equal_int(ierr_21, 0, "ierr should be 0 for valid input: commutativity 2->1")
    call assert_equal_real(angle_12, -angle_21, TOL, "Anti-commutativity of signed angles")
  end subroutine test_mathematical_properties

end module mod_test_clock_hand_angles
