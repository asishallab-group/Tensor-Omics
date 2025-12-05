!> @file mod_test_loess_debug.f90
!> Unit test suite for LOESS smoothing (tox_loess.F90)
!> @details Unit tests for LOESS smoothing, including masking and edge cases.

module mod_test_loess_debug
  use asserts
  use f42_utils
  use tox_loess
  use tox_errors, only: set_ok, is_ok
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

    ! Workspace size parameters
  integer(int32) :: liv, lv
  
  ! Test data
  integer(int32), allocatable :: iv(:) 
  real(real64), allocatable :: v(:)
  real(real64), allocatable :: x(:,:), y(:), w(:), xq(:,:), y_out(:)
  integer(int32) :: ierr

contains

  !> Get array of all available LOESS tests.
  subroutine get_all_tests(all_tests)
    type(test_case), intent(out) :: all_tests(15)
    all_tests(1) = test_case("test_loess_constant_input", test_loess_constant_input)
    all_tests(2) = test_case("test_loess_linear_trend", test_loess_linear_trend)
    all_tests(3) = test_case("test_loess_outlier_suppression", test_loess_outlier_suppression)
    all_tests(4) = test_case("test_simple_interpolation", test_simple_interpolation)
    all_tests(5) = test_case("test_loess_sparse_fallback", test_loess_sparse_fallback)
    all_tests(6) = test_case("test_loess_single_point", test_loess_single_point)
    all_tests(7) = test_case("test_loess_identical_points", test_loess_identical_points)
    all_tests(8) = test_case("test_loess_linear_interp", test_loess_linear_interp)
    all_tests(9) = test_case("test_loess_weight_decay", test_loess_weight_decay)
    all_tests(10) = test_case("test_loess_mask_exclusion", test_loess_mask_exclusion)
    all_tests(11) = test_case("test_loess_fallback", test_loess_fallback)
    all_tests(12) = test_case("test_loess_edge_query", test_loess_edge_query)
    all_tests(13) = test_case("test_loess_invalid_dimensions", test_loess_invalid_dimensions)
    all_tests(14) = test_case("test_loess_invalid_parameters", test_loess_invalid_parameters)
    all_tests(15) = test_case("test_loess_invalid_indices", test_loess_invalid_indices)
    
  end subroutine get_all_tests

  !> Run all LOESS smoothing tests.
  subroutine run_all_tests_loess_debug()
    type(test_case) :: all_tests(15)
    integer(int32) :: i
    call get_all_tests(all_tests)
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All LOESS smoothing tests passed successfully."
  end subroutine run_all_tests_loess_debug

  !> Run specific LOESS smoothing tests by name.
  subroutine run_named_tests_loess_debug(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(15)
    integer(int32) :: i, j
    logical :: found
    call get_all_tests(all_tests)
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
  end subroutine run_named_tests_loess_debug

  subroutine test_loess_constant_input()
    integer(int32), parameter :: d = 1, n = 100, n_pred = 50
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 2
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Constant data
    x = 5.0_real64
    y = 10.0_real64
    w = 1.0_real64
    xq = 5.0_real64
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Constant input error check')
    call assert_true(all(abs(y_out - 10.0_real64) < 1.0e-6_real64), &
                    'Constant input test - should return constant value')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_constant_input

    subroutine test_simple_interpolation()
        integer(int32), parameter :: d = 1, n = 5, n_pred = 1
        real(real64) :: span = 1.0_real64  ! Use all points
        integer(int32) :: degree = 1
        
        call setup_workspace(d, n)
        allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
        
        ! Simple 5 points
        x(1, 1) = 100.0_real64; y(1) = 1.0_real64
        x(1, 2) = 200.0_real64; y(2) = 2.0_real64
        x(1, 3) = 300.0_real64; y(3) = 3.0_real64
        x(1, 4) = 400.0_real64; y(4) = 4.0_real64
        x(1, 5) = 500.0_real64; y(5) = 5.0_real64
        w = 1.0_real64
        
        ! Query at 350 (between points 3 and 4)
        xq(1, 1) = 350.0_real64
        
        call tox_loess_predict(d, n, x, y, w, span, degree, &
                                iv, liv, v, lv, n_pred, xq, y_out, ierr)
        
        print *, "Simple test: Query at 350"
        print *, "Expected: ~3.5, Got: ", y_out(1)
        
        call cleanup_arrays()
    end subroutine test_simple_interpolation

  subroutine test_loess_linear_trend()
    integer(int32), parameter :: d = 1, n = 100, n_pred = 50
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 2
    integer(int32) :: i
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Linear trend: y = 0.5 * x
    do i = 1, n
      x(1, i) = real(i, real64)
      y(i) = 0.5_real64 * x(1, i)
    end do
    w = 1.0_real64
    
    ! Query points at midpoints
    do i = 1, n_pred
      xq(1, i) = real(i, real64) + 0.5_real64
    end do
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Linear trend error check')
    
    ! Check predictions are close to true values (allow some smoothing error)
    do i = 1, n_pred
      call assert_true(abs(y_out(i) - 0.5_real64 * xq(1, i)) < 0.05_real64, &
                      'Linear trend test - point ')
    end do
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_linear_trend

  subroutine test_loess_outlier_suppression()
    integer(int32), parameter :: d = 1, n = 100, n_pred = 50
    real(real64) :: span = 0.25_real64  ! Smaller span to better suppress outliers
    integer(int32) :: degree = 2
    integer(int32) :: i
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Data with outlier
    do i = 1, n-1
      x(1, i) = 10.0_real64
      y(i) = 5.0_real64
    end do
    x(1, n) = 100.0_real64
    y(n) = 99.0_real64  ! Outlier
    w = 1.0_real64
    
    ! Query at the non-outlier location
    xq = 10.0_real64
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Outlier suppression error check')
    
    ! Should suppress the outlier and give values close to 5
    do i = 1, n_pred
      call assert_true(abs(y_out(i) - 5.0_real64) < 0.5_real64, &
                      'Outlier suppression test - point ')
    end do
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_outlier_suppression

  subroutine test_loess_sparse_fallback()
    integer(int32), parameter :: d = 1, n = 100, n_pred = 50
    real(real64) :: span = 0.03_real64  ! Very small span to test fallback
    integer(int32) :: degree = 2
    integer(int32) :: i
    real(real64) :: expected
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Sparse data
    do i = 1, n
      x(1, i) = real(i, real64) * 100.0_real64
      y(i) = real(i, real64)
    end do
    w = 1.0_real64
    
    ! Query points in between training points
    do i = 1, n_pred
      xq(1, i) = real(i, real64) * 100.0_real64 + 50.0_real64
    end do
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Sparse fallback error check')
    
    ! With very small span, should get values close to nearest neighbor
    do i = 1, n_pred
        print *, 'point ', i, ' y_out: ', y_out(i)
        expected = real(i, real64) + 0.5_real64  ! Linear interpolation between i and i+1
        call assert_true(abs(y_out(i) - expected) < 1.0_real64, &
                        'Sparse fallback test - point ' // &
                        ' got ' // trim(real_to_str(y_out(i))) // &
                        ' expected ~' // trim(real_to_str(expected)))
    end do
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_sparse_fallback

  subroutine test_loess_single_point()
    integer(int32), parameter :: d = 1, n = 1, n_pred = 1
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 2
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Single point
    x(1, 1) = 0.0_real64
    y(1) = 42.0_real64
    w(1) = 1.0_real64
    xq(1, 1) = 0.0_real64
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 201, 'Single point error check')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_single_point

  subroutine test_loess_identical_points()
    integer(int32), parameter :: d = 1, n = 2, n_pred = 1
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 1
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Identical y values
    x(1, 1) = 0.0_real64
    x(1, 2) = 1.0_real64
    y(1) = 1.0_real64
    y(2) = 1.0_real64
    w = 1.0_real64
    xq(1, 1) = 0.0_real64
    
    call set_ok(ierr)
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Identical points error check')
    call assert_true(abs(y_out(1) - 1.0_real64) < 1e-6_real64, 'Identical points test')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_identical_points

  subroutine test_loess_linear_interp()
    integer(int32), parameter :: d = 1, n = 2, n_pred = 1
    real(real64) :: span = 1.0_real64  ! Use all points
    integer(int32) :: degree = 1  ! Linear interpolation
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Linear data
    x(1, 1) = 0.0_real64
    x(1, 2) = 2.0_real64
    y(1) = 0.0_real64
    y(2) = 2.0_real64
    w = 1.0_real64
    xq(1, 1) = 1.0_real64  ! Midpoint
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    print *, y_out
    call assert_equal_int(ierr, 0, 'Linear interpolation error check')
    call assert_true(abs(y_out(1) - 1.0_real64) < 0.1_real64, 'Linear interpolation test')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_linear_interp

  subroutine test_loess_weight_decay()
    integer(int32), parameter :: d = 1, n = 2, n_pred = 1
    real(real64) :: span = 1.0_real64
    integer(int32) :: degree = 1
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Linear data
    x(1, 1) = 0.0_real64
    x(1, 2) = 10.0_real64
    y(1) = 0.0_real64
    y(2) = 10.0_real64
    w = 1.0_real64
    xq(1, 1) = 0.0_real64  ! Query at left endpoint
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Weight decay error check')
    
    ! With tricube weights, point at x=0 should have more weight than point at x=10
    ! So prediction should be less than simple average (5.0)
    call assert_true(y_out(1) < 5.0_real64, 'Weight decay test')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_weight_decay

  subroutine test_loess_mask_exclusion()
    integer(int32), parameter :: d = 1, n = 3, n_pred = 1
    real(real64) :: span = 0.5_real64
    integer(int32) :: degree = 2
    integer(int32) :: i
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Linear data
    x(1, 1) = 0.0_real64
    x(1, 2) = 10.0_real64
    x(1, 3) = 20.0_real64
    y(1) = 0.0_real64
    y(2) = 10.0_real64
    y(3) = 20.0_real64
    
    ! Use weights to mask - only use middle point
    w(1) = 0.0_real64  ! Mask out first point
    w(2) = 1.0_real64  ! Use middle point
    w(3) = 0.0_real64  ! Mask out last point
    
    xq(1, 1) = 0.0_real64
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Mask exclusion error check')
    
    ! With only middle point weighted, prediction at x=0 should be close to y=10
    ! (though with very small span, might get exactly 10)
    call assert_true(abs(y_out(1) - 10.0_real64) < 1.0_real64, 'Mask exclusion test')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_mask_exclusion

  subroutine test_loess_fallback()
    integer(int32), parameter :: d = 1, n = 1, n_pred = 1
    real(real64) :: span = 0.1_real64
    integer(int32) :: degree = 2
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Single point far from query
    x(1, 1) = 0.0_real64
    y(1) = 123.0_real64
    w(1) = 1.0_real64
    xq(1, 1) = 100.0_real64  ! Far from training point
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 201, 'Fallback error check')
    !! One point is not enough for degree 2 fit
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_fallback

  subroutine test_loess_edge_query()
    integer(int32), parameter :: d = 1, n = 5, n_pred = 2
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 2
    integer(int32) :: i
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Linear data
    do i = 1, n
      x(1, i) = real(i, real64)
      y(i) = real(i, real64) * 10.0_real64
    end do
    w = 1.0_real64
    
    ! Query at edges (extrapolation)
    xq(1, 1) = 0.0_real64  ! Left of data
    xq(1, 2) = 6.0_real64  ! Right of data
    
    ! Call LOESS
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Assertions
    call assert_equal_int(ierr, 0, 'Edge query error check')
    
    ! Edge predictions should be reasonable extrapolations
    call assert_true(y_out(1) < 15.0_real64, 'Left edge query test')
    call assert_true(y_out(2) > 45.0_real64, 'Right edge query test')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_edge_query

  ! Error handling tests
  subroutine test_loess_invalid_dimensions()
    integer(int32), parameter :: d = 1, n = 5, n_pred = 2
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 2
    
    ! Setup
    call setup_workspace(d, n)
    
    ! Allocate arrays
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Valid data
    x = 1.0_real64
    y = 10.0_real64
    w = 1.0_real64
    xq = 0.0_real64
    
    ! Test with zero dimensions - should handle gracefully
    ! Note: Implementation may or may not check this
    ! We'll test that it doesn't crash
    
    ! Cleanup
    call cleanup_arrays()
    
    print *, "test_loess_invalid_dimensions: No specific error codes in this implementation"
    
  end subroutine test_loess_invalid_dimensions

  subroutine test_loess_invalid_parameters()
    integer(int32), parameter :: d = 1, n = 5, n_pred = 2
    real(real64) :: span
    integer(int32) :: degree
    
    ! Setup workspace for valid case first
    call setup_workspace(d, n)
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Valid data
    x = 1.0_real64
    y = 10.0_real64
    w = 1.0_real64
    xq = 0.0_real64
    
    ! Test with invalid span (negative)
    span = -1.0_real64
    degree = 2
    
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Should handle negative span (implementation clips to EPS)
    call assert_true(ierr == 0, 'Negative span handled')
    
    ! Test with invalid degree (should use default or handle gracefully)
    span = 0.75_real64
    degree = 5  ! Invalid degree
    
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    ! Implementation may handle this or use default
    call assert_true(ierr == 201, 'Invalid degree handled')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_invalid_parameters

  subroutine test_loess_invalid_indices()
    ! Note: This implementation doesn't use explicit index arrays like the old one
    ! Instead, it uses weights to mask points. So we test weight-based masking.
    integer(int32), parameter :: d = 1, n = 5, n_pred = 2
    real(real64) :: span = 0.75_real64
    integer(int32) :: degree = 2
    
    ! Setup
    call setup_workspace(d, n)
    allocate(x(d, n), y(n), w(n), xq(d, n_pred), y_out(n_pred))
    
    ! Valid data with some zero weights (masked points)
    x(1, 1) = 1.0_real64
    x(1, 2) = 2.0_real64
    x(1, 3) = 3.0_real64
    x(1, 4) = 4.0_real64
    x(1, 5) = 5.0_real64
    
    y = 10.0_real64
    w = 1.0_real64
    w(5) = 0.0_real64  ! Mask last point
    
    xq = 3.0_real64
    
    ! Call LOESS - should handle zero weights gracefully
    call tox_loess_predict(d, n, x, y, w, span, degree, &
                          iv, liv, v, lv, n_pred, xq, y_out, ierr)
    
    call assert_equal_int(ierr, 0, 'Zero weights handled')
    
    ! Cleanup
    call cleanup_arrays()
    
  end subroutine test_loess_invalid_indices

  ! Helper subroutines
  subroutine setup_workspace(d, n)
    integer(int32), intent(in) :: d, n
    
    ! Get required workspace sizes
    call tox_loess_required_workspace(d, n, liv, lv)
    
    ! Allocate workspace
    if (allocated(iv)) deallocate(iv)
    if (allocated(v)) deallocate(v)
    
    allocate(iv(liv), v(lv))
    
    ! Initialize workspace (optional, but good practice)
    iv = 0
    v = 0.0_real64
    
  end subroutine setup_workspace
  
  subroutine cleanup_arrays()
    if (allocated(x)) deallocate(x)
    if (allocated(y)) deallocate(y)
    if (allocated(w)) deallocate(w)
    if (allocated(xq)) deallocate(xq)
    if (allocated(y_out)) deallocate(y_out)
  end subroutine cleanup_arrays

    function real_to_str(r) result(str)
        real(real64), intent(in) :: r
        character(len=32) :: str
        write(str, '(F15.6)') r
        str = adjustl(str)
    end function real_to_str

end module mod_test_loess_debug