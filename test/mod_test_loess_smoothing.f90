!> @file mod_test_loess_smoothing.f90
!> Unit test suite for LOESS smoothing from tox_loess
!> @details Unit tests for LOESS smoothing, including edge cases and input validation.
!

module mod_test_loess_smoothing
  use asserts
  use tox_errors, only: ERR_INVALID_INPUT, ERR_NAN_INF, ERR_SIZE_MISMATCH, ERR_EMPTY_INPUT
  use tox_loess,  only: loess_alloc
  use f42_utils
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use test_suite
  implicit none
  public

contains

  !> Get array of all available LOESS tests.
  function get_all_tests_loess_smoothing() result(all_tests)
    type(test_case), allocatable :: all_tests(:)
    allocate(all_tests(22))

    all_tests(1)  = test_case("test_loess_constant_input",       test_loess_constant_input)
    all_tests(2)  = test_case("test_loess_linear_trend",         test_loess_linear_trend)
    all_tests(3)  = test_case("test_loess_outlier_suppression",  test_loess_outlier_suppression)
    all_tests(4)  = test_case("test_loess_sparse_like_spacing",  test_loess_sparse_like_spacing)
    all_tests(5)  = test_case("test_loess_single_point",         test_loess_single_point)
    all_tests(6)  = test_case("test_loess_identical_points",     test_loess_identical_points)
    all_tests(7)  = test_case("test_loess_linear_interp",        test_loess_linear_interp)
    all_tests(8)  = test_case("test_loess_edge_query",           test_loess_edge_query)
    all_tests(9) = test_case("test_loess_invalid_span",         test_loess_invalid_span)
    all_tests(10) = test_case("test_loess_degree_0_constant",    test_loess_degree_0_constant)
    all_tests(11) = test_case("test_loess_robust_same_as_plain", test_loess_robust_same_as_plain)
    all_tests(12) = test_case("test_loess_size_mismatch", test_loess_size_mismatch)
    all_tests(13) = test_case("test_loess_empty_input", test_loess_empty_input)
    all_tests(14) = test_case("test_loess_invalid_mode", test_loess_invalid_mode)
    all_tests(15) = test_case("test_loess_invalid_degree", test_loess_invalid_degree)
    all_tests(16) = test_case("test_loess_span_too_small", test_loess_span_too_small)
    all_tests(17) = test_case("test_loess_effective_points_guard", test_loess_effective_points_guard)
    all_tests(18) = test_case("test_loess_robust_requires_iters", test_loess_robust_requires_iters)
    all_tests(19) = test_case("test_loess_nan_in_x", test_loess_nan_in_x)
    all_tests(20) = test_case("test_loess_degenerate_x_range_fallback", test_loess_degenerate_x_range_fallback)
    all_tests(21) = test_case("test_loess_insufficient_unique_x_fallback", test_loess_insufficient_unique_x_fallback)
    all_tests(22) = test_case("test_loess_many_points_low_variance_y", test_loess_many_points_low_variance_y)

  end function get_all_tests_loess_smoothing

  ! ---------------------------------------------------------------------------
  ! TESTS
  ! ---------------------------------------------------------------------------

  !> Constant input: y should remain constant everywhere.
  subroutine test_loess_constant_input()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr, i

    do i = 1, n
      x_ref(i) = real(i, real64) 
    end do
    y_ref = 10.0_real64

    call loess_alloc(x_ref, y_ref, 0.5_real64, 1_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, 0, "Constant input ierr==0")
    call assert_true(all(abs(yhat - 10.0_real64) < 1.0e-6_real64), "Constant input yhat constant")
  end subroutine test_loess_constant_input

  !> Linear trend: y = 0.5 x, should be approximated well.
  subroutine test_loess_linear_trend()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr, i

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    y_ref = 0.5_real64 * x_ref

    call loess_alloc(x_ref, y_ref, 0.7_real64, 1_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, 0, "Linear trend ierr==0")
    call assert_true(all(abs(yhat - y_ref) < 0.05_real64), "Linear trend approx ok")
  end subroutine test_loess_linear_trend

  !> Robust outlier suppression: one strong outlier should be down-weighted.
  subroutine test_loess_outlier_suppression()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr,i 

    do i = 1, n
      x_ref(i) = real(i, real64)  
      y_ref(i) = 5.0_real64        
    end do

    ! One outlier
    y_ref(50) = 1000.0_real64

    call loess_alloc(x_ref, y_ref, 0.7_real64, 1_int32, yhat, 1, 3, ierr)

    call assert_equal_int(ierr, 0, "Robust outlier suppression ierr==0")
    call assert_true(all(abs(yhat - 5.0_real64) < 0.05_real64), "Outlier suppressed (near 5)")
  end subroutine test_loess_outlier_suppression

  !> Sparse-like spacing: large gaps in x should still produce finite outputs.
  subroutine test_loess_sparse_like_spacing()
    integer(int32), parameter :: n = 80
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr, i

    do i = 1, n
      x_ref(i) = 100.0_real64 * real(i, real64)
      y_ref(i) = real(i, real64)
    end do

    call loess_alloc(x_ref, y_ref, 0.7_real64, 1_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, 0, "Sparse-like spacing ierr==0")
    call assert_true(all(.not. (yhat /= yhat)), "Sparse-like spacing: no NaNs")  ! NaN check: yhat==yhat
  end subroutine test_loess_sparse_like_spacing

  !> Single point: loess can't smooth a single point, therefore it should return an error.
  subroutine test_loess_single_point()
    integer(int32), parameter :: n = 1
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr

    x_ref(1) = 0.0_real64
    y_ref(1) = 42.0_real64

    call loess_alloc(x_ref, y_ref, 1.0_real64, 1_int32, yhat, 0, 0, ierr)
  call assert_equal_int(ierr, 0, "n=1 should return OK (identity fallback)")
  call assert_equal_real(yhat(1), y_ref(1), EPS, "n=1 should return yhat=y")
  end subroutine test_loess_single_point

  !> Identical Y values: y constant, output should be constant.
  subroutine test_loess_identical_points()
    integer(int32), parameter :: n = 10
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr, i

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    
    y_ref = 1.0_real64

    call loess_alloc(x_ref, y_ref, 1.0_real64, 1_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, 0, "Constant Y ierr==0")
    ! Verificamos que yhat sea 1.0 en todos los puntos
    call assert_true(all(abs(yhat - 1.0_real64) < 1.0e-10_real64), "Output remains constant")
  end subroutine test_loess_identical_points

  !> Linear interpolation-ish sanity: with two points (0,0) and (2,2), query at 1 should be near 1.
  subroutine test_loess_linear_interp()
    integer(int32), parameter :: n = 4
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr

    x_ref = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]
    y_ref = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]

    call loess_alloc(x_ref, y_ref, 1.0_real64, 1_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, 0, "Linear interp ierr==0")
    call assert_true(abs(yhat(2) - 1.0_real64) < 1.0e-6_real64, "Linear trend preserved")
  end subroutine test_loess_linear_interp


  !> Edge queries: Expect near boundary values 
  subroutine test_loess_edge_query()
    integer(int32), parameter :: n = 10 
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr, i

    do i = 1, n
      x_ref(i) = real(i, real64)
      y_ref(i) = real(i, real64) * 10.0_real64 ! y = 10, 20, ..., 100
    end do

    call loess_alloc(x_ref, y_ref, 0.5_real64, 1_int32, yhat, 0, 0, ierr)
    call assert_equal_int(ierr, 0, "Edge query ierr==0")
    call assert_true(abs(yhat(1) - 10.0_real64) < 2.0_real64, "Low edge (x=1) near 10")
    call assert_true(abs(yhat(10) - 100.0_real64) < 2.0_real64, "High edge (x=10) near 100")
  end subroutine test_loess_edge_query

  !> Invalid span: negative span should trigger ERR_INVALID_INPUT.
  subroutine test_loess_invalid_span()
    integer(int32), parameter :: n = 5
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr

    x_ref = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
    y_ref = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]

    call loess_alloc(x_ref, y_ref, -1.0_real64, 1_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, ERR_INVALID_INPUT, "Invalid span -> ERR_INVALID_INPUT")
  end subroutine test_loess_invalid_span

  !> Degree 0 with constant y: should still return constant.
  subroutine test_loess_degree_0_constant()
    integer(int32), parameter :: n = 50
    real(real64) :: x_ref(n), y_ref(n), yhat(n)
    integer(int32) :: ierr, i

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    y_ref = 7.0_real64


    call loess_alloc(x_ref, y_ref, 0.7_real64, 0_int32, yhat, 0, 0, ierr)

    call assert_equal_int(ierr, 0, "Degree 0 constant ierr==0")
    call assert_true(all(abs(yhat - 7.0_real64) < 1.0e-6_real64), "Degree 0 constant preserved")
  end subroutine test_loess_degree_0_constant

  !> Robust equals plain in a clean dataset (no outliers): results should be close.
  subroutine test_loess_robust_same_as_plain()
    integer(int32), parameter :: n = 60
    real(real64) :: x_ref(n), y_ref(n), yhat_plain(n), yhat_robust(n)
    integer(int32) :: ierr1, ierr2, i

    do i = 1, n
      x_ref(i) = real(i, real64)
    end do
    y_ref = 2.0_real64 * x_ref + 1.0_real64

    call loess_alloc(x_ref, y_ref, 0.7_real64, 1_int32, yhat_plain, 0, 0, ierr1)
    call loess_alloc(x_ref, y_ref, 0.7_real64, 1_int32, yhat_robust, 1, 3, ierr2)

    call assert_equal_int(ierr1, 0, "Plain ierr==0")
    call assert_equal_int(ierr2, 0, "Robust ierr==0")
    call assert_true(all(abs(yhat_plain - yhat_robust) < 0.1_real64), "Robust ~ Plain (clean data)")
  end subroutine test_loess_robust_same_as_plain

  subroutine test_loess_size_mismatch()
    integer(int32), parameter :: nx=5, ny=6
    real(real64) :: x(nx), y(ny), yhat(ny)
    integer(int32) :: ierr, i

    x = [(real(i,real64), i=1,nx)]
    y = [(real(i,real64), i=1,ny)]

    call loess_alloc(x, y, 0.5_real64, 1_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_SIZE_MISMATCH, "size(x)!=size(y) should error")
  end subroutine

  subroutine test_loess_empty_input()
    real(real64), allocatable :: x(:), y(:), yhat(:)
    integer(int32) :: ierr

    allocate(x(0), y(0), yhat(0))
    call loess_alloc(x, y, 0.5_real64, 1_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_EMPTY_INPUT, "n=0 should return ERR_EMPTY_INPUT")
  end subroutine

  subroutine test_loess_invalid_mode()
    integer(int32), parameter :: n=10
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr,i

    x = [(real(i,real64), i=1,n)]
    y = [(real(i,real64), i=1,n)]

    call loess_alloc(x, y, 0.5_real64, 1_int32, yhat, 2_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "mode not in {0,1} should error")
  end subroutine

  subroutine test_loess_invalid_degree()
    integer(int32), parameter :: n=20
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr,i

    x = [(real(i,real64), i=1,n)]
    y = [(real(i,real64), i=1,n)]

    call loess_alloc(x, y, 0.5_real64, 3_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "degree>2 should error")
  end subroutine

  subroutine test_loess_span_too_small()
    integer(int32), parameter :: n=50
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr,i

    x = [(real(i,real64), i=1,n)]
    y = sin(x)

    call loess_alloc(x, y, 0.0_real64, 1_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "span<=EPS should error")
  end subroutine

  subroutine test_loess_effective_points_guard()
    integer(int32), parameter :: n=100
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr,i

    x = [(real(i,real64), i=1,n)]
    y = sin(x)

    ! degree=2 => requiere n_eff >= 5. span=0.01 => n_eff=1
    call loess_alloc(x, y, 0.01_real64, 2_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "n_eff too small should error")
  end subroutine

  subroutine test_loess_robust_requires_iters()
    integer(int32), parameter :: n=30
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr,i

    x = [(real(i,real64), i=1,n)]
    y = sin(x)

    call loess_alloc(x, y, 0.5_real64, 1_int32, yhat, 1_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_INVALID_INPUT, "robust with n_iters<1 should error")
  end subroutine

  subroutine test_loess_nan_in_x()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    integer(int32), parameter :: n=10
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr,i

    x = [(real(i,real64), i=1,n)]
    y = [(real(i,real64), i=1,n)]
    x(5) = ieee_value(x(5), ieee_quiet_nan)

    call loess_alloc(x, y, 0.5_real64, 1_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, ERR_NAN_INF, "NaN in x should return ERR_NAN_INF")
  end subroutine

  subroutine test_loess_degenerate_x_range_fallback()
    integer(int32), parameter :: n=50
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr
    integer(int32) :: i

    x = 1.0_real64
    do i=1,n
      y(i) = real(i, real64)
    end do

    call loess_alloc(x, y, 0.5_real64, 1_int32, yhat, 0_int32, 0_int32, ierr)

    call assert_equal_int(ierr, 0, "degenerate x should return OK fallback")
    call assert_true(maxval(abs(yhat - y)) == 0.0_real64, "fallback should return yhat=y")
  end subroutine

  subroutine test_loess_insufficient_unique_x_fallback()
    integer(int32), parameter :: n=60
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr, i

    do i=1,n
      select case (mod(i,3))
      case (0); x(i) = 0.0_real64
      case (1); x(i) = 1.0_real64
      case (2); x(i) = 2.0_real64
      end select
      y(i) = sin(real(i,real64))
    end do

    call loess_alloc(x, y, 0.8_real64, 2_int32, yhat, 0_int32, 0_int32, ierr)

    call assert_equal_int(ierr, 0, "insufficient unique x should OK-fallback")
    call assert_true(maxval(abs(yhat - y)) == 0.0_real64, "fallback should return yhat=y")
  end subroutine

  subroutine test_loess_many_points_low_variance_y()
    integer(int32), parameter :: n=200
    real(real64) :: x(n), y(n), yhat(n)
    integer(int32) :: ierr, i
    real(real64) :: maxdiff

    do i=1,n
      x(i) = real(i, real64) / real(n, real64)   ! [0,1]
      y(i) = 1.0_real64 + 1.0e-10_real64 * sin(100.0_real64 * x(i))
    end do

    call loess_alloc(x, y, 0.3_real64, 1_int32, yhat, 0_int32, 0_int32, ierr)
    call assert_equal_int(ierr, 0, "low-variance y should still succeed")

    ! yhat debería estar muy cerca de y (misma escala)
    maxdiff = maxval(abs(yhat - y))
    call assert_true(maxdiff < 1.0e-8_real64, "yhat should remain close for low-variance signal")
  end subroutine


end module mod_test_loess_smoothing
