!> filepath:mod_test_loess_smoothing.f90
!> Unit test suite for LOESS smoothing from f42_utils
module mod_test_loess_smoothing
  use asserts
  use f42_utils
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use mod_test_suite, only: test_case
  implicit none
  public

  
contains

  !> Get array of all available LOESS tests.
  function get_all_tests_loess_smoothing() result(all_tests)
    type(test_case),allocatable :: all_tests(:)
    allocate(all_tests(14))
    all_tests(1) = test_case("test_loess_constant_input", test_loess_constant_input)
    all_tests(2) = test_case("test_loess_linear_trend", test_loess_linear_trend)
    all_tests(3) = test_case("test_loess_outlier_suppression", test_loess_outlier_suppression)
    all_tests(4) = test_case("test_loess_sparse_fallback", test_loess_sparse_fallback)
    all_tests(5) = test_case("test_loess_single_point", test_loess_single_point)
    all_tests(6) = test_case("test_loess_identical_points", test_loess_identical_points)
    all_tests(7) = test_case("test_loess_linear_interp", test_loess_linear_interp)
    all_tests(8) = test_case("test_loess_weight_decay", test_loess_weight_decay)
    all_tests(9) = test_case("test_loess_mask_exclusion", test_loess_mask_exclusion)
    all_tests(10) = test_case("test_loess_fallback", test_loess_fallback)
    all_tests(11) = test_case("test_loess_edge_query", test_loess_edge_query)
    all_tests(12) = test_case("test_loess_invalid_dimensions", test_loess_invalid_dimensions)
    all_tests(13) = test_case("test_loess_invalid_parameters", test_loess_invalid_parameters)
    all_tests(14) = test_case("test_loess_invalid_indices", test_loess_invalid_indices)
  end function get_all_tests_loess_smoothing

  
  
  
  !> Test LOESS smoothing with constant input
  subroutine test_loess_constant_input()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(n/2), y_out(n/2)
    integer(int32) :: indices_used(n), i, ierr
    x_ref = 5.0_real64; y_ref = 10.0_real64; x_query = 5.0_real64; indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, n/2, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Constant input error check')
    call assert_true(all(abs(y_out - 10.0_real64) < 1.0e-6_real64), 'Constant input test')
  end subroutine
  
  !> Test LOESS smoothing with a linear trend 
  subroutine test_loess_linear_trend()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(n/2), y_out(n/2)
    integer(int32) :: indices_used(n), i, ierr
    x_ref = [(real(i, kind=real64), i = 1, n)]
    y_ref = 0.5_real64 * x_ref
    x_query = [(real(i, kind=real64) + 0.5_real64, i = 1, n/2)]
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, n/2, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Linear trend error check')
    call assert_true(all(abs(y_out - 0.5_real64 * x_query) < 0.05_real64), 'Linear trend test')
  end subroutine
 
  !> Test LOESS smoothing with an outlier 
  subroutine test_loess_outlier_suppression()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(50), y_out(50)
    integer(int32) :: indices_used(n), i, ierr
    x_ref(1:n-1) = 10.0_real64
    x_ref(n) = 100.0_real64
    y_ref(1:n-1) = 5.0_real64
    y_ref(n) = 99.0_real64
    x_query = 10.0_real64
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, 50, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Outlier suppression error check')
    call assert_true(all(abs(y_out - 5.0_real64) < 0.01_real64), 'Outlier suppression test')
  end subroutine
  
  !> Test LOESS smoothing with sparse data 
  subroutine test_loess_sparse_fallback()
    integer(int32), parameter :: n = 100
    real(real64) :: x_ref(n), y_ref(n), x_query(50), y_out(50)
    integer(int32) :: indices_used(n), i, ierr
    x_ref = [(real(i, kind=real64) * 100.0_real64, i = 1, n)]
    y_ref = [(real(i, kind=real64), i = 1, n)]
    x_query = [(real(i, kind=real64) * 100.0_real64 + 50.0_real64, i = 1, 50)]
    indices_used = [(i, i = 1, n)]
    call loess_smooth_2d(n, 50, x_ref, y_ref, indices_used, n, x_query, 1.0_real64, 3.0_real64, y_out, ierr)
    call assert_equal_int(ierr, 0, 'Sparse fallback error check')
    call assert_true(all(abs(y_out - y_ref(indices_used(1:50))) < 1.0e-6_real64), 'Sparse fallback test')
  end subroutine
  
  !> Test LOESS smoothing with a single point 
  subroutine test_loess_single_point()
    real(real64) :: x1(1), y1(1), xq1(1), yout1(1)
    integer(int32) :: ierr
    x1 = 0.0_real64; y1(1) = 42.0_real64; xq1 = 0.0_real64
    call loess_smooth_2d(1, 1, x1, y1, (/1/), 1, xq1, 0.1_real64, 1.0_real64, yout1, ierr)
    call assert_equal_int(ierr, 0, 'Single point error check')
    call assert_equal_real(yout1(1), y1(1), 1e-6_real64, 'Single point test')
  end subroutine
  
  !> Test LOESS smoothing with identical points 
  subroutine test_loess_identical_points()
    real(real64) :: x2(2), y2(2), xq2(1), yout2(1)
    integer(int32) :: idx2(2), ierr
    x2 = (/0.0_real64, 1.0_real64/); y2 = (/1.0_real64, 1.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth_2d(2, 1, x2, y2, idx2, 2, xq2, 0.1_real64, 1.0_real64, yout2, ierr)
    call assert_equal_int(ierr, 0, 'Identical points error check')
    call assert_equal_real(yout2(1), y2(1), 1e-6_real64, 'Identical points test')
  end subroutine

  !> Test LOESS smoothing with linear interpolation 
  subroutine test_loess_linear_interp()
    real(real64) :: x3(2), y3(2), xq3(1), yout3(1)
    integer(int32) :: idx3(2), ierr
    x3 = (/0.0_real64, 2.0_real64/); y3 = (/0.0_real64, 2.0_real64/); xq3 = 1.0_real64; idx3 = (/1,2/)
    call loess_smooth_2d(2, 1, x3, y3, idx3, 2, xq3, 0.5_real64, 3.0_real64, yout3, ierr)
    call assert_equal_int(ierr, 0, 'Linear interpolation error check')
    call assert_true(abs(yout3(1) - 1.0_real64) < 0.1_real64, 'Linear interpolation test')
  end subroutine

  !> Test LOESS smoothing with weight decay
  subroutine test_loess_weight_decay()
    real(real64) :: x2(2), y2(2), xq2(1), yout2(1)
    integer(int32) :: idx2(2), ierr
    x2 = (/0.0_real64, 10.0_real64/); y2 = (/0.0_real64, 10.0_real64/); xq2 = 0.0_real64; idx2 = (/1,2/)
    call loess_smooth_2d(2, 1, x2, y2, idx2, 2, xq2, 1.0_real64, 3.0_real64, yout2, ierr)
    call assert_equal_int(ierr, 0, 'Weight decay error check')
    call assert_true(yout2(1) < 5.0_real64, 'Weight decay test')
  end subroutine

  !> Test LOESS smoothing with mask exclusion 
  subroutine test_loess_mask_exclusion()
    real(real64) :: x3(3), y3(3), xq3(1), yout3(1)
    integer(int32) :: indices_used(1), ierr
    x3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    y3 = (/0.0_real64, 10.0_real64, 20.0_real64/)
    xq3 = 0.0_real64
    ! Pre-filter: only use index 2 (the middle point)
    indices_used = (/2/)  ! Use only the middle point
    call loess_smooth_2d(3, 1, x3, y3, indices_used, 1, xq3, 10.0_real64, 3.0_real64, yout3, ierr)
    call assert_equal_int(ierr, 0, 'Mask exclusion error check')
    call assert_true(abs(yout3(1) - 10.0_real64) < 1.0_real64, 'Mask exclusion test')
  end subroutine

  !> Test LOESS smoothing with insufficient points 
  subroutine test_loess_fallback()
    real(real64) :: x1(1), y1(1), xq1(1), yout1(1)
    integer(int32) :: ierr
    x1 = 0.0_real64; y1(1) = 123.0_real64; xq1 = 100.0_real64
    call loess_smooth_2d(1, 1, x1, y1, (/1/), 1, xq1, 0.1_real64, 1.0_real64, yout1, ierr)
    call assert_equal_int(ierr, 0, 'Fallback error check')
    call assert_equal_real(yout1(1), y1(1), 1e-6_real64, 'Fallback test')
  end subroutine
  
  !> Test edge query points (extrapolation)
  subroutine test_loess_edge_query() 
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 0, 'Edge query error check')
    call assert_true(abs(yout(1) - 10.0_real64) < 1.0 .and. abs(yout(2) - 50.0_real64) < 1.0, 'Edge query test')
  end subroutine

  !> Test LOESS smoothing with invalid dimentions
  subroutine test_loess_invalid_dimensions()
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    
    ! Test with zero dimensions
    call loess_smooth_2d(0, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 202, 'Zero n_total error check')
    
    call loess_smooth_2d(5, 0, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 202, 'Zero n_target error check')
    
    call loess_smooth_2d(5, 2, x, y, idx, 0, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 202, 'Zero n_used error check')
  end subroutine

  !> Test LOESS smoothing with Invalid parameters
  subroutine test_loess_invalid_parameters()
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    idx = (/1,2,3,4,5/)
    xq = (/0.0, 6.0/)
    
    ! Test with negative kernel_sigma
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, -1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Negative sigma error check')
    
    ! Test with negative kernel_cutoff
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, -3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Negative cutoff error check')
  end subroutine
 
  !> Test LOESS smoothing  with invalid indices
  subroutine test_loess_invalid_indices()
    real(real64) :: x(5), y(5), xq(2), yout(2)
    integer(int32) :: idx(5), ierr
    x = (/1.0, 2.0, 3.0, 4.0, 5.0/)
    y = (/10.0, 20.0, 30.0, 40.0, 50.0/)
    xq = (/0.0, 6.0/)
    
    ! Test with index out of bounds (too high)
    idx = (/1, 2, 3, 4, 6/)  ! index 6 is out of bounds for n_total=5
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Index too high error check')
    
    ! Test with index out of bounds (too low)
    idx = (/0, 2, 3, 4, 5/)  ! index 0 is out of bounds
    call loess_smooth_2d(5, 2, x, y, idx, 5, xq, 1.0_real64, 3.0_real64, yout, ierr)
    call assert_equal_int(ierr, 201, 'Index too low error check')
  end subroutine

end module mod_test_loess_smoothing
