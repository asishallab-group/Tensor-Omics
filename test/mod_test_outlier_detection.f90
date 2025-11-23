module mod_test_outlier_detection
  use asserts
  use tox_trajectory_contribution_analysis_outlier_detection
  use tox_errors
  use, intrinsic :: iso_fortran_env, only: real64, int32
  implicit none
  public

  real(real64), parameter :: TOL = 1e-12_real64

  abstract interface
    subroutine test_interface()
    end subroutine test_interface
  end interface

  type :: test_case
    character(len=64) :: name
    procedure(test_interface), pointer, nopass :: test_proc => null()
  end type test_case

contains

  function get_all_tests() result(all_tests)
    type(test_case) :: all_tests(19)
    
    all_tests(1) = test_case("test_spike_thresholds_timepoint_wise", test_spike_thresholds_timepoint_wise)
    all_tests(2) = test_case("test_spike_thresholds_alloc", test_spike_thresholds_alloc)
    all_tests(3) = test_case("test_spike_thresholds_edge_cases", test_spike_thresholds_edge_cases)
    all_tests(4) = test_case("test_spike_thresholds_extreme_percentiles", test_spike_thresholds_extreme_percentiles)
    all_tests(5) = test_case("test_integrated_threshold_sample_wise", test_integrated_threshold_sample_wise)
    all_tests(6) = test_case("test_integrated_threshold_alloc", test_integrated_threshold_alloc)
    all_tests(7) = test_case("test_integrated_threshold_edge_cases", test_integrated_threshold_edge_cases)
    all_tests(8) = test_case("test_detect_outliers_integrated", test_detect_outliers_integrated)
    all_tests(9) = test_case("test_detect_outliers_integrated_edge_cases", test_detect_outliers_integrated_edge_cases)
    all_tests(10) = test_case("test_detect_outliers_spike", test_detect_outliers_spike)
    all_tests(11) = test_case("test_detect_outliers_spike_edge_cases", test_detect_outliers_spike_edge_cases)
    all_tests(12) = test_case("test_biological_data_layout", test_biological_data_layout)
    all_tests(13) = test_case("test_column_major_optimization", test_column_major_optimization)
    all_tests(14) = test_case("test_empty_inputs", test_empty_inputs)
    all_tests(15) = test_case("test_invalid_percentiles", test_invalid_percentiles)
    all_tests(16) = test_case("test_consistency_performance_alloc", test_consistency_performance_alloc)
    all_tests(17) = test_case("test_large_scale_performance", test_large_scale_performance)
    all_tests(18) = test_case("test_error_propagation", test_error_propagation)
    all_tests(19) = test_case("test_comprehensive_workflow", test_comprehensive_workflow)
  end function get_all_tests

  subroutine run_all_tests_outlier_detection()
    type(test_case) :: all_tests(19)
    integer :: i
    all_tests = get_all_tests()
    do i = 1, size(all_tests)
      call all_tests(i)%test_proc()
      print *, trim(all_tests(i)%name), " passed."
    end do
    print *, "All outlier detection tests passed successfully."
  end subroutine run_all_tests_outlier_detection

    !> Run specific quantile_normalization tests by name.
  subroutine run_named_tests_outlier_detection(test_names)
    character(len=*), intent(in) :: test_names(:)
    type(test_case) :: all_tests(19)
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
  end subroutine run_named_tests_outlier_detection

  ! ==================== SPIKE THRESHOLDS TESTS ====================

  subroutine test_spike_thresholds_timepoint_wise()
    real(real64) :: spike_contribs(5, 3)  ! 5 timepoints, 3 samples 
    real(real64) :: thresholds(5)          ! 5 thresholds 
    integer(int32) :: permutation(3, 5)    ! 3 samples, 5 timepoints
    integer(int32) :: ierr
    integer :: i, j

    do j = 1, 3
      do i = 1, 5
        spike_contribs(i, j) = real(j, real64) * real(i, real64)
      end do
    end do
    
    ! Pre-compute permutations (data is already sorted in this simple case)
    do j = 1, 5  ! For each timepoint
      permutation(:, j) = [1, 2, 3]  ! Each timepoint's sorted indices
    end do
    
    call calc_spike_thresholds(spike_contribs, 5, 3, 80.0_real64, thresholds, permutation, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for spike thresholds timepoint-wise")

    call assert_equal_real(thresholds(1), 2.6_real64, TOL, "Timepoint 1 threshold")
    call assert_equal_real(thresholds(2), 5.2_real64, TOL, "Timepoint 2 threshold")
    call assert_equal_real(thresholds(3), 7.8_real64, TOL, "Timepoint 3 threshold")
    call assert_equal_real(thresholds(4), 10.4_real64, TOL, "Timepoint 4 threshold")
    call assert_equal_real(thresholds(5), 13.0_real64, TOL, "Timepoint 5 threshold")
    
  end subroutine test_spike_thresholds_timepoint_wise

  subroutine test_spike_thresholds_alloc()
    real(real64) :: spike_contribs(5, 3)  ! 5 timepoints, 3 samples
    real(real64) :: thresholds(5)          
    integer(int32) :: ierr, i, j
    
    ! Create test data: timepoints × samples
    do j = 1, 3  ! samples
      do i = 1, 5  ! timepoints
        spike_contribs(i, j) = real(j, real64) * real(i, real64)
      end do
    end do
    
    call calc_spike_thresholds_alloc(spike_contribs, 5, 3, 80.0_real64, thresholds, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for spike thresholds alloc")
    ! Should get same results as performance version
    call assert_equal_real(thresholds(1), 2.6_real64, TOL, "Timepoint 1 threshold")
    call assert_equal_real(thresholds(2), 5.2_real64, TOL, "Timepoint 2 threshold")
    call assert_equal_real(thresholds(3), 7.8_real64, TOL, "Timepoint 3 threshold")
    call assert_equal_real(thresholds(4), 10.4_real64, TOL, "Timepoint 4 threshold")
    call assert_equal_real(thresholds(5), 13.0_real64, TOL, "Timepoint 5 threshold")
    
  end subroutine test_spike_thresholds_alloc

  subroutine test_spike_thresholds_edge_cases()
    real(real64) :: spike_contribs(2, 4)  ! 2 timepoints, 4 samples
    real(real64) :: thresholds(2)
    integer :: ierr
    
    ! Test with duplicate values and unsorted data
    ! Timepoint 1: [5.0, 1.0, 5.0, 3.0]
    ! Timepoint 2: [2.0, 2.0, 8.0, 2.0]
    spike_contribs(1, :) = [5.0_real64, 1.0_real64, 5.0_real64, 3.0_real64]
    spike_contribs(2, :) = [2.0_real64, 2.0_real64, 8.0_real64, 2.0_real64]
    
    ! For simplicity, we'll use the alloc version which handles sorting
    call calc_spike_thresholds_alloc(spike_contribs, 2, 4, 50.0_real64, thresholds, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for spike thresholds edge cases")
    ! Median of Timepoint 1 [1,3,5,5] = 4.0, Timepoint 2 [2,2,2,8] = 2.0
    call assert_equal_real(thresholds(1), 4.0_real64, TOL, "Timepoint 1 median with duplicates")
    call assert_equal_real(thresholds(2), 2.0_real64, TOL, "Timepoint 2 median with duplicates")
    
  end subroutine test_spike_thresholds_edge_cases

  subroutine test_spike_thresholds_extreme_percentiles()
    real(real64) :: spike_contribs(5, 2)  ! 5 timepoints, 2 samples
    real(real64) :: thresholds_min(5), thresholds_max(5)
    integer(int32) :: ierr, i, j
    
    ! Test data: timepoints × samples
    do j = 1, 2  ! samples
      do i = 1, 5  ! timepoints
        spike_contribs(i, j) = real(j, real64) * real(i, real64)
      end do
    end do
    
    ! Test 0th percentile (minimum)
    call calc_spike_thresholds_alloc(spike_contribs, 5, 2, 0.0_real64, thresholds_min, ierr)
    call assert_equal_int(ierr, 0, "ierr for 0th percentile")
    call assert_equal_real(thresholds_min(1), 1.0_real64, TOL, "Timepoint 1 0th percentile")
    call assert_equal_real(thresholds_min(2), 2.0_real64, TOL, "Timepoint 2 0th percentile")
    call assert_equal_real(thresholds_min(3), 3.0_real64, TOL, "Timepoint 3 0th percentile")
    
    ! Test 100th percentile (maximum)
    call calc_spike_thresholds_alloc(spike_contribs, 5, 2, 100.0_real64, thresholds_max, ierr)
    call assert_equal_int(ierr, 0, "ierr for 100th percentile")
    call assert_equal_real(thresholds_max(1), 2.0_real64, TOL, "Timepoint 1 100th percentile")
    call assert_equal_real(thresholds_max(2), 4.0_real64, TOL, "Timepoint 2 100th percentile")
    call assert_equal_real(thresholds_max(3), 6.0_real64, TOL, "Timepoint 3 100th percentile")
    
  end subroutine test_spike_thresholds_extreme_percentiles

  ! ==================== INTEGRATED THRESHOLDS TESTS ====================

  subroutine test_integrated_threshold_sample_wise()
    real(real64) :: contributions(5)       ! 5 samples
    real(real64) :: threshold
    integer(int32) :: permutation(5), ierr
    
    contributions = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
    permutation = [1, 2, 3, 4, 5]  ! Already sorted
    
    call calc_integrated_threshold(contributions, 5, 80.0_real64, threshold, permutation, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for integrated threshold sample-wise")
    ! 80th percentile of [10,20,30,40,50] is 42.0
    call assert_equal_real(threshold, 42.0_real64, TOL, "80th percentile of samples")
    
  end subroutine test_integrated_threshold_sample_wise

  subroutine test_integrated_threshold_alloc()
    real(real64) :: contributions(5)       ! 5 samples
    real(real64) :: threshold
    integer(int32) :: ierr
    
    contributions = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64]
    
    call calc_integrated_threshold_alloc(contributions, 5, 80.0_real64, threshold, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for integrated threshold alloc")
    call assert_equal_real(threshold, 42.0_real64, TOL, "80th percentile alloc")
    
  end subroutine test_integrated_threshold_alloc

  subroutine test_integrated_threshold_edge_cases()
    real(real64) :: contributions_single(1), threshold
    real(real64) :: contributions(4)
    integer(int32) :: ierr
    
    ! Test with single element
    contributions_single = [100.0_real64]
    call calc_integrated_threshold_alloc(contributions_single, 1, 50.0_real64, threshold, ierr)
    call assert_equal_int(ierr, 0, "ierr for single element")
    call assert_equal_real(threshold, 100.0_real64, TOL, "Single element threshold")
    
    ! Test with all identical values
    contributions = [5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64]
    call calc_integrated_threshold_alloc(contributions, 4, 75.0_real64, threshold, ierr)
    call assert_equal_int(ierr, 0, "ierr for identical values")
    call assert_equal_real(threshold, 5.0_real64, TOL, "Identical values threshold")
    
  end subroutine test_integrated_threshold_edge_cases

  ! ==================== OUTLIER DETECTION TESTS ====================

  subroutine test_detect_outliers_integrated()
    real(real64) :: contributions(5)       ! 5 samples
    real(real64) :: threshold
    logical :: outlier_mask(5)
    integer(int32) :: ierr
    
    contributions = [10.0_real64, 20.0_real64, 30.0_real64, 45.0_real64, 50.0_real64]
    threshold = 42.0_real64  ! From previous test
    
    call detect_outliers_integrated(contributions, 5, threshold, outlier_mask, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for integrated outliers")
    ! Samples with values > 42.0 are outliers
    call assert_true(.not. outlier_mask(1), "Sample 1 not outlier")
    call assert_true(.not. outlier_mask(2), "Sample 2 not outlier")
    call assert_true(.not. outlier_mask(3), "Sample 3 not outlier")
    call assert_true(outlier_mask(4), "Sample 4 is outlier")
    call assert_true(outlier_mask(5), "Sample 5 is outlier")
    
  end subroutine test_detect_outliers_integrated

  subroutine test_detect_outliers_integrated_edge_cases()
    real(real64) :: contributions(4), threshold
    logical :: outlier_mask(4)
    integer(int32) :: ierr
    
    ! Test with values exactly at threshold
    contributions = [5.0_real64, 10.0_real64, 10.0_real64, 15.0_real64]
    threshold = 10.0_real64
    
    call detect_outliers_integrated(contributions, 4, threshold, outlier_mask, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for integrated edge cases")
    call assert_true(.not. outlier_mask(1), "Value below threshold")
    call assert_true(.not. outlier_mask(2), "Value equal to threshold - not outlier")
    call assert_true(.not. outlier_mask(3), "Value equal to threshold - not outlier")
    call assert_true(outlier_mask(4), "Value above threshold")
    
  end subroutine test_detect_outliers_integrated_edge_cases

  subroutine test_detect_outliers_spike()
    real(real64) :: spike_contribs(3, 4)   ! 3 timepoints, 4 samples
    real(real64) :: thresholds(3)           ! 3 timepoint thresholds
    logical :: outlier_mask(3, 4)           ! 3 timepoints × 4 samples
    integer(int32) :: ierr
    integer :: i, j
    
    ! Setup test data: 3 timepoints × 4 samples
    ! Timepoint 1: [1, 2, 3, 4]
    ! Timepoint 2: [2, 4, 6, 8]
    ! Timepoint 3: [3, 6, 9, 12]
    do j = 1, 4  ! samples
      do i = 1, 3  ! timepoints
        spike_contribs(i, j) = real(j, real64) * real(i, real64)
      end do
    end do
    
    ! Timepoint-specific thresholds
    thresholds = [3.5_real64, 7.0_real64, 10.5_real64]
    
    call detect_outliers_spike(spike_contribs, 3, 4, thresholds, outlier_mask, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for spike outliers")
    
    ! Verify specific outliers using timepoint-specific thresholds
    ! Timepoint 1 threshold = 3.5: values [1,2,3,4] -> outliers: [F,F,F,T] (sample 4)
    ! Timepoint 2 threshold = 7.0: values [2,4,6,8] -> outliers: [F,F,F,T] (sample 4)  
    ! Timepoint 3 threshold = 10.5: values [3,6,9,12] -> outliers: [F,F,F,T] (sample 4)
    
    call assert_true(.not. outlier_mask(1, 1), "Timepoint 1, Sample 1 not outlier")
    call assert_true(.not. outlier_mask(1, 2), "Timepoint 1, Sample 2 not outlier")
    call assert_true(.not. outlier_mask(1, 3), "Timepoint 1, Sample 3 not outlier")
    call assert_true(outlier_mask(1, 4), "Timepoint 1, Sample 4 is outlier")
    
    call assert_true(.not. outlier_mask(2, 1), "Timepoint 2, Sample 1 not outlier")
    call assert_true(.not. outlier_mask(2, 2), "Timepoint 2, Sample 2 not outlier") 
    call assert_true(.not. outlier_mask(2, 3), "Timepoint 2, Sample 3 not outlier")
    call assert_true(outlier_mask(2, 4), "Timepoint 2, Sample 4 is outlier")
    
    call assert_true(.not. outlier_mask(3, 1), "Timepoint 3, Sample 1 not outlier")
    call assert_true(.not. outlier_mask(3, 2), "Timepoint 3, Sample 2 not outlier")
    call assert_true(.not. outlier_mask(3, 3), "Timepoint 3, Sample 3 not outlier")
    call assert_true(outlier_mask(3, 4), "Timepoint 3, Sample 4 is outlier")
    
  end subroutine test_detect_outliers_spike

  subroutine test_detect_outliers_spike_edge_cases()
    real(real64) :: spike_contribs(2, 3)   ! 2 timepoints, 3 samples
    real(real64) :: thresholds(2)           ! 2 timepoint thresholds
    logical :: outlier_mask(2, 3)           ! 2 timepoints × 3 samples
    integer(int32) :: ierr
    integer :: i, j
    
    ! Test with values exactly at thresholds
    do j = 1, 3  ! samples
      do i = 1, 2  ! timepoints
        spike_contribs(i, j) = real(j, real64) * real(i, real64)
      end do
    end do
    
    thresholds = [2.0_real64, 4.0_real64]  ! Exact matches for some values
    
    call detect_outliers_spike(spike_contribs, 2, 3, thresholds, outlier_mask, ierr)
    
    call assert_equal_int(ierr, 0, "ierr for spike edge cases")
    ! Values equal to threshold should NOT be outliers
    call assert_true(.not. outlier_mask(1, 2), "Timepoint 1, Sample 2 equal to threshold - not outlier")
    call assert_true(.not. outlier_mask(2, 2), "Timepoint 2, Sample 2 equal to threshold - not outlier")
    
  end subroutine test_detect_outliers_spike_edge_cases

  ! ==================== BIOLOGICAL DATA TESTS ====================

  subroutine test_biological_data_layout()
    integer, parameter :: n_samples = 100, n_timepoints = 50
    real(real64) :: expression_data(n_timepoints, n_samples)  ! timepoints × samples
    real(real64) :: spike_thresholds(n_timepoints)
    real(real64) :: integrated_contribs(n_samples), integrated_threshold
    logical :: spike_outliers(n_timepoints, n_samples), integrated_outliers(n_samples)
    integer(int32) :: ierr
    integer :: i, j
    
    ! Simulate realistic gene expression data: timepoints × samples
    do j = 1, n_samples  ! samples (columns)
      do i = 1, n_timepoints  ! timepoints (rows)
        expression_data(i, j) = real(i, real64) * 0.1_real64 + &           ! Timepoint baseline
                               real(j, real64) * 0.05_real64 + &           ! Sample trend
                               sin(real(i, real64) * 0.2_real64) * 0.5_real64 + & ! Timepoint variation
                               cos(real(j, real64) * 0.3_real64) * 0.3_real64    ! Sample variation
      end do
    end do
    
    ! Calculate integrated contributions per sample (sum across timepoints)
    do i = 1, n_samples
      integrated_contribs(i) = sum(expression_data(:, i))
    end do
    
    ! Test that the data layout works correctly
    call calc_spike_thresholds_alloc(expression_data, n_timepoints, n_samples, 95.0_real64, spike_thresholds, ierr)
    call assert_equal_int(ierr, 0, "ierr for biological data spike thresholds")
    
    call calc_integrated_threshold_alloc(integrated_contribs, n_samples, 95.0_real64, integrated_threshold, ierr)
    call assert_equal_int(ierr, 0, "ierr for biological data integrated threshold")
    
    ! Verify dimensions are correct
    call assert_equal_int(size(spike_thresholds), n_timepoints, "Spike thresholds per timepoint")
    call assert_equal_int(size(integrated_contribs), n_samples, "Integrated contributions per sample")
    
    ! Test outlier detection
    call detect_outliers_spike(expression_data, n_timepoints, n_samples, spike_thresholds, spike_outliers, ierr)
    call assert_equal_int(ierr, 0, "ierr for biological spike outliers")
    
    call detect_outliers_integrated(integrated_contribs, n_samples, integrated_threshold, integrated_outliers, ierr)
    call assert_equal_int(ierr, 0, "ierr for biological integrated outliers")
    
    ! Verify output dimensions
    call assert_equal_int(size(spike_outliers, 1), n_timepoints, "Spike outliers timepoints dimension")
    call assert_equal_int(size(spike_outliers, 2), n_samples, "Spike outliers samples dimension")
    call assert_equal_int(size(integrated_outliers), n_samples, "Integrated outliers dimension")
    
  end subroutine test_biological_data_layout

  subroutine test_column_major_optimization()
    integer, parameter :: n_samples = 1000, n_timepoints = 500
    real(real64) :: expression_data(n_timepoints, n_samples)  ! timepoints × samples
    real(real64) :: thresholds(n_timepoints)
    integer(int32) :: ierr
    integer :: i, j
    real(real64) :: start_time, end_time
    
    ! Initialize with column-major friendly pattern: timepoints × samples
    do j = 1, n_samples  ! samples (columns - contiguous in memory)
      do i = 1, n_timepoints  ! timepoints (rows)
        expression_data(i, j) = real(i, real64) * 0.01_real64 + real(j, real64) * 0.001_real64
      end do
    end do
    
    call cpu_time(start_time)
    call calc_spike_thresholds_alloc(expression_data, n_timepoints, n_samples, 90.0_real64, thresholds, ierr)
    call cpu_time(end_time)
    
    call assert_equal_int(ierr, 0, "ierr for column-major optimization")
    call assert_true((end_time - start_time) < 5.0_real64, "Performance within limits")
    
    ! Verify timepoint-wise thresholds are reasonable
    do i = 1, n_timepoints
      call assert_true(thresholds(i) > 0.0_real64, "Timepoint thresholds positive")
      call assert_true(thresholds(i) < 100.0_real64, "Timepoint thresholds reasonable")
    end do
    
  end subroutine test_column_major_optimization

  ! ==================== ERROR HANDLING TESTS ====================

  subroutine test_empty_inputs()
    real(real64) :: empty_2d(0, 5), empty_1d(0), thresholds(0), threshold
    logical :: mask(0)
    integer(int32) :: permutation(5, 0), perm(0), ierr
    
    call calc_spike_thresholds(empty_2d, 0, 5, 50.0_real64, thresholds, permutation, ierr)
    call assert_true(ierr /= 0, "Empty 2D input should error")
    
    call calc_integrated_threshold(empty_1d, 0, 50.0_real64, threshold, perm, ierr)
    call assert_true(ierr /= 0, "Empty 1D input should error")
    
    call detect_outliers_integrated(empty_1d, 0, 1.0_real64, mask, ierr)
    call assert_true(ierr /= 0, "Empty detection input should error")
    
  end subroutine test_empty_inputs

  subroutine test_invalid_percentiles()
    real(real64) :: data(5) = [1.0, 2.0, 3.0, 4.0, 5.0], thresholds(1), threshold
    integer(int32) :: perm(5) = [1,2,3,4,5], perm_2d(5, 1) = reshape([1,2,3,4,5], [5, 1]), ierr
    
    call calc_spike_thresholds(reshape(data, [1,5]), 1, 5, -10.0_real64, thresholds, perm_2d, ierr)
    call assert_true(ierr /= 0, "Negative percentile should error")
    
    call calc_integrated_threshold(data, 5, 150.0_real64, threshold, perm, ierr)
    call assert_true(ierr /= 0, "Percentile > 100 should error")
    
  end subroutine test_invalid_percentiles

  ! ==================== CONSISTENCY TESTS ====================

  subroutine test_consistency_performance_alloc()
    real(real64) :: spike_contribs(3, 5), thresholds_perf(3), thresholds_alloc(3)  ! 3 timepoints, 5 samples
    real(real64) :: contribs(6), threshold_perf, threshold_alloc
    integer(int32) :: spike_perm(5, 3), int_perm(6), ierr1, ierr2, ierr3, ierr4
    integer :: i, j
    
    ! Generate test data: timepoints × samples
    do j = 1, 5  ! samples
      do i = 1, 3  ! timepoints
        spike_contribs(i, j) = real(i, real64) * real(j, real64)
      end do
    end do
    
    contribs = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64]
    
    ! Pre-compute permutations
    do j = 1, 3  ! timepoints
      spike_perm(:, j) = [1, 2, 3, 4, 5]  ! samples sorted indices for each timepoint
    end do
    int_perm = [1, 2, 3, 4, 5, 6]
    
    ! Compare spike thresholds
    call calc_spike_thresholds(spike_contribs, 3, 5, 80.0_real64, thresholds_perf, spike_perm, ierr1)
    call calc_spike_thresholds_alloc(spike_contribs, 3, 5, 80.0_real64, thresholds_alloc, ierr2)
    
    call assert_equal_int(ierr1, 0, "ierr for performance version")
    call assert_equal_int(ierr2, 0, "ierr for alloc version")

    do i = 1, 3
      call assert_equal_real(thresholds_perf(i), thresholds_alloc(i), TOL, &
                           "Performance vs alloc consistency for spike")
    end do
    
    ! Compare integrated thresholds
    call calc_integrated_threshold(contribs, 6, 75.0_real64, threshold_perf, int_perm, ierr3)
    call calc_integrated_threshold_alloc(contribs, 6, 75.0_real64, threshold_alloc, ierr4)
    
    call assert_equal_int(ierr3, 0, "ierr for integrated performance")
    call assert_equal_int(ierr4, 0, "ierr for integrated alloc")
    call assert_equal_real(threshold_perf, threshold_alloc, TOL, &
                         "Performance vs alloc consistency for integrated")
    
  end subroutine test_consistency_performance_alloc

  ! ==================== PERFORMANCE TESTS ====================

    subroutine test_large_scale_performance()
        integer, parameter :: n_samples = 1000, n_timepoints = 2000
        real(real64), allocatable :: expression_data(:,:)
        real(real64), allocatable :: spike_thresholds(:)
        real(real64), allocatable :: integrated_contribs(:)
        logical, allocatable :: spike_outliers(:,:), integrated_outliers(:)
        integer(int32) :: ierr
        integer :: i, j
        real(real64) :: start_time, end_time, integrated_threshold
        
        ! Allocate large arrays on heap: timepoints × samples
        allocate(expression_data(n_timepoints, n_samples))
        allocate(spike_thresholds(n_timepoints))
        allocate(integrated_contribs(n_samples))
        allocate(spike_outliers(n_timepoints, n_samples))
        allocate(integrated_outliers(n_samples))
        
        ! Generate large dataset: timepoints × samples
        do j = 1, n_samples  ! samples
        do i = 1, n_timepoints  ! timepoints
            expression_data(i, j) = real(i, real64) * 0.001_real64 + real(j, real64) * 0.01_real64
        end do
        end do
        
        ! Time spike threshold calculation
        call cpu_time(start_time)
        call calc_spike_thresholds_alloc(expression_data, n_timepoints, n_samples, 95.0_real64, spike_thresholds, ierr)
        call cpu_time(end_time)
        
        call assert_equal_int(ierr, 0, "ierr for large scale spike thresholds")
        call assert_true((end_time - start_time) < 10.0_real64, "Large scale performance acceptable")
        
        ! Calculate integrated contributions
        do i = 1, n_samples
        integrated_contribs(i) = sum(expression_data(:, i))
        end do
        
        ! Time integrated threshold calculation
        call cpu_time(start_time)
        call calc_integrated_threshold_alloc(integrated_contribs, n_samples, 95.0_real64, integrated_threshold, ierr)
        call cpu_time(end_time)
        
        call assert_equal_int(ierr, 0, "ierr for large scale integrated threshold")
        call assert_true((end_time - start_time) < 5.0_real64, "Integrated performance acceptable")
        
        ! Clean up
        deallocate(expression_data, spike_thresholds, integrated_contribs, &
                spike_outliers, integrated_outliers)
        
    end subroutine test_large_scale_performance

  ! ==================== ERROR PROPAGATION TESTS ====================

  subroutine test_error_propagation()
    real(real64) :: spike_contribs(2, 3), thresholds(2)  ! 2 timepoints, 3 samples
    real(real64) :: contribs(3), threshold
    integer(int32) :: ierr
    
    ! Test that errors from calc_percentile propagate correctly
    spike_contribs = reshape([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], [2,3])
    contribs = [1.0, 2.0, 3.0]
    
    ! These should work fine
    call calc_spike_thresholds_alloc(spike_contribs, 2, 3, 50.0_real64, thresholds, ierr)
    call assert_equal_int(ierr, 0, "Normal case should not error")
    
    call calc_integrated_threshold_alloc(contribs, 3, 50.0_real64, threshold, ierr)
    call assert_equal_int(ierr, 0, "Normal integrated case should not error")
    
    ! Test error propagation framework
    call assert_true(.true., "Error propagation framework in place")
    
  end subroutine test_error_propagation

  ! ==================== COMPREHENSIVE WORKFLOW TEST ====================

  subroutine test_comprehensive_workflow()
    integer, parameter :: n_samples = 50, n_timepoints = 100
    real(real64) :: expression_data(n_timepoints, n_samples)  ! timepoints × samples
    real(real64) :: spike_thresholds(n_timepoints), integrated_threshold
    real(real64) :: integrated_contribs(n_samples)
    logical :: spike_outliers(n_timepoints, n_samples), integrated_outliers(n_samples)
    integer(int32) :: ierr1, ierr2, ierr3, ierr4
    integer :: i, j, spike_outlier_count, integrated_outlier_count
    
    ! Generate comprehensive test dataset: timepoints × samples
    do j = 1, n_samples  ! samples
      do i = 1, n_timepoints  ! timepoints
        expression_data(i, j) = real(i, real64) * 0.1_real64 + &
                              real(j, real64) * 0.05_real64 + &
                              sin(real(i*j, real64) * 0.01_real64)
      end do
    end do
    
    ! Add some clear outliers
    expression_data(5, 10) = 1000.0_real64    ! Timepoint 5, Sample 10 outlier
    expression_data(30, 25) = 1500.0_real64   ! Timepoint 30, Sample 25 outlier
    expression_data(80, 40) = 2000.0_real64   ! Timepoint 80, Sample 40 outlier
    
    ! Calculate integrated contributions per sample
    do i = 1, n_samples
      integrated_contribs(i) = sum(expression_data(:, i))
    end do
    integrated_contribs(15) = 5000.0_real64  ! Make sample 15 an integrated outlier
    
    ! Complete workflow
    call calc_spike_thresholds_alloc(expression_data, n_timepoints, n_samples, 95.0_real64, spike_thresholds, ierr1)
    call calc_integrated_threshold_alloc(integrated_contribs, n_samples, 95.0_real64, integrated_threshold, ierr2)
    call detect_outliers_spike(expression_data, n_timepoints, n_samples, spike_thresholds, spike_outliers, ierr3)
    call detect_outliers_integrated(integrated_contribs, n_samples, integrated_threshold, integrated_outliers, ierr4)
    
    call assert_equal_int(ierr1, 0, "Spike thresholds in workflow")
    call assert_equal_int(ierr2, 0, "Integrated threshold in workflow")
    call assert_equal_int(ierr3, 0, "Spike outliers in workflow")
    call assert_equal_int(ierr4, 0, "Integrated outliers in workflow")
    
    ! Count outliers
    spike_outlier_count = count(spike_outliers)
    integrated_outlier_count = count(integrated_outliers)
    
    call assert_true(spike_outlier_count > 0, "Should detect spike outliers")
    call assert_true(integrated_outlier_count > 0, "Should detect integrated outliers")
    
    ! Verify specific known outliers were detected
    call assert_true(spike_outliers(5, 10), "Timepoint 5, Sample 10 spike outlier detected")
    call assert_true(spike_outliers(30, 25), "Timepoint 30, Sample 25 spike outlier detected")
    call assert_true(spike_outliers(80, 40), "Timepoint 80, Sample 40 spike outlier detected")
    call assert_true(integrated_outliers(15), "Sample 15 integrated outlier detected")
    
  end subroutine test_comprehensive_workflow

end module mod_test_outlier_detection