!> Module for trajectory explanatory contribution analysis (TCA)
!! Provides functionality for detecting outliers in trajectory and spike contributions
module tox_trajectory_contribution_analysis_outlier_detection
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_ALLOC_FAIL, set_ok, set_err, is_err, validate_dimension_size
    use f42_utils, only: calc_percentile, calc_percentile_alloc, sort_real
    
    implicit none
    
    private
    public :: calc_spike_thresholds, calc_integrated_threshold, &
              detect_outliers_integrated, detect_outliers_spike, &
              calc_spike_thresholds_alloc, calc_integrated_threshold_alloc
    
contains

    !> Calculate empirical thresholds for spike contributions (using pre-sorted permutations)
    pure subroutine calc_spike_thresholds(spike_contribs, n_timepoints, n_samples, percentile_val, thresholds, permutation, ierr)
        integer(int32), intent(in) :: n_timepoints
        !! Number of timepoints
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: spike_contribs(n_timepoints, n_samples)
        !! 2D array of spike contributions
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(n_timepoints)
        !! 1D array of thresholds for each timepoint
        integer(int32), intent(in) :: permutation(n_samples, n_timepoints)
        !! Pre-computed permutation indices
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: j
        
        ! Initialize error
        call set_ok(ierr)
        
        
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if(is_err(ierr)) return
        
        ! Calculate threshold for each timepoint using pre-sorted permutations
        do j = 1, n_timepoints
            call calc_percentile(spike_contribs(j, :), permutation(:, j), percentile_val, thresholds(j), ierr)
            if (is_err(ierr)) then
                return
            end if
        end do
        
    end subroutine calc_spike_thresholds

    !> Calculate empirical thresholds for spike contributions (with internal allocations and sorting)
    subroutine calc_spike_thresholds_alloc(spike_contribs, n_timepoints, n_samples, percentile_val, thresholds, ierr)
        integer(int32), intent(in) :: n_timepoints
        !! Number of timepoints
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: spike_contribs(n_timepoints, n_samples)
        !! 2D array of spike contributions
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: thresholds(n_timepoints)
        !! 1D array of thresholds
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) :: j, i, k
        integer(int32), allocatable :: permutation(:, :), stack_left(:), stack_right(:), temp_perm(:)
        real(real64), allocatable :: sample_data(:)
        
        ! Initialize error
        call set_ok(ierr)
        
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if(is_err(ierr)) return
        
        allocate(permutation(n_samples, n_timepoints), stat=ierr)
        allocate(temp_perm(n_timepoints), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            return
        end if

        do i = 1, n_timepoints
            permutation(:, i) = [(k, k = 1, n_samples)]
        end do
        
        allocate(stack_left(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation)
            return
        end if
        
        allocate(stack_right(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation, stack_left)
            return
        end if
        
        allocate(sample_data(n_samples), stat=ierr)
        if (is_err(ierr)) then
            call set_err(ierr, ERR_ALLOC_FAIL)
            deallocate(permutation, stack_left, stack_right)
            return
        end if
        
        ! Pre-compute permutations
        do j = 1, n_timepoints
            sample_data = spike_contribs(j, :)
            call sort_real(sample_data, permutation(:, j), stack_left, stack_right)
        end do
        
        ! Now call the main subroutine with pre-computed permutations
        call calc_spike_thresholds(spike_contribs, n_timepoints, n_samples, percentile_val, thresholds, permutation, ierr)
        
    end subroutine calc_spike_thresholds_alloc

    !> Calculate empirical threshold for integrated (trajectory-level) contributions
    !!
    !! Determines which SAMPLES are significant overall across all genes.
    !! Uses pre-computed permutation array to avoid internal sorting.
    pure subroutine calc_integrated_threshold(contributions, n_samples, percentile_val, threshold, permutation, ierr)
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: contributions(n_samples)
        !! 1D array of integrated contributions
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(in) :: permutation(n_samples)
        !! Pre-computed permutation indices for sorted contributions
        integer(int32), intent(out) :: ierr
        !! Error code
        
        ! Input validation is handled by calc_percentile
        call calc_percentile(contributions, permutation, percentile_val, threshold, ierr)
        
    end subroutine calc_integrated_threshold

    !> Calculate empirical threshold for integrated contributions (with internal allocation and sorting)
    !!
    !! Convenience version that handles workspace allocation and sorting internally.
    !! For performance-critical code, use calc_integrated_threshold with pre-computed permutation.
    subroutine calc_integrated_threshold_alloc(contributions, n_samples, percentile_val, threshold, ierr)
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: contributions(n_samples)
        !! 1D array of integrated contributions
        real(real64), intent(in) :: percentile_val
        !! Percentile value for threshold (0.0-100.0)
        real(real64), intent(out) :: threshold
        !! Scalar threshold value
        integer(int32), intent(out) :: ierr
        !! Error code
        
        ! Use the existing calc_percentile_alloc from f42_utils
        call calc_percentile_alloc(contributions, percentile_val, threshold, ierr)
        
    end subroutine calc_integrated_threshold_alloc

    !> Detect outliers in integrated (trajectory-level) contributions
    !!
    !! Identifies SAMPLES where the integrated contribution exceeds the empirical threshold.
    !! Used to find trajectories with significant overall explanatory contributions.
    pure subroutine detect_outliers_integrated(contributions, n_samples, threshold, outlier_mask, ierr)
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: contributions(n_samples)
        !! Array of integrated contributions
        real(real64), intent(in) :: threshold
        !! Scalar threshold value for outlier detection
        logical, intent(out) :: outlier_mask(n_samples)
        !! 1D logical array indicating outliers
        integer(int32), intent(out) :: ierr
        !! Error code
        integer(int32) :: i
        
        ! Initialize error
        call set_ok(ierr)
        
        ! Input validation        
        call validate_dimension_size(n_samples, ierr)
        if(is_err(ierr)) return
        
        ! Detect outlier SAMPLES
        do i = 1, n_samples
            outlier_mask(i) = contributions(i) > threshold
        end do
        
    end subroutine detect_outliers_integrated

    !> Detect outliers in spike contributions
    pure subroutine detect_outliers_spike(spike_contribs, n_timepoints, n_samples, thresholds, outlier_mask, ierr)
        integer(int32), intent(in) :: n_timepoints
        !! Number of timepoints
        integer(int32), intent(in) :: n_samples
        !! Number of samples
        real(real64), intent(in) :: spike_contribs(n_timepoints, n_samples)
        !! 2D array of spike contributions 
        real(real64), intent(in) :: thresholds(n_timepoints)
        !! Array of thresholds
        logical, intent(out) :: outlier_mask(n_timepoints, n_samples)
        !! 2D logical array indicating outliers
        integer(int32), intent(out) :: ierr
        !! Error code
        
        integer(int32) ::i, j
        
        ! Initialize error
        call set_ok(ierr)
        
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_samples, ierr)
        if(is_err(ierr)) return
        
        do j = 1, n_samples
            do i = 1, n_timepoints
                outlier_mask(i, j) = spike_contribs(i, j) > thresholds(i)
            end do
        end do
        
    end subroutine detect_outliers_spike
end module tox_trajectory_contribution_analysis_outlier_detection

!> C wrapper for calc_spike_thresholds (with pre-computed 2D permutations)
subroutine calc_spike_thresholds_expert_C(spike_contribs, n_timepoints, n_samples, &
                                    percentile_val, thresholds, permutation, ierr) &
                                    bind(C, name="calc_spike_thresholds_expert_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_contribution_analysis_outlier_detection, only: calc_spike_thresholds
    
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    integer(c_int), intent(in), value :: n_timepoints
    !! Number of genes
    real(c_double), intent(in) :: spike_contribs(n_timepoints, n_samples)
    !! 2D array of spike contributions
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: thresholds(n_timepoints)
    !! 1D array of thresholds
    integer(c_int), intent(in) :: permutation(n_samples, n_timepoints) 
    !! Pre-computed permutation indices
    integer(c_int), intent(out) :: ierr
    !! Error code

    call calc_spike_thresholds(spike_contribs, n_timepoints, n_samples, percentile_val, &
                                thresholds, permutation, ierr)
    
end subroutine calc_spike_thresholds_expert_C

!> C wrapper for calc_spike_thresholds_alloc (with internal allocations)
!! Data layout: spike_contribs(n_samples, n_genes) - rows=samples, columns=genes
subroutine calc_spike_thresholds_C(spike_contribs, n_timepoints, n_samples, &
                                        percentile_val, thresholds, ierr) &
                                        bind(C, name="calc_spike_thresholds_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_contribution_analysis_outlier_detection, only: calc_spike_thresholds_alloc
    integer(c_int), intent(in), value :: n_timepoints
    !! number of timepoints
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples    
    real(c_double), intent(in) :: spike_contribs(n_timepoints, n_samples)
    !! Array of spike contributions
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: thresholds(n_timepoints)
    !! 1D array of thresholds
    integer(c_int), intent(out) :: ierr
    !! Error code

    ! Call Fortran subroutine
    call calc_spike_thresholds_alloc(spike_contribs, n_timepoints, n_samples, percentile_val, &
                                    thresholds, ierr)

end subroutine calc_spike_thresholds_C

!> C wrapper for calc_integrated_threshold (with pre-computed permutation)
!! Input: contributions(n_samples) - one integrated contribution per sample
subroutine calc_integrated_threshold_expert_C(contributions, n_samples, percentile_val, &
                                        threshold, permutation, ierr) &
                                        bind(C, name="calc_integrated_threshold_expert_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_contribution_analysis_outlier_detection, only: calc_integrated_threshold
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in) :: contributions(n_samples)
    !! 1D array of integrated contributions [n_samples]
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: threshold
    !! Scalar threshold value
    integer(c_int), intent(in) :: permutation(n_samples)
    !! Pre-computed permutation indices for sorted contributions
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Call Fortran subroutine
    call calc_integrated_threshold(contributions, n_samples, percentile_val, &
                                threshold, permutation, ierr)
    
end subroutine calc_integrated_threshold_expert_C

!> C wrapper for calc_integrated_threshold_alloc (with internal allocations)
!! Input: contributions(n_samples) - one integrated contribution per sample
subroutine calc_integrated_threshold_C(contributions, n_samples, percentile_val, &
                                            threshold, ierr) &
                                            bind(C, name="calc_integrated_threshold_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_contribution_analysis_outlier_detection, only: calc_integrated_threshold_alloc
    real(c_double), intent(in) :: contributions(n_samples)
    !! 1D array of integrated contributions
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in), value :: percentile_val
    !! Percentile value for threshold (0.0-100.0)
    real(c_double), intent(out) :: threshold
    !! Scalar threshold value
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    ! Call Fortran subroutine
    call calc_integrated_threshold_alloc(contributions, n_samples, percentile_val, &
                                        threshold, ierr)
    
end subroutine calc_integrated_threshold_C

!> C wrapper for detect_outliers_integrated
!! Identifies outlier samples based on integrated contributions
subroutine detect_outliers_integrated_expert_C(contributions, n_samples, threshold, &
                                        outlier_mask, ierr) &
                                        bind(C, name="detect_outliers_integrated_expert_C")
    use tox_conversions, only: logical_as_c_int
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, ERR_ALLOC_FAIL, is_ok, is_err, ERR_EMPTY_INPUT
    use tox_trajectory_contribution_analysis_outlier_detection, only: detect_outliers_integrated
    real(c_double), intent(in) :: contributions(n_samples)
    !! Array of integrated contributions
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    real(c_double), intent(in), value :: threshold
    !! Scalar threshold value for outlier detection
    integer(c_int), intent(out) :: outlier_mask(n_samples)  ! c_int for logical
    !! 1D array indicating outlier samples 
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    logical, allocatable :: mask_1d(:)
    integer :: i
    
    ! Initialize error
    call set_ok(ierr)
    
    ! Check for valid dimensions
    if (n_samples <= 0) then
        call set_err(ierr, ERR_EMPTY_INPUT)
        return
    end if

    allocate(mask_1d(n_samples), stat=ierr)
    if(is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if
    
    ! Call Fortran subroutine
    call detect_outliers_integrated(contributions, n_samples, threshold, &
                                    mask_1d, ierr)
    
    if (is_ok(ierr)) then
        call logical_as_c_int(mask_1d, outlier_mask)
    end if
    
    deallocate(mask_1d)
    
end subroutine detect_outliers_integrated_expert_C

!> C wrapper for detect_outliers_spike
subroutine detect_outliers_spike_expert_C(spike_contribs, n_timepoints, n_samples, thresholds, &
                                    outlier_mask, ierr) &
                                    bind(C, name="detect_outliers_spike_expert_C")
    use tox_conversions, only: logical_as_c_int
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_errors, only: set_ok, set_err, is_ok, is_err, ERR_ALLOC_FAIL
    use tox_trajectory_contribution_analysis_outlier_detection, only: detect_outliers_spike
    integer(c_int), intent(in), value :: n_samples
    !! Number of samples
    integer(c_int), intent(in), value :: n_timepoints
    !! Number of genes
    real(c_double), intent(in) :: spike_contribs(n_timepoints, n_samples)
    !! 2D array of spike contributions
    real(c_double), intent(in) :: thresholds(n_timepoints)
    !! Array of thresholds for each gene 
    integer(c_int), intent(out) :: outlier_mask(n_timepoints, n_samples)  ! c_int for logical
    !! 2D array indicating outlier 
    integer(c_int), intent(out) :: ierr
    !! Error code
    
    logical, allocatable :: mask_2d(:,:)
    integer :: i, j
    
    ! Initialize error
    call set_ok(ierr)

    allocate(mask_2d(n_timepoints, n_samples), stat=ierr)
    if(is_err(ierr)) then
        call set_err(ierr, ERR_ALLOC_FAIL)
        return
    end if
    
    ! Call Fortran subroutine
    call detect_outliers_spike(spike_contribs, n_timepoints, n_samples, thresholds, mask_2d, ierr)
    
    if(is_ok(ierr)) then
        ! Convert Fortran logical to C int (0=false, 1=true)
        call logical_as_c_int(mask_2d, outlier_mask)
    end if
    
    deallocate(mask_2d)
    
end subroutine detect_outliers_spike_expert_C