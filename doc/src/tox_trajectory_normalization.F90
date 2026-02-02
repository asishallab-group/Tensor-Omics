#include "macros.h"
module tox_trajectory_normalization
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_NAN_INF, set_ok, set_err, is_err, validate_dimension_size, validate_all_in_range_real, ERR_DIVISION_BY_ZERO
    use safeguard
    use f42_utils, only: is_close
    implicit none
    
    private
    public :: normalize_variable_timeseries, &
              normalize_single_trajectory, &
              normalize_all_trajectories
    
contains
    
    !> Normalize a single variable across time using min-max scaling
    pure subroutine normalize_variable_timeseries(v, v_norm, n_points, ierr, status)
        integer(int32), intent(in) :: n_points
        !! Vector length (number of time points)
        real(real64), intent(in) :: v(n_points)
        !! Original time series
        real(real64), intent(out) :: v_norm(n_points)
        !! Normalized time series        
        integer(int32), intent(out) :: ierr
        !! Error code
        integer(int32), intent(out) :: status
        !! Status code for specific warnings
        
        real(real64) :: min_val, max_val, denominator
        integer(int32) :: i
        
        ! Initialize
        call set_ok(ierr)
        call set_ok(status)

        call validate_all_in_range_real(v, n_points, ierr)
        call validate_dimension_size(n_points, ierr)
        if (is_err(ierr)) return

        min_val = minval(v)
        max_val = maxval(v)
        
        ! Calculate denominator
        denominator = max_val - min_val
        
        ! Check for division by zero (min approximately equal to max)
        if (is_close(denominator, 0.0_real64)) then
            v_norm = 0.0_real64
            call set_err(status, ERR_DIVISION_BY_ZERO)
            return
        end if
        
        ! Apply min-max normalization
        do i = 1, n_points
            v_norm(i) = (v(i) - min_val) / denominator
        end do
        
    end subroutine normalize_variable_timeseries
    
    !> Normalize all factors in a single trajectory independently across time
    !! Input: trajectory(n_factors, n_timepoints) for ONE sample/entity
    pure subroutine normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, ierr, status)
        integer(int32), intent(in) :: n_factors
        !! Number of factors/variables
        integer(int32), intent(in) :: n_timepoints
        !! Number of time points
        real(real64), intent(in) :: trajectory(n_timepoints, n_factors)
        !! Original trajectory for one sample
        real(real64), intent(out) :: trajectory_norm(n_timepoints, n_factors)
        !! Normalized trajectory for one sample
        integer(int32), intent(out) :: ierr
        !! Error code
        integer(int32), intent(out) :: status
        !! Status code for specific warnings
        
        integer(int32) :: i_factor, tmp_status

        ! Initialize
        call set_ok(ierr)
        call set_ok(status)
        call set_ok(tmp_status)
        
        call validate_dimension_size(n_factors, ierr)
        if (is_err(ierr)) return
        
        ! Normalize each factor independently across time
        do i_factor = 1, n_factors
            call normalize_variable_timeseries( &
                trajectory(:, i_factor), &           ! Time series for this factor
                trajectory_norm(:, i_factor), &      ! Normalized time series
                n_timepoints, ierr, tmp_status)
            
            if (is_err(tmp_status)) status = tmp_status
            if (is_err(ierr)) return
        end do
        
    end subroutine normalize_single_trajectory
    
    !> Normalize all trajectories across multiple entities
    !! Input: trajectories(n_factors, n_samples, n_timepoints)
    !! Normalizes each factor independently across time for each sample
    pure subroutine normalize_all_trajectories(trajectories, trajectories_norm, &
                                          n_factors, n_samples, n_timepoints, ierr, status)
        integer(int32), intent(in) :: n_factors
        !! Number of factors
        integer(int32), intent(in) :: n_samples
        !! Number of samples/entities
        integer(int32), intent(in) :: n_timepoints
        !! Number of time points
        real(real64), intent(in) :: trajectories(n_factors, n_samples, n_timepoints)
        !! Original trajectories
        real(real64), intent(out) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        !! Normalized trajectories
        integer(int32), intent(out) :: ierr
        !! Error code
        integer(int32), intent(out) :: status
        !! Status code for specific warnings

        integer(int32) :: i_sample, i_factor, i_timepoint, tmp_status

        real(real64), dimension(n_timepoints) :: temp_series, temp_series_norm
        
        call set_ok(ierr)
        call set_ok(status)
        call set_ok(tmp_status)

        call validate_dimension_size(n_samples, ierr)
        if (is_err(ierr)) return
        
        ! Normalize each sample/entity independently
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                do i_timepoint = 1, n_timepoints
                    temp_series(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
                end do

                call normalize_variable_timeseries(temp_series, temp_series_norm, n_timepoints, ierr, tmp_status)

                if (is_err(tmp_status)) status = tmp_status
                if (is_err(ierr)) return
                
                trajectories_norm(i_factor, i_sample, :) = temp_series_norm
            end do
        end do
    end subroutine normalize_all_trajectories

end module tox_trajectory_normalization

pure subroutine normalize_variable_timeseries_C(v, v_norm, n_points, ierr, status) bind(C, name="normalize_variable_timeseries_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_normalization, only: normalize_variable_timeseries
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_points
    !! Vector length
    real(c_double), dimension(n_points), intent(in), target :: v
    !! Original time series
    real(c_double), dimension(n_points), intent(out), target :: v_norm
    !! Normalized time series        
    integer(c_int), intent(out), target :: ierr
    !! Error code
    integer(c_int), intent(out), target :: status
    !! Status code for specific warnings

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(v)
    M_CHECK_NON_NULL(v_norm)
    M_CHECK_NON_NULL(status)

    call normalize_variable_timeseries(v, v_norm, n_points, ierr, status)

end subroutine normalize_variable_timeseries_C

pure subroutine normalize_single_trajectory_C(trajectory, trajectory_norm, n_factors, n_timepoints, ierr, status) bind(C, name="normalize_single_trajectory_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_normalization, only: normalize_single_trajectory
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    !! Number of factors
    integer(c_int), intent(in), target :: n_timepoints
    !! Number of time points
    real(c_double), dimension(n_timepoints, n_factors), intent(in), target :: trajectory
    !! Original trajectory for one sample
    real(c_double), dimension(n_timepoints, n_factors), intent(out), target :: trajectory_norm
    !! Normalized trajectory for one sample
    integer(c_int), intent(out), target :: ierr
    !! Error code
    integer(c_int), intent(out), target :: status
    !! Status code for specific warnings

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(trajectory)
    M_CHECK_NON_NULL(trajectory_norm)
    M_CHECK_NON_NULL(status)

    call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, ierr, status)

end subroutine normalize_single_trajectory_C

pure subroutine normalize_all_trajectories_C(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr, status) bind(C, name="normalize_all_trajectories_C")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_trajectory_normalization, only: normalize_all_trajectories
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_factors
    !! Number of factors
    integer(c_int), intent(in), target :: n_samples
    !! Number of samples
    integer(c_int), intent(in), target :: n_timepoints
    !! Number of time points
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
    !! Original trajectories
    real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(out), target :: trajectories_norm
    !! Normalized trajectories
    integer(c_int), intent(out), target :: ierr
    !! Error code
    integer(c_int), intent(out), target :: status
    !! Status code for specific warnings

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_factors)
    M_CHECK_NON_NULL(n_timepoints)
    M_CHECK_NON_NULL(n_samples)
    M_CHECK_NON_NULL(trajectories)
    M_CHECK_NON_NULL(trajectories_norm)
    M_CHECK_NON_NULL(status)

    call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, ierr, status)

end subroutine normalize_all_trajectories_C