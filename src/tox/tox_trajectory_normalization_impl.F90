#include <src/macros.h>

!> Implementations for min-max normalization of factor trajectories over time.
!| Each factor's time series is independently rescaled to `[0,1]` per sample/entity. The
!| generator turns these into the validating (and, for the all-trajectories implementation, the
!| allocating) wrappers in module `tox_trajectory_normalization`.
module tox_trajectory_normalization_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: set_ok, set_err, ERR_DIVISION_BY_ZERO
    use f42_utils, only: is_close
    M_IMPLICIT_NONE

    private
    public :: normalize_variable_timeseries_impl, &
              normalize_single_trajectory_impl, &
              normalize_all_trajectories_impl

contains

    !> summary: Normalize a single variable across time using min-max scaling
    !| AUTHOR_AARON_SCHROEDER
    pure subroutine normalize_variable_timeseries_impl(v, v_norm, n_points, status)
        integer(int32), intent(in) :: n_points
            !! Vector length (number of time points)
        real(real64), intent(in) :: v(n_points)
            !! Original time series
        real(real64), intent(out) :: v_norm(n_points)
            !! Normalized time series
        integer(int32), intent(out) :: status
            !! Status code for specific warnings

        real(real64) :: min_val, max_val, denominator
        integer(int32) :: i_point

        call set_ok(status)

        min_val = minval(v)
        max_val = maxval(v)
        denominator = max_val - min_val

        ! Check for division by zero (min approximately equal to max)
        if (is_close(denominator, 0.0_real64)) then
            v_norm = 0.0_real64
            call set_err(status, ERR_DIVISION_BY_ZERO)
            return
        end if

        do concurrent (i_point = 1:n_points) shared(v_norm, v, min_val, denominator)
            v_norm(i_point) = (v(i_point) - min_val)/denominator
        end do
    end subroutine normalize_variable_timeseries_impl

    !> summary: Normalize all factors in a single trajectory independently across time
    !| AUTHOR_AARON_SCHROEDER
    !| Input: `trajectory(n_timepoints, n_factors)` for ONE sample/entity
    pure subroutine normalize_single_trajectory_impl(trajectory, trajectory_norm, n_factors, n_timepoints, status)
        integer(int32), intent(in) :: n_factors
            !! Number of factors/variables
        integer(int32), intent(in) :: n_timepoints
            !! Number of time points
        real(real64), intent(in) :: trajectory(n_timepoints, n_factors)
            !! Original trajectory for one sample
        real(real64), intent(out) :: trajectory_norm(n_timepoints, n_factors)
            !! Normalized trajectory for one sample
        integer(int32), dimension(n_factors), intent(out) :: status
            !! Status code for specific warnings, one per factor

        integer(int32) :: i_factor

        call set_ok(status)

        ! Normalize each factor independently across time
        do i_factor = 1, n_factors
            call normalize_variable_timeseries_impl( &
                trajectory(:, i_factor), &
                trajectory_norm(:, i_factor), &
                n_timepoints, status(i_factor))
        end do
    end subroutine normalize_single_trajectory_impl

    !> summary: Normalize all trajectories across multiple entities
    !| AUTHOR_AARON_SCHROEDER
    !| Input: `trajectories(n_factors, n_samples, n_timepoints)`; each factor is normalized
    !| independently across time for each sample.
    pure subroutine normalize_all_trajectories_impl(trajectories, trajectories_norm, &
                                               n_factors, n_samples, n_timepoints, &
                                               tmp_series, tmp_series_norm, status)
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
        real(real64), intent(out) :: tmp_series(n_timepoints)
            !! Work array: one factor's time series, gathered contiguously
        real(real64), intent(out) :: tmp_series_norm(n_timepoints)
            !! Work array: the normalized time series
        integer(int32), dimension(n_factors, n_samples), intent(out) :: status
            !! Status code for specific warnings, one per factor per sample

        integer(int32) :: i_sample, i_factor, i_timepoint

        call set_ok(status)

        ! Normalize each sample/entity independently
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                do concurrent (i_timepoint = 1:n_timepoints) shared(tmp_series, trajectories, i_sample, i_factor)
                    tmp_series(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
                end do

                call normalize_variable_timeseries_impl(tmp_series, tmp_series_norm, n_timepoints, status(i_factor, i_sample))

                trajectories_norm(i_factor, i_sample, :) = tmp_series_norm
            end do
        end do
    end subroutine normalize_all_trajectories_impl

end module tox_trajectory_normalization_impl
