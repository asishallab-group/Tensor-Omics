#include <src/macros.h>

!> Module for min-max normalization of factor trajectories over time.
!| Each factor's time series is independently rescaled to `[0,1]` per sample/entity.
module tox_trajectory_normalization
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_NAN_INF, set_ok, set_err, is_err, validate_dimension_size, validate_all_in_range_real, ERR_DIVISION_BY_ZERO, ERR_ALLOC_FAIL
    use f42_utils, only: is_close
    M_IMPLICIT_NONE

    private
    public :: normalize_variable_timeseries, &
              normalize_single_trajectory, &
              normalize_all_trajectories_alloc

contains

    !> M_EXPORT_C
    !| summary: Normalize a single variable across time using min-max scaling
    !| AUTHOR_AARON_SCHROEDER
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
        integer(int32) :: i_point

        ! Initialize
        call set_ok(ierr)
        call set_ok(status)

        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(v, n_points, ierr, arg_pos=1_int32)
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
        do concurrent (i_point = 1:n_points) shared(v_norm, v, min_val, denominator)
            v_norm(i_point) = (v(i_point) - min_val)/denominator
        end do

    end subroutine normalize_variable_timeseries

    !> M_EXPORT_C
    !| summary: Normalize all factors in a single trajectory independently across time
    !| AUTHOR_AARON_SCHROEDER
    !| Input: `trajectory(n_factors, n_timepoints)` for ONE sample/entity
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
        integer(int32), dimension(n_factors), intent(out) :: status
            !! Status code for specific warnings

        integer(int32) :: i_factor

        ! Initialize
        call set_ok(ierr)
        call set_ok(status)

        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        if (is_err(ierr)) return

        ! Normalize each factor independently across time
        do i_factor = 1, n_factors
            call normalize_variable_timeseries( &
                trajectory(:, i_factor), &           ! Time series for this factor
                trajectory_norm(:, i_factor), &      ! Normalized time series
                n_timepoints, ierr, status(i_factor))

            if (is_err(ierr)) return
        end do

    end subroutine normalize_single_trajectory

    !> M_EXPORT_C
    !| summary: Normalize all trajectories across multiple entities
    !| AUTHOR_AARON_SCHROEDER
    !| Input: `trajectories(n_factors, n_samples, n_timepoints)`
    !| Normalizes each factor independently across time for each sample
    pure subroutine normalize_all_trajectories_alloc(trajectories, trajectories_norm, &
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
        integer(int32), dimension(n_factors, n_samples), intent(out) :: status
            !! Status code for specific warnings

        integer(int32) :: i_sample, i_factor, i_timepoint

        real(real64), dimension(:), allocatable :: tmp_series, tmp_series_norm

        call set_ok(ierr)
        call set_ok(status)

        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_series(n_timepoints))
        M_ALLOCATE(tmp_series_norm(n_timepoints))

        ! Normalize each sample/entity independently
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                do concurrent (i_timepoint = 1:n_timepoints) shared(tmp_series, trajectories, i_sample, i_factor)
                    tmp_series(i_timepoint) = trajectories(i_factor, i_sample, i_timepoint)
                end do

                call normalize_variable_timeseries(tmp_series, tmp_series_norm, n_timepoints, ierr, status(i_factor, i_sample))

                if (is_err(ierr)) return

                trajectories_norm(i_factor, i_sample, :) = tmp_series_norm
            end do
        end do
    end subroutine normalize_all_trajectories_alloc

end module tox_trajectory_normalization
