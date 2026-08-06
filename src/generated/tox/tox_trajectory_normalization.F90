#include <src/macros.h>

!> summary: Wrappers for [[tox_trajectory_normalization_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_trajectory_normalization
    use tox_trajectory_normalization_kernel, only: normalize_all_trajectories_kernel, normalize_single_trajectory_kernel, normalize_variable_timeseries_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: normalize_variable_timeseries
    public :: normalize_single_trajectory
    public :: normalize_all_trajectories
    public :: normalize_all_trajectories_alloc

contains

    !> summary: Validates its inputs, then calls [[tox_trajectory_normalization_kernel(module):normalize_variable_timeseries_kernel]].
    subroutine normalize_variable_timeseries(&
            v,&
            v_norm,&
            n_points,&
            status,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! Vector length (number of time points)
        real(real64), dimension(n_points), intent(in) :: v
            !! Original time series
        real(real64), dimension(n_points), intent(out) :: v_norm
            !! Normalized time series
        integer(int32), intent(out) :: status
            !! Status code for specific warnings
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(v, n_points, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call normalize_variable_timeseries_kernel(&
            v = v,&
            v_norm = v_norm,&
            n_points = n_points,&
            status = status&
        )
    end subroutine normalize_variable_timeseries

    !> summary: Validates its inputs, then calls [[tox_trajectory_normalization_kernel(module):normalize_single_trajectory_kernel]].
    subroutine normalize_single_trajectory(&
            trajectory,&
            trajectory_norm,&
            n_factors,&
            n_timepoints,&
            status,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! Number of factors/variables
        integer(int32), intent(in) :: n_timepoints
            !! Number of time points
        real(real64), dimension(n_timepoints, n_factors), intent(in) :: trajectory
            !! Original trajectory for one sample
        real(real64), dimension(n_timepoints, n_factors), intent(out) :: trajectory_norm
            !! Normalized trajectory for one sample
        integer(int32), dimension(n_factors), intent(out) :: status
            !! Status code for specific warnings, one per factor
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(trajectory, n_timepoints * n_factors, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call normalize_single_trajectory_kernel(&
            trajectory = trajectory,&
            trajectory_norm = trajectory_norm,&
            n_factors = n_factors,&
            n_timepoints = n_timepoints,&
            status = status&
        )
    end subroutine normalize_single_trajectory

    !> summary: Validates its inputs, then calls [[tox_trajectory_normalization_kernel(module):normalize_all_trajectories_kernel]].
    !| independently across time for each sample.
    subroutine normalize_all_trajectories(&
            trajectories,&
            trajectories_norm,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            tmp_series,&
            tmp_series_norm,&
            status,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! Number of factors
        integer(int32), intent(in) :: n_samples
            !! Number of samples/entities
        integer(int32), intent(in) :: n_timepoints
            !! Number of time points
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Original trajectories
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(out) :: trajectories_norm
            !! Normalized trajectories
        real(real64), dimension(n_timepoints), intent(out) :: tmp_series
            !! Work array: one factor's time series, gathered contiguously
        real(real64), dimension(n_timepoints), intent(out) :: tmp_series_norm
            !! Work array: the normalized time series
        integer(int32), dimension(n_factors, n_samples), intent(out) :: status
            !! Status code for specific warnings, one per factor per sample
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call normalize_all_trajectories_kernel(&
            trajectories = trajectories,&
            trajectories_norm = trajectories_norm,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            tmp_series = tmp_series,&
            tmp_series_norm = tmp_series_norm,&
            status = status&
        )
    end subroutine normalize_all_trajectories

    !> summary: Allocates its work arrays, then calls [[tox_trajectory_normalization_kernel(module):normalize_all_trajectories_kernel]].
    !| independently across time for each sample.
    subroutine normalize_all_trajectories_alloc(&
            trajectories,&
            trajectories_norm,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            status,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! Number of factors
        integer(int32), intent(in) :: n_samples
            !! Number of samples/entities
        integer(int32), intent(in) :: n_timepoints
            !! Number of time points
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! Original trajectories
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(out) :: trajectories_norm
            !! Normalized trajectories
        integer(int32), dimension(n_factors, n_samples), intent(out) :: status
            !! Status code for specific warnings, one per factor per sample
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_series
        real(real64), dimension(:), allocatable :: tmp_series_norm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_series(n_timepoints))
        M_ALLOCATE(tmp_series_norm(n_timepoints))

        call normalize_all_trajectories_kernel(&
            trajectories = trajectories,&
            trajectories_norm = trajectories_norm,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            tmp_series = tmp_series,&
            tmp_series_norm = tmp_series_norm,&
            status = status&
        )
    end subroutine normalize_all_trajectories_alloc

end module tox_trajectory_normalization
