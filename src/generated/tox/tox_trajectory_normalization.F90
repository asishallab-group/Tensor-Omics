#include <src/macros.h>

!> Min-max normalization of factor trajectories over time.
!|
!| Each factor's time series is rescaled to `[0,1]` independently, per sample -- so trajectories
!| of very different magnitudes can be compared by shape. `normalize_single_trajectory` does one;
!| `normalize_all_trajectories` does a whole tensor of them in one call.
!|
!| Generated from [[tox_trajectory_normalization_impl(module)]]; do not edit -- regenerate instead.
module tox_trajectory_normalization
    use tox_trajectory_normalization_impl, only: normalize_all_trajectories_impl, normalize_single_trajectory_impl, normalize_variable_timeseries_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: normalize_variable_timeseries
    public :: normalize_single_trajectory
    public :: normalize_all_trajectories
    public :: normalize_all_trajectories_expert

contains

    !> summary: Validates its inputs, then calls [[tox_trajectory_normalization_impl(module):normalize_variable_timeseries_impl]].
    pure subroutine normalize_variable_timeseries(&
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

        call normalize_variable_timeseries_impl(&
            v = v,&
            v_norm = v_norm,&
            n_points = n_points,&
            status = status&
        )
    end subroutine normalize_variable_timeseries

    !> summary: Validates its inputs, then calls [[tox_trajectory_normalization_impl(module):normalize_single_trajectory_impl]].
    pure subroutine normalize_single_trajectory(&
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

        call normalize_single_trajectory_impl(&
            trajectory = trajectory,&
            trajectory_norm = trajectory_norm,&
            n_factors = n_factors,&
            n_timepoints = n_timepoints,&
            status = status&
        )
    end subroutine normalize_single_trajectory

    !> summary: Validates its inputs, prepares what [[tox_trajectory_normalization_impl(module):normalize_all_trajectories_impl]] needs, then calls it. The entry point to reach for first; see [[tox_trajectory_normalization(module):normalize_all_trajectories_expert]] to prepare it yourself.
    !| independently across time for each sample.
    pure subroutine normalize_all_trajectories(&
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

        call normalize_all_trajectories_impl(&
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

    !> summary: Validates its inputs, then calls [[tox_trajectory_normalization_impl(module):normalize_all_trajectories_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_trajectory_normalization(module):normalize_all_trajectories]] does both.
    !| independently across time for each sample.
    pure subroutine normalize_all_trajectories_expert(&
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

        call normalize_all_trajectories_impl(&
            trajectories = trajectories,&
            trajectories_norm = trajectories_norm,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            tmp_series = tmp_series,&
            tmp_series_norm = tmp_series_norm,&
            status = status&
        )
    end subroutine normalize_all_trajectories_expert

end module tox_trajectory_normalization
