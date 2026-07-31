#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_trajectory_normalization(module)]]
!| Module for min-max normalization of factor trajectories over time.
!| Each factor's time series is independently rescaled to `[0,1]` per sample/entity.
module tox_trajectory_normalization_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: normalize_variable_timeseries_c
    public :: normalize_single_trajectory_c
    public :: normalize_all_trajectories_c

contains

    !> summary: C-wrapper for [[tox_trajectory_normalization(module):normalize_variable_timeseries(subroutine)]]
    subroutine normalize_variable_timeseries_c(&
            v,&
            v_norm,&
            n_points,&
            ierr,&
            status&
        ) bind(C, name="normalize_variable_timeseries_c")
        use tox_trajectory_normalization, only: normalize_variable_timeseries

        integer(c_int), intent(in), target :: n_points
            !! Vector length (number of time points)
        real(c_double), dimension(n_points), intent(in), target :: v
            !! Original time series
        real(c_double), dimension(n_points), intent(out), target :: v_norm
            !! Normalized time series
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), intent(out), target :: status
            !! Status code for specific warnings

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(status)
        M_CHECK_ARRAY_NON_NULL(v, n_points)
        M_CHECK_ARRAY_NON_NULL(v_norm, n_points)

        call normalize_variable_timeseries(&
            v = v,&
            v_norm = v_norm,&
            n_points = n_points,&
            ierr = ierr,&
            status = status&
        )
    end subroutine normalize_variable_timeseries_c

    !> summary: C-wrapper for [[tox_trajectory_normalization(module):normalize_single_trajectory(subroutine)]]
    subroutine normalize_single_trajectory_c(&
            trajectory,&
            trajectory_norm,&
            n_factors,&
            n_timepoints,&
            ierr,&
            status&
        ) bind(C, name="normalize_single_trajectory_c")
        use tox_trajectory_normalization, only: normalize_single_trajectory

        integer(c_int), intent(in), target :: n_factors
            !! Number of factors/variables
        integer(c_int), intent(in), target :: n_timepoints
            !! Number of time points
        real(c_double), dimension(n_timepoints, n_factors), intent(in), target :: trajectory
            !! Original trajectory for one sample
        real(c_double), dimension(n_timepoints, n_factors), intent(out), target :: trajectory_norm
            !! Normalized trajectory for one sample
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), dimension(n_factors), intent(out), target :: status
            !! Status code for specific warnings

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_ARRAY_NON_NULL(trajectory, n_timepoints * n_factors)
        M_CHECK_ARRAY_NON_NULL(trajectory_norm, n_timepoints * n_factors)
        M_CHECK_ARRAY_NON_NULL(status, n_factors)

        call normalize_single_trajectory(&
            trajectory = trajectory,&
            trajectory_norm = trajectory_norm,&
            n_factors = n_factors,&
            n_timepoints = n_timepoints,&
            ierr = ierr,&
            status = status&
        )
    end subroutine normalize_single_trajectory_c

    !> summary: C-wrapper for [[tox_trajectory_normalization(module):normalize_all_trajectories_alloc(subroutine)]]
    !| Normalizes each factor independently across time for each sample
    subroutine normalize_all_trajectories_c(&
            trajectories,&
            trajectories_norm,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            ierr,&
            status&
        ) bind(C, name="normalize_all_trajectories_c")
        use tox_trajectory_normalization, only: normalize_all_trajectories_alloc

        integer(c_int), intent(in), target :: n_factors
            !! Number of factors
        integer(c_int), intent(in), target :: n_samples
            !! Number of samples/entities
        integer(c_int), intent(in), target :: n_timepoints
            !! Number of time points
        real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
            !! Original trajectories
        real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(out), target :: trajectories_norm
            !! Normalized trajectories
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), dimension(n_factors, n_samples), intent(out), target :: status
            !! Status code for specific warnings

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_ARRAY_NON_NULL(trajectories, n_factors * n_samples * n_timepoints)
        M_CHECK_ARRAY_NON_NULL(trajectories_norm, n_factors * n_samples * n_timepoints)
        M_CHECK_ARRAY_NON_NULL(status, n_factors * n_samples)

        call normalize_all_trajectories_alloc(&
            trajectories = trajectories,&
            trajectories_norm = trajectories_norm,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            ierr = ierr,&
            status = status&
        )
    end subroutine normalize_all_trajectories_c

end module tox_trajectory_normalization_c
#endif
