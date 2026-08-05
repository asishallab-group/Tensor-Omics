#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_trajectory_contribution_analysis_kernel(module)]]
!| Module for quantifying how much one trajectory (a "factor") contributes to another (a "dependent") over time.
!|
!| Contributions are computed per timepoint as the product of both series' deviations from a chosen
!| baseline, for raw expression trajectories as well as for their velocity (first difference) and
!| acceleration (second difference) derivatives. Statistical significance of an observed contribution can
!| be assessed via a permutation test that recomputes the same contribution against a randomly chosen
!| other sample.
module tox_trajectory_contribution_analysis_kernel_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_1d_as_string
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: compute_all_contributions_kernel_c
    public :: compute_baselines_factor_dependent_kernel_c
    public :: compute_velocity_trajectory_kernel_c
    public :: compute_acceleration_from_velocity_trajectory_kernel_c
    public :: compute_velocity_acceleration_contributions_kernel_c

contains

    !> summary: C-wrapper for [[tox_trajectory_contribution_analysis_kernel(module):compute_all_contributions_kernel(subroutine)]]
    subroutine compute_all_contributions_kernel_c(&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            factor_indices,&
            n_selected_factors,&
            dependent_indices,&
            n_selected_dependents,&
            baseline_mode,&
            local_contributions,&
            total_contributions,&
            tmp_factors,&
            tmp_dependent,&
            ierr&
        ) bind(C, name="compute_all_contributions_kernel_c")
        use tox_trajectory_contribution_analysis_kernel, only: compute_all_contributions_kernel
        use tox_trajectory_contribution_analysis_kernel, only: MODE_MEAN, MODE_MIN, MODE_RAW

        integer(c_int), intent(in), target :: n_factors
            !! number of factors
        integer(c_int), intent(in), target :: n_samples
            !! number of samples
        integer(c_int), intent(in), target :: n_timepoints
            !! number of timepoints
        integer(c_int), intent(in), target :: n_selected_factors
            !! number of selected factors in `factor_indices`
        integer(c_int), intent(in), target :: n_selected_dependents
            !! number of selected dependents in `dependent_indices`
        real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
            !! expression vectors across different samples over time
        integer(c_int), dimension(n_selected_factors), intent(in), target :: factor_indices
            !! indices of factors to compute the contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(c_int), dimension(n_selected_dependents), intent(in), target :: dependent_indices
            !! indices of dependents to compute the contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        character(len=1, kind=c_char), dimension(4), intent(in), target :: baseline_mode
            !! | Mode              | Value                                                                       |
            !! |-------------------|-----------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_kernel(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_kernel(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_kernel(module):MODE_MIN(variable)]]  |
        real(c_double), dimension(n_timepoints, n_selected_factors, n_selected_dependents, n_samples), intent(out), target :: local_contributions
            !! Per-timepoint contributions per sample-dependent-factor combination
        real(c_double), dimension(n_selected_factors, n_selected_dependents, n_samples), intent(out), target :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination
        real(c_double), dimension(n_timepoints, n_selected_factors), intent(out), target :: tmp_factors
            !! Working array to hold the currently handled sample's factors in contiguous memory
        real(c_double), dimension(n_timepoints), intent(out), target :: tmp_dependent
            !! Working array to hold the currently handled dependent in contiguous memory
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(int32) :: baseline_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_NON_NULL(n_selected_factors)
        M_CHECK_NON_NULL(n_selected_dependents)
        M_CHECK_ARRAY_NON_NULL(trajectories, n_factors * n_samples * n_timepoints)
        M_CHECK_ARRAY_NON_NULL(factor_indices, n_selected_factors)
        M_CHECK_ARRAY_NON_NULL(dependent_indices, n_selected_dependents)
        M_CHECK_ARRAY_NON_NULL(baseline_mode, 4)
        M_CHECK_ARRAY_NON_NULL(local_contributions, n_timepoints * n_selected_factors * n_selected_dependents * n_samples)
        M_CHECK_ARRAY_NON_NULL(total_contributions, n_selected_factors * n_selected_dependents * n_samples)
        M_CHECK_ARRAY_NON_NULL(tmp_factors, n_timepoints * n_selected_factors)
        M_CHECK_ARRAY_NON_NULL(tmp_dependent, n_timepoints)

        block
            character(len=:), allocatable :: baseline_mode_f
            call c_char_1d_as_string(baseline_mode, baseline_mode_f, ierr)
            if (is_err(ierr)) return

            select case (baseline_mode_f)
                case ("raw")
                    baseline_mode_mode_f = MODE_RAW
                case ("mean")
                    baseline_mode_mode_f = MODE_MEAN
                case ("min")
                    baseline_mode_mode_f = MODE_MIN
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call compute_all_contributions_kernel(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            factor_indices = factor_indices,&
            n_selected_factors = n_selected_factors,&
            dependent_indices = dependent_indices,&
            n_selected_dependents = n_selected_dependents,&
            baseline_mode = baseline_mode_mode_f,&
            local_contributions = local_contributions,&
            total_contributions = total_contributions,&
            tmp_factors = tmp_factors,&
            tmp_dependent = tmp_dependent,&
            ierr = ierr&
        )
    end subroutine compute_all_contributions_kernel_c

    !> summary: C-wrapper for [[tox_trajectory_contribution_analysis_kernel(module):compute_baselines_factor_dependent_kernel(subroutine)]]
    subroutine compute_baselines_factor_dependent_kernel_c(&
            n_timepoints,&
            factor,&
            dependent,&
            baseline_mode,&
            factor_baseline,&
            dependent_baseline,&
            ierr&
        ) bind(C, name="compute_baselines_factor_dependent_kernel_c")
        use tox_trajectory_contribution_analysis_kernel, only: compute_baselines_factor_dependent_kernel
        use tox_trajectory_contribution_analysis_kernel, only: MODE_MEAN, MODE_MIN, MODE_RAW

        integer(c_int), intent(in), target :: n_timepoints
            !! Number of timepoints in both factor and dependent arrays
        real(c_double), dimension(n_timepoints), intent(in), target :: factor
            !! Factor time series, length n_timepoints
        real(c_double), dimension(n_timepoints), intent(in), target :: dependent
            !! Dependent variable time series, length n_timepoints
        character(len=1, kind=c_char), dimension(4), intent(in), target :: baseline_mode
            !! | Mode              | Value                                                                       |
            !! |-------------------|-----------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_kernel(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_kernel(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_kernel(module):MODE_MIN(variable)]]  |
        real(c_double), intent(out), target :: factor_baseline
            !! Computed baseline for factor
        real(c_double), intent(out), target :: dependent_baseline
            !! Computed baseline for dependent variable
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(int32) :: baseline_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_NON_NULL(factor_baseline)
        M_CHECK_NON_NULL(dependent_baseline)
        M_CHECK_ARRAY_NON_NULL(factor, n_timepoints)
        M_CHECK_ARRAY_NON_NULL(dependent, n_timepoints)
        M_CHECK_ARRAY_NON_NULL(baseline_mode, 4)

        block
            character(len=:), allocatable :: baseline_mode_f
            call c_char_1d_as_string(baseline_mode, baseline_mode_f, ierr)
            if (is_err(ierr)) return

            select case (baseline_mode_f)
                case ("raw")
                    baseline_mode_mode_f = MODE_RAW
                case ("mean")
                    baseline_mode_mode_f = MODE_MEAN
                case ("min")
                    baseline_mode_mode_f = MODE_MIN
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call compute_baselines_factor_dependent_kernel(&
            n_timepoints = n_timepoints,&
            factor = factor,&
            dependent = dependent,&
            baseline_mode = baseline_mode_mode_f,&
            factor_baseline = factor_baseline,&
            dependent_baseline = dependent_baseline,&
            ierr = ierr&
        )
    end subroutine compute_baselines_factor_dependent_kernel_c

    !> summary: C-wrapper for [[tox_trajectory_contribution_analysis_kernel(module):compute_velocity_trajectory_kernel(subroutine)]]
    subroutine compute_velocity_trajectory_kernel_c(&
            trajectory,&
            velocity,&
            n_timepoints,&
            ierr&
        ) bind(C, name="compute_velocity_trajectory_kernel_c")
        use tox_trajectory_contribution_analysis_kernel, only: compute_velocity_trajectory_kernel

        integer(c_int), intent(in), target :: n_timepoints
            !! number of timepoints
        real(c_double), dimension(n_timepoints), intent(in), target :: trajectory
            !! input position trajectory
        real(c_double), dimension(max(0_int32, n_timepoints - 1)), intent(out), target :: velocity
            !! output velocity trajectory
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_ARRAY_NON_NULL(trajectory, n_timepoints)
        M_CHECK_ARRAY_NON_NULL(velocity, (max(0_int32, n_timepoints - 1)))

        call compute_velocity_trajectory_kernel(&
            trajectory = trajectory,&
            velocity = velocity,&
            n_timepoints = n_timepoints&
        )
    end subroutine compute_velocity_trajectory_kernel_c

    !> summary: C-wrapper for [[tox_trajectory_contribution_analysis_kernel(module):compute_acceleration_from_velocity_trajectory_kernel(subroutine)]]
    subroutine compute_acceleration_from_velocity_trajectory_kernel_c(&
            velocity,&
            acceleration,&
            n_timepoints,&
            ierr&
        ) bind(C, name="compute_acceleration_from_velocity_trajectory_kernel_c")
        use tox_trajectory_contribution_analysis_kernel, only: compute_acceleration_from_velocity_trajectory_kernel

        integer(c_int), intent(in), target :: n_timepoints
            !! number of timepoints
        real(c_double), dimension(max(0_int32, n_timepoints - 1)), intent(in), target :: velocity
            !! velocity trajectory
        real(c_double), dimension(max(0_int32, n_timepoints - 2)), intent(out), target :: acceleration
            !! acceleration trajectory
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_ARRAY_NON_NULL(velocity, (max(0_int32, n_timepoints - 1)))
        M_CHECK_ARRAY_NON_NULL(acceleration, (max(0_int32, n_timepoints - 2)))

        call compute_acceleration_from_velocity_trajectory_kernel(&
            velocity = velocity,&
            acceleration = acceleration,&
            n_timepoints = n_timepoints&
        )
    end subroutine compute_acceleration_from_velocity_trajectory_kernel_c

    !> summary: C-wrapper for [[tox_trajectory_contribution_analysis_kernel(module):compute_velocity_acceleration_contributions_kernel(subroutine)]]
    !| @note
    !| Performance layout:
    !|
    !| `trajectories` uses `(n_factors, n_samples, n_timepoints)`.
    !| Velocity and acceleration use time-first layouts:
    !|
    !| - `velocity`     -> `(max(0, n_timepoints-1), n_factors, n_samples)`
    !| - `acceleration` -> `(max(0, n_timepoints-2), n_factors, n_samples)`
    !|
    !| This keeps slices like `velocity(:, factor, sample)` contiguous,
    !| avoids expensive tmporaries, and improves cache efficiency.
    !| @endnote
    subroutine compute_velocity_acceleration_contributions_kernel_c(&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            baseline_mode,&
            tmp_factors,&
            tmp_dependent,&
            tmp_contributions,&
            contrib_velocity,&
            velocity_contribution_series,&
            contrib_acceleration,&
            acceleration_contribution_series,&
            ierr&
        ) bind(C, name="compute_velocity_acceleration_contributions_kernel_c")
        use tox_trajectory_contribution_analysis_kernel, only: compute_velocity_acceleration_contributions_kernel
        use tox_trajectory_contribution_analysis_kernel, only: MODE_MEAN, MODE_MIN, MODE_RAW

        integer(c_int), intent(in), target :: n_factors
            !! number of factors
        integer(c_int), intent(in), target :: n_samples
            !! number of samples
        integer(c_int), intent(in), target :: n_timepoints
            !! number of timepoints
        real(c_double), dimension(n_factors, n_samples, n_timepoints), intent(in), target :: trajectories
            !! input position trajectories
        character(len=1, kind=c_char), dimension(4), intent(in), target :: baseline_mode
            !! | Mode              | Value                                                                       |
            !! |-------------------|-----------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_kernel(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_kernel(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_kernel(module):MODE_MIN(variable)]]  |
        real(c_double), dimension(n_timepoints - 1, n_factors), intent(out), target :: tmp_factors
            !! workspace for factor data (used for velocity and acceleration)
        real(c_double), dimension(n_timepoints - 1), intent(out), target :: tmp_dependent
            !! workspace for dependent data (used for velocity and acceleration)
        real(c_double), dimension(n_timepoints - 1), intent(out), target :: tmp_contributions
            !! workspace for contribution calculations (used for velocity and acceleration)
        real(c_double), dimension(n_factors, n_factors, n_samples), intent(out), target :: contrib_velocity
            !! output velocity contributions
        real(c_double), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out), target :: velocity_contribution_series
            !! output velocity contribution series
        real(c_double), dimension(n_factors, n_factors, n_samples), intent(out), target :: contrib_acceleration
            !! output acceleration contributions
        real(c_double), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out), target :: acceleration_contribution_series
            !! output acceleration contribution series
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(int32) :: baseline_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_factors)
        M_CHECK_NON_NULL(n_samples)
        M_CHECK_NON_NULL(n_timepoints)
        M_CHECK_ARRAY_NON_NULL(trajectories, n_factors * n_samples * n_timepoints)
        M_CHECK_ARRAY_NON_NULL(baseline_mode, 4)
        M_CHECK_ARRAY_NON_NULL(tmp_factors, (n_timepoints - 1) * n_factors)
        M_CHECK_ARRAY_NON_NULL(tmp_dependent, (n_timepoints - 1))
        M_CHECK_ARRAY_NON_NULL(tmp_contributions, (n_timepoints - 1))
        M_CHECK_ARRAY_NON_NULL(contrib_velocity, n_factors * n_factors * n_samples)
        M_CHECK_ARRAY_NON_NULL(velocity_contribution_series, n_timepoints * n_factors * n_factors * n_samples)
        M_CHECK_ARRAY_NON_NULL(contrib_acceleration, n_factors * n_factors * n_samples)
        M_CHECK_ARRAY_NON_NULL(acceleration_contribution_series, n_timepoints * n_factors * n_factors * n_samples)

        block
            character(len=:), allocatable :: baseline_mode_f
            call c_char_1d_as_string(baseline_mode, baseline_mode_f, ierr)
            if (is_err(ierr)) return

            select case (baseline_mode_f)
                case ("raw")
                    baseline_mode_mode_f = MODE_RAW
                case ("mean")
                    baseline_mode_mode_f = MODE_MEAN
                case ("min")
                    baseline_mode_mode_f = MODE_MIN
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call compute_velocity_acceleration_contributions_kernel(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            baseline_mode = baseline_mode_mode_f,&
            tmp_factors = tmp_factors,&
            tmp_dependent = tmp_dependent,&
            tmp_contributions = tmp_contributions,&
            contrib_velocity = contrib_velocity,&
            velocity_contribution_series = velocity_contribution_series,&
            contrib_acceleration = contrib_acceleration,&
            acceleration_contribution_series = acceleration_contribution_series,&
            ierr = ierr&
        )
    end subroutine compute_velocity_acceleration_contributions_kernel_c

end module tox_trajectory_contribution_analysis_kernel_c
#endif
