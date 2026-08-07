#include <src/macros.h>

!> summary: Wrappers for [[tox_trajectory_contribution_analysis_impl(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_trajectory_contribution_analysis
    use tox_trajectory_contribution_analysis_impl, only: MODE_MEAN, MODE_MIN, MODE_RAW, compute_acceleration_from_velocity_impl
    use tox_trajectory_contribution_analysis_impl, only: compute_acceleration_from_velocity_trajectory_impl, compute_all_contributions_impl, compute_baselines_factor_dependent_impl, compute_contributions_impl
    use tox_trajectory_contribution_analysis_impl, only: compute_p_values_impl, compute_velocity_acceleration_contributions_impl, compute_velocity_trajectories_impl, compute_velocity_trajectory_impl
    use tox_trajectory_contribution_analysis_impl, only: perform_permutation_test_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_int
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: perform_permutation_test
    public :: perform_permutation_test_expert
    public :: compute_p_values
    public :: compute_contributions
    public :: compute_all_contributions
    public :: compute_all_contributions_expert
    public :: compute_baselines_factor_dependent
    public :: compute_velocity_trajectory
    public :: compute_acceleration_from_velocity_trajectory
    public :: compute_velocity_trajectories
    public :: compute_acceleration_from_velocity
    public :: compute_velocity_acceleration_contributions
    public :: compute_velocity_acceleration_contributions_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_trajectory_contribution_analysis(module):perform_permutation_test_expert]] to prepare it yourself.
    subroutine perform_permutation_test(&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            factor_idx,&
            dependent_idx,&
            sample_idx,&
            baseline_mode,&
            n_permutations,&
            local_contributions,&
            total_contributions,&
            random_seed,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), intent(in) :: factor_idx
            !! index of factor to compute the permutation contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), intent(in) :: dependent_idx
            !! index of dependent to compute the permutation contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), intent(in) :: sample_idx
            !! index of sample to compute the permutation contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_timepoints, n_permutations), intent(out) :: local_contributions
            !! Per-timepoint contributions per permutation
        real(real64), dimension(n_permutations), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per permutation
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for random number generation.
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(:), allocatable :: tmp_factor
        real(real64), dimension(:), allocatable :: tmp_dependent

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_in_range_int(factor_idx, ierr, arg_pos=5_int32, min=1_int32, max=n_factors)
        call validate_in_range_int(dependent_idx, ierr, arg_pos=6_int32, min=1_int32, max=n_factors)
        call validate_in_range_int(sample_idx, ierr, arg_pos=7_int32, min=1_int32, max=n_samples)
        call validate_dimension_size(n_permutations, ierr, arg_pos=9_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=8_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_factor(n_timepoints))
        M_ALLOCATE(tmp_dependent(n_timepoints))

        call perform_permutation_test_impl(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            factor_idx = factor_idx,&
            dependent_idx = dependent_idx,&
            sample_idx = sample_idx,&
            baseline_mode = baseline_mode,&
            n_permutations = n_permutations,&
            local_contributions = local_contributions,&
            total_contributions = total_contributions,&
            tmp_factor = tmp_factor,&
            tmp_dependent = tmp_dependent,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine perform_permutation_test

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_trajectory_contribution_analysis(module):perform_permutation_test]] does both.
    subroutine perform_permutation_test_expert(&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            factor_idx,&
            dependent_idx,&
            sample_idx,&
            baseline_mode,&
            n_permutations,&
            local_contributions,&
            total_contributions,&
            tmp_factor,&
            tmp_dependent,&
            random_seed,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), intent(in) :: factor_idx
            !! index of factor to compute the permutation contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), intent(in) :: dependent_idx
            !! index of dependent to compute the permutation contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), intent(in) :: sample_idx
            !! index of sample to compute the permutation contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_samples`.
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_timepoints, n_permutations), intent(out) :: local_contributions
            !! Per-timepoint contributions per permutation
        real(real64), dimension(n_permutations), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per permutation
        real(real64), dimension(n_timepoints), intent(out) :: tmp_factor
            !! Working array to hold the factor in contiguous memory
        real(real64), dimension(n_timepoints), intent(out) :: tmp_dependent
            !! Working array to hold the random dependent in contiguous memory
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for random number generation.
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_in_range_int(factor_idx, ierr, arg_pos=5_int32, min=1_int32, max=n_factors)
        call validate_in_range_int(dependent_idx, ierr, arg_pos=6_int32, min=1_int32, max=n_factors)
        call validate_in_range_int(sample_idx, ierr, arg_pos=7_int32, min=1_int32, max=n_samples)
        call validate_dimension_size(n_permutations, ierr, arg_pos=9_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=8_int32)
        if (is_err(ierr)) return
#endif

        call perform_permutation_test_impl(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            factor_idx = factor_idx,&
            dependent_idx = dependent_idx,&
            sample_idx = sample_idx,&
            baseline_mode = baseline_mode,&
            n_permutations = n_permutations,&
            local_contributions = local_contributions,&
            total_contributions = total_contributions,&
            tmp_factor = tmp_factor,&
            tmp_dependent = tmp_dependent,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine perform_permutation_test_expert

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_p_values_impl]].
    !| Given the permutation tests ([[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test_impl(subroutine)]]),
    !| this counts how many of the permutation contributions were at least as high as the real ones.
    pure subroutine compute_p_values(&
            local_contributions_observed,&
            total_contribution_observed,&
            local_contributions_perm_test,&
            total_contributions_perm_test,&
            n_timepoints,&
            n_permutations,&
            local_p_values,&
            total_p_value,&
            ierr&
        )
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_timepoints), intent(in) :: local_contributions_observed
            !! Per-timepoint contributions for the observed factor-dependent-sample combination
        real(real64), intent(in) :: total_contribution_observed
            !! Total contribution (`sum(local_contributions)`) for the observed factor-dependent-sample combination
        real(real64), dimension(n_timepoints, n_permutations), intent(in) :: local_contributions_perm_test
            !! Per-timepoint contributions for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_permutations), intent(in) :: total_contributions_perm_test
            !! Total contribution (`sum(local_contributions)`) for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_timepoints), intent(out) :: local_p_values
            !! calculated p values for local contributions, like: `(local_contributions_perm_test >= local_contributions_observed)/n_permutations`
        real(real64), intent(out) :: total_p_value
            !! calculated p values for total contributions, like: `(total_contributions_perm_test >= total_contribution_observed)/n_permutations`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_real(total_contribution_observed, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_permutations, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(local_contributions_observed, n_timepoints, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(local_contributions_perm_test, n_timepoints * n_permutations, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(total_contributions_perm_test, n_permutations, ierr, arg_pos=4_int32)
        if (is_err(ierr)) return
#endif

        call compute_p_values_impl(&
            local_contributions_observed = local_contributions_observed,&
            total_contribution_observed = total_contribution_observed,&
            local_contributions_perm_test = local_contributions_perm_test,&
            total_contributions_perm_test = total_contributions_perm_test,&
            n_timepoints = n_timepoints,&
            n_permutations = n_permutations,&
            local_p_values = local_p_values,&
            total_p_value = total_p_value&
        )
    end subroutine compute_p_values

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_contributions_impl]].
    pure subroutine compute_contributions(&
            factor,&
            dependent,&
            n_dims,&
            baseline_mode,&
            local_contributions,&
            total_contribution,&
            ierr&
        )
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `factor` and `dependent`
        real(real64), dimension(n_dims), intent(in) :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_dims), intent(in) :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_dims), intent(out) :: local_contributions
            !! Per-element contributions
        real(real64), intent(out) :: total_contribution
            !! Total contribution (`sum(local_contributions)`)
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dims, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(factor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(dependent, n_dims, ierr, arg_pos=2_int32)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
        if (is_err(ierr)) return
#endif

        call compute_contributions_impl(&
            factor = factor,&
            dependent = dependent,&
            n_dims = n_dims,&
            baseline_mode = baseline_mode,&
            local_contributions = local_contributions,&
            total_contribution = total_contribution,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_contributions

    !> summary: Validates its inputs, prepares what [[tox_trajectory_contribution_analysis_impl(module):compute_all_contributions_impl]] needs, then calls it. The entry point to reach for first; see [[tox_trajectory_contribution_analysis(module):compute_all_contributions_expert]] to prepare it yourself.
    pure subroutine compute_all_contributions(&
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
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_selected_factors
            !! number of selected factors in `factor_indices`
        integer(int32), intent(in) :: n_selected_dependents
            !! number of selected dependents in `dependent_indices`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), dimension(n_selected_factors), intent(in) :: factor_indices
            !! indices of factors to compute the contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), dimension(n_selected_dependents), intent(in) :: dependent_indices
            !! indices of dependents to compute the contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_timepoints, n_selected_factors, n_selected_dependents, n_samples), intent(out) :: local_contributions
            !! Per-timepoint contributions per sample-dependent-factor combination
        real(real64), dimension(n_selected_factors, n_selected_dependents, n_samples), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(:, :), allocatable :: tmp_factors
        real(real64), dimension(:), allocatable :: tmp_dependent

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_selected_factors, ierr, arg_pos=6_int32)
        call validate_dimension_size(n_selected_dependents, ierr, arg_pos=8_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(factor_indices, n_selected_factors, ierr, arg_pos=5_int32, min=1_int32, max=n_factors)
        call validate_all_in_range_int(dependent_indices, n_selected_dependents, ierr, arg_pos=7_int32, min=1_int32, max=n_factors)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=9_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_factors(n_timepoints, n_selected_factors))
        M_ALLOCATE(tmp_dependent(n_timepoints))

        call compute_all_contributions_impl(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            factor_indices = factor_indices,&
            n_selected_factors = n_selected_factors,&
            dependent_indices = dependent_indices,&
            n_selected_dependents = n_selected_dependents,&
            baseline_mode = baseline_mode,&
            local_contributions = local_contributions,&
            total_contributions = total_contributions,&
            tmp_factors = tmp_factors,&
            tmp_dependent = tmp_dependent,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_all_contributions

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_all_contributions_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_trajectory_contribution_analysis(module):compute_all_contributions]] does both.
    pure subroutine compute_all_contributions_expert(&
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
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_selected_factors
            !! number of selected factors in `factor_indices`
        integer(int32), intent(in) :: n_selected_dependents
            !! number of selected dependents in `dependent_indices`
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), dimension(n_selected_factors), intent(in) :: factor_indices
            !! indices of factors to compute the contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), dimension(n_selected_dependents), intent(in) :: dependent_indices
            !! indices of dependents to compute the contributions for
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_factors`.
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_timepoints, n_selected_factors, n_selected_dependents, n_samples), intent(out) :: local_contributions
            !! Per-timepoint contributions per sample-dependent-factor combination
        real(real64), dimension(n_selected_factors, n_selected_dependents, n_samples), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per sample-dependent-factor combination
        real(real64), dimension(n_timepoints, n_selected_factors), intent(out) :: tmp_factors
            !! Working array to hold the currently handled sample's factors in contiguous memory
        real(real64), dimension(n_timepoints), intent(out) :: tmp_dependent
            !! Working array to hold the currently handled dependent in contiguous memory
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_selected_factors, ierr, arg_pos=6_int32)
        call validate_dimension_size(n_selected_dependents, ierr, arg_pos=8_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(factor_indices, n_selected_factors, ierr, arg_pos=5_int32, min=1_int32, max=n_factors)
        call validate_all_in_range_int(dependent_indices, n_selected_dependents, ierr, arg_pos=7_int32, min=1_int32, max=n_factors)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=9_int32)
        if (is_err(ierr)) return
#endif

        call compute_all_contributions_impl(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            factor_indices = factor_indices,&
            n_selected_factors = n_selected_factors,&
            dependent_indices = dependent_indices,&
            n_selected_dependents = n_selected_dependents,&
            baseline_mode = baseline_mode,&
            local_contributions = local_contributions,&
            total_contributions = total_contributions,&
            tmp_factors = tmp_factors,&
            tmp_dependent = tmp_dependent,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_all_contributions_expert

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_baselines_factor_dependent_impl]].
    pure subroutine compute_baselines_factor_dependent(&
            n_timepoints,&
            factor,&
            dependent,&
            baseline_mode,&
            factor_baseline,&
            dependent_baseline,&
            ierr&
        )
        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in both factor and dependent arrays
        real(real64), dimension(n_timepoints), intent(in) :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_timepoints), intent(in) :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), intent(out) :: factor_baseline
            !! Computed baseline for factor
        real(real64), intent(out) :: dependent_baseline
            !! Computed baseline for dependent variable
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_timepoints, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(factor, n_timepoints, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(dependent, n_timepoints, ierr, arg_pos=3_int32)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
        if (is_err(ierr)) return
#endif

        call compute_baselines_factor_dependent_impl(&
            n_timepoints = n_timepoints,&
            factor = factor,&
            dependent = dependent,&
            baseline_mode = baseline_mode,&
            factor_baseline = factor_baseline,&
            dependent_baseline = dependent_baseline,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_baselines_factor_dependent

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_velocity_trajectory_impl]].
    pure subroutine compute_velocity_trajectory(&
            trajectory,&
            velocity,&
            n_timepoints,&
            ierr&
        )
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_timepoints), intent(in) :: trajectory
            !! input position trajectory
        real(real64), dimension(max(0_int32, n_timepoints - 1)), intent(out) :: velocity
            !! output velocity trajectory
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_timepoints, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(trajectory, n_timepoints, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call compute_velocity_trajectory_impl(&
            trajectory = trajectory,&
            velocity = velocity,&
            n_timepoints = n_timepoints&
        )
    end subroutine compute_velocity_trajectory

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_acceleration_from_velocity_trajectory_impl]].
    pure subroutine compute_acceleration_from_velocity_trajectory(&
            velocity,&
            acceleration,&
            n_timepoints,&
            ierr&
        )
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(max(0_int32, n_timepoints - 1)), intent(in) :: velocity
            !! velocity trajectory
        real(real64), dimension(max(0_int32, n_timepoints - 2)), intent(out) :: acceleration
            !! acceleration trajectory
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_all_in_range_real(velocity, (max(0_int32, n_timepoints - 1)), ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call compute_acceleration_from_velocity_trajectory_impl(&
            velocity = velocity,&
            acceleration = acceleration,&
            n_timepoints = n_timepoints&
        )
    end subroutine compute_acceleration_from_velocity_trajectory

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_velocity_trajectories_impl]].
    pure subroutine compute_velocity_trajectories(&
            trajectories,&
            velocity,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! input position trajectories
        real(real64), dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples), intent(out) :: velocity
            !! output velocity trajectories
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

        call compute_velocity_trajectories_impl(&
            trajectories = trajectories,&
            velocity = velocity,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints&
        )
    end subroutine compute_velocity_trajectories

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_acceleration_from_velocity_impl]].
    pure subroutine compute_acceleration_from_velocity(&
            velocity,&
            acceleration,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples), intent(in) :: velocity
            !! input velocity trajectories
        real(real64), dimension(max(0_int32, n_timepoints - 2), n_factors, n_samples), intent(out) :: acceleration
            !! output acceleration trajectories
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(velocity, (max(0_int32, n_timepoints - 1)) * n_factors * n_samples, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        call compute_acceleration_from_velocity_impl(&
            velocity = velocity,&
            acceleration = acceleration,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints&
        )
    end subroutine compute_acceleration_from_velocity

    !> summary: Validates its inputs, prepares what [[tox_trajectory_contribution_analysis_impl(module):compute_velocity_acceleration_contributions_impl]] needs, then calls it. The entry point to reach for first; see [[tox_trajectory_contribution_analysis(module):compute_velocity_acceleration_contributions_expert]] to prepare it yourself.
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
    pure subroutine compute_velocity_acceleration_contributions(&
            trajectories,&
            n_factors,&
            n_samples,&
            n_timepoints,&
            baseline_mode,&
            contrib_velocity,&
            velocity_contribution_series,&
            contrib_acceleration,&
            acceleration_contribution_series,&
            ierr&
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! input position trajectories
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_factors, n_factors, n_samples), intent(out) :: contrib_velocity
            !! output velocity contributions
        real(real64), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out) :: velocity_contribution_series
            !! output velocity contribution series
        real(real64), dimension(n_factors, n_factors, n_samples), intent(out) :: contrib_acceleration
            !! output acceleration contributions
        real(real64), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out) :: acceleration_contribution_series
            !! output acceleration contribution series
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(:, :), allocatable :: tmp_factors
        real(real64), dimension(:), allocatable :: tmp_dependent
        real(real64), dimension(:), allocatable :: tmp_contributions

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_factors(n_timepoints - 1, n_factors))
        M_ALLOCATE(tmp_dependent(n_timepoints - 1))
        M_ALLOCATE(tmp_contributions(n_timepoints - 1))

        call compute_velocity_acceleration_contributions_impl(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            baseline_mode = baseline_mode,&
            tmp_factors = tmp_factors,&
            tmp_dependent = tmp_dependent,&
            tmp_contributions = tmp_contributions,&
            contrib_velocity = contrib_velocity,&
            velocity_contribution_series = velocity_contribution_series,&
            contrib_acceleration = contrib_acceleration,&
            acceleration_contribution_series = acceleration_contribution_series,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_velocity_acceleration_contributions

    !> summary: Validates its inputs, then calls [[tox_trajectory_contribution_analysis_impl(module):compute_velocity_acceleration_contributions_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_trajectory_contribution_analysis(module):compute_velocity_acceleration_contributions]] does both.
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
    pure subroutine compute_velocity_acceleration_contributions_expert(&
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
        )
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! input position trajectories
        integer(int32), intent(in) :: baseline_mode
            !! | Mode              | Value                                                                     |
            !! |-------------------|---------------------------------------------------------------------------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis_impl(module):MODE_RAW(variable)]]  |
            !! | arithmetic mean   | [[tox_trajectory_contribution_analysis_impl(module):MODE_MEAN(variable)]] |
            !! | minimum value     | [[tox_trajectory_contribution_analysis_impl(module):MODE_MIN(variable)]]  |
        real(real64), dimension(n_timepoints - 1, n_factors), intent(out) :: tmp_factors
            !! workspace for factor data (used for velocity and acceleration)
        real(real64), dimension(n_timepoints - 1), intent(out) :: tmp_dependent
            !! workspace for dependent data (used for velocity and acceleration)
        real(real64), dimension(n_timepoints - 1), intent(out) :: tmp_contributions
            !! workspace for contribution calculations (used for velocity and acceleration)
        real(real64), dimension(n_factors, n_factors, n_samples), intent(out) :: contrib_velocity
            !! output velocity contributions
        real(real64), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out) :: velocity_contribution_series
            !! output velocity contribution series
        real(real64), dimension(n_factors, n_factors, n_samples), intent(out) :: contrib_acceleration
            !! output acceleration contributions
        real(real64), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out) :: acceleration_contribution_series
            !! output acceleration contribution series
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(trajectories, n_factors * n_samples * n_timepoints, ierr, arg_pos=1_int32)
        if (baseline_mode /= MODE_RAW .and. baseline_mode /= MODE_MEAN .and. baseline_mode /= MODE_MIN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (is_err(ierr)) return
#endif

        call compute_velocity_acceleration_contributions_impl(&
            trajectories = trajectories,&
            n_factors = n_factors,&
            n_samples = n_samples,&
            n_timepoints = n_timepoints,&
            baseline_mode = baseline_mode,&
            tmp_factors = tmp_factors,&
            tmp_dependent = tmp_dependent,&
            tmp_contributions = tmp_contributions,&
            contrib_velocity = contrib_velocity,&
            velocity_contribution_series = velocity_contribution_series,&
            contrib_acceleration = contrib_acceleration,&
            acceleration_contribution_series = acceleration_contribution_series,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_velocity_acceleration_contributions_expert

end module tox_trajectory_contribution_analysis
