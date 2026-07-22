#include <src/macros.h>

!> Module for quantifying how much one trajectory (a "factor") contributes to another (a "dependent") over time.
!|
!| Contributions are computed per timepoint as the product of both series' deviations from a chosen
!| baseline, for raw expression trajectories as well as for their velocity (first difference) and
!| acceleration (second difference) derivatives. Statistical significance of an observed contribution can
!| be assessed via a permutation test that recomputes the same contribution against a randomly chosen
!| other sample.
module tox_trajectory_contribution_analysis
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int, validate_in_range_real, validate_in_range_int
    use f42_utils, only: init_random, rand_range
    M_IMPLICIT_NONE

    ! Baseline computation modes
#define CM_VALIDATE_MODE_BASELINE(ARG) call validate_in_range_int(baseline_mode, ierr, min=1_int32, max=3_int32, ARG)
#define CM_MODE_BASELINE_RAW 1_int32
#define CM_MODE_BASELINE_MIN 2_int32
#define CM_MODE_BASELINE_MEAN 3_int32
    integer(int32), parameter :: BASELINE_RAW = CM_MODE_BASELINE_RAW
        !! Baseline mode code: no centering, contributions are computed from raw values.
    integer(int32), parameter :: BASELINE_MIN = CM_MODE_BASELINE_MIN
        !! Baseline mode code: contributions are centered on each series' minimum value.
    integer(int32), parameter :: BASELINE_MEAN = CM_MODE_BASELINE_MEAN
        !! Baseline mode code: contributions are centered on each series' arithmetic mean.

    ! Public aliases the interface generator maps the `baseline_mode` strings onto: an arg
    ! named `*_mode` takes the `MODE_` prefix, so "raw"/"min"/"mean" resolve to these.
    integer(int32), parameter, public :: MODE_RAW = CM_MODE_BASELINE_RAW
        !! `baseline_mode` value for no centering (raw values)
    integer(int32), parameter, public :: MODE_MIN = CM_MODE_BASELINE_MIN
        !! `baseline_mode` value for centering on each series' minimum value
    integer(int32), parameter, public :: MODE_MEAN = CM_MODE_BASELINE_MEAN
        !! `baseline_mode` value for centering on each series' arithmetic mean

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Selects a random sample different to `current_sample`. For reproducibility call [[f42_utils(module):init_random(subroutine)]] beforehand.
    subroutine select_random_sample(n_samples, current_sample, random_sample, ierr)
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: current_sample
            !! sample that won't be selected
        integer(int32), intent(out) :: random_sample
            !! selected random sample
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_in_range_int(n_samples, ierr, min=2_int32)
        call validate_in_range_int(current_sample, ierr, min=1_int32, max=n_samples)

        if (is_err(ierr)) return

        call select_random_sample_helper(n_samples, current_sample, random_sample)
    end subroutine select_random_sample

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Selects a random sample different to `current_sample`. For reproducibility call [[f42_utils(module):init_random(subroutine)]] beforehand.
    subroutine select_random_sample_helper(n_samples, current_sample, random_sample)
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: current_sample
            !! sample that won't be selected
        integer(int32), intent(out) :: random_sample
            !! selected random sample

        ! pick random sample in range [1,n_samples-1], so a set without `current_sample`
        random_sample = int(rand_range(1.0_real64, real(n_samples, real64)), kind=int32)

        ! adjust picked sample for excluded `current_sample`
        if (random_sample >= current_sample) then
            random_sample = random_sample + 1
        end if
    end subroutine select_random_sample_helper

    !> M_EXPORT_C
    !| summary: Permutation test for a factor-dependent pair using a random different sample's dependent
    !| AUTHOR_FRANZ_ERIC_SILL
    !| For a given factor-dependent pair, this subroutine calculates the contributions by taking the same dependent but from a random different sample.
    subroutine perform_permutation_test(trajectories, n_factors, n_samples, n_timepoints, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, local_contributions, total_contributions, tmp_factor, tmp_dependent, ierr, random_seed)
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), intent(in) :: factor_idx
            !! index of factor to compute the permutation contributions for
        integer(int32), intent(in) :: dependent_idx
            !! index of dependent to compute the permutation contributions for
        integer(int32), intent(in) :: sample_idx
            !! index of sample to compute the permutation contributions for
        integer(int32), intent(in) :: baseline_mode
            !! Used mode for baseline calculation (see [[tox_trajectory_contribution_analysis(module):compute_baselines_factor_dependent(subroutine)]])
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_timepoints, n_permutations), intent(out) :: local_contributions
            !! Per-timepoint contributions per permutation
        real(real64), dimension(n_permutations), intent(out) :: total_contributions
            !! Total contribution (`sum(local_contributions)`) per permutation
        real(real64), dimension(n_timepoints), intent(out) :: tmp_factor
            !! Working array to hold the factor in contiguous memory
        real(real64), dimension(n_timepoints), intent(out) :: tmp_dependent
            !! Working array to hold the random dependent in contiguous memory
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for random number generation.

        integer(int32) :: random_sample, i_perm, i_timepoint

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_permutations, ierr, arg_pos=9_int32)
        call validate_all_in_range_real(trajectories, size(trajectories, kind=int32), ierr, arg_pos=1_int32)
        call validate_in_range_int(factor_idx, ierr, min=1, max=n_factors, arg_pos=5_int32)
        call validate_in_range_int(dependent_idx, ierr, min=1, max=n_factors, arg_pos=6_int32)
        call validate_in_range_int(sample_idx, ierr, min=1, max=n_samples, arg_pos=7_int32)
        CM_VALIDATE_MODE_BASELINE(arg_pos=8_int32)

        if (is_err(ierr)) return

        call perform_permutation_test_helper(trajectories, n_factors, n_samples, n_timepoints, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, local_contributions, total_contributions, tmp_factor, tmp_dependent, random_seed, ierr)
    end subroutine perform_permutation_test

    !> AUTHOR_FRANZ_ERIC_SILL
    !| (no input validation) For a given factor-dependent pair, this subroutine calculates the contributions by taking the same dependent but from a random different sample.
    subroutine perform_permutation_test_helper(trajectories, n_factors, n_samples, n_timepoints, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, local_contributions, total_contributions, tmp_factor, tmp_dependent, random_seed, ierr)
        integer(int32), intent(in) :: n_factors
            !! number of factors
        integer(int32), intent(in) :: n_samples
            !! number of samples
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! expression vectors across different samples over time
        integer(int32), intent(in) :: factor_idx
            !! index of factor to compute the permutation contributions for
        integer(int32), intent(in) :: dependent_idx
            !! index of dependent to compute the permutation contributions for
        integer(int32), intent(in) :: sample_idx
            !! index of sample to compute the permutation contributions for
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
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

        integer(int32) :: random_sample, i_perm, i_timepoint

        call set_ok(ierr)

        if (present(random_seed)) then
            call init_random(random_seed)
        end if

        do concurrent (i_timepoint = 1:n_timepoints) shared(tmp_factor, trajectories, factor_idx, sample_idx)
            tmp_factor(i_timepoint) = trajectories(factor_idx, sample_idx, i_timepoint)
        end do

        do i_perm = 1, n_permutations
            call select_random_sample_helper(n_samples, sample_idx, random_sample)

            do concurrent (i_timepoint = 1:n_timepoints) shared(tmp_dependent, trajectories, dependent_idx, random_sample)
                tmp_dependent(i_timepoint) = trajectories(dependent_idx, random_sample, i_timepoint)
            end do

            call compute_contributions_helper(tmp_factor, tmp_dependent, n_timepoints, baseline_mode, local_contributions(:, i_perm), total_contributions(i_perm), ierr)
            if (is_err(ierr)) return
        end do
    end subroutine perform_permutation_test_helper

    !> M_EXPORT_C
    !| summary: Computes empirical p-values from the permutation-test contributions
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Once the permutation tests are calculated ([[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]),
    !| this subroutine calculates the p values for the contributions, i.e. how many of the permutation contributions were at least as high as the real contributions.
    pure subroutine compute_p_values(local_contributions_observed, total_contribution_observed, local_contributions_perm_test, total_contributions_perm_test, n_timepoints, n_permutations, local_p_values, total_p_value, ierr)
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_timepoints), intent(in) :: local_contributions_observed
            !! Per-timepoint contributions for the observed factor-dependent-sample combination
        real(real64), intent(in) :: total_contribution_observed
            !! Total contribution (`sum(local_contributions)`) for the observed factor-dependent-sample combination
        real(real64), dimension(n_timepoints, n_permutations), intent(in) :: local_contributions_perm_test
            !! Per-timepoint contributions for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_permutations), intent(in) :: total_contributions_perm_test
            !! Total contribution (`sum(local_contributions)`) for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_timepoints), intent(out) :: local_p_values
            !! calculated p values for local contributions, like: `(local_contributions_perm_test >= local_contributions_observed)/n_permutations`
        real(real64), intent(out) :: total_p_value
            !! calculated p values for total contributions, like: `(total_contributions_perm_test >= total_contribution_observed)/n_permutations`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_permutations, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(local_contributions_observed, size(local_contributions_observed, kind=int32), ierr, arg_pos=1_int32)
        call validate_in_range_real(total_contribution_observed, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(local_contributions_perm_test, size(local_contributions_perm_test, kind=int32), ierr, arg_pos=3_int32)
        call validate_all_in_range_real(total_contributions_perm_test, size(total_contributions_perm_test, kind=int32), ierr, arg_pos=4_int32)

        if (is_err(ierr)) return

        call compute_p_values_helper(local_contributions_observed, total_contribution_observed, local_contributions_perm_test, total_contributions_perm_test, n_timepoints, n_permutations, local_p_values, total_p_value)
    end subroutine compute_p_values

    !> AUTHOR_FRANZ_ERIC_SILL
    !| (no input validation) Once the permutation tests are calculated ([[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]),
    !| this subroutine calculates the p values for the contributions, i.e. how many of the permutation contributions were at least as high as the real contributions.
    pure subroutine compute_p_values_helper(local_contributions_observed, total_contribution_observed, local_contributions_perm_test, total_contributions_perm_test, n_timepoints, n_permutations, local_p_values, total_p_value)
        integer(int32), intent(in) :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: n_permutations
            !! number of permutations to perform
        real(real64), dimension(n_timepoints), intent(in) :: local_contributions_observed
            !! Per-timepoint contributions for the observed factor-dependent-sample combination
        real(real64), intent(in) :: total_contribution_observed
            !! Total contribution (`sum(local_contributions)`) for the observed factor-dependent-sample combination
        real(real64), dimension(n_timepoints, n_permutations), intent(in) :: local_contributions_perm_test
            !! Per-timepoint contributions for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_permutations), intent(in) :: total_contributions_perm_test
            !! Total contribution (`sum(local_contributions)`) for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis(module):perform_permutation_test(subroutine)]]
        real(real64), dimension(n_timepoints), intent(out) :: local_p_values
            !! calculated p values for local contributions, like: `(local_contributions_perm_test >= local_contributions_observed)/n_permutations`
        real(real64), intent(out) :: total_p_value
            !! calculated p values for total contributions, like: `(total_contributions_perm_test >= total_contribution_observed)/n_permutations`

        integer(int32) :: i_perm, i_timepoint

        total_p_value = 0.0_real64
        local_p_values = 0.0_real64

        do i_perm = 1, n_permutations
            if (total_contributions_perm_test(i_perm) >= total_contribution_observed) then
                total_p_value = total_p_value + 1.0_real64
            end if

            do concurrent (i_timepoint = 1:n_timepoints) shared(local_contributions_perm_test, i_perm, local_contributions_observed, local_p_values)
                if (local_contributions_perm_test(i_timepoint, i_perm) >= local_contributions_observed(i_timepoint)) then
                    local_p_values(i_timepoint) = local_p_values(i_timepoint) + 1.0_real64
                end if
            end do
        end do

        do concurrent (i_timepoint = 1:n_timepoints) shared(local_p_values, n_permutations)
            local_p_values(i_timepoint) = anint(local_p_values(i_timepoint), kind=real64)/real(n_permutations, kind=real64)
        end do

        total_p_value = anint(total_p_value, kind=real64)/real(n_permutations, kind=real64)
    end subroutine compute_p_values_helper

    !> M_EXPORT_C
    !| summary: Contribution analysis for a specific factor-dependent pair
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine compute_contributions(factor, dependent, n_dims, baseline_mode, local_contributions, total_contribution, ierr)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `factor` and `dependent`
        real(real64), dimension(n_dims), intent(in) :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_dims), intent(in) :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        real(real64), dimension(n_dims), intent(out) :: local_contributions
            !! Per-element contributions
        real(real64), intent(out) :: total_contribution
            !! Total contribution (`sum(local_contributions)`)
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dims, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(factor, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(dependent, n_dims, ierr, arg_pos=2_int32)
        CM_VALIDATE_MODE_BASELINE(arg_pos=4_int32)

        if (is_err(ierr)) return

        call compute_contributions_helper(factor, dependent, n_dims, baseline_mode, local_contributions, total_contribution, ierr)
    end subroutine compute_contributions

    !> AUTHOR_FRANZ_ERIC_SILL
    !| (no input validation) This routine performs contribution analysis for a specific factor–dependent pair, no input validation
    pure subroutine compute_contributions_helper(factor, dependent, n_dims, baseline_mode, local_contributions, total_contribution, ierr)
        integer(int32), intent(in) :: n_dims
            !! Number of elements in `factor` and `dependent`
        real(real64), dimension(n_dims), intent(in) :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_dims), intent(in) :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        real(real64), dimension(n_dims), intent(out) :: local_contributions
            !! Per-element contributions
        real(real64), intent(out) :: total_contribution
            !! Total contribution (`sum(local_contributions)`)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_dim
        real(real64) :: factor_baseline, dependent_baseline

        call compute_baselines_factor_dependent(n_dims, factor, dependent, baseline_mode, factor_baseline, dependent_baseline, ierr)
        if (is_err(ierr)) return

        ! Per-timepoint contribution is the product of both series' deviations from their baseline
        ! (an uncentered/instantaneous covariance term): it is positive when factor and dependent move
        ! together (both above or both below baseline) and negative when they move oppositely.
        ! Summing over all timepoints gives the total co-variation between factor and dependent.
        total_contribution = 0.0_real64
        do concurrent (i_dim = 1:n_dims) shared(local_contributions, factor, factor_baseline, dependent, dependent_baseline) reduce(+:total_contribution)
            local_contributions(i_dim) = (factor(i_dim) - factor_baseline) * (dependent(i_dim) - dependent_baseline)
            total_contribution = total_contribution + local_contributions(i_dim)
        end do
    end subroutine compute_contributions_helper

    !> M_EXPORT_C
    !| summary: Contribution analysis for every selected factor-dependent pair
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine compute_all_contributions(trajectories, n_factors, n_samples, n_timepoints, factor_indices, n_selected_factors, dependent_indices, n_selected_dependents, baseline_mode, local_contributions, total_contributions, tmp_factors, tmp_dependent, ierr)
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
        integer(int32), dimension(n_selected_dependents), intent(in) :: dependent_indices
            !! indices of dependents to compute the contributions for
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
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

        integer(int32) :: i_timepoint, i_dependent, i_factor, i_sel_factor, i_sel_dependent, i_sample

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        call validate_dimension_size(n_selected_factors, ierr)
        call validate_dimension_size(n_selected_dependents, ierr)
        call validate_all_in_range_real(trajectories, size(trajectories, kind=int32), ierr)
        call validate_all_in_range_int(factor_indices, n_selected_factors, ierr, min=1, max=n_factors)
        call validate_all_in_range_int(dependent_indices, n_selected_dependents, ierr, min=1, max=n_factors)
        CM_VALIDATE_MODE_BASELINE(arg_pos=9_int32)

        if (is_err(ierr)) return

        !TODO optimize: work across samples is embarrassingly parallel (each i_sample only reads its own slice of `trajectories` and writes its own slice of the outputs), but this loop is sequential because `tmp_factors`/`tmp_dependent` are single shared work buffers reused across iterations. Parallelizing would require per-iteration private copies of these buffers (e.g. `local()` in a do concurrent), which is worth considering given how aggressively do concurrent is used elsewhere in this file for sample loops.
        do i_sample = 1, n_samples
            ! create factor vectors for current sample
            do i_sel_factor = 1, n_selected_factors
                i_factor = factor_indices(i_sel_factor)
                do i_timepoint = 1, n_timepoints
                    tmp_factors(i_timepoint, i_sel_factor) = trajectories(i_factor, i_sample, i_timepoint)
                end do
            end do

            ! calculate contributions for each factor-dependent combination
            do i_sel_dependent = 1, n_selected_dependents
                ! create dependent vector for current sample
                i_dependent = dependent_indices(i_sel_dependent)
                do i_timepoint = 1, n_timepoints
                    tmp_dependent(i_timepoint) = trajectories(i_dependent, i_sample, i_timepoint)
                end do

                do i_sel_factor = 1, n_selected_factors
                    call compute_contributions_helper(tmp_factors(:, i_sel_factor), tmp_dependent, n_timepoints, baseline_mode, local_contributions(:, i_sel_factor, i_sel_dependent, i_sample), total_contributions(i_sel_factor, i_sel_dependent, i_sample), ierr)
                    if (is_err(ierr)) return
                end do
            end do
        end do
    end subroutine compute_all_contributions

    !> M_EXPORT_C
    !| summary: Compute scalar baselines for a factor and dependent variable time series
    !| AUTHOR_JITU_DABA
    pure subroutine compute_baselines_factor_dependent(n_timepoints, factor, dependent, baseline_mode, &
                                                       factor_baseline, dependent_baseline, ierr)

        integer(int32), intent(in) :: n_timepoints
            !! Number of timepoints in both factor and dependent arrays
        real(real64), dimension(n_timepoints), intent(in)  :: factor
            !! Factor time series, length n_timepoints
        real(real64), dimension(n_timepoints), intent(in)  :: dependent
            !! Dependent variable time series, length n_timepoints
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        real(real64), intent(out) :: factor_baseline
            !! Computed baseline for factor
        real(real64), intent(out) :: dependent_baseline
            !! Computed baseline for dependent variable
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        factor_baseline = 0.0_real64
        dependent_baseline = 0.0_real64

        ! Validate that n_timepoints > 0
        call validate_dimension_size(n_timepoints, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return

        select case (baseline_mode)

        case (BASELINE_RAW)
            ! Raw contributions: no centering
            factor_baseline = 0.0_real64
            dependent_baseline = 0.0_real64

        case (BASELINE_MIN)
            ! Minimum-centered contributions
            factor_baseline = minval(factor)
            dependent_baseline = minval(dependent)

        case (BASELINE_MEAN)
            ! Mean-centered contributions
            factor_baseline = sum(factor)/real(n_timepoints, kind=real64)
            dependent_baseline = sum(dependent)/real(n_timepoints, kind=real64)

        case default
            call set_err(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
            return
        end select

        ! Validate that baselines are finite (non-NaN, non-Inf)
        call validate_in_range_real(factor_baseline, ierr)
        call validate_in_range_real(dependent_baseline, ierr)

    end subroutine compute_baselines_factor_dependent

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Helper to map baseline baseline_mode string ("min", "mean", "raw") to integer constant
    pure subroutine get_baseline_mode(c_mode_str, baseline_mode, ierr)
        use, intrinsic :: iso_c_binding, only: c_char
        use tox_conversions, only: c_char_1d_as_string
        character(len=1, kind=c_char), dimension(4), intent(in) :: c_mode_str
            !! baseline_mode string ("min", "mean", "raw")
        integer(int32), intent(out) :: baseline_mode
            !! integer representation for the baseline_mode passed by `c_mode_str`
        integer(int32), intent(out) :: ierr
            !! Error code

        character(len=:), allocatable :: mode_str_f

        call set_ok(ierr)

        call c_char_1d_as_string(c_mode_str, mode_str_f, ierr)
        if (is_err(ierr)) return

        select case (trim(mode_str_f))
        case ("raw")
            baseline_mode = BASELINE_RAW
        case ("min")
            baseline_mode = BASELINE_MIN
        case ("mean")
            baseline_mode = BASELINE_MEAN
        case default
            call set_err(ierr, ERR_INVALID_INPUT)
        end select
    end subroutine get_baseline_mode

    !> AUTHOR_JITU_DABA
    !| This routine computes velocity trajectory from a single position trajectory, no input validation
    pure subroutine compute_velocity_trajectory_helper(trajectory, velocity, n_timepoints)

        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_timepoints), intent(in)  :: trajectory
            !! input position trajectory
        real(real64), dimension(max(0_int32, n_timepoints - 1)), intent(out) :: velocity
            !! output velocity trajectory

        integer(int32) :: i_vel, n_vel

        n_vel = size(velocity, kind=int32)
        if (n_vel == 0_int32) return

        do concurrent (i_vel = 1:n_vel) shared(velocity, trajectory)
            velocity(i_vel) = trajectory(i_vel + 1) - trajectory(i_vel)
        end do
    end subroutine compute_velocity_trajectory_helper

    !> AUTHOR_JITU_DABA
    !| This routine computes velocity trajectory from the `trajectories(n_factors, n_samples, n_timepoints)` matrix
    pure subroutine compute_factor_velocity_from_trajectories_helper(trajectories, n_factors, n_samples, n_timepoints, factor_idx, sample_idx, velocity)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! input position trajectories
        integer(int32), intent(in)  :: factor_idx
            !! factor index of the trajectory
        integer(int32), intent(in)  :: sample_idx
            !! sample index of the trajectory
        real(real64), dimension(max(0, n_timepoints - 1)), intent(out) :: velocity
            !! Calculated velocity for the selected trajectory

        integer(int32) :: n_vel, i_vel

        n_vel = size(velocity, kind=int32)
        if (n_vel == 0_int32) return

        do concurrent (i_vel = 1:n_vel) shared(velocity, trajectories, factor_idx, sample_idx)
            velocity(i_vel) = trajectories(factor_idx, sample_idx, i_vel + 1) - trajectories(factor_idx, sample_idx, i_vel)
        end do
    end subroutine compute_factor_velocity_from_trajectories_helper

    !> M_EXPORT_C
    !| summary: Compute velocity trajectory from a single position trajectory
    !| AUTHOR_JITU_DABA
    pure subroutine compute_velocity_trajectory(trajectory, velocity, n_timepoints, ierr)

        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(n_timepoints), intent(in)  :: trajectory
            !! input position trajectory
        real(real64), dimension(max(0_int32, n_timepoints - 1)), intent(out) :: velocity
            !! output velocity trajectory

        call set_ok(ierr)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=2)
        call validate_all_in_range_real(trajectory, n_timepoints, ierr, arg_pos=1)
        if (is_err(ierr)) return

        call compute_velocity_trajectory_helper(trajectory, velocity, n_timepoints)
    end subroutine compute_velocity_trajectory

    !> M_EXPORT_C
    !| summary: Compute acceleration trajectory from a single velocity trajectory
    !| AUTHOR_JITU_DABA
    pure subroutine compute_acceleration_from_velocity_trajectory(velocity, acceleration, &
                                                                  n_timepoints, ierr)
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(max(0_int32, n_timepoints - 1)), intent(in)  :: velocity
            !! velocity trajectory
        real(real64), dimension(max(0_int32, n_timepoints - 2)), intent(out) :: acceleration
            !! acceleration trajectory

        call compute_velocity_trajectory(velocity, acceleration, size(velocity, kind=int32), ierr)
    end subroutine compute_acceleration_from_velocity_trajectory

    !> AUTHOR_JITU_DABA
    !| This routine computes velocity trajectories from position trajectories, no input validation
    pure subroutine compute_velocity_trajectories_helper(trajectories, velocity, &
                                                         n_factors, n_samples, n_timepoints)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in)  :: trajectories
            !! input position trajectories
        real(real64), dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples), intent(out) ::  velocity
            !! output velocity trajectories

        integer(int32) :: i_factor, i_sample
        integer(int32) :: i_vel, n_vel

        n_vel = size(velocity, dim=1, kind=int32)
        if (n_vel == 0_int32) return

        do concurrent (i_sample = 1:n_samples) shared(trajectories, n_factors, n_timepoints, velocity)
            do concurrent (i_factor = 1:n_factors) shared(trajectories, n_samples, n_timepoints, i_sample, velocity)
                call compute_factor_velocity_from_trajectories_helper(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, velocity(:, i_factor, i_sample))
            end do
        end do
    end subroutine compute_velocity_trajectories_helper

    !> M_EXPORT_C
    !| summary: Compute velocity trajectories for all factors and samples
    !| AUTHOR_JITU_DABA
    pure subroutine compute_velocity_trajectories(trajectories, velocity, &
                                                  n_factors, n_samples, n_timepoints, ierr)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in)  :: trajectories
            !! input position trajectories
        real(real64), dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples), intent(out) :: velocity
            !! output velocity trajectories

        call set_ok(ierr)
        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(trajectories, size(trajectories, kind=int32), ierr, arg_pos=1_int32)
        if (is_err(ierr)) return

        call compute_velocity_trajectories_helper(trajectories, velocity, n_factors, n_samples, n_timepoints)
    end subroutine compute_velocity_trajectories

    !> AUTHOR_JITU_DABA
    !| This routine computes acceleration trajectories from velocity trajectories, no input validation
    pure subroutine compute_acceleration_from_velocity_helper(velocity, acceleration, &
                                                              n_factors, n_samples, n_timepoints)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        real(real64), dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples), intent(in)  ::  velocity
            !! input velocity trajectories
        real(real64), dimension(max(0_int32, n_timepoints - 2), n_factors, n_samples), intent(out) :: acceleration
            !! output acceleration trajectories

        integer(int32) :: i_factor, i_sample, n_vel

        n_vel = size(velocity, dim=1, kind=int32)
        if (n_vel < 2) return

        do concurrent (i_sample = 1:n_samples) shared(velocity, acceleration, n_factors, n_vel)
            do concurrent (i_factor = 1:n_factors) shared(velocity, acceleration, n_vel, i_sample)
                call compute_velocity_trajectory_helper(velocity(:, i_factor, i_sample), &
                                                        acceleration(:, i_factor, i_sample), &
                                                        n_vel)
            end do
        end do
    end subroutine compute_acceleration_from_velocity_helper

    !> M_EXPORT_C
    !| summary: Compute acceleration trajectories from velocity trajectories
    !| AUTHOR_JITU_DABA
    pure subroutine compute_acceleration_from_velocity(velocity, acceleration, &
                                                       n_factors, n_samples, n_timepoints, ierr)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        integer(int32), intent(out) :: ierr
            !! Error code
        real(real64), dimension(max(0_int32, n_timepoints - 1), n_factors, n_samples), intent(in)  :: velocity
            !! input velocity trajectories
        real(real64), dimension(max(0_int32, n_timepoints - 2), n_factors, n_samples), intent(out) :: acceleration
            !! output acceleration trajectories

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(velocity, size(velocity, kind=int32), ierr, arg_pos=1_int32)
        if (is_err(ierr)) return

        call compute_acceleration_from_velocity_helper(velocity, acceleration, n_factors, n_samples, n_timepoints)
    end subroutine compute_acceleration_from_velocity

    !> M_EXPORT_C
    !| summary: Compute velocity and acceleration contributions for all variable pairs (expert entry point)
    !| AUTHOR_JITU_DABA
    !| @note
    !| Performance layout:
    !|
    !| `trajectories` uses `(n_factors, n_samples, n_timepoints)`.
    !| Velocity and acceleration use time-first layouts:
    !|
    !|    - `velocity`     -> `(max(0, n_timepoints-1), n_factors, n_samples)`
    !|    - `acceleration` -> `(max(0, n_timepoints-2), n_factors, n_samples)`
    !|
    !| This keeps slices like `velocity(:, factor, sample)` contiguous,
    !| avoids expensive tmporaries, and improves cache efficiency.
    !| @endnote
    pure subroutine compute_velocity_acceleration_contributions(trajectories, n_factors, n_samples, n_timepoints, baseline_mode, &
                                                                tmp_factors, tmp_dependent, tmp_contributions, &
                                                                contrib_velocity, velocity_contribution_series, &
                                                                contrib_acceleration, acceleration_contribution_series, ierr)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! input position trajectories

        ! Workspace (preallocated by caller, reused for both velocity and acceleration)
        real(real64), dimension(n_timepoints - 1, n_factors), intent(out) :: tmp_factors
            !! workspace for factor data (used for velocity and acceleration)
        real(real64), dimension(n_timepoints - 1), intent(out) :: tmp_dependent
            !! workspace for dependent data (used for velocity and acceleration)
        real(real64), dimension(n_timepoints - 1), intent(out) :: tmp_contributions
            !! workspace for contribution calculations (used for velocity and acceleration)

        ! Outputs
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

        integer(int32) :: i_sample, i_factor, i_dependent, i_timepoint, tmp_ierr
        integer(int32) :: n_vel, n_acc

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_samples, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_timepoints, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(trajectories, size(trajectories, kind=int32), ierr, arg_pos=1_int32)
        CM_VALIDATE_MODE_BASELINE(arg_pos=5_int32)
        if (is_err(ierr)) return

        contrib_velocity = 0.0_real64
        velocity_contribution_series = 0.0_real64
        contrib_acceleration = 0.0_real64
        acceleration_contribution_series = 0.0_real64

        if (n_timepoints <= 1) return

        n_vel = n_timepoints - 1_int32
        n_acc = n_timepoints - 2_int32

        do i_sample = 1, n_samples
            ! ---- Step 1: velocity contributions ----
            ! GFORTRAN BUG: do concurrent (i_factor = 1:n_factors) shared(trajectories, n_samples, n_timepoints, i_sample, tmp_factors)
            do i_factor = 1, n_factors
                call compute_factor_velocity_from_trajectories_helper(trajectories, n_factors, n_samples, n_timepoints, i_factor, i_sample, tmp_factors(:, i_factor))
            end do

            do i_dependent = 1, n_factors
                call compute_factor_velocity_from_trajectories_helper(trajectories, n_factors, n_samples, n_timepoints, i_dependent, i_sample, tmp_dependent)

                do i_factor = 1, n_factors
                    ! velocity has one fewer timepoint than the original series (finite differences
                    ! consume one endpoint), so its contributions are written starting at output index 2
                    ! to stay aligned with the original time axis; index 1 has no velocity and is left at zero.
                    velocity_contribution_series(1, i_factor, i_dependent, i_sample) = 0.0_real64

                    call compute_contributions_helper( &
                        tmp_factors(:, i_factor), tmp_dependent, n_vel, baseline_mode, &
                        tmp_contributions, contrib_velocity(i_factor, i_dependent, i_sample), tmp_ierr)
                    if (is_err(tmp_ierr)) ierr = tmp_ierr

                    do concurrent (i_timepoint = 2:n_timepoints) shared(velocity_contribution_series, i_factor, i_dependent, i_sample, tmp_contributions)
                        velocity_contribution_series(i_timepoint, i_factor, i_dependent, i_sample) = &
                            tmp_contributions(i_timepoint - 1)
                    end do
                end do
            end do

            ! ---- Step 2: acceleration contributions (reuse same workspace) ----

            if (n_acc <= 0) cycle

            do i_factor = 1, n_factors
                call compute_velocity_trajectory_helper(tmp_factors(:, i_factor), tmp_contributions, n_vel)
                tmp_factors(1:n_acc, i_factor) = tmp_contributions(1:n_acc)
            end do

            do i_dependent = 1, n_factors
                call compute_factor_velocity_from_trajectories_helper(trajectories, n_factors, n_samples, n_timepoints, i_dependent, i_sample, tmp_contributions)
                call compute_velocity_trajectory_helper(tmp_contributions, tmp_dependent, n_vel)

                do i_factor = 1, n_factors
                    ! acceleration has two fewer timepoints than the original series (a second finite
                    ! difference), so its contributions are written starting at output index 3; the first
                    ! two indices have no acceleration and are left at zero.
                    acceleration_contribution_series(1, i_factor, i_dependent, i_sample) = 0.0_real64
                    acceleration_contribution_series(2, i_factor, i_dependent, i_sample) = 0.0_real64

                    call compute_contributions_helper( &
                        tmp_factors(1:n_acc, i_factor), tmp_dependent(1:n_acc), n_acc, baseline_mode, &
                        tmp_contributions(1:n_acc), contrib_acceleration(i_factor, i_dependent, i_sample), tmp_ierr)
                    if (is_err(tmp_ierr)) ierr = tmp_ierr

                    do concurrent (i_timepoint = 3:n_timepoints) shared(acceleration_contribution_series, i_factor, i_dependent, i_sample, tmp_contributions)
                        acceleration_contribution_series(i_timepoint, i_factor, i_dependent, i_sample) = &
                            tmp_contributions(i_timepoint - 2)
                    end do
                end do
            end do
        end do
    end subroutine compute_velocity_acceleration_contributions

    !> M_EXPORT_C
    !| summary: Compute velocity and acceleration contributions for all variable pairs (allocating entry point)
    !| AUTHOR_JITU_DABA
    !| @note
    !| Performance layout:
    !|
    !| `trajectories` uses `(n_factors, n_samples, n_timepoints)`.
    !| Velocity and acceleration use time-first layouts:
    !|
    !|    - `velocity`     -> `(max(0, n_timepoints-1), n_factors, n_samples)`
    !|    - `acceleration` -> `(max(0, n_timepoints-2), n_factors, n_samples)`
    !|
    !| This keeps slices like `velocity(:, factor, sample)` contiguous,
    !| avoids expensive tmporaries, and improves cache efficiency.
    !| @endnote
    pure subroutine compute_velocity_acceleration_contributions_alloc(trajectories, n_factors, n_samples, n_timepoints, baseline_mode, &
                                                                      contrib_velocity, velocity_contribution_series, &
                                                                      contrib_acceleration, acceleration_contribution_series, ierr)

        integer(int32), intent(in)  :: n_factors
            !! number of factors
        integer(int32), intent(in)  :: n_samples
            !! number of samples
        integer(int32), intent(in)  :: n_timepoints
            !! number of timepoints
        integer(int32), intent(in) :: baseline_mode
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Raw/zero baseline | [[tox_trajectory_contribution_analysis(module):MODE_RAW(variable)]] |
            !! | arithmetic mean | [[tox_trajectory_contribution_analysis(module):MODE_MEAN(variable)]] |
            !! | minimum value | [[tox_trajectory_contribution_analysis(module):MODE_MIN(variable)]] |
            !!
        real(real64), dimension(n_factors, n_samples, n_timepoints), intent(in) :: trajectories
            !! input position trajectories

        real(real64), dimension(n_factors, n_factors, n_samples), intent(out) :: contrib_velocity
            !! output velocity contributions
        real(real64), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out) :: velocity_contribution_series
            !! output acceleration contributions
        real(real64), dimension(n_factors, n_factors, n_samples), intent(out) :: contrib_acceleration
            !! output acceleration contributions
        real(real64), dimension(n_timepoints, n_factors, n_factors, n_samples), intent(out) :: acceleration_contribution_series
            !! output acceleration contributions
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Workspace (allocated once here)
        real(real64), allocatable :: tmp_factors(:, :), tmp_dependent(:), tmp_contributions(:)
            !! workspace arrays (reused for velocity and acceleration)

        call set_ok(ierr)

        call validate_dimension_size(n_factors, ierr)
        call validate_dimension_size(n_samples, ierr)
        call validate_dimension_size(n_timepoints, ierr)
        if (is_err(ierr)) return

        ! Allocate reusable work vectors once
        if (n_timepoints > 1) then
            M_ALLOCATE(tmp_factors(n_timepoints - 1, n_factors))
            M_ALLOCATE(tmp_dependent(n_timepoints - 1))
            M_ALLOCATE(tmp_contributions(n_timepoints - 1))
        end if

        call compute_velocity_acceleration_contributions(trajectories, n_factors, n_samples, n_timepoints, baseline_mode, &
                                                         tmp_factors, tmp_dependent, tmp_contributions, &
                                                         contrib_velocity, velocity_contribution_series, &
                                                         contrib_acceleration, acceleration_contribution_series, ierr)
    end subroutine compute_velocity_acceleration_contributions_alloc
end module tox_trajectory_contribution_analysis
