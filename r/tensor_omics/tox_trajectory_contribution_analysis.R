# Generated. Do not edit.

#' For a factor-dependent pair, calculates the contributions against the same dependent taken from a random different sample
#'
#' @param trajectories a numeric array of rank 3. expression vectors across different samples over time
#' @param factor_idx a integer scalar. index of factor to compute the permutation contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param dependent_idx a integer scalar. index of dependent to compute the permutation contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param sample_idx a integer scalar. index of sample to compute the permutation contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_samples`.
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @param n_permutations a integer scalar. number of permutations to perform
#' @param random_seed a integer scalar. Seed to use for random number generation.
#' @return a named list with elements `local_contributions`, `total_contributions`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::perform_permutation_test}.
#' @export
perform_permutation_test_expert <- function(trajectories, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, random_seed = NULL) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    factor_idx <- .tox_as_integer_scalar(factor_idx, "factor_idx")
    dependent_idx <- .tox_as_integer_scalar(dependent_idx, "dependent_idx")
    sample_idx <- .tox_as_integer_scalar(sample_idx, "sample_idx")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    if (!is.null(random_seed))
        random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    .result <- .Call("perform_permutation_test_expert_call", trajectories, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, random_seed)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_idx", "dependent_idx", "sample_idx", "baseline_mode", "n_permutations", "local_contributions", "total_contributions", "tmp_factor", "tmp_dependent", "random_seed", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        local_contributions = .result$local_contributions,
        total_contributions = .result$total_contributions
    )
}

#' For a factor-dependent pair, calculates the contributions against the same dependent taken from a random different sample
#'
#' @param trajectories a numeric array of rank 3. expression vectors across different samples over time
#' @param factor_idx a integer scalar. index of factor to compute the permutation contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param dependent_idx a integer scalar. index of dependent to compute the permutation contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param sample_idx a integer scalar. index of sample to compute the permutation contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_samples`.
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @param n_permutations a integer scalar. number of permutations to perform
#' @param random_seed a integer scalar. Seed to use for random number generation.
#' @return a named list with elements `local_contributions`, `total_contributions`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::perform_permutation_test_alloc}.
#' @export
perform_permutation_test <- function(trajectories, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, random_seed = NULL) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    factor_idx <- .tox_as_integer_scalar(factor_idx, "factor_idx")
    dependent_idx <- .tox_as_integer_scalar(dependent_idx, "dependent_idx")
    sample_idx <- .tox_as_integer_scalar(sample_idx, "sample_idx")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    if (!is.null(random_seed))
        random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    .result <- .Call("perform_permutation_test_call", trajectories, factor_idx, dependent_idx, sample_idx, baseline_mode, n_permutations, random_seed)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_idx", "dependent_idx", "sample_idx", "baseline_mode", "n_permutations", "local_contributions", "total_contributions", "random_seed", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        local_contributions = .result$local_contributions,
        total_contributions = .result$total_contributions
    )
}

#' Calculates the p values for the contributions once the permutation tests are done
#'
#' Given the permutation tests ([[tox_trajectory_contribution_analysis_kernel(module):perform_permutation_test_kernel(subroutine)]]),
#' this counts how many of the permutation contributions were at least as high as the real ones.
#'
#' @param local_contributions_observed a numeric vector. Per-timepoint contributions for the observed factor-dependent-sample combination
#' @param total_contribution_observed a numeric scalar. Total contribution (`sum(local_contributions)`) for the observed factor-dependent-sample combination
#' @param local_contributions_perm_test a numeric matrix. Per-timepoint contributions for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis_kernel(module):perform_permutation_test(subroutine)]]
#' @param total_contributions_perm_test a numeric vector. Total contribution (`sum(local_contributions)`) for the factor-dependent-random_sample combinations from [[tox_trajectory_contribution_analysis_kernel(module):perform_permutation_test(subroutine)]]
#' @return a named list with elements `local_p_values`, `total_p_value`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_p_values}.
#' @export
compute_p_values <- function(local_contributions_observed, total_contribution_observed, local_contributions_perm_test, total_contributions_perm_test) {
    local_contributions_observed <- .tox_as_double_vector(local_contributions_observed, "local_contributions_observed")
    total_contribution_observed <- .tox_as_double_scalar(total_contribution_observed, "total_contribution_observed")
    local_contributions_perm_test <- .tox_as_double_matrix(local_contributions_perm_test, "local_contributions_perm_test")
    total_contributions_perm_test <- .tox_as_double_vector(total_contributions_perm_test, "total_contributions_perm_test")
    if (dim(local_contributions_perm_test)[1] != length(local_contributions_observed))
        .tox_shape_error("local_contributions_perm_test", dim(local_contributions_perm_test)[1], "local_contributions_observed", length(local_contributions_observed))
    if (length(total_contributions_perm_test) != dim(local_contributions_perm_test)[2])
        .tox_shape_error("total_contributions_perm_test", length(total_contributions_perm_test), "local_contributions_perm_test", dim(local_contributions_perm_test)[2])

    .result <- .Call("compute_p_values_call", local_contributions_observed, total_contribution_observed, local_contributions_perm_test, total_contributions_perm_test)
    .arguments <- c("local_contributions_observed", "total_contribution_observed", "local_contributions_perm_test", "total_contributions_perm_test", "n_timepoints", "n_permutations", "local_p_values", "total_p_value", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        local_p_values = .result$local_p_values,
        total_p_value = .result$total_p_value
    )
}

#' Performs contribution analysis for a specific factor-dependent pair
#'
#' @param factor a numeric vector. Factor time series, length n_timepoints
#' @param dependent a numeric vector. Dependent variable time series, length n_timepoints
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `local_contributions`, `total_contribution`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_contributions}.
#' @export
compute_contributions <- function(factor, dependent, baseline_mode) {
    factor <- .tox_as_double_vector(factor, "factor")
    dependent <- .tox_as_double_vector(dependent, "dependent")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    if (length(dependent) != length(factor))
        .tox_shape_error("dependent", length(dependent), "factor", length(factor))

    .result <- .Call("compute_contributions_call", factor, dependent, baseline_mode)
    .arguments <- c("factor", "dependent", "n_dims", "baseline_mode", "local_contributions", "total_contribution", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        local_contributions = .result$local_contributions,
        total_contribution = .result$total_contribution
    )
}

#' Contribution analysis for every selected factor-dependent pair
#'
#' @param trajectories a numeric array of rank 3. expression vectors across different samples over time
#' @param factor_indices a integer vector. indices of factors to compute the contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param dependent_indices a integer vector. indices of dependents to compute the contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `local_contributions`, `total_contributions`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_all_contributions}.
#' @export
compute_all_contributions_expert <- function(trajectories, factor_indices, dependent_indices, baseline_mode) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    factor_indices <- .tox_as_integer_vector(factor_indices, "factor_indices")
    dependent_indices <- .tox_as_integer_vector(dependent_indices, "dependent_indices")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    .result <- .Call("compute_all_contributions_expert_call", trajectories, factor_indices, dependent_indices, baseline_mode)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_indices", "n_selected_factors", "dependent_indices", "n_selected_dependents", "baseline_mode", "local_contributions", "total_contributions", "tmp_factors", "tmp_dependent", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        local_contributions = .result$local_contributions,
        total_contributions = .result$total_contributions
    )
}

#' Contribution analysis for every selected factor-dependent pair
#'
#' @param trajectories a numeric array of rank 3. expression vectors across different samples over time
#' @param factor_indices a integer vector. indices of factors to compute the contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param dependent_indices a integer vector. indices of dependents to compute the contributions for
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_factors`.
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `local_contributions`, `total_contributions`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_all_contributions_alloc}.
#' @export
compute_all_contributions <- function(trajectories, factor_indices, dependent_indices, baseline_mode) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    factor_indices <- .tox_as_integer_vector(factor_indices, "factor_indices")
    dependent_indices <- .tox_as_integer_vector(dependent_indices, "dependent_indices")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    .result <- .Call("compute_all_contributions_call", trajectories, factor_indices, dependent_indices, baseline_mode)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_indices", "n_selected_factors", "dependent_indices", "n_selected_dependents", "baseline_mode", "local_contributions", "total_contributions", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        local_contributions = .result$local_contributions,
        total_contributions = .result$total_contributions
    )
}

#' Compute scalar baselines for a factor and dependent variable time series
#'
#' @param factor a numeric vector. Factor time series, length n_timepoints
#' @param dependent a numeric vector. Dependent variable time series, length n_timepoints
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `factor_baseline`, `dependent_baseline`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_baselines_factor_dependent}.
#' @export
compute_baselines_factor_dependent <- function(factor, dependent, baseline_mode) {
    factor <- .tox_as_double_vector(factor, "factor")
    dependent <- .tox_as_double_vector(dependent, "dependent")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    if (length(dependent) != length(factor))
        .tox_shape_error("dependent", length(dependent), "factor", length(factor))

    .result <- .Call("compute_baselines_factor_dependent_call", factor, dependent, baseline_mode)
    .arguments <- c("n_timepoints", "factor", "dependent", "baseline_mode", "factor_baseline", "dependent_baseline", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        factor_baseline = .result$factor_baseline,
        dependent_baseline = .result$dependent_baseline
    )
}

#' Compute velocity trajectory from a single position trajectory
#'
#' @param trajectory a numeric vector. input position trajectory
#' @return output velocity trajectory
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_velocity_trajectory}.
#' @export
compute_velocity_trajectory <- function(trajectory) {
    trajectory <- .tox_as_double_vector(trajectory, "trajectory")
    .result <- .Call("compute_velocity_trajectory_call", trajectory)
    .arguments <- c("trajectory", "velocity", "n_timepoints", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$velocity
}

#' Compute acceleration trajectory from a single velocity trajectory
#'
#' @param velocity a numeric vector. velocity trajectory
#' @param n_timepoints a integer scalar. number of timepoints
#' @return acceleration trajectory
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_acceleration_from_velocity_trajectory}.
#' @export
compute_acceleration_from_velocity_trajectory <- function(velocity, n_timepoints) {
    velocity <- .tox_as_double_vector(velocity, "velocity")
    n_timepoints <- .tox_as_integer_scalar(n_timepoints, "n_timepoints")
    .result <- .Call("compute_acceleration_from_velocity_trajectory_call", velocity, n_timepoints)
    .arguments <- c("velocity", "acceleration", "n_timepoints", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$acceleration
}

#' Computes velocity trajectories from position trajectories
#'
#' @param trajectories a numeric array of rank 3. input position trajectories
#' @return output velocity trajectories
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_velocity_trajectories}.
#' @export
compute_velocity_trajectories <- function(trajectories) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    .result <- .Call("compute_velocity_trajectories_call", trajectories)
    .arguments <- c("trajectories", "velocity", "n_factors", "n_samples", "n_timepoints", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$velocity
}

#' Computes acceleration trajectories from velocity trajectories
#'
#' @param velocity a numeric array of rank 3. input velocity trajectories
#' @param n_timepoints a integer scalar. number of timepoints
#' @return output acceleration trajectories
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_acceleration_from_velocity}.
#' @export
compute_acceleration_from_velocity <- function(velocity, n_timepoints) {
    velocity <- .tox_as_double_array(velocity, "velocity", 3L)
    n_timepoints <- .tox_as_integer_scalar(n_timepoints, "n_timepoints")
    .result <- .Call("compute_acceleration_from_velocity_call", velocity, n_timepoints)
    .arguments <- c("velocity", "acceleration", "n_factors", "n_samples", "n_timepoints", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$acceleration
}

#' Compute velocity and acceleration contributions for all variable pairs (expert entry point)
#'
#' @note
#' Performance layout:
#'
#' `trajectories` uses `(n_factors, n_samples, n_timepoints)`.
#' Velocity and acceleration use time-first layouts:
#'
#' - `velocity`     -> `(max(0, n_timepoints-1), n_factors, n_samples)`
#' - `acceleration` -> `(max(0, n_timepoints-2), n_factors, n_samples)`
#'
#' This keeps slices like `velocity(:, factor, sample)` contiguous,
#' avoids expensive tmporaries, and improves cache efficiency.
#' @endnote
#'
#' @param trajectories a numeric array of rank 3. input position trajectories
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `contrib_velocity`, `velocity_contribution_series`, `contrib_acceleration`, `acceleration_contribution_series`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_velocity_acceleration_contributions}.
#' @export
compute_velocity_acceleration_contributions_expert <- function(trajectories, baseline_mode) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    .result <- .Call("compute_velocity_acceleration_contributions_expert_call", trajectories, baseline_mode)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "baseline_mode", "tmp_factors", "tmp_dependent", "tmp_contributions", "contrib_velocity", "velocity_contribution_series", "contrib_acceleration", "acceleration_contribution_series", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        contrib_velocity = .result$contrib_velocity,
        velocity_contribution_series = .result$velocity_contribution_series,
        contrib_acceleration = .result$contrib_acceleration,
        acceleration_contribution_series = .result$acceleration_contribution_series
    )
}

#' Compute velocity and acceleration contributions for all variable pairs (expert entry point)
#'
#' @note
#' Performance layout:
#'
#' `trajectories` uses `(n_factors, n_samples, n_timepoints)`.
#' Velocity and acceleration use time-first layouts:
#'
#' - `velocity`     -> `(max(0, n_timepoints-1), n_factors, n_samples)`
#' - `acceleration` -> `(max(0, n_timepoints-2), n_factors, n_samples)`
#'
#' This keeps slices like `velocity(:, factor, sample)` contiguous,
#' avoids expensive tmporaries, and improves cache efficiency.
#' @endnote
#'
#' @param trajectories a numeric array of rank 3. input position trajectories
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `contrib_velocity`, `velocity_contribution_series`, `contrib_acceleration`, `acceleration_contribution_series`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_contribution_analysis::compute_velocity_acceleration_contributions_alloc}.
#' @export
compute_velocity_acceleration_contributions <- function(trajectories, baseline_mode) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    .result <- .Call("compute_velocity_acceleration_contributions_call", trajectories, baseline_mode)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "baseline_mode", "contrib_velocity", "velocity_contribution_series", "contrib_acceleration", "acceleration_contribution_series", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        contrib_velocity = .result$contrib_velocity,
        velocity_contribution_series = .result$velocity_contribution_series,
        contrib_acceleration = .result$contrib_acceleration,
        acceleration_contribution_series = .result$acceleration_contribution_series
    )
}
