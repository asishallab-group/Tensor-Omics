# Generated. Do not edit.

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
#' Generated from the Fortran module \code{tox_trajectory_contribution_analysis_kernel}.
#' @export
compute_all_contributions_kernel <- function(trajectories, factor_indices, dependent_indices, baseline_mode) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    factor_indices <- .tox_as_integer_vector(factor_indices, "factor_indices")
    dependent_indices <- .tox_as_integer_vector(dependent_indices, "dependent_indices")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    .result <- .Call("compute_all_contributions_kernel_call", trajectories, factor_indices, dependent_indices, baseline_mode)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "factor_indices", "n_selected_factors", "dependent_indices", "n_selected_dependents", "baseline_mode", "local_contributions", "total_contributions", "tmp_factors", "tmp_dependent", "ierr")
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
#' Generated from the Fortran module \code{tox_trajectory_contribution_analysis_kernel}.
#' @export
compute_baselines_factor_dependent_kernel <- function(factor, dependent, baseline_mode) {
    factor <- .tox_as_double_vector(factor, "factor")
    dependent <- .tox_as_double_vector(dependent, "dependent")
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    if (length(dependent) != length(factor))
        .tox_shape_error("dependent", length(dependent), "factor", length(factor))

    .result <- .Call("compute_baselines_factor_dependent_kernel_call", factor, dependent, baseline_mode)
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
#' Generated from the Fortran module \code{tox_trajectory_contribution_analysis_kernel}.
#' @export
compute_velocity_trajectory_kernel <- function(trajectory) {
    trajectory <- .tox_as_double_vector(trajectory, "trajectory")
    .result <- .Call("compute_velocity_trajectory_kernel_call", trajectory)
    .arguments <- c("trajectory", "velocity", "n_timepoints")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$velocity
}

#' Compute acceleration trajectory from a single velocity trajectory
#'
#' @param velocity a numeric vector. velocity trajectory
#' @param n_timepoints a integer scalar. number of timepoints
#' @return acceleration trajectory
#'
#' Generated from the Fortran module \code{tox_trajectory_contribution_analysis_kernel}.
#' @export
compute_acceleration_from_velocity_trajectory_kernel <- function(velocity, n_timepoints) {
    velocity <- .tox_as_double_vector(velocity, "velocity")
    n_timepoints <- .tox_as_integer_scalar(n_timepoints, "n_timepoints")
    .result <- .Call("compute_acceleration_from_velocity_trajectory_kernel_call", velocity, n_timepoints)
    .arguments <- c("velocity", "acceleration", "n_timepoints")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$acceleration
}

#' Compute velocity and acceleration contributions for all variable pairs
#'
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
#'
#' @param trajectories a numeric array of rank 3. input position trajectories
#' @param baseline_mode a string, one of "raw", "mean", "min"
#' @return a named list with elements `contrib_velocity`, `velocity_contribution_series`, `contrib_acceleration`, `acceleration_contribution_series`.
#'
#' Generated from the Fortran module \code{tox_trajectory_contribution_analysis_kernel}.
#' @export
compute_velocity_acceleration_contributions_kernel <- function(trajectories, baseline_mode) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    baseline_mode <- .tox_as_mode(baseline_mode, "baseline_mode", c("raw", "mean", "min"))
    .result <- .Call("compute_velocity_acceleration_contributions_kernel_call", trajectories, baseline_mode)
    .arguments <- c("trajectories", "n_factors", "n_samples", "n_timepoints", "baseline_mode", "tmp_factors", "tmp_dependent", "tmp_contributions", "contrib_velocity", "velocity_contribution_series", "contrib_acceleration", "acceleration_contribution_series", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        contrib_velocity = .result$contrib_velocity,
        velocity_contribution_series = .result$velocity_contribution_series,
        contrib_acceleration = .result$contrib_acceleration,
        acceleration_contribution_series = .result$acceleration_contribution_series
    )
}
