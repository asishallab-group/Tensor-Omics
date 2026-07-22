# Generated. Do not edit.

#' Normalize a single variable across time using min-max scaling
#'
#' @param v a numeric vector. Original time series
#' @return a named list with elements `v_norm`, `status`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_normalization::normalize_variable_timeseries}.
#' @export
normalize_variable_timeseries <- function(v) {
    v <- .tox_as_double_vector(v, "v")
    .result <- .normalize_variable_timeseries_rcpp(v)
    .arguments <- c("v", "v_norm", "n_points", "ierr", "status")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        v_norm = .result$v_norm,
        status = .result$status
    )
}

#' Normalize all factors in a single trajectory independently across time
#'
#' @param trajectory a numeric matrix. Original trajectory for one sample
#' @return a named list with elements `trajectory_norm`, `status`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_normalization::normalize_single_trajectory}.
#' @export
normalize_single_trajectory <- function(trajectory) {
    trajectory <- .tox_as_double_matrix(trajectory, "trajectory")
    .result <- .normalize_single_trajectory_rcpp(trajectory)
    .arguments <- c("trajectory", "trajectory_norm", "n_factors", "n_timepoints", "ierr", "status")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        trajectory_norm = .result$trajectory_norm,
        status = .result$status
    )
}

#' Normalize all trajectories across multiple entities
#'
#' Normalizes each factor independently across time for each sample
#'
#' @param trajectories a numeric array of rank 3. Original trajectories
#' @return a named list with elements `trajectories_norm`, `status`.
#'
#' Generated from the Fortran procedure \code{tox_trajectory_normalization::normalize_all_trajectories_alloc}.
#' @export
normalize_all_trajectories <- function(trajectories) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    .result <- .normalize_all_trajectories_rcpp(trajectories)
    .arguments <- c("trajectories", "trajectories_norm", "n_factors", "n_samples", "n_timepoints", "ierr", "status")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        trajectories_norm = .result$trajectories_norm,
        status = .result$status
    )
}
