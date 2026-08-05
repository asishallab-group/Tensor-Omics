# Generated. Do not edit.

#' Normalize a single variable across time using min-max scaling
#'
#' Generated from the Fortran module \code{tox_trajectory_normalization}.
#'
#' @param v a numeric vector. Original time series
#' @return a named list with elements:
#'   \item{v_norm}{a numeric vector. Normalized time series}
#'   \item{status}{a integer scalar. Status code for specific warnings}
#' @export
normalize_variable_timeseries <- function(v) {
    v <- .tox_as_double_vector(v, "v")
    .result <- .Call("normalize_variable_timeseries_call", v)
    .arguments <- c("v", "v_norm", "n_points", "status", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        v_norm = .result$v_norm,
        status = .result$status
    )
}

#' Normalize all factors in a single trajectory independently across time
#'
#' Generated from the Fortran module \code{tox_trajectory_normalization}.
#'
#' @param trajectory a numeric matrix. Original trajectory for one sample
#' @return a named list with elements:
#'   \item{trajectory_norm}{a numeric matrix. Normalized trajectory for one sample}
#'   \item{status}{a integer vector. Status code for specific warnings, one per factor}
#' @export
normalize_single_trajectory <- function(trajectory) {
    trajectory <- .tox_as_double_matrix(trajectory, "trajectory")
    .result <- .Call("normalize_single_trajectory_call", trajectory)
    .arguments <- c("trajectory", "trajectory_norm", "n_factors", "n_timepoints", "status", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        trajectory_norm = .result$trajectory_norm,
        status = .result$status
    )
}

#' Normalize all trajectories across multiple entities
#'
#' independently across time for each sample.
#'
#' Generated from the Fortran module \code{tox_trajectory_normalization}.
#'
#' @param trajectories a numeric array of rank 3. Original trajectories
#' @return a named list with elements:
#'   \item{trajectories_norm}{a numeric array of rank 3. Normalized trajectories}
#'   \item{status}{a integer matrix. Status code for specific warnings, one per factor per sample}
#' @export
normalize_all_trajectories_expert <- function(trajectories) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    .result <- .Call("normalize_all_trajectories_expert_call", trajectories)
    .arguments <- c("trajectories", "trajectories_norm", "n_factors", "n_samples", "n_timepoints", "tmp_series", "tmp_series_norm", "status", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        trajectories_norm = .result$trajectories_norm,
        status = .result$status
    )
}

#' Normalize all trajectories across multiple entities
#'
#' independently across time for each sample.
#'
#' Generated from the Fortran module \code{tox_trajectory_normalization}.
#'
#' @param trajectories a numeric array of rank 3. Original trajectories
#' @return a named list with elements:
#'   \item{trajectories_norm}{a numeric array of rank 3. Normalized trajectories}
#'   \item{status}{a integer matrix. Status code for specific warnings, one per factor per sample}
#' @export
normalize_all_trajectories <- function(trajectories) {
    trajectories <- .tox_as_double_array(trajectories, "trajectories", 3L)
    .result <- .Call("normalize_all_trajectories_call", trajectories)
    .arguments <- c("trajectories", "trajectories_norm", "n_factors", "n_samples", "n_timepoints", "status", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        trajectories_norm = .result$trajectories_norm,
        status = .result$status
    )
}
