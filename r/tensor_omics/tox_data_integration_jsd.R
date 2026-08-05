# Generated. Do not edit.

#' Compute the shared residual range [-R, R] from a pooled set of absolute residuals
#'
#' This takes the pool already built; `determine_study_shared_residual_range` builds it from
#' the neighborhood residuals of two studies first, if that is what is at hand.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param abs_residual_pool a numeric vector. The absolute residual values of the concatenated S1,S2 residuals
#'   NaN is permitted for this value.
#' @param abs_residual_pool_perm a integer vector. The permutation vector that sorts `abs_residual_pool`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `pool_size`.
#' @param residual_range_quantile a numeric scalar. Quantile for determining the residual range
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `95.0`.
#' @return a numeric scalar. Computed residual range (R)
#' @export
determine_shared_residual_range_expert <- function(abs_residual_pool, abs_residual_pool_perm, residual_range_quantile = 95.0) {
    abs_residual_pool <- .tox_as_double_vector(abs_residual_pool, "abs_residual_pool")
    abs_residual_pool_perm <- .tox_as_integer_vector(abs_residual_pool_perm, "abs_residual_pool_perm")
    residual_range_quantile <- .tox_as_double_scalar(residual_range_quantile, "residual_range_quantile")
    if (length(abs_residual_pool_perm) != length(abs_residual_pool))
        .tox_shape_error("abs_residual_pool_perm", length(abs_residual_pool_perm), "abs_residual_pool", length(abs_residual_pool))

    .result <- .Call("determine_shared_residual_range_expert_call", abs_residual_pool, abs_residual_pool_perm, residual_range_quantile)
    .arguments <- c("abs_residual_pool", "abs_residual_pool_perm", "pool_size", "shared_residual_range", "residual_range_quantile", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shared_residual_range
}

#' Compute the shared residual range [-R, R] from a pooled set of absolute residuals
#'
#' This takes the pool already built; `determine_study_shared_residual_range` builds it from
#' the neighborhood residuals of two studies first, if that is what is at hand.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param abs_residual_pool a numeric vector. The absolute residual values of the concatenated S1,S2 residuals
#'   NaN is permitted for this value.
#' @param residual_range_quantile a numeric scalar. Quantile for determining the residual range
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `95.0`.
#' @return a numeric scalar. Computed residual range (R)
#' @export
determine_shared_residual_range <- function(abs_residual_pool, residual_range_quantile = 95.0) {
    abs_residual_pool <- .tox_as_double_vector(abs_residual_pool, "abs_residual_pool")
    residual_range_quantile <- .tox_as_double_scalar(residual_range_quantile, "residual_range_quantile")
    .result <- .Call("determine_shared_residual_range_call", abs_residual_pool, residual_range_quantile)
    .arguments <- c("abs_residual_pool", "pool_size", "shared_residual_range", "residual_range_quantile", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shared_residual_range
}

#' Compute the shared residual range [-R, R] from the neighborhood residuals of two studies
#'
#' Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
#' as `determine_shared_residual_range` does.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param residual_range_quantile a numeric scalar. Quantile for determining the residual range
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `95.0`.
#' @return a numeric scalar. Computed residual range (R)
#' @export
determine_study_shared_residual_range_expert <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile = 95.0) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    residual_range_quantile <- .tox_as_double_scalar(residual_range_quantile, "residual_range_quantile")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("determine_study_shared_residual_range_expert_call", neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "tmp_abs_residual_pool", "tmp_abs_residual_pool_perm", "shared_residual_range", "residual_range_quantile", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shared_residual_range
}

#' Compute the shared residual range [-R, R] from the neighborhood residuals of two studies
#'
#' Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
#' as `determine_shared_residual_range` does.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param residual_range_quantile a numeric scalar. Quantile for determining the residual range
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `95.0`.
#' @return a numeric scalar. Computed residual range (R)
#' @export
determine_study_shared_residual_range <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile = 95.0) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    residual_range_quantile <- .tox_as_double_scalar(residual_range_quantile, "residual_range_quantile")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("determine_study_shared_residual_range_call", neighborhood_residuals_S1, neighborhood_residuals_S2, residual_range_quantile)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "shared_residual_range", "residual_range_quantile", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shared_residual_range
}

#' Summarize the neighborhood residuals in absolute histogram counts and probability mass functions
#'
#' The probability mass function `pmf(residual, bin)` is actually a matrix.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param neighborhood_residuals a numeric array of rank 3. Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param shared_residual_range a numeric scalar. Computed residual range (R)
#'   The minimum valid value is `0.0`.
#' @param n_bins a integer scalar. Number of equally sized histogram bins in range [-R,R]
#' @param neighbor_mask a logical matrix. Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
#' @return a named list with elements:
#'   \item{counts}{a integer matrix. Absolute counts of a residual per bin}
#'   \item{pmf}{a numeric matrix. `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`}
#'   \item{included_n_reps}{a integer vector. Stores the count of non-NaN replicates (included ones)}
#' @export
build_residual_histograms <- function(neighborhood_residuals, shared_residual_range, n_bins, neighbor_mask = NULL) {
    neighborhood_residuals <- .tox_as_double_array(neighborhood_residuals, "neighborhood_residuals", 3L)
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    if (!is.null(neighbor_mask))
        neighbor_mask <- .tox_as_logical(neighbor_mask, "neighbor_mask")
    .result <- .Call("build_residual_histograms_call", neighborhood_residuals, shared_residual_range, n_bins, neighbor_mask)
    .arguments <- c("neighborhood_residuals", "n_reps", "n_neighbors", "n_points", "shared_residual_range", "n_bins", "counts", "pmf", "included_n_reps", "neighbor_mask", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        counts = .result$counts,
        pmf = .result$pmf,
        included_n_reps = .result$included_n_reps
    )
}

#' Compute the Jensen-Shannon divergence per reference point from two histograms
#'
#' Takes the probabilities `pmf` produced by `build_residual_histograms`.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param pmf_S1 a numeric matrix. Computed normalized histogram counts for study 1
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @param pmf_S2 a numeric matrix. Computed normalized histogram counts for study 2
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a numeric vector. Jensen-Shannon divergence per reference point
#' @export
compute_divergence_per_reference_point <- function(pmf_S1, pmf_S2) {
    pmf_S1 <- .tox_as_double_matrix(pmf_S1, "pmf_S1")
    pmf_S2 <- .tox_as_double_matrix(pmf_S2, "pmf_S2")
    if (dim(pmf_S2)[1] != dim(pmf_S1)[1])
        .tox_shape_error("pmf_S2", dim(pmf_S2)[1], "pmf_S1", dim(pmf_S1)[1])
    if (dim(pmf_S2)[2] != dim(pmf_S1)[2])
        .tox_shape_error("pmf_S2", dim(pmf_S2)[2], "pmf_S1", dim(pmf_S1)[2])

    .result <- .Call("compute_divergence_per_reference_point_call", pmf_S1, pmf_S2)
    .arguments <- c("pmf_S1", "pmf_S2", "n_points", "n_bins", "js_divergences", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$js_divergences
}

#' Compute the global weighted Jensen-Shannon divergence from the per-neighbor divergences
#'
#' Takes the divergences produced by `compute_divergence_per_reference_point`.
#'
#' Generated from the Fortran module \code{tox_data_integration_jsd}.
#'
#' @param js_divergences a numeric vector. Jensen-Shannon divergence per reference point, computed for studies S1 and S2
#'   The minimum valid value is `0.0`.
#' @param included_n_reps_S1 a integer vector. Count of non-NaN residuals (included ones) in study 1
#'   The minimum valid value is `0`.
#' @param included_n_reps_S2 a integer vector. Count of non-NaN residuals (included ones) in study 2
#'   The minimum valid value is `0`.
#' @return a named list with elements:
#'   \item{global_js_divergence}{a numeric scalar. Weighted global Jensen-Shannon divergence}
#'   \item{weights}{a numeric vector. Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`}
#' @export
compute_weighted_global_divergence <- function(js_divergences, included_n_reps_S1, included_n_reps_S2) {
    js_divergences <- .tox_as_double_vector(js_divergences, "js_divergences")
    included_n_reps_S1 <- .tox_as_integer_vector(included_n_reps_S1, "included_n_reps_S1")
    included_n_reps_S2 <- .tox_as_integer_vector(included_n_reps_S2, "included_n_reps_S2")
    if (length(included_n_reps_S1) != length(js_divergences))
        .tox_shape_error("included_n_reps_S1", length(included_n_reps_S1), "js_divergences", length(js_divergences))
    if (length(included_n_reps_S2) != length(js_divergences))
        .tox_shape_error("included_n_reps_S2", length(included_n_reps_S2), "js_divergences", length(js_divergences))

    .result <- .Call("compute_weighted_global_divergence_call", js_divergences, included_n_reps_S1, included_n_reps_S2)
    .arguments <- c("js_divergences", "n_points", "included_n_reps_S1", "included_n_reps_S2", "global_js_divergence", "weights", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        global_js_divergence = .result$global_js_divergence,
        weights = .result$weights
    )
}
