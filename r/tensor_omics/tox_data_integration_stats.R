# Generated. Do not edit.

#' Estimate how likely the observed divergence is to occur by chance
#'
#' Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
#' the work copies, so the caller's own arrays are left untouched.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_stats::gjct_permutation_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param global_jsd_observed a numeric scalar. Observed global JSD value for both studies
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies
#'   The minimum valid value is `0.0`.
#' @param n_permutations a integer scalar. Number of permutations to perform
#' @param random_seed a integer scalar. Seed to use for shuffling
#' @param neighbor_mask_S1 a logical matrix. Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
#' @param neighbor_mask_S2 a logical matrix. Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
#' @return a named list with elements:
#'   \item{jsd_null}{a numeric vector. Vector of global divergence values obtained under the null hypothesis}
#'   \item{p_value}{a numeric scalar. Empirical p-value of the permutation test: \eqn{\frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1}}}
#' @export
gjct_permutation_test <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed = NULL, neighbor_mask_S1 = NULL, neighbor_mask_S2 = NULL) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    global_jsd_observed <- .tox_as_double_scalar(global_jsd_observed, "global_jsd_observed")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    if (!is.null(random_seed))
        random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    if (!is.null(neighbor_mask_S1))
        neighbor_mask_S1 <- .tox_as_logical_matrix(neighbor_mask_S1, "neighbor_mask_S1")
    if (!is.null(neighbor_mask_S2))
        neighbor_mask_S2 <- .tox_as_logical_matrix(neighbor_mask_S2, "neighbor_mask_S2")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("gjct_permutation_test_call", neighborhood_residuals_S1, neighborhood_residuals_S2, global_jsd_observed, n_bins, shared_residual_range, n_permutations, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2", "ierr")
    .sources <- c(NA_character_, NA_character_, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", NA_character_, NA_character_, NA_character_, "jsd_null", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        jsd_null = .result$jsd_null,
        p_value = .result$p_value
    )
}
