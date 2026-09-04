# Generated. Do not edit.

#' Compute the family-level compatibility score between two studies for a single gene family
#'
#' Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
#' the residual samples to the genes belonging to the family `family_idx`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_per_family::fjct_compute_jsd}, whose argument names
#' are the ones an error message reports.
#'
#' @param family_idx a integer scalar. Index of the family that should be analyzed
#'   The minimum valid value is `1`.
#' @param gene_to_family_S1 a integer vector. Mapping for study 1: Each index (gene) holds the index of its family
#'   The minimum valid value is `0`.
#' @param gene_to_family_S2 a integer vector. Mapping for study 2: Each index (gene) holds the index of its family
#'   The minimum valid value is `0`.
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighborhood_genes_S1 a integer matrix. Indices of the selected neighborhood genes of study 1
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes_S1`.
#' @param neighborhood_genes_S2 a integer matrix. Indices of the selected neighborhood genes of study 2
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes_S2`.
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies
#'   The minimum valid value is `0.0`.
#' @return a named list with elements:
#'   \item{js_divergences}{a numeric vector. Jensen-Shannon divergence per reference point, computed for studies S1 and S2}
#'   \item{included_n_reps_S1}{a integer vector. Count of non-NaN residuals (included ones) in study 1}
#'   \item{included_n_reps_S2}{a integer vector. Count of non-NaN residuals (included ones) in study 2}
#'   \item{total_included_n_reps}{a integer scalar. Total number of included replicates from both studies (\eqn{\text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2)})}
#'   \item{global_js_divergence}{a numeric scalar. Weighted global Jensen-Shannon divergence}
#'   \item{weights}{a numeric vector. Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`}
#' @export
fjct_compute_jsd <- function(family_idx, gene_to_family_S1, gene_to_family_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, neighborhood_genes_S1, neighborhood_genes_S2, n_bins, shared_residual_range) {
    family_idx <- .tox_as_integer_scalar(family_idx, "family_idx")
    gene_to_family_S1 <- .tox_as_integer_vector(gene_to_family_S1, "gene_to_family_S1")
    gene_to_family_S2 <- .tox_as_integer_vector(gene_to_family_S2, "gene_to_family_S2")
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    neighborhood_genes_S1 <- .tox_as_integer_matrix(neighborhood_genes_S1, "neighborhood_genes_S1")
    neighborhood_genes_S2 <- .tox_as_integer_matrix(neighborhood_genes_S2, "neighborhood_genes_S2")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_genes_S1)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_genes_S1", dim(neighborhood_genes_S1)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_genes_S2)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_genes_S2", dim(neighborhood_genes_S2)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighborhood_genes_S1)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_genes_S1", dim(neighborhood_genes_S1)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighborhood_genes_S2)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_genes_S2", dim(neighborhood_genes_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("fjct_compute_jsd_call", family_idx, gene_to_family_S1, gene_to_family_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, neighborhood_genes_S1, neighborhood_genes_S2, n_bins, shared_residual_range)
    .arguments <- c("family_idx", "gene_to_family_S1", "gene_to_family_S2", "n_genes_S1", "n_genes_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_genes_S1", "neighborhood_genes_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "ierr")
    .sources <- c(NA_character_, NA_character_, NA_character_, "gene_to_family_S1", "gene_to_family_S2", NA_character_, NA_character_, NA_character_, NA_character_, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        js_divergences = .result$js_divergences,
        included_n_reps_S1 = .result$included_n_reps_S1,
        included_n_reps_S2 = .result$included_n_reps_S2,
        total_included_n_reps = .result$total_included_n_reps,
        global_js_divergence = .result$global_js_divergence,
        weights = .result$weights
    )
}

#' Compute the compatibility score between two studies for a single masked sub-neighborhood
#'
#' Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
#' residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
#' those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
#' the masks for from a family index.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_per_family::fjct_compute_masked_jsd}, whose argument names
#' are the ones an error message reports.
#'
#' @param neighborhood_residuals_S1 a numeric array of rank 3. Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighborhood_residuals_S2 a numeric array of rank 3. Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
#'   NaN is permitted for this value.
#' @param neighbor_mask_S1 a logical matrix. Mask selecting the neighbors of study 1 to include
#' @param neighbor_mask_S2 a logical matrix. Mask selecting the neighbors of study 2 to include
#' @param n_bins a integer scalar. Number of equally sized histogram bins used for the studies
#' @param shared_residual_range a numeric scalar. Computed residual range for both studies
#'   The minimum valid value is `0.0`.
#' @return a named list with elements:
#'   \item{js_divergences}{a numeric vector. Jensen-Shannon divergence per reference point, computed for studies S1 and S2}
#'   \item{included_n_reps_S1}{a integer vector. Count of non-NaN residuals (included ones) in study 1}
#'   \item{included_n_reps_S2}{a integer vector. Count of non-NaN residuals (included ones) in study 2}
#'   \item{total_included_n_reps}{a integer scalar. Total number of included replicates from both studies (\eqn{\text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2)})}
#'   \item{global_js_divergence}{a numeric scalar. Weighted global Jensen-Shannon divergence}
#'   \item{weights}{a numeric vector. Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`}
#'   \item{pmf_S1}{a numeric matrix. Normalized histogram counts for study 1}
#'   \item{pmf_S2}{a numeric matrix. Normalized histogram counts for study 2}
#' @export
fjct_compute_masked_jsd <- function(neighborhood_residuals_S1, neighborhood_residuals_S2, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range) {
    neighborhood_residuals_S1 <- .tox_as_double_array(neighborhood_residuals_S1, "neighborhood_residuals_S1", 3L)
    neighborhood_residuals_S2 <- .tox_as_double_array(neighborhood_residuals_S2, "neighborhood_residuals_S2", 3L)
    neighbor_mask_S1 <- .tox_as_logical(neighbor_mask_S1, "neighbor_mask_S1")
    neighbor_mask_S2 <- .tox_as_logical(neighbor_mask_S2, "neighbor_mask_S2")
    n_bins <- .tox_as_integer_scalar(n_bins, "n_bins")
    shared_residual_range <- .tox_as_double_scalar(shared_residual_range, "shared_residual_range")
    if (dim(neighborhood_residuals_S2)[2] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighbor_mask_S1)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighbor_mask_S1", dim(neighbor_mask_S1)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighbor_mask_S2)[1] != dim(neighborhood_residuals_S1)[2])
        .tox_shape_error("neighbor_mask_S2", dim(neighbor_mask_S2)[1], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[2])
    if (dim(neighborhood_residuals_S2)[3] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighborhood_residuals_S2", dim(neighborhood_residuals_S2)[3], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighbor_mask_S1)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighbor_mask_S1", dim(neighbor_mask_S1)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])
    if (dim(neighbor_mask_S2)[2] != dim(neighborhood_residuals_S1)[3])
        .tox_shape_error("neighbor_mask_S2", dim(neighbor_mask_S2)[2], "neighborhood_residuals_S1", dim(neighborhood_residuals_S1)[3])

    .result <- .Call("fjct_compute_masked_jsd_call", neighborhood_residuals_S1, neighborhood_residuals_S2, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range)
    .arguments <- c("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "neighbor_mask_S1", "neighbor_mask_S2", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "pmf_S1", "pmf_S2", "ierr")
    .sources <- c(NA_character_, NA_character_, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", NA_character_, NA_character_, "pmf_S1", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        js_divergences = .result$js_divergences,
        included_n_reps_S1 = .result$included_n_reps_S1,
        included_n_reps_S2 = .result$included_n_reps_S2,
        total_included_n_reps = .result$total_included_n_reps,
        global_js_divergence = .result$global_js_divergence,
        weights = .result$weights,
        pmf_S1 = .result$pmf_S1,
        pmf_S2 = .result$pmf_S2
    )
}

#' Compute the per-sub-neighborhood contribution score
#'
#' Combines
#'
#' 1. how divergent the family is between the studies (`global_js_divergences`), and
#' 2. how much residual support the family has overall (`total_included_n_reps_per_f`),
#'
#' using the outputs of `fjct_compute_jsd`, collected for the analyzed sub-neighborhoods.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_per_family::fjct_compute_contribution_scores}, whose argument names
#' are the ones an error message reports.
#'
#' @param global_js_divergences a numeric vector. Per-sub-neighborhood weighted global JSD
#'   The minimum valid value is `0.0`.
#' @param total_included_n_reps_per_f a integer vector. Per-sub-neighborhood `total_included_n_reps`
#'   The minimum valid value is `0`.
#' @return a named list with elements:
#'   \item{support_weights}{a numeric vector. Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)}
#'   \item{contribution_scores}{a numeric vector. Per-sub-neighborhood calculated contribution ( \eqn{support\_weights_i * global\_js\_divergences_i} )}
#' @export
fjct_compute_contribution_scores <- function(global_js_divergences, total_included_n_reps_per_f) {
    global_js_divergences <- .tox_as_double_vector(global_js_divergences, "global_js_divergences")
    total_included_n_reps_per_f <- .tox_as_integer_vector(total_included_n_reps_per_f, "total_included_n_reps_per_f")
    if (length(total_included_n_reps_per_f) != length(global_js_divergences))
        .tox_shape_error("total_included_n_reps_per_f", length(total_included_n_reps_per_f), "global_js_divergences", length(global_js_divergences))

    .result <- .Call("fjct_compute_contribution_scores_call", global_js_divergences, total_included_n_reps_per_f)
    .arguments <- c("global_js_divergences", "total_included_n_reps_per_f", "k_families", "support_weights", "contribution_scores", "ierr")
    .sources <- c(NA_character_, NA_character_, "global_js_divergences", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        support_weights = .result$support_weights,
        contribution_scores = .result$contribution_scores
    )
}
