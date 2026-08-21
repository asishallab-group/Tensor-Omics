# Generated. Do not edit.

#' Compute family scaling factors (dscale) to normalize distances
#'
#' Uses LOESS on the median/stddev of intra-family distances for scaling, regardless of orthologs.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers::compute_family_scaling}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_families a integer scalar. Total number of gene families
#' @param distances a numeric vector. Array of Euclidean distances for each gene
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param gene_to_fam a integer vector. Mapping of each gene to its family (1-based)
#' @param span a numeric scalar. Span parameter for LOESS smoothing, passed straight to
#'   \code{\link{loess_fit_plain}}, so it is held to that
#'   procedure's own range rather than to the NaN tolerance the distance data carries.
#'   The default value is `0.7`.
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. Degree of the LOESS polynomial
#'   The default value is `2`.
#' @param mode a string, one of "plain", "robust". Mode for LOESS fitting
#'   The default value is `"robust"`.
#' @param n_iters a integer scalar. Number of iterations for robust LOESS fitting
#'   The default value is `3`.
#' @return a named list with elements:
#'   \item{dscale}{a numeric vector. Array of scaling factors per family (output)}
#'   \item{loess_x}{a numeric vector. Reference x-coordinates for LOESS smoothing}
#'   \item{loess_y}{a numeric vector. Reference y-coordinates for LOESS smoothing}
#'   \item{indices_used}{a integer vector. Indices of reference points used for smoothing}
#'   \item{low_sd_cutoff}{a numeric scalar. cutoff used to filter families with low std}
#'   \item{excluded_low_sd}{a integer vector. Mask to save those families that have low sd}
#' @export
compute_family_scaling <- function(n_families, distances, gene_to_fam, span = 0.7, degree = 2L, mode = "robust", n_iters = 3L) {
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    distances <- .tox_as_double_vector(distances, "distances")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    mode <- .tox_as_mode(mode, "mode", c("plain", "robust"))
    n_iters <- .tox_as_integer_scalar(n_iters, "n_iters")
    if (length(gene_to_fam) != length(distances))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "distances", length(distances))

    .result <- .Call("compute_family_scaling_call", n_families, distances, gene_to_fam, span, degree, mode, n_iters)
    .arguments <- c("n_genes", "n_families", "distances", "gene_to_fam", "dscale", "loess_x", "loess_y", "indices_used", "span", "degree", "mode", "n_iters", "low_sd_cutoff", "excluded_low_sd", "ierr")
    .sources <- c("distances", "dscale", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        dscale = .result$dscale,
        loess_x = .result$loess_x,
        loess_y = .result$loess_y,
        indices_used = .result$indices_used,
        low_sd_cutoff = .result$low_sd_cutoff,
        excluded_low_sd = .result$excluded_low_sd
    )
}

#' Compute the hybrid RDI (Relative Distance Index) for each gene
#'
#' RDI = Euclidean distance / family scaling factor
#'
#' Generated from the Fortran procedure \code{tox_get_outliers::compute_rdi}, whose argument names
#' are the ones an error message reports.
#'
#' @param distances a numeric vector. Array of Euclidean distances for each gene to its centroid
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param gene_to_fam a integer vector. Gene-to-family mapping (1-based indexing)
#' @param dscale a numeric vector. Array of scaling factors for each family
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @return a named list with elements:
#'   \item{rdi}{a numeric vector. Output array of RDI values for each gene}
#'   \item{sorted_rdi}{a numeric vector. Work array for sorting (dimension n_genes)}
#'   \item{perm}{a integer vector. Permutation array for sorting (dimension n_genes, should be pre-initialized with 1:n_genes)}
#' @export
compute_rdi <- function(distances, gene_to_fam, dscale) {
    distances <- .tox_as_double_vector(distances, "distances")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    dscale <- .tox_as_double_vector(dscale, "dscale")
    if (length(gene_to_fam) != length(distances))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "distances", length(distances))

    .result <- .Call("compute_rdi_call", distances, gene_to_fam, dscale)
    .arguments <- c("n_genes", "distances", "gene_to_fam", "dscale", "rdi", "sorted_rdi", "perm", "ierr")
    .sources <- c("distances", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        rdi = .result$rdi,
        sorted_rdi = .result$sorted_rdi,
        perm = .result$perm
    )
}

#' Identify gene outliers based on the top percentile of RDI values
#'
#' Expects sorted_rdi to be filtered (no negative values) and perm should be sorted in ascending order before calling.
#' If sorted_rdi contains negatives or perm is not sorted, tmp_results may be invalid.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers::identify_outliers}, whose argument names
#' are the ones an error message reports.
#'
#' @param rdi a numeric vector. Array of RDI values for each gene
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param sorted_rdi a numeric vector. Sorted RDI array (must be filtered to remove negatives and sorted in ascending order before calling)
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param perm a integer vector. Permutation array with sorted indices
#' @param percentile a numeric scalar. Percentile threshold as a fraction in [0,1] (top 5% for the default).
#'   The default value is `0.95`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a named list with elements:
#'   \item{is_outlier}{a logical vector. Output boolean array indicating outliers}
#'   \item{threshold}{a numeric scalar. Output threshold value used for detection}
#'   \item{quantile}{a numeric vector. Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
#'     observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
#'     Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
#'     upper-tail quantile is used.}
#' @export
identify_outliers <- function(rdi, sorted_rdi, perm, percentile = 0.95) {
    rdi <- .tox_as_double_vector(rdi, "rdi")
    sorted_rdi <- .tox_as_double_vector(sorted_rdi, "sorted_rdi")
    perm <- .tox_as_integer_vector(perm, "perm")
    percentile <- .tox_as_double_scalar(percentile, "percentile")
    if (length(sorted_rdi) != length(rdi))
        .tox_shape_error("sorted_rdi", length(sorted_rdi), "rdi", length(rdi))
    if (length(perm) != length(rdi))
        .tox_shape_error("perm", length(perm), "rdi", length(rdi))

    .result <- .Call("identify_outliers_call", rdi, sorted_rdi, perm, percentile)
    .arguments <- c("n_genes", "rdi", "sorted_rdi", "perm", "is_outlier", "threshold", "quantile", "percentile", "ierr")
    .sources <- c("rdi", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        is_outlier = .result$is_outlier,
        threshold = .result$threshold,
        quantile = .result$quantile
    )
}

#' Main routine to detect outliers using RDI and LOESS-based scaling
#'
#' Orchestrates the full pipeline: per-family scaling via
#' \code{\link{compute_family_scaling}}, the RDI per gene via
#' \code{\link{compute_rdi}}, then flags outliers via
#' \code{\link{identify_outliers}}.
#'
#' Generated from the Fortran procedure \code{tox_get_outliers::detect_outliers}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_families a integer scalar. Total number of gene families
#' @param distances a numeric vector. Array of Euclidean distances for each gene to its centroid
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param gene_to_fam a integer vector. Gene-to-family mapping (1-based indexing)
#' @param percentile a numeric scalar. Percentile threshold as a fraction in [0,1] for outlier detection.
#'   The default value is `0.95`.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#' @return a named list with elements:
#'   \item{is_outlier}{a logical vector. Output boolean array indicating outliers}
#'   \item{loess_x}{a numeric vector. Reference x-coordinates.}
#'   \item{loess_y}{a numeric vector. Reference y-coordinates (length n_total).}
#'   \item{loess_n}{a integer vector. Indices of reference points used for smoothing.}
#'   \item{quantile}{a numeric vector. Empirical one-sided upper-tail quantile (effect-size measure) for each gene, i.e. how extreme an
#'     observed distance is relative to all observed distances -- NOT a null-hypothesis-testing p-value.
#'     Returned in the same order as the input RDI array. Because distances are non-negative, a one-sided
#'     upper-tail quantile is used.}
#' @export
detect_outliers <- function(n_families, distances, gene_to_fam, percentile = 0.95) {
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    distances <- .tox_as_double_vector(distances, "distances")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    percentile <- .tox_as_double_scalar(percentile, "percentile")
    if (length(gene_to_fam) != length(distances))
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "distances", length(distances))

    .result <- .Call("detect_outliers_call", n_families, distances, gene_to_fam, percentile)
    .arguments <- c("n_genes", "n_families", "distances", "gene_to_fam", "is_outlier", "loess_x", "loess_y", "loess_n", "quantile", "ierr", "percentile")
    .sources <- c("distances", "loess_x", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        is_outlier = .result$is_outlier,
        loess_x = .result$loess_x,
        loess_y = .result$loess_y,
        loess_n = .result$loess_n,
        quantile = .result$quantile
    )
}
