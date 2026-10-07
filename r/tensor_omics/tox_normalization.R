# Generated. Do not edit.

#' Normalizes an input vector to unit length in-place
#'
#' Only the zero vector has no direction, and is rejected with ERR_DIVISION_BY_ZERO. Any other
#' vector becomes its unit vector, however short or long: its own length is never formed, so it
#' may underflow or overflow, as it does for entries near the smallest normal number (about
#' 2.2e-308) or near the largest number (about 1.8e308). Where subnormal numbers are flushed to
#' zero, by a fast floating-point model or by the host program, a vector made only of
#' subnormals reads as the zero vector.
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_unit_length}, whose argument names
#' are the ones an error message reports.
#'
#' @param vector a numeric vector. Vector that will be normalized to unit length
#' @return a numeric vector. Vector that will be normalized to unit length
#' @export
normalize_unit_length <- function(vector) {
    vector <- .tox_as_double_vector(vector, "vector")
    .result <- .Call("normalize_unit_length_call", vector)
    .arguments <- c("vector", "n_dims", "ierr")
    .sources <- c(NA_character_, "vector", NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$vector
}

#' Complete normalization pipeline for gene expression data.
#'
#' Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalization_pipeline}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @param reps_per_tissue a integer vector. Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
#'   e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
#' @param span a numeric scalar. LOESS span parameter.
#'   The default value is `0.7`.
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. LOESS degree parameter.
#'   The default value is `2`.
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @param use_quantile a logical scalar. Use quantile normalization.
#'   The default value is `FALSE`.
#' @return a numeric matrix. Log-transformed grouped `expr`
#' @export
normalization_pipeline <- function(expr, reps_per_tissue, span = 0.7, degree = 2L, use_quantile = FALSE) {
    expr <- .tox_as_double_matrix(expr, "expr")
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    use_quantile <- .tox_as_logical_scalar(use_quantile, "use_quantile")
    .result <- .Call("normalization_pipeline_call", expr, reps_per_tissue, span, degree, use_quantile)
    .arguments <- c("n_genes", "n_replicates", "expr", "log_transformed_expr", "reps_per_tissue", "n_tissues", "span", "degree", "use_quantile", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, "log_transformed_expr", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$log_transformed_expr
}

#' Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
#'
#' This procedure applies a global stabilization based on the relationship between
#' gene-wise mean expression and empirical standard deviation.
#' Genes whose standard deviation is exactly zero carry no information about the trend: they
#' are left out of the fit and returned unchanged. Where the fitted trend is at or below zero
#' -- a LOESS fit can dip below zero even on non-negative data -- a gene is divided by its own
#' standard deviation instead, so no gene changes sign. A spread or a fit that is merely small
#' is used as it is.
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_by_std_dev}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @param span a numeric scalar. LOESS span parameter.
#'   The default value is `0.7`.
#'   The minimum valid value is `EPS_LOESS`.
#'   The maximum valid value is `1.0`.
#' @param degree a integer scalar. LOESS degree parameter.
#'   The default value is `2`.
#'   The minimum valid value is `0`.
#'   The maximum valid value is `2`.
#' @return a numeric matrix. Normalized `expr`
#' @export
normalize_by_std_dev <- function(expr, span = 0.7, degree = 2L) {
    expr <- .tox_as_double_matrix(expr, "expr")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    .result <- .Call("normalize_by_std_dev_call", expr, span, degree)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "span", "degree", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$normalized_expr
}

#' Normalizes each gene's expression vector using `sqrt(mean(x^2))`
#'
#' across tissues (not classical standard deviation).
#' Only a gene whose values are all exactly zero has no root mean square, and is left as it is;
#' every other gene is divided by it, however small or large its values, up to the largest
#' real64 (about 1.8e308).
#'
#' Generated from the Fortran procedure \code{tox_normalization::root_mean_sq_normalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @return a numeric matrix. Normalized `expr`
#' @export
root_mean_sq_normalization <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("root_mean_sq_normalization_call", expr)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$normalized_expr
}

#' Quantile normalization of a gene expression matrix (F42-compliant).
#'
#' Computes average expression per rank across tissues.
#' Tied values within a replicate share the mean of the rank means their ranks span, so values
#' that are equal before normalization stay equal after it, as in `preprocessCore` and limma's
#' `normalizeQuantiles`. Where those rank means are all equal, the tie gets exactly that value.
#' The rank means themselves do not depend on ties.
#'
#' No mean overflows for finite input, however close to the largest real64 (about 1.8e308) the
#' values are, and a single replicate comes back exactly as it was, ties included.
#'
#' Generated from the Fortran procedure \code{tox_normalization::quantile_normalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @return a named list with elements:
#'   \item{normalized_expr}{a numeric matrix. Normalized `expr`}
#'   \item{rank_means}{a numeric vector. The mean of each rank across tissues, one per gene}
#' @export
quantile_normalization <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("quantile_normalization_call", expr)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "rank_means", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        normalized_expr = .result$normalized_expr,
        rank_means = .result$rank_means
    )
}

#' Apply `log2(x + 1)` transformation to each element of the input matrix.
#'
#' This subroutine performs element-wise `log2(x + 1)` transformation on a
#' matrix flattened in column-major order. The `log2` is computed as `log1p(x)/log(2)`:
#' `log1p` keeps the digits of a tiny `x` that forming `x + 1` would round away, and dividing
#' by `log(2)` avoids the non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
#'
#' Generated from the Fortran procedure \code{tox_normalization::log2_transformation}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix, from \code{\link{calc_tiss_avg}}
#' @return a numeric matrix. Log-transformed `expr`
#' @export
log2_transformation <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("log2_transformation_call", expr)
    .arguments <- c("n_genes", "n_tissues", "expr", "transformed_expr", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$transformed_expr
}

#' Calculate tissue averages by averaging replicates within each tissue.
#'
#' For each tissue of tissue replicates, this subroutine computes the average
#' expression per gene.
#'
#' Generated from the Fortran procedure \code{tox_normalization::calc_tiss_avg}, whose argument names
#' are the ones an error message reports.
#'
#' @param reps_per_tissue a integer vector. Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
#'   e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
#'   The minimum valid value is `1`.
#' @param expr a numeric matrix. Gene Expression matrix
#' @return a numeric matrix. Tissue averages per gene
#' @export
calc_tiss_avg <- function(reps_per_tissue, expr) {
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("calc_tiss_avg_call", reps_per_tissue, expr)
    .arguments <- c("n_genes", "n_replicates", "n_tissues", "reps_per_tissue", "expr", "tissue_averages", "ierr")
    .sources <- c("expr", "expr", "reps_per_tissue", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$tissue_averages
}

#' Calculate `log2 fold changes` between condition and control groups.
#'
#' For each control-condition pair, this subroutine computes the `log2 fold change`
#' by subtracting the expression value in the control group from the corresponding
#' value in the condition group, for all genes.
#' A difference too large for real64 -- possible only near `huge`, as in `huge - (-huge)` -- is
#' reported as ERR_NAN_INF instead of being written into the result as Inf.
#'
#' Generated from the Fortran procedure \code{tox_normalization::calc_fchange}, whose argument names
#' are the ones an error message reports.
#'
#' @param control_tissues a integer vector. Control tissue indices
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_tissues`.
#' @param condition_tissues a integer vector. Condition tissue indices
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_tissues`.
#' @param expr a numeric matrix. Gene Expression matrix, from \code{\link{calc_tiss_avg}}
#' @return a numeric matrix. Output matrix for fold changes
#' @export
calc_fchange <- function(control_tissues, condition_tissues, expr) {
    control_tissues <- .tox_as_integer_vector(control_tissues, "control_tissues")
    condition_tissues <- .tox_as_integer_vector(condition_tissues, "condition_tissues")
    expr <- .tox_as_double_matrix(expr, "expr")
    if (length(condition_tissues) != length(control_tissues))
        .tox_shape_error("condition_tissues", length(condition_tissues), "control_tissues", length(control_tissues))

    .result <- .Call("calc_fchange_call", control_tissues, condition_tissues, expr)
    .arguments <- c("n_genes", "n_tissues", "n_pairs", "control_tissues", "condition_tissues", "expr", "fold_changes", "ierr")
    .sources <- c("expr", "expr", "control_tissues", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$fold_changes
}
