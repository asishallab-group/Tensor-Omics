# Generated. Do not edit.

#' Normalizes an input vector to unit length in-place
#'
#' @param vector a numeric vector. Vector that will be normalized to unit length
#' @return Vector that will be normalized to unit length
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_unit_length}.
#' @export
normalize_unit_length <- function(vector) {
    vector <- .tox_as_double_vector(vector, "vector")
    .result <- .normalize_unit_length_rcpp(vector)
    .arguments <- c("vector", "n_dims", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$vector
}

#' Complete normalization pipeline for gene expression data.
#'
#' Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @param reps_per_tissue a integer vector. Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
#' @param span a numeric scalar. LOESS span parameter.
#' @param degree a integer scalar. LOESS degree parameter.
#' @param use_quantile a logical scalar. Use quantile normalization.
#' @return Log-transformed grouped `expr`
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalization_pipeline_alloc}.
#' @export
normalization_pipeline <- function(expr, reps_per_tissue, span = 0.7, degree = 2L, use_quantile = FALSE) {
    expr <- .tox_as_double_matrix(expr, "expr")
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    use_quantile <- .tox_as_logical(use_quantile, "use_quantile")
    .result <- .normalization_pipeline_rcpp(expr, reps_per_tissue, span, degree, use_quantile)
    .arguments <- c("n_genes", "n_replicates", "expr", "log_transformed_expr", "reps_per_tissue", "n_tissues", "span", "degree", "use_quantile", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$log_transformed_expr
}

#' Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
#'
#' This procedure applies a global stabilization based on the relationship between
#' gene-wise mean expression and empirical standard deviation.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @param span a numeric scalar. LOESS span parameter.
#' @param degree a integer scalar. LOESS degree parameter.
#' @return Normalized `expr`
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_by_std_dev_alloc}.
#' @export
normalize_by_std_dev <- function(expr, span = 0.7, degree = 2L) {
    expr <- .tox_as_double_matrix(expr, "expr")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    .result <- .normalize_by_std_dev_rcpp(expr, span, degree)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "span", "degree", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$normalized_expr
}

#' Normalizes each gene's expression vector using `sqrt(mean(x^2))`
#'
#' across tissues (not classical standard deviation).
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @return Normalized `expr`
#'
#' Generated from the Fortran procedure \code{tox_normalization::root_mean_sq_normalization}.
#' @export
root_mean_sq_normalization <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .root_mean_sq_normalization_rcpp(expr)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$normalized_expr
}

#' Quantile normalization of a gene expression matrix (F42-compliant).
#'
#' Computes average expression per rank across tissues.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#' @return a named list with elements `normalized_expr`, `rank_means`.
#'
#' Generated from the Fortran procedure \code{tox_normalization::quantile_normalization}.
#' @export
quantile_normalization <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .quantile_normalization_rcpp(expr)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "rank_means", "tmp_genes_row", "tmp_perm", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        normalized_expr = .result$normalized_expr,
        rank_means = .result$rank_means
    )
}

#' Apply `log2(x + 1)` transformation to each element of the input matrix.
#'
#' This subroutine performs element-wise `log2(x + 1)` transformation on a
#' matrix flattened in column-major order. The `log2` is computed via:
#' `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
#' non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
#'
#' @param expr a numeric matrix. Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
#' @return Log-transformed `expr`
#'
#' Generated from the Fortran procedure \code{tox_normalization::log2_transformation}.
#' @export
log2_transformation <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .log2_transformation_rcpp(expr)
    .arguments <- c("n_genes", "n_tissues", "expr", "transformed_expr", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$transformed_expr
}

#' Calculate tissue averages by averaging replicates within each tissue.
#'
#' For each tissue of tissue replicates, this subroutine computes the average
#' expression per gene.
#'
#' @param reps_per_tissue a integer vector. Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
#' @param expr a numeric matrix. Gene Expression matrix
#' @return Tissue averages per gene
#'
#' Generated from the Fortran procedure \code{tox_normalization::calc_tiss_avg}.
#' @export
calc_tiss_avg <- function(reps_per_tissue, expr) {
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .calc_tiss_avg_rcpp(reps_per_tissue, expr)
    .arguments <- c("n_genes", "n_tissues", "reps_per_tissue", "expr", "tissue_averages", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$tissue_averages
}

#' Calculate `log2 fold changes` between condition and control groups.
#'
#' For each control-condition pair, this subroutine computes the `log2 fold change`
#' by subtracting the expression value in the control group from the corresponding
#' value in the condition group, for all genes.
#'
#' @param control_tissues a integer vector. Control tissue indices
#' @param condition_tissues a integer vector. Condition tissue indices
#' @param expr a numeric matrix. Gene Expression matrix, from [[tox_normalization(module):calc_tiss_avg(subroutine)]]
#' @return Output matrix for fold changes
#'
#' Generated from the Fortran procedure \code{tox_normalization::calc_fchange}.
#' @export
calc_fchange <- function(control_tissues, condition_tissues, expr) {
    control_tissues <- .tox_as_integer_vector(control_tissues, "control_tissues")
    condition_tissues <- .tox_as_integer_vector(condition_tissues, "condition_tissues")
    expr <- .tox_as_double_matrix(expr, "expr")
    if (length(condition_tissues) != length(control_tissues))
        .tox_shape_error("condition_tissues", length(condition_tissues), "control_tissues", length(control_tissues))

    .result <- .calc_fchange_rcpp(control_tissues, condition_tissues, expr)
    .arguments <- c("n_genes", "n_tissues", "n_pairs", "control_tissues", "condition_tissues", "expr", "fold_changes", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$fold_changes
}
