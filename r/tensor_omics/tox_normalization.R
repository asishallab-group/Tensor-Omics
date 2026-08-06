# Generated. Do not edit.

#' Normalizes an input vector to unit length in-place
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_unit_length}, whose argument names
#' are the ones an error message reports.
#'
#' @param vector a numeric vector. Vector that will be normalized to unit length
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @return a numeric vector. Vector that will be normalized to unit length
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
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
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param reps_per_tissue a integer vector. Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
#'   e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
#' @param span a numeric scalar. LOESS span parameter.
#'   The default value is `0.7`.
#' @param degree a integer scalar. LOESS degree parameter.
#'   The default value is `2`.
#' @param use_quantile a logical scalar. Use quantile normalization.
#'   The default value is `FALSE`.
#' @return a numeric matrix. Log-transformed grouped `expr`
#' @export
normalization_pipeline_expert <- function(expr, reps_per_tissue, span = 0.7, degree = 2L, use_quantile = FALSE) {
    expr <- .tox_as_double_matrix(expr, "expr")
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    use_quantile <- .tox_as_logical(use_quantile, "use_quantile")
    .tox_loess_required_workspace_result <- tox_loess_required_workspace(n_dim = 1L, max_neighborhood_size = dim(expr)[2], save_factorization = FALSE)
    int_workspace_size <- .tox_loess_required_workspace_result$int_workspace_size
    real_workspace_size <- .tox_loess_required_workspace_result$real_workspace_size

    .result <- .Call("normalization_pipeline_expert_call", expr, reps_per_tissue, int_workspace_size, real_workspace_size, span, degree, use_quantile)
    .arguments <- c("n_genes", "n_replicates", "expr", "log_transformed_expr", "reps_per_tissue", "n_tissues", "tmp_expr_copy", "tmp_loess_y", "tmp_indices_used", "tmp_yhat_global", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "tmp_loess_weights", "tmp_eval_points", "tmp_robust_weights", "tmp_combined_weights", "tmp_residuals", "tmp_permutation_indices", "span", "degree", "use_quantile", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, "log_transformed_expr", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "tmp_int_workspace", NA_character_, "tmp_real_workspace", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$log_transformed_expr
}

#' Complete normalization pipeline for gene expression data.
#'
#' Final result is in log_transformed_expr. If fold change is needed, call calc_fchange separately.
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalization_pipeline_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param reps_per_tissue a integer vector. Number of replicates per tissue in `expr`. It describes, which slices in `expr` relate to which tissue,
#'   e.g. `[2,3]` means `5` total replicates per gene, the first two of which belong to the first tissue and the remaining three to the second.
#' @param span a numeric scalar. LOESS span parameter.
#'   The default value is `0.7`.
#' @param degree a integer scalar. LOESS degree parameter.
#'   The default value is `2`.
#' @param use_quantile a logical scalar. Use quantile normalization.
#'   The default value is `FALSE`.
#' @return a numeric matrix. Log-transformed grouped `expr`
#' @export
normalization_pipeline <- function(expr, reps_per_tissue, span = 0.7, degree = 2L, use_quantile = FALSE) {
    expr <- .tox_as_double_matrix(expr, "expr")
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    use_quantile <- .tox_as_logical(use_quantile, "use_quantile")
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
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_by_std_dev}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param span a numeric scalar. LOESS span parameter.
#'   The default value is `0.7`.
#' @param degree a integer scalar. LOESS degree parameter.
#'   The default value is `2`.
#' @return a numeric matrix. Normalized `expr`
#' @export
normalize_by_std_dev_expert <- function(expr, span = 0.7, degree = 2L) {
    expr <- .tox_as_double_matrix(expr, "expr")
    span <- .tox_as_double_scalar(span, "span")
    degree <- .tox_as_integer_scalar(degree, "degree")
    .tox_loess_required_workspace_result <- tox_loess_required_workspace(n_dim = 1L, max_neighborhood_size = dim(expr)[2], save_factorization = FALSE)
    int_workspace_size <- .tox_loess_required_workspace_result$int_workspace_size
    real_workspace_size <- .tox_loess_required_workspace_result$real_workspace_size

    .result <- .Call("normalize_by_std_dev_expert_call", expr, int_workspace_size, real_workspace_size, span, degree)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "tmp_loess_x", "tmp_loess_y", "tmp_indices_used", "tmp_yhat_global", "tmp_int_workspace", "int_workspace_size", "tmp_real_workspace", "real_workspace_size", "tmp_hat_diag", "tmp_loess_weights", "tmp_eval_points", "tmp_robust_weights", "tmp_combined_weights", "tmp_residuals", "tmp_permutation_indices", "span", "degree", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, "tmp_int_workspace", NA_character_, "tmp_real_workspace", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$normalized_expr
}

#' Normalizes each gene's expression vector using LOESS-stabilized standard deviation.
#'
#' This procedure applies a global stabilization based on the relationship between
#' gene-wise mean expression and empirical standard deviation.
#'
#' Generated from the Fortran procedure \code{tox_normalization::normalize_by_std_dev_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param span a numeric scalar. LOESS span parameter.
#'   The default value is `0.7`.
#' @param degree a integer scalar. LOESS degree parameter.
#'   The default value is `2`.
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
#'
#' Generated from the Fortran procedure \code{tox_normalization::root_mean_sq_normalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
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
#'
#' Generated from the Fortran procedure \code{tox_normalization::quantile_normalization}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @return a named list with elements:
#'   \item{normalized_expr}{a numeric matrix. Normalized `expr`}
#'   \item{rank_means}{a numeric vector. The mean of each rank across tissues, one per gene}
#' @export
quantile_normalization_expert <- function(expr) {
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("quantile_normalization_expert_call", expr)
    .arguments <- c("n_genes", "n_replicates", "expr", "normalized_expr", "rank_means", "tmp_genes_row", "tmp_perm", "ierr")
    .sources <- c("expr", "expr", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        normalized_expr = .result$normalized_expr,
        rank_means = .result$rank_means
    )
}

#' Quantile normalization of a gene expression matrix (F42-compliant).
#'
#' Computes average expression per rank across tissues.
#'
#' Generated from the Fortran procedure \code{tox_normalization::quantile_normalization_alloc}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
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
#' matrix flattened in column-major order. The `log2` is computed via:
#' `log(x + 1) / log(2)`, which is numerically equivalent and avoids the
#' non-portable `log2` intrinsic for compatibility with WebAssembly (WASM).
#'
#' Generated from the Fortran procedure \code{tox_normalization::log2_transformation}, whose argument names
#' are the ones an error message reports.
#'
#' @param expr a numeric matrix. Gene Expression matrix, from \code{\link{calc_tiss_avg}}
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
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
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @return a numeric matrix. Tissue averages per gene
#' @export
calc_tiss_avg <- function(reps_per_tissue, expr) {
    reps_per_tissue <- .tox_as_integer_vector(reps_per_tissue, "reps_per_tissue")
    expr <- .tox_as_double_matrix(expr, "expr")
    .result <- .Call("calc_tiss_avg_call", reps_per_tissue, expr)
    .arguments <- c("n_genes", "n_tissues", "reps_per_tissue", "expr", "tissue_averages", "ierr")
    .sources <- c("expr", "reps_per_tissue", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$tissue_averages
}

#' Calculate `log2 fold changes` between condition and control groups.
#'
#' For each control-condition pair, this subroutine computes the `log2 fold change`
#' by subtracting the expression value in the control group from the corresponding
#' value in the condition group, for all genes.
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
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
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
