# Generated. Do not edit.

#' Computes the element-wise mean for a given set of vectors.
#'
#' Generated from the Fortran procedure \code{tox_gene_centroids::mean_vector}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. The input matrix of all gene expression vectors (n_axes x n_genes).
#' @param gene_indices a integer vector. An array containing the column indices of the selected genes in 'expression_vectors'.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes`.
#' @return a numeric vector. The output vector representing the computed centroid.
#' @export
mean_vector <- function(expression_vectors, gene_indices) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_indices <- .tox_as_integer_vector(gene_indices, "gene_indices")
    .result <- .Call("mean_vector_call", expression_vectors, gene_indices)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_indices", "n_selected_genes", "centroid", "ierr")
    .sources <- c(NA_character_, "expression_vectors", "expression_vectors", NA_character_, "gene_indices", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$centroid
}

#' Iterates over families, filters gene indices, and computes centroids.
#'
#' Generated from the Fortran procedure \code{tox_gene_centroids::group_centroid_orthologs}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. The input matrix of all gene expression vectors (n_axes x n_genes).
#' @param gene_to_family a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param n_families a integer scalar. Total number of gene families to compute centroids for.
#' @param n_omp_threads a integer scalar. The most threads the parallel loop may use (OpenMP): it caps the team size.
#'   It also sizes the per-thread work arrays, where there are any, as their last extent: the plain entry point allocates that many columns, one in a build without OpenMP, and the expert one takes arrays with that many.
#'   It must be at least 1 in every build; without OpenMP the loop runs on one thread whatever is passed.
#'   Pass 1 when you parallelize around this call yourself, or call it in a forked worker process.
#'   The minimum valid value is `1`.
#' @param ortholog_set a logical vector. A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
#' @return a numeric matrix. The output matrix (n_axes x n_families) to store the computed centroids.
#' @export
group_centroid_orthologs <- function(expression_vectors, gene_to_family, n_families, n_omp_threads, ortholog_set) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_family <- .tox_as_integer_vector(gene_to_family, "gene_to_family")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    n_omp_threads <- .tox_as_integer_scalar(n_omp_threads, "n_omp_threads")
    ortholog_set <- .tox_as_logical_vector(ortholog_set, "ortholog_set")
    if (length(gene_to_family) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_family", 1L, length(gene_to_family), "'expression_vectors' implies n_genes == ", dim(expression_vectors)[2], "n_genes")
    if (length(ortholog_set) != dim(expression_vectors)[2])
        .tox_shape_error("ortholog_set", 1L, length(ortholog_set), "'expression_vectors' implies n_genes == ", dim(expression_vectors)[2], "n_genes")

    .result <- .Call("group_centroid_orthologs_call", expression_vectors, gene_to_family, n_families, n_omp_threads, ortholog_set)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "n_omp_threads", "ortholog_set", "ierr")
    .sources <- c(NA_character_, "expression_vectors", "expression_vectors", NA_character_, "centroid_matrix", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$centroid_matrix
}

#' Iterates over families, filters gene indices, and computes centroids.
#'
#' Generated from the Fortran procedure \code{tox_gene_centroids::group_centroid_all}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. The input matrix of all gene expression vectors (n_axes x n_genes).
#' @param gene_to_family a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @param n_families a integer scalar. Total number of gene families to compute centroids for.
#' @param n_omp_threads a integer scalar. The most threads the parallel loop may use (OpenMP): it caps the team size.
#'   It also sizes the per-thread work arrays, where there are any, as their last extent: the plain entry point allocates that many columns, one in a build without OpenMP, and the expert one takes arrays with that many.
#'   It must be at least 1 in every build; without OpenMP the loop runs on one thread whatever is passed.
#'   Pass 1 when you parallelize around this call yourself, or call it in a forked worker process.
#'   The minimum valid value is `1`.
#' @return a numeric matrix. The output matrix (n_axes x n_families) to store the computed centroids.
#' @export
group_centroid_all <- function(expression_vectors, gene_to_family, n_families, n_omp_threads) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_family <- .tox_as_integer_vector(gene_to_family, "gene_to_family")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    n_omp_threads <- .tox_as_integer_scalar(n_omp_threads, "n_omp_threads")
    if (length(gene_to_family) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_family", 1L, length(gene_to_family), "'expression_vectors' implies n_genes == ", dim(expression_vectors)[2], "n_genes")

    .result <- .Call("group_centroid_all_call", expression_vectors, gene_to_family, n_families, n_omp_threads)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "n_omp_threads", "ierr")
    .sources <- c(NA_character_, "expression_vectors", "expression_vectors", NA_character_, "centroid_matrix", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$centroid_matrix
}
