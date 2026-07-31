# Generated. Do not edit.

#' Computes the element-wise mean for a given set of vectors.
#'
#' @param expression_vectors a numeric matrix. The input matrix of all gene expression vectors (n_axes x n_genes).
#' @param gene_indices a integer vector. An array containing the column indices of the selected genes in 'expression_vectors'.
#' @return The output vector representing the computed centroid.
#'
#' Generated from the Fortran procedure \code{tox_gene_centroids::mean_vector}.
#' @export
mean_vector <- function(expression_vectors, gene_indices) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_indices <- .tox_as_integer_vector(gene_indices, "gene_indices")
    .result <- .Call("mean_vector_call", expression_vectors, gene_indices)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_indices", "n_selected_genes", "centroid", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$centroid
}

#' Iterates over families, filters gene indices, and computes centroids.
#'
#' @param expression_vectors a numeric matrix. The input matrix of all gene expression vectors (n_axes x n_genes).
#' @param gene_to_family a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0_int32` for unassigned genes
#' @param n_families a integer scalar. Total number of gene families to compute centroids for.
#' @param mode a string, one of "group_orthologs", "group_all". used mode for grouping
#' @param ortholog_set a logical vector. A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
#' @return The output matrix (n_axes x n_families) to store the computed centroids.
#'
#' Generated from the Fortran procedure \code{tox_gene_centroids::group_centroid}.
#' @export
group_centroid <- function(expression_vectors, gene_to_family, n_families, mode, ortholog_set = NULL) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_family <- .tox_as_integer_vector(gene_to_family, "gene_to_family")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    mode <- .tox_as_mode(mode, "mode", c("group_orthologs", "group_all"))
    if (!is.null(ortholog_set))
        ortholog_set <- .tox_as_logical(ortholog_set, "ortholog_set")
    if (length(gene_to_family) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_family", length(gene_to_family), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("group_centroid_call", expression_vectors, gene_to_family, n_families, mode, ortholog_set)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "mode", "tmp_group_indices", "ierr", "ortholog_set")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$centroid_matrix
}
