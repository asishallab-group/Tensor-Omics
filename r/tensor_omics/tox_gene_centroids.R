# Generated. Do not edit.

#' Computes the element-wise mean of the selected gene vectors.
#'
#' A selection without genes gives the zero vector.
#'
#' Generated from the Fortran procedure \code{tox_gene_centroids::mean_vector}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. The input matrix of all gene expression vectors (n_axes x n_genes).
#' @param genes_selection_mask a logical vector. `TRUE` for the genes (columns of `expression_vectors`) to average
#' @return a numeric vector. The output vector representing the computed centroid.
#' @export
mean_vector <- function(expression_vectors, genes_selection_mask) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    genes_selection_mask <- .tox_as_logical_vector(genes_selection_mask, "genes_selection_mask")
    if (length(genes_selection_mask) != dim(expression_vectors)[2])
        .tox_shape_error("genes_selection_mask", length(genes_selection_mask), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("mean_vector_call", expression_vectors, genes_selection_mask)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "genes_selection_mask", "n_selected_genes", "centroid", "ierr")
    .sources <- c(NA_character_, "expression_vectors", "expression_vectors", NA_character_, "genes_selection_mask", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$centroid
}

#' Computes one centroid per gene family, of all its genes or of its orthologs only.
#'
#' A family without selected genes gets the zero vector.
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
#' @param ortholog_set a logical vector. A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
#' @return a numeric matrix. The output matrix (n_axes x n_families) to store the computed centroids.
#' @export
group_centroid_orthologs <- function(expression_vectors, gene_to_family, n_families, ortholog_set) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_family <- .tox_as_integer_vector(gene_to_family, "gene_to_family")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    ortholog_set <- .tox_as_logical_vector(ortholog_set, "ortholog_set")
    if (length(gene_to_family) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_family", length(gene_to_family), "expression_vectors", dim(expression_vectors)[2])
    if (length(ortholog_set) != dim(expression_vectors)[2])
        .tox_shape_error("ortholog_set", length(ortholog_set), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("group_centroid_orthologs_call", expression_vectors, gene_to_family, n_families, ortholog_set)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "ortholog_set", "ierr")
    .sources <- c(NA_character_, "expression_vectors", "expression_vectors", NA_character_, "centroid_matrix", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$centroid_matrix
}

#' Computes one centroid per gene family, of all its genes or of its orthologs only.
#'
#' A family without selected genes gets the zero vector.
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
#' @return a numeric matrix. The output matrix (n_axes x n_families) to store the computed centroids.
#' @export
group_centroid_all <- function(expression_vectors, gene_to_family, n_families) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    gene_to_family <- .tox_as_integer_vector(gene_to_family, "gene_to_family")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    if (length(gene_to_family) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_family", length(gene_to_family), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("group_centroid_all_call", expression_vectors, gene_to_family, n_families)
    .arguments <- c("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "ierr")
    .sources <- c(NA_character_, "expression_vectors", "expression_vectors", NA_character_, "centroid_matrix", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$centroid_matrix
}
