# Generated. Do not edit.

#' Compute the shift vector field for all genes.
#'
#' Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
#'
#' @param expression_vectors a numeric matrix. Gene expression matrix
#' @param family_centroids a numeric matrix. Family centroid matrix
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#' @return Output, real matrix array, stores the centroid of the gene's family in `shift_vectors(:, 1, i_gene)` (zero vector if no family assigned) and the shift vectors in `shift_vectors(:, 2, i_gene)`
#'
#' Generated from the Fortran procedure \code{tox_shift_vectors::compute_shift_vector_field}.
#' @export
compute_shift_vector_field <- function(expression_vectors, family_centroids, gene_to_fam) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    if (dim(family_centroids)[1] != dim(expression_vectors)[1])
        .tox_shape_error("family_centroids", dim(family_centroids)[1], "expression_vectors", dim(expression_vectors)[1])
    if (length(gene_to_fam) != dim(expression_vectors)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "expression_vectors", dim(expression_vectors)[2])

    .result <- .Call("compute_shift_vector_field_call", expression_vectors, family_centroids, gene_to_fam)
    .arguments <- c("n_tissues", "n_genes", "n_families", "expression_vectors", "family_centroids", "gene_to_fam", "shift_vectors", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$shift_vectors
}
