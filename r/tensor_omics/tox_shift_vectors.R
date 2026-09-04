# Generated. Do not edit.

#' Compute the shift vector field for all genes.
#'
#' Computes the shift vectors by subtracting the corresponding family centroid from the expression vector.
#'
#' Generated from the Fortran procedure \code{tox_shift_vectors::compute_shift_vector_field}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression_vectors a numeric matrix. Gene expression matrix
#' @param family_centroids a numeric matrix. Family centroid matrix
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @return a numeric array of rank 3. Output, real matrix array. For each gene it holds two vectors: the centroid of the gene's family first (a zero vector if no family is assigned), then the shift vector
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
    .sources <- c("expression_vectors", "expression_vectors", "family_centroids", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$shift_vectors
}
