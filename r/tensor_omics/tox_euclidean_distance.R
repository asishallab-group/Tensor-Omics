# Generated. Do not edit.

#' Compute the Euclidean distance between two vectors.
#'
#' Calculates the L2 norm: `result = sqrt(sum((vec1_i - vec2_i)**2))`
#'
#' Generated from the Fortran procedure \code{tox_euclidean_distance::euclidean_distance}, whose argument names
#' are the ones an error message reports.
#'
#' @param vec1 a numeric vector. First expression vector
#' @param vec2 a numeric vector. Second expression vector
#' @return a numeric scalar. Output scalar distance
#' @export
euclidean_distance <- function(vec1, vec2) {
    vec1 <- .tox_as_double_vector(vec1, "vec1")
    vec2 <- .tox_as_double_vector(vec2, "vec2")
    if (length(vec2) != length(vec1))
        .tox_shape_error("vec2", length(vec2), "vec1", length(vec1))

    .result <- .Call("euclidean_distance_call", vec1, vec2)
    .arguments <- c("vec1", "vec2", "n_elements", "result", "ierr")
    .sources <- c(NA_character_, NA_character_, "vec1", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$result
}

#' Compute distance from each gene to its corresponding family centroid.
#'
#' For each gene, extracts its expression vector and the centroid of its assigned family, then computes the Euclidean distance between them.
#'
#' Generated from the Fortran procedure \code{tox_euclidean_distance::distance_to_centroid}, whose argument names
#' are the ones an error message reports.
#'
#' @param genes a numeric matrix. Gene expression matrix (n_tissues × n_genes), column-major
#' @param centroids a numeric matrix. Family centroid matrix (n_tissues × n_families), column-major
#' @param gene_to_fam a integer vector. Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0` for unassigned genes
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_families`.
#'   The value `0` is additionally accepted.
#' @return a numeric vector. Euclidean distance from each gene to its own family's centroid. A gene carrying
#'   `0` in `gene_to_fam` has no centroid to measure against and
#'   receives `-1.0` instead of a distance -- a value no Euclidean
#'   distance can take, so that summarising this array without excluding those genes is
#'   visibly wrong rather than quietly biased. Which genes those are is `gene_to_fam`,
#'   not a negative entry here.
#' @export
distance_to_centroid <- function(genes, centroids, gene_to_fam) {
    genes <- .tox_as_double_matrix(genes, "genes")
    centroids <- .tox_as_double_matrix(centroids, "centroids")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    if (length(gene_to_fam) != dim(genes)[2])
        .tox_shape_error("gene_to_fam", length(gene_to_fam), "genes", dim(genes)[2])
    if (dim(centroids)[1] != dim(genes)[1])
        .tox_shape_error("centroids", dim(centroids)[1], "genes", dim(genes)[1])

    .result <- .Call("distance_to_centroid_call", genes, centroids, gene_to_fam)
    .arguments <- c("n_genes", "n_families", "genes", "centroids", "gene_to_fam", "distances", "n_tissues", "ierr")
    .sources <- c("genes", "centroids", NA_character_, NA_character_, NA_character_, NA_character_, "genes", NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$distances
}
