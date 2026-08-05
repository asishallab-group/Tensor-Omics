# Generated. Do not edit.

#' Validate full data structure
#'
#' @param n_genes a integer scalar. Expected number of genes
#' @param n_families a integer scalar. Expected number of families
#' @param n_samples a integer scalar. Expected number of samples
#' @param gene_ids a character vector. Gene ids
#' @param gene_family_ids a character vector. Gene family ids
#' @param gene_to_fam a integer vector. gene to family mapping
#' @param expression_vectors a numeric matrix. Expression vectors
#' @param family_centroids a numeric matrix. Family centroids
#' @param shift_vectors a numeric matrix. Shift vectors
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_data_structure <- function(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors) {
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    n_samples <- .tox_as_integer_scalar(n_samples, "n_samples")
    gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    gene_family_ids <- .tox_as_character(gene_family_ids, "gene_family_ids")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    shift_vectors <- .tox_as_double_matrix(shift_vectors, "shift_vectors")
    .result <- .Call("validate_data_structure_call", n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors)
    .arguments <- c("n_genes", "n_families", "n_samples", "gene_ids", "gene_family_ids", "gene_to_fam", "expression_vectors", "family_centroids", "shift_vectors", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Validate gene to family mapping
#'
#' @param gene_to_fam a integer vector. gene to family mapping
#' @param n_families a integer scalar. number of families
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_gene_to_family_mapping <- function(gene_to_fam, n_families) {
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    .result <- .Call("validate_gene_to_family_mapping_call", gene_to_fam, n_families)
    .arguments <- c("gene_to_fam", "n_families", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Validate expresssion data
#'
#' @param expression_vectors a numeric matrix. Expression vectors
#' @param check_non_negative a logical scalar. Defines if non negative should be checked
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_expression_data <- function(expression_vectors, check_non_negative) {
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    check_non_negative <- .tox_as_logical(check_non_negative, "check_non_negative")
    .result <- .Call("validate_expression_data_call", expression_vectors, check_non_negative)
    .arguments <- c("expression_vectors", "check_non_negative", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Validate the family centroids
#'
#' @param family_centroids a numeric matrix. Family centroids array
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_family_centroids <- function(family_centroids) {
    family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    .result <- .Call("validate_family_centroids_call", family_centroids)
    .arguments <- c("family_centroids", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Validates shift vectors
#'
#' @param shift_vectors a numeric matrix. shift vectors
#' @param expression_vectors a numeric matrix. expression vectors
#' @param family_centroids a numeric matrix. family centroids
#' @param gene_to_fam a integer vector. gene to family mapping
#' @param n_samples a integer scalar. Number of samples
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_shift_vectors <- function(shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples) {
    shift_vectors <- .tox_as_double_matrix(shift_vectors, "shift_vectors")
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    n_samples <- .tox_as_integer_scalar(n_samples, "n_samples")
    .result <- .Call("validate_shift_vectors_call", shift_vectors, expression_vectors, family_centroids, gene_to_fam, n_samples)
    .arguments <- c("shift_vectors", "expression_vectors", "family_centroids", "gene_to_fam", "n_samples", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Validate that no string appears more than once
#'
#' @param str_arr a character vector. string array
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_string_array_uniqueness <- function(str_arr) {
    str_arr <- .tox_as_character(str_arr, "str_arr")
    .result <- .Call("validate_string_array_uniqueness_call", str_arr)
    .arguments <- c("str_arr", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Comprehensive validation routine, combining all checks
#'
#' @param n_genes a integer scalar. Number of genes
#' @param n_families a integer scalar. Number of families
#' @param n_samples a integer scalar. Number of samples
#' @param gene_ids a character vector. Gene ids array
#' @param gene_family_ids a character vector. gene family ids
#' @param gene_to_fam a integer vector. gene to family mapping
#' @param expression_vectors a numeric matrix. Expression vectors
#' @param family_centroids a numeric matrix. family centroids
#' @param shift_vectors a numeric matrix. shift vectors
#' @param check_uniqueness a logical scalar. Check ID arrays for uniqueness.
#'   The default value is `TRUE`.
#' @param check_shift_consistency a logical scalar. Check consitency of shift array.
#'   The default value is `TRUE`.
#' @return invisibly `NULL`; called for its effect.
#'
#' Generated from the Fortran module \code{tox_data_validation}.
#' @export
validate_all_data <- function(n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors, check_uniqueness = TRUE, check_shift_consistency = TRUE) {
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    n_families <- .tox_as_integer_scalar(n_families, "n_families")
    n_samples <- .tox_as_integer_scalar(n_samples, "n_samples")
    gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    gene_family_ids <- .tox_as_character(gene_family_ids, "gene_family_ids")
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    shift_vectors <- .tox_as_double_matrix(shift_vectors, "shift_vectors")
    check_uniqueness <- .tox_as_logical(check_uniqueness, "check_uniqueness")
    check_shift_consistency <- .tox_as_logical(check_shift_consistency, "check_shift_consistency")
    .result <- .Call("validate_all_data_call", n_genes, n_families, n_samples, gene_ids, gene_family_ids, gene_to_fam, expression_vectors, family_centroids, shift_vectors, check_uniqueness, check_shift_consistency)
    .arguments <- c("n_genes", "n_families", "n_samples", "gene_ids", "gene_family_ids", "gene_to_fam", "expression_vectors", "family_centroids", "shift_vectors", "ierr", "check_uniqueness", "check_shift_consistency")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}
