# Generated. Do not edit.

#' Creates a zip archive with generic file lists
#'
#' Generated from the Fortran procedure \code{tox_data_archive::create_zip_archive}, whose argument names
#' are the ones an error message reports.
#'
#' @param zip_filename a string. Name of the zip file to create
#' @param keys a character vector. Array of keys for manifest entries
#' @param filenames a character vector. Array of filenames to add to zip
#' @return invisibly `NULL`; called for its effect.
#' @export
create_zip_archive <- function(zip_filename, keys, filenames) {
    zip_filename <- .tox_as_character(zip_filename, "zip_filename")
    keys <- .tox_as_character(keys, "keys")
    filenames <- .tox_as_character(filenames, "filenames")
    .result <- .Call("create_zip_archive_call", zip_filename, keys, filenames)
    .arguments <- c("zip_filename", "keys", "filenames", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Save standard tox data
#'
#' Generated from the Fortran procedure \code{tox_data_archive::save_tox_data}, whose argument names
#' are the ones an error message reports.
#'
#' @param zip_filename a string. Zip filename
#' @param gene_ids a character vector. Gene ids array, will be saved if provided
#' @param gene_ids_file a string. Name of the gene ids file
#' @param expression a numeric matrix. Expression vectors array, will be saved if provided
#' @param expression_file a string. Name of the expression file
#' @param gene_to_family a integer vector. Gene to family mapping array, will be saved if provided
#' @param gene_to_family_file a string. Name of the gene to family mapping file
#' @param family_ids a character vector. Family ids array, will be saved if provided
#' @param family_ids_file a string. Name of the family ids file
#' @param family_centroids a numeric matrix. Family centroids array, will be saved if provided
#' @param family_centroids_file a string. Name of the family centroids file
#' @param shift_vectors a numeric matrix. Shift vectors array, will be saved if provided
#' @param shift_vectors_file a string. Name of the shift vectors file
#' @return invisibly `NULL`; called for its effect.
#' @export
save_tox_data <- function(zip_filename, gene_ids = NULL, gene_ids_file = NULL, expression = NULL, expression_file = NULL, gene_to_family = NULL, gene_to_family_file = NULL, family_ids = NULL, family_ids_file = NULL, family_centroids = NULL, family_centroids_file = NULL, shift_vectors = NULL, shift_vectors_file = NULL) {
    zip_filename <- .tox_as_character(zip_filename, "zip_filename")
    if (!is.null(gene_ids))
        gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    if (!is.null(gene_ids_file))
        gene_ids_file <- .tox_as_character(gene_ids_file, "gene_ids_file")
    if (!is.null(expression))
        expression <- .tox_as_double_matrix(expression, "expression")
    if (!is.null(expression_file))
        expression_file <- .tox_as_character(expression_file, "expression_file")
    if (!is.null(gene_to_family))
        gene_to_family <- .tox_as_integer_vector(gene_to_family, "gene_to_family")
    if (!is.null(gene_to_family_file))
        gene_to_family_file <- .tox_as_character(gene_to_family_file, "gene_to_family_file")
    if (!is.null(family_ids))
        family_ids <- .tox_as_character(family_ids, "family_ids")
    if (!is.null(family_ids_file))
        family_ids_file <- .tox_as_character(family_ids_file, "family_ids_file")
    if (!is.null(family_centroids))
        family_centroids <- .tox_as_double_matrix(family_centroids, "family_centroids")
    if (!is.null(family_centroids_file))
        family_centroids_file <- .tox_as_character(family_centroids_file, "family_centroids_file")
    if (!is.null(shift_vectors))
        shift_vectors <- .tox_as_double_matrix(shift_vectors, "shift_vectors")
    if (!is.null(shift_vectors_file))
        shift_vectors_file <- .tox_as_character(shift_vectors_file, "shift_vectors_file")
    .result <- .Call("save_tox_data_call", zip_filename, gene_ids, gene_ids_file, expression, expression_file, gene_to_family, gene_to_family_file, family_ids, family_ids_file, family_centroids, family_centroids_file, shift_vectors, shift_vectors_file)
    .arguments <- c("zip_filename", "ierr", "gene_ids", "gene_ids_file", "expression", "expression_file", "gene_to_family", "gene_to_family_file", "family_ids", "family_ids_file", "family_centroids", "family_centroids_file", "shift_vectors", "shift_vectors_file")
    .status <- check_err_code(.result$ierr, .arguments)

    invisible(NULL)
}

#' Report the shape of every member of a tox data archive
#'
#' Each count (and each string length) is 0 when the corresponding member is absent, so a
#' caller can size all six output buffers up front. Character members report both an element
#' count and a per-element string length. Pairs with
#' [[tox_data_archive(module):read_tox_data_into(subroutine)]].
#'
#' Generated from the Fortran procedure \code{tox_data_archive::get_tox_data_dims}, whose argument names
#' are the ones an error message reports.
#'
#' @param zip_filename a string. Name of the zip file
#' @return a named list with elements:
#'   \item{n_gene_ids}{a integer scalar. Number of gene ids, 0 if absent}
#'   \item{gene_id_len}{a integer scalar. String length of each gene id, 0 if absent}
#'   \item{n_expression_rows}{a integer scalar. Rows (samples) of the expression matrix, 0 if absent}
#'   \item{n_expression_cols}{a integer scalar. Columns (genes) of the expression matrix, 0 if absent}
#'   \item{n_gene_to_family}{a integer scalar. Number of gene-to-family entries, 0 if absent}
#'   \item{n_family_ids}{a integer scalar. Number of family ids, 0 if absent}
#'   \item{family_id_len}{a integer scalar. String length of each family id, 0 if absent}
#'   \item{n_family_centroids_rows}{a integer scalar. Rows (samples) of the family centroids matrix, 0 if absent}
#'   \item{n_family_centroids_cols}{a integer scalar. Columns (families) of the family centroids matrix, 0 if absent}
#'   \item{n_shift_vectors_rows}{a integer scalar. Rows of the shift vectors matrix, 0 if absent}
#'   \item{n_shift_vectors_cols}{a integer scalar. Columns of the shift vectors matrix, 0 if absent}
#' @export
get_tox_data_dims <- function(zip_filename) {
    zip_filename <- .tox_as_character(zip_filename, "zip_filename")
    .result <- .Call("get_tox_data_dims_call", zip_filename)
    .arguments <- c("zip_filename", "n_gene_ids", "gene_id_len", "n_expression_rows", "n_expression_cols", "n_gene_to_family", "n_family_ids", "family_id_len", "n_family_centroids_rows", "n_family_centroids_cols", "n_shift_vectors_rows", "n_shift_vectors_cols", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        n_gene_ids = .result$n_gene_ids,
        gene_id_len = .result$gene_id_len,
        n_expression_rows = .result$n_expression_rows,
        n_expression_cols = .result$n_expression_cols,
        n_gene_to_family = .result$n_gene_to_family,
        n_family_ids = .result$n_family_ids,
        family_id_len = .result$family_id_len,
        n_family_centroids_rows = .result$n_family_centroids_rows,
        n_family_centroids_cols = .result$n_family_centroids_cols,
        n_shift_vectors_rows = .result$n_shift_vectors_rows,
        n_shift_vectors_cols = .result$n_shift_vectors_cols
    )
}

#' Read a tox data archive into caller-provided buffers
#'
#' Fills every buffer from the archive; size them from
#' [[tox_data_archive(module):get_tox_data_dims(subroutine)]] first. A member that is absent
#' has a zero extent and is left untouched.
#'
#' Generated from the Fortran procedure \code{tox_data_archive::read_tox_data_into}, whose argument names
#' are the ones an error message reports.
#'
#' @param zip_filename a string. Name of the zip file
#' @return a named list with elements:
#'   \item{gene_ids}{a character vector. Gene ids}
#'   \item{expression}{a numeric matrix. Expression vectors}
#'   \item{gene_to_family}{a integer vector. Gene to family mapping}
#'   \item{family_ids}{a character vector. Family ids}
#'   \item{family_centroids}{a numeric matrix. Family centroids}
#'   \item{shift_vectors}{a numeric matrix. Shift vectors}
#' @export
read_tox_data_into <- function(zip_filename) {
    zip_filename <- .tox_as_character(zip_filename, "zip_filename")
    .get_tox_data_dims_result <- get_tox_data_dims(zip_filename = zip_filename)
    n_gene_ids <- .get_tox_data_dims_result$n_gene_ids
    gene_id_len <- .get_tox_data_dims_result$gene_id_len
    n_expression_rows <- .get_tox_data_dims_result$n_expression_rows
    n_expression_cols <- .get_tox_data_dims_result$n_expression_cols
    n_gene_to_family <- .get_tox_data_dims_result$n_gene_to_family
    n_family_ids <- .get_tox_data_dims_result$n_family_ids
    family_id_len <- .get_tox_data_dims_result$family_id_len
    n_family_centroids_rows <- .get_tox_data_dims_result$n_family_centroids_rows
    n_family_centroids_cols <- .get_tox_data_dims_result$n_family_centroids_cols
    n_shift_vectors_rows <- .get_tox_data_dims_result$n_shift_vectors_rows
    n_shift_vectors_cols <- .get_tox_data_dims_result$n_shift_vectors_cols

    .result <- .Call("read_tox_data_into_call", zip_filename, n_gene_ids, gene_id_len, n_expression_rows, n_expression_cols, n_gene_to_family, n_family_ids, family_id_len, n_family_centroids_rows, n_family_centroids_cols, n_shift_vectors_rows, n_shift_vectors_cols)
    .arguments <- c("zip_filename", "n_gene_ids", "gene_id_len", "gene_ids", "n_expression_rows", "n_expression_cols", "expression", "n_gene_to_family", "gene_to_family", "n_family_ids", "family_id_len", "family_ids", "n_family_centroids_rows", "n_family_centroids_cols", "family_centroids", "n_shift_vectors_rows", "n_shift_vectors_cols", "shift_vectors", "ierr")
    .sources <- c(NA_character_, "gene_ids", NA_character_, NA_character_, "expression", "expression", NA_character_, "gene_to_family", NA_character_, "family_ids", NA_character_, NA_character_, "family_centroids", "family_centroids", NA_character_, "shift_vectors", "shift_vectors", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        gene_ids = .result$gene_ids,
        expression = .result$expression,
        gene_to_family = .result$gene_to_family,
        family_ids = .result$family_ids,
        family_centroids = .result$family_centroids,
        shift_vectors = .result$shift_vectors
    )
}
