# Generated. Do not edit.

#' Read expression vectors from csv/tsv files
#'
#' @param file_list a character vector. List of files to read from
#' @param gene_ids a character vector. Array of gene IDS
#' @param expression_vectors a numeric matrix. Array of expression vectors
#' @param n_header_rows a integer scalar. Number of header rows to skip
#' @param gene_col a integer scalar. Index of column with gene_ids
#' @param value_cols a integer vector. Indicies of columns containing values
#' @param start_row a integer scalar. Row in the expression vectors to start in
#' @param delimiter a string. optional delimiter
#'   The default value is `char(9)`.
#' @return Array of expression vectors
#'
#' Generated from the Fortran procedure \code{tox_data_tools::read_expression_vectors_tsv}.
#' @export
read_expression_vectors_tsv <- function(file_list, gene_ids, expression_vectors, n_header_rows, gene_col, value_cols, start_row, delimiter = "	") {
    file_list <- .tox_as_character(file_list, "file_list")
    gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    expression_vectors <- .tox_as_double_matrix(expression_vectors, "expression_vectors")
    n_header_rows <- .tox_as_integer_scalar(n_header_rows, "n_header_rows")
    gene_col <- .tox_as_integer_scalar(gene_col, "gene_col")
    value_cols <- .tox_as_integer_vector(value_cols, "value_cols")
    start_row <- .tox_as_integer_scalar(start_row, "start_row")
    delimiter <- .tox_as_character(delimiter, "delimiter")
    .result <- .Call("read_expression_vectors_tsv_call", file_list, gene_ids, expression_vectors, n_header_rows, gene_col, value_cols, start_row, delimiter)
    .arguments <- c("file_list", "gene_ids", "expression_vectors", "n_header_rows", "gene_col", "value_cols", "start_row", "ierr", "delimiter")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$expression_vectors
}

#' Only read the gene ids from a tsv file
#'
#' @param filename a string. Name of the file
#' @param gene_ids_strlen a integer scalar. length of the strings in `gene_ids`
#' @param n_gene_ids_elements a integer scalar. number of elements in `gene_ids`
#' @param n_header_rows a integer scalar. number of headers to skip
#' @param gene_col a integer scalar. Index of the column containing gene ids
#' @return gene ids array
#'
#' Generated from the Fortran procedure \code{tox_data_tools::read_gene_ids_from_tsv_file}.
#' @export
read_gene_ids_from_tsv_file <- function(filename, gene_ids_strlen, n_gene_ids_elements, n_header_rows, gene_col) {
    filename <- .tox_as_character(filename, "filename")
    gene_ids_strlen <- .tox_as_integer_scalar(gene_ids_strlen, "gene_ids_strlen")
    n_gene_ids_elements <- .tox_as_integer_scalar(n_gene_ids_elements, "n_gene_ids_elements")
    n_header_rows <- .tox_as_integer_scalar(n_header_rows, "n_header_rows")
    gene_col <- .tox_as_integer_scalar(gene_col, "gene_col")
    .result <- .Call("read_gene_ids_from_tsv_file_call", filename, gene_ids_strlen, n_gene_ids_elements, n_header_rows, gene_col)
    .arguments <- c("filename", "gene_ids", "n_header_rows", "gene_col", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$gene_ids
}

#' Read a family file (Orthofinder)
#'
#' @param filename a string. Name of the file
#' @param gene_ids a character vector. gene ids array
#' @param family_ids_strlen a integer scalar. length of the strings in `family_ids`
#' @param n_family_ids_elements a integer scalar. number of elements in `family_ids`
#' @param n_gene_to_fam_elements a integer scalar. number of elements in `gene_to_fam`
#' @return a named list with elements `family_ids`, `gene_to_fam`.
#'
#' Generated from the Fortran procedure \code{tox_data_tools::read_orthofinder_file}.
#' @export
read_orthofinder_file <- function(filename, gene_ids, family_ids_strlen, n_family_ids_elements, n_gene_to_fam_elements) {
    filename <- .tox_as_character(filename, "filename")
    gene_ids <- .tox_as_character(gene_ids, "gene_ids")
    family_ids_strlen <- .tox_as_integer_scalar(family_ids_strlen, "family_ids_strlen")
    n_family_ids_elements <- .tox_as_integer_scalar(n_family_ids_elements, "n_family_ids_elements")
    n_gene_to_fam_elements <- .tox_as_integer_scalar(n_gene_to_fam_elements, "n_gene_to_fam_elements")
    .result <- .Call("read_orthofinder_file_call", filename, gene_ids, family_ids_strlen, n_family_ids_elements, n_gene_to_fam_elements)
    .arguments <- c("filename", "gene_ids", "family_ids", "gene_to_fam", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        family_ids = .result$family_ids,
        gene_to_fam = .result$gene_to_fam
    )
}

#' Helper to create a mask of genes that are unassigned
#'
#' @param gene_to_fam a integer vector. gene to family mapping
#' @return a named list with elements `mask`, `n_genes_kept`.
#'
#' Generated from the Fortran procedure \code{tox_data_tools::get_unassigned_mask}.
#' @export
get_unassigned_mask <- function(gene_to_fam) {
    gene_to_fam <- .tox_as_integer_vector(gene_to_fam, "gene_to_fam")
    .result <- .Call("get_unassigned_mask_call", gene_to_fam)
    .arguments <- c("gene_to_fam", "mask", "n_genes_kept")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        mask = .result$mask,
        n_genes_kept = .result$n_genes_kept
    )
}
