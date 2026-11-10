# Generated. Do not edit.

#' Checks the state of a bit/paralog in `bit_mask` -> TRUE if 1 else FALSE
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis_impl::mask_check_state}, whose argument names
#' are the ones an error message reports.
#'
#' @param bit_mask a integer vector. chunked mask to mark active paralogs
#' @param i_gene a integer scalar. index of paralog to be marked active
#' @return a logical scalar. check result
#' @export
mask_check_state <- function(bit_mask, i_gene) {
    bit_mask <- .tox_as_integer_vector(bit_mask, "bit_mask")
    i_gene <- .tox_as_integer_scalar(i_gene, "i_gene")
    .result <- .Call("mask_check_state_call", bit_mask, i_gene)
    .arguments <- c("bit_mask", "i_gene")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$state
}

#' Determines the needed chunk count for subset bit masks (an integer has only 32 bits)
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis_impl::mask_chunk_count}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_genes a integer scalar. number of genes
#' @return a integer scalar. number of 32 bit chunks a mask needs to encode `n_genes` genes
#'
#'   Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32.0)` and represents the number of chunks
#' @export
mask_chunk_count <- function(n_genes) {
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    .result <- .Call("mask_chunk_count_call", n_genes)
    .arguments <- c("n_genes", "count")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$count
}
