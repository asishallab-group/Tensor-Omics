# Generated. Do not edit.

#' Checks the state of a bit/paralog in `bit_mask` -> TRUE if 1 else FALSE
#'
#' @param bit_mask a integer vector. chunked mask to mark active paralogs
#' @param i_gene a integer scalar. index of paralog to be marked active
#' @return check result
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis_kernel::mask_check_state}.
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
#' @param n_genes a integer scalar. number of genes
#' @return number of 32 bit chunks a mask needs to encode `n_genes` genes
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis_kernel::mask_chunk_count}.
#' @export
mask_chunk_count <- function(n_genes) {
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    .result <- .Call("mask_chunk_count_call", n_genes)
    .arguments <- c("n_genes", "count")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$count
}

#' Calculates the needed size for the paralog-subsets work array
#'
#' The `detect_*` subroutines need a work array for the to be tested subsets.
#' In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
#' This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
#'
#' This subroutine calculates the needed size for the work array.
#'
#' @param max_subset_size a integer scalar. maximum size that a subset must not exceed. Zero is in range and means there is
#'   nothing to size a work array for, which is reported back as a size of zero.
#'   The minimum valid value is `0`.
#'   If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.
#'
#'   Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
#' @param n_genes a integer scalar. number of genes
#' @param filtered_paralogs_mask a integer vector. Output mask with all genes disabled that did not pass the filter
#' @return a named list with elements `max_subset_size`, `work_array_size`.
#'
#' Generated from the Fortran procedure \code{tox_paralog_analysis_kernel::calc_work_arr_paralog_subsets_size}.
#' @export
calc_work_arr_paralog_subsets_size <- function(max_subset_size, n_genes, filtered_paralogs_mask) {
    max_subset_size <- .tox_as_integer_scalar(max_subset_size, "max_subset_size")
    n_genes <- .tox_as_integer_scalar(n_genes, "n_genes")
    filtered_paralogs_mask <- .tox_as_integer_vector(filtered_paralogs_mask, "filtered_paralogs_mask")
    .result <- .Call("calc_work_arr_paralog_subsets_size_call", max_subset_size, n_genes, filtered_paralogs_mask)
    .arguments <- c("max_subset_size", "n_genes", "work_array_size", "filtered_paralogs_mask", "n_mask_chunks", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_subset_size = .result$max_subset_size,
        work_array_size = .result$work_array_size
    )
}
