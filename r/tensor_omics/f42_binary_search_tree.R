# Generated. Do not edit.

#' Build the BST index by sorting indices using values in x
#'
#' Generated from the Fortran procedure \code{f42_binary_search_tree::build_bst_index}, whose argument names
#' are the ones an error message reports.
#'
#' @param values a numeric vector. Input real array to be indexed
#' @return a integer vector. Output permutation index
#' @export
build_bst_index <- function(values) {
    values <- .tox_as_double_vector(values, "values")
    .result <- .Call("build_bst_index_call", values)
    .arguments <- c("values", "n_values", "sorted_indices", "ierr")
    .sources <- c(NA_character_, "values", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$sorted_indices
}

#' Perform a 1D range query over the sorted index
#'
#' Generated from the Fortran procedure \code{f42_binary_search_tree::bst_range_query}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{values_perm} and sorts it by \code{values}.
#' Call \code{bst_range_query_expert} to do that yourself.
#'
#' @param values a numeric vector. Input real array
#' @param lower_bound a numeric scalar. Lower bound of range (inclusive)
#'   The maximum valid value is `upper_bound`.
#' @param upper_bound a numeric scalar. Upper bound of range (inclusive)
#' @return a integer vector. Output array of matching indices.
#'   The first `n_matches` elements will hold the results.
#' @export
bst_range_query <- function(values, lower_bound, upper_bound) {
    values <- .tox_as_double_vector(values, "values")
    lower_bound <- .tox_as_double_scalar(lower_bound, "lower_bound")
    upper_bound <- .tox_as_double_scalar(upper_bound, "upper_bound")
    .result <- .Call("bst_range_query_call", values, lower_bound, upper_bound)
    .arguments <- c("values", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr")
    .sources <- c(NA_character_, "values", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    utils::head(.result$output_indices, .result$n_matches)
}

#' Perform a 1D range query over the sorted index
#'
#' Generated from the Fortran procedure \code{f42_binary_search_tree::bst_range_query_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{values_perm} yourself.
#' \code{bst_range_query} seeds \code{values_perm} and sorts it by \code{values}.
#'
#' @param values a numeric vector. Input real array
#' @param values_perm a integer vector. Permutation of `values` in ascending order -- the BST index. The allocating entry
#'   point builds and heapsorts it for you; the expert one takes whatever order you
#'   supply, so a caller that already holds one from
#'   \code{\link{build_bst_index}} can reuse it across queries.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_values`.
#' @param lower_bound a numeric scalar. Lower bound of range (inclusive)
#'   The maximum valid value is `upper_bound`.
#' @param upper_bound a numeric scalar. Upper bound of range (inclusive)
#' @return a integer vector. Output array of matching indices.
#'   The first `n_matches` elements will hold the results.
#' @export
bst_range_query_expert <- function(values, values_perm, lower_bound, upper_bound) {
    values <- .tox_as_double_vector(values, "values")
    values_perm <- .tox_as_integer_vector(values_perm, "values_perm")
    lower_bound <- .tox_as_double_scalar(lower_bound, "lower_bound")
    upper_bound <- .tox_as_double_scalar(upper_bound, "upper_bound")
    if (length(values_perm) != length(values))
        .tox_shape_error("values_perm", length(values_perm), "values", length(values))

    .result <- .Call("bst_range_query_expert_call", values, values_perm, lower_bound, upper_bound)
    .arguments <- c("values", "values_perm", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr")
    .sources <- c(NA_character_, NA_character_, "values", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    utils::head(.result$output_indices, .result$n_matches)
}
