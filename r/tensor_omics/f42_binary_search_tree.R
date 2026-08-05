# Generated. Do not edit.

#' Build the BST index by sorting indices using values in x
#'
#' @param values a numeric vector. Input real array to be indexed
#' @return Output permutation index
#'
#' Generated from the Fortran module \code{f42_binary_search_tree}.
#' @export
build_bst_index <- function(values) {
    values <- .tox_as_double_vector(values, "values")
    .result <- .Call("build_bst_index_call", values)
    .arguments <- c("values", "n_values", "sorted_indices", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$sorted_indices
}

#' Perform a 1D range query over the sorted index
#'
#' @param values a numeric vector. Input real array
#' @param sorted_indices a integer vector. Permutation index array (sorted)
#' @param lower_bound a numeric scalar. Lower bound of range (inclusive)
#' @param upper_bound a numeric scalar. Upper bound of range (inclusive)
#' @return Output array of matching indices.
#'
#' Generated from the Fortran module \code{f42_binary_search_tree}.
#' @export
bst_range_query <- function(values, sorted_indices, lower_bound, upper_bound) {
    values <- .tox_as_double_vector(values, "values")
    sorted_indices <- .tox_as_integer_vector(sorted_indices, "sorted_indices")
    lower_bound <- .tox_as_double_scalar(lower_bound, "lower_bound")
    upper_bound <- .tox_as_double_scalar(upper_bound, "upper_bound")
    if (length(sorted_indices) != length(values))
        .tox_shape_error("sorted_indices", length(sorted_indices), "values", length(values))

    .result <- .Call("bst_range_query_call", values, sorted_indices, lower_bound, upper_bound)
    .arguments <- c("values", "sorted_indices", "n_values", "lower_bound", "upper_bound", "output_indices", "n_matches", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    utils::head(.result$output_indices, .result$n_matches)
}
