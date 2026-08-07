# Generated. Do not edit.

#' Compute the Empirical Distribution Function (EDF) from a sorted permutation
#'
#' Returns the sorted unique values and their cumulative frequencies in [0,1].
#' The number of unique values can be determined by finding the last non-zero cdf_value.
#'
#' Generated from the Fortran procedure \code{f42_stats::compute_edf}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{values_perm} and sorts it by \code{values}.
#' Call \code{compute_edf_expert} to do that yourself.
#'
#' @param values a numeric vector. Array of observed data values (e.g., contributions or spikes).
#' @return a named list with elements:
#'   \item{unique_values}{a numeric vector. Sorted unique data values.
#'     The first `n_unique` elements will hold the results.}
#'   \item{cdf_values}{a numeric vector. Corresponding cumulative frequencies between 0 and 1.
#'     The first `n_unique` elements will hold the results.}
#' @export
compute_edf <- function(values) {
    values <- .tox_as_double_vector(values, "values")
    .result <- .Call("compute_edf_call", values)
    .arguments <- c("values", "n_values", "unique_values", "cdf_values", "n_unique", "ierr")
    .sources <- c(NA_character_, "values", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        unique_values = utils::head(.result$unique_values, .result$n_unique),
        cdf_values = utils::head(.result$cdf_values, .result$n_unique)
    )
}

#' Compute the Empirical Distribution Function (EDF) from a sorted permutation
#'
#' Returns the sorted unique values and their cumulative frequencies in [0,1].
#' The number of unique values can be determined by finding the last non-zero cdf_value.
#'
#' Generated from the Fortran procedure \code{f42_stats::compute_edf_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{values_perm} yourself.
#' \code{compute_edf} seeds \code{values_perm} and sorts it by \code{values}.
#'
#' @param values a numeric vector. Array of observed data values (e.g., contributions or spikes).
#' @param values_perm a integer vector. Permutation of `values` in ascending order. The allocating entry point builds
#'   and heapsorts it for you; the expert one takes whatever order you supply.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_values`.
#' @return a named list with elements:
#'   \item{unique_values}{a numeric vector. Sorted unique data values.
#'     The first `n_unique` elements will hold the results.}
#'   \item{cdf_values}{a numeric vector. Corresponding cumulative frequencies between 0 and 1.
#'     The first `n_unique` elements will hold the results.}
#' @export
compute_edf_expert <- function(values, values_perm) {
    values <- .tox_as_double_vector(values, "values")
    values_perm <- .tox_as_integer_vector(values_perm, "values_perm")
    if (length(values_perm) != length(values))
        .tox_shape_error("values_perm", length(values_perm), "values", length(values))

    .result <- .Call("compute_edf_expert_call", values, values_perm)
    .arguments <- c("values", "n_values", "values_perm", "unique_values", "cdf_values", "n_unique", "ierr")
    .sources <- c(NA_character_, "values", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        unique_values = utils::head(.result$unique_values, .result$n_unique),
        cdf_values = utils::head(.result$cdf_values, .result$n_unique)
    )
}

#' Calculate the percentile of an array given a sorted permutation
#'
#' Uses linear interpolation between adjacent values.
#'
#' Generated from the Fortran procedure \code{f42_stats::calc_percentile}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{array_perm} and sorts it by \code{array}.
#' Call \code{calc_percentile_expert} to do that yourself.
#'
#' @param array a numeric vector. input array
#' @param percentile a numeric scalar. desired percentile (0-100)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#' @param n_considered a integer scalar. How many leading entries of `array_perm` the percentile is taken over, for a
#'   percentile of a subset -- the trailing entries are ignored rather than sliced
#'   off, so the permutation stays the shape the sort produced. Zero, the default,
#'   considers all `n_array` of them.
#'   The default value is `0`.
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_array`.
#' @return a numeric scalar. output percentile value
#' @export
calc_percentile <- function(array, percentile, n_considered = 0L) {
    array <- .tox_as_double_vector(array, "array")
    percentile <- .tox_as_double_scalar(percentile, "percentile")
    n_considered <- .tox_as_integer_scalar(n_considered, "n_considered")
    .result <- .Call("calc_percentile_call", array, percentile, n_considered)
    .arguments <- c("array", "n_array", "percentile", "value", "n_considered", "ierr")
    .sources <- c(NA_character_, "array", NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$value
}

#' Calculate the percentile of an array given a sorted permutation
#'
#' Uses linear interpolation between adjacent values.
#'
#' Generated from the Fortran procedure \code{f42_stats::calc_percentile_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{array_perm} yourself.
#' \code{calc_percentile} seeds \code{array_perm} and sorts it by \code{array}.
#'
#' @param array a numeric vector. input array
#' @param array_perm a integer vector. Permutation of `array` in ascending order. The allocating entry point builds and
#'   heapsorts it for you; the expert one takes whatever order you supply.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_array`.
#' @param percentile a numeric scalar. desired percentile (0-100)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#' @param n_considered a integer scalar. How many leading entries of `array_perm` the percentile is taken over, for a
#'   percentile of a subset -- the trailing entries are ignored rather than sliced
#'   off, so the permutation stays the shape the sort produced. Zero, the default,
#'   considers all `n_array` of them.
#'   The default value is `0`.
#'   The minimum valid value is `0`.
#'   The maximum valid value is `n_array`.
#' @return a numeric scalar. output percentile value
#' @export
calc_percentile_expert <- function(array, array_perm, percentile, n_considered = 0L) {
    array <- .tox_as_double_vector(array, "array")
    array_perm <- .tox_as_integer_vector(array_perm, "array_perm")
    percentile <- .tox_as_double_scalar(percentile, "percentile")
    n_considered <- .tox_as_integer_scalar(n_considered, "n_considered")
    if (length(array_perm) != length(array))
        .tox_shape_error("array_perm", length(array_perm), "array", length(array))

    .result <- .Call("calc_percentile_expert_call", array, array_perm, percentile, n_considered)
    .arguments <- c("array", "n_array", "array_perm", "percentile", "value", "n_considered", "ierr")
    .sources <- c(NA_character_, "array", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$value
}
