# Generated. Do not edit.

#' Performs LOESS smoothing on a set of data points
#'
#' Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
#' The user must pre-filter data and provide only valid indices in indices_used.
#'
#' Generated from the Fortran procedure \code{f42_stats::loess_smooth_2d}, whose argument names
#' are the ones an error message reports.
#'
#' @param x_ref a numeric vector. Reference x-coordinates.
#' @param y_ref a numeric vector. Reference y-coordinates (length n_total).
#' @param indices_used a integer vector. Indices of reference points used for smoothing (only valid indices).
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_total`.
#' @param x_query a numeric vector. Target x-coordinates to smooth.
#' @param kernel_sigma a numeric scalar. Bandwidth parameter for the kernel.
#'   The minimum valid value is `0.0`.
#' @param kernel_cutoff a numeric scalar. Cutoff for the kernel, not used if zero
#'   The minimum valid value is `0.0`.
#' @return a numeric vector. Output smoothed values (length n_target).
#' @export
loess_smooth_2d <- function(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff) {
    x_ref <- .tox_as_double_vector(x_ref, "x_ref")
    y_ref <- .tox_as_double_vector(y_ref, "y_ref")
    indices_used <- .tox_as_integer_vector(indices_used, "indices_used")
    x_query <- .tox_as_double_vector(x_query, "x_query")
    kernel_sigma <- .tox_as_double_scalar(kernel_sigma, "kernel_sigma")
    kernel_cutoff <- .tox_as_double_scalar(kernel_cutoff, "kernel_cutoff")
    if (length(y_ref) != length(x_ref))
        .tox_shape_error("y_ref", length(y_ref), "x_ref", length(x_ref))

    .result <- .Call("loess_smooth_2d_call", x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    .arguments <- c("n_total", "n_target", "x_ref", "y_ref", "indices_used", "n_used", "x_query", "kernel_sigma", "kernel_cutoff", "y_out", "ierr")
    .sources <- c("x_ref", "x_query", NA_character_, NA_character_, NA_character_, "indices_used", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$y_out
}

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
#' @param percentile a numeric scalar. desired percentile as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
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
#' @param percentile a numeric scalar. desired percentile as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
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

#' Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)
#'
#' This is NOT a null-hypothesis-testing p-value: each distance is compared against the
#' observed distribution it was drawn from, not an independently generated null distribution.
#' It instead measures how extreme an observed distance is relative to all observed distances.
#'
#' Implements:
#' Q(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
#'
#' Because distances are non-negative, a one-sided upper-tail quantile is used.
#'
#' Assumptions / preconditions:
#' - sorted_rdi(1:n_genes) contains the empirical distribution D.
#' - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
#'
#' Generated from the Fortran procedure \code{f42_stats::compute_scaled_distance_quantile}, whose argument names
#' are the ones an error message reports.
#'
#' This entry point seeds \code{sorted_rdi_perm} and sorts it by \code{sorted_rdi}.
#' Call \code{compute_scaled_distance_quantile_expert} to do that yourself.
#'
#' @param rdi a numeric vector. empirical distribution D
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param sorted_rdi a numeric vector. empirical distribution D with non negative values
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param c_const a numeric scalar. Constant used in the computation, typically 1
#' @return a numeric vector. Output array to store the computed quantile for each gene.
#' @export
compute_scaled_distance_quantile <- function(rdi, sorted_rdi, c_const) {
    rdi <- .tox_as_double_vector(rdi, "rdi")
    sorted_rdi <- .tox_as_double_vector(sorted_rdi, "sorted_rdi")
    c_const <- .tox_as_double_scalar(c_const, "c_const")
    if (length(sorted_rdi) != length(rdi))
        .tox_shape_error("sorted_rdi", length(sorted_rdi), "rdi", length(rdi))

    .result <- .Call("compute_scaled_distance_quantile_call", rdi, sorted_rdi, c_const)
    .arguments <- c("n_genes", "rdi", "sorted_rdi", "quantile", "c_const", "ierr")
    .sources <- c("rdi", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$quantile
}

#' Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)
#'
#' This is NOT a null-hypothesis-testing p-value: each distance is compared against the
#' observed distribution it was drawn from, not an independently generated null distribution.
#' It instead measures how extreme an observed distance is relative to all observed distances.
#'
#' Implements:
#' Q(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
#'
#' Because distances are non-negative, a one-sided upper-tail quantile is used.
#'
#' Assumptions / preconditions:
#' - sorted_rdi(1:n_genes) contains the empirical distribution D.
#' - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
#'
#' Generated from the Fortran procedure \code{f42_stats::compute_scaled_distance_quantile_expert}, whose argument names
#' are the ones an error message reports.
#'
#' The expert entry point: you supply \code{sorted_rdi_perm} yourself.
#' \code{compute_scaled_distance_quantile} seeds \code{sorted_rdi_perm} and sorts it by \code{sorted_rdi}.
#'
#' @param rdi a numeric vector. empirical distribution D
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param sorted_rdi a numeric vector. empirical distribution D with non negative values
#'   NaN is permitted for this value.
#'   Infinite values are permitted for this value.
#' @param sorted_rdi_perm a integer vector. Permutation of `sorted_rdi` in ascending order. The allocating entry point builds
#'   and heapsorts it for you; the expert one takes whatever order you supply.
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_genes`.
#' @param c_const a numeric scalar. Constant used in the computation, typically 1
#' @return a numeric vector. Output array to store the computed quantile for each gene.
#' @export
compute_scaled_distance_quantile_expert <- function(rdi, sorted_rdi, sorted_rdi_perm, c_const) {
    rdi <- .tox_as_double_vector(rdi, "rdi")
    sorted_rdi <- .tox_as_double_vector(sorted_rdi, "sorted_rdi")
    sorted_rdi_perm <- .tox_as_integer_vector(sorted_rdi_perm, "sorted_rdi_perm")
    c_const <- .tox_as_double_scalar(c_const, "c_const")
    if (length(sorted_rdi) != length(rdi))
        .tox_shape_error("sorted_rdi", length(sorted_rdi), "rdi", length(rdi))
    if (length(sorted_rdi_perm) != length(rdi))
        .tox_shape_error("sorted_rdi_perm", length(sorted_rdi_perm), "rdi", length(rdi))

    .result <- .Call("compute_scaled_distance_quantile_expert_call", rdi, sorted_rdi, sorted_rdi_perm, c_const)
    .arguments <- c("n_genes", "rdi", "sorted_rdi", "sorted_rdi_perm", "quantile", "c_const", "ierr")
    .sources <- c("rdi", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$quantile
}
