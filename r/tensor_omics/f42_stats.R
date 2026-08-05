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
#' @param x_query a numeric vector. Target x-coordinates to smooth.
#' @param kernel_sigma a numeric scalar. Bandwidth parameter for the kernel.
#' @param kernel_cutoff a numeric scalar. Cutoff for the kernel, not used if zero
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

#' Compute the Empirical Distribution Function (EDF) from pre-sorted permutation
#'
#' Returns the sorted unique values and their cumulative frequencies in [0,1].
#' Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
#' The number of unique values can be determined by finding the last non-zero cdf_value.
#'
#' Generated from the Fortran procedure \code{f42_stats::compute_edf}, whose argument names
#' are the ones an error message reports.
#'
#' @param values a numeric vector. Array of observed data values (e.g., contributions or spikes).
#' @param perm a integer vector. Pre-sorted permutation indices (must be sorted by values[perm]).
#' @return a named list with elements:
#'   \item{unique_values}{a numeric vector. Sorted unique data values.
#'     The first `n_unique` elements will hold the results.}
#'   \item{cdf_values}{a numeric vector. Corresponding cumulative frequencies between 0 and 1.
#'     The first `n_unique` elements will hold the results.}
#' @export
compute_edf_expert <- function(values, perm) {
    values <- .tox_as_double_vector(values, "values")
    perm <- .tox_as_integer_vector(perm, "perm")
    if (length(perm) != length(values))
        .tox_shape_error("perm", length(perm), "values", length(values))

    .result <- .Call("compute_edf_expert_call", values, perm)
    .arguments <- c("values", "n_values", "perm", "unique_values", "cdf_values", "n_unique", "ierr")
    .sources <- c(NA_character_, "values", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        unique_values = utils::head(.result$unique_values, .result$n_unique),
        cdf_values = utils::head(.result$cdf_values, .result$n_unique)
    )
}

#' Sorts the values and computes the Empirical Distribution Function (EDF)
#'
#' Allocates workspace internally and performs sorting before computing EDF.
#' Use this for convenience; use compute_edf directly for custom sorting.
#'
#' Generated from the Fortran procedure \code{f42_stats::compute_edf_alloc}, whose argument names
#' are the ones an error message reports.
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
#' @param rdi a numeric vector. empirical distribution D
#' @param sorted_rdi a numeric vector. empirical distribution D with non negative values
#' @param perm a integer vector. Permutation array with sorted indices for sorted_rdi
#' @param c_const a numeric scalar. Constant used in the computation, typically 1
#' @return a numeric vector. Output array to store the computed quantile for each gene.
#' @export
compute_scaled_distance_quantile <- function(rdi, sorted_rdi, perm, c_const) {
    rdi <- .tox_as_double_vector(rdi, "rdi")
    sorted_rdi <- .tox_as_double_vector(sorted_rdi, "sorted_rdi")
    perm <- .tox_as_integer_vector(perm, "perm")
    c_const <- .tox_as_double_scalar(c_const, "c_const")
    if (length(sorted_rdi) != length(rdi))
        .tox_shape_error("sorted_rdi", length(sorted_rdi), "rdi", length(rdi))
    if (length(perm) != length(rdi))
        .tox_shape_error("perm", length(perm), "rdi", length(rdi))

    .result <- .Call("compute_scaled_distance_quantile_call", rdi, sorted_rdi, perm, c_const)
    .arguments <- c("n_genes", "rdi", "sorted_rdi", "perm", "quantile", "c_const")
    .sources <- c("rdi", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$quantile
}
