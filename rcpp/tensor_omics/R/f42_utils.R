# Generated. Do not edit.

#' Performs LOESS smoothing on a set of data points
#'
#' Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
#' The user must pre-filter data and provide only valid indices in indices_used.
#'
#' @param x_ref a numeric vector. Reference x-coordinates.
#' @param y_ref a numeric vector. Reference y-coordinates (length n_total).
#' @param indices_used a integer vector. Indices of reference points used for smoothing (only valid indices).
#' @param x_query a numeric vector. Target x-coordinates to smooth.
#' @param kernel_sigma a numeric scalar. Bandwidth parameter for the kernel.
#' @param kernel_cutoff a numeric scalar. Cutoff for the kernel, not used if zero
#' @return Output smoothed values (length n_target).
#'
#' Generated from the Fortran procedure \code{f42_utils::loess_smooth_2d}.
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

    .result <- .loess_smooth_2d_rcpp(x_ref, y_ref, indices_used, x_query, kernel_sigma, kernel_cutoff)
    .arguments <- c("n_total", "n_target", "x_ref", "y_ref", "indices_used", "n_used", "x_query", "kernel_sigma", "kernel_cutoff", "y_out", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$y_out
}

#' Compute the Empirical Distribution Function (EDF) from pre-sorted permutation
#'
#' Returns the sorted unique values and their cumulative frequencies in [0,1].
#' Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
#' The number of unique values can be determined by finding the last non-zero cdf_value.
#'
#' @param values a numeric vector. Array of observed data values (e.g., contributions or spikes).
#' @param perm a integer vector. Pre-sorted permutation indices (must be sorted by values[perm]).
#' @return a named list with elements `unique_values`, `cdf_values`.
#'
#' Generated from the Fortran procedure \code{f42_utils::compute_edf}.
#' @export
compute_edf_expert <- function(values, perm) {
    values <- .tox_as_double_vector(values, "values")
    perm <- .tox_as_integer_vector(perm, "perm")
    if (length(perm) != length(values))
        .tox_shape_error("perm", length(perm), "values", length(values))

    .result <- .compute_edf_expert_rcpp(values, perm)
    .arguments <- c("values", "n_values", "perm", "unique_values", "cdf_values", "n_unique", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

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
#' @param values a numeric vector. Array of observed data values (e.g., contributions or spikes).
#' @return a named list with elements `unique_values`, `cdf_values`.
#'
#' Generated from the Fortran procedure \code{f42_utils::compute_edf_alloc}.
#' @export
compute_edf <- function(values) {
    values <- .tox_as_double_vector(values, "values")
    .result <- .compute_edf_rcpp(values)
    .arguments <- c("values", "n_values", "unique_values", "cdf_values", "n_unique", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        unique_values = utils::head(.result$unique_values, .result$n_unique),
        cdf_values = utils::head(.result$cdf_values, .result$n_unique)
    )
}

#' Calculate empirical p-values for scaled expression distances (RDI)
#'
#' Implements:
#' P(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
#'
#' Because distances are non-negative, a one-sided upper-tail empirical p-value is used.
#'
#' Assumptions / preconditions:
#' - sorted_rdi(1:n_genes) contains the empirical distribution D.
#' - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
#'
#' @param rdi a numeric vector. empirical distribution D
#' @param sorted_rdi a numeric vector. empirical distribution D with non negative values
#' @param perm a integer vector. Permutation array with sorted indices for sorted_rdi
#' @param c_const a numeric scalar. Constant used in the computation, typically 1
#' @return Output array to store the computed p-values for each gene.
#'
#' Generated from the Fortran procedure \code{f42_utils::compute_empirical_p_values}.
#' @export
compute_empirical_p_values <- function(rdi, sorted_rdi, perm, c_const) {
    rdi <- .tox_as_double_vector(rdi, "rdi")
    sorted_rdi <- .tox_as_double_vector(sorted_rdi, "sorted_rdi")
    perm <- .tox_as_integer_vector(perm, "perm")
    c_const <- .tox_as_double_scalar(c_const, "c_const")
    if (length(sorted_rdi) != length(rdi))
        .tox_shape_error("sorted_rdi", length(sorted_rdi), "rdi", length(rdi))
    if (length(perm) != length(rdi))
        .tox_shape_error("perm", length(perm), "rdi", length(rdi))

    .result <- .compute_empirical_p_values_rcpp(rdi, sorted_rdi, perm, c_const)
    .arguments <- c("n_genes", "rdi", "sorted_rdi", "perm", "p_values", "c_const")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$p_values
}
