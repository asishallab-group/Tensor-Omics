# Generated. Do not edit.

#' Calculate the number of neighbors to be used for constructing neighborhoods
#'
#' The `desired_size` works as upper limit, as the actual neighborhood size might be lower
#' due to few genes with non-NaN mean.
#'
#' @param n_pool a integer scalar. Total number of pooled mean-expression values across both studies
#' @param n_points a integer scalar. Number of reference points
#' @param mean_S a numeric vector. Per-gene mean expression values
#'   NaN is permitted for this value.
#' @param desired_size a integer scalar. Optional desired neighborhood size
#'   The default value is `1000`.
#' @return Calculated neighborhood size
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing_kernel::calc_neighborhood_size}.
#' @export
calc_neighborhood_size <- function(n_pool, n_points, mean_S, desired_size = 1000L) {
    n_pool <- .tox_as_integer_scalar(n_pool, "n_pool")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    desired_size <- .tox_as_integer_scalar(desired_size, "desired_size")
    .result <- .Call("calc_neighborhood_size_call", n_pool, n_points, mean_S, desired_size)
    .arguments <- c("n_pool", "n_points", "n_genes_S", "mean_S", "desired_size")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$n_neighbors
}
