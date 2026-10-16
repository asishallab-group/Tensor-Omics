# Generated. Do not edit.

#' Calculate the number of neighbors to be used for constructing neighborhoods
#'
#' The `desired_size` works as upper limit, as the actual neighborhood size might be lower
#' due to few genes with non-NaN mean.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing_impl::calc_neighborhood_size}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_pool a integer scalar. Total number of pooled mean-expression values across both studies
#' @param n_points a integer scalar. Number of reference points
#' @param mean_S a numeric vector. Per-gene mean expression values
#'   NaN is permitted for this value.
#' @param desired_size a integer scalar. Optional desired neighborhood size
#'   The default value is `1000`.
#' @return a integer scalar. Calculated neighborhood size
#' @export
calc_neighborhood_size <- function(n_pool, n_points, mean_S, desired_size = 1000L) {
    n_pool <- .tox_as_integer_scalar(n_pool, "n_pool")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    mean_S <- .tox_as_double_vector(mean_S, "mean_S")
    desired_size <- .tox_as_integer_scalar(desired_size, "desired_size")
    .result <- .Call("calc_neighborhood_size_call", n_pool, n_points, mean_S, desired_size)
    .arguments <- c("n_pool", "n_points", "n_genes_S", "mean_S", "desired_size")
    .sources <- c(NA_character_, NA_character_, "mean_S", NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$n_neighbors
}

#' Pool one reference point's residuals across every neighbor and every study
#'
#' Gathers one reference point's residuals from every neighbor gene of every study into one flat
#' array, given that point's per-study neighbor gene indices (one column of the
#' `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
#' \code{\link{construct_neighborhoods_ranged}}).
#' Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
#' `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
#' \code{\link{run_js_comp_test}} and
#' \code{\link{run_js_comp_test_parameter_search}}
#' feed to the bin-count search, so a caller can pass the same input to
#' \code{\link{determine_bin_count_occupancy_exhaustive}}
#' on real data.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_preprocessing_impl::gather_pooled_neighborhood_residuals}, whose argument names
#' are the ones an error message reports.
#'
#' @param residuals a numeric array of rank 3. Matrix of signed residuals per study, NaN allowed for missing values; not validated
#' @param neighborhood_indices_point a integer matrix. Gene indices of one reference point's neighborhood, per study; each in
#'   [1, max_n_genes_all_studies]
#' @return a numeric vector. This reference point's pooled residuals, laid out as a
#'   (max_n_reps_all_studies, n_neighbors, n_studies) array
#' @export
gather_pooled_neighborhood_residuals <- function(residuals, neighborhood_indices_point) {
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    neighborhood_indices_point <- .tox_as_integer_matrix(neighborhood_indices_point, "neighborhood_indices_point")
    if (dim(neighborhood_indices_point)[2] != dim(residuals)[3])
        .tox_shape_error("neighborhood_indices_point", dim(neighborhood_indices_point)[2], "residuals", dim(residuals)[3])

    .result <- .Call("gather_pooled_neighborhood_residuals_call", residuals, neighborhood_indices_point)
    .arguments <- c("residuals", "max_n_reps_all_studies", "max_n_genes_all_studies", "n_neighbors", "n_studies", "neighborhood_indices_point", "pooled_residuals", "ierr")
    .sources <- c(NA_character_, "residuals", "residuals", "neighborhood_indices_point", "residuals", NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    .result$pooled_residuals
}
