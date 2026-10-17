# Generated. Do not edit.

#' Recommend upper bounds for the js-comp-test candidate-grid work arrays
#'
#' For the candidate grid of
#' \code{\link{generate_js_comp_test_candidates}}:
#' `max_n_points_candidate` is exactly its largest `n_points` candidate, and
#' `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`,
#' in every `neighborhood_size_mode` (the bound uses the large-neighborhood factor `KX_FACTORS(1)`)
#' and every `grid_floor_mode` (the floor candidate of the clamp sits at the floor this bound
#' already uses). Pass the same `grid_sample_size_mode`, `grid_points_factor` and
#' `grid_min_points` as the parameter search that is sized with the result.
#'
#' Rejects a non-positive `max_n_genes_all_studies` or `n_studies`, an unknown
#' `grid_sample_size_mode`, a non-positive `grid_points_factor` and a `grid_min_points` outside
#' `[1, MAX_POINTS]` with invalid input (a non-finite factor with `ERR_NAN_INF`); both outputs
#' are then left undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_fixed_k_search_impl::calc_js_comp_test_candidate_bounds}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_studies a integer scalar. Number of studies; used only by the pooled grid sample size
#'   The minimum valid value is `1`.
#' @param grid_sample_size_mode a string, one of "grid_sample_size_per_study", "grid_sample_size_pooled". Sample size `N` in the grid's first reference-point count; must match the parameter
#'   search's own. An experiment switch
#'
#'   The default value is `"grid_sample_size_per_study"`.
#' @param grid_points_factor a numeric scalar. Factor in the grid's first reference-point count `ceiling(grid_points_factor * sqrt(N))`;
#'   must match the parameter search's own. An experiment switch
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `4.0`.
#' @param grid_min_points a integer scalar. The grid's floor; must match the parameter search's own. An experiment switch: 300
#'   (`MIN_POINTS`) is the floor before it existed
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_POINTS`.
#'   The default value is `200`.
#' @return a named list with elements:
#'   \item{max_n_points_candidate}{a integer scalar. Exact upper bound on the `n_points` candidate the grid ever produces}
#'   \item{max_n_neighbors_candidate}{a integer scalar. Safe upper bound on the `n_neighbors` candidate the grid ever produces}
#' @export
calc_js_comp_test_candidate_bounds <- function(max_n_genes_all_studies, n_studies, grid_sample_size_mode = "grid_sample_size_per_study", grid_points_factor = 4.0, grid_min_points = 200L) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    grid_sample_size_mode <- .tox_as_mode(grid_sample_size_mode, "grid_sample_size_mode", c("grid_sample_size_per_study", "grid_sample_size_pooled"))
    grid_points_factor <- .tox_as_double_scalar(grid_points_factor, "grid_points_factor")
    grid_min_points <- .tox_as_integer_scalar(grid_min_points, "grid_min_points")
    .result <- .Call("calc_js_comp_test_candidate_bounds_call", max_n_genes_all_studies, n_studies, grid_sample_size_mode, grid_points_factor, grid_min_points)
    .arguments <- c("max_n_genes_all_studies", "n_studies", "max_n_points_candidate", "max_n_neighbors_candidate", "grid_sample_size_mode", "grid_points_factor", "grid_min_points", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_n_points_candidate = .result$max_n_points_candidate,
        max_n_neighbors_candidate = .result$max_n_neighbors_candidate
    )
}
