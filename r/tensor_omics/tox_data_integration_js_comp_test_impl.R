# Generated. Do not edit.

#' Recommend upper bounds for the js-comp-test candidate-grid work arrays
#'
#' For the candidate grid of
#' \code{\link{generate_js_comp_test_candidates}}:
#' `max_n_points_candidate` is exactly its largest `n_points` candidate, and
#' `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`.
#'
#' Rejects a non-positive `max_n_genes_all_studies` with invalid input; both outputs are then
#' left undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_impl::calc_js_comp_test_candidate_bounds}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @return a named list with elements:
#'   \item{max_n_points_candidate}{a integer scalar. Exact upper bound on the `n_points` candidate the grid ever produces}
#'   \item{max_n_neighbors_candidate}{a integer scalar. Safe upper bound on the `n_neighbors` candidate the grid ever produces}
#' @export
calc_js_comp_test_candidate_bounds <- function(max_n_genes_all_studies) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    .result <- .Call("calc_js_comp_test_candidate_bounds_call", max_n_genes_all_studies)
    .arguments <- c("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_n_points_candidate = .result$max_n_points_candidate,
        max_n_neighbors_candidate = .result$max_n_neighbors_candidate
    )
}

#' Recommend a reference-point capacity for the adaptive js-comp-test candidate sequence
#'
#' Sizes the per-point work arrays of the adaptive parameter search over the candidates of
#' \code{\link{generate_adaptive_js_comp_test_candidates}}.
#' A candidate's reference-point count is only known after its neighborhoods are grown, so this
#' is a **practical capacity, not a proven bound**:
#' `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
#' \code{ADAPTIVE_POINT_CAPACITY_FACTOR},
#' `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
#' candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.
#'
#' Rejects, with invalid input naming neither argument when the bound is joint, a
#' `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
#' first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
#' `huge(1)/4 = 536870911` (so that its `k_max = 4*k_start_1` fits); the output is then
#' left undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_impl::calc_adaptive_js_comp_test_bounds}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @return a integer scalar. Practical upper bound on the number of reference points any adaptive candidate is
#'   allowed to emerge with
#' @export
calc_adaptive_js_comp_test_bounds <- function(max_n_genes_all_studies, n_studies) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    .result <- .Call("calc_adaptive_js_comp_test_bounds_call", max_n_genes_all_studies, n_studies)
    .arguments <- c("max_n_genes_all_studies", "n_studies", "max_n_points_candidate", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$max_n_points_candidate
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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_impl::gather_pooled_neighborhood_residuals}, whose argument names
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

#' Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
#'
#' Sizes the top-k/bottom-k heaps
#' \code{\link{bootstrap_histogram}} uses
#' for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
#' the default 2.5 and `n_bootstraps=1000`.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_impl::calc_js_comp_test_n_top_k_jsds}, whose argument names
#' are the ones an error message reports.
#'
#' @param n_bootstraps a integer scalar. Number of bootstrap resamples that will be performed
#'   The minimum valid value is `1`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level, as a percentage reserved at each end of the bootstrap
#'   distribution (e.g. 2.5 for a 95% two-sided interval)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @return a integer scalar. Number of elements to keep at each end of the bootstrap distribution, at least 1
#' @export
calc_js_comp_test_n_top_k_jsds <- function(n_bootstraps, two_sided_bootstrapping_significance_level = 2.5) {
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    .result <- .Call("calc_js_comp_test_n_top_k_jsds_call", n_bootstraps, two_sided_bootstrapping_significance_level)
    .arguments <- c("n_bootstraps", "two_sided_bootstrapping_significance_level", "n_top_k")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$n_top_k
}
