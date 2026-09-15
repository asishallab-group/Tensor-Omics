# Generated. Do not edit.

#' Recommend upper bounds for the js-comp-test candidate-grid work arrays
#'
#' Closed-form from the GAMMA-decay candidate-grid formula
#' \code{\link{generate_js_comp_test_candidates}}
#' uses: `max_n_points_candidate` is exactly the grid's first (largest) `n_points` candidate.
#' `max_n_neighbors_candidate` is a safe, not necessarily tight, upper bound on the grid's
#' largest `n_neighbors` candidate -- reached with the smallest `KX_FACTORS` entry at the
#' smallest `n_points_high` the grid loop ever uses, which by construction never drops below
#' `n_points_low`.
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
    .arguments <- c("max_n_genes_all_studies", "max_n_points_candidate", "max_n_neighbors_candidate")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_n_points_candidate = .result$max_n_points_candidate,
        max_n_neighbors_candidate = .result$max_n_neighbors_candidate
    )
}

#' Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
#'
#' Ported from 125-stabilize-jscomp's inline `n_bootstrapping_top_k_jsds` computation in
#' `determine_js_comp_test_n_points_n_neighbors_alloc`: sizes the top-k/bottom-k heaps
#' \code{\link{bootstrap_histogram}} uses
#' to track a two-sided bootstrap confidence interval's endpoints. E.g. for a 95% two-sided
#' interval (2.5% reserved at each end, the default) with `n_bootstraps=1000`,
#' `n_top_k = max(1, floor(0.025*1000)) = 25`.
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
