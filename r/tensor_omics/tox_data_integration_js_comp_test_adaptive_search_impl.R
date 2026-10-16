# Generated. Do not edit.

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
#' `floor(huge(1)/ADAPTIVE_K_MAX_FACTOR)` (\code{ADAPTIVE_K_MAX_FACTOR}), so that its
#' `k_max = ceiling(ADAPTIVE_K_MAX_FACTOR*k_start_1)` fits; the output is then left undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_adaptive_search_impl::calc_adaptive_js_comp_test_bounds}, whose argument names
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
