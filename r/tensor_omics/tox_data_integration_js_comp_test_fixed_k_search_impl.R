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
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_fixed_k_search_impl::calc_js_comp_test_candidate_bounds}, whose argument names
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
