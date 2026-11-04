# Generated. Do not edit.

#' Recommend the work-array capacities of the k-NN study-label Permutation-Test
#'
#' Sizes the per-point and per-neighbor work arrays of
#' \code{\link{run_js_comp_test_on_sample_assignment}}
#' and
#' \code{\link{run_js_comp_test_permutation}}
#' so that both the fixed `(n_points, n_neighbors)` run and a rerun of the parameter search fit:
#' `max_n_points_candidate = max(n_points, P)` and `max_n_neighbors_candidate = max(n_neighbors,
#' K)`, with `P` and `K` the bounds of
#' \code{\link{calc_js_comp_test_candidate_bounds}}
#' for the same `grid_points_factor` and `grid_min_points`.
#' `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
#' \code{\link{calc_js_comp_test_n_top_k_jsds}} computes it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation_bounds::calc_js_comp_test_permutation_bounds}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_points a integer scalar. Number of reference points of the fixed-parameter run
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Number of neighbors per neighborhood of the fixed-parameter run
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param n_bootstraps a integer scalar. Number of bootstraps of the parameter search, when it is rerun
#'   The minimum valid value is `1`.
#'   The default value is `1000`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param grid_points_factor a numeric scalar. The rerun parameter search's `grid_points_factor`
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `4.0`.
#' @param grid_min_points a integer scalar. The rerun parameter search's `grid_min_points`
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_POINTS`.
#'   The default value is `200`.
#' @return a named list with elements:
#'   \item{max_n_points_candidate}{a integer scalar. Capacity for reference points}
#'   \item{max_n_neighbors_candidate}{a integer scalar. Capacity for neighbors per neighborhood}
#'   \item{n_bootstrapping_top_k_jsds}{a integer scalar. Bootstrap heap size (values kept at each end of the bootstrap distribution)}
#' @export
calc_js_comp_test_permutation_bounds <- function(max_n_genes_all_studies, n_points, n_neighbors, n_bootstraps = 1000L, two_sided_bootstrapping_significance_level = 2.5, grid_points_factor = 4.0, grid_min_points = 200L) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    grid_points_factor <- .tox_as_double_scalar(grid_points_factor, "grid_points_factor")
    grid_min_points <- .tox_as_integer_scalar(grid_min_points, "grid_min_points")
    .result <- .Call("calc_js_comp_test_permutation_bounds_call", max_n_genes_all_studies, n_points, n_neighbors, n_bootstraps, two_sided_bootstrapping_significance_level, grid_points_factor, grid_min_points)
    .arguments <- c("max_n_genes_all_studies", "n_points", "n_neighbors", "n_bootstraps", "two_sided_bootstrapping_significance_level", "grid_points_factor", "grid_min_points", "max_n_points_candidate", "max_n_neighbors_candidate", "n_bootstrapping_top_k_jsds", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_n_points_candidate = .result$max_n_points_candidate,
        max_n_neighbors_candidate = .result$max_n_neighbors_candidate,
        n_bootstrapping_top_k_jsds = .result$n_bootstrapping_top_k_jsds
    )
}

#' Recommend the work-array capacity of the adaptive study-label Permutation-Test
#'
#' Sizes the per-point work arrays of
#' \code{\link{run_js_comp_test_adaptive_on_sample_assignment}}
#' and
#' \code{\link{run_js_comp_test_adaptive_permutation}}:
#' `max_n_points_candidate = max(A, min(N, ceiling(F * N / k_start) + 1))`, with `N =
#' max_n_genes_all_studies * n_studies` (the padded pool size), `F` the factor
#' \code{ADAPTIVE_POINT_CAPACITY_FACTOR}
#' and `A` the capacity of the adaptive parameter search,
#' \code{\link{calc_adaptive_js_comp_test_bounds}}.
#' As there, this is a **practical capacity, not a proven bound**: how many reference points a
#' construction builds is only known once it has run, and a shuffled data set whose
#' construction builds more fails with invalid input (it is then counted as a failed
#' permutation). `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
#' \code{\link{calc_js_comp_test_n_top_k_jsds}} computes it.
#'
#' Reports invalid input, naming neither argument, for the joint bounds on
#' `max_n_genes_all_studies` and `n_studies` that the adaptive search's capacity rejects (a pool
#' size or a first candidate's `k_max` beyond the 32-bit range); the outputs are then left
#' undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation_bounds::calc_adaptive_js_comp_test_permutation_bounds}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @param k_start a integer scalar. Pooled entries every neighborhood of the fixed-knob run takes unconditionally
#'   The minimum valid value is `1`.
#' @param n_bootstraps a integer scalar. Number of bootstraps of the parameter search, when it is rerun
#'   The minimum valid value is `1`.
#'   The default value is `1000`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @return a named list with elements:
#'   \item{max_n_points_candidate}{a integer scalar. Practical capacity for reference points}
#'   \item{n_bootstrapping_top_k_jsds}{a integer scalar. Bootstrap heap size (values kept at each end of the bootstrap distribution)}
#' @export
calc_adaptive_js_comp_test_permutation_bounds <- function(max_n_genes_all_studies, n_studies, k_start, n_bootstraps = 1000L, two_sided_bootstrapping_significance_level = 2.5) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    k_start <- .tox_as_integer_scalar(k_start, "k_start")
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    .result <- .Call("calc_adaptive_js_comp_test_permutation_bounds_call", max_n_genes_all_studies, n_studies, k_start, n_bootstraps, two_sided_bootstrapping_significance_level)
    .arguments <- c("max_n_genes_all_studies", "n_studies", "k_start", "n_bootstraps", "two_sided_bootstrapping_significance_level", "max_n_points_candidate", "n_bootstrapping_top_k_jsds", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    list(
        max_n_points_candidate = .result$max_n_points_candidate,
        n_bootstrapping_top_k_jsds = .result$n_bootstrapping_top_k_jsds
    )
}
