# Generated. Do not edit.

#' Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
#'
#' Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
#' `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
#' `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
#' `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
#' Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
#' candidate. Bin counts are not part of a candidate; they are chosen per reference point later.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_fixed_k_search::generate_js_comp_test_candidates}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @return a integer matrix. Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
#'   The first `n_candidates` elements will hold the results.
#' @export
generate_js_comp_test_candidates <- function(max_n_genes_all_studies) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    .result <- .Call("generate_js_comp_test_candidates_call", max_n_genes_all_studies)
    .arguments <- c("max_n_genes_all_studies", "candidates_n_points_n_neighbors", "n_candidates", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$candidates_n_points_n_neighbors[, seq_len(.result$n_candidates), drop = FALSE]
}

#' Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting
#'
#' Walks the candidates of
#' \code{\link{generate_js_comp_test_candidates}}
#' from finest to coarsest resolution. For each candidate:
#'
#' 1. Reference points from the pooled gene means
#' (\code{\link{pool_means}}) and every
#' study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
#' (\code{\link{check_neighborhood_overlaps}}).
#' 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
#' residuals, then every study's histograms.
#' 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
#' `min_residuals_per_bin`
#' (\code{\link{check_mean_pmf_min_counts}}).
#' A candidate that passes both gates is *admissible*.
#' 4. The observed global JSD per study, its bootstrap confidence interval
#' (\code{\link{bootstrap_histogram}}),
#' and both plateau criteria: CI overlap with the best candidate so far
#' (\code{\link{check_plateau_condition}})
#' and relative effect-size stability
#' (\code{\link{check_effect_size_plateau_condition}}).
#' Both are always computed and traced; `plateau_mode` decides which one stops the search. A
#' stop by the effect-size criterion returns the candidate that triggered it.
#'
#' The returned candidate:
#'
#' 1. A plateau was found, or the grid has a single candidate and it is admissible: the best
#' candidate, `plateau_established = TRUE`.
#' 2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
#' median confidence-interval width across studies, `plateau_established = FALSE`.
#' 3. No candidate was admissible: the first (finest) candidate, with
#' `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
#' ranges `0.0`, and `plateau_established = FALSE`.
#'
#' Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
#' for its meaningful part. A failure to create the bootstrap's random number generator is
#' returned in `ierr` and stops the search; the outputs are then undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_fixed_k_search::run_js_comp_test_parameter_search}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Per-gene mean expression values for all studies
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Matrix of signed residuals per study
#'   NaN is permitted for this value.
#' @param n_bootstraps a integer scalar. Number of bootstraps to perform for a candidate pair
#'   The minimum valid value is `1`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the studies' confidence-interval overlaps are joined into the plateau condition
#' @param max_n_points_candidate a integer scalar. Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
#'   at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
#'   It is recommended to compute this argument from the `max_n_points_candidate` output produced by \code{\link{calc_js_comp_test_candidate_bounds}}.
#'   The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
#' @param max_n_neighbors_candidate a integer scalar. Upper bound on the grid's `n_neighbors` candidates; at least what
#'   `calc_js_comp_test_candidate_bounds` returns, else invalid input
#'   It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by \code{\link{calc_js_comp_test_candidate_bounds}}.
#'   The minimum valid value is `max(1, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2*real(max(MIN_POINTS, min(ceiling(4.0*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
#' @param min_residuals_per_bin a integer scalar. Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
#'   one threshold for both, so a point the search marks as failed also fails the gate
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap of two consecutive neighborhoods (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Minimum fractional overlap of a candidate's confidence interval with the best one so
#'   far, per `join_method`. `0` makes the first admissible candidate plateau under the
#'   CI-overlap criterion (for baseline sweeps)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Which plateau criterion decides when the search stops
#'
#'   The default value is `"plateau_ci_overlap"`.
#' @param delta_median_threshold a numeric scalar. Upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Floor of the denominators of the relative JSD change and relative CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Number of consecutive qualifying transitions for an effect-size plateau
#'   The minimum valid value is `1`.
#'   The default value is `2`.
#' @param m_min a integer scalar. Smallest bin count the bin-count search tests (M_min)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `3`.
#' @param m_max a integer scalar. Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_N_BINS`.
#'   The default value is `120`.
#' @param gamma_occupancy a numeric scalar. Geometric growth factor of the occupancy search; unused in the Stone modes
#'   The minimum valid value is `above(1.0)`.
#'   The default value is `1.25`.
#' @param lower_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's lower residual-range bound
#'   Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
#'   with invalid input reported on upper_residual_range_quantile.
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.05`.
#' @param upper_residual_range_quantile a numeric scalar. Quantile in [0,1] of each reference point's upper residual-range bound
#'   The minimum valid value is `lower_residual_range_quantile`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.95`.
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed for the GSL random number generator
#'   The default value is `42`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen, for every candidate; the
#'   Stone modes test every bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{n_points}{a integer scalar. The finally chosen candidate's `n_points`}
#'   \item{n_neighbors}{a integer scalar. The finally chosen candidate's `n_neighbors`}
#'   \item{n_bins_per_point}{a integer vector. The chosen candidate's bin count per reference point; the first `n_points` entries
#'     are meaningful, the rest are 0}
#'   \item{shared_residual_range_low}{a numeric vector. The chosen candidate's lower residual-range bound (R_low) per reference point; the
#'     first `n_points` entries are meaningful, `0.0` when no candidate was admissible; the
#'     rest are `0.0`}
#'   \item{shared_residual_range_high}{a numeric vector. The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`}
#'   \item{best_candidate_pair_confidence_interval}{a numeric matrix. The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
#'     was admissible}
#'   \item{plateau_established}{a logical scalar. `TRUE` when a plateau was found (by the criterion `plateau_mode` selects) or the
#'     only candidate was admissible; `FALSE` otherwise, also when the smallest-uncertainty
#'     fallback returns a real candidate and interval}
#'   \item{trace_n_points}{a integer vector. Per-admissible-candidate `n_points`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_neighbors}{a integer vector. Per-admissible-candidate `n_neighbors`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_global_js_divergence}{a numeric matrix. Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_lower}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
#'     (`L_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_upper}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
#'     (`U_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width}{a numeric matrix. Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
#'     L_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width_relative}{a numeric matrix. Per-admissible-candidate, per-study relative confidence-interval width
#'     (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_overlap}{a numeric matrix. Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
#'     interval with the best candidate's so far, relative to its own width (the values
#'     `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta}{a numeric matrix. Per-admissible-candidate, per-study relative JSD change from the previous admissible
#'     candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_median}{a numeric vector. Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta_max}{a numeric vector. Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_selected_n_bins}{a integer matrix. Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
#'     `1:trace_n_points(t)` of column `t` are meaningful (the rest are 0), so a Python/R
#'     caller slices `[:trace_n_points[t], t]`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_occupancy_failed}{a logical matrix. Per-admissible-candidate, per-reference-point flag: `TRUE` iff the bin-count search
#'     failed and the point used `m_min` bins. Jagged like trace_selected_n_bins (`FALSE`
#'     beyond)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_pooled_residuals}{a integer matrix. Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
#'     like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_min_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_mean_bin_occupancy}{a numeric matrix. Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_max_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_sturges_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_fd_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
#'     Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_low}{a numeric matrix. Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
#'     Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_high}{a numeric matrix. Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
#'     Jagged like trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#' @export
run_js_comp_test_parameter_search <- function(gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, max_n_neighbors_candidate, min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_ci_overlap", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bin_count_mode = "bin_count_occupancy") {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
    max_n_points_candidate <- .tox_as_integer_scalar(max_n_points_candidate, "max_n_points_candidate")
    max_n_neighbors_candidate <- .tox_as_integer_scalar(max_n_neighbors_candidate, "max_n_neighbors_candidate")
    min_residuals_per_bin <- .tox_as_integer_scalar(min_residuals_per_bin, "min_residuals_per_bin")
    min_neighbor_overlap <- .tox_as_double_scalar(min_neighbor_overlap, "min_neighbor_overlap")
    succeeding_ci_overlap <- .tox_as_double_scalar(succeeding_ci_overlap, "succeeding_ci_overlap")
    plateau_mode <- .tox_as_mode(plateau_mode, "plateau_mode", c("plateau_ci_overlap", "plateau_effect_size", "plateau_both"))
    delta_median_threshold <- .tox_as_double_scalar(delta_median_threshold, "delta_median_threshold")
    delta_max_threshold <- .tox_as_double_scalar(delta_max_threshold, "delta_max_threshold")
    delta_epsilon <- .tox_as_double_scalar(delta_epsilon, "delta_epsilon")
    delta_min_consecutive_transitions <- .tox_as_integer_scalar(delta_min_consecutive_transitions, "delta_min_consecutive_transitions")
    m_min <- .tox_as_integer_scalar(m_min, "m_min")
    m_max <- .tox_as_integer_scalar(m_max, "m_max")
    gamma_occupancy <- .tox_as_double_scalar(gamma_occupancy, "gamma_occupancy")
    lower_residual_range_quantile <- .tox_as_double_scalar(lower_residual_range_quantile, "lower_residual_range_quantile")
    upper_residual_range_quantile <- .tox_as_double_scalar(upper_residual_range_quantile, "upper_residual_range_quantile")
    two_sided_bootstrapping_significance_level <- .tox_as_double_scalar(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed <- .tox_as_integer_scalar(random_seed, "random_seed")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_parameter_search_call", gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, max_n_neighbors_candidate, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "max_n_neighbors_candidate", "n_points", "n_neighbors", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_pair_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_n_points", "trace_n_neighbors", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, "n_bins_per_point", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        n_points = .result$n_points,
        n_neighbors = .result$n_neighbors,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        best_candidate_pair_confidence_interval = .result$best_candidate_pair_confidence_interval,
        plateau_established = .result$plateau_established,
        trace_n_points = utils::head(.result$trace_n_points, .result$n_admissible_evaluated),
        trace_n_neighbors = utils::head(.result$trace_n_neighbors, .result$n_admissible_evaluated),
        trace_global_js_divergence = .result$trace_global_js_divergence[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_lower = .result$trace_ci_lower[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_upper = .result$trace_ci_upper[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_width = .result$trace_ci_width[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_width_relative = .result$trace_ci_width_relative[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_ci_overlap = .result$trace_ci_overlap[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_delta = .result$trace_delta[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_delta_median = utils::head(.result$trace_delta_median, .result$n_admissible_evaluated),
        trace_delta_max = utils::head(.result$trace_delta_max, .result$n_admissible_evaluated),
        trace_selected_n_bins = .result$trace_selected_n_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_occupancy_failed = .result$trace_occupancy_failed[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_n_pooled_residuals = .result$trace_n_pooled_residuals[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_min_bin_occupancy = .result$trace_min_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_mean_bin_occupancy = .result$trace_mean_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_max_bin_occupancy = .result$trace_max_bin_occupancy[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_sturges_bins = .result$trace_sturges_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_fd_bins = .result$trace_fd_bins[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_shared_residual_range_low = .result$trace_shared_residual_range_low[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        trace_shared_residual_range_high = .result$trace_shared_residual_range_high[, seq_len(.result$n_admissible_evaluated), drop = FALSE]
    )
}
