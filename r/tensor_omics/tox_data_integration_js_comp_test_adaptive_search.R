# Generated. Do not edit.

#' Generate the ascending adaptive (k_start, k_step, k_max) candidate sequence
#'
#' The adaptive counterpart of
#' \code{\link{generate_js_comp_test_candidates}}.
#' An adaptive candidate fixes only its growth knobs (each neighborhood grows from `k_start`
#' pooled entries in rounds of `k_step` up to `k_max`); its number of reference points emerges.
#' `k_start` shrinks from one candidate to the next, so the candidates **ascend** in the number
#' of reference points they are expected to produce.
#'
#' With `N = max_n_genes_all_studies * n_studies` (the padded pool size) and
#' `floor_k = max(10, 2*n_studies)`:
#'
#' - `k_start_1 = max(floor_k, ceiling(0.02 * N))`
#' - `k_start_t = floor(k_start_1 * 0.8**(t-1))`
#' - stop once `k_start_t < floor_k`, or once 16 candidates are accepted
#' - `k_step_t = ceiling(0.25 * k_start_t)`
#' - `k_max_t = ceiling(4 * k_start_t)`
#'
#' The constants are \code{ADAPTIVE_START_FRACTION} (0.02),
#' \code{ADAPTIVE_GAMMA} (0.8),
#' \code{ADAPTIVE_K_STEP_FRACTION} (0.25),
#' \code{ADAPTIVE_K_MAX_FACTOR} (4),
#' \code{ADAPTIVE_K_START_MIN_ABS} (10) and
#' \code{MAX_CANDIDATE_PAIRS} (16).
#' The values strictly descend without deduplication. At least one candidate is always
#' returned, even for data smaller than `k_start_1`, whose construction then fails with a
#' too-few-means status in the search.
#'
#' Not checked here: the values fit 32-bit integers only while
#' `k_start_1 <= floor(huge(1)/ADAPTIVE_K_MAX_FACTOR)` (\code{ADAPTIVE_K_MAX_FACTOR});
#' \code{\link{calc_adaptive_js_comp_test_bounds}}
#' rejects inputs beyond that.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_adaptive_search::generate_adaptive_js_comp_test_candidates}, whose argument names
#' are the ones an error message reports.
#'
#' @param max_n_genes_all_studies a integer scalar. Maximum number of genes across all studies
#'   The minimum valid value is `1`.
#' @param n_studies a integer scalar. Number of studies
#'   The minimum valid value is `1`.
#' @return a integer matrix. Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
#'   The first `n_candidates` elements will hold the results.
#' @export
generate_adaptive_js_comp_test_candidates <- function(max_n_genes_all_studies, n_studies) {
    max_n_genes_all_studies <- .tox_as_integer_scalar(max_n_genes_all_studies, "max_n_genes_all_studies")
    n_studies <- .tox_as_integer_scalar(n_studies, "n_studies")
    .result <- .Call("generate_adaptive_js_comp_test_candidates_call", max_n_genes_all_studies, n_studies)
    .arguments <- c("max_n_genes_all_studies", "n_studies", "candidates_k_start_k_step_k_max", "n_candidates", "ierr")
    .status <- check_err_code(.result$ierr, .arguments)

    .result$candidates_k_start_k_step_k_max[, seq_len(.result$n_candidates), drop = FALSE]
}

#' Search the ascending adaptive (k_start, k_step, k_max) candidate sequence for a stable JSD parameter setting (Issue #217)
#'
#' The adaptive-neighborhood counterpart of
#' \code{\link{run_js_comp_test_parameter_search}}.
#' Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
#' \code{\link{generate_adaptive_js_comp_test_candidates}}
#' in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
#' The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.
#'
#' **Per candidate**, in order:
#'
#' 1. Grow its neighborhoods exactly as
#' \code{\link{construct_adaptive_neighborhoods}}
#' does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
#' `min_study_neighbors`.
#' 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
#' reference points emerged.
#' 3. Gate 1:
#' \code{\link{check_neighborhood_overlaps}}
#' on the pooled ranges.
#' 4. The pipeline of
#' \code{\link{run_js_comp_test_adaptive}}
#' up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
#' (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
#' count
#' (\code{\link{check_mean_pmf_min_counts}}).
#' 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
#' first plateau `plateau_mode` selects.
#'
#' **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
#' `candidate_n_points` and `candidate_status`, until the search stops. The status is the first
#' failing check, in the order above:
#' \code{ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED}
#' (too few means or residuals),
#' \code{ADAPTIVE_CANDIDATE_EMPTY_STUDY}
#' (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
#' some point),
#' \code{ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED},
#' \code{ADAPTIVE_CANDIDATE_OVERLAP_FAILED},
#' \code{ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED},
#' or else \code{ADAPTIVE_CANDIDATE_EVALUATED}
#' (admissible, with a trace column).
#'
#' **Result**, one of three cases, as in the fixed-k search:
#'
#' 1. A plateau was found, or the sequence had a single candidate and it passed both gates:
#' that candidate's knobs, `plateau_established = TRUE`.
#' 2. No plateau but at least one admissible candidate: the admissible candidate with the
#' smallest median confidence-interval width, `plateau_established = FALSE`.
#' 3. Nothing admissible, including a single candidate that was rejected: candidate 1's knobs,
#' confidence interval `-1.0`, `plateau_established = FALSE`.
#'
#' `n_points` is the selected candidate's logged point count; in case 3 the bin counts are
#' `m_min` and the ranges `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
#' `construct_adaptive_neighborhoods` with the returned knobs and the same construction
#' options, then `run_js_comp_test_adaptive` with the same bin-count options
#' (`bin_count_mode` included): its `global_js_divergence` equals the selected candidate's
#' `trace_global_js_divergence` column.
#'
#' **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
#' place of `trace_n_neighbors`; the per-point traces are jagged in the same way.
#'
#' **Memory.** The per-point arrays are sized by `max_n_points_candidate` (from
#' \code{\link{calc_adaptive_js_comp_test_bounds}});
#' the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
#' `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
#' `max_n_reps_all_studies * N` values each, the size of `residuals`.
#'
#' A failure to create the bootstrap's random number generator is returned in `ierr` and
#' stops the search; the outputs are then undefined.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_adaptive_search::run_js_comp_test_adaptive_parameter_search}, whose argument names
#' are the ones an error message reports.
#'
#' @param gene_means a numeric matrix. Mean expression of every gene in every study, NaN for a missing gene
#'   NaN is permitted for this value.
#' @param residuals a numeric array of rank 3. Signed residuals of every replicate of every gene in every study, NaN for a missing value
#'   NaN is permitted for this value.
#' @param n_bootstraps a integer scalar. Number of bootstraps to perform for a candidate
#'   The minimum valid value is `1`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the studies' confidence-interval overlaps are joined into the plateau condition
#' @param max_n_points_candidate a integer scalar. Most reference points a candidate may emerge with; a candidate with more is rejected
#'   with the capacity status. The per-point outputs and traces below are sized by it, so
#'   the caller has to know it up front
#'   It is recommended to compute this argument from the `max_n_points_candidate` output produced by \code{\link{calc_adaptive_js_comp_test_bounds}}.
#'   The minimum valid value is `1`.
#' @param min_residuals_per_bin a integer scalar. Minimum count every bin must reach, both in the bin-count search (the floor in every
#'   `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
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
#'   Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
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
#' @param tau a numeric scalar. Largest relative dispersion increase a growth round may cause and still be committed
#'   The minimum valid value is `0.0`.
#'   The default value is `0.1`.
#' @param mad_distance_factor a numeric scalar. Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
#'   The minimum valid value is `0.0`.
#'   The default value is `1.0`.
#' @param max_pooled_residuals a integer scalar. Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
#'   The minimum valid value is `0`.
#'   The default value is `0`.
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have for a candidate to be
#'   admissible
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen, for every candidate; the
#'   Stone modes test every bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @return a named list with elements:
#'   \item{k_start}{a integer scalar. The finally chosen candidate's `k_start`}
#'   \item{k_step}{a integer scalar. The finally chosen candidate's `k_step`}
#'   \item{k_max}{a integer scalar. The finally chosen candidate's `k_max`}
#'   \item{n_points}{a integer scalar. Number of reference points the finally chosen candidate's construction produced}
#'   \item{n_bins_per_point}{a integer vector. The chosen candidate's bin count per reference point; the first `n_points` entries
#'     are meaningful. Without an admissible candidate, the first
#'     `min(n_points, max_n_points_candidate)` entries are `m_min`. The rest are 0}
#'   \item{shared_residual_range_low}{a numeric vector. The chosen candidate's lower residual-range bound (R_low) per reference point; the
#'     first `n_points` entries are meaningful, `0.0` without an admissible candidate; the
#'     rest are `0.0`}
#'   \item{shared_residual_range_high}{a numeric vector. The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`}
#'   \item{best_candidate_confidence_interval}{a numeric matrix. The bootstrapped JSD confidence interval of the finally chosen candidate;
#'     `-1.0` throughout where no admissible candidate backs the result}
#'   \item{plateau_established}{a logical scalar. `TRUE` when a plateau was found (by the criterion `plateau_mode` selects) or the
#'     only candidate was admissible; `FALSE` otherwise}
#'   \item{trace_k_start}{a integer vector. Per-admissible-candidate `k_start`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_k_step}{a integer vector. Per-admissible-candidate `k_step`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_k_max}{a integer vector. Per-admissible-candidate `k_max`
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_n_points}{a integer vector. Per-admissible-candidate number of reference points
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_global_js_divergence}{a numeric matrix. Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_lower}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_upper}{a numeric matrix. Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width}{a numeric matrix. Per-admissible-candidate, per-study confidence-interval width
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_width_relative}{a numeric matrix. Per-admissible-candidate, per-study relative confidence-interval width
#'     (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_ci_overlap}{a numeric matrix. Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
#'     interval with the best candidate's so far, relative to its own width; `-1.0` for the
#'     first admissible candidate
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_delta}{a numeric matrix. Per-admissible-candidate, per-study relative JSD change from the previous admissible
#'     candidate; `-1.0` for the first
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
#'   \item{trace_n_pooled_residuals}{a integer matrix. Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_min_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_mean_bin_occupancy}{a numeric matrix. Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_max_bin_occupancy}{a integer matrix. Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_sturges_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
#'     trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_fd_bins}{a integer matrix. Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
#'     as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_low}{a numeric matrix. Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{trace_shared_residual_range_high}{a numeric matrix. Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
#'     jagged as trace_selected_n_bins
#'     The first `n_admissible_evaluated` elements will hold the results.}
#'   \item{candidate_k_start}{a integer vector. `k_start` of every candidate tried, in sequence order
#'     The first `n_candidates_tried` elements will hold the results.}
#'   \item{candidate_n_points}{a integer vector. Number of reference points each tried candidate's construction produced (0 for too
#'     few means)
#'     The first `n_candidates_tried` elements will hold the results.}
#'   \item{candidate_status}{a integer vector. What happened to every tried candidate: `0` evaluated,
#'     `1` construction failed,
#'     `2` a study neighborhood below `min_study_neighbors`
#'     or a study without a non-NaN residual at some point,
#'     `3` capacity exceeded,
#'     `4` overlap gate failed,
#'     `5` consensus-pmf count gate failed
#'     The first `n_candidates_tried` elements will hold the results.}
#' @export
run_js_comp_test_adaptive_parameter_search <- function(gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_ci_overlap", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L, bin_count_mode = "bin_count_occupancy") {
    gene_means <- .tox_as_double_matrix(gene_means, "gene_means")
    residuals <- .tox_as_double_array(residuals, "residuals", 3L)
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
    max_n_points_candidate <- .tox_as_integer_scalar(max_n_points_candidate, "max_n_points_candidate")
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
    tau <- .tox_as_double_scalar(tau, "tau")
    mad_distance_factor <- .tox_as_double_scalar(mad_distance_factor, "mad_distance_factor")
    max_pooled_residuals <- .tox_as_integer_scalar(max_pooled_residuals, "max_pooled_residuals")
    min_study_neighbors <- .tox_as_integer_scalar(min_study_neighbors, "min_study_neighbors")
    bin_count_mode <- .tox_as_mode(bin_count_mode, "bin_count_mode", c("bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled"))
    if (dim(residuals)[3] != dim(gene_means)[2])
        .tox_shape_error("residuals", dim(residuals)[3], "gene_means", dim(gene_means)[2])
    if (dim(residuals)[2] != dim(gene_means)[1])
        .tox_shape_error("residuals", dim(residuals)[2], "gene_means", dim(gene_means)[1])

    .result <- .Call("run_js_comp_test_adaptive_parameter_search_call", gene_means, residuals, n_bootstraps, join_method, max_n_points_candidate, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors, bin_count_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "gene_means", "residuals", "n_bootstraps", "join_method", "max_n_points_candidate", "k_start", "k_step", "k_max", "n_points", "n_bins_per_point", "shared_residual_range_low", "shared_residual_range_high", "best_candidate_confidence_interval", "plateau_established", "n_admissible_evaluated", "trace_k_start", "trace_k_step", "trace_k_max", "trace_n_points", "trace_global_js_divergence", "trace_ci_lower", "trace_ci_upper", "trace_ci_width", "trace_ci_width_relative", "trace_ci_overlap", "trace_delta", "trace_delta_median", "trace_delta_max", "trace_selected_n_bins", "trace_occupancy_failed", "trace_n_pooled_residuals", "trace_min_bin_occupancy", "trace_mean_bin_occupancy", "trace_max_bin_occupancy", "trace_sturges_bins", "trace_fd_bins", "trace_shared_residual_range_low", "trace_shared_residual_range_high", "n_candidates_tried", "candidate_k_start", "candidate_n_points", "candidate_status", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "ierr")
    .sources <- c("gene_means", "gene_means", "residuals", NA_character_, NA_character_, NA_character_, NA_character_, "n_bins_per_point", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        k_start = .result$k_start,
        k_step = .result$k_step,
        k_max = .result$k_max,
        n_points = .result$n_points,
        n_bins_per_point = .result$n_bins_per_point,
        shared_residual_range_low = .result$shared_residual_range_low,
        shared_residual_range_high = .result$shared_residual_range_high,
        best_candidate_confidence_interval = .result$best_candidate_confidence_interval,
        plateau_established = .result$plateau_established,
        trace_k_start = utils::head(.result$trace_k_start, .result$n_admissible_evaluated),
        trace_k_step = utils::head(.result$trace_k_step, .result$n_admissible_evaluated),
        trace_k_max = utils::head(.result$trace_k_max, .result$n_admissible_evaluated),
        trace_n_points = utils::head(.result$trace_n_points, .result$n_admissible_evaluated),
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
        trace_shared_residual_range_high = .result$trace_shared_residual_range_high[, seq_len(.result$n_admissible_evaluated), drop = FALSE],
        candidate_k_start = utils::head(.result$candidate_k_start, .result$n_candidates_tried),
        candidate_n_points = utils::head(.result$candidate_n_points, .result$n_candidates_tried),
        candidate_status = utils::head(.result$candidate_status, .result$n_candidates_tried)
    )
}
