# Generated. Do not edit.

#' Run the k-NN JS-Comp-Test on one assignment of samples to studies
#'
#' One building block of the study-label Permutation-Test: relabels the samples by
#' `sample_assignment` and runs the whole JS-Comp-Test on the relabelled data set.
#'
#' **Input layout.** `expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)`
#' holds every study: gene `g` is the same gene in every study, a gene missing in a study is NaN
#' there, and a replicate slot whose values are all NaN is padding (a study with fewer samples).
#' Every other slot is a sample, numbered in slot order: replicate fastest, then study (sample 1
#' is the first non-padding slot of study 1).
#'
#' **Assignment.** `sample_assignment` is a permutation of `1..n_samples`: destination sample
#' `k` (the `k`-th non-padding slot) receives all genes of source sample `sample_assignment(k)`.
#' So every study keeps its slots, its sample count and its padding; only which samples sit in
#' them changes. The identity `1..n_samples` is the observed data.
#'
#' **Steps.** Per study, the gene means and residuals of the relabelled data
#' (\code{\link{compute_gene_means}},
#' \code{\link{compute_residuals}}); the
#' reference points from the pooled means
#' (\code{\link{pool_means}} on their
#' ascending, NaN-last order); then
#' \code{\link{run_js_comp_test}} with
#' `n_permutations = 0` (its own residual-level permutation test does not run) and the bin-count
#' options as given. On the identity assignment the result is exactly that of a caller doing
#' these steps on the observed data with `run_js_comp_test`.
#'
#' **Rerunning the parameter search.** With `rerun_parameter_search = TRUE`, the relabelled
#' data set first runs
#' \code{\link{run_js_comp_test_parameter_search}}
#' with the search options below, and the run then uses the `(n_points, n_neighbors)` it returns
#' (also its fallback when no candidate was admissible); `n_points` and `n_neighbors` are then
#' not used for the run, only to size the work arrays. `used_n_points` and `used_n_neighbors`
#' report what the run used. Run time: every call runs the whole JS-Comp-Test; with the rerun
#' option it runs the whole parameter search first.
#'
#' **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
#' the samples, so after relabelling that study may have values for it.
#'
#' Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
#' no sample at all, when `n_samples` is not the number of samples in `expression`, or when
#' `sample_assignment` is not a permutation of `1..n_samples`. (Fortran expert tier only, where
#' the caller supplies the work-array capacities: with the rerun option it also fails when a
#' capacity is below the search's own bound, `calc_js_comp_test_candidate_bounds`.) Errors of
#' the search and the run are returned with their code (for example `ERR_INVALID_INPUT` when some
#' study has no residual at some reference point of the relabelled data). `global_js_divergence`
#' is then NaN.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation::run_js_comp_test_on_sample_assignment}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression a numeric array of rank 3. Expression of every sample slot of every gene in every study (see above); NaN for a
#'   missing value, a missing gene or a padding slot
#'   NaN is permitted for this value.
#' @param sample_assignment a integer vector. A permutation of `1..n_samples`: destination sample `k` receives source sample
#'   `sample_assignment(k)`; the identity is the observed data
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_samples`.
#' @param n_points a integer scalar. Number of reference points of the run (sizes the work arrays only, with the rerun option)
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
#'   rerun option)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param rerun_parameter_search a logical scalar. `TRUE`: run the fixed-k parameter search on the relabelled data first and use its
#'   `(n_points, n_neighbors)`; `FALSE`: use `n_points` and `n_neighbors`
#'   The default value is `FALSE`.
#' @param n_bootstraps a integer scalar. Number of bootstraps per candidate of the rerun parameter search
#'   The minimum valid value is `1`.
#'   The default value is `1000`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the rerun search joins the studies' confidence-interval overlaps
#'
#'   The default value is `"join_min"`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach in the bin-count search (and,
#'   in the rerun search, in every consensus bin)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Rerun search: which plateau criterion stops the search; by default, as in the
#'   searches, the change of the observed JSD values from one admissible candidate to the next
#'
#'   The default value is `"plateau_effect_size"`.
#' @param delta_median_threshold a numeric scalar. Rerun search: upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Rerun search: upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Rerun search: floor of the denominators of the relative JSD change and CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Rerun search: two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Rerun search: seed of the bootstrap's random number generator
#'   The default value is `42`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @param neighborhood_size_mode a string, one of "neighborhood_size_large", "neighborhood_size_small", "neighborhood_size_both". Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
#'   `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
#'   by default: a higher JSCompTest value from small neighborhoods only shows that they are
#'   too small
#'
#'   The default value is `"neighborhood_size_large"`.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0`.
#' @return a named list with elements:
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus, for this assignment; NaN on error}
#'   \item{used_n_points}{a integer scalar. Number of reference points the run used (`n_points`, or the search's choice)}
#'   \item{used_n_neighbors}{a integer scalar. Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)}
#' @export
run_js_comp_test_on_sample_assignment <- function(expression, sample_assignment, n_points, n_neighbors, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bin_count_mode = "bin_count_occupancy", neighborhood_size_mode = "neighborhood_size_large", delta_abs_tolerance = 0.0) {
    expression <- .tox_as_double_array(expression, "expression", 3L)
    sample_assignment <- .tox_as_integer_vector(sample_assignment, "sample_assignment")
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    rerun_parameter_search <- .tox_as_logical_scalar(rerun_parameter_search, "rerun_parameter_search")
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
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
    neighborhood_size_mode <- .tox_as_mode(neighborhood_size_mode, "neighborhood_size_mode", c("neighborhood_size_large", "neighborhood_size_small", "neighborhood_size_both"))
    delta_abs_tolerance <- .tox_as_double_scalar(delta_abs_tolerance, "delta_abs_tolerance")
    .result <- .Call("run_js_comp_test_on_sample_assignment_call", expression, sample_assignment, n_points, n_neighbors, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, bin_count_mode, neighborhood_size_mode, delta_abs_tolerance)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "n_points", "n_neighbors", "global_js_divergence", "used_n_points", "used_n_neighbors", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, "sample_assignment", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence = .result$global_js_divergence,
        used_n_points = .result$used_n_points,
        used_n_neighbors = .result$used_n_neighbors
    )
}

#' Run the adaptive JS-Comp-Test on one assignment of samples to studies
#'
#' The adaptive counterpart of
#' \code{\link{run_js_comp_test_on_sample_assignment}},
#' with the same input layout and assignment convention: `expression` holds every study (NaN for
#' a missing gene, all-NaN replicate slots are padding, samples numbered replicate fastest, then
#' study), and destination sample `k` receives all genes of source sample `sample_assignment(k)`;
#' the identity is the observed data.
#'
#' **Steps.** Per study, the gene means and residuals of the relabelled data; the adaptive
#' neighborhoods for the growth knobs `k_start`, `k_step`, `k_max` and the construction options
#' (\code{\link{construct_adaptive_neighborhoods}});
#' then \code{\link{run_js_comp_test_adaptive}}
#' on its ranges with `n_permutations = 0`. `n_points` is the number of neighborhoods the
#' construction built. On the identity assignment the result is exactly that of a caller doing
#' these steps on the observed data.
#'
#' **Rerunning the parameter search.** With `rerun_parameter_search = TRUE`, the relabelled
#' data first runs
#' \code{\link{run_js_comp_test_adaptive_parameter_search}}
#' (with its own capacity, `calc_adaptive_js_comp_test_bounds`, so it rejects the same candidates
#' as a direct call), and the construction and run then use the knobs it returns (also its
#' fallback when no candidate was admissible); `k_start`, `k_step` and `k_max` are then not used,
#' except `k_start` to size the work arrays. `used_k_start`, `used_k_step` and `used_k_max` report
#' what was used. Run time: every call runs the whole JS-Comp-Test; with the rerun option it runs
#' the whole parameter search first.
#'
#' **Known limitation.** A gene missing in a whole study (NaN in all its samples) travels with
#' the samples, so after relabelling that study may have values for it.
#'
#' Errors (`ierr`): the call fails with invalid input (`ERR_INVALID_INPUT`) when some study has
#' no sample at all, when `n_samples` is not the number of samples in `expression`, when
#' `sample_assignment` is not a permutation of `1..n_samples`, or when the construction builds no
#' neighborhood or more neighborhoods than the work arrays hold. (Fortran expert tier only, where
#' the caller supplies the work-array capacity: with the rerun option it also fails when the
#' capacity is below the search's own.) Errors of the search and the run are returned with their
#' code. `global_js_divergence` is then NaN.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation::run_js_comp_test_adaptive_on_sample_assignment}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression a numeric array of rank 3. Expression of every sample slot of every gene in every study (see above); NaN for a
#'   missing value, a missing gene or a padding slot
#'   NaN is permitted for this value.
#' @param sample_assignment a integer vector. A permutation of `1..n_samples`: destination sample `k` receives source sample
#'   `sample_assignment(k)`; the identity is the observed data
#'   The minimum valid value is `1`.
#'   The maximum valid value is `n_samples`.
#' @param k_start a integer scalar. Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
#'   The minimum valid value is `1`.
#' @param k_step a integer scalar. Pooled entries staged per adaptive growth round
#'   The minimum valid value is `1`.
#' @param k_max a integer scalar. Largest number of pooled entries in one neighborhood
#'   The minimum valid value is `k_start`.
#' @param rerun_parameter_search a logical scalar. `TRUE`: run the adaptive parameter search on the relabelled data first and use its
#'   knobs; `FALSE`: use `k_start`, `k_step` and `k_max`
#'   The default value is `FALSE`.
#' @param n_bootstraps a integer scalar. Number of bootstraps per candidate of the rerun parameter search
#'   The minimum valid value is `1`.
#'   The default value is `1000`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the rerun search joins the studies' confidence-interval overlaps
#'
#'   The default value is `"join_min"`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach in the bin-count search (and,
#'   in the rerun search, in every consensus bin)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Rerun search: which plateau criterion stops the search; by default, as in the
#'   searches, the change of the observed JSD values from one admissible candidate to the next
#'
#'   The default value is `"plateau_effect_size"`.
#' @param delta_median_threshold a numeric scalar. Rerun search: upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Rerun search: upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Rerun search: floor of the denominators of the relative JSD change and CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Rerun search: two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Rerun search: seed of the bootstrap's random number generator
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
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have (construction status; in
#'   the rerun search, candidate admissibility)
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0`.
#' @return a named list with elements:
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study against the consensus, for this assignment; NaN on error}
#'   \item{used_k_start}{a integer scalar. `k_start` the construction used (`k_start`, or the search's choice)}
#'   \item{used_k_step}{a integer scalar. `k_step` the construction used (`k_step`, or the search's choice)}
#'   \item{used_k_max}{a integer scalar. `k_max` the construction used (`k_max`, or the search's choice)}
#'   \item{n_points}{a integer scalar. Number of neighborhoods the construction built (0 when it built none)}
#' @export
run_js_comp_test_adaptive_on_sample_assignment <- function(expression, sample_assignment, k_start, k_step, k_max, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L, bin_count_mode = "bin_count_occupancy", delta_abs_tolerance = 0.0) {
    expression <- .tox_as_double_array(expression, "expression", 3L)
    sample_assignment <- .tox_as_integer_vector(sample_assignment, "sample_assignment")
    k_start <- .tox_as_integer_scalar(k_start, "k_start")
    k_step <- .tox_as_integer_scalar(k_step, "k_step")
    k_max <- .tox_as_integer_scalar(k_max, "k_max")
    rerun_parameter_search <- .tox_as_logical_scalar(rerun_parameter_search, "rerun_parameter_search")
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
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
    delta_abs_tolerance <- .tox_as_double_scalar(delta_abs_tolerance, "delta_abs_tolerance")
    .result <- .Call("run_js_comp_test_adaptive_on_sample_assignment_call", expression, sample_assignment, k_start, k_step, k_max, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors, bin_count_mode, delta_abs_tolerance)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "k_start", "k_step", "k_max", "global_js_divergence", "used_k_start", "used_k_step", "used_k_max", "n_points", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "delta_abs_tolerance", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, "sample_assignment", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence = .result$global_js_divergence,
        used_k_start = .result$used_k_start,
        used_k_step = .result$used_k_step,
        used_k_max = .result$used_k_max,
        n_points = .result$n_points
    )
}

#' Study-label Permutation-Test of the k-NN JS-Comp-Test
#'
#' Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value
#' (null hypothesis: it does not). The same data, parameters and options as
#' \code{\link{run_js_comp_test_on_sample_assignment}},
#' whose input layout (all studies in `expression`, NaN for a missing gene, all-NaN replicate
#' slots are padding) and assignment convention apply here too.
#'
#' 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
#' 2. `n_permutations` times, the study labels of the samples are shuffled without replacement (a
#' uniformly random permutation, unbiased Fisher-Yates, from one random number generator seeded
#' with `random_seed`); every study keeps its sample count. The whole JS-Comp-Test runs on the
#' shuffled data; row `i` of `jsd_null` holds its values.
#' 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`,
#' the Laplace correction of the residual-level test, so it is never 0 and at least
#' `1/(n_permutations + 1)`.
#'
#' **Rerunning the parameter search** (`rerun_parameter_search = TRUE`): every data set, the
#' observed one included, first runs the fixed-k parameter search and then the run with the
#' `(n_points, n_neighbors)` it selects; `n_points` and `n_neighbors` then only size the work
#' arrays. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
#' runs the whole parameter search first. The routine is sequential: a caller who wants
#' permutations in parallel calls
#' \code{\link{run_js_comp_test_on_sample_assignment}}
#' itself with its own assignments (e.g. from R), one call per assignment.
#'
#' **Failures.** A permutation whose JS-Comp-Test fails with invalid input (for example a study
#' left without residuals at some reference point) gets a NaN row in `jsd_null`, is counted in
#' `n_failed_permutations` and does not count in `n_valid`. No valid permutation gives NaN
#' p-values (a direct call with `n_permutations = 0` too; the published entry points reject 0,
#' since it is the extent of `jsd_null`). Any other error ends the test and is returned, with NaN
#' p-values; the rows of the permutations run so far are kept. An error on the observed data is
#' returned at once (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study
#' without any sample (the call fails with invalid input). In the Fortran expert tier, where the
#' caller supplies the work-array capacities, a capacity below the search's own bound fails the
#' same way with the rerun option.
#'
#' **Known limitation.** A gene missing in a whole study travels with the samples, so after
#' shuffling that study may have values for it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation::run_js_comp_test_label_permutation_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression a numeric array of rank 3. Expression of every sample slot of every gene in every study (see above); NaN for a
#'   missing value, a missing gene or a padding slot
#'   NaN is permitted for this value.
#' @param n_points a integer scalar. Number of reference points of the run (sizes the work arrays only, with the rerun option)
#'   The minimum valid value is `1`.
#' @param n_neighbors a integer scalar. Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
#'   rerun option)
#'   The minimum valid value is `1`.
#'   The maximum valid value is `max_n_genes_all_studies`.
#' @param n_permutations a integer scalar. Number of random assignments (rows of `jsd_null`)
#'   The minimum valid value is `1`.
#' @param rerun_parameter_search a logical scalar. `TRUE`: run the fixed-k parameter search on the relabelled data first and use its
#'   `(n_points, n_neighbors)`; `FALSE`: use `n_points` and `n_neighbors`
#'   The default value is `FALSE`.
#' @param n_bootstraps a integer scalar. Number of bootstraps per candidate of the rerun parameter search
#'   The minimum valid value is `1`.
#'   The default value is `1000`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the rerun search joins the studies' confidence-interval overlaps
#'
#'   The default value is `"join_min"`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach in the bin-count search (and,
#'   in the rerun search, in every consensus bin)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Rerun search: which plateau criterion stops the search; by default, as in the
#'   searches, the change of the observed JSD values from one admissible candidate to the next
#'
#'   The default value is `"plateau_effect_size"`.
#' @param delta_median_threshold a numeric scalar. Rerun search: upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Rerun search: upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Rerun search: floor of the denominators of the relative JSD change and CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Rerun search: two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed of the random number generator of the label draws, and of the bootstrap of
#'   every rerun parameter search (a separate generator)
#'   The default value is `42`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @param neighborhood_size_mode a string, one of "neighborhood_size_large", "neighborhood_size_small", "neighborhood_size_both". Rerun search: which neighborhood sizes the k-NN candidate grid pairs with each
#'   `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are left out
#'   by default: a higher JSCompTest value from small neighborhoods only shows that they are
#'   too small
#'
#'   The default value is `"neighborhood_size_large"`.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0`.
#' @return a named list with elements:
#'   \item{global_js_divergence_observed}{a numeric vector. Weighted global JSD of each study on the observed data (the identity assignment)}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
#'     whose value reached the observed one; NaN when none is valid or on error}
#'   \item{jsd_null}{a numeric matrix. Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation}
#'   \item{n_failed_permutations}{a integer scalar. Number of permutations whose JS-Comp-Test failed with invalid input}
#' @export
run_js_comp_test_label_permutation_test <- function(expression, n_points, n_neighbors, n_permutations, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bin_count_mode = "bin_count_occupancy", neighborhood_size_mode = "neighborhood_size_large", delta_abs_tolerance = 0.0) {
    expression <- .tox_as_double_array(expression, "expression", 3L)
    n_points <- .tox_as_integer_scalar(n_points, "n_points")
    n_neighbors <- .tox_as_integer_scalar(n_neighbors, "n_neighbors")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    rerun_parameter_search <- .tox_as_logical_scalar(rerun_parameter_search, "rerun_parameter_search")
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
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
    neighborhood_size_mode <- .tox_as_mode(neighborhood_size_mode, "neighborhood_size_mode", c("neighborhood_size_large", "neighborhood_size_small", "neighborhood_size_both"))
    delta_abs_tolerance <- .tox_as_double_scalar(delta_abs_tolerance, "delta_abs_tolerance")
    .result <- .Call("run_js_comp_test_label_permutation_test_call", expression, n_points, n_neighbors, n_permutations, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, bin_count_mode, neighborhood_size_mode, delta_abs_tolerance)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_points", "n_neighbors", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, NA_character_, NA_character_, "jsd_null", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence_observed = .result$global_js_divergence_observed,
        p_values = .result$p_values,
        jsd_null = .result$jsd_null,
        n_failed_permutations = .result$n_failed_permutations
    )
}

#' Study-label Permutation-Test of the adaptive JS-Comp-Test
#'
#' The adaptive counterpart of
#' \code{\link{run_js_comp_test_label_permutation_test}}:
#' the same test (observed value, `n_permutations` shuffles of the study labels without
#' replacement, p-value `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`), each
#' data set run by
#' \code{\link{run_js_comp_test_adaptive_on_sample_assignment}}
#' with the same data, knobs and options; its input layout and assignment convention apply here too.
#'
#' **Rerunning the parameter search** (`rerun_parameter_search = TRUE`): every data set first
#' runs the adaptive parameter search and then the construction and run with the knobs it
#' selects. **Run time:** each permutation runs the whole JS-Comp-Test; with the rerun option it
#' runs the whole parameter search first. The routine is sequential: a caller who wants
#' permutations in parallel calls
#' \code{\link{run_js_comp_test_adaptive_on_sample_assignment}}
#' itself with its own assignments (e.g. from R), one call per assignment.
#'
#' **Failures.** A permutation that fails with invalid input (a study left without residuals at
#' some reference point, or a construction that builds no neighborhood or more neighborhoods
#' than the work arrays hold) gets a NaN row in `jsd_null`, is counted in `n_failed_permutations`
#' and does not count in `n_valid`. No valid permutation gives NaN p-values (a direct call with
#' `n_permutations = 0` too; the published entry points reject 0, since it is the extent of
#' `jsd_null`). Any other error ends the test and is returned, with NaN p-values; the rows of the
#' permutations run so far are kept. An error on the observed data is returned at once
#' (`global_js_divergence_observed`, `p_values` and `jsd_null` NaN), as is a study without any
#' sample (the call fails with invalid input). In the Fortran expert tier, where the caller
#' supplies the work-array capacity, a capacity below the search's own fails the same way with
#' the rerun option.
#'
#' **Known limitation.** A gene missing in a whole study travels with the samples, so after
#' shuffling that study may have values for it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation::run_js_comp_test_adaptive_label_permutation_test}, whose argument names
#' are the ones an error message reports.
#'
#' @param expression a numeric array of rank 3. Expression of every sample slot of every gene in every study (see above); NaN for a
#'   missing value, a missing gene or a padding slot
#'   NaN is permitted for this value.
#' @param k_start a integer scalar. Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
#'   The minimum valid value is `1`.
#' @param k_step a integer scalar. Pooled entries staged per adaptive growth round
#'   The minimum valid value is `1`.
#' @param k_max a integer scalar. Largest number of pooled entries in one neighborhood
#'   The minimum valid value is `k_start`.
#' @param n_permutations a integer scalar. Number of random assignments (rows of `jsd_null`)
#'   The minimum valid value is `1`.
#' @param rerun_parameter_search a logical scalar. `TRUE`: run the adaptive parameter search on the relabelled data first and use its
#'   knobs; `FALSE`: use `k_start`, `k_step` and `k_max`
#'   The default value is `FALSE`.
#' @param n_bootstraps a integer scalar. Number of bootstraps per candidate of the rerun parameter search
#'   The minimum valid value is `1`.
#'   The default value is `1000`.
#' @param join_method a string, one of "join_min", "join_max", "join_median". How the rerun search joins the studies' confidence-interval overlaps
#'
#'   The default value is `"join_min"`.
#' @param min_residuals_per_bin a integer scalar. Minimum number of pooled residuals every bin must reach in the bin-count search (and,
#'   in the rerun search, in every consensus bin)
#'   The minimum valid value is `0`.
#'   The default value is `10`.
#' @param min_neighbor_overlap a numeric scalar. Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.1`.
#' @param succeeding_ci_overlap a numeric scalar. Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `1.0`.
#'   The default value is `0.9`.
#' @param plateau_mode a string, one of "plateau_ci_overlap", "plateau_effect_size", "plateau_both". Rerun search: which plateau criterion stops the search; by default, as in the
#'   searches, the change of the observed JSD values from one admissible candidate to the next
#'
#'   The default value is `"plateau_effect_size"`.
#' @param delta_median_threshold a numeric scalar. Rerun search: upper bound on the median relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.05`.
#' @param delta_max_threshold a numeric scalar. Rerun search: upper bound on the largest relative JSD change for a transition to qualify
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `0.10`.
#' @param delta_epsilon a numeric scalar. Rerun search: floor of the denominators of the relative JSD change and CI width
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `1.0e-10`.
#' @param delta_min_consecutive_transitions a integer scalar. Rerun search: number of consecutive qualifying transitions for an effect-size plateau
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
#' @param two_sided_bootstrapping_significance_level a numeric scalar. Rerun search: two-sided significance level of the bootstrap interval, in percent per end
#'   The minimum valid value is `0.0`.
#'   The maximum valid value is `100.0`.
#'   The default value is `2.5`.
#' @param random_seed a integer scalar. Seed of the random number generator of the label draws, and of the bootstrap of
#'   every rerun parameter search (a separate generator)
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
#' @param min_study_neighbors a integer scalar. Fewest entries of each study every neighborhood must have (construction status; in
#'   the rerun search, candidate admissibility)
#'   The minimum valid value is `1`.
#'   The default value is `1`.
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_occupancy"`.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0`.
#' @return a named list with elements:
#'   \item{global_js_divergence_observed}{a numeric vector. Weighted global JSD of each study on the observed data (the identity assignment)}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
#'     whose value reached the observed one; NaN when none is valid or on error}
#'   \item{jsd_null}{a numeric matrix. Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation}
#'   \item{n_failed_permutations}{a integer scalar. Number of permutations whose JS-Comp-Test failed with invalid input}
#' @export
run_js_comp_test_adaptive_label_permutation_test <- function(expression, k_start, k_step, k_max, n_permutations, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L, bin_count_mode = "bin_count_occupancy", delta_abs_tolerance = 0.0) {
    expression <- .tox_as_double_array(expression, "expression", 3L)
    k_start <- .tox_as_integer_scalar(k_start, "k_start")
    k_step <- .tox_as_integer_scalar(k_step, "k_step")
    k_max <- .tox_as_integer_scalar(k_max, "k_max")
    n_permutations <- .tox_as_integer_scalar(n_permutations, "n_permutations")
    rerun_parameter_search <- .tox_as_logical_scalar(rerun_parameter_search, "rerun_parameter_search")
    n_bootstraps <- .tox_as_integer_scalar(n_bootstraps, "n_bootstraps")
    join_method <- .tox_as_mode(join_method, "join_method", c("join_min", "join_max", "join_median"))
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
    delta_abs_tolerance <- .tox_as_double_scalar(delta_abs_tolerance, "delta_abs_tolerance")
    .result <- .Call("run_js_comp_test_adaptive_label_permutation_test_call", expression, k_start, k_step, k_max, n_permutations, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors, bin_count_mode, delta_abs_tolerance)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "k_start", "k_step", "k_max", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "delta_abs_tolerance", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, NA_character_, NA_character_, NA_character_, "jsd_null", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence_observed = .result$global_js_divergence_observed,
        p_values = .result$p_values,
        jsd_null = .result$jsd_null,
        n_failed_permutations = .result$n_failed_permutations
    )
}
