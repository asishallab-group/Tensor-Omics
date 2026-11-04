# Generated. Do not edit.

#' Run the adaptive JS-Comp-Test on one assignment of samples to studies
#'
#' The adaptive counterpart of
#' \code{\link{run_js_comp_test_on_sample_assignment}},
#' with the same input layout and assignment convention (the identity is the observed data).
#'
#' Steps: per study, the gene means and residuals of the relabelled data; the adaptive neighborhoods for the
#' growth knobs `k_start`, `k_step`, `k_max` and the construction options
#' (\code{\link{construct_adaptive_neighborhoods}}); then
#' \code{\link{run_js_comp_test_adaptive}} on its ranges, with `divergence_mode` as given. `n_points` is the number
#' of neighborhoods the construction built.
#'
#' With `rerun_parameter_search = TRUE` the relabelled data first run
#' \code{\link{run_js_comp_test_adaptive_parameter_search}}
#' (with its own capacity, so it rejects the same candidates as a direct call), and the construction and run use the knobs
#' it returns (also its fallback when no candidate was admissible). `k_start`, `k_step` and `k_max` are then not used, except
#' `k_start` to size the work arrays, and `used_k_start`, `used_k_step`, `used_k_max` report what was used. Every call runs the
#' whole JS-Comp-Test, with the rerun option the whole parameter search first.
#'
#' Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
#' number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. Also when the construction does not succeed: its status is not ok (too few means, a
#' neighborhood with too few residuals or with fewer than `min_study_neighbors` entries of some study), it builds no neighborhood, or more
#' neighborhoods than the work arrays hold. In the Fortran expert tier, with the rerun option, also when the capacity is below
#' the search's own.
#' Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
#' Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
#' values for it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation_adaptive::run_js_comp_test_adaptive_on_sample_assignment}, whose argument names
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
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
#'   The default value is `"bin_count_stone_equal"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0001`.
#' @param divergence_mode a string, one of "divergence_symmetric", "divergence_pairwise_k2". How each study's per-point divergence is computed, in the run (observed and every
#'   shuffled data set) and in the rerun search
#'
#'   The default value is `"divergence_pairwise_k2"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @return a named list with elements:
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
#'     two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
#'     error}
#'   \item{used_k_start}{a integer scalar. `k_start` the construction used (`k_start`, or the search's choice)}
#'   \item{used_k_step}{a integer scalar. `k_step` the construction used (`k_step`, or the search's choice)}
#'   \item{used_k_max}{a integer scalar. `k_max` the construction used (`k_max`, or the search's choice)}
#'   \item{n_points}{a integer scalar. Number of neighborhoods the construction built (0 when it built none)}
#' @export
run_js_comp_test_adaptive_on_sample_assignment <- function(expression, sample_assignment, k_start, k_step, k_max, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L, bin_count_mode = "bin_count_stone_equal", delta_abs_tolerance = 0.0001, divergence_mode = "divergence_pairwise_k2") {
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
    divergence_mode <- .tox_as_mode(divergence_mode, "divergence_mode", c("divergence_symmetric", "divergence_pairwise_k2"))
    .result <- .Call("run_js_comp_test_adaptive_on_sample_assignment_call", expression, sample_assignment, k_start, k_step, k_max, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors, bin_count_mode, delta_abs_tolerance, divergence_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "k_start", "k_step", "k_max", "global_js_divergence", "used_k_start", "used_k_step", "used_k_max", "n_points", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "delta_abs_tolerance", "divergence_mode", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, "sample_assignment", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence = .result$global_js_divergence,
        used_k_start = .result$used_k_start,
        used_k_step = .result$used_k_step,
        used_k_max = .result$used_k_max,
        n_points = .result$n_points
    )
}

#' Permutation-Test of the adaptive JS-Comp-Test that permutes the samples' study labels
#'
#' The adaptive counterpart of
#' \code{\link{run_js_comp_test_permutation}}:
#' the same test, each data set run by
#' \code{\link{run_js_comp_test_adaptive_on_sample_assignment}}
#' with the same data, knobs and options.
#'
#' 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
#' 2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
#' `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
#' shuffled data. Row `i` of `jsd_null` holds its values.
#' 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.
#'
#' With `rerun_parameter_search = TRUE` every data set, the observed one included, first runs the adaptive parameter search. `random_seed` then
#' also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
#' draws. The routine is sequential: a caller who wants parallel permutations calls
#' \code{\link{run_js_comp_test_adaptive_on_sample_assignment}}
#' itself, once per assignment.
#'
#' Failures: a permutation that fails with invalid input (a study left without residuals at some reference point, or a construction that does not
#' succeed: a status other than ok, no neighborhood, or more neighborhoods than the work arrays hold) gets a NaN row in `jsd_null`, is
#' counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
#' `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
#' returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
#' is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own fails the same way with the rerun option.
#' Known limitation: a gene missing in a whole study travels with the samples.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation_adaptive::run_js_comp_test_adaptive_permutation}, whose argument names
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
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
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
#'   The default value is `"bin_count_stone_equal"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0001`.
#' @param divergence_mode a string, one of "divergence_symmetric", "divergence_pairwise_k2". How each study's per-point divergence is computed, in the run (observed and every
#'   shuffled data set) and in the rerun search
#'
#'   The default value is `"divergence_pairwise_k2"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @return a named list with elements:
#'   \item{global_js_divergence_observed}{a numeric vector. Weighted global JSD of each study on the observed data (the identity assignment)}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
#'     whose value reached the observed one; NaN when none is valid or on error}
#'   \item{jsd_null}{a numeric matrix. Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation}
#'   \item{n_failed_permutations}{a integer scalar. Number of permutations whose JS-Comp-Test failed with invalid input}
#' @export
run_js_comp_test_adaptive_permutation <- function(expression, k_start, k_step, k_max, n_permutations, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, tau = 0.1, mad_distance_factor = 1.0, max_pooled_residuals = 0L, min_study_neighbors = 1L, bin_count_mode = "bin_count_stone_equal", delta_abs_tolerance = 0.0001, divergence_mode = "divergence_pairwise_k2") {
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
    divergence_mode <- .tox_as_mode(divergence_mode, "divergence_mode", c("divergence_symmetric", "divergence_pairwise_k2"))
    .result <- .Call("run_js_comp_test_adaptive_permutation_call", expression, k_start, k_step, k_max, n_permutations, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, tau, mad_distance_factor, max_pooled_residuals, min_study_neighbors, bin_count_mode, delta_abs_tolerance, divergence_mode)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "k_start", "k_step", "k_max", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "tau", "mad_distance_factor", "max_pooled_residuals", "min_study_neighbors", "bin_count_mode", "delta_abs_tolerance", "divergence_mode", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, NA_character_, NA_character_, NA_character_, "jsd_null", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence_observed = .result$global_js_divergence_observed,
        p_values = .result$p_values,
        jsd_null = .result$jsd_null,
        n_failed_permutations = .result$n_failed_permutations
    )
}
