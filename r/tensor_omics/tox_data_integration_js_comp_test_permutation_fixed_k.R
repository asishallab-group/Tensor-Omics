# Generated. Do not edit.

#' Run the k-NN JS-Comp-Test on one assignment of samples to studies
#'
#' One building block of the study-label Permutation-Test: relabels the samples by `sample_assignment` and runs
#' the whole JS-Comp-Test on the relabelled data set. `expression` holds every study (input layout as in the module
#' description). `sample_assignment` is a permutation of `1..n_samples`: destination sample `k` receives all genes of source
#' sample `sample_assignment(k)`, so every study keeps its slots, sample count and padding. The identity is the observed
#' data, and its result is exactly that of a caller doing the steps below on the observed data.
#'
#' Steps: per study, the gene means and residuals of the relabelled data
#' (\code{\link{compute_gene_means}}, \code{\link{compute_residuals}}); the reference points from
#' the pooled means (\code{\link{pool_means}}); then \code{\link{run_js_comp_test}} with the
#' bin-count options and `divergence_mode` as given.
#'
#' With `rerun_parameter_search = TRUE` the relabelled data first run
#' \code{\link{run_js_comp_test_parameter_search}}
#' with the search options below, and the run uses the `(n_points, n_neighbors)` it returns (also its fallback when no
#' candidate was admissible). `n_points` and `n_neighbors` then only size the work arrays, and `used_n_points` and
#' `used_n_neighbors` report what the run used. Every call runs the whole JS-Comp-Test, with the rerun option the
#' whole parameter search first.
#'
#' Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
#' number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. In the Fortran expert tier, with the rerun option, also when a capacity is below the
#' search's own bound (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`).
#' Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
#' Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
#' values for it.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation_fixed_k::run_js_comp_test_on_sample_assignment}, whose argument names
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
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_stone_equal"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param neighborhood_size_mode a string, one of "neighborhood_size_large", "neighborhood_size_small". Rerun search: which neighborhood size the k-NN candidate grid pairs with each
#'   `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
#'   default: a higher JSCompTest value from small neighborhoods only shows that they are
#'   too small
#'
#'   The default value is `"neighborhood_size_large"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0001`.
#' @param divergence_mode a string, one of "divergence_symmetric", "divergence_pairwise_k2". How each study's per-point divergence is computed, in the run (observed and every
#'   shuffled data set) and in the rerun search
#'
#'   The default value is `"divergence_pairwise_k2"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param grid_points_factor a numeric scalar. Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
#'   grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `4.0`.
#' @param grid_min_points a integer scalar. Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
#'   floor, its last and coarsest `n_points` candidate
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_POINTS`.
#'   The default value is `200`.
#' @return a named list with elements:
#'   \item{global_js_divergence}{a numeric vector. Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
#'     two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
#'     error}
#'   \item{used_n_points}{a integer scalar. Number of reference points the run used (`n_points`, or the search's choice)}
#'   \item{used_n_neighbors}{a integer scalar. Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)}
#' @export
run_js_comp_test_on_sample_assignment <- function(expression, sample_assignment, n_points, n_neighbors, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bin_count_mode = "bin_count_stone_equal", neighborhood_size_mode = "neighborhood_size_large", delta_abs_tolerance = 0.0001, divergence_mode = "divergence_pairwise_k2", grid_points_factor = 4.0, grid_min_points = 200L) {
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
    neighborhood_size_mode <- .tox_as_mode(neighborhood_size_mode, "neighborhood_size_mode", c("neighborhood_size_large", "neighborhood_size_small"))
    delta_abs_tolerance <- .tox_as_double_scalar(delta_abs_tolerance, "delta_abs_tolerance")
    divergence_mode <- .tox_as_mode(divergence_mode, "divergence_mode", c("divergence_symmetric", "divergence_pairwise_k2"))
    grid_points_factor <- .tox_as_double_scalar(grid_points_factor, "grid_points_factor")
    grid_min_points <- .tox_as_integer_scalar(grid_min_points, "grid_min_points")
    .result <- .Call("run_js_comp_test_on_sample_assignment_call", expression, sample_assignment, n_points, n_neighbors, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, bin_count_mode, neighborhood_size_mode, delta_abs_tolerance, divergence_mode, grid_points_factor, grid_min_points)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_samples", "sample_assignment", "n_points", "n_neighbors", "global_js_divergence", "used_n_points", "used_n_neighbors", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "divergence_mode", "grid_points_factor", "grid_min_points", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, "sample_assignment", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence = .result$global_js_divergence,
        used_n_points = .result$used_n_points,
        used_n_neighbors = .result$used_n_neighbors
    )
}

#' Permutation-Test of the k-NN JS-Comp-Test that permutes the samples' study labels
#'
#' Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value (null
#' hypothesis: it does not), with the data, parameters and options of
#' \code{\link{run_js_comp_test_on_sample_assignment}}.
#'
#' 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
#' 2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
#' `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
#' shuffled data. Row `i` of `jsd_null` holds its values.
#' 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.
#'
#' With `rerun_parameter_search = TRUE` every data set, the observed one included, first runs the fixed-k parameter search. `random_seed` then
#' also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
#' draws. The routine is sequential: a caller who wants parallel permutations calls
#' \code{\link{run_js_comp_test_on_sample_assignment}}
#' itself, once per assignment.
#'
#' Failures: a permutation that fails with invalid input (a safeguard: the relabelling keeps every study's sample count, so the k-NN test does not normally trigger it) gets a NaN row in `jsd_null`, is
#' counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
#' `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
#' returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
#' is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own bound
#' (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`) fails the same way with the rerun option.
#' Known limitation: a gene missing in a whole study travels with the samples.
#'
#' Generated from the Fortran procedure \code{tox_data_integration_js_comp_test_permutation_fixed_k::run_js_comp_test_permutation}, whose argument names
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
#' @param bin_count_mode a string, one of "bin_count_occupancy", "bin_count_stone_equal", "bin_count_stone_size_weighted", "bin_count_stone_pooled". How each reference point's histogram bin count is chosen; the Stone modes test every
#'   bin count from `m_min` to `m_max`
#'   (\code{\link{determine_bin_count_stone}}),
#'   with `min_residuals_per_bin` as the floor
#'
#'   The default value is `"bin_count_stone_equal"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param neighborhood_size_mode a string, one of "neighborhood_size_large", "neighborhood_size_small". Rerun search: which neighborhood size the k-NN candidate grid pairs with each
#'   `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
#'   default: a higher JSCompTest value from small neighborhoods only shows that they are
#'   too small
#'
#'   The default value is `"neighborhood_size_large"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param delta_abs_tolerance a numeric scalar. Rerun search (used only with `rerun_parameter_search`): Changes of the JSD of a study up to this absolute amount always count as no change; above it the relative thresholds apply (0: purely relative change).
#'   The minimum valid value is `0.0`.
#'   The default value is `0.0001`.
#' @param divergence_mode a string, one of "divergence_symmetric", "divergence_pairwise_k2". How each study's per-point divergence is computed, in the run (observed and every
#'   shuffled data set) and in the rerun search
#'
#'   The default value is `"divergence_pairwise_k2"`.
#'   Trailing blanks are not part of the value, as in Fortran: they are dropped.
#' @param grid_points_factor a numeric scalar. Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
#'   grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
#'   The minimum valid value is `above(0.0)`.
#'   The default value is `4.0`.
#' @param grid_min_points a integer scalar. Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
#'   floor, its last and coarsest `n_points` candidate
#'   The minimum valid value is `1`.
#'   The maximum valid value is `MAX_POINTS`.
#'   The default value is `200`.
#' @return a named list with elements:
#'   \item{global_js_divergence_observed}{a numeric vector. Weighted global JSD of each study on the observed data (the identity assignment)}
#'   \item{p_values}{a numeric vector. Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
#'     whose value reached the observed one; NaN when none is valid or on error}
#'   \item{jsd_null}{a numeric matrix. Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation}
#'   \item{n_failed_permutations}{a integer scalar. Number of permutations whose JS-Comp-Test failed with invalid input}
#' @export
run_js_comp_test_permutation <- function(expression, n_points, n_neighbors, n_permutations, rerun_parameter_search = FALSE, n_bootstraps = 1000L, join_method = "join_min", min_residuals_per_bin = 10L, min_neighbor_overlap = 0.1, succeeding_ci_overlap = 0.9, plateau_mode = "plateau_effect_size", delta_median_threshold = 0.05, delta_max_threshold = 0.1, delta_epsilon = 1e-10, delta_min_consecutive_transitions = 2L, m_min = 3L, m_max = 120L, gamma_occupancy = 1.25, lower_residual_range_quantile = 0.05, upper_residual_range_quantile = 0.95, two_sided_bootstrapping_significance_level = 2.5, random_seed = 42L, bin_count_mode = "bin_count_stone_equal", neighborhood_size_mode = "neighborhood_size_large", delta_abs_tolerance = 0.0001, divergence_mode = "divergence_pairwise_k2", grid_points_factor = 4.0, grid_min_points = 200L) {
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
    neighborhood_size_mode <- .tox_as_mode(neighborhood_size_mode, "neighborhood_size_mode", c("neighborhood_size_large", "neighborhood_size_small"))
    delta_abs_tolerance <- .tox_as_double_scalar(delta_abs_tolerance, "delta_abs_tolerance")
    divergence_mode <- .tox_as_mode(divergence_mode, "divergence_mode", c("divergence_symmetric", "divergence_pairwise_k2"))
    grid_points_factor <- .tox_as_double_scalar(grid_points_factor, "grid_points_factor")
    grid_min_points <- .tox_as_integer_scalar(grid_min_points, "grid_min_points")
    .result <- .Call("run_js_comp_test_permutation_call", expression, n_points, n_neighbors, n_permutations, rerun_parameter_search, n_bootstraps, join_method, min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, plateau_mode, delta_median_threshold, delta_max_threshold, delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, lower_residual_range_quantile, upper_residual_range_quantile, two_sided_bootstrapping_significance_level, random_seed, bin_count_mode, neighborhood_size_mode, delta_abs_tolerance, divergence_mode, grid_points_factor, grid_min_points)
    .arguments <- c("n_studies", "max_n_genes_all_studies", "max_n_reps_all_studies", "expression", "n_points", "n_neighbors", "n_permutations", "global_js_divergence_observed", "p_values", "jsd_null", "n_failed_permutations", "rerun_parameter_search", "n_bootstraps", "join_method", "min_residuals_per_bin", "min_neighbor_overlap", "succeeding_ci_overlap", "plateau_mode", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "m_min", "m_max", "gamma_occupancy", "lower_residual_range_quantile", "upper_residual_range_quantile", "two_sided_bootstrapping_significance_level", "random_seed", "bin_count_mode", "neighborhood_size_mode", "delta_abs_tolerance", "divergence_mode", "grid_points_factor", "grid_min_points", "ierr")
    .sources <- c("expression", "expression", "expression", NA_character_, NA_character_, NA_character_, "jsd_null", NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_, NA_character_)
    .status <- check_err_code(.result$ierr, .arguments, .sources)

    list(
        global_js_divergence_observed = .result$global_js_divergence_observed,
        p_values = .result$p_values,
        jsd_null = .result$jsd_null,
        n_failed_permutations = .result$n_failed_permutations
    )
}
