#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_fixed_k_search(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Fixed-k Parameter Search
!|
!| The data-driven choice of the JS-Comp-Test's fixed-k neighborhoods, `n_points` reference
!| points of `n_neighbors` genes each.
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| walks the `(n_points, n_neighbors)` candidate grid of
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| from finest to coarsest and stops at the first admissible candidate whose bootstrapped JSDs
!| have reached a plateau, by the criteria of
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_fixed_k_search_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: generate_js_comp_test_candidates_c
    public :: run_js_comp_test_parameter_search_c
    public :: run_js_comp_test_parameter_search_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_fixed_k_search(module):generate_js_comp_test_candidates(subroutine)]]
    !| Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
    !| `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
    !| `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
    !| `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
    !| Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
    !| candidate. Bin counts are not part of a candidate; they are chosen per reference point later.
    subroutine generate_js_comp_test_candidates_c(&
            max_n_genes_all_studies,&
            candidates_n_points_n_neighbors,&
            n_candidates,&
            ierr&
        ) bind(C, name="generate_js_comp_test_candidates_c")
        use tox_data_integration_js_comp_test_fixed_k_search, only: generate_js_comp_test_candidates

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(2, 16), intent(out), target :: candidates_n_points_n_neighbors
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates
            !! Number of candidate pairs actually filled (at most 16)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_candidates)
        M_CHECK_ARRAY_NON_NULL(candidates_n_points_n_neighbors, 2 * 16)

        call generate_js_comp_test_candidates(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            candidates_n_points_n_neighbors = candidates_n_points_n_neighbors,&
            n_candidates = n_candidates,&
            ierr = ierr&
        )
    end subroutine generate_js_comp_test_candidates_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_fixed_k_search(module):run_js_comp_test_parameter_search(subroutine)]]
    !| Walks the candidates of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !| study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !| ([[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !| residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !| `min_residuals_per_bin`
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !| ([[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]]),
    !| and both plateau criteria: CI overlap with the best candidate so far
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]])
    !| and relative effect-size stability
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
    !| Both are always computed and traced; `plateau_mode` decides which one stops the search. A
    !| stop by the effect-size criterion returns the candidate that triggered it.
    !|
    !| The returned candidate:
    !|
    !| 1. A plateau was found, or the grid has a single candidate and it is admissible: the best
    !| candidate, `plateau_established = .true.`.
    !| 2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
    !| median confidence-interval width across studies, `plateau_established = .false.`.
    !| 3. No candidate was admissible: the first (finest) candidate, with
    !| `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
    !| ranges `0.0`, and `plateau_established = .false.`.
    !|
    !| Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
    !| for its meaningful part. A failure to create the bootstrap's random number generator is
    !| returned in `ierr` and stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_parameter_search_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            n_bootstraps,&
            join_method,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_points,&
            n_neighbors,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            best_candidate_pair_confidence_interval,&
            plateau_established,&
            n_admissible_evaluated,&
            trace_n_points,&
            trace_n_neighbors,&
            trace_global_js_divergence,&
            trace_ci_lower,&
            trace_ci_upper,&
            trace_ci_width,&
            trace_ci_width_relative,&
            trace_ci_overlap,&
            trace_delta,&
            trace_delta_median,&
            trace_delta_max,&
            trace_selected_n_bins,&
            trace_occupancy_failed,&
            trace_n_pooled_residuals,&
            trace_min_bin_occupancy,&
            trace_mean_bin_occupancy,&
            trace_max_bin_occupancy,&
            trace_sturges_bins,&
            trace_fd_bins,&
            trace_shared_residual_range_low,&
            trace_shared_residual_range_high,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            bin_count_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_parameter_search_c")
        use tox_data_integration_js_comp_test_fixed_k_search, only: run_js_comp_test_parameter_search
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_points_candidate
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
        integer(c_int), intent(out), target :: n_points
            !! The finally chosen candidate's `n_points`
        integer(c_int), intent(out), target :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` when no candidate was admissible
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
            !! The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
            !! was admissible
        logical(c_bool), intent(out), target :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise, also when the smallest-uncertainty
            !! fallback returns a real candidate and interval
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array (an inadmissible candidate has no entry)
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate `n_points`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width (the values
            !! `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out), target :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            !! like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
            !! one threshold for both, so a point the search marks as failed also fails the gate
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                | Value                                                                                        |
            !! |-------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_ci_overlap'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                                   |
            !! |----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_occupancy'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f
        integer(int32) :: bin_count_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(plateau_established)
        M_CHECK_NON_NULL(n_admissible_evaluated)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(best_candidate_pair_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(trace_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_neighbors, 16)
        M_CHECK_ARRAY_NON_NULL(trace_global_js_divergence, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_lower, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_upper, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_width, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_width_relative, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_overlap, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_median, 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_max, 16)
        M_CHECK_ARRAY_NON_NULL(trace_selected_n_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_occupancy_failed, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_pooled_residuals, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_min_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_mean_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_max_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_sturges_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_fd_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_shared_residual_range_low, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_shared_residual_range_high, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

        block
            character(len=:), pointer :: join_method_f
            join_method_f => c_char_as_view(join_method)

            select case (join_method_f)
                case ("join_min")
                    join_method_mode_f = METHOD_JOIN_MIN
                case ("join_max")
                    join_method_mode_f = METHOD_JOIN_MAX
                case ("join_median")
                    join_method_mode_f = METHOD_JOIN_MEDIAN
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: plateau_mode_f
            plateau_mode_f => c_char_as_view(plateau_mode)

            select case (plateau_mode_f)
                case ("plateau_ci_overlap")
                    plateau_mode_mode_f = MODE_PLATEAU_CI_OVERLAP
                case ("plateau_effect_size")
                    plateau_mode_mode_f = MODE_PLATEAU_EFFECT_SIZE
                case ("plateau_both")
                    plateau_mode_mode_f = MODE_PLATEAU_BOTH
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
                case ("bin_count_occupancy")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_OCCUPANCY
                case ("bin_count_stone_equal")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_EQUAL
                case ("bin_count_stone_size_weighted")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
                case ("bin_count_stone_pooled")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_parameter_search(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            plateau_established = plateau_established,&
            n_admissible_evaluated = n_admissible_evaluated,&
            trace_n_points = trace_n_points,&
            trace_n_neighbors = trace_n_neighbors,&
            trace_global_js_divergence = trace_global_js_divergence,&
            trace_ci_lower = trace_ci_lower,&
            trace_ci_upper = trace_ci_upper,&
            trace_ci_width = trace_ci_width,&
            trace_ci_width_relative = trace_ci_width_relative,&
            trace_ci_overlap = trace_ci_overlap,&
            trace_delta = trace_delta,&
            trace_delta_median = trace_delta_median,&
            trace_delta_max = trace_delta_max,&
            trace_selected_n_bins = trace_selected_n_bins,&
            trace_occupancy_failed = trace_occupancy_failed,&
            trace_n_pooled_residuals = trace_n_pooled_residuals,&
            trace_min_bin_occupancy = trace_min_bin_occupancy,&
            trace_mean_bin_occupancy = trace_mean_bin_occupancy,&
            trace_max_bin_occupancy = trace_max_bin_occupancy,&
            trace_sturges_bins = trace_sturges_bins,&
            trace_fd_bins = trace_fd_bins,&
            trace_shared_residual_range_low = trace_shared_residual_range_low,&
            trace_shared_residual_range_high = trace_shared_residual_range_high,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_parameter_search_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_fixed_k_search(module):run_js_comp_test_parameter_search_expert(subroutine)]]
    !| Walks the candidates of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !| study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !| ([[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !| residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !| `min_residuals_per_bin`
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !| ([[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]]),
    !| and both plateau criteria: CI overlap with the best candidate so far
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]])
    !| and relative effect-size stability
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
    !| Both are always computed and traced; `plateau_mode` decides which one stops the search. A
    !| stop by the effect-size criterion returns the candidate that triggered it.
    !|
    !| The returned candidate:
    !|
    !| 1. A plateau was found, or the grid has a single candidate and it is admissible: the best
    !| candidate, `plateau_established = .true.`.
    !| 2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
    !| median confidence-interval width across studies, `plateau_established = .false.`.
    !| 3. No candidate was admissible: the first (finest) candidate, with
    !| `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
    !| ranges `0.0`, and `plateau_established = .false.`.
    !|
    !| Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
    !| for its meaningful part. A failure to create the bootstrap's random number generator is
    !| returned in `ierr` and stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_parameter_search_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            n_bootstraps,&
            join_method,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds,&
            n_points,&
            n_neighbors,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            best_candidate_pair_confidence_interval,&
            plateau_established,&
            n_admissible_evaluated,&
            trace_n_points,&
            trace_n_neighbors,&
            trace_global_js_divergence,&
            trace_ci_lower,&
            trace_ci_upper,&
            trace_ci_width,&
            trace_ci_width_relative,&
            trace_ci_overlap,&
            trace_delta,&
            trace_delta_median,&
            trace_delta_max,&
            trace_selected_n_bins,&
            trace_occupancy_failed,&
            trace_n_pooled_residuals,&
            trace_min_bin_occupancy,&
            trace_mean_bin_occupancy,&
            trace_max_bin_occupancy,&
            trace_sturges_bins,&
            trace_fd_bins,&
            trace_shared_residual_range_low,&
            trace_shared_residual_range_high,&
            tmp_gene_means_perms,&
            tmp_gene_means_perm_all,&
            tmp_x_star,&
            tmp_neighborhood_indices_all_studies,&
            tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered,&
            tmp_point_pool,&
            tmp_point_draw,&
            tmp_n_bins_per_point,&
            tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high,&
            tmp_pmfs,&
            tmp_counts,&
            tmp_included_n_reps,&
            tmp_mean_pmf,&
            tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence,&
            tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_occupancy_failed,&
            tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy,&
            tmp_sturges_bins,&
            tmp_fd_bins,&
            tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high,&
            min_residuals_per_bin,&
            min_neighbor_overlap,&
            succeeding_ci_overlap,&
            plateau_mode,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            bin_count_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_parameter_search_expert_c")
        use tox_data_integration_js_comp_test_fixed_k_search, only: run_js_comp_test_parameter_search_expert
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_points_candidate
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(out), target :: n_points
            !! The finally chosen candidate's `n_points`
        integer(c_int), intent(out), target :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` when no candidate was admissible
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
            !! The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
            !! was admissible
        logical(c_bool), intent(out), target :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise, also when the smallest-uncertainty
            !! fallback returns a real candidate and interval
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array (an inadmissible candidate has no entry)
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate `n_points`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width (the values
            !! `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out), target :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            !! like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_x_star
            !! Working array: the current candidate's reference points
        integer(c_int), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out), target :: tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate
        integer(c_int), dimension(2, max_n_points_candidate), intent(out), target :: tmp_neighborhood_range
            !! Working array: one study's neighborhood spans, reused per study
        real(c_double), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        integer(c_int), dimension(256), intent(out), target :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(c_int), dimension(256), intent(out), target :: tmp_point_draw
            !! Working array of the bootstrap: one study's draw at one point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_bins_per_point
            !! Working array: the current candidate's bin count per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_low
            !! Working array: the current candidate's R_low per reference point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_high
            !! Working array: the current candidate's R_high per reference point
        real(c_double), dimension(max_n_points_candidate, 256, n_studies), intent(out), target :: tmp_pmfs
            !! Working array: every study's pmf for the current candidate
        integer(c_int), dimension(max_n_points_candidate, 256, n_studies), intent(out), target :: tmp_counts
            !! Working array: every study's histogram counts for the current candidate
        integer(c_int), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(c_double), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf
        integer(c_int), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, then bootstrap scratch
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_weights
            !! Working array: the current candidate's per-point weights, then bootstrap scratch
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD, then bootstrap scratch
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, bootstrapped in place
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array: the bootstrap's bottom-k/top-k heaps
        real(c_double), dimension(n_studies), intent(out), target :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_delta_perm
            !! Working array: sorting permutation for the per-study medians
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_best_uncertainty_confidence_interval
            !! Working array: confidence interval of the smallest-uncertainty admissible candidate
        real(c_double), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        logical(c_bool), dimension(max_n_points_candidate), intent(out), target :: tmp_occupancy_failed
            !! Working array: the current candidate's occupancy_failed flag per reference point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_pooled_residuals
            !! Working array: the current candidate's pooled residual count (N_j) per point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_min_bin_occupancy
            !! Working array: the current candidate's minimum bin occupancy per point
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_bin_occupancy
            !! Working array: the current candidate's mean bin occupancy per point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_max_bin_occupancy
            !! Working array: the current candidate's maximum bin occupancy per point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_sturges_bins
            !! Working array: the current candidate's Sturges' rule diagnostic per point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_fd_bins
            !! Working array: the current candidate's Freedman-Diaconis diagnostic per point
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_n_bins_per_point
            !! Working array: bin counts per point of the current best candidate
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: bin counts per point of the smallest-uncertainty admissible candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_low
            !! Working array: R_low per point of the current best candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_high
            !! Working array: R_high per point of the current best candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: R_low per point of the smallest-uncertainty admissible candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: R_high per point of the smallest-uncertainty admissible candidate
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
            !! one threshold for both, so a point the search marks as failed also fails the gate
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                | Value                                                                                        |
            !! |-------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_ci_overlap'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                                   |
            !! |----------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_occupancy'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f
        integer(int32) :: bin_count_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(plateau_established)
        M_CHECK_NON_NULL(n_admissible_evaluated)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(best_candidate_pair_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(trace_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_neighbors, 16)
        M_CHECK_ARRAY_NON_NULL(trace_global_js_divergence, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_lower, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_upper, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_width, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_width_relative, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_ci_overlap, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta, n_studies * 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_median, 16)
        M_CHECK_ARRAY_NON_NULL(trace_delta_max, 16)
        M_CHECK_ARRAY_NON_NULL(trace_selected_n_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_occupancy_failed, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_pooled_residuals, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_min_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_mean_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_max_bin_occupancy, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_sturges_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_fd_bins, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_shared_residual_range_low, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(trace_shared_residual_range_high, max_n_points_candidate * 16)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_x_star, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_indices_all_studies, max_n_neighbors_candidate * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_range, 2 * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * max_n_neighbors_candidate * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, max_n_points_candidate * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, max_n_points_candidate * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_included_n_reps, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_delta_perm, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_occupancy_failed, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_n_pooled_residuals, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_min_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_max_bin_occupancy, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_sturges_bins, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_fd_bins, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(plateau_mode, 19)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

        block
            character(len=:), pointer :: join_method_f
            join_method_f => c_char_as_view(join_method)

            select case (join_method_f)
                case ("join_min")
                    join_method_mode_f = METHOD_JOIN_MIN
                case ("join_max")
                    join_method_mode_f = METHOD_JOIN_MAX
                case ("join_median")
                    join_method_mode_f = METHOD_JOIN_MEDIAN
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: plateau_mode_f
            plateau_mode_f => c_char_as_view(plateau_mode)

            select case (plateau_mode_f)
                case ("plateau_ci_overlap")
                    plateau_mode_mode_f = MODE_PLATEAU_CI_OVERLAP
                case ("plateau_effect_size")
                    plateau_mode_mode_f = MODE_PLATEAU_EFFECT_SIZE
                case ("plateau_both")
                    plateau_mode_mode_f = MODE_PLATEAU_BOTH
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block
        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
                case ("bin_count_occupancy")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_OCCUPANCY
                case ("bin_count_stone_equal")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_EQUAL
                case ("bin_count_stone_size_weighted")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
                case ("bin_count_stone_pooled")
                    bin_count_mode_mode_f = MODE_BIN_COUNT_STONE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call run_js_comp_test_parameter_search_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            plateau_established = plateau_established,&
            n_admissible_evaluated = n_admissible_evaluated,&
            trace_n_points = trace_n_points,&
            trace_n_neighbors = trace_n_neighbors,&
            trace_global_js_divergence = trace_global_js_divergence,&
            trace_ci_lower = trace_ci_lower,&
            trace_ci_upper = trace_ci_upper,&
            trace_ci_width = trace_ci_width,&
            trace_ci_width_relative = trace_ci_width_relative,&
            trace_ci_overlap = trace_ci_overlap,&
            trace_delta = trace_delta,&
            trace_delta_median = trace_delta_median,&
            trace_delta_max = trace_delta_max,&
            trace_selected_n_bins = trace_selected_n_bins,&
            trace_occupancy_failed = trace_occupancy_failed,&
            trace_n_pooled_residuals = trace_n_pooled_residuals,&
            trace_min_bin_occupancy = trace_min_bin_occupancy,&
            trace_mean_bin_occupancy = trace_mean_bin_occupancy,&
            trace_max_bin_occupancy = trace_max_bin_occupancy,&
            trace_sturges_bins = trace_sturges_bins,&
            trace_fd_bins = trace_fd_bins,&
            trace_shared_residual_range_low = trace_shared_residual_range_low,&
            trace_shared_residual_range_high = trace_shared_residual_range_high,&
            tmp_gene_means_perms = tmp_gene_means_perms,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_neighborhood_indices_all_studies = tmp_neighborhood_indices_all_studies,&
            tmp_neighborhood_range = tmp_neighborhood_range,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            tmp_pmfs = tmp_pmfs,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps = tmp_included_n_reps,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_mean_pmf_included_n_reps = tmp_mean_pmf_included_n_reps,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_confidence_interval = tmp_confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_prev_global_js_divergence = tmp_prev_global_js_divergence,&
            tmp_delta_perm = tmp_delta_perm,&
            tmp_best_uncertainty_confidence_interval = tmp_best_uncertainty_confidence_interval,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_occupancy_failed = tmp_occupancy_failed,&
            tmp_n_pooled_residuals = tmp_n_pooled_residuals,&
            tmp_min_bin_occupancy = tmp_min_bin_occupancy,&
            tmp_mean_bin_occupancy = tmp_mean_bin_occupancy,&
            tmp_max_bin_occupancy = tmp_max_bin_occupancy,&
            tmp_sturges_bins = tmp_sturges_bins,&
            tmp_fd_bins = tmp_fd_bins,&
            tmp_best_n_bins_per_point = tmp_best_n_bins_per_point,&
            tmp_best_uncertainty_n_bins_per_point = tmp_best_uncertainty_n_bins_per_point,&
            tmp_best_shared_residual_range_low = tmp_best_shared_residual_range_low,&
            tmp_best_shared_residual_range_high = tmp_best_shared_residual_range_high,&
            tmp_best_uncertainty_shared_residual_range_low = tmp_best_uncertainty_shared_residual_range_low,&
            tmp_best_uncertainty_shared_residual_range_high = tmp_best_uncertainty_shared_residual_range_high,&
            min_residuals_per_bin = min_residuals_per_bin,&
            min_neighbor_overlap = min_neighbor_overlap,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_mode = plateau_mode_mode_f,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_parameter_search_expert_c

end module tox_data_integration_js_comp_test_fixed_k_search_c
#endif
