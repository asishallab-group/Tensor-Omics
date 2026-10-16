#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Adaptive Parameter Search
!|
!| The data-driven choice of the JS-Comp-Test's adaptive neighborhoods, each grown from `k_start`
!| pooled entries in rounds of `k_step` up to `k_max` while the residual dispersion stays stable.
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| walks the `(k_start, k_step, k_max)` candidates of
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
!| and stops at the first admissible candidate whose bootstrapped JSDs have reached a plateau, by
!| the criteria of
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]];
!| every candidate it tries is logged with an `ADAPTIVE_CANDIDATE_*` status.
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| builds the neighborhoods for the chosen knobs. `calc_adaptive_js_comp_test_bounds` sizes the
!| work arrays for a caller that allocates its own.
!|
!| Generated from [[tox_data_integration_js_comp_test_adaptive_search_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test_adaptive_search
    use f42_safeguard
    use tox_data_integration_js_comp_test_adaptive_search_impl, only: generate_adaptive_js_comp_test_candidates_impl, run_js_comp_test_adaptive_parameter_search_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use tox_data_integration_js_comp_test_binning_impl, only: MAX_N_BINS, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED
    use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
    use tox_data_integration_js_comp_test_search_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH
    use tox_data_integration_js_comp_test_search_impl, only: MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE
    use tox_data_integration_stats_impl, only: calc_js_comp_test_n_top_k_jsds
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_real
    use tox_errors, only: validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: generate_adaptive_js_comp_test_candidates
    public :: run_js_comp_test_adaptive_parameter_search
    public :: run_js_comp_test_adaptive_parameter_search_expert

contains

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl]].
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]].
    !| An adaptive candidate fixes only its growth knobs (each neighborhood grows from `k_start`
    !| pooled entries in rounds of `k_step` up to `k_max`); its number of reference points emerges.
    !| `k_start` shrinks from one candidate to the next, so the candidates **ascend** in the number
    !| of reference points they are expected to produce.
    !|
    !| With `N = max_n_genes_all_studies * n_studies` (the padded pool size) and
    !| `floor_k = max(10, 2*n_studies)`:
    !|
    !| - `k_start_1 = max(floor_k, ceiling(0.02 * N))`
    !| - `k_start_t = floor(k_start_1 * 0.8**(t-1))`
    !| - stop once `k_start_t < floor_k`, or once 16 candidates are accepted
    !| - `k_step_t = ceiling(0.25 * k_start_t)`
    !| - `k_max_t = ceiling(4 * k_start_t)`
    !|
    !| The constants are [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_START_FRACTION(variable)]] (0.02),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_GAMMA(variable)]] (0.8),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_STEP_FRACTION(variable)]] (0.25),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]] (4),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_START_MIN_ABS(variable)]] (10) and
    !| [[tox_data_integration_js_comp_test_search_impl(module):MAX_CANDIDATE_PAIRS(variable)]] (16).
    !| The values strictly descend without deduplication. At least one candidate is always
    !| returned, even for data smaller than `k_start_1`, whose construction then fails with a
    !| too-few-means status in the search.
    !|
    !| Not checked here: the values fit 32-bit integers only while
    !| `k_start_1 <= floor(huge(1_int32)/ADAPTIVE_K_MAX_FACTOR)` ([[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]]);
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]
    !| rejects inputs beyond that.
    pure subroutine generate_adaptive_js_comp_test_candidates(&
            max_n_genes_all_studies,&
            n_studies,&
            candidates_k_start_k_step_k_max,&
            n_candidates,&
            ierr&
        )
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(3, 16), intent(out) :: candidates_k_start_k_step_k_max
            !! Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
            !! The first `n_candidates` elements will hold the results.
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate triples actually filled (at least 1, at most 16)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=2_int32, min=1_int32)
        if (is_err(ierr)) return
#endif

        call generate_adaptive_js_comp_test_candidates_impl(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            candidates_k_start_k_step_k_max = candidates_k_start_k_step_k_max,&
            n_candidates = n_candidates&
        )
    end subroutine generate_adaptive_js_comp_test_candidates

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test_adaptive_search(module):run_js_comp_test_adaptive_parameter_search_expert]] to prepare it yourself.
    !| The adaptive-neighborhood counterpart of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !| Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
    !| in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
    !| The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.
    !|
    !| **Per candidate**, in order:
    !|
    !| 1. Grow its neighborhoods exactly as
    !| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    !| `min_study_neighbors`.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !| reference points emerged.
    !| 3. Gate 1:
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !| on the pooled ranges.
    !| 4. The pipeline of
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !| up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    !| (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    !| count
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    !| first plateau `plateau_mode` selects.
    !|
    !| **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
    !| `candidate_n_points` and `candidate_status`, until the search stops. The status is the first
    !| failing check, in the order above:
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED(variable)]]
    !| (too few means or residuals),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_EMPTY_STUDY(variable)]]
    !| (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
    !| some point),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED(variable)]],
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_OVERLAP_FAILED(variable)]],
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED(variable)]],
    !| or else [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_EVALUATED(variable)]]
    !| (admissible, with a trace column).
    !|
    !| **Result**, one of three cases, as in the fixed-k search:
    !|
    !| 1. A plateau was found, or the sequence had a single candidate and it passed both gates:
    !| that candidate's knobs, `plateau_established = .true.`.
    !| 2. No plateau but at least one admissible candidate: the admissible candidate with the
    !| smallest median confidence-interval width, `plateau_established = .false.`.
    !| 3. Nothing admissible, including a single candidate that was rejected: candidate 1's knobs,
    !| confidence interval `-1.0`, `plateau_established = .false.`.
    !|
    !| `n_points` is the selected candidate's logged point count; in case 3 the bin counts are
    !| `m_min` and the ranges `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
    !| `construct_adaptive_neighborhoods` with the returned knobs and the same construction
    !| options, then `run_js_comp_test_adaptive` with the same bin-count options
    !| (`bin_count_mode` included): its `global_js_divergence` equals the selected candidate's
    !| `trace_global_js_divergence` column.
    !|
    !| **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
    !| place of `trace_n_neighbors`; the per-point traces are jagged in the same way.
    !|
    !| **Memory.** The per-point arrays are sized by `max_n_points_candidate` (from
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]);
    !| the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
    !| `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
    !| `max_n_reps_all_studies * N` values each, the size of `residuals`.
    !|
    !| A failure to create the bootstrap's random number generator is returned in `ierr` and
    !| stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_adaptive_parameter_search(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            n_bootstraps,&
            join_method,&
            max_n_points_candidate,&
            k_start,&
            k_step,&
            k_max,&
            n_points,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            best_candidate_confidence_interval,&
            plateau_established,&
            n_admissible_evaluated,&
            trace_k_start,&
            trace_k_step,&
            trace_k_max,&
            trace_n_points,&
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
            n_candidates_tried,&
            candidate_k_start,&
            candidate_n_points,&
            candidate_status,&
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
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            bin_count_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_points_candidate
            !! Most reference points a candidate may emerge with; a candidate with more is rejected
            !! with the capacity status. The per-point outputs and traces below are sized by it, so
            !! the caller has to know it up front
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds]].
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(out) :: k_start
            !! The finally chosen candidate's `k_start`
        integer(int32), intent(out) :: k_step
            !! The finally chosen candidate's `k_step`
        integer(int32), intent(out) :: k_max
            !! The finally chosen candidate's `k_max`
        integer(int32), intent(out) :: n_points
            !! Number of reference points the finally chosen candidate's construction produced
        integer(int32), dimension(max_n_points_candidate), intent(out) :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful. Without an admissible candidate, the first
            !! `min(n_points, max_n_points_candidate)` entries are `m_min`. The rest are 0
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` without an admissible candidate; the
            !! rest are `0.0`
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_confidence_interval
            !! The bootstrapped JSD confidence interval of the finally chosen candidate;
            !! `-1.0_real64` throughout where no admissible candidate backs the result
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise
        integer(int32), intent(out) :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array
        integer(int32), dimension(16), intent(out) :: trace_k_start
            !! Per-admissible-candidate `k_start`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width; `-1.0` for the
            !! first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful (the rest are 0), so a Python/R
            !! caller slices `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins (`.false.`
            !! beyond)
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            !! as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), intent(out) :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(int32), dimension(16), intent(out) :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! The first `n_candidates_tried` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: candidate_n_points
            !! Number of reference points each tried candidate's construction produced (0 for too
            !! few means)
            !! The first `n_candidates_tried` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: candidate_status
            !! What happened to every tried candidate: `0_int32` evaluated,
            !! `1_int32` construction failed,
            !! `2_int32` a study neighborhood below `min_study_neighbors`
            !! or a study without a non-NaN residual at some point,
            !! `3_int32` capacity exceeded,
            !! `4_int32` overlap gate failed,
            !! `5_int32` consensus-pmf count gate failed
            !! The first `n_candidates_tried` elements will hold the results.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count every bin must reach, both in the bin-count search (the floor in every
            !! `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                | Value                                                                                        |
            !! |-------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :), allocatable :: tmp_pooled_neighborhood_range
        integer(int32), dimension(:, :), allocatable :: tmp_n_neighbors_per_point
        integer(int32), dimension(:), allocatable :: tmp_stop_reason
        real(real64), dimension(:), allocatable :: tmp_neighborhood_dispersion
        real(real64), dimension(:), allocatable :: tmp_neighborhood_mad
        integer(int32), dimension(:, :), allocatable :: tmp_point_neighborhood_indices
        real(real64), dimension(:, :), allocatable :: tmp_neighbor_residuals
        integer(int32), dimension(:), allocatable :: tmp_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_high
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        integer(int32), dimension(:, :, :), allocatable :: tmp_counts
        integer(int32), dimension(:, :), allocatable :: tmp_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:), allocatable :: tmp_mean_pmf_included_n_reps
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw
        real(real64), dimension(:, :), allocatable :: tmp_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:), allocatable :: tmp_prev_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_delta_perm
        real(real64), dimension(:, :), allocatable :: tmp_best_uncertainty_confidence_interval
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        logical(c_bool), dimension(:), allocatable :: tmp_occupancy_failed
        integer(int32), dimension(:), allocatable :: tmp_n_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_min_bin_occupancy
        real(real64), dimension(:), allocatable :: tmp_mean_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_max_bin_occupancy
        integer(int32), dimension(:), allocatable :: tmp_sturges_bins
        integer(int32), dimension(:), allocatable :: tmp_fd_bins
        integer(int32), dimension(:), allocatable :: tmp_best_n_bins_per_point
        integer(int32), dimension(:), allocatable :: tmp_best_uncertainty_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_shared_residual_range_high
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_best_uncertainty_shared_residual_range_high

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=46_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=47_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=48_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=50_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=51_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=52_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=53_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=54_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=55_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=56_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=57_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=58_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=59_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(tau, ierr, arg_pos=61_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=62_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=63_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=64_int32, min=1_int32)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=49_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=65_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_n_top_k_jsds(&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            n_top_k = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_x_star(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_stop_reason(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        M_ALLOCATE(tmp_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_counts(max_n_points_candidate, 256, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_counts(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_point_pool(256))
        M_ALLOCATE(tmp_point_draw(256))
        M_ALLOCATE(tmp_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_delta_perm(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_occupancy_failed(max_n_points_candidate))
        M_ALLOCATE(tmp_n_pooled_residuals(max_n_points_candidate))
        M_ALLOCATE(tmp_min_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_mean_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_max_bin_occupancy(max_n_points_candidate))
        M_ALLOCATE(tmp_sturges_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_fd_bins(max_n_points_candidate))
        M_ALLOCATE(tmp_best_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))

        call run_js_comp_test_adaptive_parameter_search_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            n_points = n_points,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            best_candidate_confidence_interval = best_candidate_confidence_interval,&
            plateau_established = plateau_established,&
            n_admissible_evaluated = n_admissible_evaluated,&
            trace_k_start = trace_k_start,&
            trace_k_step = trace_k_step,&
            trace_k_max = trace_k_max,&
            trace_n_points = trace_n_points,&
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
            n_candidates_tried = n_candidates_tried,&
            candidate_k_start = candidate_k_start,&
            candidate_n_points = candidate_n_points,&
            candidate_status = candidate_status,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_pooled_neighborhood_range = tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point = tmp_n_neighbors_per_point,&
            tmp_stop_reason = tmp_stop_reason,&
            tmp_neighborhood_dispersion = tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad = tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
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
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
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
            plateau_mode = plateau_mode,&
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
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_parameter_search

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test_adaptive_search(module):run_js_comp_test_adaptive_parameter_search]] does both.
    !| The adaptive-neighborhood counterpart of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !| Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
    !| in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
    !| The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.
    !|
    !| **Per candidate**, in order:
    !|
    !| 1. Grow its neighborhoods exactly as
    !| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    !| `min_study_neighbors`.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !| reference points emerged.
    !| 3. Gate 1:
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !| on the pooled ranges.
    !| 4. The pipeline of
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !| up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    !| (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    !| count
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    !| first plateau `plateau_mode` selects.
    !|
    !| **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
    !| `candidate_n_points` and `candidate_status`, until the search stops. The status is the first
    !| failing check, in the order above:
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED(variable)]]
    !| (too few means or residuals),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_EMPTY_STUDY(variable)]]
    !| (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
    !| some point),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED(variable)]],
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_OVERLAP_FAILED(variable)]],
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED(variable)]],
    !| or else [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_CANDIDATE_EVALUATED(variable)]]
    !| (admissible, with a trace column).
    !|
    !| **Result**, one of three cases, as in the fixed-k search:
    !|
    !| 1. A plateau was found, or the sequence had a single candidate and it passed both gates:
    !| that candidate's knobs, `plateau_established = .true.`.
    !| 2. No plateau but at least one admissible candidate: the admissible candidate with the
    !| smallest median confidence-interval width, `plateau_established = .false.`.
    !| 3. Nothing admissible, including a single candidate that was rejected: candidate 1's knobs,
    !| confidence interval `-1.0`, `plateau_established = .false.`.
    !|
    !| `n_points` is the selected candidate's logged point count; in case 3 the bin counts are
    !| `m_min` and the ranges `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
    !| `construct_adaptive_neighborhoods` with the returned knobs and the same construction
    !| options, then `run_js_comp_test_adaptive` with the same bin-count options
    !| (`bin_count_mode` included): its `global_js_divergence` equals the selected candidate's
    !| `trace_global_js_divergence` column.
    !|
    !| **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
    !| place of `trace_n_neighbors`; the per-point traces are jagged in the same way.
    !|
    !| **Memory.** The per-point arrays are sized by `max_n_points_candidate` (from
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]);
    !| the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
    !| `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
    !| `max_n_reps_all_studies * N` values each, the size of `residuals`.
    !|
    !| A failure to create the bootstrap's random number generator is returned in `ierr` and
    !| stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_adaptive_parameter_search_expert(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            n_bootstraps,&
            join_method,&
            max_n_points_candidate,&
            n_bootstrapping_top_k_jsds,&
            k_start,&
            k_step,&
            k_max,&
            n_points,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            best_candidate_confidence_interval,&
            plateau_established,&
            n_admissible_evaluated,&
            trace_k_start,&
            trace_k_step,&
            trace_k_max,&
            trace_n_points,&
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
            n_candidates_tried,&
            candidate_k_start,&
            candidate_n_points,&
            candidate_status,&
            tmp_gene_means_perm_all,&
            tmp_x_star,&
            tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point,&
            tmp_stop_reason,&
            tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals,&
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
            tmp_point_pool,&
            tmp_point_draw,&
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
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            bin_count_mode,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_points_candidate
            !! Most reference points a candidate may emerge with; a candidate with more is rejected
            !! with the capacity status. The per-point outputs and traces below are sized by it, so
            !! the caller has to know it up front
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                                  |
            !! |-----------------------------------------|----------------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(out) :: k_start
            !! The finally chosen candidate's `k_start`
        integer(int32), intent(out) :: k_step
            !! The finally chosen candidate's `k_step`
        integer(int32), intent(out) :: k_max
            !! The finally chosen candidate's `k_max`
        integer(int32), intent(out) :: n_points
            !! Number of reference points the finally chosen candidate's construction produced
        integer(int32), dimension(max_n_points_candidate), intent(out) :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful. Without an admissible candidate, the first
            !! `min(n_points, max_n_points_candidate)` entries are `m_min`. The rest are 0
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` without an admissible candidate; the
            !! rest are `0.0`
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_confidence_interval
            !! The bootstrapped JSD confidence interval of the finally chosen candidate;
            !! `-1.0_real64` throughout where no admissible candidate backs the result
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise
        integer(int32), intent(out) :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array
        integer(int32), dimension(16), intent(out) :: trace_k_start
            !! Per-admissible-candidate `k_start`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width; `-1.0` for the
            !! first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful (the rest are 0), so a Python/R
            !! caller slices `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins (`.false.`
            !! beyond)
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            !! as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), intent(out) :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(int32), dimension(16), intent(out) :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! The first `n_candidates_tried` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: candidate_n_points
            !! Number of reference points each tried candidate's construction produced (0 for too
            !! few means)
            !! The first `n_candidates_tried` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: candidate_status
            !! What happened to every tried candidate: `0_int32` evaluated,
            !! `1_int32` construction failed,
            !! `2_int32` a study neighborhood below `min_study_neighbors`
            !! or a study without a non-NaN residual at some point,
            !! `3_int32` capacity exceeded,
            !! `4_int32` overlap gate failed,
            !! `5_int32` consensus-pmf count gate failed
            !! The first `n_candidates_tried` elements will hold the results.
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_x_star
            !! Working array: the current candidate's reference points
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_neighborhood_range
            !! Working array: the current candidate's `[first, last]` pooled position per neighborhood
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: tmp_n_neighbors_per_point
            !! Working array: the current candidate's per-study entry count per neighborhood
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_stop_reason
            !! Working array: the current candidate's stop reason per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_dispersion
            !! Working array: the current candidate's dispersion per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_mad
            !! Working array: the current candidate's MAD of the pooled means per neighborhood
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array: the current candidate's per-point bin count
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array: the current candidate's per-point R_low
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array: the current candidate's per-point R_high
        real(real64), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's pmf for the current candidate
        integer(int32), dimension(max_n_points_candidate, 256, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's histogram counts for the current candidate
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, then bootstrap scratch
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-point weights, then bootstrap scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD, then bootstrap scratch
        integer(int32), dimension(256), intent(out) :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_point_draw
            !! Working array of the bootstrap: one study's draw at one point
        real(real64), dimension(2, n_studies), intent(out) :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, bootstrapped in place
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array: sorting permutation for the per-study medians
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array: confidence interval of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array: the current candidate's per-point occupancy_failed flag
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array: the current candidate's per-point pooled residual count
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array: the current candidate's per-point minimum bin occupancy
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array: the current candidate's per-point mean bin occupancy
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array: the current candidate's per-point maximum bin occupancy
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array: the current candidate's per-point Sturges' rule diagnostic
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array: the current candidate's per-point Freedman-Diaconis rule diagnostic
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array: snapshot of tmp_n_bins_per_point for the current best candidate
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: snapshot of tmp_n_bins_per_point for the smallest-uncertainty candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array: snapshot of tmp_shared_residual_range_low for the current best candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array: snapshot of tmp_shared_residual_range_high for the current best candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: R_low per point of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: R_high per point of the smallest-uncertainty admissible candidate
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count every bin must reach, both in the bin-count search (the floor in every
            !! `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                | Value                                                                                        |
            !! |-------------------------------------|----------------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=9_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=94_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=95_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=96_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=98_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=99_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=100_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=101_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=102_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=103_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=104_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=105_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=106_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=107_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(tau, ierr, arg_pos=109_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=110_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=111_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=112_int32, min=1_int32)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=97_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=113_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_adaptive_parameter_search_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            n_points = n_points,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            best_candidate_confidence_interval = best_candidate_confidence_interval,&
            plateau_established = plateau_established,&
            n_admissible_evaluated = n_admissible_evaluated,&
            trace_k_start = trace_k_start,&
            trace_k_step = trace_k_step,&
            trace_k_max = trace_k_max,&
            trace_n_points = trace_n_points,&
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
            n_candidates_tried = n_candidates_tried,&
            candidate_k_start = candidate_k_start,&
            candidate_n_points = candidate_n_points,&
            candidate_status = candidate_status,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_x_star = tmp_x_star,&
            tmp_pooled_neighborhood_range = tmp_pooled_neighborhood_range,&
            tmp_n_neighbors_per_point = tmp_n_neighbors_per_point,&
            tmp_stop_reason = tmp_stop_reason,&
            tmp_neighborhood_dispersion = tmp_neighborhood_dispersion,&
            tmp_neighborhood_mad = tmp_neighborhood_mad,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
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
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
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
            plateau_mode = plateau_mode,&
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
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_parameter_search_expert

end module tox_data_integration_js_comp_test_adaptive_search
