#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

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
module tox_data_integration_js_comp_test_adaptive_search_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err_once, is_err, validate_in_range_int, ERR_INVALID_INPUT
    use tox_data_integration_preprocessing_impl, only: grow_adaptive_neighborhoods, ADAPTIVE_STATUS_OK, &
                                                       ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD, sort_pooled_gene_means
    use tox_data_integration_js_comp_test_binning_impl, only: determine_adaptive_point_bin_counts
    use tox_data_integration_js_comp_test_histograms_impl, only: build_adaptive_histograms
    use tox_data_integration_js_comp_test_search_impl, only: check_neighborhood_overlaps_impl, score_candidate_histograms, &
                                                             evaluate_admissible_candidate, select_final_candidate, &
                                                             MAX_CANDIDATE_PAIRS, NOT_COMPUTED
    M_IMPLICIT_NONE
    private
    public :: calc_adaptive_js_comp_test_bounds, generate_adaptive_js_comp_test_candidates_impl, &
             run_js_comp_test_adaptive_parameter_search_impl, &
             ADAPTIVE_START_FRACTION, ADAPTIVE_GAMMA, ADAPTIVE_K_STEP_FRACTION, ADAPTIVE_K_MAX_FACTOR, &
             ADAPTIVE_K_START_MIN_ABS, ADAPTIVE_POINT_CAPACITY_FACTOR, &
             ADAPTIVE_CANDIDATE_EVALUATED, ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED, ADAPTIVE_CANDIDATE_EMPTY_STUDY, &
             ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED, ADAPTIVE_CANDIDATE_OVERLAP_FAILED, &
             ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED

    ! Tuning constants of the adaptive candidate sequence (generate_adaptive_js_comp_test_candidates_impl)
    ! and its capacity (calc_adaptive_js_comp_test_bounds), relative to the padded pool size
    ! `N = max_n_genes_all_studies * n_studies`. Provisional, pending real-data validation.
    real(real64), parameter :: ADAPTIVE_START_FRACTION = 0.02_real64
        !! First (largest) candidate `k_start`, as a fraction of the padded pool size `N`
    real(real64), parameter :: ADAPTIVE_GAMMA = 0.8_real64
        !! Decay factor of `k_start` between successive adaptive candidates
    real(real64), parameter :: ADAPTIVE_K_STEP_FRACTION = 0.25_real64
        !! A candidate's `k_step`, as a fraction of its own `k_start`
    real(real64), parameter :: ADAPTIVE_K_MAX_FACTOR = 4.0_real64
        !! A candidate's `k_max`, as a multiple of its own `k_start`
    integer(int32), parameter :: ADAPTIVE_K_START_MIN_ABS = 10_int32
        !! Absolute lower bound on any candidate's `k_start` (raised to `2*n_studies` when that
        !! is larger)
    real(real64), parameter :: ADAPTIVE_POINT_CAPACITY_FACTOR = 2.0_real64
        !! Safety factor of the practical reference-point capacity over the naive estimate
        !! `N / k_start`, in calc_adaptive_js_comp_test_bounds

    ! Per-candidate status codes of run_js_comp_test_adaptive_parameter_search_impl's log
    ! (output codes, CM_ macros first like construct_adaptive_neighborhoods_impl's).
#define CM_ADAPTIVE_CANDIDATE_EVALUATED 0_int32
#define CM_ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED 1_int32
#define CM_ADAPTIVE_CANDIDATE_EMPTY_STUDY 2_int32
#define CM_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED 3_int32
#define CM_ADAPTIVE_CANDIDATE_OVERLAP_FAILED 4_int32
#define CM_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED 5_int32
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_EVALUATED = CM_ADAPTIVE_CANDIDATE_EVALUATED
        !! Candidate status: passed both admissibility gates and was bootstrapped
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED = CM_ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED
        !! Candidate status: the construction reported too few means or too few residuals
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_EMPTY_STUDY = CM_ADAPTIVE_CANDIDATE_EMPTY_STUDY
        !! Candidate status: some neighborhood has fewer than `min_study_neighbors` entries of some
        !! study, or some study has no non-NaN residual in some neighborhood
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED = CM_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED
        !! Candidate status: more reference points emerged than `max_n_points_candidate`
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_OVERLAP_FAILED = CM_ADAPTIVE_CANDIDATE_OVERLAP_FAILED
        !! Candidate status: two consecutive neighborhoods overlap by less than `min_neighbor_overlap`
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED = CM_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED
        !! Candidate status: some bin of the consensus pmf has fewer than `min_residuals_per_bin` residuals

contains

    !> M_EXPORT_C
    !| summary: Recommend a reference-point capacity for the adaptive js-comp-test candidate sequence
    !| AUTHOR_LASZLO_LANG
    !| Sizes the per-point work arrays of the adaptive parameter search over the candidates of
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| A candidate's reference-point count is only known after its neighborhoods are grown, so this
    !| is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
    !| `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
    !| candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.
    !|
    !| Rejects, with invalid input naming neither argument when the bound is joint, a
    !| `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
    !| first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
    !| `floor(huge(1_int32)/ADAPTIVE_K_MAX_FACTOR)` ([[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]]), so that its
    !| `k_max = ceiling(ADAPTIVE_K_MAX_FACTOR*k_start_1)` fits; the output is then left undefined.
    pure subroutine calc_adaptive_js_comp_test_bounds(max_n_genes_all_studies, n_studies, max_n_points_candidate, ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Practical upper bound on the number of reference points any adaptive candidate is
            !! allowed to emerge with
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: k_starts(MAX_CANDIDATE_PAIRS), n_k_starts
        integer(int64) :: n_pool, capacity

        call set_ok(ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=2_int32, min=1_int32)
        if (is_err(ierr)) return
        n_pool = int(max_n_genes_all_studies, int64)*int(n_studies, int64)
        ! Joint bounds on (G, S), so no arg_pos: the pool size must fit a 32-bit integer (array
        ! extents derive from it), and so must the first candidate's
        ! k_max = ceiling(ADAPTIVE_K_MAX_FACTOR*k_start_1).
        if (n_pool > int(huge(1_int32), int64)) call set_err_once(ierr, ERR_INVALID_INPUT)
        if (adaptive_k_start_1(max_n_genes_all_studies, n_studies) > &
            floor(real(huge(1_int32), real64)/ADAPTIVE_K_MAX_FACTOR, kind=int64)) &
            call set_err_once(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        call compute_adaptive_k_starts(max_n_genes_all_studies, n_studies, k_starts, n_k_starts)

        capacity = ceiling(ADAPTIVE_POINT_CAPACITY_FACTOR*real(n_pool, real64)/ &
                           real(k_starts(n_k_starts), real64), kind=int64) + 1_int64
        max_n_points_candidate = int(min(n_pool, capacity), kind=int32)
    end subroutine calc_adaptive_js_comp_test_bounds

    !> summary: Generate the ascending adaptive (k_start, k_step, k_max) candidate sequence
    !| AUTHOR_LASZLO_LANG
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
    !| - stop once `k_start_t < floor_k`, or once M_MAX_CANDIDATE_PAIRS candidates are accepted
    !| - `k_step_t = ceiling(0.25 * k_start_t)`
    !| - `k_max_t = ceiling(4 * k_start_t)`
    !|
    !| The constants are [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_START_FRACTION(variable)]] (0.02),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_GAMMA(variable)]] (0.8),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_STEP_FRACTION(variable)]] (0.25),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]] (4),
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_START_MIN_ABS(variable)]] (10) and
    !| [[tox_data_integration_js_comp_test_search_impl(module):MAX_CANDIDATE_PAIRS(variable)]] (M_MAX_CANDIDATE_PAIRS).
    !| The values strictly descend without deduplication. At least one candidate is always
    !| returned, even for data smaller than `k_start_1`, whose construction then fails with a
    !| too-few-means status in the search.
    !|
    !| Not checked here: the values fit 32-bit integers only while
    !| `k_start_1 <= floor(huge(1_int32)/ADAPTIVE_K_MAX_FACTOR)` ([[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]]);
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]
    !| rejects inputs beyond that.
    pure subroutine generate_adaptive_js_comp_test_candidates_impl(max_n_genes_all_studies, n_studies, &
                                                                   candidates_k_start_k_step_k_max, n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: candidates_k_start_k_step_k_max(3, M_MAX_CANDIDATE_PAIRS)
            !! Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate triples actually filled (at least 1, at most M_MAX_CANDIDATE_PAIRS)

        integer(int32) :: k_starts(MAX_CANDIDATE_PAIRS), i_candidate, k_start

        call compute_adaptive_k_starts(max_n_genes_all_studies, n_studies, k_starts, n_candidates)

        do concurrent(i_candidate=1:n_candidates) local(k_start) shared(k_starts, candidates_k_start_k_step_k_max)
            k_start = k_starts(i_candidate)
            candidates_k_start_k_step_k_max(1, i_candidate) = k_start
            candidates_k_start_k_step_k_max(2, i_candidate) = &
                ceiling(ADAPTIVE_K_STEP_FRACTION*real(k_start, real64), kind=int32)
            candidates_k_start_k_step_k_max(3, i_candidate) = &
                ceiling(ADAPTIVE_K_MAX_FACTOR*real(k_start, real64), kind=int32)
        end do
    end subroutine generate_adaptive_js_comp_test_candidates_impl

    !> summary: The adaptive candidate sequence's k_start values -- its one definition
    !| AUTHOR_LASZLO_LANG
    !| Shared by
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
    !| and `calc_adaptive_js_comp_test_bounds`, so both see the same smallest `k_start`. No
    !| deduplication is needed while `(1 - ADAPTIVE_GAMMA) * ADAPTIVE_K_START_MIN_ABS >= 1`: one
    !| decay step then lowers any `k_start >= floor_k` by at least 1
    !| (`test_adaptive_constants_preclude_dead_clamps` guards this). Each value is computed from
    !| `k_start_1` directly, so rounding does not accumulate.
    pure subroutine compute_adaptive_k_starts(max_n_genes_all_studies, n_studies, k_starts, n_k_starts)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies, at least 1
        integer(int32), intent(in) :: n_studies
            !! Number of studies, at least 1
        integer(int32), intent(out) :: k_starts(MAX_CANDIDATE_PAIRS)
            !! Accepted `k_start` values, strictly descending; only the first `n_k_starts` are set
        integer(int32), intent(out) :: n_k_starts
            !! Number of accepted `k_start` values, between 1 and MAX_CANDIDATE_PAIRS

        integer(int64) :: floor_k, k_start_1, k_start

        floor_k = adaptive_floor_k(n_studies)
        k_start_1 = adaptive_k_start_1(max_n_genes_all_studies, n_studies)

        ! k_start_1 >= floor_k, so the first value is always accepted.
        n_k_starts = 0_int32
        do while (n_k_starts < MAX_CANDIDATE_PAIRS)
            k_start = floor(real(k_start_1, real64)*ADAPTIVE_GAMMA**n_k_starts, kind=int64)
            if (k_start < floor_k) exit
            n_k_starts = n_k_starts + 1_int32
            k_starts(n_k_starts) = int(k_start, int32)
        end do
    end subroutine compute_adaptive_k_starts

    !> summary: The adaptive candidate sequence's lower bound on `k_start`, `max(10, 2*n_studies)`
    !| AUTHOR_LASZLO_LANG
    pure integer(int64) function adaptive_floor_k(n_studies) result(floor_k)
        integer(int32), intent(in) :: n_studies
            !! Number of studies, at least 1

        floor_k = max(int(ADAPTIVE_K_START_MIN_ABS, int64), 2_int64*int(n_studies, int64))
    end function adaptive_floor_k

    !> summary: The adaptive candidate sequence's first `k_start`, in 64-bit so it cannot overflow
    !| AUTHOR_LASZLO_LANG
    !| `max(floor_k, ceiling(0.02 * G * S))`; shared by the sequence and by the representability
    !| check in `calc_adaptive_js_comp_test_bounds`.
    pure integer(int64) function adaptive_k_start_1(max_n_genes_all_studies, n_studies) result(k_start_1)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies, at least 1
        integer(int32), intent(in) :: n_studies
            !! Number of studies, at least 1

        integer(int64) :: n_pool

        n_pool = int(max_n_genes_all_studies, int64)*int(n_studies, int64)
        k_start_1 = max(adaptive_floor_k(n_studies), ceiling(ADAPTIVE_START_FRACTION*real(n_pool, real64), kind=int64))
    end function adaptive_k_start_1

    !> summary: Search the ascending adaptive (k_start, k_step, k_max) candidate sequence for a stable JSD parameter setting (Issue #217)
    !| AUTHOR_LASZLO_LANG
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
    !|    [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !|    does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    !|    `min_study_neighbors`.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !|    reference points emerged.
    !| 3. Gate 1:
    !|    [[tox_data_integration_js_comp_test_search_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !|    on the pooled ranges.
    !| 4. The pipeline of
    !|    [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !|    up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    !|    (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    !|    count
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    !|    first plateau `plateau_mode` selects.
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
    !|    that candidate's knobs, `plateau_established = .true.`.
    !| 2. No plateau but at least one admissible candidate: the admissible candidate with the
    !|    smallest median confidence-interval width, `plateau_established = .false.`.
    !| 3. Nothing admissible, including a single candidate that was rejected: candidate 1's knobs,
    !|    confidence interval `-1.0`, `plateau_established = .false.`.
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
    subroutine run_js_comp_test_adaptive_parameter_search_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                                gene_means, residuals, n_bootstraps, join_method, &
                                                                max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                                                k_start, k_step, k_max, n_points, n_bins_per_point, &
                                                                shared_residual_range_low, shared_residual_range_high, &
                                                                best_candidate_confidence_interval, plateau_established, &
                                                                n_admissible_evaluated, trace_k_start, trace_k_step, &
                                                                trace_k_max, trace_n_points, trace_global_js_divergence, &
                                                                trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                                                trace_ci_width_relative, trace_ci_overlap, trace_delta, &
                                                                trace_delta_median, trace_delta_max, trace_selected_n_bins, &
                                                                trace_occupancy_failed, &
                                                                trace_n_pooled_residuals, trace_min_bin_occupancy, &
                                                                trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                                                trace_sturges_bins, trace_fd_bins, &
                                                                trace_shared_residual_range_low, &
                                                                trace_shared_residual_range_high, n_candidates_tried, &
                                                                candidate_k_start, candidate_n_points, candidate_status, &
                                                                tmp_gene_means_perm_all, tmp_x_star, &
                                                                tmp_pooled_neighborhood_range, tmp_n_neighbors_per_point, &
                                                                tmp_stop_reason, tmp_neighborhood_dispersion, &
                                                                tmp_neighborhood_mad, tmp_point_neighborhood_indices, &
                                                                tmp_neighbor_residuals, tmp_n_bins_per_point, &
                                                                tmp_shared_residual_range_low, tmp_shared_residual_range_high, &
                                                                tmp_pmfs, tmp_counts, tmp_included_n_reps, tmp_mean_pmf, &
                                                                tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                                                tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                                                tmp_point_pool, tmp_point_draw, &
                                                                tmp_confidence_interval, tmp_bootstrapping_top_k_jsds, &
                                                                tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                                tmp_best_uncertainty_confidence_interval, &
                                                                tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                                tmp_bin_counts_search, tmp_residual_to_study_perm, &
                                                                tmp_study_block_ends, tmp_study_block_sizes, tmp_occupancy_failed, &
                                                                tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                                tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                                                tmp_sturges_bins, tmp_fd_bins, tmp_best_n_bins_per_point, &
                                                                tmp_best_uncertainty_n_bins_per_point, &
                                                                tmp_best_shared_residual_range_low, &
                                                                tmp_best_shared_residual_range_high, &
                                                                tmp_best_uncertainty_shared_residual_range_low, &
                                                                tmp_best_uncertainty_shared_residual_range_high, &
                                                                min_residuals_per_bin, min_neighbor_overlap, &
                                                                succeeding_ci_overlap, plateau_mode, delta_median_threshold, &
                                                                delta_max_threshold, delta_epsilon, &
                                                                delta_min_consecutive_transitions, m_min, m_max, &
                                                                gamma_occupancy, lower_residual_range_quantile, &
                                                                upper_residual_range_quantile, &
                                                                two_sided_bootstrapping_significance_level, random_seed, &
                                                                tau, mad_distance_factor, max_pooled_residuals, &
                                                                min_study_neighbors, bin_count_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_points_candidate
            !! Most reference points a candidate may emerge with; a candidate with more is rejected
            !! with the capacity status. The per-point outputs and traces below are sized by it, so
            !! the caller has to know it up front
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_bounds, tox_data_integration_js_comp_test_adaptive_search_impl, JUST_INFO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_stats_impl, AUTO)
            !! DM_MIN(1_int32)
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
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_k_start
            !! Per-admissible-candidate `k_start`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width; `-1.0` for the
            !! first admissible candidate
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful (the rest are 0), so a Python/R
            !! caller slices `[:trace_n_points[t], t]`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins (`.false.`
            !! beyond)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            !! as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), intent(out) :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: candidate_n_points
            !! Number of reference points each tried candidate's construction produced (0 for too
            !! few means)
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: candidate_status
            !! What happened to every tried candidate: `CM_ADAPTIVE_CANDIDATE_EVALUATED` evaluated,
            !! `CM_ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED` construction failed,
            !! `CM_ADAPTIVE_CANDIDATE_EMPTY_STUDY` a study neighborhood below `min_study_neighbors`
            !! or a study without a non-NaN residual at some point,
            !! `CM_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED` capacity exceeded,
            !! `CM_ADAPTIVE_CANDIDATE_OVERLAP_FAILED` overlap gate failed,
            !! `CM_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED` consensus-pmf count gate failed
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
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
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's pmf for the current candidate
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's histogram counts for the current candidate
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, then bootstrap scratch
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-point weights, then bootstrap scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD, then bootstrap scratch
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
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
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
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
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_CI_OVERLAP)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(M_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! M_LOWER_RESIDUAL_RANGE_QUANTILE_ORDER_DOC
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_OCCUPANCY)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created

        integer(int32) :: candidates_k_start_k_step_k_max(3, MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, max_n_bins
        integer(int32) :: best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: actual_min_residuals_per_bin, selected_candidate_index
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index, actual_m_min
        integer(int32) :: actual_max_pooled_residuals, actual_min_study_neighbors, construction_status, max_n_neighbors
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value, actual_tau, actual_mad_distance_factor
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found, has_empty_study

        call set_ok(ierr)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        M_DEFAULT_VAL(plateau_mode, actual_plateau_mode, M_MODE_PLATEAU_CI_OVERLAP)
        M_DEFAULT_VAL(delta_median_threshold, actual_delta_median_threshold, M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_max_threshold, actual_delta_max_threshold, M_DELTA_MAX_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_epsilon, actual_delta_epsilon, M_DELTA_EPSILON_DEFAULT)
        M_DEFAULT_VAL(delta_min_consecutive_transitions, actual_delta_min_consecutive_transitions, M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        M_DEFAULT_VAL(m_min, actual_m_min, M_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(tau, actual_tau, M_ADAPTIVE_TAU_DEFAULT)
        M_DEFAULT_VAL(mad_distance_factor, actual_mad_distance_factor, M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        M_DEFAULT_VAL(max_pooled_residuals, actual_max_pooled_residuals, M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        M_DEFAULT_VAL(min_study_neighbors, actual_min_study_neighbors, M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)

        call generate_adaptive_js_comp_test_candidates_impl(max_n_genes_all_studies, n_studies, &
                                                            candidates_k_start_k_step_k_max, n_candidates)

        ! The construction's own sort (construct_adaptive_neighborhoods_impl), done once: every
        ! candidate grows on it, and a pooled position names the same entry as there.
        call sort_pooled_gene_means(gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perm_all)

        best_candidate_confidence_interval = NOT_COMPUTED
        best_candidate_index = 1_int32
        best_exceeded_ci_overlap_count = 0_int32
        best_uncertainty_value = huge(1.0_real64)
        best_uncertainty_candidate_index = 1_int32
        plateau_found = logical(.false., kind=c_bool)
        n_admissible_evaluated = 0_int32
        n_candidates_tried = 0_int32
        n_consecutive_effect_size_ok = 0_int32

        ! Test candidates, fewest emerging points first, and stop at the first JSD plateau.
        do i_candidate = 1, n_candidates
            k_start = candidates_k_start_k_step_k_max(1, i_candidate)
            k_step = candidates_k_start_k_step_k_max(2, i_candidate)
            k_max = candidates_k_start_k_step_k_max(3, i_candidate)

            call grow_adaptive_neighborhoods(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                             residuals, tmp_gene_means_perm_all, k_start, k_step, k_max, actual_tau, &
                                             actual_mad_distance_factor, actual_max_pooled_residuals, &
                                             actual_min_study_neighbors, n_points, tmp_x_star, &
                                             tmp_pooled_neighborhood_range, tmp_n_neighbors_per_point, tmp_stop_reason, &
                                             tmp_neighborhood_dispersion, tmp_neighborhood_mad, max_n_neighbors, &
                                             construction_status)

            n_candidates_tried = i_candidate
            candidate_k_start(i_candidate) = k_start
            candidate_n_points(i_candidate) = n_points
            candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_EVALUATED

            ! Construction status first: an empty study neighborhood leaves that study's pmf
            ! empty, and too few means/residuals means the growth did not run as specified.
            if (construction_status == ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_EMPTY_STUDY
                cycle
            else if (construction_status /= ADAPTIVE_STATUS_OK) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED
                cycle
            end if
            if (n_points > max_n_points_candidate) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED
                cycle
            end if

            ! First admissibility gate, on the pooled ranges the construction returns.
            call check_neighborhood_overlaps_impl(tmp_pooled_neighborhood_range(1:2, 1:n_points), n_points, &
                                                   actual_min_neighbor_overlap, all_have_min_neighbor_overlap)
            if (.not. all_have_min_neighbor_overlap) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_OVERLAP_FAILED
                cycle
            end if

            ! ===== PASS B (per point, sequential: one pooling buffer serves every point), the same
            ! helper as run_js_comp_test_adaptive_impl.
            call determine_adaptive_point_bin_counts(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                     n_studies, n_points, tmp_gene_means_perm_all, &
                                                     tmp_pooled_neighborhood_range, tmp_n_neighbors_per_point, &
                                                     tmp_n_bins_per_point, tmp_occupancy_failed, &
                                                     tmp_shared_residual_range_low, tmp_shared_residual_range_high, &
                                                     tmp_n_pooled_residuals, tmp_min_bin_occupancy, tmp_mean_bin_occupancy, &
                                                     tmp_max_bin_occupancy, tmp_sturges_bins, tmp_fd_bins, &
                                                     tmp_point_neighborhood_indices, tmp_pooled_residuals, &
                                                     tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                                     tmp_residual_to_study_perm, tmp_study_block_ends, tmp_study_block_sizes, &
                                                     m_min=m_min, m_max=m_max, &
                                                     min_residuals_per_bin=actual_min_residuals_per_bin, &
                                                     gamma_occupancy=gamma_occupancy, &
                                                     lower_residual_range_quantile=lower_residual_range_quantile, &
                                                     upper_residual_range_quantile=upper_residual_range_quantile, &
                                                     bin_count_mode=bin_count_mode)

            ! ===== PASS C (per point, per study, sequential: one gather buffer serves every pair).
            ! From here on the work arrays sized by an upper bound (max_n_points_candidate,
            ! M_MAX_N_BINS) are passed whole to explicit-shape dummies declared with the
            ! candidate's own sizes, so every routine sees their leading elements packed
            ! (sequence association) and nothing is copied.
            max_n_bins = maxval(tmp_n_bins_per_point(1:n_points))
            call build_adaptive_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                           n_points, max_n_bins, max_n_bins, tmp_gene_means_perm_all, &
                                           tmp_pooled_neighborhood_range, tmp_n_neighbors_per_point, &
                                           tmp_shared_residual_range_low, tmp_shared_residual_range_high, &
                                           tmp_n_bins_per_point, tmp_counts, tmp_pmfs, tmp_included_n_reps, &
                                           tmp_point_neighborhood_indices, tmp_neighbor_residuals)

            ! Second admissibility gate. A study whose genes at a point have only NaN residuals passes
            ! the construction's gene count but not this; it is logged as an empty-study candidate.
            call score_candidate_histograms(n_points, max_n_bins, n_studies, tmp_n_bins_per_point, &
                                            actual_min_residuals_per_bin, tmp_pmfs, tmp_counts, tmp_included_n_reps, &
                                            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                            tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                            tmp_confidence_interval, has_empty_study, all_bins_have_min_count)
            if (has_empty_study) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_EMPTY_STUDY
                cycle
            end if
            if (.not. all_bins_have_min_count) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED
                cycle
            end if

            n_admissible_evaluated = n_admissible_evaluated + 1_int32
            trace_k_start(n_admissible_evaluated) = k_start
            trace_k_step(n_admissible_evaluated) = k_step
            trace_k_max(n_admissible_evaluated) = k_max
            trace_n_points(n_admissible_evaluated) = n_points
            call evaluate_admissible_candidate(n_points, max_n_bins, n_studies, max_n_points_candidate, i_candidate, &
                                               n_admissible_evaluated, n_bootstraps, n_bootstrapping_top_k_jsds, &
                                               join_method, actual_plateau_mode, actual_succeeding_ci_overlap, &
                                               actual_delta_median_threshold, actual_delta_max_threshold, &
                                               actual_delta_epsilon, actual_delta_min_consecutive_transitions, &
                                               tmp_n_bins_per_point, tmp_occupancy_failed, tmp_n_pooled_residuals, &
                                               tmp_min_bin_occupancy, tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                               tmp_sturges_bins, tmp_fd_bins, tmp_shared_residual_range_low, &
                                               tmp_shared_residual_range_high, tmp_included_n_reps, tmp_mean_pmf_counts, &
                                               tmp_mean_pmf_included_n_reps, tmp_counts, tmp_pmfs, tmp_mean_pmf, &
                                               tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                               tmp_confidence_interval, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_ci_overlap, &
                                               trace_delta, trace_delta_median, trace_delta_max, trace_selected_n_bins, &
                                               trace_occupancy_failed, trace_n_pooled_residuals, trace_min_bin_occupancy, &
                                               trace_mean_bin_occupancy, trace_max_bin_occupancy, trace_sturges_bins, &
                                               trace_fd_bins, trace_shared_residual_range_low, &
                                               trace_shared_residual_range_high, tmp_prev_global_js_divergence, &
                                               n_consecutive_effect_size_ok, best_candidate_confidence_interval, &
                                               best_candidate_index, best_exceeded_ci_overlap_count, &
                                               tmp_best_n_bins_per_point, tmp_best_shared_residual_range_low, &
                                               tmp_best_shared_residual_range_high, best_uncertainty_value, &
                                               best_uncertainty_candidate_index, tmp_best_uncertainty_confidence_interval, &
                                               tmp_best_uncertainty_n_bins_per_point, &
                                               tmp_best_uncertainty_shared_residual_range_low, &
                                               tmp_best_uncertainty_shared_residual_range_high, plateau_found, &
                                               tmp_bootstrapping_top_k_jsds, tmp_point_pool, tmp_point_draw, &
                                               tmp_delta_perm, ierr, two_sided_bootstrapping_significance_level, &
                                               random_seed)
            if (is_err(ierr)) return
            if (plateau_found) exit
        end do

        ! The point count is the chosen candidate's logged one; the knobs come from the sequence. A
        ! single candidate counts as established only when it was admissible; a rejected one falls
        ! through to case 3, which returns that same candidate.
        call select_final_candidate(max_n_points_candidate, n_studies, n_candidates, n_admissible_evaluated, plateau_found, &
                                    actual_m_min, best_candidate_index, candidate_n_points(best_candidate_index), &
                                    tmp_best_n_bins_per_point, tmp_best_shared_residual_range_low, &
                                    tmp_best_shared_residual_range_high, best_uncertainty_candidate_index, &
                                    candidate_n_points(best_uncertainty_candidate_index), &
                                    tmp_best_uncertainty_confidence_interval, tmp_best_uncertainty_n_bins_per_point, &
                                    tmp_best_uncertainty_shared_residual_range_low, &
                                    tmp_best_uncertainty_shared_residual_range_high, candidate_n_points(1), &
                                    selected_candidate_index, n_points, n_bins_per_point, shared_residual_range_low, &
                                    shared_residual_range_high, best_candidate_confidence_interval, plateau_established)
        k_start = candidates_k_start_k_step_k_max(1, selected_candidate_index)
        k_step = candidates_k_start_k_step_k_max(2, selected_candidate_index)
        k_max = candidates_k_start_k_step_k_max(3, selected_candidate_index)
    end subroutine run_js_comp_test_adaptive_parameter_search_impl

end module tox_data_integration_js_comp_test_adaptive_search_impl
