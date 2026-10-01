#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
!| before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
!| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| generates candidate `(n_points, n_neighbors)` pairs only -- each candidate's real
!| per-neighborhood histogram bin count is decided later, per reference point, by
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]),
!| two admissibility gates a candidate must pass before it is bootstrapped
!| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]),
!| and the plateau check that decides when the search has converged
!| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
!| is a brute-force reference implementation of the same per-point bin-count search, exhaustively
!| testing every candidate `M` instead of the fast geometric-search-then-refinement the production
!| routine uses, for validating that the fast search's own result is correct.
!| `calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
!| allocates its own. For adaptive (heteroscedastic, Issue #217) neighborhood construction,
!| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
!| generates the ascending `(k_start, k_step, k_max)` growth-knob candidates instead, and
!| `calc_adaptive_js_comp_test_bounds` recommends the reference-point capacity their search's
!| per-point arrays are sized by;
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| is that search, the adaptive counterpart of the fixed-k parameter search, with the same
!| plateau criteria. For the selected knobs,
!| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| rebuilds the reference points and their neighborhoods, pooled across all studies, and
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]]
!| runs the final test on them. Together they reproduce exactly what the search traced for that
!| candidate. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| For fixed-k neighborhoods, the final test itself (histograms, JSDs, weighting and the K-study
!| permutation test) is
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]].
module tox_data_integration_js_comp_test_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_math_impl, only: clamp, is_close, LOG_2
    use f42_stats_impl, only: calc_percentile_impl
    use f42_sort_impl, only: init_perm, sort_array_heapsort, sort_real_heapsort_expl_size, binary_search_insertion
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_multinomial
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, validate_dimension_size, &
                          validate_in_range_int, ERR_INVALID_INPUT
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, &
                                             build_residual_histograms_impl, calc_pmf_impl
    use tox_data_integration_preprocessing_impl, only: construct_neighborhoods_ranged_impl, pool_means_impl
    use tox_data_integration_stats_impl, only: gjct_permutation_test_impl
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_candidate_bounds, gather_pooled_neighborhood_residuals, &
             determine_point_bin_count, build_point_study_histogram, &
             estimate_bin_count_impl, determine_bin_count_occupancy_impl, &
             determine_bin_count_occupancy_exhaustive_impl, &
             generate_js_comp_test_candidates_impl, &
             check_neighborhood_overlaps_impl, check_mean_pmf_min_counts_impl, check_plateau_condition_impl, &
             check_effect_size_plateau_condition_impl, create_mean_pmf_impl, create_mean_pmf_only_impl, &
             calc_js_comp_test_n_top_k_jsds, bootstrap_histogram_impl, run_js_comp_test_impl, &
             run_js_comp_test_parameter_search_impl, METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, &
             MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, MODE_PLATEAU_BOTH, KX_FACTORS, MAX_POINTS, MIN_POINTS, &
             GAMMA, MAX_POINT_CANDIDATES, MAX_CANDIDATE_PAIRS, MAX_N_BINS, &
             generate_adaptive_js_comp_test_candidates_impl, calc_adaptive_js_comp_test_bounds, &
             ADAPTIVE_START_FRACTION, ADAPTIVE_GAMMA, ADAPTIVE_K_STEP_FRACTION, ADAPTIVE_K_MAX_FACTOR, &
             ADAPTIVE_K_START_MIN_ABS, ADAPTIVE_POINT_CAPACITY_FACTOR, &
             calc_sorted_slice_mad, construct_adaptive_neighborhoods_impl, ADAPTIVE_MIN_VALID_RESIDUALS, &
             materialize_pooled_neighborhood, run_js_comp_test_adaptive_impl, &
             ADAPTIVE_STOP_TAU, ADAPTIVE_STOP_K_MAX, ADAPTIVE_STOP_RESIDUAL_CAP, ADAPTIVE_STOP_EXHAUSTED, &
             ADAPTIVE_STOP_ZERO_DISPERSION, ADAPTIVE_STOP_TOO_FEW_RESIDUALS, ADAPTIVE_STATUS_OK, &
             ADAPTIVE_STATUS_TOO_FEW_MEANS, ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD, ADAPTIVE_STATUS_TOO_FEW_RESIDUALS, &
             run_js_comp_test_adaptive_parameter_search_impl, ADAPTIVE_CANDIDATE_EVALUATED, &
             ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED, ADAPTIVE_CANDIDATE_EMPTY_STUDY, &
             ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED, ADAPTIVE_CANDIDATE_OVERLAP_FAILED, &
             ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED

    ! `join_method`'s mode table (see check_plateau_condition_impl below). The generator derives
    ! a mode argument's required parameter prefix from the argument's own name -- `join_method`
    ! ends in `_method`, so the prefix is `METHOD_`, not `JOIN_` (verified against
    ! helper/codegen/ir/roles.py's `mode_alias_of`/`_mode_values`); `JOIN_` is kept as part of
    ! each parameter's own name, matching the codebase's existing `baseline_mode` ->
    ! `MODE_BASELINE_RAW` precedent for the same reason.
    integer(int32), parameter :: METHOD_JOIN_MIN = 0_int32
        !! Join method: succeeds only once every study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MAX = 1_int32
        !! Join method: succeeds once any one study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MEDIAN = 2_int32
        !! Join method: succeeds once a majority of studies' confidence-interval overlaps exceed
        !! the threshold, in [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]

    ! Module-private tuning constants, ported verbatim from 125-stabilize-jscomp's
    ! `determine_js_comp_test_n_points_n_neighbors_alloc`. No argument, no directive -- these
    ! never cross a wrapper boundary; where an array would need one of them as a dimension (e.g.
    ! MAX_CANDIDATE_PAIRS below), the dimension is written as the literal value instead, since a
    ! dimension expression that names a symbol the wrapper module never `use`s would not compile
    ! there (confirmed against how `emit/fortran_wrapper.py` renders `argument.dimension.extents`
    ! verbatim, with no import path for a bare dimension token).
    real(real64), parameter :: KX_FACTORS(2) = [0.25_real64, 0.5_real64]
        !! Ascending, as part of the neighbor-candidate denominator: first factor gives the
        !! larger neighbor-candidate value, second the smaller
    integer(int32), parameter :: MAX_POINTS = 1500_int32
    integer(int32), parameter :: MIN_POINTS = 300_int32
    real(real64), parameter :: GAMMA = 0.8_real64
    integer(int32), parameter :: MAX_POINT_CANDIDATES = 8_int32
        !! How often n_points_high may be multiplied by GAMMA before the grid gives up
    integer(int32), parameter :: MAX_CANDIDATE_PAIRS = 16_int32
        !! = size(KX_FACTORS) * MAX_POINT_CANDIDATES = 2 * 8. Kept a literal (rather than the
        !! product expression) because it is also, independently, the literal array-dimension
        !! bound written on candidates_n_points_n_neighbors below -- the two
        !! must agree, and a compile-time PARAMETER expression cannot itself be documented as
        !! equal to a dimension literal the generator reads as plain text.
    integer(int32), parameter :: MAX_N_BINS = 256_int32
        !! Fixed ceiling for a candidate's estimated histogram bin count (see
        !! estimate_bin_count_impl below): later stages of this pipeline size their histogram
        !! work arrays to this fixed bound rather than to a data-dependent one.

    ! Issue #217: tuning constants of the adaptive (heteroscedastic) neighborhood construction's
    ! ascending candidate sequence (generate_adaptive_js_comp_test_candidates_impl below) and of
    ! its capacity producer (calc_adaptive_js_comp_test_bounds). Plain parameters, like GAMMA /
    ! KX_FACTORS above: none of them is an argument default, so none needs a CM_ macro. Every one
    ! of them is relative to the padded pool size `N = max_n_genes_all_studies * n_studies`, so
    ! the sequence scales with the data instead of hitting a fixed floor the way the fixed-k
    ! grid's `n_points_low = max(MIN_POINTS, ...)` does. The candidate cap is
    ! MAX_CANDIDATE_PAIRS (16), shared with the fixed-k grid so both searches' trace arrays have
    ! the same candidate extent. Provisional: the values are to be settled by real-data
    ! validation of the adaptive search (Issue #217).
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

    ! Issue #217: the adaptive neighborhood construction itself
    ! (construct_adaptive_neighborhoods_impl below). ADAPTIVE_MIN_VALID_RESIDUALS is Aaron
    ! Schroeder's literal 10 from his noise model's `gather_residuals_helper`. The stop reasons and
    ! construction statuses are output codes, not mode arguments, so they are plain parameters
    ! documented in prose on the routine rather than in a mode table. Each value is a CM_ macro
    ! first, so that the parameter and the documentation that quotes it share one literal.
#define CM_ADAPTIVE_MIN_VALID_RESIDUALS 10_int32
#define CM_ADAPTIVE_STOP_TAU 1_int32
#define CM_ADAPTIVE_STOP_K_MAX 2_int32
#define CM_ADAPTIVE_STOP_RESIDUAL_CAP 3_int32
#define CM_ADAPTIVE_STOP_EXHAUSTED 4_int32
#define CM_ADAPTIVE_STOP_ZERO_DISPERSION 5_int32
#define CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS 6_int32
#define CM_ADAPTIVE_STATUS_OK 0_int32
#define CM_ADAPTIVE_STATUS_TOO_FEW_MEANS 1_int32
#define CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD 2_int32
#define CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS 3_int32
    integer(int32), parameter :: ADAPTIVE_MIN_VALID_RESIDUALS = CM_ADAPTIVE_MIN_VALID_RESIDUALS
        !! Minimum number of non-NaN residuals a neighborhood must pool after its unconditional
        !! `k_start` phase for its dispersion-controlled growth to run
    integer(int32), parameter :: ADAPTIVE_STOP_TAU = CM_ADAPTIVE_STOP_TAU
        !! Stop reason: the next round would have raised the dispersion by more than `tau`
    integer(int32), parameter :: ADAPTIVE_STOP_K_MAX = CM_ADAPTIVE_STOP_K_MAX
        !! Stop reason: the neighborhood reached `k_max` pooled entries
    integer(int32), parameter :: ADAPTIVE_STOP_RESIDUAL_CAP = CM_ADAPTIVE_STOP_RESIDUAL_CAP
        !! Stop reason: the next entry would have pushed the pooled residuals past `max_pooled_residuals`
    integer(int32), parameter :: ADAPTIVE_STOP_EXHAUSTED = CM_ADAPTIVE_STOP_EXHAUSTED
        !! Stop reason: no pooled entry is left on either side of the neighborhood
    integer(int32), parameter :: ADAPTIVE_STOP_ZERO_DISPERSION = CM_ADAPTIVE_STOP_ZERO_DISPERSION
        !! Stop reason: the `k_start` entries' residuals are all zero, so no relative change is defined
    integer(int32), parameter :: ADAPTIVE_STOP_TOO_FEW_RESIDUALS = CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS
        !! Stop reason: fewer than ADAPTIVE_MIN_VALID_RESIDUALS non-NaN residuals after the
        !! `k_start` phase, so the adaptive phase never ran
    integer(int32), parameter :: ADAPTIVE_STATUS_OK = CM_ADAPTIVE_STATUS_OK
        !! Construction status: every neighborhood was built without a failure
    integer(int32), parameter :: ADAPTIVE_STATUS_TOO_FEW_MEANS = CM_ADAPTIVE_STATUS_TOO_FEW_MEANS
        !! Construction status: fewer non-NaN pooled means than `k_start`; no neighborhood was built
    integer(int32), parameter :: ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD = CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD
        !! Construction status: some neighborhood has fewer than `min_study_neighbors` entries of some study
    integer(int32), parameter :: ADAPTIVE_STATUS_TOO_FEW_RESIDUALS = CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS
        !! Construction status: some neighborhood stopped with ADAPTIVE_STOP_TOO_FEW_RESIDUALS

    ! Defaults of construct_adaptive_neighborhoods_impl's optional knobs, as CM_ macros so that
    ! DM_DEFAULT and M_DEFAULT_VAL share one literal each.
#define CM_ADAPTIVE_TAU_DEFAULT 0.1_real64
#define CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT 1.0_real64
#define CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT 0_int32
#define CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT 1_int32

    ! Issue #217: what the adaptive parameter search
    ! (run_js_comp_test_adaptive_parameter_search_impl below) did with each candidate it tried,
    ! returned in its per-candidate log. Output codes, not a mode argument, so plain parameters,
    ! each a CM_ macro first so the parameter and the documentation quoting it share one literal.
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
        !! Candidate status: some neighborhood has fewer than `min_study_neighbors` entries of some study
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED = CM_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED
        !! Candidate status: more reference points emerged than `max_n_points_candidate`
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_OVERLAP_FAILED = CM_ADAPTIVE_CANDIDATE_OVERLAP_FAILED
        !! Candidate status: two consecutive neighborhoods overlap by less than `min_neighbor_overlap`
    integer(int32), parameter :: ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED = CM_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED
        !! Candidate status: some bin of the consensus pmf has fewer than `min_residuals_per_bin` residuals

    ! `plateau_mode`'s mode table (see run_js_comp_test_parameter_search_impl below). Issue #178:
    ! CI overlap can be too strict a plateau criterion once bootstrap confidence intervals are very
    ! narrow, so MODE_PLATEAU_EFFECT_SIZE/MODE_PLATEAU_BOTH let a caller declare a plateau from
    ! relative-effect-size stability (see check_effect_size_plateau_condition_impl) instead of, or
    ! in addition to, CI overlap. CM_MODE_PLATEAU_CI_OVERLAP is also plateau_mode's DM_DEFAULT: a
    ! genuine Fortran parameter can't be named inside a DM_DEFAULT doc-macro, so the value is
    ! defined once as a preprocessor macro and both the parameter and the default reference it --
    ! same pattern as tox_get_outliers_impl.F90's CM_FAMILY_MODE_DEFAULT.
#define CM_MODE_PLATEAU_CI_OVERLAP 0_int32
    integer(int32), parameter :: MODE_PLATEAU_CI_OVERLAP = CM_MODE_PLATEAU_CI_OVERLAP
        !! Plateau mode: CI overlap only, exactly the pre-Issue-#178 behavior
    integer(int32), parameter :: MODE_PLATEAU_EFFECT_SIZE = 1_int32
        !! Plateau mode: relative-effect-size stability only
    integer(int32), parameter :: MODE_PLATEAU_BOTH = 2_int32
        !! Plateau mode: either CI overlap or relative-effect-size stability

    ! Issue #178's suggested defaults for the relative-effect-size plateau criterion. Explicitly
    ! provisional: the issue itself says "the exact thresholds should be validated empirically",
    ! and that validation is blocked on 2 separate open issues (KX_FACTORS default and the
    ! Freedman-Diaconis bin-count overestimate) -- see the project's JSD-Comp-Test follow-up issue.
    ! (A third candidate blocker, a one-sided-vs-symmetric JSD formula question, was raised in the
    ! same follow-up issue but confirmed by the issue's own author to be a mistake in the issue
    ! text, not a real discrepancy -- the code's symmetric formula is correct as written.) All four
    ! thresholds are exposed as optional arguments precisely so a caller can override them once the
    ! remaining validation happens, without a code change. Defined as CM_ macros, not just
    ! parameters, so DM_DEFAULT and M_DEFAULT_VAL below share one literal each instead of
    ! duplicating it by hand.
#define CM_DELTA_MEDIAN_THRESHOLD_DEFAULT 0.05_real64
#define CM_DELTA_MAX_THRESHOLD_DEFAULT 0.10_real64
#define CM_DELTA_EPSILON_DEFAULT 1.0e-10_real64
#define CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT 2_int32

    ! Issue #187's suggested defaults for occupancy-constrained per-neighborhood histogram binning
    ! (determine_bin_count_occupancy_impl below). CM_OCCUPANCY_GAMMA_DEFAULT is deliberately a
    ! distinct constant from the module's own GAMMA=0.8 above: that one is the point/neighbor
    ! candidate-grid's decay factor (shrinks n_points_high toward n_points_low); this one is the
    ! occupancy bin-count search's own growth factor (grows a candidate M toward m_max), an
    ! unrelated tuning knob that happens to reuse the letter "gamma" in the issue text.
#define CM_OCCUPANCY_M_MIN_DEFAULT 3_int32
#define CM_OCCUPANCY_M_MAX_DEFAULT 120_int32
#define CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT 10_int32
#define CM_OCCUPANCY_GAMMA_DEFAULT 1.25_real64
#define CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.05_real64
#define CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT 0.95_real64

contains

    !> M_EXPORT_C
    !| summary: Recommend upper bounds for the js-comp-test candidate-grid work arrays
    !| AUTHOR_LASZLO_LANG
    !| Closed-form from the GAMMA-decay candidate-grid formula
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| uses: `max_n_points_candidate` is exactly the grid's first (largest) `n_points` candidate.
    !| `max_n_neighbors_candidate` is a safe, not necessarily tight, upper bound on the grid's
    !| largest `n_neighbors` candidate -- reached with the smallest `KX_FACTORS` entry at the
    !| smallest `n_points_high` the grid loop ever uses, which by construction never drops below
    !| `n_points_low`.
    !|
    !| Rejects a non-positive `max_n_genes_all_studies` with invalid input; both outputs are then
    !| left undefined.
    pure subroutine calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, &
                                                        max_n_neighbors_candidate, ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Exact upper bound on the `n_points` candidate the grid ever produces
        integer(int32), intent(out) :: max_n_neighbors_candidate
            !! Safe upper bound on the `n_neighbors` candidate the grid ever produces
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: n_points_low

        call set_ok(ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        if (is_err(ierr)) return

        max_n_points_candidate = clamp(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), &
                                       min_val=MIN_POINTS, max_val=MAX_POINTS)
        n_points_low = max(MIN_POINTS, ceiling(0.2_real64*real(max_n_points_candidate, real64), kind=int32))
        max_n_neighbors_candidate = max(1_int32, floor(real(max_n_genes_all_studies, real64)/ &
                                                        (KX_FACTORS(1)*real(n_points_low, real64)), kind=int32))
    end subroutine calc_js_comp_test_candidate_bounds

    !> M_EXPORT_C
    !| summary: Recommend a reference-point capacity for the adaptive js-comp-test candidate sequence
    !| AUTHOR_LASZLO_LANG
    !| Sizes the per-point work arrays of the adaptive (Issue #217) parameter search, whose
    !| candidates come from
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| An adaptive candidate's reference-point count is not known before its neighborhoods are
    !| grown, so this is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the safety factor
    !| 2 being [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
    !| `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the
    !| smallest `k_start` the candidate sequence contains (the candidate expected to emerge with
    !| the most reference points). The only provable bound is `N` itself, because every new
    !| reference point's seed strictly advances through the pool; sizing every per-point array by
    !| `N` would be prohibitive, so a candidate that emerges with more points than this capacity is
    !| rejected by the adaptive search with a capacity status rather than stored.
    !|
    !| Computed in 64-bit integer / double precision throughout, so `N` itself never overflows.
    !| Enforces the candidate generator's representability bound, which the generator itself
    !| cannot check: the first candidate's `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))`
    !| must not exceed `huge(1_int32)/4 = 536870911`, so that its `k_max = 4*k_start_1` fits a
    !| 32-bit integer; otherwise this routine reports invalid input. The bound is joint in both
    !| arguments, so the error names neither.
    !|
    !| Rejects a `max_n_genes_all_studies` or `n_studies` below 1, or a pair beyond that
    !| representability bound, with invalid input; the output is then left undefined.
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

        integer(int32) :: k_starts(16), n_k_starts
        integer(int64) :: n_pool, capacity

        call set_ok(ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=2_int32, min=1_int32)
        if (is_err(ierr)) return
        ! Joint bound on (G, S): no arg_pos, because tox_errors' position 0 means "not one
        ! argument", and blaming either one alone would name the wrong argument half the time.
        if (adaptive_k_start_1(max_n_genes_all_studies, n_studies) > int(huge(1_int32)/4_int32, int64)) &
            call set_err_once(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        call compute_adaptive_k_starts(max_n_genes_all_studies, n_studies, k_starts, n_k_starts)

        n_pool = int(max_n_genes_all_studies, int64)*int(n_studies, int64)
        capacity = ceiling(ADAPTIVE_POINT_CAPACITY_FACTOR*real(n_pool, real64)/ &
                           real(k_starts(n_k_starts), real64), kind=int64) + 1_int64
        max_n_points_candidate = int(min(n_pool, capacity), kind=int32)
    end subroutine calc_adaptive_js_comp_test_bounds

    !> summary: Estimate the histogram bin count for one (n_points, n_neighbors) candidate
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    !| the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    !| is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    !| range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    !| independently, each clamped on its own to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins and returned as
    !| `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    !| interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    !| skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
    pure subroutine estimate_bin_count_impl(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, n_neighbors, &
                                            shared_residual_range, n_bins, sturges_bins, fd_bins)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), intent(in) :: residuals(n_residuals)
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: residuals_perm(n_residuals)
            !! Sorting permutation for `residuals`, ascending, NaN last
            !! DM_MIN(1_int32)
            !! DM_MAX(n_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! DM_MIN(0.0_real64)
        integer(int32), intent(out) :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            !! to the (clamped) sturges_bins when the interquartile range is too close to zero to
            !! divide by (see the is_close guard below)

        integer(int32) :: i_pool, n_pool, sturges_raw, fd_raw
        real(real64) :: half_bin_width, quartile_25, quartile_75, n_reps_neighborhood

        ! NaN sorts last under the ascending permutation -- find the last non-NaN position by
        ! scanning back from the end, exactly as the already-shipped
        ! determine_shared_residual_range_impl does for the same reason.
        n_pool = n_residuals
        do i_pool = n_pool, 1, -1
            if (ieee_is_nan(residuals(residuals_perm(i_pool)))) then
                n_pool = n_pool - 1
            else
                exit
            end if
        end do

        if (n_pool == 0) then
            n_bins = 1_int32
            sturges_bins = 1_int32
            fd_bins = 1_int32
            return
        end if

        n_reps_neighborhood = real(max_n_reps_all_studies*n_neighbors, kind=real64)

        ! Sturges
        sturges_raw = 1_int32 + nint(log(n_reps_neighborhood)/LOG_2, kind=int32)
        sturges_bins = max(1_int32, min(sturges_raw, MAX_N_BINS))

        ! Freedman-Diaconis
        call calc_percentile_impl(residuals, n_residuals, residuals_perm, 0.25_real64, quartile_25, n_considered=n_pool)
        call calc_percentile_impl(residuals, n_residuals, residuals_perm, 0.75_real64, quartile_75, n_considered=n_pool)

        half_bin_width = (quartile_75 - quartile_25)/(n_reps_neighborhood**(1.0_real64/3.0_real64))
        if (.not. is_close(half_bin_width, 0.0_real64)) then
            fd_raw = nint(shared_residual_range/half_bin_width, kind=int32)
        else
            fd_raw = sturges_raw
        end if
        fd_bins = max(1_int32, min(fd_raw, MAX_N_BINS))

        n_bins = max(sturges_bins, fd_bins)
    end subroutine estimate_bin_count_impl

    !> summary: Determine one neighborhood's occupancy-constrained histogram bin count (Issue #187)
    !| AUTHOR_LASZLO_LANG
    !| Implements Issue #187's two-stage geometric-search-then-local-refinement algorithm for one
    !| neighborhood's pooled residuals (`pooled_residuals`, across all its neighbors and all
    !| studies): find the largest bin count `M` in `[m_min, m_max]` whose equal-width histogram
    !| over `[shared_residual_range_low, shared_residual_range_high]` -- this neighborhood's own
    !| asymmetric range, the `lower_residual_range_quantile`/`upper_residual_range_quantile`
    !| percentiles of its own pooled signed residuals, rather than a single dataset-wide symmetric
    !| range -- has every bin at or above `min_residuals_per_bin` (the occupancy criterion), rather
    !| than the generic
    !| Sturges/Freedman-Diaconis rule
    !| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]] alone
    !| applies, which is why that routine is still called here too -- purely for the
    !| `sturges_bins`/`fd_bins` diagnostic outputs, never for the decision itself.
    !|
    !| Stage 1 grows a candidate `trial_m` geometrically from `m_min` (`next_m =
    !| ceiling(gamma_occupancy*trial_m)`, guaranteed to advance by at least 1 via
    !| `max(trial_m + 1, next_m)`, clamped to `m_max`) until occupancy first fails, recording the
    !| largest admissible `m_valid` and the first inadmissible `m_invalid`; reaching `m_max` while
    !| still admissible returns it immediately, skipping stage 2 entirely. Stage 2 then tests every
    !| integer strictly between `m_valid` and `m_invalid` -- not a binary search, since equal-width
    !| bin boundaries are recomputed for every candidate `M` and occupancy is therefore not
    !| guaranteed monotonic in `M` (Issue #187 is explicit about this) -- and keeps the largest one
    !| that still passes. Both stages reuse the occupancy diagnostics (`min`/`max_bin_occupancy`)
    !| computed for whichever `M` ends up selected, rather than recomputing them a second time
    !| afterward; `mean_bin_occupancy` needs no such bookkeeping, since every non-NaN pooled
    !| residual lands in exactly one bin at any `M`, so it is always exactly
    !| `n_pooled_residuals / selected_n_bins`.
    !|
    !| Per the issue's own FAILURE policy, `occupancy_failed = .true.` (even `m_min` bins could not
    !| satisfy the occupancy criterion, or every pooled residual is NaN) means the neighborhood
    !| should be rejected by the caller rather than built from `selected_n_bins` -- which is still
    !| set to `m_min` in that case, purely so a caller ignoring `occupancy_failed` has *a* value to
    !| build with, never as an indication the search actually found `m_min` admissible.
    !|
    !| The outer geometric search is a genuine sequential `do`/`exit` state machine, not `do
    !| concurrent`: `m_valid`/`m_invalid` accumulate across iterations and each iteration's
    !| continuation depends on the previous one's outcome, exactly the "loops with data-dependent
    !| control flow across iterations" case Fortran_Coding_Guides.pdf Sec 10 carves out as the
    !| deliberate exception to `do concurrent`.
    pure subroutine determine_bin_count_occupancy_impl(pooled_residuals, pooled_residuals_perm, n_residuals, &
                                                        max_n_reps_all_studies, n_neighbors, &
                                                        selected_n_bins, occupancy_failed, shared_residual_range_low, &
                                                        shared_residual_range_high, n_pooled_residuals, &
                                                        min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, &
                                                        sturges_bins, fd_bins, tmp_bin_counts, m_min, m_max, &
                                                        min_residuals_per_bin, gamma_occupancy, &
                                                        lower_residual_range_quantile, upper_residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! DM_MIN(1_int32)
            !! DM_MAX(n_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: selected_n_bins
            !! The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
            !! every bin at or above min_residuals_per_bin; m_min when occupancy_failed
        logical(c_bool), intent(out) :: occupancy_failed
            !! `.true.` iff even m_min bins could not satisfy the occupancy criterion (including
            !! the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
            !! caller should reject this neighborhood rather than build a histogram from
            !! selected_n_bins
        real(real64), intent(out) :: shared_residual_range_low
            !! This neighborhood's own lower residual-range bound (R_low): the
            !! lower_residual_range_quantile percentile of its own pooled signed residuals --
            !! replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            !! occupancy_failed because every pooled residual is NaN
        real(real64), intent(out) :: shared_residual_range_high
            !! This neighborhood's own upper residual-range bound (R_high): the
            !! upper_residual_range_quantile percentile of its own pooled signed residuals --
            !! replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            !! occupancy_failed because every pooled residual is NaN
        integer(int32), intent(out) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count at selected_n_bins; 0 when occupancy_failed
        real(real64), intent(out) :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            !! occupancy_failed
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count at selected_n_bins; 0 when occupancy_failed
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule estimate for this neighborhood's own pooled residuals, from
            !! estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
            !! estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        integer(int32), intent(out) :: tmp_bin_counts(256)
            !! Working array: per-bin counts of whichever candidate bin count is currently being
            !! tested, reused throughout the search (256 = MAX_N_BINS, written as a literal since a
            !! generated wrapper's dummy dimension cannot reference a module parameter)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); if a caller passes
            !! `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
            !! relying on an unconfirmed generator capability to bound one optional argument by
            !! another
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        real(real64) :: actual_gamma_occupancy, actual_lower_residual_range_quantile, actual_upper_residual_range_quantile
        integer(int32) :: i_pool, trial_m, next_m, m_valid, m_invalid, best_m, min_occ, best_min_occ, best_max_occ
        integer(int32) :: combined_n_bins_unused
        real(real64) :: half_span
        logical(c_bool) :: found_valid

        M_DEFAULT_VAL(m_min, actual_m_min, CM_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, CM_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(gamma_occupancy, actual_gamma_occupancy, CM_OCCUPANCY_GAMMA_DEFAULT)
        M_DEFAULT_VAL(lower_residual_range_quantile, actual_lower_residual_range_quantile, CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        M_DEFAULT_VAL(upper_residual_range_quantile, actual_upper_residual_range_quantile, CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        actual_m_max = max(actual_m_min, actual_m_max)

        ! NaN sorts last under the ascending permutation -- same scan-back convention as
        ! estimate_bin_count_impl / determine_shared_residual_range_impl. Moved ahead of the
        ! estimate_bin_count_impl diagnostic call below so n_pooled_residuals is already known when
        ! shared_residual_range_low/high are derived from it.
        n_pooled_residuals = n_residuals
        do i_pool = n_pooled_residuals, 1, -1
            if (ieee_is_nan(pooled_residuals(pooled_residuals_perm(i_pool)))) then
                n_pooled_residuals = n_pooled_residuals - 1_int32
            else
                exit
            end if
        end do

        ! This neighborhood's own asymmetric range, replacing the old dataset-wide symmetric
        ! shared_residual_range -- 0.0/0.0 when every pooled residual is NaN, matching this
        ! routine's existing degenerate-fallback convention (min_bin_occupancy etc. below).
        if (n_pooled_residuals == 0_int32) then
            shared_residual_range_low = 0.0_real64
            shared_residual_range_high = 0.0_real64
        else
            call calc_percentile_impl(pooled_residuals, n_residuals, pooled_residuals_perm, &
                                      actual_lower_residual_range_quantile, shared_residual_range_low, &
                                      n_considered=n_pooled_residuals)
            call calc_percentile_impl(pooled_residuals, n_residuals, pooled_residuals_perm, &
                                      actual_upper_residual_range_quantile, shared_residual_range_high, &
                                      n_considered=n_pooled_residuals)
        end if

        ! Diagnostics only -- estimate_bin_count_impl's own all-NaN branch (n_bins=sturges_bins=
        ! fd_bins=1) already handles n_pool==0 correctly regardless of what half_span is, so this
        ! call stays unconditional: it must run before any early return, so sturges_bins/fd_bins --
        ! both intent(out) -- are always defined. Fed half the new span, not the full span:
        ! estimate_bin_count_impl's own Freedman-Diaconis arithmetic expects a one-sided range (R)
        ! the way shared_residual_range always was before this neighborhood-local range existed.
        half_span = (shared_residual_range_high - shared_residual_range_low)/2.0_real64
        call estimate_bin_count_impl(pooled_residuals, pooled_residuals_perm, n_residuals, max_n_reps_all_studies, &
                                     n_neighbors, half_span, combined_n_bins_unused, sturges_bins, fd_bins)

        if (n_pooled_residuals == 0_int32) then
            occupancy_failed = .true.
            selected_n_bins = actual_m_min
            min_bin_occupancy = 0_int32
            mean_bin_occupancy = 0.0_real64
            max_bin_occupancy = 0_int32
            return
        end if

        found_valid = .false.
        m_valid = -1_int32
        m_invalid = -1_int32
        best_min_occ = 0_int32
        best_max_occ = 0_int32
        trial_m = actual_m_min

        ! Stage 1: geometric coarse search. Sequential by construction -- see the doc block above.
        do
            call histogram_bin_counts(pooled_residuals, pooled_residuals_perm, n_residuals, n_pooled_residuals, &
                                      shared_residual_range_low, shared_residual_range_high, trial_m, &
                                      tmp_bin_counts(1:trial_m), min_occ)

            if (min_occ >= actual_min_residuals_per_bin) then
                m_valid = trial_m
                found_valid = .true.
                best_min_occ = min_occ
                best_max_occ = maxval(tmp_bin_counts(1:trial_m))
            else
                m_invalid = trial_m
                exit
            end if

            if (trial_m == actual_m_max) then
                ! Early return: still admissible at m_max, so stage 2's local refinement has
                ! nothing left to search.
                selected_n_bins = actual_m_max
                occupancy_failed = .false.
                min_bin_occupancy = best_min_occ
                max_bin_occupancy = best_max_occ
                mean_bin_occupancy = real(n_pooled_residuals, real64)/real(actual_m_max, real64)
                return
            end if

            next_m = ceiling(actual_gamma_occupancy*real(trial_m, real64), kind=int32)
            next_m = max(trial_m + 1_int32, next_m) ! guarantee progress after integer rounding
            trial_m = min(next_m, actual_m_max)
        end do

        if (.not. found_valid) then
            ! Even m_min itself is unsupported -- FAILURE per Issue #187's own pseudocode.
            occupancy_failed = .true.
            selected_n_bins = actual_m_min
            min_bin_occupancy = 0_int32
            mean_bin_occupancy = 0.0_real64
            max_bin_occupancy = 0_int32
            return
        end if

        ! Stage 2: exhaustive local refinement over (m_valid, m_invalid). Every integer is tested,
        ! deliberately not a binary search -- occupancy is not guaranteed monotonic in M once bin
        ! boundaries are recomputed per candidate (Issue #187's own warning).
        best_m = m_valid
        do trial_m = m_valid + 1_int32, m_invalid - 1_int32
            call histogram_bin_counts(pooled_residuals, pooled_residuals_perm, n_residuals, n_pooled_residuals, &
                                      shared_residual_range_low, shared_residual_range_high, trial_m, &
                                      tmp_bin_counts(1:trial_m), min_occ)
            if (min_occ >= actual_min_residuals_per_bin) then
                best_m = trial_m
                best_min_occ = min_occ
                best_max_occ = maxval(tmp_bin_counts(1:trial_m))
            end if
        end do

        selected_n_bins = best_m
        occupancy_failed = .false.
        min_bin_occupancy = best_min_occ
        max_bin_occupancy = best_max_occ
        mean_bin_occupancy = real(n_pooled_residuals, real64)/real(best_m, real64)
    end subroutine determine_bin_count_occupancy_impl

    !> Bins already-permutation-sorted pooled residuals into `trial_m` equal-width bins over
    !| `[shared_residual_range_low, shared_residual_range_high]`, this neighborhood's own
    !| asymmetric range, mirroring
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]'s own
    !| clamp-then-bin-index formula, but for a flat 1-D pooled array instead of a 3-D
    !| per-study-per-point one. Not published: a private helper for
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
    !| geometric search and local-refinement loops, called O(log(m_max/m_min)) + refinement times
    !| per neighborhood -- kept separate from the heavier per-study/per-point routine above rather
    !| than reusing it, since that one is called far less often and has the wrong shape for this.
    !| `pure`, allocates nothing, explicit-shape dummies only: `bin_counts`'s own extent is
    !| `trial_m`, which is fine for an internal helper that is never itself wrapped (only `tmp_`
    !| dummies of a wrapped `_impl` routine need a literal dimension).
    !|
    !| Only `pooled_residuals_perm(1:n_pool)` is scanned -- the caller has already established that
    !| any trailing NaN entries sort last and are excluded, exactly as
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| does for its own `n_pooled_residuals`.
    !|
    !| The accumulation loop is a plain sequential `do`, not `do concurrent`: distinct `i_pool`
    !| iterations can increment the very same `bin_counts(bin_idx)`, exactly the same data race
    !| `build_residual_histograms_impl` avoids the same way for its own per-replicate loop.
    pure subroutine histogram_bin_counts(pooled_residuals, pooled_residuals_perm, n_residuals, n_pool, &
                                         shared_residual_range_low, shared_residual_range_high, trial_m, &
                                         bin_counts, min_occupancy)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        integer(int32), intent(in) :: n_pool
            !! Number of non-NaN pooled residuals -- only pooled_residuals_perm(1:n_pool) is used
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! Pooled signed residuals for one neighborhood
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
        real(real64), intent(in) :: shared_residual_range_low
            !! This neighborhood's own lower residual-range bound (R_low)
        real(real64), intent(in) :: shared_residual_range_high
            !! This neighborhood's own upper residual-range bound (R_high)
        integer(int32), intent(in) :: trial_m
            !! Candidate number of equally sized histogram bins in range [R_low, R_high]
        integer(int32), intent(out) :: bin_counts(trial_m)
            !! Absolute count of pooled residuals per bin
        integer(int32), intent(out) :: min_occupancy
            !! minval(bin_counts) -- the occupancy criterion's own test statistic

        real(real64) :: bin_width, clamped_residual
        integer(int32) :: i_pool, bin_idx

        ! Same zero-range guard as build_residual_histograms_impl: a degenerate range (e.g. every
        ! pooled residual is 0, or the 0.0/0.0 all-NaN fallback) falls back to a fixed bin width
        ! instead of dividing by zero or a negative width.
        if (shared_residual_range_high - shared_residual_range_low <= 0.0_real64) then
            bin_width = 1.0_real64
        else
            bin_width = (shared_residual_range_high - shared_residual_range_low)/real(trial_m, real64)
        end if

        bin_counts = 0_int32
        do i_pool = 1, n_pool
            clamped_residual = clamp(pooled_residuals(pooled_residuals_perm(i_pool)), min_val=shared_residual_range_low, &
                                     max_val=shared_residual_range_high)
            bin_idx = min(trial_m, int((clamped_residual - shared_residual_range_low)/bin_width) + 1_int32)
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1_int32
        end do

        min_occupancy = minval(bin_counts)
    end subroutine histogram_bin_counts

    !> summary: Exhaustive brute-force reference implementation of Issue #187's occupancy search
    !| AUTHOR_LASZLO_LANG
    !| Tests every candidate bin count `M` in `[m_min, m_max]` independently and keeps the largest
    !| one whose pooled histogram satisfies the occupancy criterion, instead of
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
    !| own fast geometric-search-then-refinement. Exists purely to validate that routine's result:
    !| occupancy is not guaranteed monotonic in `M` once bin boundaries are recomputed per candidate
    !| (Issue #187 is explicit about this), so a search that stops at the first failure can in
    !| principle miss a larger, independently-admissible `M` the geometric ladder never tries. Takes
    !| `shared_residual_range_low`/`shared_residual_range_high`/`n_pooled_residuals` as direct
    !| inputs, already produced by a prior
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| call, rather than re-deriving them -- this isolates the comparison to just the `M`-selection
    !| algorithm, uncontaminated by a second independent percentile computation.
    !|
    !| Every candidate is independent (no early exit, no state carried between iterations, unlike
    !| the production routine's own Stage 1/Stage 2), so this is a genuine `do concurrent` with
    !| `reduce(max:...)`, not a sequential search: `candidate_bin_counts` is declared local to the
    !| loop (an ordinary MAX_N_BINS-sized local, not a `tmp_` dummy -- precedented by
    !| `run_js_comp_test_impl`'s own `candidates_n_points_n_neighbors` local) so each concurrent
    !| iteration gets its own private scratch instead of racing on a shared buffer sliced by
    !| `trial_m`. `reduce(max:...)` has no "argmax" form, so the winning `M`'s own
    !| min/mean/max_bin_occupancy are recovered with one extra, ordinary (non-concurrent) call to
    !| `histogram_bin_counts` for `best_m` alone once the reduction is done -- trivial cost next to
    !| the search itself.
    !|
    !| `n_pooled_residuals == 0` and a degenerate zero-width range need no special-case branch here:
    !| `histogram_bin_counts` already guards the degenerate range internally, and an all-zero pool
    !| naturally resolves to `occupancy_failed` (or a trivial pass at `min_residuals_per_bin=0`,
    !| with `mean_bin_occupancy=0.0`, no divide-by-zero) -- intentional, not an oversight.
    pure subroutine determine_bin_count_occupancy_exhaustive_impl(pooled_residuals, pooled_residuals_perm, &
                                                                   n_residuals, n_pooled_residuals, &
                                                                   shared_residual_range_low, shared_residual_range_high, &
                                                                   selected_n_bins, occupancy_failed, min_bin_occupancy, &
                                                                   mean_bin_occupancy, max_bin_occupancy, m_min, m_max, &
                                                                   min_residuals_per_bin)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! DM_MIN(1_int32)
            !! DM_MAX(n_residuals)
        integer(int32), intent(in) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! DM_MIN(0_int32)
            !! DM_MAX(n_residuals)
        real(real64), intent(in) :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(real64), intent(in) :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! DM_MIN(shared_residual_range_low)
        integer(int32), intent(out) :: selected_n_bins
            !! The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            !! above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        logical(c_bool), intent(out) :: occupancy_failed
            !! `.true.` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
            !! criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
            !! is not itself 0)
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count at selected_n_bins; 0 when occupancy_failed
        real(real64), intent(out) :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins; 0 when occupancy_failed
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count at selected_n_bins; 0 when occupancy_failed
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
            !! implementation clamps it up to `m_min` internally rather than relying on an
            !! unconfirmed generator capability to bound one optional argument by another
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        integer(int32) :: trial_m, best_m, min_occ
        integer(int32) :: candidate_bin_counts(MAX_N_BINS)

        M_DEFAULT_VAL(m_min, actual_m_min, CM_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, CM_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        actual_m_max = max(actual_m_min, actual_m_max)

        best_m = actual_m_min - 1_int32 ! sentinel: stays below actual_m_min if nothing passes
        do concurrent(trial_m=actual_m_min:actual_m_max) local(candidate_bin_counts, min_occ) &
                shared(pooled_residuals, pooled_residuals_perm, n_residuals, n_pooled_residuals, &
                       shared_residual_range_low, shared_residual_range_high, actual_min_residuals_per_bin) &
                reduce(max:best_m)
            call histogram_bin_counts(pooled_residuals, pooled_residuals_perm, n_residuals, n_pooled_residuals, &
                                      shared_residual_range_low, shared_residual_range_high, trial_m, &
                                      candidate_bin_counts(1:trial_m), min_occ)
            if (min_occ >= actual_min_residuals_per_bin) best_m = max(best_m, trial_m)
        end do

        if (best_m < actual_m_min) then
            occupancy_failed = .true.
            selected_n_bins = actual_m_min
            min_bin_occupancy = 0_int32
            mean_bin_occupancy = 0.0_real64
            max_bin_occupancy = 0_int32
            return
        end if

        selected_n_bins = best_m
        occupancy_failed = .false.
        call histogram_bin_counts(pooled_residuals, pooled_residuals_perm, n_residuals, n_pooled_residuals, &
                                  shared_residual_range_low, shared_residual_range_high, best_m, &
                                  candidate_bin_counts(1:best_m), min_occ)
        min_bin_occupancy = min_occ
        max_bin_occupancy = maxval(candidate_bin_counts(1:best_m))
        mean_bin_occupancy = real(n_pooled_residuals, real64)/real(best_m, real64)
    end subroutine determine_bin_count_occupancy_exhaustive_impl

    !> M_EXPORT_C
    !| summary: Pool one reference point's residuals across every neighbor and every study
    !| AUTHOR_LASZLO_LANG
    !| Given one reference point's own per-study neighbor gene indices (one column of a larger
    !| `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)`, as produced by
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]),
    !| gathers that point's residual values from every neighbor gene, across every study, into one
    !| flat pooled array. This is the exact same pooling
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| perform internally, per reference point, before handing the result to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
    !| own occupancy search -- published so a caller can reconstruct that exact same input directly
    !| on real data and feed it to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
    !| (or to `determine_bin_count_occupancy` itself), to check whether the fast search and the
    !| exhaustive reference ever actually disagree in practice, not just on a synthetic fixture.
    pure subroutine gather_pooled_neighborhood_residuals(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                         n_neighbors, n_studies, neighborhood_indices_point, &
                                                         pooled_residuals, ierr)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN explicitly allowed for missing values;
            !! unvalidated on this export path, so any value passes through as-is
        integer(int32), intent(in) :: neighborhood_indices_point(n_neighbors, n_studies)
            !! Gene indices of one reference point's neighborhood, per study -- one column of a
            !! larger neighborhood_indices_all_studies(n_neighbors, n_points, n_studies), as sliced
            !! by the caller
        real(real64), intent(out) :: pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies)
            !! The pooled residual values for this reference point, across every neighbor and every
            !! study, laid out exactly as a (max_n_reps_all_studies, n_neighbors, n_studies) array
            !! would be
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: i_study, i_neighbor, i_rep, gene_idx, n_predecessors

        call set_ok(ierr)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return

        ! No size(...)-based cross-array shape check here, unlike validate_shift_vectors's pattern:
        ! residuals/neighborhood_indices_point/pooled_residuals are explicit-shape dummies sized by
        ! this routine's own scalar arguments, so size(residuals, 1) is tautologically equal to
        ! max_n_reps_all_studies -- such a check could never fire. Real cross-argument shape
        ! protection for Python/R callers lives one layer up, in the generated wrapper, which
        ! derives all 4 extents from the arrays' own shapes and cross-checks them independently.

        ! Every value of neighborhood_indices_point is a gene index used directly to subscript
        ! residuals(i_rep, gene_idx, i_study) below -- a mismatched/stale index (a real mistake a
        ! caller can make once this is R/Python-callable, not just a synthetic worry) would be a
        ! genuine out-of-bounds read, not merely a wrong-answer bug. A plain sequential `do`, NOT
        ! `do concurrent`: multiple iterations calling set_err_once on the same shared ierr scalar
        ! -- a read-then-conditionally-write, not a plain write -- is not standard-conforming
        ! do-concurrent semantics, unlike the gather loop below it (which only ever writes).
        do i_study = 1, n_studies
            do i_neighbor = 1, n_neighbors
                gene_idx = neighborhood_indices_point(i_neighbor, i_study)
                if (gene_idx < 1 .or. gene_idx > max_n_genes_all_studies) then
                    call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)
                end if
            end do
        end do
        if (is_err(ierr)) return

        do concurrent(i_study=1:n_studies)
            do concurrent(i_neighbor=1:n_neighbors) local(gene_idx, n_predecessors) &
                    shared(neighborhood_indices_point, max_n_reps_all_studies, n_neighbors, i_study)
                gene_idx = neighborhood_indices_point(i_neighbor, i_study)
                n_predecessors = ((i_study - 1)*n_neighbors + (i_neighbor - 1))*max_n_reps_all_studies
                do concurrent(i_rep=1:max_n_reps_all_studies) shared(pooled_residuals, residuals, n_predecessors, &
                                                                     gene_idx, i_study)
                    pooled_residuals(n_predecessors + i_rep) = residuals(i_rep, gene_idx, i_study)
                end do
            end do
        end do
    end subroutine gather_pooled_neighborhood_residuals

    !> summary: Pool one reference point's residuals across all studies and pick its histogram bin count
    !| AUTHOR_LASZLO_LANG
    !| Pass B of the JSD-Comp-Test for ONE reference point. Each study may have its own neighbor
    !| count `point_n_neighbors(i_study)`. This is an internal helper of this module, not
    !| published (like `jct_compute_jsd_pipeline_helper` in `tox_data_integration_jsd_impl`).
    !|
    !| It calls `gather_pooled_neighborhood_residuals` once per study, with `n_studies = 1` and
    !| `n_neighbors = point_n_neighbors(i_study)`, and writes each study's block right after the
    !| previous one in `tmp_pooled_residuals`. With the same neighbor count `k` for every study
    !| this gives exactly the layout of a single multi-study call: neighbor `i_neighbor` of study
    !| `i_study` starts at `((i_study - 1)*k + (i_neighbor - 1))*max_n_reps_all_studies`. The
    !| filled part, `1:max_n_reps_all_studies*sum(point_n_neighbors)`, is then sorted (NaN last) and
    !| handed to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]].
    !|
    !| That routine uses its `n_neighbors` only for the Sturges/Freedman-Diaconis diagnostics
    !| (`sturges_bins`, `fd_bins`). It gets the rounded mean per-study neighbor count
    !| `(sum(point_n_neighbors) + n_studies/2)/n_studies` (integer division), which is exactly `k`
    !| when every study has `k` neighbors.
    !|
    !| Precondition, not checked here: every `point_n_neighbors(i_study)` is at least 1 and at most
    !| `max_n_neighbors`. A gene index outside `[1, max_n_genes_all_studies]` is reported by
    !| `gather_pooled_neighborhood_residuals` itself; its error code (without its argument position,
    !| which numbers that routine's own arguments) is returned in `ierr` and this helper returns at
    !| once, leaving every other output undefined.
    pure subroutine determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                              max_n_neighbors, point_neighborhood_indices, point_n_neighbors, &
                                              selected_n_bins, occupancy_failed, shared_residual_range_low, &
                                              shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                              mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                              tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts, &
                                              ierr, m_min, m_max, min_residuals_per_bin, gamma_occupancy, &
                                              lower_residual_range_quantile, upper_residual_range_quantile)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_neighbors
            !! Leading extent of `point_neighborhood_indices`: the largest neighbor count any study
            !! may have for this reference point
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: point_neighborhood_indices(max_n_neighbors, n_studies)
            !! Gene indices of this reference point's neighborhood, per study. Only the first
            !! `point_n_neighbors(i_study)` entries of column `i_study` are read
        integer(int32), intent(in) :: point_n_neighbors(n_studies)
            !! Number of neighbors of this reference point in each study, each in
            !! `[1, max_n_neighbors]`
        integer(int32), intent(out) :: selected_n_bins
            !! This reference point's bin count, from determine_bin_count_occupancy_impl
        logical(c_bool), intent(out) :: occupancy_failed
            !! `.true.` iff even `m_min` bins could not satisfy the occupancy criterion, from
            !! determine_bin_count_occupancy_impl
        real(real64), intent(out) :: shared_residual_range_low
            !! This reference point's lower residual-range bound (R_low), from
            !! determine_bin_count_occupancy_impl
        real(real64), intent(out) :: shared_residual_range_high
            !! This reference point's upper residual-range bound (R_high), from
            !! determine_bin_count_occupancy_impl
        integer(int32), intent(out) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from determine_bin_count_occupancy_impl
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count at `selected_n_bins`, from determine_bin_count_occupancy_impl
        real(real64), intent(out) :: mean_bin_occupancy
            !! Mean bin count at `selected_n_bins`, from determine_bin_count_occupancy_impl
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count at `selected_n_bins`, from determine_bin_count_occupancy_impl
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule diagnostic, computed with the rounded mean per-study neighbor count
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule diagnostic, computed with the rounded mean per-study neighbor
            !! count
        real(real64), intent(out) :: tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors*n_studies)
            !! Working array: this reference point's pooled residuals; only the leading
            !! `max_n_reps_all_studies*sum(point_n_neighbors)` entries are used
        integer(int32), intent(out) :: tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors*n_studies)
            !! Working array: sorting permutation for the used part of `tmp_pooled_residuals`
        integer(int32), intent(out) :: tmp_bin_counts(MAX_N_BINS)
            !! Working array forwarded to determine_bin_count_occupancy_impl's `tmp_bin_counts`
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT if a gene index is out of range (see above)
        integer(int32), intent(in), optional :: m_min
            !! Forwarded to determine_bin_count_occupancy_impl
        integer(int32), intent(in), optional :: m_max
            !! Forwarded to determine_bin_count_occupancy_impl
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Forwarded to determine_bin_count_occupancy_impl
        real(real64), intent(in), optional :: gamma_occupancy
            !! Forwarded to determine_bin_count_occupancy_impl
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Forwarded to determine_bin_count_occupancy_impl
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Forwarded to determine_bin_count_occupancy_impl

        integer(int32) :: i_study, n_study_neighbors, pool_offset, gather_ierr, diagnostic_n_neighbors

        call set_ok(ierr)

        ! Sequential: each study's block starts where the previous one ended.
        pool_offset = 0_int32
        do i_study = 1, n_studies
            n_study_neighbors = point_n_neighbors(i_study)
            call gather_pooled_neighborhood_residuals(residuals(:, :, i_study:i_study), max_n_reps_all_studies, &
                                                      max_n_genes_all_studies, n_study_neighbors, 1_int32, &
                                                      point_neighborhood_indices(1:n_study_neighbors, i_study:i_study), &
                                                      tmp_pooled_residuals(pool_offset + 1: &
                                                                           pool_offset + max_n_reps_all_studies*n_study_neighbors), &
                                                      gather_ierr)
            if (is_err(gather_ierr)) then
                call set_err_once(ierr, get_err_code(gather_ierr))
                return
            end if
            pool_offset = pool_offset + max_n_reps_all_studies*n_study_neighbors
        end do

        call init_perm(tmp_pooled_residuals_perm(1:pool_offset))
        call sort_array_heapsort(tmp_pooled_residuals(1:pool_offset), tmp_pooled_residuals_perm(1:pool_offset))

        diagnostic_n_neighbors = (sum(point_n_neighbors) + n_studies/2_int32)/n_studies

        call determine_bin_count_occupancy_impl(tmp_pooled_residuals(1:pool_offset), tmp_pooled_residuals_perm(1:pool_offset), &
                                                pool_offset, max_n_reps_all_studies, diagnostic_n_neighbors, &
                                                selected_n_bins, occupancy_failed, shared_residual_range_low, &
                                                shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                                mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                                tmp_bin_counts, m_min, m_max, min_residuals_per_bin, gamma_occupancy, &
                                                lower_residual_range_quantile, upper_residual_range_quantile)
    end subroutine determine_point_bin_count

    !> summary: Build one study's residual histogram for one reference point
    !| AUTHOR_LASZLO_LANG
    !| Pass C of the JSD-Comp-Test for ONE (reference point, study) pair. This is an internal
    !| helper of this module, not published. It copies the residuals of the `n_neighbors`
    !| neighbor genes into `tmp_neighbor_residuals`, then calls
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] with a
    !| single reference point. That routine bins every reference point independently, so the
    !| result is exactly the matching column of a call over all points at once.
    !|
    !| The single-point call writes into local `(1, MAX_N_BINS)` buffers, and their first
    !| `max_n_bins` entries are then copied into `counts`/`pmf`. Entries beyond `n_bins` are zero,
    !| as build_residual_histograms_impl guarantees. NaN residuals are skipped and not counted in
    !| `included_n_reps`.
    !|
    !| Precondition, not checked here: every entry of `neighbor_indices` is in
    !| `[1, max_n_genes_all_studies]`, and `1 <= n_bins <= max_n_bins <= MAX_N_BINS`.
    pure subroutine build_point_study_histogram(residuals_study, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                n_neighbors, neighbor_indices, shared_residual_range_low, &
                                                shared_residual_range_high, n_bins, max_n_bins, counts, pmf, &
                                                included_n_reps, tmp_neighbor_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors of this reference point in this study
        real(real64), intent(in) :: residuals_study(max_n_reps_all_studies, max_n_genes_all_studies)
            !! This study's signed residuals, NaN allowed for missing values
        integer(int32), intent(in) :: neighbor_indices(n_neighbors)
            !! Gene indices of this reference point's neighbors in this study
        real(real64), intent(in) :: shared_residual_range_low
            !! Lower bound of this reference point's histogram range (R_low)
        real(real64), intent(in) :: shared_residual_range_high
            !! Upper bound of this reference point's histogram range (R_high)
        integer(int32), intent(in) :: n_bins
            !! This reference point's bin count
        integer(int32), intent(in) :: max_n_bins
            !! Length of `counts`/`pmf`: the widest bin count of any reference point
        integer(int32), intent(out) :: counts(max_n_bins)
            !! Absolute count of residuals per bin; zero beyond `n_bins`
        real(real64), intent(out) :: pmf(max_n_bins)
            !! `counts` divided by `included_n_reps` (all zero when `included_n_reps` is 0)
        integer(int32), intent(out) :: included_n_reps
            !! Number of non-NaN residuals that were binned
        real(real64), intent(out) :: tmp_neighbor_residuals(max_n_reps_all_studies, n_neighbors)
            !! Working array: the gathered residuals of this reference point's neighbors

        integer(int32) :: i_neighbor
        integer(int32) :: point_n_bins(1), point_included_n_reps(1)
        integer(int32) :: point_counts(1, MAX_N_BINS)
        real(real64) :: point_range_low(1), point_range_high(1)
        real(real64) :: point_pmf(1, MAX_N_BINS)

        do concurrent(i_neighbor=1:n_neighbors) shared(tmp_neighbor_residuals, residuals_study, neighbor_indices)
            tmp_neighbor_residuals(:, i_neighbor) = residuals_study(:, neighbor_indices(i_neighbor))
        end do

        point_n_bins(1) = n_bins
        point_range_low(1) = shared_residual_range_low
        point_range_high(1) = shared_residual_range_high

        call build_residual_histograms_impl(tmp_neighbor_residuals, max_n_reps_all_studies, n_neighbors, 1_int32, &
                                            point_range_low, point_range_high, max_n_bins, point_n_bins, &
                                            point_counts(1:1, 1:max_n_bins), point_pmf(1:1, 1:max_n_bins), &
                                            point_included_n_reps)

        counts = point_counts(1, 1:max_n_bins)
        pmf = point_pmf(1, 1:max_n_bins)
        included_n_reps = point_included_n_reps(1)
    end subroutine build_point_study_histogram

    !> summary: Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
    !| AUTHOR_LASZLO_LANG
    !| Ported from the grid-building half of 125-stabilize-jscomp's
    !| `determine_js_comp_test_n_points_n_neighbors_helper`: starting from an initial
    !| `n_points_high` (clamped between MIN_POINTS and MAX_POINTS), repeatedly multiplies by
    !| GAMMA until it would drop below `n_points_low`, and for each distinct resulting
    !| `n_points` candidate pairs it with up to `size(KX_FACTORS)` distinct `n_neighbors`
    !| candidates. A duplicate `n_points` or `n_neighbors` value (from clamping or
    !| integer rounding) collapses rather than repeating -- this is real, derived behavior the
    !| grid depends on to avoid redundant candidates at small `max_n_genes_all_studies`, not a
    !| bug: a small enough `max_n_genes_all_studies` collapses the whole grid down to exactly one
    !| candidate.
    !|
    !| Issue #187: this routine no longer estimates a per-candidate histogram bin count as a side
    !| effect -- both
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| now compute real per-neighborhood bin counts via
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| once a candidate has passed admissibility, superseding the old global-pool
    !| `estimate_bin_count_impl` estimate this routine used to produce for every candidate
    !| regardless of admissibility.
    pure subroutine generate_js_comp_test_candidates_impl(max_n_genes_all_studies, candidates_n_points_n_neighbors, &
                                                          n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: candidates_n_points_n_neighbors(2, 16)
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate pairs actually filled (at most MAX_CANDIDATE_PAIRS = 16)

        real(real64) :: n_points_high, n_points_low
        integer(int32) :: i_point_candidate, point_candidate, prev_point_candidate
        integer(int32) :: i_neighbor_candidate, neighbor_candidate, prev_neighbor_candidate

        n_points_high = real(clamp(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), &
                                   min_val=MIN_POINTS, max_val=MAX_POINTS), kind=real64)
        n_points_low = real(max(MIN_POINTS, ceiling(0.2_real64*n_points_high, kind=int32)), kind=real64)

        prev_point_candidate = -1_int32
        n_candidates = 0_int32

        ! Sequential by construction: each iteration's dedup against the previous one requires
        ! the loop to run in order, so this is a plain `do`, not `do concurrent`.
        do i_point_candidate = 1, MAX_POINT_CANDIDATES
            if (n_points_high < n_points_low) exit

            point_candidate = nint(n_points_high, kind=int32)
            if (point_candidate /= prev_point_candidate) then
                prev_point_candidate = point_candidate
                prev_neighbor_candidate = -1_int32

                do i_neighbor_candidate = 1, size(KX_FACTORS, kind=int32)
                    neighbor_candidate = max(1_int32, floor(real(max_n_genes_all_studies, real64)/ &
                                                            (KX_FACTORS(i_neighbor_candidate)*n_points_high), kind=int32))

                    if (neighbor_candidate /= prev_neighbor_candidate) then
                        prev_neighbor_candidate = neighbor_candidate

                        n_candidates = n_candidates + 1_int32
                        candidates_n_points_n_neighbors(1, n_candidates) = point_candidate
                        candidates_n_points_n_neighbors(2, n_candidates) = neighbor_candidate
                    end if
                end do
            end if
            n_points_high = n_points_high*GAMMA
        end do
    end subroutine generate_js_comp_test_candidates_impl

    !> summary: Generate the ascending adaptive (k_start, k_step, k_max) candidate sequence
    !| AUTHOR_LASZLO_LANG
    !| The adaptive-neighborhood counterpart (Issue #217) of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]].
    !| An adaptive candidate does not fix its number of reference points: each neighborhood is
    !| grown from `k_start` pooled entries in rounds of `k_step` up to `k_max` while its residual
    !| dispersion stays stable, and the next reference point is placed beyond it, so the point
    !| count emerges from the growth knobs. This routine only generates those knobs. Because
    !| `k_start` shrinks from one candidate to the next, the candidates **ascend** in the number of
    !| reference points they are expected to produce.
    !|
    !| With `N = max_n_genes_all_studies * n_studies` (the padded pool size, not the count of
    !| non-NaN means, so that the capacity producer
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]
    !| and the search agree on the same number) and `floor_k = max(10, 2*n_studies)`:
    !|
    !| - `k_start_1 = max(floor_k, ceiling(0.02 * N))`
    !| - `k_start_t = floor(k_start_1 * 0.8**(t-1))`
    !| - stop once `k_start_t < floor_k`, or once 16 candidates are accepted
    !| - `k_step_t = ceiling(0.25 * k_start_t)`
    !| - `k_max_t = ceiling(4 * k_start_t)`
    !|
    !| The constants are [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_START_FRACTION(variable)]] (0.02),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_GAMMA(variable)]] (0.8),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_K_STEP_FRACTION(variable)]] (0.25),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]] (4),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_K_START_MIN_ABS(variable)]] (10) and
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_CANDIDATE_PAIRS(variable)]] (16).
    !| With them, one decay step lowers any `k_start >= 10` by at least 2, so the values are
    !| strictly descending without any deduplication, every `k_step` is at least 1 and every
    !| `k_max` at least `k_start`.
    !|
    !| Each `k_start_t` is computed in double precision directly from `k_start_1`, not by
    !| rounding the previous candidate again, so rounding never accumulates along the sequence.
    !|
    !| Every knob scales with `N`, so the sequence does not collapse for small data the way the
    !| fixed-k grid does below about 8,743 genes, where its absolute `n_points_low` floor leaves
    !| room for only a single `n_points` value. The first candidate is always accepted, so even a
    !| data set smaller than `k_start_1` (for example 3 genes in 1 study, where `k_start_1 = 10`)
    !| yields exactly one candidate; growing its neighborhoods then fails with a too-few-means
    !| status, which is the adaptive search's to report, not this routine's.
    !|
    !| `N` is formed in 64-bit integer arithmetic and never overflows. The generated values are
    !| representable as long as `k_start_1` does not exceed `huge(1_int32)/4 = 536870911` (so
    !| that `k_max_1` fits a 32-bit integer), i.e. for `N` up to about 2.7e10 and `n_studies` up
    !| to about 2.7e8 -- far beyond any real data set. This routine does not check that bound;
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]],
    !| which the adaptive search needs for its sizing anyway, rejects inputs beyond it.
    pure subroutine generate_adaptive_js_comp_test_candidates_impl(max_n_genes_all_studies, n_studies, &
                                                                   candidates_k_start_k_step_k_max, n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: candidates_k_start_k_step_k_max(3, 16)
            !! Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending (so the
            !! emerging number of reference points ascends)
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate triples actually filled (at least 1, at most 16)

        integer(int32) :: k_starts(16), i_candidate, k_start

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
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
    !| and its capacity producer `calc_adaptive_js_comp_test_bounds`, so the smallest `k_start`
    !| the producer sizes for is by construction the last one the generator emits. The rule is
    !| documented on the generator. Every value is appended without deduplication: since
    !| `(1 - ADAPTIVE_GAMMA) * ADAPTIVE_K_START_MIN_ABS >= 1` (`0.2 * 10 = 2`), one decay step lowers
    !| any `k_start >= floor_k` by at least 1, so consecutive floors always differ. A retuning that
    !| breaks that condition fails `test_adaptive_constants_preclude_dead_clamps` rather than
    !| silently emitting a repeated candidate.
    pure subroutine compute_adaptive_k_starts(max_n_genes_all_studies, n_studies, k_starts, n_k_starts)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies, at least 1
        integer(int32), intent(in) :: n_studies
            !! Number of studies, at least 1
        integer(int32), intent(out) :: k_starts(16)
            !! Accepted `k_start` values, strictly descending; only the first `n_k_starts` are set
        integer(int32), intent(out) :: n_k_starts
            !! Number of accepted `k_start` values, between 1 and MAX_CANDIDATE_PAIRS

        integer(int64) :: floor_k, k_start_1, k_start

        floor_k = adaptive_floor_k(n_studies)
        k_start_1 = adaptive_k_start_1(max_n_genes_all_studies, n_studies)

        ! A plain `do`: the loop stops at the first value below floor_k, whose position is not
        ! known in advance. k_start_1 >= floor_k, so the first value is always accepted.
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

    !> summary: Median and raw median absolute deviation of an already-sorted slice, in O(k)
    !| AUTHOR_LASZLO_LANG
    !| Computes the median and the raw (unscaled) median absolute deviation (MAD) of the `k =
    !| last_pos - first_pos + 1` values `values(perm(first_pos:last_pos))`, which `perm` already
    !| orders ascending, with no work buffer and no sort. Both medians follow
    !| [[f42_stats_impl(module):calc_percentile_impl(interface)]]'s convention: the value at rank
    !| `0.5*(k-1) + 1`, linearly interpolated between its two neighbors (so the median of an even
    !| count is the mean of the middle pair). The absolute deviations are never stored: those of
    !| the entries up to the median's lower rank grow leftward and those after it grow rightward,
    !| so merging the two runs outward from the median (two pointers, ties to the left) visits them
    !| in ascending order, and the merge stops at the rank it needs. A single value has itself as
    !| median and a MAD of 0.
    !|
    !| Used by
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| on a neighborhood's pooled means. Precondition, not checked: `1 <= first_pos <= last_pos <=
    !| n_values`, and the slice contains no NaN.
    pure subroutine calc_sorted_slice_mad(values, n_values, perm, first_pos, last_pos, median, mad)
        integer(int32), intent(in) :: n_values
            !! Number of elements in `values`
        real(real64), intent(in) :: values(n_values)
            !! Values the slice is taken from
        integer(int32), intent(in) :: perm(n_values)
            !! Permutation ordering `values` ascending at least over `first_pos:last_pos`
        integer(int32), intent(in) :: first_pos
            !! First position of the slice in `perm`
        integer(int32), intent(in) :: last_pos
            !! Last position of the slice in `perm`
        real(real64), intent(out) :: median
            !! Median of the slice
        real(real64), intent(out) :: mad
            !! Raw median absolute deviation of the slice from `median`

        integer(int32) :: n_slice, lower_rank, left_pos, right_pos, i_rank
        real(real64) :: rank, fraction, left_dev, right_dev, dev, lower_dev

        n_slice = last_pos - first_pos + 1
        if (n_slice == 1_int32) then
            median = values(perm(first_pos))
            mad = 0.0_real64
            return
        end if

        ! calc_percentile_impl's rank at percentile 0.5. For n_slice >= 2 it is below n_slice,
        ! so both interpolation neighbors lie inside the slice.
        rank = 0.5_real64*real(n_slice - 1_int32, real64) + 1.0_real64
        lower_rank = floor(rank, kind=int32)
        fraction = rank - real(lower_rank, real64)
        left_pos = first_pos + lower_rank - 1_int32
        right_pos = left_pos + 1_int32
        median = values(perm(left_pos)) + fraction*(values(perm(right_pos)) - values(perm(left_pos)))

        ! Sequential two-pointer merge: each step depends on which side the previous one took.
        ! The deviations at ranks lower_rank and lower_rank + 1 are interpolated like the median.
        lower_dev = 0.0_real64
        dev = 0.0_real64
        do i_rank = 1_int32, lower_rank + 1_int32
            if (left_pos >= first_pos .and. right_pos <= last_pos) then
                left_dev = abs(values(perm(left_pos)) - median)
                right_dev = abs(values(perm(right_pos)) - median)
                if (left_dev <= right_dev) then
                    dev = left_dev
                    left_pos = left_pos - 1_int32
                else
                    dev = right_dev
                    right_pos = right_pos + 1_int32
                end if
            else if (left_pos >= first_pos) then
                dev = abs(values(perm(left_pos)) - median)
                left_pos = left_pos - 1_int32
            else
                dev = abs(values(perm(right_pos)) - median)
                right_pos = right_pos + 1_int32
            end if
            if (i_rank == lower_rank) lower_dev = dev
        end do
        mad = lower_dev + fraction*(dev - lower_dev)
    end subroutine calc_sorted_slice_mad

    !> summary: Grow Issue #217's adaptive neighborhoods over an already-sorted pool of gene means
    !| AUTHOR_LASZLO_LANG
    !| The construction of
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]],
    !| documented there, on a permutation the caller has already sorted ascending with NaN last.
    !| Kept separate so the adaptive parameter search can grow every candidate on one shared sort.
    !| `pooled_means` and `pooled_residuals` are `gene_means(G, S)` and `residuals(R, G, S)` seen
    !| through sequence association as `G*S` flat entries, entry `(s - 1)*G + g` being gene `g` of
    !| study `s`.
    pure subroutine grow_adaptive_neighborhoods(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                pooled_means, pooled_residuals, pooled_means_perm, k_start, k_step, &
                                                k_max, tau, mad_distance_factor, max_pooled_residuals, &
                                                min_study_neighbors, n_points, x_star, pooled_neighborhood_range, &
                                                n_neighbors_per_point, stop_reason, neighborhood_dispersion, &
                                                neighborhood_mad, max_n_neighbors, construction_status)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: pooled_means(max_n_genes_all_studies*n_studies)
            !! Every study's gene means, flattened, NaN allowed
        real(real64), intent(in) :: pooled_residuals(max_n_reps_all_studies, max_n_genes_all_studies*n_studies)
            !! Every study's residuals, one column per flat entry of `pooled_means`, NaN allowed
        integer(int32), intent(in) :: pooled_means_perm(max_n_genes_all_studies*n_studies)
            !! Permutation sorting `pooled_means` ascending, NaN last
        integer(int32), intent(in) :: k_start
            !! Entries every neighborhood takes unconditionally, at least 1
        integer(int32), intent(in) :: k_step
            !! Entries staged per adaptive round, at least 1
        integer(int32), intent(in) :: k_max
            !! Largest neighborhood, at least `k_start`
        real(real64), intent(in) :: tau
            !! Largest accepted relative dispersion increase per round, at least 0
        real(real64), intent(in) :: mad_distance_factor
            !! Multiple of a neighborhood's MAD the next target lies beyond it, at least 0
        integer(int32), intent(in) :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals, 0 for none
        integer(int32), intent(in) :: min_study_neighbors
            !! Fewest entries of any one study a neighborhood may have, at least 1
        integer(int32), intent(out) :: n_points
            !! Number of neighborhoods built
        real(real64), intent(out) :: x_star(max_n_genes_all_studies*n_studies)
            !! Reference point of each neighborhood
        integer(int32), intent(out) :: pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies)
            !! `[first, last]` sorted pooled position of each neighborhood
        integer(int32), intent(out) :: n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies)
            !! Entries of each study in each neighborhood
        integer(int32), intent(out) :: stop_reason(max_n_genes_all_studies*n_studies)
            !! Why each neighborhood stopped growing
        real(real64), intent(out) :: neighborhood_dispersion(max_n_genes_all_studies*n_studies)
            !! Mean absolute residual of each final neighborhood, NaN without a non-NaN residual
        real(real64), intent(out) :: neighborhood_mad(max_n_genes_all_studies*n_studies)
            !! Raw MAD of each final neighborhood's pooled means
        integer(int32), intent(out) :: max_n_neighbors
            !! Largest entry count of one study in one neighborhood, 0 without a neighborhood
        integer(int32), intent(out) :: construction_status
            !! Overall outcome of the construction

        integer(int32) :: pool_size, n_valid_means, i_pos, seed_pos, first_pos, last_pos, left_pos, right_pos
        integer(int32) :: n_entries, n_valid, entry_pos, entry_n_valid, n_planned, n_staged, round_n_valid
        integer(int32) :: stage_left, stage_right, stage_first, stage_last, stop_code, target_pos, study_idx
        real(real64) :: x_star_val, sum_abs, entry_sum_abs, round_sum_abs, s_old, s_new, median, mad, target
        logical(c_bool) :: capped, has_too_few_residuals, has_empty_study

        pool_size = max_n_genes_all_studies*n_studies

        ! NaN sorts last: every position after n_valid_means holds a NaN mean.
        n_valid_means = pool_size
        do i_pos = pool_size, 1_int32, -1_int32
            if (ieee_is_nan(pooled_means(pooled_means_perm(i_pos)))) then
                n_valid_means = n_valid_means - 1_int32
            else
                exit
            end if
        end do

        n_points = 0_int32
        max_n_neighbors = 0_int32
        if (n_valid_means < k_start) then
            construction_status = ADAPTIVE_STATUS_TOO_FEW_MEANS
            return
        end if

        has_too_few_residuals = .false.
        has_empty_study = .false.
        seed_pos = 1_int32
        ! Sequential by construction: every seed is placed beyond the previous neighborhood.
        do
            n_points = n_points + 1_int32
            x_star_val = pooled_means(pooled_means_perm(seed_pos))
            first_pos = seed_pos
            last_pos = seed_pos
            left_pos = seed_pos - 1_int32
            right_pos = seed_pos + 1_int32
            n_entries = 1_int32
            call pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, &
                                             pooled_means_perm(seed_pos), sum_abs, n_valid)

            ! Phase 1: the first k_start entries, nearest first, unconditionally. k_start <=
            ! n_valid_means, so a frontier is always left while this loop runs.
            do while (n_entries < k_start)
                entry_pos = nearest_pooled_position(pooled_means, pool_size, pooled_means_perm, n_valid_means, &
                                                    x_star_val, left_pos, right_pos)
                if (entry_pos == left_pos) then
                    left_pos = left_pos - 1_int32
                    first_pos = entry_pos
                else
                    right_pos = right_pos + 1_int32
                    last_pos = entry_pos
                end if
                call pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, &
                                                 pooled_means_perm(entry_pos), entry_sum_abs, entry_n_valid)
                sum_abs = sum_abs + entry_sum_abs
                n_valid = n_valid + entry_n_valid
                n_entries = n_entries + 1_int32
            end do

            if (n_valid < ADAPTIVE_MIN_VALID_RESIDUALS) then
                stop_code = ADAPTIVE_STOP_TOO_FEW_RESIDUALS
                has_too_few_residuals = .true.
            else if (sum_abs == 0.0_real64) then
                stop_code = ADAPTIVE_STOP_ZERO_DISPERSION
            else
                ! Phase 2: rounds of up to k_step entries, each judged against the previous
                ! committed round's dispersion. n_valid >= ADAPTIVE_MIN_VALID_RESIDUALS >= 1 and
                ! sum_abs > 0 here, and a committed round never lowers sum_abs, so s_old > 0.
                s_old = sum_abs/real(n_valid, real64)
                do
                    if (n_entries == k_max) then
                        stop_code = ADAPTIVE_STOP_K_MAX
                        exit
                    end if
                    if (left_pos < 1_int32 .and. right_pos > n_valid_means) then
                        stop_code = ADAPTIVE_STOP_EXHAUSTED
                        exit
                    end if

                    n_planned = min(k_step, k_max - n_entries)
                    stage_left = left_pos
                    stage_right = right_pos
                    stage_first = first_pos
                    stage_last = last_pos
                    round_sum_abs = 0.0_real64
                    round_n_valid = 0_int32
                    n_staged = 0_int32
                    capped = .false.
                    do while (n_staged < n_planned)
                        if (stage_left < 1_int32 .and. stage_right > n_valid_means) exit
                        entry_pos = nearest_pooled_position(pooled_means, pool_size, pooled_means_perm, n_valid_means, &
                                                            x_star_val, stage_left, stage_right)
                        call pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, &
                                                         pooled_means_perm(entry_pos), entry_sum_abs, entry_n_valid)
                        ! Whole entries only: one whose residuals would overflow the cap ends the
                        ! round before it, because nearest-first order may not skip it.
                        if (max_pooled_residuals > 0_int32) then
                            if (n_valid + round_n_valid + entry_n_valid > max_pooled_residuals) then
                                capped = .true.
                                exit
                            end if
                        end if
                        if (entry_pos == stage_left) then
                            stage_left = stage_left - 1_int32
                            stage_first = entry_pos
                        else
                            stage_right = stage_right + 1_int32
                            stage_last = entry_pos
                        end if
                        round_sum_abs = round_sum_abs + entry_sum_abs
                        round_n_valid = round_n_valid + entry_n_valid
                        n_staged = n_staged + 1_int32
                    end do

                    ! A frontier was left at the top of this round, so an empty round can only
                    ! come from the cap.
                    if (n_staged == 0_int32) then
                        stop_code = ADAPTIVE_STOP_RESIDUAL_CAP
                        exit
                    end if

                    ! A round without a non-NaN residual gives exactly s_old back, and commits.
                    s_new = (sum_abs + round_sum_abs)/real(n_valid + round_n_valid, real64)
                    if ((s_new - s_old)/s_old > tau) then
                        stop_code = ADAPTIVE_STOP_TAU
                        exit
                    end if

                    left_pos = stage_left
                    right_pos = stage_right
                    first_pos = stage_first
                    last_pos = stage_last
                    sum_abs = sum_abs + round_sum_abs
                    n_valid = n_valid + round_n_valid
                    n_entries = n_entries + n_staged
                    s_old = s_new
                    if (capped) then
                        stop_code = ADAPTIVE_STOP_RESIDUAL_CAP
                        exit
                    end if
                end do
            end if

            x_star(n_points) = x_star_val
            stop_reason(n_points) = stop_code
            pooled_neighborhood_range(1, n_points) = first_pos
            pooled_neighborhood_range(2, n_points) = last_pos

            ! Decode each pooled position to its study; the studies partition the neighborhood.
            n_neighbors_per_point(:, n_points) = 0_int32
            do i_pos = first_pos, last_pos
                study_idx = (pooled_means_perm(i_pos) - 1_int32)/max_n_genes_all_studies + 1_int32
                n_neighbors_per_point(study_idx, n_points) = n_neighbors_per_point(study_idx, n_points) + 1_int32
            end do
            if (any(n_neighbors_per_point(:, n_points) < min_study_neighbors)) has_empty_study = .true.
            max_n_neighbors = max(max_n_neighbors, maxval(n_neighbors_per_point(:, n_points)))

            if (n_valid > 0_int32) then
                neighborhood_dispersion(n_points) = sum_abs/real(n_valid, real64)
            else
                neighborhood_dispersion(n_points) = M_NAN
            end if
            call calc_sorted_slice_mad(pooled_means, pool_size, pooled_means_perm, first_pos, last_pos, median, mad)
            neighborhood_mad(n_points) = mad

            if (last_pos == n_valid_means) exit

            ! Next seed: the first position after last_pos whose mean is nearest the target, ties
            ! to the lower mean; every branch takes the first position of a run of equal means. Seeds strictly advance, so there are at most n_valid_means points.
            target = pooled_means(pooled_means_perm(last_pos)) + mad_distance_factor*mad
            target_pos = binary_search_insertion(pooled_means, pooled_means_perm, target, &
                                                 lower_idx=last_pos + 1_int32, upper_idx=n_valid_means)
            if (target_pos > n_valid_means) then
                ! Beyond the largest mean: the first position of the largest mean's run.
                seed_pos = binary_search_insertion(pooled_means, pooled_means_perm, &
                                                   pooled_means(pooled_means_perm(n_valid_means)), &
                                                   lower_idx=last_pos + 1_int32, upper_idx=n_valid_means)
            else if (target_pos == last_pos + 1_int32) then
                seed_pos = target_pos
            else if (pooled_means(pooled_means_perm(target_pos)) - target < &
                     target - pooled_means(pooled_means_perm(target_pos - 1_int32))) then
                seed_pos = target_pos
            else
                ! The lower candidate wins: take the first position of its run of equal means,
                ! the lowest position at that distance.
                seed_pos = binary_search_insertion(pooled_means, pooled_means_perm, &
                                                   pooled_means(pooled_means_perm(target_pos - 1_int32)), &
                                                   lower_idx=last_pos + 1_int32, upper_idx=target_pos - 1_int32)
            end if
        end do

        ! Fixed priority, independent of which neighborhood failed first: too few residuals means
        ! growth itself ran without its dispersion rule, an empty study only affects membership.
        if (has_too_few_residuals) then
            construction_status = ADAPTIVE_STATUS_TOO_FEW_RESIDUALS
        else if (has_empty_study) then
            construction_status = ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD
        else
            construction_status = ADAPTIVE_STATUS_OK
        end if
    end subroutine grow_adaptive_neighborhoods

    !> summary: The unvisited pooled position nearest a reference point, ties to the lower mean
    !| AUTHOR_LASZLO_LANG
    !| Chooses between the two frontiers `left_pos` and `right_pos` of a neighborhood growing
    !| around `x_star_val` (a frontier is gone once `left_pos < 1` or `right_pos > n_valid_means`).
    !| Does not move either frontier; the caller advances the one returned. Precondition, not
    !| checked: at least one frontier is left.
    pure integer(int32) function nearest_pooled_position(pooled_means, pool_size, pooled_means_perm, n_valid_means, &
                                                         x_star_val, left_pos, right_pos) result(entry_pos)
        integer(int32), intent(in) :: pool_size
            !! Number of pooled entries
        real(real64), intent(in) :: pooled_means(pool_size)
            !! Pooled gene means
        integer(int32), intent(in) :: pooled_means_perm(pool_size)
            !! Permutation sorting `pooled_means` ascending, NaN last
        integer(int32), intent(in) :: n_valid_means
            !! Number of non-NaN pooled means
        real(real64), intent(in) :: x_star_val
            !! Reference point the neighborhood grows around
        integer(int32), intent(in) :: left_pos
            !! Next candidate position below the neighborhood
        integer(int32), intent(in) :: right_pos
            !! Next candidate position above the neighborhood

        if (left_pos >= 1_int32 .and. right_pos <= n_valid_means) then
            if (abs(pooled_means(pooled_means_perm(left_pos)) - x_star_val) <= &
                abs(pooled_means(pooled_means_perm(right_pos)) - x_star_val)) then
                entry_pos = left_pos
            else
                entry_pos = right_pos
            end if
        else if (left_pos >= 1_int32) then
            entry_pos = left_pos
        else
            entry_pos = right_pos
        end if
    end function nearest_pooled_position

    !> summary: Sum of absolute values and count of one pooled entry's non-NaN residuals
    !| AUTHOR_LASZLO_LANG
    pure subroutine pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, flat_idx, &
                                                entry_sum_abs, entry_n_valid)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: pool_size
            !! Number of pooled entries
        real(real64), intent(in) :: pooled_residuals(max_n_reps_all_studies, pool_size)
            !! Residuals, one column per pooled entry, NaN allowed
        integer(int32), intent(in) :: flat_idx
            !! Flat entry whose residuals are summarized
        real(real64), intent(out) :: entry_sum_abs
            !! Sum of the absolute non-NaN residuals
        integer(int32), intent(out) :: entry_n_valid
            !! Number of non-NaN residuals

        integer(int32) :: i_rep
        real(real64) :: residual

        entry_sum_abs = 0.0_real64
        entry_n_valid = 0_int32
        do concurrent(i_rep=1:max_n_reps_all_studies) local(residual) shared(pooled_residuals, flat_idx) &
            reduce(+:entry_sum_abs, entry_n_valid)
            residual = pooled_residuals(i_rep, flat_idx)
            if (.not. ieee_is_nan(residual)) then
                entry_sum_abs = entry_sum_abs + abs(residual)
                entry_n_valid = entry_n_valid + 1_int32
            end if
        end do
    end subroutine pooled_entry_residual_stats

    !> summary: Grow adaptive, dispersion-controlled neighborhoods over all studies' pooled gene means (Issue #217)
    !| AUTHOR_LASZLO_LANG
    !| The heteroscedastic neighborhood construction of Issue #217, after Aaron Schroeder's
    !| noise-model growth: instead of one fixed neighbor count, each neighborhood grows while the
    !| mean absolute residual of its pool stays stable, and the next reference point is placed
    !| beyond it. The number of reference points emerges from the growth knobs `k_start`, `k_step`
    !| and `k_max`, all counted in pooled entries (one entry is one gene of one study).
    !|
    !| **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
    !| n_studies` flat entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted
    !| ascending; NaN means sort last and take no part. A neighborhood is a contiguous run
    !| `a..b` of that sorted order, returned as `pooled_neighborhood_range`, so its per-study
    !| membership follows by decoding every position back to its study, and the studies' counts
    !| `n_neighbors_per_point` add up to `b - a + 1` exactly. These ranges are what
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !| takes; unlike the fixed-k ranges they are not extended over tied means, because a range is
    !| the exact membership here.
    !|
    !| **Growth of one neighborhood** around its reference point `x_star`, the mean at its seed
    !| position. Entries join nearest first by `|mean - x_star|`, a tie going to the lower mean.
    !| The dispersion `S` is the mean absolute value of every non-NaN residual of every entry in
    !| the pool.
    !|
    !| 1. The first `k_start` entries join unconditionally.
    !| 2. If they hold fewer than
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !|    (`CM_ADAPTIVE_MIN_VALID_RESIDUALS`) non-NaN residuals, the neighborhood stops there with stop reason
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !|    (`CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS`). If their residuals are all zero, it stops with
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !|    (`CM_ADAPTIVE_STOP_ZERO_DISPERSION`), since no relative change is defined.
    !| 3. Otherwise it grows in rounds: the next `min(k_step, k_max - count)` nearest entries are
    !|    staged and `S_new` of the pool with them is compared with the previous round's `S_old`.
    !|    If `(S_new - S_old) / S_old > tau` the round is discarded and growth stops
    !|    ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TAU(variable)]], `CM_ADAPTIVE_STOP_TAU`);
    !|    otherwise the round is committed and `S_old = S_new`. A round without any non-NaN
    !|    residual leaves `S` unchanged and is committed.
    !| 4. Growth also stops once the neighborhood holds `k_max` entries
    !|    ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `CM_ADAPTIVE_STOP_K_MAX`),
    !|    or once no entry is left on either side
    !|    ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `CM_ADAPTIVE_STOP_EXHAUSTED`);
    !|    `k_max` is checked first. A round cut short by the pool's end is still evaluated.
    !| 5. With `max_pooled_residuals > 0`, an adaptive round stages whole entries only while the
    !|    pool's non-NaN residuals stay at or below that cap; NaN residuals do not count toward
    !|    it, and the `k_start` phase is not capped. A round cut short by the cap is evaluated as
    !|    usual and, if committed, ends growth; a round that could stage nothing ends it at once.
    !|    Both stop with
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`CM_ADAPTIVE_STOP_RESIDUAL_CAP`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After neighborhood `i`
    !| ends at sorted position `b_i`, the next target is `mean(b_i) + mad_distance_factor *
    !| MAD_i`, `MAD_i` being the raw median absolute deviation of neighborhood `i`'s pooled means
    !| (see [[tox_data_integration_js_comp_test_impl(module):calc_sorted_slice_mad(subroutine)]]),
    !| and the next seed is the position after `b_i` whose mean is nearest that target (a tie
    !| going to the lower mean; the largest pooled mean if the target lies beyond it); where
    !| several positions share that mean, the seed is the first of them. Seeds therefore
    !| strictly advance, which bounds the number of reference points by `N`; the
    !| construction ends with the neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`CM_ADAPTIVE_STATUS_TOO_FEW_MEANS`) when fewer than `k_start` pooled means are non-NaN (all-NaN means included): no
    !| neighborhood is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction
    !| always completes, so that it can be inspected, and the status reports
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS`) if any neighborhood stopped with too few residuals, else
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD`) if any neighborhood has fewer than `min_study_neighbors` entries of some study, else
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`CM_ADAPTIVE_STATUS_OK`). This
    !| priority is fixed, whichever neighborhood failed first.
    !|
    !| All outputs sized by `N` are filled only for their first `n_points` reference points.
    pure subroutine construct_adaptive_neighborhoods_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                          gene_means, residuals, k_start, k_step, k_max, n_points, &
                                                          x_star, pooled_neighborhood_range, n_neighbors_per_point, &
                                                          stop_reason, neighborhood_dispersion, neighborhood_mad, &
                                                          max_n_neighbors, construction_status, &
                                                          tmp_gene_means_perm_all, tau, mad_distance_factor, &
                                                          max_pooled_residuals, min_study_neighbors)
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
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! DM_MIN(k_start)
        integer(int32), intent(out) :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: stop_reason
            !! Why each neighborhood stopped growing: `CM_ADAPTIVE_STOP_TAU` tau exceeded,
            !! `CM_ADAPTIVE_STOP_K_MAX` `k_max` reached, `CM_ADAPTIVE_STOP_RESIDUAL_CAP` residual cap
            !! reached, `CM_ADAPTIVE_STOP_EXHAUSTED` pool exhausted, `CM_ADAPTIVE_STOP_ZERO_DISPERSION`
            !! zero dispersion, `CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS` too few residuals
            !! DM_RESULT_SIZE_IS(n_points)
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood (the last committed
            !! dispersion); NaN when the neighborhood holds no non-NaN residual at all
            !! DM_RESULT_SIZE_IS(n_points)
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), intent(out) :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(int32), intent(out) :: construction_status
            !! Overall outcome: `CM_ADAPTIVE_STATUS_OK` ok, `CM_ADAPTIVE_STATUS_TOO_FEW_MEANS` too few
            !! means, `CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD` a study neighborhood below
            !! `min_study_neighbors`, `CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS` a neighborhood with too
            !! few residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, seeded and sorted here
        real(real64), intent(in), optional :: tau
            !! Largest relative increase of the dispersion an adaptive round may cause and still
            !! be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's median absolute deviation that the next reference
            !! point's target lies beyond the neighborhood's largest mean
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on the non-NaN residuals adaptive rounds may grow a neighborhood's pool to; 0
            !! for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)

        real(real64) :: actual_tau, actual_mad_distance_factor
        integer(int32) :: actual_max_pooled_residuals, actual_min_study_neighbors

        M_DEFAULT_VAL(tau, actual_tau, CM_ADAPTIVE_TAU_DEFAULT)
        M_DEFAULT_VAL(mad_distance_factor, actual_mad_distance_factor, CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        M_DEFAULT_VAL(max_pooled_residuals, actual_max_pooled_residuals, CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        M_DEFAULT_VAL(min_study_neighbors, actual_min_study_neighbors, CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)

        ! The same pooled sort as the fixed-k search, so an adaptive search can share this order.
        call init_perm(tmp_gene_means_perm_all)
        call sort_real_heapsort_expl_size(gene_means, tmp_gene_means_perm_all, max_n_genes_all_studies*n_studies)

        call grow_adaptive_neighborhoods(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                         residuals, tmp_gene_means_perm_all, k_start, k_step, k_max, actual_tau, &
                                         actual_mad_distance_factor, actual_max_pooled_residuals, &
                                         actual_min_study_neighbors, n_points, x_star, pooled_neighborhood_range, &
                                         n_neighbors_per_point, stop_reason, neighborhood_dispersion, &
                                         neighborhood_mad, max_n_neighbors, construction_status)
    end subroutine construct_adaptive_neighborhoods_impl

    !> summary: Decode one adaptive neighborhood's pooled positions into per-study gene lists
    !| AUTHOR_LASZLO_LANG
    !| The per-study membership of one neighborhood built by
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]:
    !| every sorted pooled position `first_pos..last_pos` names the flat entry `flat =
    !| perm_all(pos)` of `gene_means(G, S)`, which is gene `mod(flat - 1, G) + 1` of study
    !| `(flat - 1)/G + 1`. Each study's genes are appended to its column of
    !| `point_neighborhood_indices` in pooled order (ascending mean), and `point_n_neighbors`
    !| counts them, so `sum(point_n_neighbors) = last_pos - first_pos + 1`. A study with no entry
    !| in the range gets a count of 0. This is an internal helper of this module, not published.
    !|
    !| Precondition, not checked here: `1 <= first_pos <= last_pos <= max_n_genes_all_studies *
    !| n_studies`, and `max_n_neighbors` is at least the largest per-study count, which
    !| `max_n_neighbors = max_n_genes_all_studies` always satisfies (a study has at most `G`
    !| distinct entries). Only the first `point_n_neighbors(i_study)` entries of column `i_study`
    !| are defined on return.
    pure subroutine materialize_pooled_neighborhood(perm_all, max_n_genes_all_studies, n_studies, first_pos, last_pos, &
                                                    max_n_neighbors, point_neighborhood_indices, point_n_neighbors)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: perm_all(max_n_genes_all_studies*n_studies)
            !! Permutation sorting the flattened `gene_means(G, S)` ascending, NaN last
        integer(int32), intent(in) :: first_pos
            !! First sorted pooled position of the neighborhood
        integer(int32), intent(in) :: last_pos
            !! Last sorted pooled position of the neighborhood
        integer(int32), intent(in) :: max_n_neighbors
            !! Leading extent of `point_neighborhood_indices`
        integer(int32), intent(out) :: point_neighborhood_indices(max_n_neighbors, n_studies)
            !! Gene indices of the neighborhood, per study, in pooled order; only the first
            !! `point_n_neighbors(i_study)` entries of column `i_study` are defined
        integer(int32), intent(out) :: point_n_neighbors(n_studies)
            !! Number of the neighborhood's entries that belong to each study

        integer(int32) :: i_pos, flat_idx, study_idx

        point_n_neighbors = 0_int32
        ! Sequential: each study's next free slot depends on the entries decoded before it.
        do i_pos = first_pos, last_pos
            flat_idx = perm_all(i_pos) - 1_int32
            study_idx = flat_idx/max_n_genes_all_studies + 1_int32
            point_n_neighbors(study_idx) = point_n_neighbors(study_idx) + 1_int32
            point_neighborhood_indices(point_n_neighbors(study_idx), study_idx) = mod(flat_idx, max_n_genes_all_studies) + 1_int32
        end do
    end subroutine materialize_pooled_neighborhood

    !> summary: Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `test_neighborhood_overlaps_helper`: the first
    !| admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, using the
    !| `[min_idx, max_idx]` neighborhood spans
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| produces. Named `check_*` rather than 125's `test_*`, a deliberate deviation from the
    !| verbatim port: the generated R binding is published under the Fortran name, and the R test
    !| harness (`r/test_helpers.R`'s `run_all_tests`) discovers every `test_`-prefixed name in the
    !| environment as a test case to run with no arguments -- a `test_`-prefixed export would be
    !| swept up and fail every R suite that sources the package, not just this module's own.
    pure subroutine check_neighborhood_overlaps_impl(neighborhood_range, n_points, min_neighbor_overlap, &
                                                      all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: neighborhood_range(2, n_points)
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span, as produced
            !! by construct_neighborhoods_ranged_impl
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`

        integer(int32) :: i_point
        real(real64) :: overlap
        logical(c_bool) :: all_pass

        all_pass = .true.
        do concurrent(i_point=1:n_points - 1) local(overlap) shared(neighborhood_range, min_neighbor_overlap) &
                reduce(.and.:all_pass)
            overlap = compute_fractional_overlap( &
                     real(neighborhood_range(1, i_point), real64), real(neighborhood_range(2, i_point), real64), &
                     real(neighborhood_range(1, i_point + 1), real64), real(neighborhood_range(2, i_point + 1), real64))
            all_pass = all_pass .and. (overlap >= min_neighbor_overlap)
        end do
        all_have_min_neighbor_overlap = logical(all_pass, kind=c_bool)
    end subroutine check_neighborhood_overlaps_impl

    !> summary: Test whether every bin of a mean pmf reaches a minimum absolute count
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `test_mean_pmf_min_counts_helper`: the second
    !| admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, checked once the
    !| first gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]])
    !| has already passed. Named `check_*` rather than 125's `test_*` for the same reason as its
    !| sibling above: a `test_`-prefixed R export collides with the R test harness's own
    !| test-discovery convention.
    !| Issue #187: `n_bins_per_point` scopes the reduction to each point's own valid bin range, so
    !| a point whose own bin count is narrower than `n_bins` (this array's second extent, i.e. the
    !| widest bin count any point uses) has its legitimate zero-padded columns skipped rather than
    !| mistaken for an occupancy failure.
    pure subroutine check_mean_pmf_min_counts_impl(mean_pmf_counts, n_bins, n_bins_per_point, n_points, min_count, &
                                                   all_bins_have_min_count)
        integer(int32), intent(in) :: n_bins
            !! Second extent of `mean_pmf_counts` -- the widest bin count any point uses (max_n_bins),
            !! not necessarily every point's own bin count
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! This point's own bin count -- only `mean_pmf_counts(1:n_bins_per_point(i_point), i_point)`
            !! is inspected; columns beyond it are legitimate zero-padding, not failures
            !! DM_MIN(1_int32)
            !! DM_MAX(n_bins)
        integer(int32), intent(in) :: mean_pmf_counts(n_bins, n_points)
            !! Absolute counts of a residual per bin for the mean pmf
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! DM_MIN(0_int32)
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` if every bin within each reference point's own `n_bins_per_point`, at every
            !! reference point, reaches at least `min_count`

        integer(int32) :: i_point, i_bin
        logical(c_bool) :: all_pass

        all_pass = .true.
        do concurrent(i_point=1:n_points, i_bin=1:n_bins) shared(mean_pmf_counts, min_count, n_bins_per_point) &
                reduce(.and.:all_pass)
            if (i_bin <= n_bins_per_point(i_point)) then
                all_pass = all_pass .and. (mean_pmf_counts(i_bin, i_point) >= min_count)
            end if
        end do
        all_bins_have_min_count = logical(all_pass, kind=c_bool)
    end subroutine check_mean_pmf_min_counts_impl

    !> Calculates the fractional overlap between two intervals: the fraction of the overlap over
    !| the total range of the first interval, so for `a_min < b_min`,
    !| \[ \frac{\min(\texttt{a\_max}, \texttt{b\_max}) - \texttt{b\_min}}{\texttt{a\_max} - \texttt{a\_min}} \]
    !| Not published: a private helper for
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]],
    !| ported from 125-stabilize-jscomp's `compute_fractional_overlap_helper`.
    pure real(real64) function compute_fractional_overlap(a_min, a_max, b_min, b_max) result(overlap_percent)
        real(real64), intent(in) :: a_min
            !! Lower bound of the first interval
        real(real64), intent(in) :: a_max
            !! Upper bound of the first interval, assumed >= a_min
        real(real64), intent(in) :: b_min
            !! Lower bound of the second interval
        real(real64), intent(in) :: b_max
            !! Upper bound of the second interval, assumed >= b_min

        real(real64) :: left_max, right_min

        if (a_max == a_min) then
            ! A zero-width first interval: maximum overlap if it sits inside the second interval
            if (b_min <= a_max .and. b_max >= a_max) then
                overlap_percent = 1.0_real64
            else
                overlap_percent = 0.0_real64
            end if
        else
            left_max = min(a_max, b_max)
            right_min = max(a_min, b_min)
            ! Negative exactly when a_min < a_max < b_min < b_max, i.e. no overlap at all
            overlap_percent = max(0.0_real64, (left_max - right_min)/(a_max - a_min))
        end if
    end function compute_fractional_overlap

    !> summary: Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `check_plateau_condition_helper`. A search over
    !| candidates (finest resolution to coarsest) stops -- "plateaus" -- either when a new
    !| candidate is no better than the previous best (the short-circuit below: keep the previous
    !| best and stop searching), or once the new candidate's confidence-interval overlap with the
    !| previous best meets the condition `join_method` names. `join_method` replaces
    !| 125-stabilize-jscomp's hand-rolled join-method-range validation macro entirely: the mode
    !| table below is itself the validation, checked against exactly the values it names.
    pure subroutine check_plateau_condition_impl(confidence_interval, best_candidate_pair_confidence_interval, n_studies, &
                                                 best_candidate_index, best_exceeded_ci_overlap_count, candidate_index, &
                                                 join_method, succeeding_ci_overlap, plateau_found)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: confidence_interval(2, n_studies)
            !! JSD confidence interval `[lower, upper]` from bootstrapping, for the candidate
            !! pair under test
        real(real64), intent(inout) :: best_candidate_pair_confidence_interval(2, n_studies)
            !! JSD confidence intervals for the current best candidate pair; overwritten with
            !! `confidence_interval` unless the new candidate is worse
        integer(int32), intent(inout) :: best_candidate_index
            !! Candidate-grid index of the current best candidate pair; overwritten with
            !! `candidate_index` unless the new candidate is worse
        integer(int32), intent(inout) :: best_exceeded_ci_overlap_count
            !! Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            !! best candidate pair; overwritten unless the new candidate is worse
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: candidate_index
            !! Candidate-grid index of the candidate pair that produced `confidence_interval`
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition: METHOD_JOIN_MIN requires every study's overlap to exceed
            !! `succeeding_ci_overlap`, METHOD_JOIN_MAX requires only one study's overlap to
            !! exceed it, and METHOD_JOIN_MEDIAN requires a majority
            !! (`count > (n_studies - 1) / 2`) to exceed it
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum fractional overlap an interval in `confidence_interval` must have with its
            !! respective interval in `best_candidate_pair_confidence_interval` to count as
            !! "exceeded"
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate

        integer(int32) :: exceeds_min_ci_overlap, i_study
        real(real64) :: overlap
        logical(c_bool) :: plateau

        exceeds_min_ci_overlap = 0_int32
        do concurrent(i_study=1:n_studies) local(overlap) shared(confidence_interval, &
                                                                  best_candidate_pair_confidence_interval, &
                                                                  succeeding_ci_overlap) reduce(+:exceeds_min_ci_overlap)
            ! IMPORTANT: the current best interval is the *second* pair passed to
            ! compute_fractional_overlap, since its range is the denominator -- a candidate whose
            ! interval sits inside the best one so far scores 1.0 (no worse), a wider one scores
            ! below the threshold.
            overlap = compute_fractional_overlap(confidence_interval(1, i_study), confidence_interval(2, i_study), &
                                                 best_candidate_pair_confidence_interval(1, i_study), &
                                                 best_candidate_pair_confidence_interval(2, i_study))
            if (overlap >= succeeding_ci_overlap) exceeds_min_ci_overlap = exceeds_min_ci_overlap + 1_int32
        end do

        if (exceeds_min_ci_overlap < best_exceeded_ci_overlap_count) then
            plateau = .true.
        else
            best_candidate_index = candidate_index
            best_exceeded_ci_overlap_count = exceeds_min_ci_overlap
            best_candidate_pair_confidence_interval = confidence_interval

            select case (join_method)
            case (METHOD_JOIN_MIN)
                plateau = exceeds_min_ci_overlap == n_studies
            case (METHOD_JOIN_MAX)
                plateau = exceeds_min_ci_overlap > 0_int32
            case default ! METHOD_JOIN_MEDIAN
                plateau = exceeds_min_ci_overlap >= (n_studies - 1_int32)/2_int32 + 1_int32
            end select
        end if
        plateau_found = logical(plateau, kind=c_bool)
    end subroutine check_plateau_condition_impl

    !> summary: Test one candidate's per-study JSD against the previous admissible candidate for a relative-effect-size plateau
    !| AUTHOR_LASZLO_LANG
    !| Implements Issue #178's relative-effect-size plateau criterion, complementary to
    !| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]'s
    !| CI-overlap one: for each study `i`, the relative change in observed JSD between successive
    !| ADMISSIBLE parameter settings (both admissibility gates already passed),
    !| `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    !| delta_epsilon)`, summarized across studies by its median (`delta_median`, via the
    !| already-shipped
    !| [[f42_stats_impl(module):calc_percentile_impl(interface)]]) and maximum (`delta_max`). A
    !| plateau is declared once both stay under their respective thresholds for
    !| `delta_min_consecutive_transitions` consecutive transitions in a row -- tracked across calls
    !| via `n_consecutive_ok`, reset the moment either threshold is missed.
    !|
    !| No transition exists for the very first admissible candidate a caller ever passes in
    !| (`has_previous = .false.`): `delta`/`delta_median`/`delta_max` are all set to
    !| `-1.0_real64` -- the same not-yet-computed sentinel
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| already uses for its own confidence-interval fallback, usable here for the same reason:
    !| every quantity this routine tracks is structurally non-negative.
    !|
    !| The 0.05/0.10 defaults `run_js_comp_test_parameter_search_impl` passes for
    !| `delta_median_threshold`/`delta_max_threshold` are Issue #178's own suggested starting
    !| point, explicitly not yet empirically validated -- see that routine's doc comment.
    pure subroutine check_effect_size_plateau_condition_impl(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                             has_previous, delta_median_threshold, delta_max_threshold, &
                                                             delta_epsilon, delta_min_consecutive_transitions, &
                                                             n_consecutive_ok, delta, delta_median, delta_max, &
                                                             plateau_found, tmp_delta_perm)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: global_js_divergence(n_studies)
            !! Current admissible candidate's observed global JSD per study
            !! DM_MIN(0.0_real64)
        real(real64), intent(in) :: prev_global_js_divergence(n_studies)
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! DM_MIN(0.0_real64)
        logical(c_bool), intent(in) :: has_previous
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! DM_MIN(above(0.0_real64))
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! DM_MIN(1_int32)
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
            !! DM_MIN(0_int32)
        real(real64), intent(out) :: delta(n_studies)
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(real64), intent(out) :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(real64), intent(out) :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(int32), intent(out) :: tmp_delta_perm(n_studies)
            !! Working array: sorting permutation for `delta`, used to compute `delta_median`

        integer(int32) :: i_study

        if (.not. has_previous) then
            delta = -1.0_real64
            delta_median = -1.0_real64
            delta_max = -1.0_real64
            n_consecutive_ok = 0_int32
            plateau_found = logical(.false., kind=c_bool)
            return
        end if

        do concurrent(i_study=1:n_studies) shared(delta, global_js_divergence, prev_global_js_divergence, delta_epsilon)
            delta(i_study) = abs(global_js_divergence(i_study) - prev_global_js_divergence(i_study)) &
                            /max(prev_global_js_divergence(i_study), delta_epsilon)
        end do

        call init_perm(tmp_delta_perm)
        call sort_array_heapsort(delta, tmp_delta_perm)
        call calc_percentile_impl(delta, n_studies, tmp_delta_perm, 0.5_real64, delta_median)
        delta_max = maxval(delta)

        if (delta_median < delta_median_threshold .and. delta_max < delta_max_threshold) then
            n_consecutive_ok = n_consecutive_ok + 1_int32
        else
            n_consecutive_ok = 0_int32
        end if

        plateau_found = logical(n_consecutive_ok >= delta_min_consecutive_transitions, kind=c_bool)
    end subroutine check_effect_size_plateau_condition_impl

    !> summary: Build the consensus pmf and its histogram counts from all studies' pmfs
    !| AUTHOR_LASZLO_LANG
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_helper`.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    pure subroutine create_mean_pmf_impl(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                                         mean_pmf_included_n_reps, mean_pmf_counts)
        integer(int32), intent(in) :: n_bins
            !! The array's first extent for `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` -- the
            !! widest histogram bin count of any reference point (`max_n_bins`, for pmfs/counts
            !! built by `build_residual_histograms_impl` from its own per-point
            !! `n_bins_per_point`). Every study must share the same per-point bin count for this
            !! to be safe: their zero-padded columns beyond that count then coincide across all
            !! `n_studies`, so the averaged/summed `mean_pmf`/`mean_pmf_counts` are zero there too
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_bins, n_points, n_studies), intent(in) :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), dimension(n_bins, n_points, n_studies), intent(in) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
            !! DM_MIN(0_int32)
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for `mean_pmf`,
            !! summed across all n_studies
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts, dim=3)`

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64
        mean_pmf_included_n_reps = 0_int32
        mean_pmf_counts = 0_int32

        ! Sequential over i_study/i_point: each accumulates into the same (i_bin, i_point) element
        ! across studies, so only the innermost i_bin loop -- which writes distinct elements -- can
        ! run concurrently, exactly matching build_residual_histograms_impl's own accumulation
        ! pattern above.
        do i_study = 1, n_studies
            do i_point = 1, n_points
                do concurrent(i_bin=1:n_bins) shared(mean_pmf, pmfs, i_point, i_study, n_studies, mean_pmf_counts, counts)
                    mean_pmf(i_bin, i_point) = mean_pmf(i_bin, i_point) + pmfs(i_bin, i_point, i_study)/real(n_studies, real64)
                    mean_pmf_counts(i_bin, i_point) = mean_pmf_counts(i_bin, i_point) + counts(i_bin, i_point, i_study)
                end do
                mean_pmf_included_n_reps(i_point) = mean_pmf_included_n_reps(i_point) + included_n_reps(i_point, i_study)
            end do
        end do
    end subroutine create_mean_pmf_impl

    !> summary: Build only the consensus pmf from all studies' pmfs, without its histogram counts
    !| AUTHOR_LASZLO_LANG
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_only_helper`: useful where the
    !| mean pmf's own counts don't matter, e.g. for the bootstrap confidence interval a later
    !| stage of this port adds.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    pure subroutine create_mean_pmf_only_impl(pmfs, n_bins, n_points, n_studies, mean_pmf)
        integer(int32), intent(in) :: n_bins
            !! The array's first extent for `pmfs`/`mean_pmf` -- the widest histogram bin count of
            !! any reference point (`max_n_bins`, for pmfs built by
            !! `build_residual_histograms_impl` from its own per-point `n_bins_per_point`). Every
            !! study must share the same per-point bin count for this to be safe: their
            !! zero-padded columns beyond that count then coincide across all `n_studies`, so the
            !! averaged `mean_pmf` is zero there too
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_bins, n_points, n_studies), intent(in) :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64

        do i_study = 1, n_studies
            do i_point = 1, n_points
                do concurrent(i_bin=1:n_bins) shared(mean_pmf, pmfs, i_point, i_study, n_studies)
                    mean_pmf(i_bin, i_point) = mean_pmf(i_bin, i_point) + pmfs(i_bin, i_point, i_study)/real(n_studies, real64)
                end do
            end do
        end do
    end subroutine create_mean_pmf_only_impl

    !> M_EXPORT_C
    !| summary: Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's inline `n_bootstrapping_top_k_jsds` computation in
    !| `determine_js_comp_test_n_points_n_neighbors_alloc`: sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| to track a two-sided bootstrap confidence interval's endpoints. E.g. for a 95% two-sided
    !| interval (2.5% reserved at each end, the default) with `n_bootstraps=1000`,
    !| `n_top_k = max(1, floor(0.025*1000)) = 25`.
    pure subroutine calc_js_comp_test_n_top_k_jsds(n_bootstraps, two_sided_bootstrapping_significance_level, n_top_k)
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples that will be performed
            !! DM_MIN(1_int32)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level, as a percentage reserved at each end of the bootstrap
            !! distribution (e.g. 2.5 for a 95% two-sided interval)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(out) :: n_top_k
            !! Number of elements to keep at each end of the bootstrap distribution, at least 1

        real(real64) :: sig_level

        M_DEFAULT_VAL(two_sided_bootstrapping_significance_level, sig_level, 2.5_real64)

        n_top_k = max(1_int32, floor(sig_level/100.0_real64*real(n_bootstraps, real64), kind=int32))
    end subroutine calc_js_comp_test_n_top_k_jsds

    !> summary: Bootstrap a confidence interval for each study's global JSD by resampling the pooled consensus histogram
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `bootstrap_histogram_helper`. Resamples from the
    !| POOLED/consensus histogram counts (`mean_pmf_counts`), not from each study's own histogram
    !| -- this is 125's deliberate design, preserved as-is (see the plan for this port). Draws
    !| random numbers via
    !| [[f42_random_gsl(module):random_multinomial(subroutine)]], so this implementation is
    !| deliberately impure, matching the project's existing precedent for the other
    !| permutation-test module's purity
    !| ([[tox_trajectory_contribution_analysis_impl(module):perform_permutation_test_impl(interface)]]).
    !|
    !| `confidence_interval` is both an input and an output: its incoming `[lower, upper]` values
    !| seed every slot of the top-k/bottom-k heaps (`tmp_bootstrapping_top_k_jsds`) -- ported
    !| verbatim from 125, which fills the whole heap with the *same* incoming reference value
    !| rather than the usual plus/minus-infinity heap initialization, so the observed
    !| (pre-bootstrap) value can only be displaced by a strictly more extreme bootstrap draw. On
    !| return it holds `[largest of the n_bootstrapping_top_k_jsds smallest bootstrap draws,
    !| smallest of the n_bootstrapping_top_k_jsds largest bootstrap draws]`.
    !|
    !| `mean_pmf_counts`/`tmp_pmfs`/`tmp_mean_pmf` are laid out bin-major (`(n_bins, n_points[,
    !| n_studies])`), matching
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] and its
    !| own `create_mean_pmf_only_impl` above -- both ported from 125, whose own convention this
    !| is. The already-shipped
    !| [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]
    !| predates 125's own port and instead takes its pmf arguments point-major
    !| (`(n_points, n_bins)`); the two calls below bridge the two conventions with an explicit
    !| `transpose`, rather than picking one shape and silently reinterpreting the other's memory
    !| under it (which would scramble every non-square `(n_bins, n_points)` histogram).
    !|
    !| A GSL allocation failure in `create_rng` is a genuine runtime error no input check could
    !| have foreseen (codegen_guide.md Sec 5.14): every work array and `confidence_interval` are
    !| then left untouched (arrays not yet written to keep their caller-visible defined state) and
    !| `ierr` reports `ERR_ALLOC_FAIL`. A `random_multinomial` draw failing (which validated,
    !| internally-consistent inputs should never trigger) is likewise folded into `ierr`, first
    !| failure only, without stopping the resampling already in flight.
    subroutine bootstrap_histogram_impl(n_bootstraps, n_bins, n_points, n_studies, mean_pmf_counts, &
                                        mean_pmf_included_n_reps, included_n_reps, n_bootstrapping_top_k_jsds, &
                                        confidence_interval, tmp_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, &
                                        tmp_mean_pmf, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                        two_sided_bootstrapping_significance_level, random_seed, ierr)
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_js_comp_test_impl, AUTO)
            !! DM_MIN(1_int32)
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array used as top-k/bottom-k heaps for efficient percentile detection in
            !! the bootstrapped values -- `(:, 1, :)` the bottom-k (lower bound), `(:, 2, :)` the
            !! top-k (upper bound)
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_counts
            !! Working array that holds one bootstrap's resampled histogram counts, one study at a time
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: tmp_pmfs
            !! Working array that holds one bootstrap's resampled pmfs, all studies
        real(real64), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf
            !! Working array for one bootstrap's own (resampled) mean pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Working array for one bootstrap's per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Working array for one bootstrap's per-point global-JSD weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array for one bootstrap's global weighted JSD values
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
            !! otherwise used here
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        integer(int32) :: i_bootstrap, i_study, i_point, i_bin, draw_ierr
        type(rng_t) :: rng

        call set_ok(ierr)
        rng = create_rng(ierr, random_seed)
        if (is_err(ierr)) then
            ! A genuine runtime failure -- leave every work array in a defined (zeroed) state and
            ! the caller's confidence_interval untouched, rather than bootstrap with an
            ! uninitialized generator.
            tmp_bootstrapping_top_k_jsds = 0.0_real64
            tmp_counts = 0_int32
            tmp_pmfs = 0.0_real64
            tmp_mean_pmf = 0.0_real64
            tmp_js_divergences = 0.0_real64
            tmp_weights = 0.0_real64
            tmp_global_js_divergence = 0.0_real64
            return
        end if

        do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, confidence_interval)
            tmp_bootstrapping_top_k_jsds(:, 1, i_study) = confidence_interval(1, i_study)
            tmp_bootstrapping_top_k_jsds(:, 2, i_study) = confidence_interval(2, i_study)
        end do

        do i_bootstrap = 1, n_bootstraps
            ! Resample: each reference point's pooled/consensus histogram is a pool with
            ! replacement, drawn from separately for each study -- sequential over i_study/i_point,
            ! since random_multinomial mutates the (impure) rng state.
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multinomial(rng, n_bins, mean_pmf_counts(:, i_point), mean_pmf_included_n_reps(i_point), &
                                            included_n_reps(i_point, i_study), tmp_counts(:, i_point), draw_ierr)
                    if (is_err(draw_ierr)) call set_err_once(ierr, get_err_code(draw_ierr))
                end do

                do concurrent(i_point=1:n_points, i_bin=1:n_bins) shared(tmp_pmfs, tmp_counts, included_n_reps, i_study)
                    if (included_n_reps(i_point, i_study) == 0_int32) then
                        tmp_pmfs(i_bin, i_point, i_study) = 0.0_real64
                    else
                        tmp_pmfs(i_bin, i_point, i_study) = real(tmp_counts(i_bin, i_point), real64) &
                                                            /real(included_n_reps(i_point, i_study), real64)
                    end if
                end do
            end do

            ! As resampling was with replacement, the mean pmf changed -- very likely.
            call create_mean_pmf_only_impl(tmp_pmfs, n_bins, n_points, n_studies, tmp_mean_pmf)

            do concurrent(i_study=1:n_studies) shared(tmp_pmfs, tmp_mean_pmf, n_points, n_bins, tmp_js_divergences, &
                                                      included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, &
                                                      tmp_weights, tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds)
                call compute_divergence_per_reference_point_impl(transpose(tmp_pmfs(:, :, i_study)), transpose(tmp_mean_pmf), &
                                                                 n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_impl(tmp_js_divergences(:, i_study), n_points, &
                                                              included_n_reps(:, i_study), mean_pmf_included_n_reps, &
                                                              tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                ! Keep the highest and lowest values seen so far.
                call bottom_k_heap_push(tmp_bootstrapping_top_k_jsds(:, 1, i_study), n_bootstrapping_top_k_jsds, &
                                        tmp_global_js_divergence(i_study))
                call top_k_heap_push(tmp_bootstrapping_top_k_jsds(:, 2, i_study), n_bootstrapping_top_k_jsds, &
                                     tmp_global_js_divergence(i_study))
            end do
        end do

        ! Set confidence interval -> [largest small value, smallest large value]
        do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, confidence_interval)
            confidence_interval(1, i_study) = tmp_bootstrapping_top_k_jsds(1, 1, i_study)
            confidence_interval(2, i_study) = tmp_bootstrapping_top_k_jsds(1, 2, i_study)
        end do

        call destroy_rng(rng)
    end subroutine bootstrap_histogram_impl

    !> summary: Push a value onto a top-k heap (a min-heap tracking the k largest values pushed)
    !| AUTHOR_LASZLO_LANG
    !| Private module-internal helper for
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]],
    !| ported from 125-stabilize-jscomp's `f42_heaps.F90` (`top_k_heap_push`) -- kept
    !| module-internal here rather than promoted to its own module, per the plan for this port.
    !| Unlike its upstream counterpart's usual pairing with `init_top_k_heap` (which fills `heap`
    !| with `-Inf`), `bootstrap_histogram_impl` seeds every slot with its own initial reference
    !| value instead, so this routine never needs (and does not provide) that initialization.
    pure subroutine top_k_heap_push(heap, heap_size, value)
        integer(int32), intent(in) :: heap_size
            !! Size of heap
        real(real64), intent(inout) :: heap(heap_size)
            !! The top-k heap array (min-heap of the k largest values pushed so far)
        real(real64), intent(in) :: value
            !! Value to push

        if (heap_size <= 0_int32) return

        if (value > heap(1)) then
            heap(1) = value
            call minheap_sift_down(heap, heap_size)
        end if
    end subroutine top_k_heap_push

    !> summary: Push a value onto a bottom-k heap (a max-heap tracking the k smallest values pushed)
    !| AUTHOR_LASZLO_LANG
    !| Private module-internal helper for
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]],
    !| ported from 125-stabilize-jscomp's `f42_heaps.F90` (`bottom_k_heap_push`) -- see
    !| [[tox_data_integration_js_comp_test_impl(module):top_k_heap_push(interface)]]'s doc comment
    !| for why no separate initialization routine is ported alongside it.
    pure subroutine bottom_k_heap_push(heap, heap_size, value)
        integer(int32), intent(in) :: heap_size
            !! Size of heap
        real(real64), intent(inout) :: heap(heap_size)
            !! The bottom-k heap array (max-heap of the k smallest values pushed so far)
        real(real64), intent(in) :: value
            !! Value to push

        if (heap_size <= 0_int32) return

        if (value < heap(1)) then
            heap(1) = value
            call maxheap_sift_down(heap, heap_size)
        end if
    end subroutine bottom_k_heap_push

    !> Top-down sift for a min-heap: swaps the root with its smaller child until the heap
    !| condition holds. Not published: a private helper of
    !| [[tox_data_integration_js_comp_test_impl(module):top_k_heap_push(interface)]], ported from
    !| 125-stabilize-jscomp's `f42_heaps.F90` (`minheap_sift_down`).
    pure subroutine minheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
            !! Current heap size
        real(real64), intent(inout) :: heap(n)
            !! The min-heap array

        integer(int32) :: i, left, right, smallest
        real(real64) :: temp

        i = 1_int32
        do
            left = 2_int32*i
            right = left + 1_int32
            if (left > n) exit

            smallest = left
            if (right <= n) then
                if (heap(right) < heap(smallest)) smallest = right
            end if

            if (heap(smallest) < heap(i)) then
                temp = heap(i)
                heap(i) = heap(smallest)
                heap(smallest) = temp
                i = smallest
            else
                exit
            end if
        end do
    end subroutine minheap_sift_down

    !> Top-down sift for a max-heap: swaps the root with its larger child until the heap
    !| condition holds. Not published: a private helper of
    !| [[tox_data_integration_js_comp_test_impl(module):bottom_k_heap_push(interface)]], ported
    !| from 125-stabilize-jscomp's `f42_heaps.F90` (`maxheap_sift_down`).
    pure subroutine maxheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
            !! Current heap size
        real(real64), intent(inout) :: heap(n)
            !! The max-heap array

        integer(int32) :: i, left, right, largest
        real(real64) :: temp

        i = 1_int32
        do
            left = 2_int32*i
            right = left + 1_int32
            if (left > n) exit

            largest = left
            if (right <= n) then
                if (heap(right) > heap(largest)) largest = right
            end if

            if (heap(largest) > heap(i)) then
                temp = heap(i)
                heap(i) = heap(largest)
                heap(largest) = temp
                i = largest
            else
                exit
            end if
        end do
    end subroutine maxheap_sift_down

    !> The part of run_js_comp_test_impl that runs after every study's histograms are built:
    !| the consensus pmf (create_mean_pmf_impl), each study's JSD against it and the weighted
    !| global JSD, the permutation test (gjct_permutation_test_impl), and the final re-derivation
    !| of each study's pmf/JSD/weights/global JSD from its untouched `counts`. Moved here unchanged
    !| from run_js_comp_test_impl, so a later driver that builds its neighborhoods differently can
    !| share it. Private, and not `pure`: gjct_permutation_test_impl is impure. If that routine
    !| fails, its error code is returned in `ierr` and this routine returns at once.
    subroutine finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, mean_pmf, &
                                   mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                   global_js_divergence, p_values, p_values_observed_consensus, &
                                   tmp_counts_point_major, tmp_pmf_point_major, tmp_permutation_mean_pmf_counts, &
                                   tmp_permutation_counts, tmp_permutation_pmfs, &
                                   tmp_permutation_js_divergences, tmp_permutation_weights, &
                                   tmp_permutation_global_js_divergence, n_permutations, random_seed, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: max_n_bins_per_point
            !! Widest bin count of any reference point; only rows `1:max_n_bins_per_point` of the
            !! bin-sized arrays below are read or written
        real(real64), dimension(MAX_N_BINS, n_points, n_studies), intent(inout) :: pmfs
            !! Per-study pmfs; read for the consensus, then re-derived from `counts`
        integer(int32), dimension(MAX_N_BINS, n_points, n_studies), intent(in) :: counts
            !! Per-study histogram counts
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(MAX_N_BINS, n_points), intent(out) :: mean_pmf
            !! The consensus pmf
        integer(int32), dimension(MAX_N_BINS, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts per bin for the consensus pmf
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study from gjct_permutation_test_impl, with each permuted
            !! study compared against the observed consensus `mean_pmf`
        integer(int32), dimension(n_points, MAX_N_BINS), intent(out) :: tmp_counts_point_major
            !! Working array forwarded to gjct_permutation_test_impl
        real(real64), dimension(n_points, MAX_N_BINS), intent(out) :: tmp_pmf_point_major
            !! Working array forwarded to gjct_permutation_test_impl and reused for calc_pmf_impl
        integer(int32), dimension(MAX_N_BINS, n_points), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to gjct_permutation_test_impl
        integer(int32), dimension(MAX_N_BINS, n_points), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to gjct_permutation_test_impl
        real(real64), dimension(MAX_N_BINS, n_points, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to gjct_permutation_test_impl
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array forwarded to gjct_permutation_test_impl
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array forwarded to gjct_permutation_test_impl
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array forwarded to gjct_permutation_test_impl
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations, already resolved from its default by the caller
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator, forwarded as-is
        integer(int32), intent(out) :: ierr
            !! Error code; the error code of gjct_permutation_test_impl if it fails (its argument
            !! position cleared), in which case this routine returns at once and the re-derived
            !! pmfs/JSDs/weights/global JSD are not computed

        integer(int32) :: i_study, permutation_ierr

        call set_ok(ierr)

        call create_mean_pmf_impl(pmfs(1:max_n_bins_per_point, :, :), counts(1:max_n_bins_per_point, :, :), max_n_bins_per_point, n_points, n_studies, &
                                  included_n_reps, mean_pmf(1:max_n_bins_per_point, :), mean_pmf_included_n_reps, &
                                  mean_pmf_counts(1:max_n_bins_per_point, :))

        do concurrent(i_study=1:n_studies) shared(pmfs, mean_pmf, n_points, max_n_bins_per_point, js_divergences, included_n_reps, &
                                                  mean_pmf_included_n_reps, global_js_divergence, weights)
            call compute_divergence_per_reference_point_impl(transpose(pmfs(1:max_n_bins_per_point, :, i_study)), &
                                                             transpose(mean_pmf(1:max_n_bins_per_point, :)), n_points, &
                                                             max_n_bins_per_point, js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
        end do

        call gjct_permutation_test_impl(n_permutations, max_n_bins_per_point, n_points, n_studies, &
                                        mean_pmf_counts(1:max_n_bins_per_point, :), mean_pmf(1:max_n_bins_per_point, :), &
                                        mean_pmf_included_n_reps, included_n_reps, global_js_divergence, p_values, &
                                        p_values_observed_consensus, &
                                        tmp_permutation_mean_pmf_counts(1:max_n_bins_per_point, :), &
                                        tmp_permutation_counts(1:max_n_bins_per_point, :), tmp_permutation_pmfs(1:max_n_bins_per_point, :, :), &
                                        tmp_permutation_js_divergences, tmp_permutation_weights, &
                                        tmp_permutation_global_js_divergence, tmp_pmf_point_major(:, 1:max_n_bins_per_point), &
                                        tmp_counts_point_major(:, 1:max_n_bins_per_point), random_seed, permutation_ierr)
        if (is_err(permutation_ierr)) then
            call set_err_once(ierr, get_err_code(permutation_ierr))
            return
        end if

        ! Re-derive each study's own pmf/JSD/weights/global JSD from its UNTOUCHED counts --
        ! mean_pmf/mean_pmf_counts are NOT re-derived, since the permutation test above only
        ! perturbs its own scratch copies (tmp_permutation_*), never mean_pmf_counts itself.
        do i_study = 1, n_studies
            call calc_pmf_impl(transpose(counts(1:max_n_bins_per_point, :, i_study)), included_n_reps(:, i_study), n_points, &
                               max_n_bins_per_point, tmp_pmf_point_major(:, 1:max_n_bins_per_point))
            pmfs(1:max_n_bins_per_point, :, i_study) = transpose(tmp_pmf_point_major(:, 1:max_n_bins_per_point))
            call compute_divergence_per_reference_point_impl(transpose(pmfs(1:max_n_bins_per_point, :, i_study)), &
                                                             transpose(mean_pmf(1:max_n_bins_per_point, :)), n_points, &
                                                             max_n_bins_per_point, js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
        end do
    end subroutine finish_js_comp_test

    !> summary: Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `js_comp_test_helper`, restructured for Issue #187's
    !| occupancy-constrained per-neighborhood histogram binning into three passes, mirroring
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]'s
    !| own Pass A/B/C split (Issue #187's own Steps 2.5/2.6):
    !|
    !| - Pass A (per study): builds every study's neighborhoods
    !|   ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]),
    !|   writing into `neighborhood_indices`/`neighborhood_range`, which already retain every
    !|   study's own values simultaneously (both are real `intent(out)` arguments sized
    !|   `(..., n_points, n_studies)` -- unlike `run_js_comp_test_parameter_search_impl`, no new
    !|   buffer was needed for this). Unlike that routine, there is no admissibility gate here, so
    !|   Pass A always runs to completion for every study.
    !| - Pass B (per point, sequential -- see the implementation body's own comment for why): pools
    !|   every study's residuals for one reference point at a time
    !|   ([[tox_data_integration_js_comp_test_impl(module):gather_pooled_neighborhood_residuals(interface)]],
    !|   now published as its own entry point) and runs Issue #187's occupancy-constrained
    !|   bin-count search on the pooled result
    !|   ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]),
    !|   deciding `n_bins_per_point(i_point)` independently for every reference point, plus the
    !|   `occupancy_failed`/`n_pooled_residuals`/`min_bin_occupancy`/`mean_bin_occupancy`/
    !|   `max_bin_occupancy`/`sturges_bins`/`fd_bins` diagnostics. `max_n_bins_per_point`
    !|   (`maxval(n_bins_per_point(1:n_points))`) is derived once after Pass B and replaces the old
    !|   caller-supplied scalar `n_bins` everywhere downstream.
    !| - Pass C (per study): re-gathers this study's residual values from the neighbor indices Pass
    !|   A already computed, then builds its residual histograms at the real per-point bin counts
    !|   ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !|
    !| `n_neighbors` may not exceed `max_n_genes_all_studies`, so every neighbor gene index Pass A
    !| produces lies in `[1, max_n_genes_all_studies]`. Pass B still checks each index and, should
    !| one ever fall outside, reports `ERR_INVALID_INPUT` and returns before Pass C could read
    !| `residuals` out of bounds -- a safeguard that valid input cannot reach.
    !|
    !| After Pass C, the pipeline continues exactly as before: pools the per-study pmfs into the
    !| consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), computes
    !| each study's observed JSD against that consensus
    !| ([[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]/[[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl(interface)]],
    !| called with the consensus pmf as the second argument), runs the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]) -- which
    !| returns the primary `p_values` (each permuted study against the consensus of the permuted
    !| studies) and, from the same permutations, `p_values_observed_consensus` (each permuted study
    !| against the observed consensus `mean_pmf`, kept for comparison) -- and
    !| finally re-derives each study's pmf/JSD/weights/global JSD from its own UNTOUCHED `counts`
    !| via
    !| [[tox_data_integration_jsd_impl(module):calc_pmf_impl(interface)]] -- `mean_pmf`/`mean_pmf_counts`
    !| are NOT re-derived, since the permutation test never modifies the observed consensus (it
    !| resamples only its own scratch copies, never `mean_pmf_counts` itself, and builds each
    !| permutation's consensus in its own scratch too), exactly as 125 relies on.
    !|
    !| **Behavioral asymmetry vs.
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| -- read before using this entry point where inadequately-supported neighborhoods must be
    !| rejected:** unlike that routine, THIS one has NO multi-candidate fallback and NO
    !| admissibility gate at all (no
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]],
    !| no [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]],
    !| no early exit). A reference point whose Pass B occupancy search fails even at `m_min`
    !| (`occupancy_failed(i_point) = .true._c_bool`) still gets a real histogram built, at
    !| `n_bins_per_point(i_point) == m_min`, and that point still contributes to
    !| `global_js_divergence` exactly like every other point -- its contribution is down-weighted
    !| only by `included_n_reps` (an orthogonal quantity: how many non-NaN replicates it has), never
    !| by bin sparsity. A caller that needs inadequately-supported neighborhoods rejected outright
    !| should use `run_js_comp_test_parameter_search_impl` instead, which gates on exactly this via
    !| `check_mean_pmf_min_counts_impl`.
    !|
    !| **A real, deliberate change to this routine's public array shapes (Issue #187):** the old
    !| mandatory scalar input `n_bins` is gone -- there is no way for a caller to know the right bin
    !| count in advance, since it is now genuinely computed inside this routine by Pass B's
    !| occupancy search, independently per reference point. Every array whose bin-sized dimension
    !| used to be sized by that input (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`,
    !| `tmp_counts_point_major`, `tmp_pmf_point_major`, `tmp_permutation_mean_pmf_counts`,
    !| `tmp_permutation_counts`, `tmp_permutation_pmfs`) is now sized to the fixed compile-time
    !| ceiling [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] (`256`)
    !| instead, exactly mirroring how `run_js_comp_test_parameter_search_impl`'s own
    !| `tmp_counts_point_major`/`tmp_pmf_point_major`/`tmp_pmfs`/`tmp_counts` etc. have been sized
    !| since Issue #187's earlier steps. The new `max_n_bins_per_point` output tells a caller how many of the
    !| LEADING bins/rows of each of those arrays are actually meaningful
    !| (`maxval(n_bins_per_point(1:n_points))`); the rest is unused padding. The generator's own
    !| result-size trimming directive cannot express this trim, because it only ever trims an
    !| array's LAST declared extent, and bins is the FIRST declared extent of every one of those
    !| arrays -- so a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves, exactly as a
    !| caller of `run_js_comp_test_parameter_search_impl`'s own jagged `trace_*` arrays already has
    !| to.
    !|
    !| `x_star` is an ordinary input here, not computed by this routine -- 125's own
    !| `js_comp_test_helper` takes it the same way, since a caller running several studies/several
    !| parameter settings is expected to compute the reference points once
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and reuse
    !| them consistently.
    !|
    !| `construct_neighborhoods_ranged_impl` reports neighbor gene INDICES, not gathered residual
    !| values (unlike its distance-sort sibling
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl(interface)]]),
    !| so Pass C gathers each neighbor's actual residual values from `residuals` itself
    !| (one slice of `tmp_neighborhood_residuals_gathered` per reference point, reused for every
    !| study) before building that point's histogram with `build_residual_histograms_impl` for
    !| that single point, whose bins are copied straight into the matching column of
    !| `counts`/`pmfs` -- no transpose is needed there. `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts`
    !| here are BIN-major (`(256, n_points, n_studies)`) to match
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]'s own
    !| convention, while `compute_divergence_per_reference_point_impl` and `calc_pmf_impl` are
    !| POINT-major (`(n_points, max_n_bins_per_point)`): the calls to those two after Pass C bridge
    !| with an explicit `transpose`, exactly as
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] and
    !| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] already do.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is returned
    !| in `ierr`, and the routine returns right there: `pmfs`, `js_divergences`, `weights` and
    !| `global_js_divergence` then hold the pre-permutation values, not the final re-derived ones.
    subroutine run_js_comp_test_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                                     gene_means, gene_means_perms, residuals, x_star, &
                                     neighborhood_indices, neighborhood_range, n_bins_per_point, &
                                     shared_residual_range_low, shared_residual_range_high, max_n_bins_per_point, &
                                     occupancy_failed, n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                     max_bin_occupancy, sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, &
                                     mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                     global_js_divergence, p_values, p_values_observed_consensus, &
                                     tmp_neighborhood_residuals_gathered, &
                                     tmp_counts_point_major, tmp_pmf_point_major, tmp_pooled_residuals, &
                                     tmp_pooled_residuals_perm, tmp_bin_counts_search, tmp_point_n_neighbors, &
                                     tmp_permutation_mean_pmf_counts, tmp_permutation_counts, tmp_permutation_pmfs, &
                                     tmp_permutation_js_divergences, tmp_permutation_weights, &
                                     tmp_permutation_global_js_divergence, n_permutations, &
                                     random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, &
                                     lower_residual_range_quantile, upper_residual_range_quantile, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! DM_ALLOW_NAN
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! DM_ALLOW_NAN
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
            !! DM_ALLOW_NAN
        integer(int32), dimension(n_neighbors, n_points, n_studies), intent(out) :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study (Pass A)
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            !! produced by construct_neighborhoods_ranged_impl (Pass A)
        integer(int32), dimension(n_points), intent(out) :: n_bins_per_point
            !! This reference point's own selected histogram bin count (Issue #187's `M_j`), from
            !! Pass B's occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood
            !! may use a different bin count
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_low
            !! This reference point's own lower residual-range bound (R_low), from Pass B's
            !! occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            !! different, asymmetric range (Step 3)
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_high
            !! This reference point's own upper residual-range bound (R_high), from Pass B's
            !! occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            !! different, asymmetric range (Step 3)
        integer(int32), intent(out) :: max_n_bins_per_point
            !! The widest `n_bins_per_point` value across all `n_points` reference points
            !! (`maxval(n_bins_per_point(1:n_points))`), derived once after Pass B. The number of
            !! leading, meaningful bins/rows in `pmfs`, `counts`, `mean_pmf` and `mean_pmf_counts`
            !! below (and of the bin-sized work arrays) -- those are all declared
            !! with a fixed 256-bin ceiling (MAX_N_BINS) rather than a caller-supplied bin count,
            !! since `n_bins_per_point` can no longer be known by a caller in advance. A Python/R
            !! caller must slice `[:max_n_bins_per_point, ...]` themselves: the generator's own result-size
            !! trimming directive cannot express this trim, because it only ever trims an array's
            !! LAST declared extent, and bins is the FIRST declared extent of every one of those
            !! arrays
        logical(c_bool), dimension(n_points), intent(out) :: occupancy_failed
            !! This reference point's `occupancy_failed` flag from Pass B
            !! (determine_bin_count_occupancy_impl) -- `.true.` iff even `m_min` bins could not
            !! satisfy the occupancy criterion for it. See this routine's own doc block above for
            !! the behavioral asymmetry this implies vs. run_js_comp_test_parameter_search_impl: a
            !! `.true.` point here still gets a real histogram and still contributes to
            !! `global_js_divergence`, it is never rejected
        integer(int32), dimension(n_points), intent(out) :: n_pooled_residuals
            !! This reference point's pooled residual count (N_j) from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(int32), dimension(n_points), intent(out) :: min_bin_occupancy
            !! This reference point's minimum bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        real(real64), dimension(n_points), intent(out) :: mean_bin_occupancy
            !! This reference point's mean bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(int32), dimension(n_points), intent(out) :: max_bin_occupancy
            !! This reference point's maximum bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(int32), dimension(n_points), intent(out) :: sturges_bins
            !! This reference point's Sturges' rule bin-count diagnostic from Pass B
            !! (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            !! decision
        integer(int32), dimension(n_points), intent(out) :: fd_bins
            !! This reference point's Freedman-Diaconis rule bin-count diagnostic from Pass B
            !! (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            !! decision
        real(real64), dimension(256, n_points, n_studies), intent(out) :: pmfs
            !! `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`. `256`
            !! = MAX_N_BINS, a fixed ceiling (see `max_n_bins_per_point` above) -- only rows `1:max_n_bins_per_point`
            !! are meaningful; a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves
        integer(int32), dimension(256, n_points, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
        real(real64), dimension(256, n_points), intent(out) :: mean_pmf
            !! The consensus pmf, from create_mean_pmf_impl. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(int32), dimension(256, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS;
            !! only rows `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl, each permuted study
            !! compared against the consensus of the permuted studies -- the primary p-value
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study from the same permutations, but with each permuted study
            !! compared against the fixed observed consensus `mean_pmf`, as the permutation test did
            !! before it recomputed the consensus per permutation. Kept for comparison only;
            !! `p_values` is the primary result
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: &
            tmp_neighborhood_residuals_gathered
            !! Working array: gathered neighborhood residual values, one slice per reference
            !! point, reused for every study (Pass C)
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_counts_point_major
            !! Working array forwarded to gjct_permutation_test_impl's own point-major counts
            !! scratch. `256` = MAX_N_BINS, written as a literal because a generated wrapper's
            !! dummy dimension cannot reference a module parameter
        real(real64), dimension(n_points, 256), intent(out) :: tmp_pmf_point_major
            !! Working array: forwarded to gjct_permutation_test_impl's own point-major pmf
            !! scratch, then reused for calc_pmf_impl's re-derived per-study pmf. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        real(real64), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals across every neighbor and
            !! every study (Pass B), reused per point -- one small buffer, not one per point, since
            !! Pass B is a deliberate sequential loop (see the implementation body's own comment)
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array forwarded to determine_bin_count_occupancy_impl's own per-bin-count
            !! search scratch, reused per point. `256` = MAX_N_BINS, matching
            !! determine_bin_count_occupancy_impl's own tmp_bin_counts dummy -- written as a literal
            !! because a generated wrapper's dummy dimension cannot reference a module parameter
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array: every study's neighbor count for the reference point in Pass B --
            !! `n_neighbors` for every study, since this routine uses one neighbor count throughout
        integer(int32), dimension(256, n_points), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to gjct_permutation_test_impl's own resampling pool. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(int32), dimension(256, n_points), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to gjct_permutation_test_impl's own per-study resampled
            !! counts. `256` = MAX_N_BINS, see tmp_counts_point_major above
        real(real64), dimension(256, n_points, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to gjct_permutation_test_impl's own resampled pmfs. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array forwarded to gjct_permutation_test_impl's own per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array forwarded to gjct_permutation_test_impl's own per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array forwarded to gjct_permutation_test_impl's own resampled global JSD values
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(1000_int32)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible in Pass B's occupancy search, forwarded to
            !! determine_bin_count_occupancy_impl
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT if a neighbor gene index is out of range (a safeguard
            !! valid input cannot reach, see above), ERR_ALLOC_FAIL if GSL could not allocate the random number generator for the
            !! permutation test

        integer(int32) :: i_study, i_point, actual_n_permutations, point_ierr, finish_ierr

        call set_ok(ierr)
        M_DEFAULT_VAL(n_permutations, actual_n_permutations, 1000_int32)

        ! ===== PASS A (per study): neighborhood construction only. Unlike
        ! run_js_comp_test_parameter_search_impl, this routine has no admissibility gate -- it
        ! always runs to completion for every study. The residual-gather and histogram-building
        ! tail that used to live in this loop moves to Pass C below, after Pass B has decided real
        ! per-point bin counts.
        do i_study = 1, n_studies
            call construct_neighborhoods_ranged_impl(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), &
                                                     gene_means_perms(:, i_study), n_neighbors, &
                                                     neighborhood_indices(:, :, i_study), neighborhood_range(:, :, i_study))
        end do

        ! ===== PASS B (per point): Issue #187's occupancy-constrained bin-count search, pooling
        ! every study's residuals for one reference point at a time (determine_point_bin_count).
        ! Deliberately a plain sequential `do`, NOT `do concurrent`: a single small
        ! tmp_pooled_residuals buffer is reused across points here rather than allocating one such
        ! buffer per point (up to n_points of them, which `do concurrent`'s `local()` would
        ! require) -- the memory-conscious choice given each iteration's work (gather, heapsort,
        ! one occupancy search) is already substantial on its own, mirroring
        ! run_js_comp_test_parameter_search_impl's own Pass B. Every study has the same
        ! n_neighbors here.
        tmp_point_n_neighbors = n_neighbors
        do i_point = 1, n_points
            call determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                           n_neighbors, neighborhood_indices(:, i_point, :), tmp_point_n_neighbors, &
                                           n_bins_per_point(i_point), occupancy_failed(i_point), &
                                           shared_residual_range_low(i_point), shared_residual_range_high(i_point), &
                                           n_pooled_residuals(i_point), min_bin_occupancy(i_point), &
                                           mean_bin_occupancy(i_point), max_bin_occupancy(i_point), &
                                           sturges_bins(i_point), fd_bins(i_point), tmp_pooled_residuals, &
                                           tmp_pooled_residuals_perm, tmp_bin_counts_search, point_ierr, &
                                           m_min=m_min, m_max=m_max, min_residuals_per_bin=min_residuals_per_bin, &
                                           gamma_occupancy=gamma_occupancy, &
                                           lower_residual_range_quantile=lower_residual_range_quantile, &
                                           upper_residual_range_quantile=upper_residual_range_quantile)
            ! Defensive: determine_point_bin_count only fails on an out-of-range gene index. Pass A
            ! is called with n_genes_S = max_n_genes_all_studies, the wrapper bounds n_neighbors by
            ! that and gene_means_perms to [1, max_n_genes_all_studies], so every index is in range
            ! and this branch cannot be reached through the wrapper. It is kept for direct callers,
            ! since Pass C would otherwise read residuals out of bounds, and has no test.
            if (is_err(point_ierr)) then
                call set_err_once(ierr, get_err_code(point_ierr))
                return
            end if
        end do

        ! ===== PASS C (per study, per point): build each (point, study) histogram from the
        ! indices Pass A already computed, now that Pass B has decided real per-point bin counts
        ! (build_point_study_histogram). max_n_bins_per_point (the widest bin count any point uses)
        ! replaces the old caller-supplied scalar n_bins for every remaining use below. Points are
        ! independent, and each gathers into its own slice of tmp_neighborhood_residuals_gathered.
        max_n_bins_per_point = maxval(n_bins_per_point(1:n_points))
        do i_study = 1, n_studies
            do concurrent(i_point=1:n_points) shared(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                     n_neighbors, neighborhood_indices, shared_residual_range_low, &
                                                     shared_residual_range_high, n_bins_per_point, max_n_bins_per_point, &
                                                     counts, pmfs, included_n_reps, &
                                                     tmp_neighborhood_residuals_gathered, i_study)
                call build_point_study_histogram(residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, &
                                                 n_neighbors, neighborhood_indices(:, i_point, i_study), &
                                                 shared_residual_range_low(i_point), shared_residual_range_high(i_point), &
                                                 n_bins_per_point(i_point), max_n_bins_per_point, &
                                                 counts(1:max_n_bins_per_point, i_point, i_study), &
                                                 pmfs(1:max_n_bins_per_point, i_point, i_study), &
                                                 included_n_reps(i_point, i_study), &
                                                 tmp_neighborhood_residuals_gathered(:, :, i_point))
            end do
        end do

        call finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, mean_pmf, &
                                 mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                 global_js_divergence, p_values, p_values_observed_consensus, &
                                 tmp_counts_point_major, tmp_pmf_point_major, tmp_permutation_mean_pmf_counts, &
                                 tmp_permutation_counts, tmp_permutation_pmfs, &
                                 tmp_permutation_js_divergences, tmp_permutation_weights, &
                                 tmp_permutation_global_js_divergence, actual_n_permutations, random_seed, &
                                 finish_ierr)
        if (is_err(finish_ierr)) then
            call set_err_once(ierr, get_err_code(finish_ierr))
            return
        end if
    end subroutine run_js_comp_test_impl

    !> summary: Run the JSD-Comp-Test pipeline on adaptive (Issue #217) neighborhoods
    !| AUTHOR_LASZLO_LANG
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]]: the
    !| same pipeline and the same outputs, on neighborhoods whose size varies per reference point
    !| and per study. Call
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| first, then pass its `pooled_neighborhood_range` (trimmed to its `n_points` reference
    !| points) here together with the same `gene_means` and `residuals`.
    !|
    !| **Membership.** The pooled gene means are sorted again exactly as the construction sorts them
    !| (the same routine on the same input), so a pooled position means the same entry here as it
    !| did there. Every range `a..b` is decoded back to per-study gene lists by
    !| [[tox_data_integration_js_comp_test_impl(module):materialize_pooled_neighborhood(subroutine)]];
    !| the per-study counts are returned as `n_neighbors_per_point`, in the construction's own
    !| `(n_studies, n_points)` orientation, and add up to `b - a + 1`.
    !|
    !| **Pipeline, per reference point.** Pass B pools the point's residuals across all its genes
    !| of all studies and picks its bin count and residual range with Issue #187's occupancy
    !| search
    !| ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]),
    !| its Sturges/Freedman-Diaconis diagnostics computed with the rounded mean per-study neighbor
    !| count. Pass C builds each study's histogram from that study's own genes of the point. Then,
    !| exactly as in `run_js_comp_test_impl`: the consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), each
    !| study's JSD against it and the weighted global JSD, the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]], with
    !| both its `p_values` and its comparison `p_values_observed_consensus`) and the
    !| final re-derivation from the untouched `counts`. Neighborhoods that decode to the same gene
    !| sets as a fixed-k run's therefore give bit-identical results. Every point is weighted by its
    !| non-NaN residual count, which under adaptive growth genuinely differs between points.
    !| Like `run_js_comp_test_impl`, there is no admissibility gate: overlap and occupancy are the
    !| parameter search's concern.
    !|
    !| **Runtime errors.** Before any histogram work, every range is checked in point order; the
    !| first failing one sets `ERR_INVALID_INPUT` and the routine returns at once, every output
    !| undefined. For each point, in this order: its first position lies after its last; its last
    !| position lies beyond the non-NaN pooled means (the entries a construction can ever use); a
    !| study has no gene in it, which would leave that study's pmf empty and the consensus pmf, an
    !| average over all studies, meaningless. A range a construction produced passes the first two
    !| checks by construction, and the third whenever its status was not an empty-study one.
    !|
    !| **Memory.** The work arrays are sized by the per-study gene bound `max_n_genes_all_studies`,
    !| the largest number of genes one study can have in one neighborhood, not by a caller-supplied
    !| neighbor count. The pooling buffers `tmp_pooled_residuals`/`tmp_pooled_residuals_perm`
    !| therefore hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies` values each,
    !| the size of `residuals` itself; Pass B and Pass C are sequential over points so one such
    !| buffer serves every point.
    !|
    !| The bin-sized arrays have the same fixed 256-bin
    !| ([[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]]) leading extent as in
    !| `run_js_comp_test_impl`, of which only the first `max_n_bins_per_point` rows are
    !| meaningful; a Python/R caller slices `[:max_n_bins_per_point, ...]` themselves.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is returned
    !| in `ierr`, and the routine returns right there: `pmfs`, `js_divergences`, `weights` and
    !| `global_js_divergence` then hold the pre-permutation values, not the final re-derived ones.
    subroutine run_js_comp_test_adaptive_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, &
                                              gene_means, residuals, pooled_neighborhood_range, n_neighbors_per_point, &
                                              n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                              max_n_bins_per_point, occupancy_failed, n_pooled_residuals, &
                                              min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, sturges_bins, &
                                              fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, &
                                              p_values, p_values_observed_consensus, tmp_gene_means_perm_all, &
                                              tmp_point_neighborhood_indices, &
                                              tmp_neighbor_residuals, tmp_counts_point_major, tmp_pmf_point_major, &
                                              tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                              tmp_permutation_mean_pmf_counts, tmp_permutation_counts, &
                                              tmp_permutation_pmfs, tmp_permutation_js_divergences, &
                                              tmp_permutation_weights, tmp_permutation_global_js_divergence, &
                                              n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, &
                                              gamma_occupancy, lower_residual_range_quantile, &
                                              upper_residual_range_quantile, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the same
            !! array the neighborhoods were constructed from
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! DM_ALLOW_NAN
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means, as returned by construct_adaptive_neighborhoods
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies*n_studies)
        integer(int32), dimension(n_studies, n_points), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(int32), dimension(n_points), intent(out) :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the occupancy search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the occupancy search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the occupancy search
        integer(int32), intent(out) :: max_n_bins_per_point
            !! The widest `n_bins_per_point` value, `maxval(n_bins_per_point)`: the number of leading,
            !! meaningful rows of `pmfs`, `counts`, `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out) :: occupancy_failed
            !! `.true.` iff even `m_min` bins could not satisfy the occupancy criterion for this
            !! reference point; such a point still gets a histogram at `m_min` bins and still
            !! contributes to `global_js_divergence`
        integer(int32), dimension(n_points), intent(out) :: n_pooled_residuals
            !! Each reference point's pooled non-NaN residual count (N_j), across all its genes of
            !! all studies
        integer(int32), dimension(n_points), intent(out) :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(real64), dimension(n_points), intent(out) :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic, with the rounded mean per-study
            !! neighbor count; never part of the occupancy search's decision
        integer(int32), dimension(n_points), intent(out) :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic, with the rounded mean
            !! per-study neighbor count; never part of the occupancy search's decision
        real(real64), dimension(256, n_points, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(256, n_points, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(256, n_points), intent(out) :: mean_pmf
            !! The consensus pmf. `256` = MAX_N_BINS; only rows `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(256, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS; only
            !! rows `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study from the permutation test, each permuted study compared
            !! against the consensus of the permuted studies -- the primary p-value
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study from the same permutations, but with each permuted study
            !! compared against the fixed observed consensus `mean_pmf`, as the permutation test did
            !! before it recomputed the consensus per permutation. Kept for comparison only;
            !! `p_values` is the primary result
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, seeded and sorted here
            !! exactly as construct_adaptive_neighborhoods sorts it
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused (Pass C)
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_counts_point_major
            !! Working array forwarded to the permutation test. `256` = MAX_N_BINS
        real(real64), dimension(n_points, 256), intent(out) :: tmp_pmf_point_major
            !! Working array forwarded to the permutation test and reused for the re-derived pmfs.
            !! `256` = MAX_N_BINS
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals across all its genes of all
            !! studies (Pass B), reused per point; as large as `residuals`
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array forwarded to the occupancy search, reused per point. `256` = MAX_N_BINS
        integer(int32), dimension(256, n_points), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to the permutation test. `256` = MAX_N_BINS
        integer(int32), dimension(256, n_points), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to the permutation test. `256` = MAX_N_BINS
        real(real64), dimension(256, n_points, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to the permutation test. `256` = MAX_N_BINS
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array forwarded to the permutation test
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(1000_int32)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible in the occupancy search
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the occupancy search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the occupancy search tests (M_max); raised to `m_min`
            !! when smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search's coarse stage
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's lower residual-range bound
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's upper residual-range bound
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for a range that is reversed, reaches past the non-NaN
            !! pooled means or leaves a study without a gene (see above), ERR_ALLOC_FAIL if GSL
            !! could not allocate the random number generator for the permutation test

        integer(int32) :: i_point, i_study, n_valid_means, n_study_neighbors, actual_n_permutations, point_ierr, &
                          finish_ierr

        call set_ok(ierr)
        M_DEFAULT_VAL(n_permutations, actual_n_permutations, 1000_int32)

        ! The construction's own sort, so every pooled position names the same entry as there.
        call init_perm(tmp_gene_means_perm_all)
        call sort_real_heapsort_expl_size(gene_means, tmp_gene_means_perm_all, max_n_genes_all_studies*n_studies)
        ! NaN sorts last, so the non-NaN means occupy exactly the first n_valid_means positions.
        n_valid_means = count(.not. ieee_is_nan(gene_means), kind=int32)

        ! Runtime checks, all before any histogram work. Sequential with an early return: the
        ! first failing point ends the routine.
        do i_point = 1, n_points
            if (pooled_neighborhood_range(1, i_point) > pooled_neighborhood_range(2, i_point) .or. &
                pooled_neighborhood_range(2, i_point) > n_valid_means) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
                return
            end if
            call materialize_pooled_neighborhood(tmp_gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                 pooled_neighborhood_range(1, i_point), &
                                                 pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                 tmp_point_neighborhood_indices, n_neighbors_per_point(:, i_point))
            if (any(n_neighbors_per_point(:, i_point) == 0_int32)) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
                return
            end if
        end do

        ! ===== PASS B (per point, sequential: one pooling buffer serves every point).
        do i_point = 1, n_points
            call materialize_pooled_neighborhood(tmp_gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                 pooled_neighborhood_range(1, i_point), &
                                                 pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                 tmp_point_neighborhood_indices, n_neighbors_per_point(:, i_point))
            call determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                           max_n_genes_all_studies, tmp_point_neighborhood_indices, &
                                           n_neighbors_per_point(:, i_point), n_bins_per_point(i_point), &
                                           occupancy_failed(i_point), shared_residual_range_low(i_point), &
                                           shared_residual_range_high(i_point), n_pooled_residuals(i_point), &
                                           min_bin_occupancy(i_point), mean_bin_occupancy(i_point), &
                                           max_bin_occupancy(i_point), sturges_bins(i_point), fd_bins(i_point), &
                                           tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                           point_ierr, m_min=m_min, m_max=m_max, &
                                           min_residuals_per_bin=min_residuals_per_bin, gamma_occupancy=gamma_occupancy, &
                                           lower_residual_range_quantile=lower_residual_range_quantile, &
                                           upper_residual_range_quantile=upper_residual_range_quantile)
            ! Defensive: determine_point_bin_count only fails on an out-of-range gene index, and
            ! every decoded index lies in [1, max_n_genes_all_studies], so this branch cannot be
            ! reached under the checks above and has no test.
            if (is_err(point_ierr)) then
                call set_err_once(ierr, get_err_code(point_ierr))
                return
            end if
        end do

        ! ===== PASS C (per point, per study, sequential: one gather buffer serves every pair).
        max_n_bins_per_point = maxval(n_bins_per_point)
        do i_point = 1, n_points
            call materialize_pooled_neighborhood(tmp_gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                 pooled_neighborhood_range(1, i_point), &
                                                 pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                 tmp_point_neighborhood_indices, n_neighbors_per_point(:, i_point))
            do i_study = 1, n_studies
                n_study_neighbors = n_neighbors_per_point(i_study, i_point)
                call build_point_study_histogram(residuals(:, :, i_study), max_n_reps_all_studies, &
                                                 max_n_genes_all_studies, n_study_neighbors, &
                                                 tmp_point_neighborhood_indices(1:n_study_neighbors, i_study), &
                                                 shared_residual_range_low(i_point), shared_residual_range_high(i_point), &
                                                 n_bins_per_point(i_point), max_n_bins_per_point, &
                                                 counts(1:max_n_bins_per_point, i_point, i_study), &
                                                 pmfs(1:max_n_bins_per_point, i_point, i_study), &
                                                 included_n_reps(i_point, i_study), &
                                                 tmp_neighbor_residuals(:, 1:n_study_neighbors))
            end do
        end do

        call finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, mean_pmf, &
                                 mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                 global_js_divergence, p_values, p_values_observed_consensus, &
                                 tmp_counts_point_major, tmp_pmf_point_major, tmp_permutation_mean_pmf_counts, &
                                 tmp_permutation_counts, tmp_permutation_pmfs, &
                                 tmp_permutation_js_divergences, tmp_permutation_weights, &
                                 tmp_permutation_global_js_divergence, actual_n_permutations, random_seed, &
                                 finish_ierr)
        if (is_err(finish_ierr)) then
            call set_err_once(ierr, get_err_code(finish_ierr))
            return
        end if
    end subroutine run_js_comp_test_adaptive_impl

    !> summary: Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `determine_js_comp_test_n_points_n_neighbors_helper` and
    !| `_alloc`, merged into one implementation now that the new `_impl` rules leave no separate
    !| hand-written allocation layer. Pools all studies' gene means, sorts them once
    !| ([[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]), generates the candidate
    !| grid from the gene count alone
    !| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| -- Issue #187, Step 2.7: this no longer needs the pooled residuals, since real per-neighborhood
    !| bin counts are decided later, in Pass B below),
    !| then walks it from finest to coarsest resolution: for each candidate, builds every study's
    !| neighborhoods and checks the first admissibility gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]);
    !| once every study passes, pools the consensus pmf and checks the second gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]);
    !| once that passes too, seeds a confidence interval with the observed JSD, bootstraps it
    !| ([[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]), and
    !| tests it against the running best candidate for a plateau
    !| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
    !| `plateau_mode` picks which of that CI-overlap criterion and Issue #178's complementary
    !| relative-effect-size one
    !| ([[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]])
    !| governs the stop condition; both are always computed and traced (`trace_*` below) once a
    !| candidate is admissible, regardless of `plateau_mode`, so a caller can compare what either
    !| criterion would have decided. The search stops (`exit`) the moment the SELECTED criterion's
    !| plateau is found -- see `plateau_mode`'s own mode table below for the accepted values.
    !|
    !| Issue #178 also names 2 blocking dependencies for validating the effect-size thresholds
    !| empirically -- the KX_FACTORS default and the Freedman-Diaconis bin-count overestimate --
    !| both deliberately left as-is here; see the project's JSD-Comp-Test follow-up issue. (A third
    !| candidate blocker, a one-sided-vs-symmetric JSD formula question, was raised in the same
    !| follow-up issue but confirmed by the issue's own author to be a mistake in the issue text,
    !| not a real discrepancy -- the code's symmetric formula is correct as written.)
    !| `delta_median_threshold`/`delta_max_threshold` default to the issue's own suggested (not yet
    !| validated) 0.05/0.10.
    !|
    !| When `plateau_mode` selects the effect-size criterion (`MODE_PLATEAU_EFFECT_SIZE` or
    !| `MODE_PLATEAU_BOTH`) and it plateaus independently of the CI-overlap criterion's own running
    !| "best candidate" bookkeeping, `best_candidate_index`/`best_candidate_pair_confidence_interval`
    !| are overridden to the candidate that actually triggered the effect-size plateau, so the
    !| candidate this routine returns is always the one that stopped the search.
    !|
    !| The finally returned candidate is decided by one of three cases:
    !|
    !| 1. A plateau was found, or the grid collapsed to a single candidate (see
    !|    [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]'s
    !|    own small-N collapse note): the best candidate is returned and `plateau_established` is
    !|    `.true.`. A lone candidate is used regardless of whether it plateaued or even passed
    !|    either gate -- the plateau machinery is bypassed entirely, exactly as 125 does. If that
    !|    lone candidate never passed both gates, no real per-point values exist for it, so every
    !|    point's bin count is set to `m_min`, both residual ranges to `0.0`, and
    !|    `best_candidate_pair_confidence_interval` stays `-1.0`.
    !| 2. No plateau, but at least one candidate was admissible: the admissible candidate with the
    !|    smallest bootstrapped uncertainty (median confidence-interval width across studies) is
    !|    returned, with its real confidence interval, and `plateau_established` is `.false.`. This
    !|    applies to every `plateau_mode`.
    !| 3. No plateau, and no candidate was ever admissible: the search falls back to the FIRST
    !|    (finest-resolution) candidate, resets `best_candidate_pair_confidence_interval` to
    !|    `-1.0`, sets every point's bin count to `m_min` and both residual ranges to `0.0`, and
    !|    `plateau_established` is `.false.`.
    !|
    !| Per the plan's work-array translation for this routine specifically: `max_n_bins_all_candidates`
    !| (data-dependent, not cheaply closed-form in 125) is replaced by the fixed
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] ceiling, so every
    !| bin-dimensioned work array below is sized to MAX_N_BINS and sliced `(1:max_n_bins, ...)` per
    !| candidate, rather than carrying a separate recommend-sized dimension argument for it.
    !| `max_n_bins` is the widest per-point bin count Issue #187's occupancy search (Pass B below)
    !| chose for the current candidate, `maxval(tmp_n_bins_per_point(1:n_points))` -- it replaces
    !| the old single scalar `n_bins` that used to come from the global-pool Sturges/FD estimate.
    !| `gene_means` is passed to
    !| [[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]
    !| as its own multi-dimensional self -- that callee declares its matching dummy with an
    !| explicit shape, so standard Fortran sequence association reinterprets the contiguous actual
    !| argument as the flat 1-D array it expects, exactly as 125's own `_alloc` layer did for the
    !| same call.
    !|
    !| Impure: calls the impure
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]. A GSL
    !| failure it reports is folded into `ierr` (first failure only) without aborting the search,
    !| matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_parameter_search_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                       gene_means, residuals, n_bootstraps, &
                                                       join_method, max_n_points_candidate, max_n_neighbors_candidate, &
                                                       n_bootstrapping_top_k_jsds, n_points, n_neighbors, n_bins_per_point, &
                                                       shared_residual_range_low, shared_residual_range_high, &
                                                       best_candidate_pair_confidence_interval, plateau_established, &
                                                       n_admissible_evaluated, &
                                                       trace_n_points, trace_n_neighbors, trace_global_js_divergence, &
                                                       trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                                       trace_ci_width_relative, trace_delta, trace_delta_median, &
                                                       trace_delta_max, trace_selected_n_bins, trace_occupancy_failed, &
                                                       trace_n_pooled_residuals, trace_min_bin_occupancy, &
                                                       trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                                       trace_sturges_bins, trace_fd_bins, &
                                                       trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                                       tmp_gene_means_perms, &
                                                       tmp_gene_means_perm_all, &
                                                       tmp_x_star, tmp_neighborhood_indices_all_studies, &
                                                       tmp_neighborhood_range, tmp_neighborhood_residuals_gathered, &
                                                       tmp_counts_point_major, tmp_pmf_point_major, tmp_n_bins_per_point, &
                                                       tmp_shared_residual_range_low, tmp_shared_residual_range_high, &
                                                       tmp_pmfs, tmp_counts, &
                                                       tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                       tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, &
                                                       tmp_global_js_divergence, tmp_confidence_interval, &
                                                       tmp_bootstrapping_top_k_jsds, tmp_prev_global_js_divergence, &
                                                       tmp_delta_perm, tmp_best_uncertainty_confidence_interval, &
                                                       tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                       tmp_bin_counts_search, tmp_occupancy_failed, &
                                                       tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                       tmp_mean_bin_occupancy, tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                       tmp_fd_bins, tmp_best_n_bins_per_point, &
                                                       tmp_best_uncertainty_n_bins_per_point, &
                                                       tmp_best_shared_residual_range_low, tmp_best_shared_residual_range_high, &
                                                       tmp_best_uncertainty_shared_residual_range_low, &
                                                       tmp_best_uncertainty_shared_residual_range_high, &
                                                       min_residuals_per_bin, min_neighbor_overlap, &
                                                       succeeding_ci_overlap, plateau_mode, delta_median_threshold, &
                                                       delta_max_threshold, delta_epsilon, &
                                                       delta_min_consecutive_transitions, m_min, m_max, gamma_occupancy, &
                                                       lower_residual_range_quantile, upper_residual_range_quantile, &
                                                       two_sided_bootstrapping_significance_level, random_seed, ierr)
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
            !! Per-gene mean expression values for all studies
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_points_candidate
            !! Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
            !! per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
            !! jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
            !! is no longer purely an internal sizing detail the plain wrapper can compute and
            !! hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
            !! rather than AUTO here now
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_impl, JUST_INFO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
            !! because anything returned is sized by it (nothing is), but because it comes from the
            !! same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
            !! now that that call can no longer run automatically inside this wrapper, splitting
            !! this one back into an AUTO call would just be a second, redundant call to the same
            !! routine for no benefit
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_impl, JUST_INFO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_js_comp_test_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate's `n_points`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(int32), dimension(max_n_points_candidate), intent(out) :: n_bins_per_point
            !! The finally chosen candidate's per-point histogram bin count, one per reference
            !! point (Issue #187: every neighborhood may use a different bin count). Only the
            !! leading `n_points` entries are meaningful, mirroring how `n_points`/`n_neighbors`
            !! above are the finally chosen candidate's own values
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The finally chosen candidate's per-point lower residual-range bound (R_low), one per
            !! reference point (Step 3: every neighborhood may use a different, asymmetric range).
            !! Only the leading `n_points` entries are meaningful, mirroring `n_bins_per_point`
            !! above; `0.0_real64` throughout in the two genuinely-degenerate cases where Pass B
            !! never ran for the returned candidate (cases 1 and 3 of the routine description)
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The finally chosen candidate's per-point upper residual-range bound (R_high),
            !! mirroring `shared_residual_range_low` above in every respect
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            !! `-1.0_real64` throughout when no candidate was ever admissible: either the lone
            !! candidate of a collapsed grid never passed both gates (`plateau_established` is
            !! `.true.`), or no plateau was found and no smallest-bootstrap-uncertainty candidate
            !! could be substituted (`plateau_established` is `.false.`)
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a real plateau was found (by whichever criterion
            !! `plateau_mode` selected) or the candidate grid never had more than one candidate to
            !! begin with. `.false.` when the search exhausted every admissible candidate
            !! without ever finding one -- Issue #178's own "report that parameter stability could
            !! not be established". When `.false.` and at least one candidate was admissible, the
            !! routine still returns a real (non-`-1.0`) candidate and confidence interval: the
            !! admissible candidate with the smallest bootstrapped uncertainty, per the issue's own
            !! fallback recommendation -- `plateau_established` is what distinguishes that case from
            !! an actual plateau, not the confidence interval's sentinel value
        integer(int32), intent(out) :: n_admissible_evaluated
            !! Number of candidates that passed both admissibility gates and got a JSD/confidence
            !! interval computed before the search stopped (by plateau or grid exhaustion) -- the
            !! number of leading, valid columns/elements in every `trace_*` array below. This is
            !! Issue #178's own index `t` domain: "position in the ordered sequence of ADMISSIBLE
            !! parameter pairs" -- a candidate that failed either gate has no `trace_*` entry at
            !! all, rather than a zero-filled one
        integer(int32), dimension(16), intent(out) :: trace_n_points
            !! Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
            !! arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
            !! candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(16), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
            !! formula omits this floor, but the same near-zero-JSD instability that motivates
            !! `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
            !! destabilizes any ratio that divides by it, whichever candidate's `J` it is)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point selected histogram bin count
            !! (Issue #187's `M_j`), from determine_bin_count_occupancy_impl. Unlike every OTHER
            !! trace_* array above, whose first extent is a fixed thing like `n_studies`, this
            !! array's first extent is `max_n_points_candidate`, NOT `n_points`, because `n_points`
            !! itself varies per candidate (that is why `trace_n_points(16)` exists as its own
            !! array): this array is genuinely JAGGED per candidate column `t` -- only rows
            !! `1:trace_n_points(t)` are meaningful for that column, rows beyond that are undefined
            !! padding. The result-size directive below only trims the LAST extent (candidates, via
            !! `n_admissible_evaluated`), not this row dimension, so a Python/R caller must
            !! additionally slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
            !! from determine_bin_count_occupancy_impl (never part of the occupancy search's own
            !! decision). Jagged per candidate column exactly as trace_selected_n_bins above --
            !! only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
            !! must slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
            !! diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
            !! search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
            !! above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
            !! caller must slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation for the flattened, all-studies-pooled `gene_means`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array: reference points for the candidate whose `n_points` is current,
            !! recomputed only when `n_points` changes between candidates
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate,
            !! retained simultaneously (unlike the single-study-reused buffer this replaces) so
            !! Pass B below can pool residuals across studies for one reference point at a time,
            !! and Pass C can re-gather each study's own residuals from the already-known indices
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_range
            !! Working array: one study's `[min_idx, max_idx]` neighborhood spans for the current
            !! candidate, reused per study
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values for the current
            !! candidate, reused per study
        integer(int32), dimension(max_n_points_candidate, 256), intent(out) :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts for the current candidate,
            !! reused per study. The `256` is
            !! [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]], written as a
            !! literal because a dimension naming a module parameter the generated wrapper never
            !! `use`s would not compile there -- the same reason MAX_CANDIDATE_PAIRS is written as
            !! `16` on generate_js_comp_test_candidates_impl's own dummy arguments
        real(real64), dimension(max_n_points_candidate, 256), intent(out) :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf for the current candidate, reused per
            !! study. `256` = MAX_N_BINS, see tmp_counts_point_major above
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array: this candidate's per-point histogram bin count, decided by Pass B's
            !! occupancy search
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]])
            !! for each reference point independently
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array: this candidate's per-point lower residual-range bound (R_low), from
            !! Pass B's occupancy search, mirroring tmp_n_bins_per_point above
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array: this candidate's per-point upper residual-range bound (R_high),
            !! mirroring tmp_shared_residual_range_low above
        real(real64), dimension(256, max_n_points_candidate, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's bin-major pmf for the current candidate. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(int32), dimension(256, max_n_points_candidate, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's bin-major histogram counts for the current candidate.
            !! `256` = MAX_N_BINS, see tmp_counts_point_major above
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(256, max_n_points_candidate), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf. `256` = MAX_N_BINS, see
            !! tmp_counts_point_major above
        integer(int32), dimension(256, max_n_points_candidate), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-reference-point JSD values, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-reference-point weights, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD per study, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(real64), dimension(2, n_studies), intent(out) :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, seeded with the observed
            !! global JSD and then bootstrapped in place
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array forwarded to bootstrap_histogram_impl's own top-k/bottom-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study,
            !! forwarded to check_effect_size_plateau_condition_impl
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array forwarded to check_effect_size_plateau_condition_impl's own median
            !! computation
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array: the confidence interval of the admissible candidate with the smallest
            !! bootstrapped uncertainty seen so far, used only internally by the no-plateau fallback
            !! (see plateau_established)
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals across every neighbor and
            !! every study (Pass B), reused per point -- one small buffer, not one per point, since
            !! Pass B is a deliberate sequential loop (see the module-internal doc comment on the
            !! implementation body)
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array forwarded to determine_bin_count_occupancy_impl's own per-bin-count
            !! search scratch, reused per point. `256` = MAX_N_BINS, matching
            !! determine_bin_count_occupancy_impl's own tmp_bin_counts dummy -- written as a
            !! literal because a generated wrapper's dummy dimension cannot reference a module
            !! parameter
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array: this candidate's per-point occupancy_failed flag from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array: this candidate's per-point pooled residual count (N_j) from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array: this candidate's per-point minimum bin occupancy from Pass B
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array: this candidate's per-point mean bin occupancy from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array: this candidate's per-point maximum bin occupancy from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array: this candidate's per-point Sturges' rule diagnostic from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array: this candidate's per-point Freedman-Diaconis rule diagnostic from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array: a snapshot of tmp_n_bins_per_point for whichever candidate is
            !! currently best_candidate_index, updated in lockstep with
            !! check_plateau_condition_impl's own best-candidate bookkeeping (and with the
            !! effect-size override below) -- never read from tmp_n_bins_per_point at loop exit,
            !! which would be stale once the loop has moved on to a later, non-best candidate
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: a snapshot of tmp_n_bins_per_point for whichever admissible candidate
            !! currently has the smallest bootstrapped uncertainty, mirroring how
            !! tmp_best_uncertainty_confidence_interval already works
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array: a snapshot of tmp_shared_residual_range_low for whichever candidate is
            !! currently best_candidate_index, updated at the exact same sites as
            !! tmp_best_n_bins_per_point above
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array: a snapshot of tmp_shared_residual_range_high, mirroring
            !! tmp_best_shared_residual_range_low above
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: a snapshot of tmp_shared_residual_range_low for whichever admissible
            !! candidate currently has the smallest bootstrapped uncertainty, updated at the exact
            !! same site as tmp_best_uncertainty_n_bins_per_point above
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: a snapshot of tmp_shared_residual_range_high, mirroring
            !! tmp_best_uncertainty_shared_residual_range_low above
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate. Reuses Issue #187's occupancy-search default rather than an
            !! independently-tunable threshold of its own: once
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! wires real per-neighborhood bin counts in, a separate laxer threshold here would
            !! silently let a candidate the occupancy search already marked `occupancy_failed`
            !! pass this gate anyway, defeating the FAILURE-detection mechanism
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.1_real64)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.9_real64)
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(CM_MODE_PLATEAU_CI_OVERLAP)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports

        integer(int32) :: candidates_n_points_n_neighbors(2, MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, i_study, i_point, i_neighbor, gene_idx, n_pool, max_n_bins, gather_ierr
        integer(int32) :: prev_n_points, best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: pool_size, bootstrap_ierr, actual_min_residuals_per_bin
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index, actual_m_min
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value, median_ci_width
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found
        logical(c_bool) :: ci_plateau_found, effect_size_plateau_found, has_previous_admissible

        call set_ok(ierr)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, 0.1_real64)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, 0.9_real64)
        M_DEFAULT_VAL(plateau_mode, actual_plateau_mode, CM_MODE_PLATEAU_CI_OVERLAP)
        M_DEFAULT_VAL(delta_median_threshold, actual_delta_median_threshold, CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_max_threshold, actual_delta_max_threshold, CM_DELTA_MAX_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_epsilon, actual_delta_epsilon, CM_DELTA_EPSILON_DEFAULT)
        M_DEFAULT_VAL(delta_min_consecutive_transitions, actual_delta_min_consecutive_transitions, CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        M_DEFAULT_VAL(m_min, actual_m_min, CM_OCCUPANCY_M_MIN_DEFAULT)

        pool_size = max_n_genes_all_studies*n_studies

        call generate_js_comp_test_candidates_impl(max_n_genes_all_studies, candidates_n_points_n_neighbors, n_candidates)

        ! Sort each study's own gene means (for construct_neighborhoods_ranged_impl) and the
        ! flattened, all-studies-pooled gene means (for pool_means_impl's x_star).
        do i_study = 1, n_studies
            call init_perm(tmp_gene_means_perms(:, i_study))
            call sort_array_heapsort(gene_means(:, i_study), tmp_gene_means_perms(:, i_study))
        end do
        call init_perm(tmp_gene_means_perm_all)
        call sort_real_heapsort_expl_size(gene_means, tmp_gene_means_perm_all, pool_size)

        best_candidate_pair_confidence_interval = -1.0_real64
        best_candidate_index = 1_int32
        best_exceeded_ci_overlap_count = 0_int32
        best_uncertainty_value = huge(1.0_real64)
        best_uncertainty_candidate_index = 1_int32
        plateau_found = logical(.false., kind=c_bool)
        prev_n_points = -1_int32
        n_admissible_evaluated = 0_int32
        n_consecutive_effect_size_ok = 0_int32
        has_previous_admissible = logical(.false., kind=c_bool)

        ! Test candidate pairs, finest resolution to coarsest, and stop at the first JSD plateau.
        do i_candidate = 1, n_candidates
            n_points = candidates_n_points_n_neighbors(1, i_candidate)
            n_neighbors = candidates_n_points_n_neighbors(2, i_candidate)

            if (prev_n_points /= n_points) then
                prev_n_points = n_points
                call pool_means_impl(gene_means, tmp_gene_means_perm_all, pool_size, n_points, n_pool, &
                                     tmp_x_star(1:n_points))
            end if

            ! ===== PASS A (per study): neighborhood construction + first admissibility gate.
            ! Short-circuit semantics unchanged from before this restructuring -- a study that
            ! fails the gate stops the remaining studies for this candidate -- but now writes into
            ! tmp_neighborhood_indices_all_studies, which retains EVERY study's indices
            ! simultaneously (rather than reusing one single-study buffer), so Pass B below can pool
            ! residuals across studies for one reference point at a time. Residual gathering and
            ! histogram construction no longer happen here: they move to Pass C below, after Pass B
            ! has decided real per-point bin counts.
            all_have_min_neighbor_overlap = logical(.true., kind=c_bool)
            do i_study = 1, n_studies
                call construct_neighborhoods_ranged_impl(n_points, tmp_x_star(1:n_points), max_n_genes_all_studies, &
                                                         gene_means(:, i_study), tmp_gene_means_perms(:, i_study), &
                                                         n_neighbors, &
                                                         tmp_neighborhood_indices_all_studies(1:n_neighbors, 1:n_points, &
                                                                                               i_study), &
                                                         tmp_neighborhood_range(1:2, 1:n_points))

                call check_neighborhood_overlaps_impl(tmp_neighborhood_range(1:2, 1:n_points), n_points, &
                                                       actual_min_neighbor_overlap, all_have_min_neighbor_overlap)
                if (.not. all_have_min_neighbor_overlap) exit
            end do
            if (.not. all_have_min_neighbor_overlap) cycle

            ! ===== PASS B (per point): Issue #187's occupancy-constrained bin-count search, pooling
            ! every study's residuals for one reference point at a time. Deliberately a plain
            ! sequential `do`, NOT `do concurrent`: this is not the "loops with data-dependent
            ! control flow across iterations" exception from Fortran_Coding_Guides.pdf Sec 10 (each
            ! iteration's work -- gather, heapsort, one occupancy search -- is fully independent of
            ! every other point), but a single small tmp_pooled_residuals buffer is reused across
            ! points here rather than allocating one such buffer per point (up to
            ! max_n_points_candidate of them, which `do concurrent`'s `local()` would require) -- the
            ! memory-conscious choice given each iteration's work is already substantial on its own.
            do i_point = 1, n_points
                call gather_pooled_neighborhood_residuals(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                    n_neighbors, n_studies, tmp_neighborhood_indices_all_studies(1:n_neighbors, i_point, 1:n_studies), &
                    tmp_pooled_residuals(1:max_n_reps_all_studies*n_neighbors*n_studies), gather_ierr)
                if (is_err(gather_ierr)) call set_err_once(ierr, get_err_code(gather_ierr))

                call init_perm(tmp_pooled_residuals_perm(1:max_n_reps_all_studies*n_neighbors*n_studies))
                call sort_array_heapsort(tmp_pooled_residuals(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                                         tmp_pooled_residuals_perm(1:max_n_reps_all_studies*n_neighbors*n_studies))

                call determine_bin_count_occupancy_impl( &
                    tmp_pooled_residuals(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                    tmp_pooled_residuals_perm(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                    max_n_reps_all_studies*n_neighbors*n_studies, max_n_reps_all_studies, n_neighbors, &
                    tmp_n_bins_per_point(i_point), tmp_occupancy_failed(i_point), &
                    tmp_shared_residual_range_low(i_point), tmp_shared_residual_range_high(i_point), &
                    tmp_n_pooled_residuals(i_point), tmp_min_bin_occupancy(i_point), tmp_mean_bin_occupancy(i_point), &
                    tmp_max_bin_occupancy(i_point), tmp_sturges_bins(i_point), tmp_fd_bins(i_point), &
                    tmp_bin_counts_search, m_min, m_max, actual_min_residuals_per_bin, gamma_occupancy, &
                    lower_residual_range_quantile, upper_residual_range_quantile)
            end do

            ! ===== PASS C (per study): re-gather this study's residuals from the indices Pass A
            ! already computed, now that Pass B has decided real per-point bin counts. max_n_bins
            ! (the widest bin count any point in this candidate uses) replaces the old scalar
            ! n_bins for every remaining use in this candidate's processing.
            max_n_bins = maxval(tmp_n_bins_per_point(1:n_points))
            do i_study = 1, n_studies
                do concurrent(i_point=1:n_points, i_neighbor=1:n_neighbors) local(gene_idx) &
                        shared(tmp_neighborhood_residuals_gathered, residuals, tmp_neighborhood_indices_all_studies, i_study)
                    gene_idx = tmp_neighborhood_indices_all_studies(i_neighbor, i_point, i_study)
                    tmp_neighborhood_residuals_gathered(:, i_neighbor, i_point) = residuals(:, gene_idx, i_study)
                end do

                call build_residual_histograms_impl(tmp_neighborhood_residuals_gathered(:, 1:n_neighbors, 1:n_points), &
                                                    max_n_reps_all_studies, n_neighbors, n_points, &
                                                    tmp_shared_residual_range_low(1:n_points), &
                                                    tmp_shared_residual_range_high(1:n_points), &
                                                    max_n_bins, tmp_n_bins_per_point(1:n_points), &
                                                    tmp_counts_point_major(1:n_points, 1:max_n_bins), &
                                                    tmp_pmf_point_major(1:n_points, 1:max_n_bins), &
                                                    tmp_included_n_reps(1:n_points, i_study))
                tmp_counts(1:max_n_bins, 1:n_points, i_study) = transpose(tmp_counts_point_major(1:n_points, 1:max_n_bins))
                tmp_pmfs(1:max_n_bins, 1:n_points, i_study) = transpose(tmp_pmf_point_major(1:n_points, 1:max_n_bins))
            end do

            ! Second admissibility gate: every bin of the consensus pmf must reach the minimum
            ! count, scoped to each point's own tmp_n_bins_per_point.
            call create_mean_pmf_impl(tmp_pmfs(1:max_n_bins, 1:n_points, 1:n_studies), &
                                      tmp_counts(1:max_n_bins, 1:n_points, 1:n_studies), &
                                      max_n_bins, n_points, n_studies, tmp_included_n_reps(1:n_points, 1:n_studies), &
                                      tmp_mean_pmf(1:max_n_bins, 1:n_points), tmp_mean_pmf_included_n_reps(1:n_points), &
                                      tmp_mean_pmf_counts(1:max_n_bins, 1:n_points))
            call check_mean_pmf_min_counts_impl(tmp_mean_pmf_counts(1:max_n_bins, 1:n_points), max_n_bins, &
                                                tmp_n_bins_per_point(1:n_points), n_points, &
                                                actual_min_residuals_per_bin, all_bins_have_min_count)
            if (.not. all_bins_have_min_count) cycle

            ! Observed JSD per study, seeding the confidence interval bootstrap_histogram_impl bootstraps in place.
            do i_study = 1, n_studies
                call compute_divergence_per_reference_point_impl(transpose(tmp_pmfs(1:max_n_bins, 1:n_points, i_study)), &
                                                                 transpose(tmp_mean_pmf(1:max_n_bins, 1:n_points)), n_points, &
                                                                 max_n_bins, tmp_js_divergences(1:n_points, i_study))
                call compute_weighted_global_divergence_impl(tmp_js_divergences(1:n_points, i_study), n_points, &
                                                              tmp_included_n_reps(1:n_points, i_study), &
                                                              tmp_mean_pmf_included_n_reps(1:n_points), &
                                                              tmp_global_js_divergence(i_study), tmp_weights(1:n_points, i_study))
                tmp_confidence_interval(1, i_study) = tmp_global_js_divergence(i_study)
                tmp_confidence_interval(2, i_study) = tmp_global_js_divergence(i_study)
            end do

            ! Record this admissible candidate's diagnostics and test the relative-effect-size
            ! plateau criterion (Issue #178) -- always, regardless of plateau_mode, so a caller can
            ! compare both criteria; only the final select-case below decides which one governs
            ! `exit`. Must happen here, using the true observed tmp_global_js_divergence, BEFORE
            ! bootstrap_histogram_impl reuses that same array as its own scratch below -- reading it
            ! afterward would silently pick up the last bootstrap replicate's value instead of the
            ! observed one.
            n_admissible_evaluated = n_admissible_evaluated + 1_int32
            trace_n_points(n_admissible_evaluated) = n_points
            trace_n_neighbors(n_admissible_evaluated) = n_neighbors
            trace_global_js_divergence(1:n_studies, n_admissible_evaluated) = tmp_global_js_divergence
            trace_selected_n_bins(1:n_points, n_admissible_evaluated) = tmp_n_bins_per_point(1:n_points)
            trace_occupancy_failed(1:n_points, n_admissible_evaluated) = tmp_occupancy_failed(1:n_points)
            trace_n_pooled_residuals(1:n_points, n_admissible_evaluated) = tmp_n_pooled_residuals(1:n_points)
            trace_min_bin_occupancy(1:n_points, n_admissible_evaluated) = tmp_min_bin_occupancy(1:n_points)
            trace_mean_bin_occupancy(1:n_points, n_admissible_evaluated) = tmp_mean_bin_occupancy(1:n_points)
            trace_max_bin_occupancy(1:n_points, n_admissible_evaluated) = tmp_max_bin_occupancy(1:n_points)
            trace_sturges_bins(1:n_points, n_admissible_evaluated) = tmp_sturges_bins(1:n_points)
            trace_fd_bins(1:n_points, n_admissible_evaluated) = tmp_fd_bins(1:n_points)
            trace_shared_residual_range_low(1:n_points, n_admissible_evaluated) = tmp_shared_residual_range_low(1:n_points)
            trace_shared_residual_range_high(1:n_points, n_admissible_evaluated) = tmp_shared_residual_range_high(1:n_points)

            call check_effect_size_plateau_condition_impl(tmp_global_js_divergence, tmp_prev_global_js_divergence, &
                                                          n_studies, has_previous_admissible, actual_delta_median_threshold, &
                                                          actual_delta_max_threshold, actual_delta_epsilon, &
                                                          actual_delta_min_consecutive_transitions, &
                                                          n_consecutive_effect_size_ok, &
                                                          trace_delta(1:n_studies, n_admissible_evaluated), &
                                                          trace_delta_median(n_admissible_evaluated), &
                                                          trace_delta_max(n_admissible_evaluated), &
                                                          effect_size_plateau_found, tmp_delta_perm)
            tmp_prev_global_js_divergence = tmp_global_js_divergence
            has_previous_admissible = logical(.true., kind=c_bool)

            ! Bootstrap the confidence interval -- same n_bootstraps/random_seed for every
            ! candidate, for comparability. The observed-value scratch above is no longer needed,
            ! so it doubles as bootstrap_histogram_impl's own work arrays, exactly as 125 does.
            call bootstrap_histogram_impl(n_bootstraps, max_n_bins, n_points, n_studies, &
                                          tmp_mean_pmf_counts(1:max_n_bins, 1:n_points), &
                                          tmp_mean_pmf_included_n_reps(1:n_points), tmp_included_n_reps(1:n_points, 1:n_studies), &
                                          n_bootstrapping_top_k_jsds, tmp_confidence_interval, &
                                          tmp_bootstrapping_top_k_jsds, tmp_counts(1:max_n_bins, 1:n_points, 1), &
                                          tmp_pmfs(1:max_n_bins, 1:n_points, 1:n_studies), tmp_mean_pmf(1:max_n_bins, 1:n_points), &
                                          tmp_js_divergences(1:n_points, 1:n_studies), tmp_weights(1:n_points, 1:n_studies), &
                                          tmp_global_js_divergence, two_sided_bootstrapping_significance_level, random_seed, &
                                          bootstrap_ierr)
            if (is_err(bootstrap_ierr)) call set_err_once(ierr, get_err_code(bootstrap_ierr))

            trace_ci_lower(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(1, 1:n_studies)
            trace_ci_upper(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(2, 1:n_studies)
            trace_ci_width(1:n_studies, n_admissible_evaluated) = &
                trace_ci_upper(1:n_studies, n_admissible_evaluated) - trace_ci_lower(1:n_studies, n_admissible_evaluated)
            trace_ci_width_relative(1:n_studies, n_admissible_evaluated) = &
                trace_ci_width(1:n_studies, n_admissible_evaluated) &
                / max(trace_global_js_divergence(1:n_studies, n_admissible_evaluated), actual_delta_epsilon)

            ! Track the admissible candidate with the smallest bootstrapped uncertainty
            ! (median confidence-interval width across studies), regardless of plateau_mode -- the
            ! ranking itself doesn't depend on which plateau criterion is selected, only whether the
            ! fallback branch below actually uses it. tmp_delta_perm was just used above to sort
            ! `delta`; it must be re-seeded and re-sorted here for trace_ci_width's own order before
            ! calc_percentile_impl can use it, since calc_percentile_impl trusts whatever order it's
            ! handed rather than sorting internally.
            call init_perm(tmp_delta_perm)
            call sort_array_heapsort(trace_ci_width(1:n_studies, n_admissible_evaluated), tmp_delta_perm)
            call calc_percentile_impl(trace_ci_width(1:n_studies, n_admissible_evaluated), n_studies, tmp_delta_perm, &
                                      0.5_real64, median_ci_width)
            if (median_ci_width < best_uncertainty_value) then
                best_uncertainty_value = median_ci_width
                best_uncertainty_candidate_index = i_candidate
                tmp_best_uncertainty_confidence_interval = tmp_confidence_interval
                tmp_best_uncertainty_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_uncertainty_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_uncertainty_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            call check_plateau_condition_impl(tmp_confidence_interval, best_candidate_pair_confidence_interval, n_studies, &
                                              best_candidate_index, best_exceeded_ci_overlap_count, i_candidate, join_method, &
                                              actual_succeeding_ci_overlap, ci_plateau_found)

            ! check_plateau_condition_impl mutates best_candidate_index internally whenever THIS
            ! candidate becomes the new best (its own doc comment: "overwritten ... unless the new
            ! candidate is worse"); it is out of scope to modify, so detect that here and snapshot
            ! this candidate's real per-point bin counts before the loop moves on to another one.
            if (best_candidate_index == i_candidate) then
                tmp_best_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            select case (actual_plateau_mode)
            case (MODE_PLATEAU_CI_OVERLAP)
                plateau_found = ci_plateau_found
            case (MODE_PLATEAU_EFFECT_SIZE)
                plateau_found = effect_size_plateau_found
            case default ! MODE_PLATEAU_BOTH
                plateau_found = logical(ci_plateau_found .or. effect_size_plateau_found, kind=c_bool)
            end select

            ! check_plateau_condition_impl's own best-candidate bookkeeping tracks its OWN
            ! (CI-overlap) criterion only. When the effect-size criterion is the one that actually
            ! plateaued (and plateau_mode gives it a say), override to the candidate that triggered
            ! it, so the candidate this routine returns is always the one that stopped the search.
            if (effect_size_plateau_found .and. actual_plateau_mode /= MODE_PLATEAU_CI_OVERLAP) then
                best_candidate_index = i_candidate
                best_candidate_pair_confidence_interval = tmp_confidence_interval
                tmp_best_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            if (plateau_found) exit
        end do

        ! Final candidate pair, one of three cases:
        ! 1. A plateau was reached, or the grid never had more than one candidate to begin with (the
        !    single-candidate/collapsed-grid path bypasses the plateau machinery entirely, exactly
        !    as 125 does): use the best candidate found, plateau_established = .true.
        ! 2. No plateau, and at least one candidate was ever admissible: Issue #178's own fallback --
        !    the admissible candidate with the smallest bootstrapped uncertainty, a real (non--1.0)
        !    confidence interval, plateau_established = .false. Applies to every plateau_mode -- the
        !    ranking computation itself doesn't depend on which plateau criterion was selected (see
        !    the comment where it's computed, a few dozen lines above the loop's own exit), only
        !    which mode's own plateau-detection condition it was gated behind, and that gate is gone
        !    now that every mode has been validated end-to-end against real Kidney data (Step 3's
        !    own m_max sweep, all 3 modes).
        ! 3. No plateau, and zero candidates were ever admissible (n_candidates counts the raw grid
        !    including gate-failed candidates, so this is distinct from case 2's "at least one
        !    admissible" -- nothing exists for case 2 to select from here): fall back to the
        !    finest-resolution (first) candidate and reset the confidence interval to -1.0, exactly
        !    as before this change, plateau_established = .false.
        if (plateau_found .or. n_candidates < 2_int32) then
            n_points = candidates_n_points_n_neighbors(1, best_candidate_index)
            n_neighbors = candidates_n_points_n_neighbors(2, best_candidate_index)
            if (n_admissible_evaluated >= 1_int32) then
                ! The normal case: best_candidate_index became admissible at some point during the
                ! loop, so tmp_best_n_bins_per_point was genuinely snapshotted for it.
                n_bins_per_point(1:n_points) = tmp_best_n_bins_per_point(1:n_points)
                shared_residual_range_low(1:n_points) = tmp_best_shared_residual_range_low(1:n_points)
                shared_residual_range_high(1:n_points) = tmp_best_shared_residual_range_high(1:n_points)
            else
                ! Degenerate: n_candidates < 2 and that lone candidate never passed gate 1, so
                ! Pass B never ran for it and tmp_best_n_bins_per_point was never snapshotted.
                ! generate_js_comp_test_candidates_impl no longer computes any bin estimate at all
                ! (Issue #187, Step 2.7) -- there is no real per-neighborhood value available for
                ! this candidate, so fall back to actual_m_min, the smallest bin count the
                ! occupancy search would ever try, as the most conservative available choice. This
                ! mirrors determine_bin_count_occupancy_impl's own FAILURE case, which uses m_min
                ! as its fallback for exactly the same "we don't have enough information" reason.
                ! The range has no equivalent fallback value -- Pass B never ran, so there is no
                ! data to derive a percentile from -- so it is zeroed, matching
                ! determine_bin_count_occupancy_impl's own all-NaN degenerate convention.
                n_bins_per_point(1:n_points) = actual_m_min
                shared_residual_range_low(1:n_points) = 0.0_real64
                shared_residual_range_high(1:n_points) = 0.0_real64
            end if
            plateau_established = logical(.true., kind=c_bool)
        else if (n_admissible_evaluated >= 1_int32) then
            n_points = candidates_n_points_n_neighbors(1, best_uncertainty_candidate_index)
            n_neighbors = candidates_n_points_n_neighbors(2, best_uncertainty_candidate_index)
            ! This branch's own condition already guarantees n_admissible_evaluated >= 1, so
            ! tmp_best_uncertainty_n_bins_per_point was genuinely snapshotted -- no fallback needed.
            n_bins_per_point(1:n_points) = tmp_best_uncertainty_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = tmp_best_uncertainty_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = tmp_best_uncertainty_shared_residual_range_high(1:n_points)
            best_candidate_pair_confidence_interval = tmp_best_uncertainty_confidence_interval
            plateau_established = logical(.false., kind=c_bool)
        else
            n_points = candidates_n_points_n_neighbors(1, 1)
            n_neighbors = candidates_n_points_n_neighbors(2, 1)
            ! Case 3, Issue #178's own "fall back to the finest-resolution (first) candidate":
            ! candidate 1 may never have passed gate 1 (this branch has no cheap way to know
            ! without re-running Pass A/B, and it is not worth it for an already-degenerate,
            ! CI-reset-to--1.0 fallback path), so there is no real per-neighborhood bin count
            ! available for it either. generate_js_comp_test_candidates_impl no longer computes any
            ! bin estimate at all (Issue #187, Step 2.7), so fall back to actual_m_min -- the
            ! smallest bin count the occupancy search would ever try -- exactly as the other
            ! degenerate branch above does, and consistent with
            ! determine_bin_count_occupancy_impl's own FAILURE case using m_min for the same reason.
            ! The range is zeroed for the same reason as the other genuinely-degenerate branch above.
            n_bins_per_point(1:n_points) = actual_m_min
            shared_residual_range_low(1:n_points) = 0.0_real64
            shared_residual_range_high(1:n_points) = 0.0_real64
            best_candidate_pair_confidence_interval = -1.0_real64
            plateau_established = logical(.false., kind=c_bool)
        end if
    end subroutine run_js_comp_test_parameter_search_impl

    !> summary: Search the ascending adaptive (k_start, k_step, k_max) candidate sequence for a stable JSD parameter setting (Issue #217)
    !| AUTHOR_LASZLO_LANG
    !| The adaptive-neighborhood counterpart of
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !| Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
    !| in order. Their `k_start` shrinks, so the number of reference points that emerges ascends.
    !| The plateau criteria, `plateau_mode` and the fallbacks are exactly the fixed-k search's.
    !|
    !| **Per candidate**, in order:
    !|
    !| 1. Grow its neighborhoods over the pooled, sorted gene means. The pooled means are sorted
    !|    once for the whole search, with the same routine on the same input as
    !|    [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]],
    !|    and the growth is that routine's own, with the same `tau`, `mad_distance_factor`,
    !|    `max_pooled_residuals` and `min_study_neighbors`. So a candidate's neighborhoods are
    !|    exactly what `construct_adaptive_neighborhoods` returns for its knobs.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !|    reference points emerged.
    !| 3. First admissibility gate:
    !|    [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !|    on the pooled ranges `[a_i, b_i]`.
    !| 4. Per point, decode the range into per-study gene lists and run Issue #187's occupancy
    !|    search on the pooled residuals (Pass B). Then build each study's histogram per point
    !|    (Pass C), the consensus pmf, and apply the second gate,
    !|    [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]].
    !|    This is the pipeline of
    !|    [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]].
    !| 5. From here on everything is the fixed-k search's post-gate bookkeeping, keyed on the
    !|    candidate index: the observed JSD and the traces, the effect-size check (before the
    !|    bootstrap, which reuses the observed-JSD buffer as scratch), the bootstrap, the
    !|    smallest-uncertainty tracking, the CI-overlap check, the `plateau_mode` selection, the
    !|    effect-size override, and the stop at the first plateau.
    !|
    !| **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
    !| `candidate_n_points` and `candidate_status`, in order, until the search stops. The status is
    !|
    !| - [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_EVALUATED(variable)]]
    !|   (`CM_ADAPTIVE_CANDIDATE_EVALUATED`): passed both gates and has a trace column;
    !| - [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED(variable)]]
    !|   (`CM_ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED`): construction status too few means or too
    !|   few residuals, so the growth could not run as specified;
    !| - [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_EMPTY_STUDY(variable)]]
    !|   (`CM_ADAPTIVE_CANDIDATE_EMPTY_STUDY`): construction status empty study neighborhood;
    !| - [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED(variable)]]
    !|   (`CM_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED`): more than `max_n_points_candidate` points;
    !| - [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_OVERLAP_FAILED(variable)]]
    !|   (`CM_ADAPTIVE_CANDIDATE_OVERLAP_FAILED`): the first gate failed;
    !| - [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED(variable)]]
    !|   (`CM_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED`): the second gate failed.
    !|
    !| The checks run in that order and the first failing one is logged. `candidate_n_points` is the
    !| number of reference points the construction produced (0 for too few means), whatever the
    !| status.
    !|
    !| **Result**, one of three cases, as in the fixed-k search:
    !|
    !| 1. A plateau was found, or the sequence had a single candidate: that candidate's knobs,
    !|    `plateau_established = .true.`.
    !| 2. No plateau but at least one admissible candidate: the admissible candidate with the
    !|    smallest median confidence-interval width, `plateau_established = .false.`.
    !| 3. Nothing admissible: candidate 1's knobs, confidence interval `-1.0`,
    !|    `plateau_established = .false.`.
    !|
    !| `n_points` is the selected candidate's own logged point count. Where no admissible candidate
    !| backs the result (case 3, or case 1 with a single, rejected candidate), every
    !| `n_bins_per_point` entry up to `min(n_points, max_n_points_candidate)` is `m_min` and both
    !| ranges are `0.0`, as in the fixed-k search. To rebuild the selected neighborhoods, call
    !| `construct_adaptive_neighborhoods` with the returned knobs and the same construction
    !| options, then `run_js_comp_test_adaptive` with the same occupancy options: its
    !| `global_js_divergence` equals the selected candidate's `trace_global_js_divergence` column.
    !|
    !| **Traces** are the fixed-k search's, with `trace_k_start`/`trace_k_step`/`trace_k_max` in
    !| place of `trace_n_neighbors`. The per-point traces are jagged per column exactly as there.
    !|
    !| **Memory.** Beyond the fixed-k search's per-point arrays, which are sized by
    !| `max_n_points_candidate` (practical capacity from
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]),
    !| the construction needs buffers over all `N = max_n_genes_all_studies * n_studies` pooled
    !| entries, the largest being `tmp_n_neighbors_per_point` with `n_studies * N` integers. The
    !| pooling buffers `tmp_pooled_residuals`/`tmp_pooled_residuals_perm` hold
    !| `max_n_reps_all_studies * N` values each, the size of `residuals`, because one study may
    !| have all its genes in one neighborhood. The histogram buffers hold `256 * max_n_points_candidate
    !| * n_studies` values each.
    !|
    !| Impure: calls the impure
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]. A GSL
    !| failure it reports is folded into `ierr` (first failure only) without aborting the search,
    !| as in the fixed-k search.
    subroutine run_js_comp_test_adaptive_parameter_search_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                                gene_means, residuals, n_bootstraps, join_method, &
                                                                max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                                                k_start, k_step, k_max, n_points, n_bins_per_point, &
                                                                shared_residual_range_low, shared_residual_range_high, &
                                                                best_candidate_confidence_interval, plateau_established, &
                                                                n_admissible_evaluated, trace_k_start, trace_k_step, &
                                                                trace_k_max, trace_n_points, trace_global_js_divergence, &
                                                                trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                                                trace_ci_width_relative, trace_delta, trace_delta_median, &
                                                                trace_delta_max, trace_selected_n_bins, trace_occupancy_failed, &
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
                                                                tmp_confidence_interval, tmp_bootstrapping_top_k_jsds, &
                                                                tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                                tmp_best_uncertainty_confidence_interval, &
                                                                tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                                tmp_bin_counts_search, tmp_occupancy_failed, &
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
                                                                min_study_neighbors, ierr)
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
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_points_candidate
            !! Most reference points a candidate may emerge with; a candidate with more is rejected
            !! with the capacity status. The per-point outputs and traces below are sized by it, so
            !! the caller has to know it up front
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_bounds, tox_data_integration_js_comp_test_impl, JUST_INFO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_js_comp_test_impl, AUTO)
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
            !! The finally chosen candidate's per-point histogram bin count; only the leading
            !! `n_points` entries are meaningful. Where no admissible candidate backs the result,
            !! the leading `min(n_points, max_n_points_candidate)` entries are `m_min`, since that
            !! candidate's point count may exceed the capacity
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The finally chosen candidate's per-point lower residual-range bound (R_low); only the
            !! leading `n_points` entries are meaningful, `0.0_real64` where no admissible candidate
            !! backs the result
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The finally chosen candidate's per-point upper residual-range bound (R_high),
            !! mirroring `shared_residual_range_low`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_confidence_interval
            !! The bootstrapped JSD confidence interval of the finally chosen candidate;
            !! `-1.0_real64` throughout where no admissible candidate backs the result
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! sequence had a single candidate; `.false.` when the search exhausted the sequence
            !! without one, in which case the smallest-uncertainty admissible candidate is
            !! returned if there is one
        integer(int32), intent(out) :: n_admissible_evaluated
            !! Number of candidates that passed both admissibility gates and got a JSD and
            !! confidence interval before the search stopped: the number of leading, valid
            !! columns/elements of every `trace_*` array
        integer(int32), dimension(16), intent(out) :: trace_k_start
            !! Per-admissible-candidate `k_start`. `16` = MAX_CANDIDATE_PAIRS, the adaptive
            !! sequence's cap
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(16), intent(out) :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(16), intent(out) :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(16), intent(out) :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study confidence-interval width divided by the observed
            !! global JSD, the denominator floored at `delta_epsilon`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0_real64` throughout at the first one
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0_real64` at the
            !! first one
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0_real64` at the
            !! first one
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point selected histogram bin count. Jagged:
            !! only rows `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller
            !! must additionally slice `[:trace_n_points[t], t]` themselves
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point `occupancy_failed` flag; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic, with the
            !! rounded mean per-study neighbor count; jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic, with
            !! the rounded mean per-study neighbor count; jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), intent(out) :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(int32), dimension(16), intent(out) :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(16), intent(out) :: candidate_n_points
            !! Number of reference points every tried candidate's construction produced (0 for too
            !! few means), whether or not the candidate was admissible
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(16), intent(out) :: candidate_status
            !! What happened to every tried candidate: `CM_ADAPTIVE_CANDIDATE_EVALUATED` evaluated,
            !! `CM_ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED` construction failed,
            !! `CM_ADAPTIVE_CANDIDATE_EMPTY_STUDY` a study neighborhood below `min_study_neighbors`,
            !! `CM_ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED` capacity exceeded,
            !! `CM_ADAPTIVE_CANDIDATE_OVERLAP_FAILED` overlap gate failed,
            !! `CM_ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED` consensus-pmf count gate failed
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, sorted once, exactly as
            !! construct_adaptive_neighborhoods sorts it
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
            !! Working array: one (reference point, study) pair's gathered residuals, reused (Pass C)
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array: the current candidate's per-point bin count from Pass B
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array: the current candidate's per-point R_low from Pass B
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array: the current candidate's per-point R_high from Pass B
        real(real64), dimension(256, max_n_points_candidate, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's bin-major pmf for the current candidate. `256` = MAX_N_BINS
        integer(int32), dimension(256, max_n_points_candidate, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's bin-major histogram counts for the current candidate.
            !! `256` = MAX_N_BINS
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(256, max_n_points_candidate), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf. `256` = MAX_N_BINS
        integer(int32), dimension(256, max_n_points_candidate), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts. `256` = MAX_N_BINS
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, reused as
            !! bootstrap_histogram_impl's scratch once consumed
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-point weights, reused as
            !! bootstrap_histogram_impl's scratch once consumed
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD per study, reused as
            !! bootstrap_histogram_impl's scratch once consumed
        real(real64), dimension(2, n_studies), intent(out) :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, seeded with the observed
            !! global JSD and then bootstrapped in place
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array forwarded to bootstrap_histogram_impl's top-k/bottom-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array: sorting permutation for the per-study medians
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array: the confidence interval of the admissible candidate with the smallest
            !! bootstrapped uncertainty seen so far
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals (Pass B), reused per point; as
            !! large as `residuals`
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array forwarded to the occupancy search, reused per point. `256` = MAX_N_BINS
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array: the current candidate's per-point occupancy_failed flag from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array: the current candidate's per-point pooled residual count from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array: the current candidate's per-point minimum bin occupancy from Pass B
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array: the current candidate's per-point mean bin occupancy from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array: the current candidate's per-point maximum bin occupancy from Pass B
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array: the current candidate's per-point Sturges' rule diagnostic from Pass B
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
            !! Working array: snapshot of tmp_shared_residual_range_low for the smallest-uncertainty
            !! candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: snapshot of tmp_shared_residual_range_high for the smallest-uncertainty
            !! candidate
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count every bin must reach, both in the occupancy search and in the second
            !! admissibility gate (the consensus pmf)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods' pooled ranges must have to
            !! pass the first admissibility gate
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.1_real64)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.9_real64)
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(CM_MODE_PLATEAU_CI_OVERLAP)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible JSD was
            !! zero
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the occupancy search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the occupancy search tests (M_max); raised to `m_min`
            !! when smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search's coarse stage
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's lower residual-range bound
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's upper residual-range bound
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        real(real64), intent(in), optional :: tau
            !! Largest relative increase of the dispersion an adaptive growth round may cause and
            !! still be committed, as in construct_adaptive_neighborhoods
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's median absolute deviation that the next reference
            !! point's target lies beyond it, as in construct_adaptive_neighborhoods
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on the non-NaN residuals adaptive rounds may grow a neighborhood's pool to; 0
            !! for no cap, as in construct_adaptive_neighborhoods
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible, as in construct_adaptive_neighborhoods
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(out) :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports

        integer(int32) :: candidates_k_start_k_step_k_max(3, MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, i_study, i_point, max_n_bins, n_study_neighbors, point_ierr
        integer(int32) :: best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: pool_size, bootstrap_ierr, actual_min_residuals_per_bin
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index, actual_m_min
        integer(int32) :: actual_max_pooled_residuals, actual_min_study_neighbors, construction_status, max_n_neighbors
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value, median_ci_width, actual_tau, actual_mad_distance_factor
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found
        logical(c_bool) :: ci_plateau_found, effect_size_plateau_found, has_previous_admissible

        call set_ok(ierr)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, 0.1_real64)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, 0.9_real64)
        M_DEFAULT_VAL(plateau_mode, actual_plateau_mode, CM_MODE_PLATEAU_CI_OVERLAP)
        M_DEFAULT_VAL(delta_median_threshold, actual_delta_median_threshold, CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_max_threshold, actual_delta_max_threshold, CM_DELTA_MAX_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_epsilon, actual_delta_epsilon, CM_DELTA_EPSILON_DEFAULT)
        M_DEFAULT_VAL(delta_min_consecutive_transitions, actual_delta_min_consecutive_transitions, CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        M_DEFAULT_VAL(m_min, actual_m_min, CM_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(tau, actual_tau, CM_ADAPTIVE_TAU_DEFAULT)
        M_DEFAULT_VAL(mad_distance_factor, actual_mad_distance_factor, CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        M_DEFAULT_VAL(max_pooled_residuals, actual_max_pooled_residuals, CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        M_DEFAULT_VAL(min_study_neighbors, actual_min_study_neighbors, CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)

        pool_size = max_n_genes_all_studies*n_studies

        call generate_adaptive_js_comp_test_candidates_impl(max_n_genes_all_studies, n_studies, &
                                                            candidates_k_start_k_step_k_max, n_candidates)

        ! The construction's own sort (construct_adaptive_neighborhoods_impl), done once: every
        ! candidate grows on it, and a pooled position names the same entry as there.
        call init_perm(tmp_gene_means_perm_all)
        call sort_real_heapsort_expl_size(gene_means, tmp_gene_means_perm_all, pool_size)

        best_candidate_confidence_interval = -1.0_real64
        best_candidate_index = 1_int32
        best_exceeded_ci_overlap_count = 0_int32
        best_uncertainty_value = huge(1.0_real64)
        best_uncertainty_candidate_index = 1_int32
        plateau_found = logical(.false., kind=c_bool)
        n_admissible_evaluated = 0_int32
        n_candidates_tried = 0_int32
        n_consecutive_effect_size_ok = 0_int32
        has_previous_admissible = logical(.false., kind=c_bool)

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

            ! ===== PASS B (per point, sequential: one pooling buffer serves every point), exactly
            ! as run_js_comp_test_adaptive_impl.
            do i_point = 1, n_points
                call materialize_pooled_neighborhood(tmp_gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                     tmp_pooled_neighborhood_range(1, i_point), &
                                                     tmp_pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                     tmp_point_neighborhood_indices, tmp_n_neighbors_per_point(:, i_point))
                call determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                               max_n_genes_all_studies, tmp_point_neighborhood_indices, &
                                               tmp_n_neighbors_per_point(:, i_point), tmp_n_bins_per_point(i_point), &
                                               tmp_occupancy_failed(i_point), tmp_shared_residual_range_low(i_point), &
                                               tmp_shared_residual_range_high(i_point), tmp_n_pooled_residuals(i_point), &
                                               tmp_min_bin_occupancy(i_point), tmp_mean_bin_occupancy(i_point), &
                                               tmp_max_bin_occupancy(i_point), tmp_sturges_bins(i_point), &
                                               tmp_fd_bins(i_point), tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                               tmp_bin_counts_search, point_ierr, m_min=m_min, m_max=m_max, &
                                               min_residuals_per_bin=actual_min_residuals_per_bin, &
                                               gamma_occupancy=gamma_occupancy, &
                                               lower_residual_range_quantile=lower_residual_range_quantile, &
                                               upper_residual_range_quantile=upper_residual_range_quantile)
                ! Defensive: determine_point_bin_count only fails on an out-of-range gene index, and
                ! every decoded index lies in [1, max_n_genes_all_studies], so this branch cannot be
                ! reached and has no test.
                if (is_err(point_ierr)) then
                    call set_err_once(ierr, get_err_code(point_ierr))
                    return
                end if
            end do

            ! ===== PASS C (per point, per study, sequential: one gather buffer serves every pair).
            max_n_bins = maxval(tmp_n_bins_per_point(1:n_points))
            do i_point = 1, n_points
                call materialize_pooled_neighborhood(tmp_gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                     tmp_pooled_neighborhood_range(1, i_point), &
                                                     tmp_pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                     tmp_point_neighborhood_indices, tmp_n_neighbors_per_point(:, i_point))
                do i_study = 1, n_studies
                    n_study_neighbors = tmp_n_neighbors_per_point(i_study, i_point)
                    call build_point_study_histogram(residuals(:, :, i_study), max_n_reps_all_studies, &
                                                     max_n_genes_all_studies, n_study_neighbors, &
                                                     tmp_point_neighborhood_indices(1:n_study_neighbors, i_study), &
                                                     tmp_shared_residual_range_low(i_point), &
                                                     tmp_shared_residual_range_high(i_point), &
                                                     tmp_n_bins_per_point(i_point), max_n_bins, &
                                                     tmp_counts(1:max_n_bins, i_point, i_study), &
                                                     tmp_pmfs(1:max_n_bins, i_point, i_study), &
                                                     tmp_included_n_reps(i_point, i_study), &
                                                     tmp_neighbor_residuals(:, 1:n_study_neighbors))
                end do
            end do

            ! Second admissibility gate: every bin of the consensus pmf must reach the minimum
            ! count, scoped to each point's own tmp_n_bins_per_point.
            call create_mean_pmf_impl(tmp_pmfs(1:max_n_bins, 1:n_points, 1:n_studies), &
                                      tmp_counts(1:max_n_bins, 1:n_points, 1:n_studies), &
                                      max_n_bins, n_points, n_studies, tmp_included_n_reps(1:n_points, 1:n_studies), &
                                      tmp_mean_pmf(1:max_n_bins, 1:n_points), tmp_mean_pmf_included_n_reps(1:n_points), &
                                      tmp_mean_pmf_counts(1:max_n_bins, 1:n_points))
            call check_mean_pmf_min_counts_impl(tmp_mean_pmf_counts(1:max_n_bins, 1:n_points), max_n_bins, &
                                                tmp_n_bins_per_point(1:n_points), n_points, &
                                                actual_min_residuals_per_bin, all_bins_have_min_count)
            if (.not. all_bins_have_min_count) then
                candidate_status(i_candidate) = ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED
                cycle
            end if

            ! From here on: the fixed-k search's post-gate bookkeeping, unchanged in logic, keyed on
            ! i_candidate, with the adaptive knobs traced in place of n_neighbors.

            ! Observed JSD per study, seeding the confidence interval bootstrap_histogram_impl bootstraps in place.
            do i_study = 1, n_studies
                call compute_divergence_per_reference_point_impl(transpose(tmp_pmfs(1:max_n_bins, 1:n_points, i_study)), &
                                                                 transpose(tmp_mean_pmf(1:max_n_bins, 1:n_points)), n_points, &
                                                                 max_n_bins, tmp_js_divergences(1:n_points, i_study))
                call compute_weighted_global_divergence_impl(tmp_js_divergences(1:n_points, i_study), n_points, &
                                                              tmp_included_n_reps(1:n_points, i_study), &
                                                              tmp_mean_pmf_included_n_reps(1:n_points), &
                                                              tmp_global_js_divergence(i_study), tmp_weights(1:n_points, i_study))
                tmp_confidence_interval(1, i_study) = tmp_global_js_divergence(i_study)
                tmp_confidence_interval(2, i_study) = tmp_global_js_divergence(i_study)
            end do

            ! Record the diagnostics and test the effect-size criterion BEFORE the bootstrap, which
            ! reuses tmp_global_js_divergence as its own scratch.
            n_admissible_evaluated = n_admissible_evaluated + 1_int32
            trace_k_start(n_admissible_evaluated) = k_start
            trace_k_step(n_admissible_evaluated) = k_step
            trace_k_max(n_admissible_evaluated) = k_max
            trace_n_points(n_admissible_evaluated) = n_points
            trace_global_js_divergence(1:n_studies, n_admissible_evaluated) = tmp_global_js_divergence
            trace_selected_n_bins(1:n_points, n_admissible_evaluated) = tmp_n_bins_per_point(1:n_points)
            trace_occupancy_failed(1:n_points, n_admissible_evaluated) = tmp_occupancy_failed(1:n_points)
            trace_n_pooled_residuals(1:n_points, n_admissible_evaluated) = tmp_n_pooled_residuals(1:n_points)
            trace_min_bin_occupancy(1:n_points, n_admissible_evaluated) = tmp_min_bin_occupancy(1:n_points)
            trace_mean_bin_occupancy(1:n_points, n_admissible_evaluated) = tmp_mean_bin_occupancy(1:n_points)
            trace_max_bin_occupancy(1:n_points, n_admissible_evaluated) = tmp_max_bin_occupancy(1:n_points)
            trace_sturges_bins(1:n_points, n_admissible_evaluated) = tmp_sturges_bins(1:n_points)
            trace_fd_bins(1:n_points, n_admissible_evaluated) = tmp_fd_bins(1:n_points)
            trace_shared_residual_range_low(1:n_points, n_admissible_evaluated) = tmp_shared_residual_range_low(1:n_points)
            trace_shared_residual_range_high(1:n_points, n_admissible_evaluated) = tmp_shared_residual_range_high(1:n_points)

            call check_effect_size_plateau_condition_impl(tmp_global_js_divergence, tmp_prev_global_js_divergence, &
                                                          n_studies, has_previous_admissible, actual_delta_median_threshold, &
                                                          actual_delta_max_threshold, actual_delta_epsilon, &
                                                          actual_delta_min_consecutive_transitions, &
                                                          n_consecutive_effect_size_ok, &
                                                          trace_delta(1:n_studies, n_admissible_evaluated), &
                                                          trace_delta_median(n_admissible_evaluated), &
                                                          trace_delta_max(n_admissible_evaluated), &
                                                          effect_size_plateau_found, tmp_delta_perm)
            tmp_prev_global_js_divergence = tmp_global_js_divergence
            has_previous_admissible = logical(.true., kind=c_bool)

            ! Bootstrap the confidence interval -- same n_bootstraps/random_seed for every candidate.
            call bootstrap_histogram_impl(n_bootstraps, max_n_bins, n_points, n_studies, &
                                          tmp_mean_pmf_counts(1:max_n_bins, 1:n_points), &
                                          tmp_mean_pmf_included_n_reps(1:n_points), tmp_included_n_reps(1:n_points, 1:n_studies), &
                                          n_bootstrapping_top_k_jsds, tmp_confidence_interval, &
                                          tmp_bootstrapping_top_k_jsds, tmp_counts(1:max_n_bins, 1:n_points, 1), &
                                          tmp_pmfs(1:max_n_bins, 1:n_points, 1:n_studies), tmp_mean_pmf(1:max_n_bins, 1:n_points), &
                                          tmp_js_divergences(1:n_points, 1:n_studies), tmp_weights(1:n_points, 1:n_studies), &
                                          tmp_global_js_divergence, two_sided_bootstrapping_significance_level, random_seed, &
                                          bootstrap_ierr)
            if (is_err(bootstrap_ierr)) call set_err_once(ierr, get_err_code(bootstrap_ierr))

            trace_ci_lower(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(1, 1:n_studies)
            trace_ci_upper(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(2, 1:n_studies)
            trace_ci_width(1:n_studies, n_admissible_evaluated) = &
                trace_ci_upper(1:n_studies, n_admissible_evaluated) - trace_ci_lower(1:n_studies, n_admissible_evaluated)
            trace_ci_width_relative(1:n_studies, n_admissible_evaluated) = &
                trace_ci_width(1:n_studies, n_admissible_evaluated) &
                / max(trace_global_js_divergence(1:n_studies, n_admissible_evaluated), actual_delta_epsilon)

            ! Smallest bootstrapped uncertainty (median CI width across studies), regardless of
            ! plateau_mode. tmp_delta_perm is re-seeded for trace_ci_width's own order.
            call init_perm(tmp_delta_perm)
            call sort_array_heapsort(trace_ci_width(1:n_studies, n_admissible_evaluated), tmp_delta_perm)
            call calc_percentile_impl(trace_ci_width(1:n_studies, n_admissible_evaluated), n_studies, tmp_delta_perm, &
                                      0.5_real64, median_ci_width)
            if (median_ci_width < best_uncertainty_value) then
                best_uncertainty_value = median_ci_width
                best_uncertainty_candidate_index = i_candidate
                tmp_best_uncertainty_confidence_interval = tmp_confidence_interval
                tmp_best_uncertainty_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_uncertainty_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_uncertainty_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            call check_plateau_condition_impl(tmp_confidence_interval, best_candidate_confidence_interval, n_studies, &
                                              best_candidate_index, best_exceeded_ci_overlap_count, i_candidate, join_method, &
                                              actual_succeeding_ci_overlap, ci_plateau_found)

            ! Snapshot this candidate's per-point values whenever it became the CI-overlap best.
            if (best_candidate_index == i_candidate) then
                tmp_best_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            select case (actual_plateau_mode)
            case (MODE_PLATEAU_CI_OVERLAP)
                plateau_found = ci_plateau_found
            case (MODE_PLATEAU_EFFECT_SIZE)
                plateau_found = effect_size_plateau_found
            case default ! MODE_PLATEAU_BOTH
                plateau_found = logical(ci_plateau_found .or. effect_size_plateau_found, kind=c_bool)
            end select

            ! When the effect-size criterion plateaued and plateau_mode gives it a say, the
            ! candidate that triggered it is the one returned.
            if (effect_size_plateau_found .and. actual_plateau_mode /= MODE_PLATEAU_CI_OVERLAP) then
                best_candidate_index = i_candidate
                best_candidate_confidence_interval = tmp_confidence_interval
                tmp_best_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            if (plateau_found) exit
        end do

        ! Final candidate, the fixed-k search's three cases (see the doc comment). The point count
        ! is the chosen candidate's logged one; the knobs come from the sequence.
        if (plateau_found .or. n_candidates < 2_int32) then
            k_start = candidates_k_start_k_step_k_max(1, best_candidate_index)
            k_step = candidates_k_start_k_step_k_max(2, best_candidate_index)
            k_max = candidates_k_start_k_step_k_max(3, best_candidate_index)
            n_points = candidate_n_points(best_candidate_index)
            if (n_admissible_evaluated >= 1_int32) then
                n_bins_per_point(1:n_points) = tmp_best_n_bins_per_point(1:n_points)
                shared_residual_range_low(1:n_points) = tmp_best_shared_residual_range_low(1:n_points)
                shared_residual_range_high(1:n_points) = tmp_best_shared_residual_range_high(1:n_points)
            else
                ! A single candidate that was rejected: Pass B never ran for it (or ran without
                ! passing the second gate), so m_min and a zero range, as in the fixed-k search.
                ! Its point count may exceed the capacity, hence the clamp.
                n_bins_per_point(1:min(n_points, max_n_points_candidate)) = actual_m_min
                shared_residual_range_low(1:min(n_points, max_n_points_candidate)) = 0.0_real64
                shared_residual_range_high(1:min(n_points, max_n_points_candidate)) = 0.0_real64
            end if
            plateau_established = logical(.true., kind=c_bool)
        else if (n_admissible_evaluated >= 1_int32) then
            k_start = candidates_k_start_k_step_k_max(1, best_uncertainty_candidate_index)
            k_step = candidates_k_start_k_step_k_max(2, best_uncertainty_candidate_index)
            k_max = candidates_k_start_k_step_k_max(3, best_uncertainty_candidate_index)
            n_points = candidate_n_points(best_uncertainty_candidate_index)
            n_bins_per_point(1:n_points) = tmp_best_uncertainty_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = tmp_best_uncertainty_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = tmp_best_uncertainty_shared_residual_range_high(1:n_points)
            best_candidate_confidence_interval = tmp_best_uncertainty_confidence_interval
            plateau_established = logical(.false., kind=c_bool)
        else
            ! Nothing admissible: candidate 1, whose point count may exceed the capacity.
            k_start = candidates_k_start_k_step_k_max(1, 1)
            k_step = candidates_k_start_k_step_k_max(2, 1)
            k_max = candidates_k_start_k_step_k_max(3, 1)
            n_points = candidate_n_points(1)
            n_bins_per_point(1:min(n_points, max_n_points_candidate)) = actual_m_min
            shared_residual_range_low(1:min(n_points, max_n_points_candidate)) = 0.0_real64
            shared_residual_range_high(1:min(n_points, max_n_points_candidate)) = 0.0_real64
            best_candidate_confidence_interval = -1.0_real64
            plateau_established = logical(.false., kind=c_bool)
        end if
    end subroutine run_js_comp_test_adaptive_parameter_search_impl

end module tox_data_integration_js_comp_test_impl
