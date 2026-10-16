#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The global JS-Comp-Test of K studies, with a data-driven choice of its neighborhood parameters.
!|
!| Fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each):
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| searches the `(n_points, n_neighbors)` candidates of
!| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| for a bootstrap plateau, and
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] runs the
!| final test (histograms, JSDs, weighting, permutation test) for the chosen pair.
!|
!| Adaptive neighborhoods, grown per reference point while the residual dispersion stays stable:
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| searches the `(k_start, k_step, k_max)` candidates of
!| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| builds the neighborhoods for the chosen knobs, and
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]] runs
!| the final test on them.
!|
!| Each reference point's histogram bin count is chosen by the occupancy search
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
!| (default) or by Stone's cross-validation
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]
!| (`bin_count_mode`). The building blocks are published too: the admissibility gates
!| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]], the
!| plateau checks
!| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]],
!| the bootstrap
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]], and the
!| brute-force reference
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds`, `calc_adaptive_js_comp_test_bounds` and
!| `calc_js_comp_test_n_top_k_jsds` size the work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan
    use f42_math_impl, only: clamp, is_close, LOG_2
    use f42_stats_impl, only: calc_quantile_impl, calc_sorted_slice_mad, histogram_bin_counts, sorted_bin_counts
    use f42_sort_impl, only: init_perm, sort_array_heapsort, sort_real_heapsort_expl_size, binary_search_insertion, &
                             block_containing_index, top_k_heap_push, bottom_k_heap_push
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
             determine_bin_count_stone_impl, group_sorted_perm_by_study, &
             MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, &
             MODE_BIN_COUNT_STONE_POOLED, &
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
             construct_adaptive_neighborhoods_impl, ADAPTIVE_MIN_VALID_RESIDUALS, &
             materialize_pooled_neighborhood, run_js_comp_test_adaptive_impl, &
             ADAPTIVE_STOP_TAU, ADAPTIVE_STOP_K_MAX, ADAPTIVE_STOP_RESIDUAL_CAP, ADAPTIVE_STOP_EXHAUSTED, &
             ADAPTIVE_STOP_ZERO_DISPERSION, ADAPTIVE_STOP_TOO_FEW_RESIDUALS, ADAPTIVE_STATUS_OK, &
             ADAPTIVE_STATUS_TOO_FEW_MEANS, ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD, ADAPTIVE_STATUS_TOO_FEW_RESIDUALS, &
             run_js_comp_test_adaptive_parameter_search_impl, ADAPTIVE_CANDIDATE_EVALUATED, &
             ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED, ADAPTIVE_CANDIDATE_EMPTY_STUDY, &
             ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED, ADAPTIVE_CANDIDATE_OVERLAP_FAILED, &
             ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED

    ! `join_method`'s mode table; the `METHOD_` prefix follows from the argument name `join_method`.
    integer(int32), parameter :: METHOD_JOIN_MIN = 0_int32
        !! Join method: succeeds only once every study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MAX = 1_int32
        !! Join method: succeeds once any one study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MEDIAN = 2_int32
        !! Join method: succeeds once a majority of studies' confidence-interval overlaps exceed
        !! the threshold, in [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]

    ! Tuning constants of the fixed-k candidate grid (generate_js_comp_test_candidates_impl).
    real(real64), parameter :: KX_FACTORS(2) = [0.25_real64, 0.5_real64]
        !! Neighbor-count factors, ascending: the first gives the larger `n_neighbors` candidate
    integer(int32), parameter :: MAX_POINTS = 1500_int32
    integer(int32), parameter :: MIN_POINTS = 300_int32
    real(real64), parameter :: GAMMA = 0.8_real64
    integer(int32), parameter :: MAX_POINT_CANDIDATES = 8_int32
        !! How often n_points_high may be multiplied by GAMMA before the grid gives up
    ! CM_MAX_CANDIDATE_PAIRS and CM_MAX_N_BINS are the single source of these two values, used
    ! both for the parameters and for the dummy extents of wrapped routines: a parameter cannot be
    ! such an extent, because the generated C and R bindings need the number itself.
#define CM_MAX_CANDIDATE_PAIRS 16
#define CM_MAX_N_BINS 256
    integer(int32), parameter :: MAX_CANDIDATE_PAIRS = CM_MAX_CANDIDATE_PAIRS
        !! = size(KX_FACTORS) * MAX_POINT_CANDIDATES = 2 * 8
    integer(int32), parameter :: MAX_N_BINS = CM_MAX_N_BINS
        !! Fixed ceiling on any histogram bin count; histogram work arrays are sized by it

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

    ! Output codes of construct_adaptive_neighborhoods_impl (not mode arguments). Each is a CM_
    ! macro first, so the parameter and the documentation quoting it share one literal.
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

    ! Defaults of construct_adaptive_neighborhoods_impl's optional knobs (shared by DM_DEFAULT and
    ! M_DEFAULT_VAL).
#define CM_ADAPTIVE_TAU_DEFAULT 0.1_real64
#define CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT 1.0_real64
#define CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT 0_int32
#define CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT 1_int32

    ! Per-candidate status codes of run_js_comp_test_adaptive_parameter_search_impl's log
    ! (output codes, CM_ macros first like the ones above).
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

    ! `plateau_mode`'s mode table. CM_MODE_PLATEAU_CI_OVERLAP is also its DM_DEFAULT, which cannot
    ! name a Fortran parameter.
#define CM_MODE_PLATEAU_CI_OVERLAP 0_int32
    integer(int32), parameter :: MODE_PLATEAU_CI_OVERLAP = CM_MODE_PLATEAU_CI_OVERLAP
        !! Plateau mode: CI overlap only
    integer(int32), parameter :: MODE_PLATEAU_EFFECT_SIZE = 1_int32
        !! Plateau mode: relative-effect-size stability only
    integer(int32), parameter :: MODE_PLATEAU_BOTH = 2_int32
        !! Plateau mode: either CI overlap or relative-effect-size stability

    ! `bin_count_mode`'s values: how a reference point's histogram bin count is chosen. Occupancy
    ! is the default of all four pipeline entry points; determine_bin_count_stone_impl rejects it.
    ! CM_ macros first, so that a DM_DEFAULT can name them.
#define CM_MODE_BIN_COUNT_OCCUPANCY 0_int32
#define CM_MODE_BIN_COUNT_STONE_EQUAL 1_int32
#define CM_MODE_BIN_COUNT_STONE_SIZE_WEIGHTED 2_int32
#define CM_MODE_BIN_COUNT_STONE_POOLED 3_int32
    integer(int32), parameter :: MODE_BIN_COUNT_OCCUPANCY = CM_MODE_BIN_COUNT_OCCUPANCY
        !! Bin-count mode: the largest bin count whose pooled histogram has every bin at or above
        !! the occupancy floor
    integer(int32), parameter :: MODE_BIN_COUNT_STONE_EQUAL = CM_MODE_BIN_COUNT_STONE_EQUAL
        !! Bin-count mode: Stone's cross-validation score of every study, summed with equal weights
    integer(int32), parameter :: MODE_BIN_COUNT_STONE_SIZE_WEIGHTED = CM_MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        !! Bin-count mode: Stone's cross-validation score of every study, weighted by the study's
        !! residual count inside the range
    integer(int32), parameter :: MODE_BIN_COUNT_STONE_POOLED = CM_MODE_BIN_COUNT_STONE_POOLED
        !! Bin-count mode: Stone's cross-validation score of the pooled residuals of all studies

    ! Defaults of the relative-effect-size plateau criterion (provisional, to be validated
    ! empirically), shared by DM_DEFAULT and M_DEFAULT_VAL.
#define CM_DELTA_MEDIAN_THRESHOLD_DEFAULT 0.05_real64
#define CM_DELTA_MAX_THRESHOLD_DEFAULT 0.10_real64
#define CM_DELTA_EPSILON_DEFAULT 1.0e-10_real64
#define CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT 2_int32

    ! Defaults of the occupancy bin-count search (determine_bin_count_occupancy_impl).
    ! CM_OCCUPANCY_GAMMA_DEFAULT (growth of a candidate bin count) is unrelated to GAMMA (decay of
    ! the fixed-k grid's n_points).
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
    !| For the candidate grid of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]:
    !| `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    !| `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`.
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
    !| Sizes the per-point work arrays of the adaptive parameter search over the candidates of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| A candidate's reference-point count is only known after its neighborhoods are grown, so this
    !| is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
    !| `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
    !| candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.
    !|
    !| Rejects, with invalid input naming neither argument when the bound is joint, a
    !| `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
    !| first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
    !| `huge(1_int32)/4 = 536870911` (so that its `k_max = 4*k_start_1` fits); the output is then
    !| left undefined.
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
        ! extents derive from it), and so must the first candidate's k_max = 4*k_start_1.
        if (n_pool > int(huge(1_int32), int64)) call set_err_once(ierr, ERR_INVALID_INPUT)
        if (adaptive_k_start_1(max_n_genes_all_studies, n_studies) > int(huge(1_int32)/4_int32, int64)) &
            call set_err_once(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        call compute_adaptive_k_starts(max_n_genes_all_studies, n_studies, k_starts, n_k_starts)

        capacity = ceiling(ADAPTIVE_POINT_CAPACITY_FACTOR*real(n_pool, real64)/ &
                           real(k_starts(n_k_starts), real64), kind=int64) + 1_int64
        max_n_points_candidate = int(min(n_pool, capacity), kind=int32)
    end subroutine calc_adaptive_js_comp_test_bounds

    !> summary: Estimate the histogram bin count for one (n_points, n_neighbors) candidate
    !| AUTHOR_LASZLO_LANG
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] and returned as
    !| `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    !| `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    !| from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    !| bin count.
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
            !! to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN

        integer(int32) :: i_pool, n_pool, sturges_raw, fd_raw
        real(real64) :: half_bin_width, quartile_25, quartile_75, n_reps_neighborhood, fd_ratio

        ! NaN sorts last under the ascending permutation: scan back for the last non-NaN position.
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

        ! Sturges, with 1 + nint(log2 n) rather than the textbook ceil(1 + log2 n) (diagnostic only).
        sturges_raw = 1_int32 + nint(log(n_reps_neighborhood)/LOG_2, kind=int32)
        sturges_bins = max(1_int32, min(sturges_raw, MAX_N_BINS))

        ! Freedman-Diaconis, with the width not doubled: R is half of the [-R, R] range already.
        call calc_quantile_impl(residuals, n_residuals, residuals_perm, 0.25_real64, quartile_25, n_considered=n_pool)
        call calc_quantile_impl(residuals, n_residuals, residuals_perm, 0.75_real64, quartile_75, n_considered=n_pool)

        half_bin_width = (quartile_75 - quartile_25)/(n_reps_neighborhood**(1.0_real64/3.0_real64))
        if (.not. is_close(half_bin_width, 0.0_real64)) then
            ! A tiny interquartile range pushes the ratio past the int32 range, where nint is
            ! processor-dependent: clamp the real ratio first (+Infinity too). A NaN ratio, which min
            ! would pass through unpredictably, falls back to Sturges.
            fd_ratio = shared_residual_range/half_bin_width
            if (ieee_is_nan(fd_ratio)) then
                fd_raw = sturges_raw
            else
                fd_raw = nint(min(fd_ratio, real(MAX_N_BINS, kind=real64)), kind=int32)
            end if
        else
            fd_raw = sturges_raw
        end if
        fd_bins = max(1_int32, min(fd_raw, MAX_N_BINS))

        n_bins = max(sturges_bins, fd_bins)
    end subroutine estimate_bin_count_impl

    !> summary: Determine one neighborhood's occupancy-constrained histogram bin count (Issue #187)
    !| AUTHOR_LASZLO_LANG
    !| For one neighborhood's pooled residuals (all its neighbors, all studies), finds the largest
    !| bin count `M` in `[m_min, m_max]` whose equal-width histogram over the neighborhood's own range
    !| `[shared_residual_range_low, shared_residual_range_high]` (the
    !| `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of these residuals)
    !| has every bin at or above `min_residuals_per_bin`. A geometric search (factor
    !| `gamma_occupancy`) brackets the first failing `M`, then every integer below it is tested,
    !| because occupancy is not monotonic in `M`.
    !|
    !| `occupancy_failed = .true.` (even `m_min` bins fail, or every residual is NaN) means the caller
    !| should reject the neighborhood; `selected_n_bins` is then `m_min` only so that it is defined.
    !| `sturges_bins`/`fd_bins` come from
    !| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]] and are
    !| diagnostics only.
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
            !! `.true.` iff even m_min bins fail the occupancy criterion (or every pooled residual is
            !! NaN); the caller should then reject this neighborhood
        real(real64), intent(out) :: shared_residual_range_low
            !! This neighborhood's lower residual-range bound (R_low): the
            !! lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        real(real64), intent(out) :: shared_residual_range_high
            !! This neighborhood's upper residual-range bound (R_high): the
            !! upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
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
            !! Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(int32), intent(out) :: tmp_bin_counts(CM_MAX_N_BINS)
            !! Working array: per-bin counts of the candidate bin count under test
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
            !! smaller
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! DM_MIN(lower_residual_range_quantile)
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

        ! NaN sorts last under the ascending permutation: scan back for the last non-NaN position.
        n_pooled_residuals = n_residuals
        do i_pool = n_pooled_residuals, 1, -1
            if (ieee_is_nan(pooled_residuals(pooled_residuals_perm(i_pool)))) then
                n_pooled_residuals = n_pooled_residuals - 1_int32
            else
                exit
            end if
        end do

        ! This neighborhood's own range; 0.0/0.0 when every pooled residual is NaN.
        if (n_pooled_residuals == 0_int32) then
            shared_residual_range_low = 0.0_real64
            shared_residual_range_high = 0.0_real64
        else
            call calc_quantile_impl(pooled_residuals, n_residuals, pooled_residuals_perm, &
                                      actual_lower_residual_range_quantile, shared_residual_range_low, &
                                      n_considered=n_pooled_residuals)
            call calc_quantile_impl(pooled_residuals, n_residuals, pooled_residuals_perm, &
                                      actual_upper_residual_range_quantile, shared_residual_range_high, &
                                      n_considered=n_pooled_residuals)
        end if

        ! Diagnostics only; before any early return so sturges_bins/fd_bins are always defined
        ! (estimate_bin_count_impl handles an all-NaN pool). It expects the one-sided R, so half the span.
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

        ! Stage 1: geometric coarse search, recording the largest admissible m_valid and the first
        ! inadmissible m_invalid. Sequential: each step depends on the previous outcome.
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
                ! Still admissible at m_max: nothing left for stage 2.
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
            ! Even m_min fails.
            occupancy_failed = .true.
            selected_n_bins = actual_m_min
            min_bin_occupancy = 0_int32
            mean_bin_occupancy = 0.0_real64
            max_bin_occupancy = 0_int32
            return
        end if

        ! Stage 2: test every integer in (m_valid, m_invalid); not a binary search, because
        ! occupancy is not monotonic in M once bin boundaries are recomputed per candidate.
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

    !> summary: Exhaustive brute-force reference implementation of Issue #187's occupancy search
    !| AUTHOR_LASZLO_LANG
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]],
    !| whose faster search could in principle miss a larger admissible `M`. Takes the range and
    !| `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
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
            !! `.true.` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion
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
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
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
        ! Candidates are independent, so concurrent; reduce(max) has no argmax, so best_m's
        ! occupancy is recomputed once afterwards.
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

    !> Stone's leave-one-out cross-validation score (Duke STA 293.3 lecture notes, ch. 6,
    !| Theorem 2) of one histogram of `n_inside` values over `trial_m` equal bins on a range of
    !| width `range_width`:
    !| `J = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2`, with `h = range_width/trial_m` and
    !| `p_j = bin_counts(j)/n`, evaluated as `J = trial_m*(2n^2 - (n+1)*Q) / (range_width*(n-1)*n^2)`
    !| with the integer `Q = sum_j bin_counts(j)^2`: while `n^3 < 2^53` the numerator is exact, so
    !| mathematically equal scores are bit-equal and the tie rule does not depend on rounding.
    !| Precondition, not checked: `n_inside >= 2` and `range_width > 0`. Private helper of
    !| `stone_bin_count_core`.
    pure real(real64) function stone_risk_estimate(bin_counts, trial_m, n_inside, range_width) result(risk)
        integer(int32), intent(in) :: trial_m
            !! Number of bins
        integer(int32), intent(in) :: bin_counts(trial_m)
            !! Count of values per bin
        integer(int32), intent(in) :: n_inside
            !! Number of values, the sum of `bin_counts`
        real(real64), intent(in) :: range_width
            !! Width of the binned range, R_high - R_low

        integer(int64) :: sum_sq
        integer(int32) :: i_bin
        real(real64) :: n_real

        sum_sq = 0_int64
        do concurrent(i_bin=1:trial_m) shared(bin_counts) reduce(+:sum_sq)
            sum_sq = sum_sq + int(bin_counts(i_bin), int64)*int(bin_counts(i_bin), int64)
        end do
        n_real = real(n_inside, real64)
        risk = real(trial_m, real64)*(2.0_real64*n_real*n_real - (n_real + 1.0_real64)*real(sum_sq, real64)) &
               /(range_width*(n_real - 1.0_real64)*n_real*n_real)
    end function stone_risk_estimate

    !> The Stone bin-count search on pooled residuals that are already grouped by study.
    !|
    !| Private core of
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]],
    !| also used by callers whose pool already holds every study as one contiguous block. Inputs:
    !|
    !| - `pooled_residuals_perm` orders all of `pooled_residuals` ascending, NaN last.
    !| - `residual_to_study_perm` holds indices into `pooled_residuals`, grouped by study (study `k`
    !|   fills positions `sum(study_block_sizes(1:k-1)) + 1` to `sum(study_block_sizes(1:k))`), each
    !|   block ascending with NaN last; `group_sorted_perm_by_study` builds it.
    !|
    !| The range is as in
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]].
    !| Every bin count `m` in `[m_min, m_max]` is tested in ascending order (the score has many
    !| local minima):
    !|
    !| 1. Each study's residuals inside the range are counted into `m` bins
    !|    (`sorted_bin_counts`). Values outside the range are left out of the score: with them
    !|    clamped into the edge bins, Stone's score keeps falling with `m` and always picks `m_max`.
    !| 2. The score depends on `bin_count_mode`: the sum of every study's
    !|    `stone_risk_estimate` (`MODE_BIN_COUNT_STONE_EQUAL`), the same sum with each study's
    !|    score multiplied by its inside count (`MODE_BIN_COUNT_STONE_SIZE_WEIGHTED`), or one score
    !|    of the summed inside counts (`MODE_BIN_COUNT_STONE_POOLED`). Any other value is treated as
    !|    equal weights. The per-study modes skip a study with fewer than 2 inside residuals.
    !| 3. `m` is admissible only if the pooled clamped histogram (the summed inside counts, plus the
    !|    pooled residuals below the range in bin 1 and those above it in bin `m`) has every bin at
    !|    or above `min_residuals_per_bin`. That histogram is the one `histogram_bin_counts` builds.
    !| 4. The admissible `m` with the smallest score wins; on a tie the smaller `m` stays.
    !|
    !| The search fails (`stone_failed = .true.`, `selected_n_bins = m_min`, occupancy diagnostics
    !| 0) when every pooled residual is NaN, the range has width 0 or less, no study has 2 inside
    !| residuals (pooled mode: fewer than 2 inside residuals in total), or no `m` is admissible.
    !| Otherwise the occupancy diagnostics describe the pooled clamped histogram at the chosen `m`.
    pure subroutine stone_bin_count_core(pooled_residuals, pooled_residuals_perm, n_residuals, n_studies, &
                                         residual_to_study_perm, study_block_sizes, max_n_reps_all_studies, n_neighbors, &
                                         bin_count_mode, selected_n_bins, stone_failed, shared_residual_range_low, &
                                         shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                         mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, tmp_bin_counts, &
                                         m_min, m_max, min_residuals_per_bin, lower_residual_range_quantile, &
                                         upper_residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! Pooled signed residuals of one neighborhood, all studies
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: residual_to_study_perm(n_residuals)
            !! Indices into `pooled_residuals`, grouped by study and sorted within each study (see above)
        integer(int32), intent(in) :: study_block_sizes(n_studies)
            !! Number of residuals of each study, NaN included; sums to `n_residuals`
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies, for the Sturges/FD diagnostics
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size, for the Sturges/FD diagnostics
        integer(int32), intent(in) :: bin_count_mode
            !! MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED or MODE_BIN_COUNT_STONE_POOLED
        integer(int32), intent(out) :: selected_n_bins
            !! The chosen bin count; m_min when stone_failed
        logical(c_bool), intent(out) :: stone_failed
            !! `.true.` iff the search failed (see above)
        real(real64), intent(out) :: shared_residual_range_low
            !! Lower residual-range bound (R_low); 0.0 when every pooled residual is NaN
        real(real64), intent(out) :: shared_residual_range_high
            !! Upper residual-range bound (R_high); 0.0 when every pooled residual is NaN
        integer(int32), intent(out) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        real(real64), intent(out) :: mean_bin_occupancy
            !! n_pooled_residuals / selected_n_bins; 0 when stone_failed
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule diagnostic from estimate_bin_count_impl
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule diagnostic from estimate_bin_count_impl
        integer(int32), intent(out) :: tmp_bin_counts(MAX_N_BINS)
            !! Working array: one study's inside counts at the bin count under test
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count tested; default CM_OCCUPANCY_M_MIN_DEFAULT
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count tested, raised to m_min if smaller; default CM_OCCUPANCY_M_MAX_DEFAULT
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Floor on every bin of the pooled clamped histogram; default
            !! CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile of R_low; default CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile of R_high; default CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        real(real64) :: actual_lower_residual_range_quantile, actual_upper_residual_range_quantile
        integer(int32) :: i_pool, i_study, block_start, block_end, trial_m, best_m, best_min_occ, best_max_occ
        integer(int32) :: n_eligible_studies, n_total_inside, n_study_below, n_study_inside
        integer(int32) :: pooled_n_below, pooled_n_inside, pooled_n_above, min_occ, max_occ, combined_n_bins_unused
        integer(int32) :: pooled_counts(MAX_N_BINS)
        real(real64) :: half_span, range_width, score, best_score
        logical(c_bool) :: pooled_scoring, size_weighted, found_admissible

        M_DEFAULT_VAL(m_min, actual_m_min, CM_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, CM_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(lower_residual_range_quantile, actual_lower_residual_range_quantile, CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        M_DEFAULT_VAL(upper_residual_range_quantile, actual_upper_residual_range_quantile, CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        actual_m_max = max(actual_m_min, actual_m_max)
        pooled_scoring = bin_count_mode == MODE_BIN_COUNT_STONE_POOLED
        size_weighted = bin_count_mode == MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

        ! NaN sorts last -- same scan-back as determine_bin_count_occupancy_impl.
        n_pooled_residuals = n_residuals
        do i_pool = n_pooled_residuals, 1, -1
            if (ieee_is_nan(pooled_residuals(pooled_residuals_perm(i_pool)))) then
                n_pooled_residuals = n_pooled_residuals - 1_int32
            else
                exit
            end if
        end do

        ! The range, exactly as in determine_bin_count_occupancy_impl (0.0/0.0 for an all-NaN pool).
        if (n_pooled_residuals == 0_int32) then
            shared_residual_range_low = 0.0_real64
            shared_residual_range_high = 0.0_real64
        else
            call calc_quantile_impl(pooled_residuals, n_residuals, pooled_residuals_perm, &
                                      actual_lower_residual_range_quantile, shared_residual_range_low, &
                                      n_considered=n_pooled_residuals)
            call calc_quantile_impl(pooled_residuals, n_residuals, pooled_residuals_perm, &
                                      actual_upper_residual_range_quantile, shared_residual_range_high, &
                                      n_considered=n_pooled_residuals)
        end if

        ! Diagnostics only, before any early return so sturges_bins/fd_bins are always defined;
        ! half the span, as in determine_bin_count_occupancy_impl.
        half_span = (shared_residual_range_high - shared_residual_range_low)/2.0_real64
        call estimate_bin_count_impl(pooled_residuals, pooled_residuals_perm, n_residuals, max_n_reps_all_studies, &
                                     n_neighbors, half_span, combined_n_bins_unused, sturges_bins, fd_bins)

        ! The failure outputs; overwritten below once an admissible bin count is found.
        stone_failed = .true.
        selected_n_bins = actual_m_min
        min_bin_occupancy = 0_int32
        mean_bin_occupancy = 0.0_real64
        max_bin_occupancy = 0_int32

        range_width = shared_residual_range_high - shared_residual_range_low
        if (n_pooled_residuals == 0_int32 .or. range_width <= 0.0_real64) return

        ! The pooled residuals below and above the range, which clamping adds to the edge bins.
        call sorted_bin_counts(pooled_residuals, n_residuals, pooled_residuals_perm(1:n_pooled_residuals), &
                               n_pooled_residuals, shared_residual_range_low, shared_residual_range_high, 1_int32, &
                               tmp_bin_counts(1:1), pooled_n_below, pooled_n_inside)
        pooled_n_above = n_pooled_residuals - pooled_n_below - pooled_n_inside

        ! Which studies can be scored: each needs 2 inside residuals (pooled mode: 2 in total).
        n_eligible_studies = 0_int32
        n_total_inside = 0_int32
        block_end = 0_int32
        do i_study = 1, n_studies
            block_start = block_end + 1_int32
            block_end = block_end + study_block_sizes(i_study)
            ! An empty study has nothing inside; skipping it also keeps an empty slice that starts
            ! past the array's end away from ifx's bounds check, which rejects it.
            if (study_block_sizes(i_study) == 0_int32) cycle
            call sorted_bin_counts(pooled_residuals, n_residuals, residual_to_study_perm(block_start:block_end), &
                                   study_block_sizes(i_study), shared_residual_range_low, shared_residual_range_high, &
                                   1_int32, tmp_bin_counts(1:1), n_study_below, n_study_inside)
            n_total_inside = n_total_inside + n_study_inside
            if (n_study_inside >= 2_int32) n_eligible_studies = n_eligible_studies + 1_int32
        end do
        if (pooled_scoring) then
            if (n_total_inside < 2_int32) return
        else
            if (n_eligible_studies == 0_int32) return
        end if

        found_admissible = .false.
        best_m = actual_m_min
        best_score = 0.0_real64
        best_min_occ = 0_int32
        best_max_occ = 0_int32
        ! Sequential: a running argmin whose tie rule depends on the order.
        do trial_m = actual_m_min, actual_m_max
            pooled_counts(1:trial_m) = 0_int32
            score = 0.0_real64
            block_end = 0_int32
            do i_study = 1, n_studies
                block_start = block_end + 1_int32
                block_end = block_end + study_block_sizes(i_study)
                ! An empty study adds no counts and no score (see the eligibility loop above).
                if (study_block_sizes(i_study) == 0_int32) cycle
                call sorted_bin_counts(pooled_residuals, n_residuals, residual_to_study_perm(block_start:block_end), &
                                       study_block_sizes(i_study), shared_residual_range_low, &
                                       shared_residual_range_high, trial_m, tmp_bin_counts(1:trial_m), &
                                       n_study_below, n_study_inside)
                pooled_counts(1:trial_m) = pooled_counts(1:trial_m) + tmp_bin_counts(1:trial_m)
                if (.not. pooled_scoring .and. n_study_inside >= 2_int32) then
                    if (size_weighted) then
                        score = score + real(n_study_inside, real64)* &
                                stone_risk_estimate(tmp_bin_counts(1:trial_m), trial_m, n_study_inside, range_width)
                    else
                        score = score + stone_risk_estimate(tmp_bin_counts(1:trial_m), trial_m, n_study_inside, &
                                                            range_width)
                    end if
                end if
            end do
            if (pooled_scoring) score = stone_risk_estimate(pooled_counts(1:trial_m), trial_m, n_total_inside, range_width)

            ! The floor, on the pooled clamped histogram.
            pooled_counts(1) = pooled_counts(1) + pooled_n_below
            pooled_counts(trial_m) = pooled_counts(trial_m) + pooled_n_above
            min_occ = minval(pooled_counts(1:trial_m))
            max_occ = maxval(pooled_counts(1:trial_m))
            if (min_occ < actual_min_residuals_per_bin) cycle

            ! Strictly smaller only, so on a tie the smaller bin count (tested first) stays.
            if (.not. found_admissible .or. score < best_score) then
                found_admissible = .true.
                best_score = score
                best_m = trial_m
                best_min_occ = min_occ
                best_max_occ = max_occ
            end if
        end do

        if (.not. found_admissible) return

        stone_failed = .false.
        selected_n_bins = best_m
        min_bin_occupancy = best_min_occ
        max_bin_occupancy = best_max_occ
        mean_bin_occupancy = real(n_pooled_residuals, real64)/real(best_m, real64)
    end subroutine stone_bin_count_core

    !> summary: Group an ascending permutation of pooled residuals by study in one pass, each study's block sorted
    !| AUTHOR_LASZLO_LANG
    !| Walks `pooled_residuals_perm` (ascending, NaN last) once and appends each index to its
    !| study's block of `residual_to_study_perm`, so every block comes out ascending, NaN last,
    !| without a sort per block (equal residuals may come out in another order, which changes no
    !| bin count). The study of an index comes from exactly one of the two optional arguments:
    !|
    !| - `residual_to_study(idx)`, a label per pooled residual; or
    !| - `study_block_ends`, for a pool whose residuals are laid out as contiguous per-study blocks:
    !|   index `idx` belongs to the first study `s` with `study_block_ends(s) >= idx`. An empty block
    !|   has the same end as the block before it and receives no index.
    !|
    !| `study_block_cursors` is the write cursor of every block: on entry the last position before
    !| study `s`'s block (0 for study 1), on exit the last position of that block. With both
    !| optional arguments present the labels are used; with neither, nothing is grouped:
    !| `residual_to_study_perm` is all 0 and the cursors are unchanged. Precondition, not
    !| checked: one of the optional arguments is present, every label is in
    !| `[1, n_studies]` or `study_block_ends` ascends to `n_residuals`, and the cursors on entry
    !| match the block sizes. Not published; public only so the tests can compare it with a sort
    !| per block.
    pure subroutine group_sorted_perm_by_study(pooled_residuals_perm, n_residuals, n_studies, study_block_cursors, &
                                               residual_to_study_perm, residual_to_study, study_block_ends)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Indices of the pooled residuals, ordered ascending by residual, NaN last
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(inout) :: study_block_cursors(n_studies)
            !! On entry the last position before each study's block, on exit the last position of it
        integer(int32), intent(out) :: residual_to_study_perm(n_residuals)
            !! The indices of `pooled_residuals_perm`, grouped by study, each block in its order
        integer(int32), intent(in), optional :: residual_to_study(n_residuals)
            !! Study label of each pooled residual
        integer(int32), intent(in), optional :: study_block_ends(n_studies)
            !! Last pool index of each study's contiguous block

        integer(int32) :: i_pos, residual_idx, study_idx

        ! Sequential: residuals of one study share a cursor, and the walk order is the block order.
        if (present(residual_to_study)) then
            do i_pos = 1, n_residuals
                residual_idx = pooled_residuals_perm(i_pos)
                study_idx = residual_to_study(residual_idx)
                study_block_cursors(study_idx) = study_block_cursors(study_idx) + 1_int32
                residual_to_study_perm(study_block_cursors(study_idx)) = residual_idx
            end do
        else if (present(study_block_ends)) then
            do i_pos = 1, n_residuals
                residual_idx = pooled_residuals_perm(i_pos)
                study_idx = block_containing_index(study_block_ends, n_studies, residual_idx)
                study_block_cursors(study_idx) = study_block_cursors(study_idx) + 1_int32
                residual_to_study_perm(study_block_cursors(study_idx)) = residual_idx
            end do
        else
            ! Neither is present (a broken precondition): nothing is placed. The cursors stay as
            ! they were; the intent(out) residual_to_study_perm is defined as all 0.
            residual_to_study_perm = 0_int32
        end if
    end subroutine group_sorted_perm_by_study

    !> Groups pooled residuals by study and sorts each study's residuals: fills
    !| `residual_to_study_perm` with the indices of the pooled residuals grouped by
    !| `residual_to_study` (study 1 first), each study's indices ordered ascending by residual with
    !| NaN last, and `study_block_sizes` with each study's residual count. Private helper of
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]].
    !| Precondition, not checked: every label is in `[1, n_studies]`, and `pooled_residuals_perm`
    !| orders the pooled residuals ascending, NaN last.
    pure subroutine build_residual_to_study_perm(pooled_residuals_perm, n_residuals, residual_to_study, n_studies, &
                                                 residual_to_study_perm, study_block_sizes)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Indices of the pooled residuals, ordered ascending by residual, NaN last
        integer(int32), intent(in) :: residual_to_study(n_residuals)
            !! Study label of each pooled residual, in `[1, n_studies]`
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(out) :: residual_to_study_perm(n_residuals)
            !! Indices of the pooled residuals, grouped by study and sorted within each study
        integer(int32), intent(out) :: study_block_sizes(n_studies)
            !! Number of pooled residuals of each study, NaN included

        integer(int32) :: i_residual, i_study, study_idx, n_in_study, block_end

        ! Count per study label. Sequential: distinct residuals can share a counter.
        study_block_sizes = 0_int32
        do i_residual = 1, n_residuals
            study_idx = residual_to_study(i_residual)
            study_block_sizes(study_idx) = study_block_sizes(study_idx) + 1_int32
        end do
        ! Each counter becomes the position just before its study's block (the write cursor) ...
        block_end = 0_int32
        do i_study = 1, n_studies
            n_in_study = study_block_sizes(i_study)
            study_block_sizes(i_study) = block_end
            block_end = block_end + n_in_study
        end do
        ! ... and, once every index is placed, the last position of that block.
        call group_sorted_perm_by_study(pooled_residuals_perm, n_residuals, n_studies, study_block_sizes, &
                                        residual_to_study_perm, residual_to_study=residual_to_study)

        ! Block ends back to block sizes. Sequential, from the last study down, so every
        ! study still reads its predecessor's end.
        do i_study = n_studies, 2, -1
            study_block_sizes(i_study) = study_block_sizes(i_study) - study_block_sizes(i_study - 1_int32)
        end do
    end subroutine build_residual_to_study_perm

    !> summary: Choose one neighborhood's histogram bin count by Stone's cross-validation rule
    !| AUTHOR_LASZLO_LANG
    !| Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
    !| Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
    !| `n` values over `m` bins of width `h`, the bin count minimising
    !| `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
    !| chosen. An alternative to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]],
    !| with the same range, floor and diagnostics.
    !|
    !| The range is the `lower_residual_range_quantile`/`upper_residual_range_quantile` quantiles of
    !| all non-NaN pooled residuals, so `h = (R_high - R_low)/m`. The score counts only the residuals
    !| **inside** the range (`n` is their number). Every `m` in `[m_min, m_max]` is tested, ties
    !| going to the smaller `m`. A single score (pooled mode, or one study) is bit-exact for up to
    !| about 200,000 residuals, so ties are exact; a sum over studies may differ in the last bit.
    !|
    !| `residual_to_study` says which study each residual comes from, and `bin_count_mode` how the
    !| studies are combined into one shared bin count:
    !| [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] (the
    !| default) minimises the sum of every study's own score `J_k(m)`,
    !| [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]]
    !| the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
    !| [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] the
    !| score `J(m)` of all studies' residuals together. The two per-study modes leave out a study
    !| with fewer than 2 residuals inside the range.
    !|
    !| `min_residuals_per_bin` stays a hard floor, with the occupancy search's definition: a bin
    !| count is admissible only if the pooled histogram with the outside residuals clamped into the
    !| edge bins has every bin at or above it (`0` turns the floor off). The minimum is taken over
    !| the admissible bin counts only.
    !|
    !| `stone_failed = .true.` means no bin count could be chosen: every pooled residual is NaN, the
    !| range has width 0 or less, no study has 2 residuals inside the range (pooled mode: fewer than
    !| 2 inside in total), or no bin count passes the floor. `selected_n_bins` is then `m_min`, only
    !| so a caller ignoring the flag has a value to build with, and the occupancy diagnostics are 0.
    !| Otherwise the diagnostics describe the pooled clamped histogram at `selected_n_bins`.
    pure subroutine determine_bin_count_stone_impl(pooled_residuals, pooled_residuals_perm, n_residuals, &
                                                    residual_to_study, n_studies, max_n_reps_all_studies, n_neighbors, &
                                                    selected_n_bins, stone_failed, shared_residual_range_low, &
                                                    shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                                    mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                                    tmp_bin_counts, tmp_residual_to_study_perm, tmp_study_block_sizes, &
                                                    bin_count_mode, m_min, m_max, min_residuals_per_bin, &
                                                    lower_residual_range_quantile, &
                                                    upper_residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! DM_MIN(1_int32)
            !! DM_MAX(n_residuals)
        integer(int32), intent(in) :: residual_to_study(n_residuals)
            !! Study each pooled residual comes from
            !! DM_MIN(1_int32)
            !! DM_MAX(n_studies)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: selected_n_bins
            !! The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
            !! Stone score; m_min when stone_failed
        logical(c_bool), intent(out) :: stone_failed
            !! `.true.` iff no bin count could be chosen (see above); the caller should then reject
            !! this neighborhood rather than build a histogram from selected_n_bins
        real(real64), intent(out) :: shared_residual_range_low
            !! This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN
        real(real64), intent(out) :: shared_residual_range_high
            !! This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN
        integer(int32), intent(out) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        real(real64), intent(out) :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            !! stone_failed
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(int32), intent(out) :: tmp_bin_counts(CM_MAX_N_BINS)
            !! Working array: one study's bin counts at the bin count under test
        integer(int32), intent(out) :: tmp_residual_to_study_perm(n_residuals)
            !! Working array: indices of `pooled_residuals` grouped by study, each study sorted
        integer(int32), intent(out) :: tmp_study_block_sizes(n_studies)
            !! Working array: number of pooled residuals of each study
        integer(int32), intent(in), optional :: bin_count_mode
            !! How the studies' scores are combined into one bin count
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Minimise the sum of the per-study scores (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(CM_MODE_BIN_COUNT_STONE_EQUAL)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin of the pooled clamped histogram must
            !! reach for a candidate bin count to be admissible (0 turns the floor off)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)

        integer(int32) :: actual_bin_count_mode

        M_DEFAULT_VAL(bin_count_mode, actual_bin_count_mode, CM_MODE_BIN_COUNT_STONE_EQUAL)

        call build_residual_to_study_perm(pooled_residuals_perm, n_residuals, residual_to_study, n_studies, &
                                          tmp_residual_to_study_perm, tmp_study_block_sizes)

        call stone_bin_count_core(pooled_residuals, pooled_residuals_perm, n_residuals, n_studies, &
                                  tmp_residual_to_study_perm, tmp_study_block_sizes, max_n_reps_all_studies, n_neighbors, &
                                  actual_bin_count_mode, selected_n_bins, stone_failed, shared_residual_range_low, &
                                  shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                  max_bin_occupancy, sturges_bins, fd_bins, tmp_bin_counts, m_min, m_max, &
                                  min_residuals_per_bin, lower_residual_range_quantile, upper_residual_range_quantile)
    end subroutine determine_bin_count_stone_impl

    !> M_EXPORT_C
    !| summary: Pool one reference point's residuals across every neighbor and every study
    !| AUTHOR_LASZLO_LANG
    !| Gathers one reference point's residuals from every neighbor gene of every study into one flat
    !| array, given that point's per-study neighbor gene indices (one column of the
    !| `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
    !| `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| feed to the bin-count search, so a caller can pass the same input to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
    !| on real data.
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
            !! Matrix of signed residuals per study, NaN allowed for missing values; not validated
        integer(int32), intent(in) :: neighborhood_indices_point(n_neighbors, n_studies)
            !! Gene indices of one reference point's neighborhood, per study; each in
            !! [1, max_n_genes_all_studies]
        real(real64), intent(out) :: pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies)
            !! This reference point's pooled residuals, laid out as a
            !! (max_n_reps_all_studies, n_neighbors, n_studies) array
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: i_study, i_neighbor, i_rep, gene_idx, n_predecessors

        call set_ok(ierr)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return

        ! Every gene index subscripts residuals directly, so a bad one would read out of bounds.
        ! Sequential: set_err_once on the shared ierr is a read-then-write, not allowed in do concurrent.
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

    !> Chooses one reference point's bin count from its gathered, sorted pool with the rule
    !| `bin_count_mode` names.
    !|
    !| The pool `pooled_residuals(1:n_residuals)` holds every study's residuals as one contiguous
    !| block, study after study: study `s` fills pool indices `study_block_ends(s-1) + 1` to
    !| `study_block_ends(s)` (from 1 for study 1). `pooled_residuals_perm` orders the whole pool
    !| ascending, NaN last.
    !|
    !| - [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]:
    !|   calls
    !|   [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !|   with exactly these arguments, so the result is the occupancy search's, bit for bit. The
    !|   two per-study work arrays are not touched.
    !| - The three Stone modes: `group_sorted_perm_by_study` groups the sorted pool by block into
    !|   `tmp_residual_to_study_perm`, the block sizes go into `tmp_study_block_sizes`, and
    !|   `stone_bin_count_core` chooses the bin count. Its `stone_failed` is returned as
    !|   `bin_count_failed`. `gamma_occupancy` is not used.
    !|
    !| Both rules get `n_neighbors` only for the Sturges/Freedman-Diaconis diagnostics. Private
    !| helper of `determine_point_bin_count` and of the fixed-k parameter search.
    !| Precondition, not checked: `bin_count_mode` is one of
    !| the four `MODE_BIN_COUNT_*` values (any other is treated as a Stone mode with equal weights),
    !| and `study_block_ends` ascends to `n_residuals`.
    pure subroutine select_point_bin_count(pooled_residuals, pooled_residuals_perm, n_residuals, n_studies, &
                                           study_block_ends, max_n_reps_all_studies, n_neighbors, bin_count_mode, &
                                           selected_n_bins, bin_count_failed, shared_residual_range_low, &
                                           shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                           mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, tmp_bin_counts, &
                                           tmp_residual_to_study_perm, tmp_study_block_sizes, m_min, m_max, &
                                           min_residuals_per_bin, gamma_occupancy, lower_residual_range_quantile, &
                                           upper_residual_range_quantile)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals, NaN included
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! One reference point's pooled residuals, one contiguous block per study
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: study_block_ends(n_studies)
            !! Last pool index of each study's block, ascending, the last one `n_residuals`
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies, for the Sturges/FD diagnostics
        integer(int32), intent(in) :: n_neighbors
            !! Neighbor count for the Sturges/FD diagnostics
        integer(int32), intent(in) :: bin_count_mode
            !! One of the four `MODE_BIN_COUNT_*` values
        integer(int32), intent(out) :: selected_n_bins
            !! The chosen bin count; m_min when bin_count_failed
        logical(c_bool), intent(out) :: bin_count_failed
            !! `.true.` iff the bin-count search failed (occupancy or Stone, per `bin_count_mode`)
        real(real64), intent(out) :: shared_residual_range_low
            !! Lower residual-range bound (R_low)
        real(real64), intent(out) :: shared_residual_range_high
            !! Upper residual-range bound (R_high)
        integer(int32), intent(out) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count of the pooled histogram at selected_n_bins
        real(real64), intent(out) :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count of the pooled histogram at selected_n_bins
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule diagnostic
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule diagnostic
        integer(int32), intent(out) :: tmp_bin_counts(MAX_N_BINS)
            !! Working array: bin counts at the bin count under test
        integer(int32), intent(out) :: tmp_residual_to_study_perm(n_residuals)
            !! Working array (Stone modes): the pool indices grouped by study, each block sorted
        integer(int32), intent(out) :: tmp_study_block_sizes(n_studies)
            !! Working array (Stone modes): the write cursor of each block, then its size
        integer(int32), intent(in), optional :: m_min
            !! Forwarded to the bin-count search
        integer(int32), intent(in), optional :: m_max
            !! Forwarded to the bin-count search
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Forwarded to the bin-count search; the floor in every mode
        real(real64), intent(in), optional :: gamma_occupancy
            !! Forwarded to the occupancy search; unused in the Stone modes
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Forwarded to the bin-count search
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Forwarded to the bin-count search

        integer(int32) :: i_study

        if (bin_count_mode == MODE_BIN_COUNT_OCCUPANCY) then
            call determine_bin_count_occupancy_impl(pooled_residuals, pooled_residuals_perm, n_residuals, &
                                                    max_n_reps_all_studies, n_neighbors, selected_n_bins, &
                                                    bin_count_failed, shared_residual_range_low, &
                                                    shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                                    mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                                    tmp_bin_counts, m_min, m_max, min_residuals_per_bin, gamma_occupancy, &
                                                    lower_residual_range_quantile, upper_residual_range_quantile)
            return
        end if

        ! The write cursor of each block starts just before it.
        tmp_study_block_sizes(1) = 0_int32
        do concurrent(i_study=2:n_studies) shared(tmp_study_block_sizes, study_block_ends)
            tmp_study_block_sizes(i_study) = study_block_ends(i_study - 1_int32)
        end do
        call group_sorted_perm_by_study(pooled_residuals_perm, n_residuals, n_studies, tmp_study_block_sizes, &
                                        tmp_residual_to_study_perm, study_block_ends=study_block_ends)
        ! The cursors now stand on the block ends; the core wants the block sizes.
        tmp_study_block_sizes(1) = study_block_ends(1)
        do concurrent(i_study=2:n_studies) shared(tmp_study_block_sizes, study_block_ends)
            tmp_study_block_sizes(i_study) = study_block_ends(i_study) - study_block_ends(i_study - 1_int32)
        end do

        call stone_bin_count_core(pooled_residuals, pooled_residuals_perm, n_residuals, n_studies, &
                                  tmp_residual_to_study_perm, tmp_study_block_sizes, max_n_reps_all_studies, n_neighbors, &
                                  bin_count_mode, selected_n_bins, bin_count_failed, shared_residual_range_low, &
                                  shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                  max_bin_occupancy, sturges_bins, fd_bins, tmp_bin_counts, m_min, m_max, &
                                  min_residuals_per_bin, lower_residual_range_quantile, upper_residual_range_quantile)
    end subroutine select_point_bin_count

    !> summary: Pool one reference point's residuals across all studies and pick its histogram bin count
    !| AUTHOR_LASZLO_LANG
    !| For ONE reference point whose studies may have different neighbor counts
    !| `point_n_neighbors`: pools the point's residuals study after study into
    !| `tmp_pooled_residuals` (the layout of `gather_pooled_neighborhood_residuals` when all counts
    !| are equal), sorts them (NaN last) and runs the bin-count search `bin_count_mode` names
    !| (`select_point_bin_count`). The Sturges/Freedman-Diaconis diagnostics get the rounded mean
    !| per-study neighbor count `(sum(point_n_neighbors) + n_studies/2)/n_studies`. Internal helper,
    !| not published.
    !|
    !| Precondition, not checked here: every `point_n_neighbors(i_study)` is at least 1 and at most
    !| `max_n_neighbors`, and `bin_count_mode` is one of the four `MODE_BIN_COUNT_*` values. A gene
    !| index outside `[1, max_n_genes_all_studies]` is reported by
    !| `gather_pooled_neighborhood_residuals` itself; its error code (without its argument position,
    !| which numbers that routine's own arguments) is returned in `ierr` and this helper returns at
    !| once, leaving every other output undefined.
    pure subroutine determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                              max_n_neighbors, n_index_points, neighborhood_indices, i_index_point, &
                                              point_n_neighbors, selected_n_bins, occupancy_failed, shared_residual_range_low, &
                                              shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, &
                                              mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                              tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts, &
                                              tmp_residual_to_study_perm, tmp_study_block_ends, tmp_study_block_sizes, &
                                              ierr, m_min, m_max, min_residuals_per_bin, gamma_occupancy, &
                                              lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_neighbors
            !! Leading extent of `neighborhood_indices`: the largest neighbor count of any study
        integer(int32), intent(in) :: n_index_points
            !! Second extent of `neighborhood_indices`: the number of reference points it holds
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices(max_n_neighbors, n_index_points, n_studies)
            !! Gene indices of the neighborhoods, per reference point, per study; only the first
            !! `point_n_neighbors(i_study)` entries of point `i_index_point` are read
        integer(int32), intent(in) :: i_index_point
            !! Which reference point of `neighborhood_indices` this call handles, in `[1, n_index_points]`
        integer(int32), intent(in) :: point_n_neighbors(n_studies)
            !! Number of neighbors of this reference point in each study, in `[1, max_n_neighbors]`
        integer(int32), intent(out) :: selected_n_bins
            !! This reference point's bin count, from the bin-count search
        logical(c_bool), intent(out) :: occupancy_failed
            !! `.true.` iff the bin-count search failed (occupancy or Stone, per `bin_count_mode`);
            !! `selected_n_bins` is then `m_min`
        real(real64), intent(out) :: shared_residual_range_low
            !! This reference point's lower residual-range bound (R_low), from the bin-count search
        real(real64), intent(out) :: shared_residual_range_high
            !! This reference point's upper residual-range bound (R_high), from the bin-count search
        integer(int32), intent(out) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from the bin-count search
        integer(int32), intent(out) :: min_bin_occupancy
            !! Minimum bin count at `selected_n_bins`, from the bin-count search
        real(real64), intent(out) :: mean_bin_occupancy
            !! Mean bin count at `selected_n_bins`, from the bin-count search
        integer(int32), intent(out) :: max_bin_occupancy
            !! Maximum bin count at `selected_n_bins`, from the bin-count search
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule diagnostic, computed with the rounded mean per-study neighbor count
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule diagnostic, with the rounded mean per-study neighbor count
        real(real64), intent(out) :: tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors*n_studies)
            !! Working array: this reference point's pooled residuals (leading part used)
        integer(int32), intent(out) :: tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors*n_studies)
            !! Working array: sorting permutation for the used part of `tmp_pooled_residuals`
        integer(int32), intent(out) :: tmp_bin_counts(MAX_N_BINS)
            !! Working array forwarded to the bin-count search's `tmp_bin_counts`
        integer(int32), intent(out) :: tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors*n_studies)
            !! Working array (Stone modes): pool indices grouped by study, each study sorted
        integer(int32), intent(out) :: tmp_study_block_ends(n_studies)
            !! Working array: the last pool index of each study's block
        integer(int32), intent(out) :: tmp_study_block_sizes(n_studies)
            !! Working array (Stone modes): the number of pooled residuals of each study
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT if a gene index is out of range (see above)
        integer(int32), intent(in), optional :: m_min
            !! Forwarded to the bin-count search
        integer(int32), intent(in), optional :: m_max
            !! Forwarded to the bin-count search
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Forwarded to the bin-count search, the occupancy floor in every mode
        real(real64), intent(in), optional :: gamma_occupancy
            !! Forwarded to the occupancy search; unused in the Stone modes
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Forwarded to the bin-count search
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Forwarded to the bin-count search
        integer(int32), intent(in), optional :: bin_count_mode
            !! Which bin-count search runs, one of the four `MODE_BIN_COUNT_*` values; default
            !! CM_MODE_BIN_COUNT_OCCUPANCY

        integer(int32) :: i_study, n_study_neighbors, pool_offset, gather_ierr, diagnostic_n_neighbors
        integer(int32) :: actual_bin_count_mode

        call set_ok(ierr)
        M_DEFAULT_VAL(bin_count_mode, actual_bin_count_mode, CM_MODE_BIN_COUNT_OCCUPANCY)

        ! Sequential: each study's block starts where the previous one ended. The index column is
        ! passed by its first element, so nothing is copied.
        pool_offset = 0_int32
        do i_study = 1, n_studies
            n_study_neighbors = point_n_neighbors(i_study)
            call gather_pooled_neighborhood_residuals(residuals(:, :, i_study:i_study), max_n_reps_all_studies, &
                                                      max_n_genes_all_studies, n_study_neighbors, 1_int32, &
                                                      neighborhood_indices(1, i_index_point, i_study), &
                                                      tmp_pooled_residuals(pool_offset + 1: &
                                                                           pool_offset + max_n_reps_all_studies*n_study_neighbors), &
                                                      gather_ierr)
            if (is_err(gather_ierr)) then
                call set_err_once(ierr, get_err_code(gather_ierr))
                return
            end if
            pool_offset = pool_offset + max_n_reps_all_studies*n_study_neighbors
            tmp_study_block_ends(i_study) = pool_offset
        end do

        call init_perm(tmp_pooled_residuals_perm(1:pool_offset))
        call sort_array_heapsort(tmp_pooled_residuals(1:pool_offset), tmp_pooled_residuals_perm(1:pool_offset))

        diagnostic_n_neighbors = (sum(point_n_neighbors) + n_studies/2_int32)/n_studies

        call select_point_bin_count(tmp_pooled_residuals(1:pool_offset), tmp_pooled_residuals_perm(1:pool_offset), &
                                    pool_offset, n_studies, tmp_study_block_ends, max_n_reps_all_studies, &
                                    diagnostic_n_neighbors, actual_bin_count_mode, selected_n_bins, occupancy_failed, &
                                    shared_residual_range_low, shared_residual_range_high, n_pooled_residuals, &
                                    min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, sturges_bins, fd_bins, &
                                    tmp_bin_counts, tmp_residual_to_study_perm(1:pool_offset), tmp_study_block_sizes, &
                                    m_min, m_max, min_residuals_per_bin, gamma_occupancy, &
                                    lower_residual_range_quantile, upper_residual_range_quantile)
    end subroutine determine_point_bin_count

    !> summary: Build one study's residual histogram for one reference point
    !| AUTHOR_LASZLO_LANG
    !| For ONE (reference point, study) pair: gathers the neighbors' residuals and bins them with
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] for a
    !| single point, which gives exactly the matching row of a call over all points. Only row
    !| `i_point` of `counts`/`pmf` is written, zero beyond `n_bins`; NaN residuals are skipped.
    !| Internal helper, not published.
    !|
    !| Precondition, not checked here: every entry of `neighbor_indices` is in
    !| `[1, max_n_genes_all_studies]`, `1 <= n_bins <= max_n_bins <= MAX_N_BINS`, and
    !| `1 <= i_point <= n_points`.
    pure subroutine build_point_study_histogram(residuals_study, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                n_neighbors, neighbor_indices, shared_residual_range_low, &
                                                shared_residual_range_high, n_bins, n_points, max_n_bins, i_point, &
                                                counts, pmf, included_n_reps, tmp_neighbor_residuals)
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
        integer(int32), intent(in) :: n_points
            !! First extent of `counts`/`pmf`: the number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! Second extent of `counts`/`pmf`: the widest bin count of any reference point
        integer(int32), intent(in) :: i_point
            !! The row of `counts`/`pmf` this reference point's histogram is written to
        integer(int32), intent(inout) :: counts(n_points, max_n_bins)
            !! The study's histogram counts; row `i_point` is overwritten, zero beyond `n_bins`
        real(real64), intent(inout) :: pmf(n_points, max_n_bins)
            !! The study's pmfs; row `i_point` is overwritten (all zero when `included_n_reps` is 0)
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

        counts(i_point, :) = point_counts(1, 1:max_n_bins)
        pmf(i_point, :) = point_pmf(1, 1:max_n_bins)
        included_n_reps = point_included_n_reps(1)
    end subroutine build_point_study_histogram

    !> Pass C of the fixed-k pipeline for ONE study: gathers every reference point's neighbor
    !| residuals into `tmp_neighborhood_residuals_gathered`, then builds all of the study's
    !| histograms with one
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] call at
    !| the per-point bin counts. Private helper of `run_js_comp_test_impl` and the fixed-k
    !| parameter search. Precondition, not checked here: every index is in
    !| `[1, max_n_genes_all_studies]` and `1 <= n_bins_per_point <= max_n_bins`.
    pure subroutine build_fixed_k_study_histograms(residuals_study, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                   n_neighbors, n_points, neighborhood_indices_study, &
                                                   shared_residual_range_low, shared_residual_range_high, &
                                                   n_bins_per_point, max_n_bins, counts, pmfs, included_n_reps, &
                                                   tmp_neighborhood_residuals_gathered)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per reference point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        real(real64), intent(in) :: residuals_study(max_n_reps_all_studies, max_n_genes_all_studies)
            !! This study's signed residuals, NaN allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices_study(n_neighbors, n_points)
            !! This study's neighbor gene indices per reference point
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! Lower histogram range bound (R_low) per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! Upper histogram range bound (R_high) per reference point
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(in) :: max_n_bins
            !! Second extent of `counts`/`pmfs`: the widest bin count of any reference point
        integer(int32), intent(out) :: counts(n_points, max_n_bins)
            !! Residual count per bin, per reference point; zero beyond each point's bin count
        real(real64), intent(out) :: pmfs(n_points, max_n_bins)
            !! `counts` normalized per reference point
        integer(int32), intent(out) :: included_n_reps(n_points)
            !! Number of non-NaN residuals binned per reference point
        real(real64), intent(out) :: tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points)
            !! Working array: the gathered residuals of every reference point's neighbors

        integer(int32) :: i_point, i_neighbor, gene_idx

        do concurrent(i_point=1:n_points, i_neighbor=1:n_neighbors) local(gene_idx) &
                shared(tmp_neighborhood_residuals_gathered, residuals_study, neighborhood_indices_study)
            gene_idx = neighborhood_indices_study(i_neighbor, i_point)
            tmp_neighborhood_residuals_gathered(:, i_neighbor, i_point) = residuals_study(:, gene_idx)
        end do

        call build_residual_histograms_impl(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies, n_neighbors, &
                                            n_points, shared_residual_range_low, shared_residual_range_high, max_n_bins, &
                                            n_bins_per_point, counts, pmfs, included_n_reps)
    end subroutine build_fixed_k_study_histograms

    !> Pass C on adaptive neighborhoods: decodes every reference point's pooled range into
    !| per-study gene lists
    !| ([[tox_data_integration_js_comp_test_impl(module):materialize_pooled_neighborhood(subroutine)]])
    !| and builds each (point, study) histogram with `build_point_study_histogram` into bins
    !| `1:max_n_bins` of `counts`/`pmfs`. Sequential, so one gather buffer serves every pair.
    !| Private helper of `run_js_comp_test_adaptive_impl` (`n_bins_extent = CM_MAX_N_BINS`) and the
    !| adaptive parameter search (`n_bins_extent = max_n_bins`, its buffers viewed packed).
    pure subroutine build_adaptive_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                              n_points, max_n_bins, n_bins_extent, gene_means_perm_all, &
                                              pooled_neighborhood_range, n_neighbors_per_point, &
                                              shared_residual_range_low, shared_residual_range_high, n_bins_per_point, &
                                              counts, pmfs, included_n_reps, tmp_point_neighborhood_indices, &
                                              tmp_neighbor_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! Widest per-point bin count; only bins `1:max_n_bins` of `counts`/`pmfs` are written
        integer(int32), intent(in) :: n_bins_extent
            !! Second extent of `counts`/`pmfs`, at least `max_n_bins`
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: gene_means_perm_all(max_n_genes_all_studies*n_studies)
            !! Sorting permutation of the pooled gene means, as the construction sorts them
        integer(int32), intent(in) :: pooled_neighborhood_range(2, n_points)
            !! Each reference point's `[first, last]` pooled position
        integer(int32), intent(out) :: n_neighbors_per_point(n_studies, n_points)
            !! Each reference point's per-study gene count, as decoded
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! Lower histogram range bound (R_low) per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! Upper histogram range bound (R_high) per reference point
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(inout) :: counts(n_points, n_bins_extent, n_studies)
            !! Every study's histogram counts; bins `1:max_n_bins` are overwritten
        real(real64), intent(inout) :: pmfs(n_points, n_bins_extent, n_studies)
            !! Every study's pmfs; bins `1:max_n_bins` are overwritten
        integer(int32), intent(out) :: included_n_reps(n_points, n_studies)
            !! Number of non-NaN residuals binned per reference point, per study
        integer(int32), intent(out) :: tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies)
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), intent(out) :: tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies)
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair

        integer(int32) :: i_point, i_study, n_study_neighbors

        do i_point = 1, n_points
            call materialize_pooled_neighborhood(gene_means_perm_all, max_n_genes_all_studies, n_studies, &
                                                 pooled_neighborhood_range(1, i_point), &
                                                 pooled_neighborhood_range(2, i_point), max_n_genes_all_studies, &
                                                 tmp_point_neighborhood_indices, n_neighbors_per_point(:, i_point))
            do i_study = 1, n_studies
                n_study_neighbors = n_neighbors_per_point(i_study, i_point)
                call build_point_study_histogram(residuals(:, :, i_study), max_n_reps_all_studies, &
                                                 max_n_genes_all_studies, n_study_neighbors, &
                                                 tmp_point_neighborhood_indices(1:n_study_neighbors, i_study), &
                                                 shared_residual_range_low(i_point), shared_residual_range_high(i_point), &
                                                 n_bins_per_point(i_point), n_points, max_n_bins, i_point, &
                                                 counts(:, 1:max_n_bins, i_study), pmfs(:, 1:max_n_bins, i_study), &
                                                 included_n_reps(i_point, i_study), &
                                                 tmp_neighbor_residuals(:, 1:n_study_neighbors))
            end do
        end do
    end subroutine build_adaptive_histograms

    !> summary: Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
    !| AUTHOR_LASZLO_LANG
    !| Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
    !| `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
    !| `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
    !| `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
    !| Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
    !| candidate. Bin counts are not part of a candidate; they are chosen per reference point later.
    pure subroutine generate_js_comp_test_candidates_impl(max_n_genes_all_studies, candidates_n_points_n_neighbors, &
                                                          n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: candidates_n_points_n_neighbors(2, CM_MAX_CANDIDATE_PAIRS)
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate pairs actually filled (at most CM_MAX_CANDIDATE_PAIRS)

        real(real64) :: n_points_high, n_points_low
        integer(int32) :: i_point_candidate, point_candidate, prev_point_candidate
        integer(int32) :: i_neighbor_candidate, neighbor_candidate, prev_neighbor_candidate

        n_points_high = real(clamp(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), &
                                   min_val=MIN_POINTS, max_val=MAX_POINTS), kind=real64)
        n_points_low = real(max(MIN_POINTS, ceiling(0.2_real64*n_points_high, kind=int32)), kind=real64)

        prev_point_candidate = -1_int32
        n_candidates = 0_int32

        ! Sequential: each iteration dedups against the previous one.
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
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]].
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
    !| - stop once `k_start_t < floor_k`, or once CM_MAX_CANDIDATE_PAIRS candidates are accepted
    !| - `k_step_t = ceiling(0.25 * k_start_t)`
    !| - `k_max_t = ceiling(4 * k_start_t)`
    !|
    !| The constants are [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_START_FRACTION(variable)]] (0.02),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_GAMMA(variable)]] (0.8),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_K_STEP_FRACTION(variable)]] (0.25),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]] (4),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_K_START_MIN_ABS(variable)]] (10) and
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_CANDIDATE_PAIRS(variable)]] (CM_MAX_CANDIDATE_PAIRS).
    !| The values strictly descend without deduplication. At least one candidate is always
    !| returned, even for data smaller than `k_start_1`, whose construction then fails with a
    !| too-few-means status in the search.
    !|
    !| Not checked here: the values fit 32-bit integers only while `k_start_1 <= huge(1_int32)/4`;
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]
    !| rejects inputs beyond that.
    pure subroutine generate_adaptive_js_comp_test_candidates_impl(max_n_genes_all_studies, n_studies, &
                                                                   candidates_k_start_k_step_k_max, n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: candidates_k_start_k_step_k_max(3, CM_MAX_CANDIDATE_PAIRS)
            !! Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate triples actually filled (at least 1, at most CM_MAX_CANDIDATE_PAIRS)

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
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
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

    !> summary: Grow Issue #217's adaptive neighborhoods over an already-sorted pool of gene means
    !| AUTHOR_LASZLO_LANG
    !| The construction of
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]],
    !| documented there, on an already sorted permutation (ascending, NaN last), so the adaptive
    !| parameter search can grow every candidate on one shared sort.
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
            ! to the lower mean; every branch takes the first position of a run of equal means.
            ! Seeds strictly advance, so there are at most n_valid_means points.
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
    !| Instead of one fixed neighbor count, each neighborhood grows while the mean absolute residual
    !| of its pool stays stable, and the next reference point is placed beyond it; the number of
    !| reference points emerges. `k_start`, `k_step` and `k_max` count pooled entries (one entry is
    !| one gene of one study).
    !|
    !| **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
    !| n_studies` entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted ascending; NaN
    !| means take no part. A neighborhood is a contiguous run `a..b` of that order
    !| (`pooled_neighborhood_range`, the input of
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]),
    !| and its per-study counts `n_neighbors_per_point` add up to `b - a + 1`.
    !|
    !| **Growth** around the reference point `x_star` (the mean at the seed position): entries join
    !| nearest first, ties to the lower mean; the dispersion `S` is the mean absolute non-NaN
    !| residual of the pool.
    !|
    !| 1. The first `k_start` entries join unconditionally. With fewer than
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !|    (`CM_ADAPTIVE_MIN_VALID_RESIDUALS`) non-NaN residuals the neighborhood stops with
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !|    (`CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS`); with all residuals zero, with
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !|    (`CM_ADAPTIVE_STOP_ZERO_DISPERSION`).
    !| 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    !|    `(S_new - S_old) / S_old > tau` is discarded and growth stops
    !|    ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TAU(variable)]], `CM_ADAPTIVE_STOP_TAU`);
    !|    otherwise it is committed. Growth also stops at `k_max` entries
    !|    ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `CM_ADAPTIVE_STOP_K_MAX`,
    !|    checked first) or when no entry is left
    !|    ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `CM_ADAPTIVE_STOP_EXHAUSTED`).
    !| 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    !|    non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    !|    is still evaluated, and growth then stops with
    !|    [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`CM_ADAPTIVE_STOP_RESIDUAL_CAP`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    !| sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    !| `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    !| absolute deviation of the neighborhood's pooled means. The construction ends with the
    !| neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`CM_ADAPTIVE_STATUS_TOO_FEW_MEANS`) when fewer than `k_start` pooled means are non-NaN: nothing
    !| is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    !| and reports, in this priority,
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS`) if any neighborhood stopped with too few residuals,
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD`) if any has fewer than `min_study_neighbors`
    !| entries of some study, else
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`CM_ADAPTIVE_STATUS_OK`).
    !|
    !| Outputs sized by `N` are filled only for their first `n_points` reference points.
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
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
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
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
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
    !| For one neighborhood built by
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]:
    !| sorted pooled position `pos` names the flat entry `flat = perm_all(pos)`, gene
    !| `mod(flat - 1, G) + 1` of study `(flat - 1)/G + 1`. Each study's genes are appended to its
    !| column of `point_neighborhood_indices` in pooled order and counted in `point_n_neighbors`.
    !| Internal helper, not published.
    !|
    !| Precondition, not checked here: `1 <= first_pos <= last_pos <= max_n_genes_all_studies *
    !| n_studies`, and `max_n_neighbors` is at least the largest per-study count (`G` always is).
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
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
    pure subroutine check_neighborhood_overlaps_impl(neighborhood_range, n_points, min_neighbor_overlap, &
                                                      all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: neighborhood_range(2, n_points)
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
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
            overlap = compute_index_span_overlap(neighborhood_range(1, i_point), neighborhood_range(2, i_point), &
                                                 neighborhood_range(1, i_point + 1), neighborhood_range(2, i_point + 1))
            all_pass = all_pass .and. (overlap >= min_neighbor_overlap)
        end do
        all_have_min_neighbor_overlap = logical(all_pass, kind=c_bool)
    end subroutine check_neighborhood_overlaps_impl

    !> summary: Test whether every bin of a mean pmf reaches a minimum absolute count
    !| AUTHOR_LASZLO_LANG
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]].
    !| Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    !| them does not count as a failure.
    pure subroutine check_mean_pmf_min_counts_impl(mean_pmf_counts, n_bins, n_bins_per_point, n_points, min_count, &
                                                   all_bins_have_min_count)
        integer(int32), intent(in) :: n_bins
            !! Second extent of `mean_pmf_counts`: the widest bin count of any point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
            !! is inspected
            !! DM_MIN(1_int32)
            !! DM_MAX(n_bins)
        integer(int32), intent(in) :: mean_pmf_counts(n_points, n_bins)
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
        do concurrent(i_bin=1:n_bins, i_point=1:n_points) shared(mean_pmf_counts, min_count, n_bins_per_point) &
                reduce(.and.:all_pass)
            if (i_bin <= n_bins_per_point(i_point)) then
                all_pass = all_pass .and. (mean_pmf_counts(i_point, i_bin) >= min_count)
            end if
        end do
        all_bins_have_min_count = logical(all_pass, kind=c_bool)
    end subroutine check_mean_pmf_min_counts_impl

    !> Calculates the fractional overlap between two intervals: the fraction of the overlap over
    !| the total range of the first interval, so for `a_min < b_min`,
    !| \[ \frac{\min(\texttt{a\_max}, \texttt{b\_max}) - \texttt{b\_min}}{\texttt{a\_max} - \texttt{a\_min}} \]
    !| For continuous intervals; inclusive integer index spans use `compute_index_span_overlap`.
    !| Private helper of
    !| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]].
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

    !> Counts the positions two inclusive integer runs `[a_min, a_max]` and `[b_min, b_max]` share,
    !| relative to the number of positions in the first run:
    !| \[ \frac{\max(0, \min(\texttt{a\_max}, \texttt{b\_max}) - \max(\texttt{a\_min}, \texttt{b\_min}) + 1)}{\texttt{a\_max} - \texttt{a\_min} + 1} \]
    !| For example `[5,6]` and `[6,7]` give 0.5. A reversed run (`max < min`) is empty, so the
    !| result is 0 rather than a division by zero. Private helper of
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]].
    pure real(real64) function compute_index_span_overlap(a_min, a_max, b_min, b_max) result(overlap_fraction)
        integer(int32), intent(in) :: a_min
            !! First position of the first run
        integer(int32), intent(in) :: a_max
            !! Last position of the first run; the run is empty if `a_max < a_min`
        integer(int32), intent(in) :: b_min
            !! First position of the second run
        integer(int32), intent(in) :: b_max
            !! Last position of the second run; the run is empty if `b_max < b_min`

        integer(int32) :: n_shared

        ! An empty first run shares no position
        if (a_max < a_min) then
            overlap_fraction = 0.0_real64
            return
        end if

        n_shared = max(0_int32, min(a_max, b_max) - max(a_min, b_min) + 1_int32)
        overlap_fraction = real(n_shared, real64)/real(a_max - a_min + 1_int32, real64)
    end function compute_index_span_overlap

    !> summary: Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau
    !| AUTHOR_LASZLO_LANG
    !| For each study, the overlap of the new candidate's interval with the best one so far is
    !| measured as a fraction of the new interval's width, and counted when it reaches
    !| `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    !| best (the best is then kept), or, after the new candidate became the best, when the count
    !| meets the condition `join_method` names.
    pure subroutine check_plateau_condition_impl(confidence_interval, best_candidate_pair_confidence_interval, n_studies, &
                                                 best_candidate_index, best_exceeded_ci_overlap_count, candidate_index, &
                                                 join_method, succeeding_ci_overlap, ci_overlap, plateau_found)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: confidence_interval(2, n_studies)
            !! Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
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
            !! How the studies' overlaps are joined into the plateau condition; the median method
            !! requires `count > (n_studies - 1) / 2`
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum fractional overlap with the best candidate's interval for a study to pass
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(out) :: ci_overlap(n_studies)
            !! Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            !! it was on entry, relative to the new interval's width; always filled. Meaningless
            !! for the very first candidate, whose best interval is only a seed
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate

        integer(int32) :: exceeds_min_ci_overlap, i_study
        real(real64) :: overlap
        logical(c_bool) :: plateau

        exceeds_min_ci_overlap = 0_int32
        do concurrent(i_study=1:n_studies) local(overlap) shared(confidence_interval, &
                                                                  best_candidate_pair_confidence_interval, &
                                                                  succeeding_ci_overlap, ci_overlap) &
                                                           reduce(+:exceeds_min_ci_overlap)
            ! The new candidate goes first: its own width is the denominator, so an interval inside
            ! the best one scores 1.0 and a wider one less.
            overlap = compute_fractional_overlap(confidence_interval(1, i_study), confidence_interval(2, i_study), &
                                                 best_candidate_pair_confidence_interval(1, i_study), &
                                                 best_candidate_pair_confidence_interval(2, i_study))
            ci_overlap(i_study) = overlap
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
    !| The relative-effect-size alternative to
    !| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]:
    !| for each study `i`, the relative change of the observed global JSD between successive
    !| admissible candidates,
    !| `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    !| delta_epsilon)`, is summarized by its median (`delta_median`) and maximum (`delta_max`). A
    !| plateau is found once both stay below their thresholds for
    !| `delta_min_consecutive_transitions` consecutive transitions, counted across calls in
    !| `n_consecutive_ok`.
    !|
    !| For the first admissible candidate (`has_previous = .false.`) there is no transition:
    !| `delta`, `delta_median` and `delta_max` are `-1.0` ("not computed").
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
            !! `.false.` for the first admissible candidate, which has no transition
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! DM_MIN(above(0.0_real64))
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! DM_MIN(1_int32)
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
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
        call calc_quantile_impl(delta, n_studies, tmp_delta_perm, 0.5_real64, delta_median)
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
    !| `mean_pmf` is the equal-weight mean of all studies' pmfs, `mean_pmf_counts` and
    !| `mean_pmf_included_n_reps` the sums of their counts. Known limitation: every study is
    !| compared against a consensus that includes itself, not a leave-one-out background.
    pure subroutine create_mean_pmf_impl(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                                         mean_pmf_included_n_reps, mean_pmf_counts)
        integer(int32), intent(in) :: n_bins
            !! Second extent of the bin-sized arrays: the widest bin count of any reference point.
            !! All studies must share each point's bin count, so their zero padding coincides
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_points, n_bins, n_studies), intent(in) :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), dimension(n_points, n_bins, n_studies), intent(in) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study; every study needs at
            !! least one at every reference point
            !! DM_MIN(1_int32)
        real(real64), dimension(n_points, n_bins), intent(out) :: mean_pmf
            !! The consensus pmf, `sum(pmfs, dim=3) / n_studies`
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point, summed over the studies
        integer(int32), dimension(n_points, n_bins), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus, `sum(counts, dim=3)`

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64
        mean_pmf_included_n_reps = 0_int32
        mean_pmf_counts = 0_int32

        ! Sequential over i_study, outermost, so every element accumulates its studies in order.
        ! For some study counts (e.g. 9, 11, 18, 20) the sum exceeds 1.0 by a few ulps when every
        ! study has all its mass in one bin, which gjct_permutation_test's DM_MAX(1.0) on mean_pmf
        ! would reject. Deliberately not clamped, so outputs stay bit-identical.
        do i_study = 1, n_studies
            do i_bin = 1, n_bins
                do concurrent(i_point=1:n_points) shared(mean_pmf, pmfs, i_bin, i_study, n_studies, mean_pmf_counts, counts)
                    mean_pmf(i_point, i_bin) = mean_pmf(i_point, i_bin) + pmfs(i_point, i_bin, i_study)/real(n_studies, real64)
                    mean_pmf_counts(i_point, i_bin) = mean_pmf_counts(i_point, i_bin) + counts(i_point, i_bin, i_study)
                end do
            end do
            do concurrent(i_point=1:n_points) shared(mean_pmf_included_n_reps, included_n_reps, i_study)
                mean_pmf_included_n_reps(i_point) = mean_pmf_included_n_reps(i_point) + included_n_reps(i_point, i_study)
            end do
        end do
    end subroutine create_mean_pmf_impl

    !> summary: Build only the consensus pmf from all studies' pmfs, without its histogram counts
    !| AUTHOR_LASZLO_LANG
    !| As [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]], for
    !| callers that need no counts (the bootstrap). Same known limitation: the consensus includes
    !| the study compared against it.
    pure subroutine create_mean_pmf_only_impl(pmfs, n_bins, n_points, n_studies, mean_pmf)
        integer(int32), intent(in) :: n_bins
            !! Second extent of the bin-sized arrays: the widest bin count of any reference point.
            !! All studies must share each point's bin count, so their zero padding coincides
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_points, n_bins, n_studies), intent(in) :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_points, n_bins), intent(out) :: mean_pmf
            !! The consensus pmf, `sum(pmfs, dim=3) / n_studies`

        integer(int32) :: i_study, i_point, i_bin

        mean_pmf = 0.0_real64

        ! Same order and the same few-ulp excess over 1.0 as create_mean_pmf_impl (not clamped).
        do i_study = 1, n_studies
            do i_bin = 1, n_bins
                do concurrent(i_point=1:n_points) shared(mean_pmf, pmfs, i_bin, i_study, n_studies)
                    mean_pmf(i_point, i_bin) = mean_pmf(i_point, i_bin) + pmfs(i_point, i_bin, i_study)/real(n_studies, real64)
                end do
            end do
        end do
    end subroutine create_mean_pmf_only_impl

    !> M_EXPORT_C
    !| summary: Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
    !| AUTHOR_LASZLO_LANG
    !| Sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    !| the default 2.5 and `n_bootstraps=1000`.
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
    !| In every bootstrap, each study's histogram at each reference point is drawn with replacement
    !| (multinomial) from the pooled consensus counts `mean_pmf_counts`, with as many residuals as
    !| the study observed; the consensus pmf, each study's JSD against it and its weighted global
    !| JSD are recomputed from the draw.
    !|
    !| `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
    !| slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
    !| values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
    !| smallest of the n_bootstrapping_top_k_jsds largest values]`.
    !|
    !| Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
    !| draw, when at some reference point the studies' counts do not add up to
    !| `mean_pmf_included_n_reps` or that total is not the sum of `mean_pmf_counts` (as
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] produces
    !| them); `ERR_ALLOC_FAIL` when the random number generator cannot be created. A failed draw
    !| later on is recorded (first failure only) without stopping the bootstrap.
    subroutine bootstrap_histogram_impl(n_bootstraps, n_bins, n_points, n_studies, mean_pmf_counts, &
                                        mean_pmf_included_n_reps, included_n_reps, n_bootstrapping_top_k_jsds, &
                                        confidence_interval, tmp_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, &
                                        tmp_mean_pmf, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                        tmp_point_pool, tmp_point_draw, two_sided_bootstrapping_significance_level, &
                                        random_seed, ierr)
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
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_js_comp_test_impl, AUTO)
            !! DM_MIN(1_int32)
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array: per study the bottom-k heap `(:, 1, :)` and the top-k heap `(:, 2, :)`
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array: one study's resampled histogram counts
        real(real64), dimension(n_points, n_bins, n_studies), intent(out) :: tmp_pmfs
            !! Working array: one bootstrap's resampled pmfs, all studies
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_mean_pmf
            !! Working array: one bootstrap's consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: one bootstrap's per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Working array: one bootstrap's per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: one bootstrap's weighted global JSD values
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_pool
            !! Working array: one reference point's pooled counts, the pool of one draw
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        integer(int32) :: i_bootstrap, i_study, i_point, draw_ierr
        type(rng_t) :: rng

        call set_ok(ierr)
        ! The pool's total is the multinomial total and the weights' denominator, so it must match
        ! the studies' draws and its own bin counts (also covers drawing from an empty pool). Not
        ! one argument's fault, so no position.
        do i_point = 1, n_points
            if (sum(included_n_reps(i_point, :)) /= mean_pmf_included_n_reps(i_point) .or. &
                mean_pmf_included_n_reps(i_point) /= sum(mean_pmf_counts(i_point, :))) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                exit
            end if
        end do
        if (.not. is_err(ierr)) rng = create_rng(ierr, random_seed)
        if (is_err(ierr)) then
            ! Leave the work arrays defined (zeroed) and confidence_interval untouched.
            tmp_bootstrapping_top_k_jsds = 0.0_real64
            tmp_counts = 0_int32
            tmp_pmfs = 0.0_real64
            tmp_mean_pmf = 0.0_real64
            tmp_js_divergences = 0.0_real64
            tmp_weights = 0.0_real64
            tmp_global_js_divergence = 0.0_real64
            tmp_point_pool = 0_int32
            tmp_point_draw = 0_int32
            return
        end if

        do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, confidence_interval)
            tmp_bootstrapping_top_k_jsds(:, 1, i_study) = confidence_interval(1, i_study)
            tmp_bootstrapping_top_k_jsds(:, 2, i_study) = confidence_interval(2, i_study)
        end do

        do i_bootstrap = 1, n_bootstraps
            ! Sequential: every draw advances the rng state. A point's bins are a strided row here,
            ! so each draw goes through the contiguous tmp_point_pool/tmp_point_draw.
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    tmp_point_pool = mean_pmf_counts(i_point, :)
                    call random_multinomial(rng, n_bins, tmp_point_pool, mean_pmf_included_n_reps(i_point), &
                                            included_n_reps(i_point, i_study), tmp_point_draw, draw_ierr)
                    if (is_err(draw_ierr)) call set_err_once(ierr, get_err_code(draw_ierr))
                    tmp_counts(i_point, :) = tmp_point_draw
                end do

                ! A study without residuals at a point is rejected by the wrappers (DM_MIN(1));
                ! calc_pmf_impl gives such a point an all-zero pmf instead of dividing by 0.
                call calc_pmf_impl(tmp_counts, included_n_reps(:, i_study), n_points, n_bins, tmp_pmfs(:, :, i_study))
            end do

            ! Drawn with replacement, so the consensus changes too.
            call create_mean_pmf_only_impl(tmp_pmfs, n_bins, n_points, n_studies, tmp_mean_pmf)

            do concurrent(i_study=1:n_studies) shared(tmp_pmfs, tmp_mean_pmf, n_points, n_bins, tmp_js_divergences, &
                                                      included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, &
                                                      tmp_weights, tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds)
                call compute_divergence_per_reference_point_impl(tmp_pmfs(:, :, i_study), tmp_mean_pmf, &
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

        ! The heap roots: [largest of the small values, smallest of the large values].
        do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, confidence_interval)
            confidence_interval(1, i_study) = tmp_bootstrapping_top_k_jsds(1, 1, i_study)
            confidence_interval(2, i_study) = tmp_bootstrapping_top_k_jsds(1, 2, i_study)
        end do

        call destroy_rng(rng)
    end subroutine bootstrap_histogram_impl

    !> The shared tail of both final runs, after every study's histograms are built: the consensus
    !| pmf (create_mean_pmf_impl), each study's JSD against it and the weighted global JSD, the
    !| permutation test (gjct_permutation_test_impl), and the re-derivation of each study's
    !| pmf/JSD/weights/global JSD from its untouched `counts`. Private; impure because the
    !| permutation test is. If that test fails, its error code is returned in `ierr` at once.
    !|
    !| The consensus is built over the full `MAX_N_BINS` extent, so `pmfs`/`counts` must be zero
    !| beyond `max_n_bins_per_point`; `mean_pmf`/`mean_pmf_counts` are then zero there too. The
    !| `tmp_permutation_*` work arrays are handed to gjct_permutation_test_impl whole, which views
    !| their leading elements packed at `max_n_bins_per_point` bins (sequence association, no copy).
    subroutine finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, mean_pmf, &
                                   mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                   global_js_divergence, p_values, p_values_observed_consensus, &
                                   tmp_permutation_mean_pmf, tmp_permutation_point_pool, tmp_permutation_point_draw, &
                                   tmp_permutation_mean_pmf_counts, tmp_permutation_counts, tmp_permutation_pmfs, &
                                   tmp_permutation_js_divergences, tmp_permutation_weights, &
                                   tmp_permutation_global_js_divergence, n_permutations, random_seed, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: max_n_bins_per_point
            !! Widest bin count of any reference point; only bins `1:max_n_bins_per_point` of the
            !! bin-sized arrays below are meaningful
        real(real64), dimension(n_points, MAX_N_BINS, n_studies), intent(inout) :: pmfs
            !! Per-study pmfs, zero beyond `max_n_bins_per_point`; read for the consensus, then
            !! re-derived from `counts`
        integer(int32), dimension(n_points, MAX_N_BINS, n_studies), intent(in) :: counts
            !! Per-study histogram counts, zero beyond `max_n_bins_per_point`
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, MAX_N_BINS), intent(out) :: mean_pmf
            !! The consensus pmf
        integer(int32), dimension(n_points, MAX_N_BINS), intent(out) :: mean_pmf_counts
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
        real(real64), dimension(n_points, MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array forwarded to gjct_permutation_test_impl
        integer(int32), dimension(MAX_N_BINS), intent(out) :: tmp_permutation_point_pool
            !! Working array forwarded to gjct_permutation_test_impl
        integer(int32), dimension(MAX_N_BINS), intent(out) :: tmp_permutation_point_draw
            !! Working array forwarded to gjct_permutation_test_impl
        integer(int32), dimension(n_points, MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to gjct_permutation_test_impl
        integer(int32), dimension(n_points, MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to gjct_permutation_test_impl
        real(real64), dimension(n_points, MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
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
            !! Error code of gjct_permutation_test_impl if it fails (argument position cleared); the
            !! re-derived pmfs/JSDs/weights/global JSD are then not computed

        integer(int32) :: i_study, permutation_ierr

        call set_ok(ierr)

        ! Over the full extent: the studies' bins 1:max_n_bins_per_point are not one contiguous
        ! block, and the zero padding sums to zero.
        call create_mean_pmf_impl(pmfs, counts, MAX_N_BINS, n_points, n_studies, included_n_reps, mean_pmf, &
                                  mean_pmf_included_n_reps, mean_pmf_counts)

        do concurrent(i_study=1:n_studies) shared(pmfs, mean_pmf, n_points, max_n_bins_per_point, js_divergences, included_n_reps, &
                                                  mean_pmf_included_n_reps, global_js_divergence, weights)
            call compute_divergence_per_reference_point_impl(pmfs(:, 1:max_n_bins_per_point, i_study), &
                                                             mean_pmf(:, 1:max_n_bins_per_point), n_points, &
                                                             max_n_bins_per_point, js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
        end do

        call gjct_permutation_test_impl(n_permutations, max_n_bins_per_point, n_points, n_studies, &
                                        mean_pmf_counts(:, 1:max_n_bins_per_point), mean_pmf(:, 1:max_n_bins_per_point), &
                                        mean_pmf_included_n_reps, included_n_reps, global_js_divergence, p_values, &
                                        p_values_observed_consensus, tmp_permutation_mean_pmf_counts, &
                                        tmp_permutation_counts, tmp_permutation_pmfs, &
                                        tmp_permutation_js_divergences, tmp_permutation_weights, &
                                        tmp_permutation_global_js_divergence, tmp_permutation_mean_pmf, &
                                        tmp_permutation_point_pool, tmp_permutation_point_draw, random_seed, permutation_ierr)
        if (is_err(permutation_ierr)) then
            call set_err_once(ierr, get_err_code(permutation_ierr))
            return
        end if

        ! Re-derive each study's pmf/JSD/weights/global JSD from its untouched counts. mean_pmf and
        ! mean_pmf_counts need no re-derivation: the permutation test only touches its scratch copies.
        do i_study = 1, n_studies
            call calc_pmf_impl(counts(:, 1:max_n_bins_per_point, i_study), included_n_reps(:, i_study), n_points, &
                               max_n_bins_per_point, pmfs(:, 1:max_n_bins_per_point, i_study))
            call compute_divergence_per_reference_point_impl(pmfs(:, 1:max_n_bins_per_point, i_study), &
                                                             mean_pmf(:, 1:max_n_bins_per_point), n_points, &
                                                             max_n_bins_per_point, js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
        end do
    end subroutine finish_js_comp_test

    !> summary: Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors) parameter setting
    !| AUTHOR_LASZLO_LANG
    !| The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
    !|
    !| 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    !|    the caller, e.g. with
    !|    [[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) form its
    !|    neighborhoods
    !|    ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| 2. Per reference point, all studies' residuals are pooled and the bin-count search
    !|    `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    !|    [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !|    by default, or Stone's rule
    !|    [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]).
    !| 3. Each study's histogram and pmf is built per point
    !|    ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !|    ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), each
    !|    study's per-point JSD against it and its weighted global JSD, and the permutation test
    !|    ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]) give
    !|    `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
    !|    `n_permutations = 0` both p-values are NaN.
    !|
    !| There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
    !| still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
    !| caller that needs such neighborhoods rejected should use
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !|
    !| The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
    !| second extent [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]]
    !| (`CM_MAX_N_BINS`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
    !| zero. A Python/R caller slices `[:, :max_n_bins_per_point, ...]` (R:
    !| `[, seq_len(max_n_bins_per_point), ]`) itself.
    !|
    !| Errors (`ierr`): `ERR_INVALID_INPUT` when some study has no non-NaN residual at some
    !| reference point (for example a study whose gene means are all missing), since the
    !| equal-weight consensus needs every study; the consensus, JSD, weight and p-value outputs are
    !| then not filled. `ERR_ALLOC_FAIL` when the permutation test cannot create its random number
    !| generator; `pmfs`, `js_divergences`, `weights` and `global_js_divergence` then hold their
    !| values from before the permutation test. An out-of-range neighbor index is also reported
    !| as `ERR_INVALID_INPUT`, a safeguard valid input cannot reach.
    subroutine run_js_comp_test_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                                     gene_means, gene_means_perms, residuals, x_star, &
                                     neighborhood_indices, neighborhood_range, n_bins_per_point, &
                                     shared_residual_range_low, shared_residual_range_high, max_n_bins_per_point, &
                                     occupancy_failed, n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                     max_bin_occupancy, sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, &
                                     mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                     global_js_divergence, p_values, p_values_observed_consensus, &
                                     tmp_neighborhood_residuals_gathered, &
                                     tmp_permutation_mean_pmf, tmp_permutation_point_pool, &
                                     tmp_permutation_point_draw, tmp_pooled_residuals, &
                                     tmp_pooled_residuals_perm, tmp_bin_counts_search, tmp_residual_to_study_perm, &
                                     tmp_study_block_ends, tmp_study_block_sizes, tmp_point_n_neighbors, &
                                     tmp_permutation_mean_pmf_counts, tmp_permutation_counts, tmp_permutation_pmfs, &
                                     tmp_permutation_js_divergences, tmp_permutation_weights, &
                                     tmp_permutation_global_js_divergence, n_permutations, &
                                     random_seed, min_residuals_per_bin, m_min, m_max, gamma_occupancy, &
                                     lower_residual_range_quantile, upper_residual_range_quantile, bin_count_mode, ierr)
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
        integer(int32), dimension(n_neighbors, n_points, n_studies), intent(out) :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span
        integer(int32), dimension(n_points), intent(out) :: n_bins_per_point
            !! Each reference point's histogram bin count `M_j`, from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(int32), intent(out) :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts` (see above)
        logical(c_bool), dimension(n_points), intent(out) :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes (see above)
        integer(int32), dimension(n_points), intent(out) :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(int32), dimension(n_points), intent(out) :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(real64), dimension(n_points), intent(out) :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(int32), dimension(n_points), intent(out) :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(real64), dimension(n_points, CM_MAX_N_BINS, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, CM_MAX_N_BINS, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, CM_MAX_N_BINS), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, CM_MAX_N_BINS), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: &
            tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        real(real64), dimension(n_points, CM_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the permutation test: consensus pmf of the permuted studies
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_permutation_point_pool
            !! Working array of the permutation test: one point's pool
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_permutation_point_draw
            !! Working array of the permutation test: one study's draw at one point
        real(real64), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array: every study's neighbor count at a reference point (all `n_neighbors`)
        integer(int32), dimension(n_points, CM_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the permutation test: the remaining pool
        integer(int32), dimension(n_points, CM_MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array of the permutation test: one study's resampled counts
        real(real64), dimension(n_points, CM_MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the permutation test: resampled pmfs
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the permutation test: per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the permutation test: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array of the permutation test: global JSD values
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(1000_int32)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(CM_MODE_BIN_COUNT_OCCUPANCY)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)

        integer(int32) :: i_study, i_point, actual_n_permutations, point_ierr, finish_ierr

        call set_ok(ierr)
        M_DEFAULT_VAL(n_permutations, actual_n_permutations, 1000_int32)

        ! ===== PASS A (per study): neighborhood construction.
        do i_study = 1, n_studies
            call construct_neighborhoods_ranged_impl(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), &
                                                     gene_means_perms(:, i_study), n_neighbors, &
                                                     neighborhood_indices(:, :, i_study), neighborhood_range(:, :, i_study))
        end do

        ! ===== PASS B (per point): pool all studies' residuals and run the bin-count search.
        ! Sequential, so one tmp_pooled_residuals buffer serves every point (do concurrent would
        ! need one per point).
        tmp_point_n_neighbors = n_neighbors
        do i_point = 1, n_points
            call determine_point_bin_count(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                           n_neighbors, n_points, neighborhood_indices, i_point, tmp_point_n_neighbors, &
                                           n_bins_per_point(i_point), occupancy_failed(i_point), &
                                           shared_residual_range_low(i_point), shared_residual_range_high(i_point), &
                                           n_pooled_residuals(i_point), min_bin_occupancy(i_point), &
                                           mean_bin_occupancy(i_point), max_bin_occupancy(i_point), &
                                           sturges_bins(i_point), fd_bins(i_point), tmp_pooled_residuals, &
                                           tmp_pooled_residuals_perm, tmp_bin_counts_search, tmp_residual_to_study_perm, &
                                           tmp_study_block_ends, tmp_study_block_sizes, point_ierr, &
                                           m_min=m_min, m_max=m_max, min_residuals_per_bin=min_residuals_per_bin, &
                                           gamma_occupancy=gamma_occupancy, &
                                           lower_residual_range_quantile=lower_residual_range_quantile, &
                                           upper_residual_range_quantile=upper_residual_range_quantile, &
                                           bin_count_mode=bin_count_mode)
            ! Defensive: only an out-of-range gene index fails here, which the wrapper's bounds
            ! exclude. Kept for direct callers, since Pass C would read out of bounds; untested.
            if (is_err(point_ierr)) then
                call set_err_once(ierr, get_err_code(point_ierr))
                return
            end if
        end do

        ! ===== PASS C (per study): build every (point, study) histogram at the chosen bin counts.
        ! The bins beyond max_n_bins_per_point are zeroed, as finish_js_comp_test requires.
        max_n_bins_per_point = maxval(n_bins_per_point(1:n_points))
        do i_study = 1, n_studies
            call build_fixed_k_study_histograms(residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, &
                                                n_neighbors, n_points, neighborhood_indices(:, :, i_study), &
                                                shared_residual_range_low, shared_residual_range_high, n_bins_per_point, &
                                                max_n_bins_per_point, counts(:, 1:max_n_bins_per_point, i_study), &
                                                pmfs(:, 1:max_n_bins_per_point, i_study), included_n_reps(:, i_study), &
                                                tmp_neighborhood_residuals_gathered)
        end do
        counts(:, max_n_bins_per_point + 1:, :) = 0_int32
        pmfs(:, max_n_bins_per_point + 1:, :) = 0.0_real64

        ! A study with no included residual at a point would leave that point's equal-weight
        ! consensus summing to less than 1. It depends on gene_means and residuals together, so
        ! no arg_pos.
        if (any(included_n_reps(1:n_points, 1:n_studies) == 0_int32)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, mean_pmf, &
                                 mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                 global_js_divergence, p_values, p_values_observed_consensus, &
                                 tmp_permutation_mean_pmf, tmp_permutation_point_pool, tmp_permutation_point_draw, &
                                 tmp_permutation_mean_pmf_counts, tmp_permutation_counts, tmp_permutation_pmfs, &
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
    !| same pipeline and outputs, on neighborhoods whose size varies per reference point and per
    !| study. Call
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| first and pass its `pooled_neighborhood_range` (its first `n_points` columns) with the same
    !| `gene_means` and `residuals`. Each range is decoded back to per-study gene lists; their sizes
    !| are returned as `n_neighbors_per_point`. Neighborhoods that decode to the same gene sets as a
    !| fixed-k run give bit-identical results. There is no admissibility gate.
    !|
    !| Errors (`ierr`): `ERR_INVALID_INPUT` before any histogram work for the first range that is
    !| reversed, reaches past the non-NaN pooled means, or leaves a study without a gene (every
    !| output then undefined); otherwise as in `run_js_comp_test_impl` (a study without a non-NaN
    !| residual at some point, or `ERR_ALLOC_FAIL` from the permutation test).
    !|
    !| The pooling work arrays hold `max_n_reps_all_studies * max_n_genes_all_studies * n_studies`
    !| values each, the size of `residuals`, since one study may contribute all its genes to one
    !| neighborhood. The bin-sized arrays have the fixed second extent
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]]; only their first
    !| `max_n_bins_per_point` bins are meaningful (Python `[:, :max_n_bins_per_point, ...]`, R
    !| `[, seq_len(max_n_bins_per_point), ]`).
    subroutine run_js_comp_test_adaptive_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, &
                                              gene_means, residuals, pooled_neighborhood_range, n_neighbors_per_point, &
                                              n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                              max_n_bins_per_point, occupancy_failed, n_pooled_residuals, &
                                              min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, sturges_bins, &
                                              fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, &
                                              p_values, p_values_observed_consensus, tmp_gene_means_perm_all, &
                                              tmp_point_neighborhood_indices, &
                                              tmp_neighbor_residuals, tmp_permutation_mean_pmf, &
                                              tmp_permutation_point_pool, tmp_permutation_point_draw, tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                              tmp_residual_to_study_perm, tmp_study_block_ends, tmp_study_block_sizes, &
                                              tmp_permutation_mean_pmf_counts, tmp_permutation_counts, &
                                              tmp_permutation_pmfs, tmp_permutation_js_divergences, &
                                              tmp_permutation_weights, tmp_permutation_global_js_divergence, &
                                              n_permutations, random_seed, min_residuals_per_bin, m_min, m_max, &
                                              gamma_occupancy, lower_residual_range_quantile, &
                                              upper_residual_range_quantile, bin_count_mode, ierr)
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
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! DM_ALLOW_NAN
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies*n_studies)
        integer(int32), dimension(n_studies, n_points), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(int32), dimension(n_points), intent(out) :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(real64), dimension(n_points), intent(out) :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(int32), intent(out) :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out) :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes
        integer(int32), dimension(n_points), intent(out) :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(int32), dimension(n_points), intent(out) :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(real64), dimension(n_points), intent(out) :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(int32), dimension(n_points), intent(out) :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(int32), dimension(n_points), intent(out) :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(real64), dimension(n_points, CM_MAX_N_BINS, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, CM_MAX_N_BINS, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(n_points, CM_MAX_N_BINS), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, CM_MAX_N_BINS), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        real(real64), dimension(n_points, CM_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf
            !! Working array forwarded to the permutation test
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_permutation_point_pool
            !! Working array forwarded to the permutation test
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_permutation_point_draw
            !! Working array forwarded to the permutation test
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        integer(int32), dimension(n_points, CM_MAX_N_BINS), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to the permutation test
        integer(int32), dimension(n_points, CM_MAX_N_BINS), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, CM_MAX_N_BINS, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array forwarded to the permutation test
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(1000_int32)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(CM_MODE_BIN_COUNT_OCCUPANCY)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)

        integer(int32) :: i_point, n_valid_means, actual_n_permutations, point_ierr, finish_ierr

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
                                           max_n_genes_all_studies, 1_int32, tmp_point_neighborhood_indices, 1_int32, &
                                           n_neighbors_per_point(:, i_point), n_bins_per_point(i_point), &
                                           occupancy_failed(i_point), shared_residual_range_low(i_point), &
                                           shared_residual_range_high(i_point), n_pooled_residuals(i_point), &
                                           min_bin_occupancy(i_point), mean_bin_occupancy(i_point), &
                                           max_bin_occupancy(i_point), sturges_bins(i_point), fd_bins(i_point), &
                                           tmp_pooled_residuals, tmp_pooled_residuals_perm, tmp_bin_counts_search, &
                                           tmp_residual_to_study_perm, tmp_study_block_ends, tmp_study_block_sizes, &
                                           point_ierr, m_min=m_min, m_max=m_max, &
                                           min_residuals_per_bin=min_residuals_per_bin, gamma_occupancy=gamma_occupancy, &
                                           lower_residual_range_quantile=lower_residual_range_quantile, &
                                           upper_residual_range_quantile=upper_residual_range_quantile, &
                                           bin_count_mode=bin_count_mode)
            ! Defensive: only an out-of-range gene index fails here, and every decoded index is in
            ! range; unreachable and untested.
            if (is_err(point_ierr)) then
                call set_err_once(ierr, get_err_code(point_ierr))
                return
            end if
        end do

        ! ===== PASS C (per point, per study, sequential: one gather buffer serves every pair).
        ! The bins beyond max_n_bins_per_point are zeroed, as finish_js_comp_test requires.
        max_n_bins_per_point = maxval(n_bins_per_point)
        counts(:, max_n_bins_per_point + 1:, :) = 0_int32
        pmfs(:, max_n_bins_per_point + 1:, :) = 0.0_real64
        call build_adaptive_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, n_points, &
                                       max_n_bins_per_point, CM_MAX_N_BINS, tmp_gene_means_perm_all, &
                                       pooled_neighborhood_range, n_neighbors_per_point, shared_residual_range_low, &
                                       shared_residual_range_high, n_bins_per_point, counts, pmfs, included_n_reps, &
                                       tmp_point_neighborhood_indices, tmp_neighbor_residuals)

        ! The gene-count check cannot see a study whose genes at a point have only NaN residuals.
        ! Such a study would leave that point's equal-weight consensus summing to less than 1.
        ! It depends on gene_means and residuals together, so no arg_pos.
        if (any(included_n_reps(1:n_points, 1:n_studies) == 0_int32)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call finish_js_comp_test(n_studies, n_points, max_n_bins_per_point, pmfs, counts, included_n_reps, mean_pmf, &
                                 mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                 global_js_divergence, p_values, p_values_observed_consensus, &
                                 tmp_permutation_mean_pmf, tmp_permutation_point_pool, tmp_permutation_point_draw, &
                                 tmp_permutation_mean_pmf_counts, tmp_permutation_counts, tmp_permutation_pmfs, &
                                 tmp_permutation_js_divergences, tmp_permutation_weights, &
                                 tmp_permutation_global_js_divergence, actual_n_permutations, random_seed, &
                                 finish_ierr)
        if (is_err(finish_ierr)) then
            call set_err_once(ierr, get_err_code(finish_ierr))
            return
        end if
    end subroutine run_js_comp_test_adaptive_impl

    !> Pass A of the fixed-k parameter search for one candidate: builds every study's
    !| neighborhoods into `neighborhood_indices` and checks the first admissibility gate
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]])
    !| study by study, stopping at the first study that fails it. Private helper; its
    !| explicit-shape dummies take the search's upper-bound-sized buffers whole and view their
    !| leading elements packed at the candidate's own sizes (see the parameter search).
    pure subroutine construct_fixed_k_candidate_neighborhoods(n_studies, max_n_genes_all_studies, n_points, n_neighbors, &
                                                              x_star, gene_means, gene_means_perms, min_neighbor_overlap, &
                                                              neighborhood_indices, neighborhood_range, &
                                                              all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        integer(int32), intent(in) :: n_neighbors
            !! The candidate's number of neighbors per reference point
        real(real64), intent(in) :: x_star(n_points)
            !! The candidate's reference points
        real(real64), intent(in) :: gene_means(max_n_genes_all_studies, n_studies)
            !! Per-gene mean expression values for all studies
        integer(int32), intent(in) :: gene_means_perms(max_n_genes_all_studies, n_studies)
            !! Each study's sorting permutation for `gene_means`
        real(real64), intent(in) :: min_neighbor_overlap
            !! Forwarded to check_neighborhood_overlaps_impl
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points, n_studies)
            !! Every study's neighbor gene indices; only the studies up to the first failing one
            !! are filled
        integer(int32), intent(out) :: neighborhood_range(2, n_points)
            !! One study's neighborhood spans, reused per study
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` iff every study passed the first admissibility gate

        integer(int32) :: i_study

        all_have_min_neighbor_overlap = logical(.true., kind=c_bool)
        do i_study = 1, n_studies
            call construct_neighborhoods_ranged_impl(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), &
                                                     gene_means_perms(:, i_study), n_neighbors, &
                                                     neighborhood_indices(:, :, i_study), neighborhood_range)
            call check_neighborhood_overlaps_impl(neighborhood_range, n_points, min_neighbor_overlap, &
                                                   all_have_min_neighbor_overlap)
            if (.not. all_have_min_neighbor_overlap) exit
        end do
    end subroutine construct_fixed_k_candidate_neighborhoods

    !> Pools one reference point's residuals across all its neighbors of all studies, for Pass B
    !| of the fixed-k parameter search: one
    !| [[tox_data_integration_js_comp_test_impl(module):gather_pooled_neighborhood_residuals(interface)]]
    !| call per study into that study's block of `pooled_residuals`, which gives exactly the layout
    !| of a single multi-study call. Private helper; `neighborhood_indices` is the search's buffer
    !| viewed packed at the candidate's sizes. A gene index out of range returns that routine's
    !| error code in `ierr` (its argument position cleared) and stops at once.
    pure subroutine gather_fixed_k_point_pool(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                              n_neighbors, n_points, neighborhood_indices, i_point, pooled_residuals, ierr)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_neighbors
            !! The candidate's number of neighbors per reference point
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices(n_neighbors, n_points, n_studies)
            !! Every study's neighbor gene indices
        integer(int32), intent(in) :: i_point
            !! The reference point to pool
        real(real64), intent(out) :: pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies)
            !! The point's pooled residuals, one block of `max_n_reps_all_studies*n_neighbors` per study
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for a gene index out of range

        integer(int32) :: i_study, block_size, gather_ierr

        call set_ok(ierr)
        block_size = max_n_reps_all_studies*n_neighbors
        do i_study = 1, n_studies
            call gather_pooled_neighborhood_residuals(residuals(:, :, i_study:i_study), max_n_reps_all_studies, &
                                                      max_n_genes_all_studies, n_neighbors, 1_int32, &
                                                      neighborhood_indices(:, i_point, i_study), &
                                                      pooled_residuals((i_study - 1)*block_size + 1:i_study*block_size), &
                                                      gather_ierr)
            if (is_err(gather_ierr)) then
                call set_err_once(ierr, get_err_code(gather_ierr))
                return
            end if
        end do
    end subroutine gather_fixed_k_point_pool

    !> Pass C of the fixed-k parameter search for one candidate: every study's histograms, by
    !| `build_fixed_k_study_histograms`. Private helper; every array dummy is a search buffer
    !| viewed packed at the candidate's sizes.
    pure subroutine build_fixed_k_candidate_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                       n_studies, n_neighbors, n_points, max_n_bins, neighborhood_indices, &
                                                       shared_residual_range_low, shared_residual_range_high, &
                                                       n_bins_per_point, counts, pmfs, included_n_reps, &
                                                       tmp_neighborhood_residuals_gathered)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_neighbors
            !! The candidate's number of neighbors per reference point
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! The candidate's widest per-point bin count
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices(n_neighbors, n_points, n_studies)
            !! Every study's neighbor gene indices
        real(real64), intent(in) :: shared_residual_range_low(n_points)
            !! Lower histogram range bound (R_low) per reference point
        real(real64), intent(in) :: shared_residual_range_high(n_points)
            !! Upper histogram range bound (R_high) per reference point
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(out) :: counts(n_points, max_n_bins, n_studies)
            !! Every study's histogram counts
        real(real64), intent(out) :: pmfs(n_points, max_n_bins, n_studies)
            !! Every study's pmfs
        integer(int32), intent(out) :: included_n_reps(n_points, n_studies)
            !! Every study's included-replicate counts
        real(real64), intent(out) :: tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points)
            !! Working array: one study's gathered residuals, reused per study

        integer(int32) :: i_study

        do i_study = 1, n_studies
            call build_fixed_k_study_histograms(residuals(:, :, i_study), max_n_reps_all_studies, max_n_genes_all_studies, &
                                                n_neighbors, n_points, neighborhood_indices(:, :, i_study), &
                                                shared_residual_range_low, shared_residual_range_high, n_bins_per_point, &
                                                max_n_bins, counts(:, :, i_study), pmfs(:, :, i_study), &
                                                included_n_reps(:, i_study), tmp_neighborhood_residuals_gathered)
        end do
    end subroutine build_fixed_k_candidate_histograms

    !> The second admissibility gate and the observed JSD of one parameter-search candidate,
    !| shared by both searches. First, every study must include at least one residual at every
    !| reference point (`has_empty_study` otherwise, and nothing else is computed): the
    !| equal-weight consensus divides by `n_studies`, so an empty study would leave it
    !| unnormalised. Then the consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]) must
    !| reach `min_count` in every bin, scoped to each point's own bin count
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| Only if both pass are each study's per-point JSD, weights and observed global JSD computed,
    !| and the confidence interval seeded with the observed value. Private helper; every array
    !| dummy is a search buffer viewed packed at the candidate's sizes.
    pure subroutine score_candidate_histograms(n_points, max_n_bins, n_studies, n_bins_per_point, min_count, pmfs, counts, &
                                               included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps, &
                                               js_divergences, weights, global_js_divergence, confidence_interval, &
                                               has_empty_study, all_bins_have_min_count)
        integer(int32), intent(in) :: n_points
            !! The candidate's number of reference points
        integer(int32), intent(in) :: max_n_bins
            !! The candidate's widest per-point bin count
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: n_bins_per_point(n_points)
            !! Bin count per reference point
        integer(int32), intent(in) :: min_count
            !! Minimum count every consensus bin must reach
        real(real64), intent(in) :: pmfs(n_points, max_n_bins, n_studies)
            !! Every study's pmfs
        integer(int32), intent(in) :: counts(n_points, max_n_bins, n_studies)
            !! Every study's histogram counts
        integer(int32), intent(in) :: included_n_reps(n_points, n_studies)
            !! Every study's included-replicate counts
        real(real64), intent(out) :: mean_pmf(n_points, max_n_bins)
            !! The consensus pmf
        integer(int32), intent(out) :: mean_pmf_counts(n_points, max_n_bins)
            !! The consensus histogram counts
        integer(int32), intent(out) :: mean_pmf_included_n_reps(n_points)
            !! The consensus included-replicate counts
        real(real64), intent(out) :: js_divergences(n_points, n_studies)
            !! Each study's per-point JSD against the consensus
        real(real64), intent(out) :: weights(n_points, n_studies)
            !! Each study's per-point weights
        real(real64), intent(out) :: global_js_divergence(n_studies)
            !! Each study's observed weighted global JSD
        real(real64), intent(out) :: confidence_interval(2, n_studies)
            !! `[J, J]` per study, the observed global JSD, to seed the bootstrap
        logical(c_bool), intent(out) :: has_empty_study
            !! `.true.` iff some study has no included residual at some reference point
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` iff every consensus bin reaches `min_count`; `.false.` when `has_empty_study`

        integer(int32) :: i_study

        has_empty_study = logical(any(included_n_reps == 0_int32), kind=c_bool)
        all_bins_have_min_count = logical(.false., kind=c_bool)
        if (has_empty_study) return

        call create_mean_pmf_impl(pmfs, counts, max_n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                                  mean_pmf_included_n_reps, mean_pmf_counts)
        call check_mean_pmf_min_counts_impl(mean_pmf_counts, max_n_bins, n_bins_per_point, n_points, min_count, &
                                            all_bins_have_min_count)
        if (.not. all_bins_have_min_count) return

        do i_study = 1, n_studies
            call compute_divergence_per_reference_point_impl(pmfs(:, :, i_study), mean_pmf, n_points, max_n_bins, &
                                                             js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
            confidence_interval(1, i_study) = global_js_divergence(i_study)
            confidence_interval(2, i_study) = global_js_divergence(i_study)
        end do
    end subroutine score_candidate_histograms

    !> summary: Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting
    !| AUTHOR_LASZLO_LANG
    !| Walks the candidates of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !|    ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !|    study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !|    ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !|    residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !|    `min_residuals_per_bin`
    !|    ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !|    A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !|    ([[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]),
    !|    and both plateau criteria: CI overlap with the best candidate so far
    !|    ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]])
    !|    and relative effect-size stability
    !|    ([[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
    !|    Both are always computed and traced; `plateau_mode` decides which one stops the search. A
    !|    stop by the effect-size criterion returns the candidate that triggered it.
    !|
    !| The returned candidate:
    !|
    !| 1. A plateau was found, or the grid has a single candidate and it is admissible: the best
    !|    candidate, `plateau_established = .true.`.
    !| 2. No plateau, but some candidate was admissible: the admissible candidate with the smallest
    !|    median confidence-interval width across studies, `plateau_established = .false.`.
    !| 3. No candidate was admissible: the first (finest) candidate, with
    !|    `best_candidate_pair_confidence_interval = -1.0`, every bin count `m_min`, both residual
    !|    ranges `0.0`, and `plateau_established = .false.`.
    !|
    !| Per-point outputs and `trace_*` arrays are sized by the candidate bounds; see each argument
    !| for its meaningful part. A failure to create the bootstrap's random number generator is
    !| returned in `ierr` and stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_parameter_search_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                       gene_means, residuals, n_bootstraps, &
                                                       join_method, max_n_points_candidate, max_n_neighbors_candidate, &
                                                       n_bootstrapping_top_k_jsds, n_points, n_neighbors, n_bins_per_point, &
                                                       shared_residual_range_low, shared_residual_range_high, &
                                                       best_candidate_pair_confidence_interval, plateau_established, &
                                                       n_admissible_evaluated, &
                                                       trace_n_points, trace_n_neighbors, trace_global_js_divergence, &
                                                       trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                                       trace_ci_width_relative, trace_ci_overlap, trace_delta, &
                                                       trace_delta_median, trace_delta_max, trace_selected_n_bins, &
                                                       trace_occupancy_failed, &
                                                       trace_n_pooled_residuals, trace_min_bin_occupancy, &
                                                       trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                                       trace_sturges_bins, trace_fd_bins, &
                                                       trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                                       tmp_gene_means_perms, &
                                                       tmp_gene_means_perm_all, &
                                                       tmp_x_star, tmp_neighborhood_indices_all_studies, &
                                                       tmp_neighborhood_range, tmp_neighborhood_residuals_gathered, &
                                                       tmp_point_pool, tmp_point_draw, tmp_n_bins_per_point, &
                                                       tmp_shared_residual_range_low, tmp_shared_residual_range_high, &
                                                       tmp_pmfs, tmp_counts, &
                                                       tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                       tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, &
                                                       tmp_global_js_divergence, tmp_confidence_interval, &
                                                       tmp_bootstrapping_top_k_jsds, tmp_prev_global_js_divergence, &
                                                       tmp_delta_perm, tmp_best_uncertainty_confidence_interval, &
                                                       tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                       tmp_bin_counts_search, tmp_residual_to_study_perm, &
                                                       tmp_study_block_ends, tmp_study_block_sizes, tmp_occupancy_failed, &
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
                                                       two_sided_bootstrapping_significance_level, random_seed, &
                                                       bin_count_mode, ierr)
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
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_points_candidate
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_impl, JUST_INFO)
            !! DM_MIN(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)))
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_impl, JUST_INFO)
            !! DM_MIN(max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32)))
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_js_comp_test_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: n_points
            !! The finally chosen candidate's `n_points`
        integer(int32), intent(out) :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(int32), dimension(max_n_points_candidate), intent(out) :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` when no candidate was admissible
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The chosen candidate's bootstrapped JSD confidence interval; `-1.0` when no candidate
            !! was admissible
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise, also when the smallest-uncertainty
            !! fallback returns a real candidate and interval
        integer(int32), intent(out) :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array (an inadmissible candidate has no entry)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_points
            !! Per-admissible-candidate `n_points`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width (the values
            !! `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        logical(c_bool), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            !! like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            !! Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            !! Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            !! Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array: the current candidate's reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_range
            !! Working array: one study's neighborhood spans, reused per study
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_point_draw
            !! Working array of the bootstrap: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array: the current candidate's bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array: the current candidate's R_low per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array: the current candidate's R_high per reference point
        real(real64), dimension(max_n_points_candidate, CM_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's pmf for the current candidate
        integer(int32), dimension(max_n_points_candidate, CM_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's histogram counts for the current candidate
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(max_n_points_candidate, CM_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf
        integer(int32), dimension(max_n_points_candidate, CM_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, then bootstrap scratch
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-point weights, then bootstrap scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD, then bootstrap scratch
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
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array: the current candidate's occupancy_failed flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array: the current candidate's pooled residual count (N_j) per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array: the current candidate's minimum bin occupancy per point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array: the current candidate's mean bin occupancy per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array: the current candidate's maximum bin occupancy per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array: the current candidate's Sturges' rule diagnostic per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array: the current candidate's Freedman-Diaconis diagnostic per point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array: bin counts per point of the current best candidate
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: bin counts per point of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array: R_low per point of the current best candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array: R_high per point of the current best candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: R_low per point of the smallest-uncertainty admissible candidate
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: R_high per point of the smallest-uncertainty admissible candidate
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
            !! one threshold for both, so a point the search marks as failed also fails the gate
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.1_real64)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.9_real64)
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(CM_MODE_PLATEAU_CI_OVERLAP)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(CM_MODE_BIN_COUNT_OCCUPANCY)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created

        integer(int32) :: candidates_n_points_n_neighbors(2, MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, i_study, i_point, n_pool, max_n_bins, gather_ierr
        integer(int32) :: prev_n_points, best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: pool_size, bootstrap_ierr, actual_min_residuals_per_bin
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index, actual_m_min, actual_bin_count_mode
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value, median_ci_width
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found, has_empty_study
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
        M_DEFAULT_VAL(bin_count_mode, actual_bin_count_mode, CM_MODE_BIN_COUNT_OCCUPANCY)

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
            !
            ! The work arrays sized by an upper bound (max_n_neighbors_candidate,
            ! max_n_points_candidate, CM_MAX_N_BINS) are used at the candidate's own sizes. They are
            ! always passed whole to explicit-shape dummies declared with the candidate's sizes, so
            ! every routine sees their leading elements packed (sequence association) and nothing is
            ! copied; every use of such a buffer goes through such a dummy, so all agree on the
            ! layout. This is why Passes A to C and the second gate live in private helpers.
            call construct_fixed_k_candidate_neighborhoods(n_studies, max_n_genes_all_studies, n_points, n_neighbors, &
                                                           tmp_x_star, gene_means, tmp_gene_means_perms, &
                                                           actual_min_neighbor_overlap, &
                                                           tmp_neighborhood_indices_all_studies, tmp_neighborhood_range, &
                                                           all_have_min_neighbor_overlap)
            if (.not. all_have_min_neighbor_overlap) cycle

            ! ===== PASS B (per point): pool all studies' residuals and run the bin-count search.
            ! Sequential although the points are independent, so one tmp_pooled_residuals buffer
            ! serves every point (do concurrent would need one per point). Every study's block has
            ! max_n_reps_all_studies*n_neighbors values, so the block ends are the same for every point.
            do concurrent(i_study=1:n_studies) shared(tmp_study_block_ends, max_n_reps_all_studies, n_neighbors)
                tmp_study_block_ends(i_study) = i_study*max_n_reps_all_studies*n_neighbors
            end do
            do i_point = 1, n_points
                call gather_fixed_k_point_pool(residuals, max_n_reps_all_studies, max_n_genes_all_studies, n_studies, &
                                               n_neighbors, n_points, tmp_neighborhood_indices_all_studies, i_point, &
                                               tmp_pooled_residuals, gather_ierr)
                ! Defensive: only an out-of-range gene index fails here, and every index comes from
                ! permutations this routine builds. Kept for direct callers, since Pass C would read
                ! out of bounds; untested.
                if (is_err(gather_ierr)) then
                    call set_err_once(ierr, get_err_code(gather_ierr))
                    return
                end if

                call init_perm(tmp_pooled_residuals_perm(1:max_n_reps_all_studies*n_neighbors*n_studies))
                call sort_array_heapsort(tmp_pooled_residuals(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                                         tmp_pooled_residuals_perm(1:max_n_reps_all_studies*n_neighbors*n_studies))

                call select_point_bin_count( &
                    tmp_pooled_residuals(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                    tmp_pooled_residuals_perm(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                    max_n_reps_all_studies*n_neighbors*n_studies, n_studies, tmp_study_block_ends, &
                    max_n_reps_all_studies, n_neighbors, actual_bin_count_mode, &
                    tmp_n_bins_per_point(i_point), tmp_occupancy_failed(i_point), &
                    tmp_shared_residual_range_low(i_point), tmp_shared_residual_range_high(i_point), &
                    tmp_n_pooled_residuals(i_point), tmp_min_bin_occupancy(i_point), tmp_mean_bin_occupancy(i_point), &
                    tmp_max_bin_occupancy(i_point), tmp_sturges_bins(i_point), tmp_fd_bins(i_point), &
                    tmp_bin_counts_search, tmp_residual_to_study_perm(1:max_n_reps_all_studies*n_neighbors*n_studies), &
                    tmp_study_block_sizes, m_min, m_max, actual_min_residuals_per_bin, gamma_occupancy, &
                    lower_residual_range_quantile, upper_residual_range_quantile)
            end do

            ! ===== PASS C (per study): every study's histograms at the chosen bin counts.
            max_n_bins = maxval(tmp_n_bins_per_point(1:n_points))
            call build_fixed_k_candidate_histograms(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                    n_studies, n_neighbors, n_points, max_n_bins, &
                                                    tmp_neighborhood_indices_all_studies, tmp_shared_residual_range_low, &
                                                    tmp_shared_residual_range_high, tmp_n_bins_per_point, tmp_counts, &
                                                    tmp_pmfs, tmp_included_n_reps, tmp_neighborhood_residuals_gathered)

            ! Second admissibility gate; an empty study (an error in run_js_comp_test_impl) only
            ! makes the candidate inadmissible here.
            call score_candidate_histograms(n_points, max_n_bins, n_studies, tmp_n_bins_per_point, &
                                            actual_min_residuals_per_bin, tmp_pmfs, tmp_counts, tmp_included_n_reps, &
                                            tmp_mean_pmf, tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                            tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                            tmp_confidence_interval, has_empty_study, all_bins_have_min_count)
            if (has_empty_study .or. .not. all_bins_have_min_count) cycle

            ! Trace the candidate and test the effect-size criterion, whatever plateau_mode is. Must
            ! come before the bootstrap, which reuses tmp_global_js_divergence as scratch.
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

            ! Same n_bootstraps/random_seed for every candidate, for comparability. The observed-value
            ! arrays above double as the bootstrap's work arrays, passed whole and viewed packed.
            call bootstrap_histogram_impl(n_bootstraps, max_n_bins, n_points, n_studies, tmp_mean_pmf_counts, &
                                          tmp_mean_pmf_included_n_reps, tmp_included_n_reps, &
                                          n_bootstrapping_top_k_jsds, tmp_confidence_interval, &
                                          tmp_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, tmp_mean_pmf, &
                                          tmp_js_divergences, tmp_weights, tmp_global_js_divergence, tmp_point_pool, &
                                          tmp_point_draw, two_sided_bootstrapping_significance_level, random_seed, &
                                          bootstrap_ierr)
            ! Defensive: the counts add up by construction, so only a failed rng allocation can
            ! land here (untested). Without the return the interval would stay at its zero-width
            ! [J, J] seed and could win the smallest-uncertainty ranking.
            if (is_err(bootstrap_ierr)) then
                call set_err_once(ierr, get_err_code(bootstrap_ierr))
                return
            end if

            trace_ci_lower(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(1, 1:n_studies)
            trace_ci_upper(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(2, 1:n_studies)
            trace_ci_width(1:n_studies, n_admissible_evaluated) = &
                trace_ci_upper(1:n_studies, n_admissible_evaluated) - trace_ci_lower(1:n_studies, n_admissible_evaluated)
            trace_ci_width_relative(1:n_studies, n_admissible_evaluated) = &
                trace_ci_width(1:n_studies, n_admissible_evaluated) &
                / max(trace_global_js_divergence(1:n_studies, n_admissible_evaluated), actual_delta_epsilon)

            ! Track the admissible candidate with the smallest median CI width (the no-plateau
            ! fallback). tmp_delta_perm must be re-sorted for trace_ci_width: calc_quantile_impl
            ! trusts the order it is handed.
            call init_perm(tmp_delta_perm)
            call sort_array_heapsort(trace_ci_width(1:n_studies, n_admissible_evaluated), tmp_delta_perm)
            call calc_quantile_impl(trace_ci_width(1:n_studies, n_admissible_evaluated), n_studies, tmp_delta_perm, &
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
                                              actual_succeeding_ci_overlap, &
                                              trace_ci_overlap(1:n_studies, n_admissible_evaluated), ci_plateau_found)
            ! The first admissible candidate was compared against the -1.0 seed: record -1.0.
            if (n_admissible_evaluated == 1_int32) trace_ci_overlap(1:n_studies, 1) = -1.0_real64

            ! This candidate became the new best: snapshot its per-point values.
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

            ! The best-candidate bookkeeping follows the CI-overlap criterion only; an effect-size
            ! plateau returns the candidate that triggered it.
            if (effect_size_plateau_found .and. actual_plateau_mode /= MODE_PLATEAU_CI_OVERLAP) then
                best_candidate_index = i_candidate
                best_candidate_pair_confidence_interval = tmp_confidence_interval
                tmp_best_n_bins_per_point(1:n_points) = tmp_n_bins_per_point(1:n_points)
                tmp_best_shared_residual_range_low(1:n_points) = tmp_shared_residual_range_low(1:n_points)
                tmp_best_shared_residual_range_high(1:n_points) = tmp_shared_residual_range_high(1:n_points)
            end if

            if (plateau_found) exit
        end do

        ! The three cases of the routine description. In case 1 some candidate was admissible, so
        ! the best-candidate snapshot exists; n_candidates counts gate-failed candidates too.
        if (plateau_found .or. (n_candidates < 2_int32 .and. n_admissible_evaluated >= 1_int32)) then
            n_points = candidates_n_points_n_neighbors(1, best_candidate_index)
            n_neighbors = candidates_n_points_n_neighbors(2, best_candidate_index)
            n_bins_per_point(1:n_points) = tmp_best_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = tmp_best_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = tmp_best_shared_residual_range_high(1:n_points)
            plateau_established = logical(.true., kind=c_bool)
        else if (n_admissible_evaluated >= 1_int32) then
            n_points = candidates_n_points_n_neighbors(1, best_uncertainty_candidate_index)
            n_neighbors = candidates_n_points_n_neighbors(2, best_uncertainty_candidate_index)
            n_bins_per_point(1:n_points) = tmp_best_uncertainty_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = tmp_best_uncertainty_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = tmp_best_uncertainty_shared_residual_range_high(1:n_points)
            best_candidate_pair_confidence_interval = tmp_best_uncertainty_confidence_interval
            plateau_established = logical(.false., kind=c_bool)
        else
            n_points = candidates_n_points_n_neighbors(1, 1)
            n_neighbors = candidates_n_points_n_neighbors(2, 1)
            ! No bin count was kept for candidate 1, so m_min (as in a failed bin-count search) and
            ! a zero range (as for an all-NaN pool).
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
    !| 1. Grow its neighborhoods exactly as
    !|    [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !|    does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    !|    `min_study_neighbors`.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !|    reference points emerged.
    !| 3. Gate 1:
    !|    [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !|    on the pooled ranges.
    !| 4. The pipeline of
    !|    [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !|    up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    !|    (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    !|    count
    !|    ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    !|    first plateau `plateau_mode` selects.
    !|
    !| **Candidate log.** Every candidate tried gets an entry in `candidate_k_start`,
    !| `candidate_n_points` and `candidate_status`, until the search stops. The status is the first
    !| failing check, in the order above:
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED(variable)]]
    !| (too few means or residuals),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_EMPTY_STUDY(variable)]]
    !| (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
    !| some point),
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED(variable)]],
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_OVERLAP_FAILED(variable)]],
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED(variable)]],
    !| or else [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_CANDIDATE_EVALUATED(variable)]]
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
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]);
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
            !! Number of values kept at each end of the bootstrap distribution (heap size)
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
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful. Without an admissible candidate, the first
            !! `min(n_points, max_n_points_candidate)` entries are `m_min`
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` without an admissible candidate
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
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_k_start
            !! Per-admissible-candidate `k_start`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width; `-1.0` for the
            !! first admissible candidate
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(n_studies, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        logical(c_bool), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            !! trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            !! as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        real(real64), dimension(max_n_points_candidate, CM_MAX_CANDIDATE_PAIRS), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! DM_RESULT_SIZE_IS(n_admissible_evaluated)
        integer(int32), intent(out) :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: candidate_n_points
            !! Number of reference points each tried candidate's construction produced (0 for too
            !! few means)
            !! DM_RESULT_SIZE_IS(n_candidates_tried)
        integer(int32), dimension(CM_MAX_CANDIDATE_PAIRS), intent(out) :: candidate_status
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
        real(real64), dimension(max_n_points_candidate, CM_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array: every study's pmf for the current candidate
        integer(int32), dimension(max_n_points_candidate, CM_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array: every study's histogram counts for the current candidate
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(real64), dimension(max_n_points_candidate, CM_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf
        integer(int32), dimension(max_n_points_candidate, CM_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: the current candidate's per-point JSD values, then bootstrap scratch
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array: the current candidate's per-point weights, then bootstrap scratch
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD, then bootstrap scratch
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_point_draw
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
        integer(int32), dimension(CM_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
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
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.1_real64)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap of a candidate's confidence interval with the best one so
            !! far, per `join_method`. `0` makes the first admissible candidate plateau under the
            !! CI-overlap criterion (for baseline sweeps)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(0.9_real64)
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(CM_MODE_PLATEAU_CI_OVERLAP)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Floor of the denominators of the relative JSD change and relative CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(CM_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions for an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(CM_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(CM_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(CM_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(CM_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(CM_MODE_BIN_COUNT_OCCUPANCY)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created

        integer(int32) :: candidates_k_start_k_step_k_max(3, MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, i_point, max_n_bins, point_ierr
        integer(int32) :: best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: pool_size, bootstrap_ierr, actual_min_residuals_per_bin
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index, actual_m_min
        integer(int32) :: actual_max_pooled_residuals, actual_min_study_neighbors, construction_status, max_n_neighbors
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value, median_ci_width, actual_tau, actual_mad_distance_factor
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found, has_empty_study
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
                                               max_n_genes_all_studies, 1_int32, tmp_point_neighborhood_indices, 1_int32, &
                                               tmp_n_neighbors_per_point(:, i_point), tmp_n_bins_per_point(i_point), &
                                               tmp_occupancy_failed(i_point), tmp_shared_residual_range_low(i_point), &
                                               tmp_shared_residual_range_high(i_point), tmp_n_pooled_residuals(i_point), &
                                               tmp_min_bin_occupancy(i_point), tmp_mean_bin_occupancy(i_point), &
                                               tmp_max_bin_occupancy(i_point), tmp_sturges_bins(i_point), &
                                               tmp_fd_bins(i_point), tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                               tmp_bin_counts_search, tmp_residual_to_study_perm, tmp_study_block_ends, &
                                               tmp_study_block_sizes, point_ierr, m_min=m_min, m_max=m_max, &
                                               min_residuals_per_bin=actual_min_residuals_per_bin, &
                                               gamma_occupancy=gamma_occupancy, &
                                               lower_residual_range_quantile=lower_residual_range_quantile, &
                                               upper_residual_range_quantile=upper_residual_range_quantile, &
                                               bin_count_mode=bin_count_mode)
                ! Defensive: determine_point_bin_count only fails on an out-of-range gene index, and
                ! every decoded index lies in [1, max_n_genes_all_studies], so this branch cannot be
                ! reached and has no test.
                if (is_err(point_ierr)) then
                    call set_err_once(ierr, get_err_code(point_ierr))
                    return
                end if
            end do

            ! ===== PASS C (per point, per study, sequential: one gather buffer serves every pair).
            ! From here on the work arrays sized by an upper bound (max_n_points_candidate,
            ! CM_MAX_N_BINS) are passed whole to explicit-shape dummies declared with the
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

            ! From here on as in the fixed-k search. Trace and test the effect-size criterion before
            ! the bootstrap, which reuses tmp_global_js_divergence as scratch.
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
            call bootstrap_histogram_impl(n_bootstraps, max_n_bins, n_points, n_studies, tmp_mean_pmf_counts, &
                                          tmp_mean_pmf_included_n_reps, tmp_included_n_reps, &
                                          n_bootstrapping_top_k_jsds, tmp_confidence_interval, &
                                          tmp_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, tmp_mean_pmf, &
                                          tmp_js_divergences, tmp_weights, tmp_global_js_divergence, tmp_point_pool, &
                                          tmp_point_draw, two_sided_bootstrapping_significance_level, random_seed, &
                                          bootstrap_ierr)
            ! Defensive: the counts add up by construction, so only a failed rng allocation can
            ! land here (untested). Without the return the interval would stay at its zero-width
            ! [J, J] seed and could win the smallest-uncertainty ranking.
            if (is_err(bootstrap_ierr)) then
                call set_err_once(ierr, get_err_code(bootstrap_ierr))
                return
            end if

            trace_ci_lower(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(1, 1:n_studies)
            trace_ci_upper(1:n_studies, n_admissible_evaluated) = tmp_confidence_interval(2, 1:n_studies)
            trace_ci_width(1:n_studies, n_admissible_evaluated) = &
                trace_ci_upper(1:n_studies, n_admissible_evaluated) - trace_ci_lower(1:n_studies, n_admissible_evaluated)
            trace_ci_width_relative(1:n_studies, n_admissible_evaluated) = &
                trace_ci_width(1:n_studies, n_admissible_evaluated) &
                / max(trace_global_js_divergence(1:n_studies, n_admissible_evaluated), actual_delta_epsilon)

            ! Track the admissible candidate with the smallest median CI width (the no-plateau
            ! fallback). tmp_delta_perm must be re-sorted for trace_ci_width.
            call init_perm(tmp_delta_perm)
            call sort_array_heapsort(trace_ci_width(1:n_studies, n_admissible_evaluated), tmp_delta_perm)
            call calc_quantile_impl(trace_ci_width(1:n_studies, n_admissible_evaluated), n_studies, tmp_delta_perm, &
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
                                              actual_succeeding_ci_overlap, &
                                              trace_ci_overlap(1:n_studies, n_admissible_evaluated), ci_plateau_found)
            ! The first admissible candidate was compared against the -1.0 seed: record -1.0.
            if (n_admissible_evaluated == 1_int32) trace_ci_overlap(1:n_studies, 1) = -1.0_real64

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
        ! is the chosen candidate's logged one; the knobs come from the sequence. A single
        ! candidate counts as established only when it was admissible; a rejected one falls
        ! through to case 3, which returns that same candidate.
        if (plateau_found .or. (n_candidates < 2_int32 .and. n_admissible_evaluated >= 1_int32)) then
            k_start = candidates_k_start_k_step_k_max(1, best_candidate_index)
            k_step = candidates_k_start_k_step_k_max(2, best_candidate_index)
            k_max = candidates_k_start_k_step_k_max(3, best_candidate_index)
            n_points = candidate_n_points(best_candidate_index)
            n_bins_per_point(1:n_points) = tmp_best_n_bins_per_point(1:n_points)
            shared_residual_range_low(1:n_points) = tmp_best_shared_residual_range_low(1:n_points)
            shared_residual_range_high(1:n_points) = tmp_best_shared_residual_range_high(1:n_points)
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
            ! Nothing admissible (including a single, rejected candidate): candidate 1, with m_min
            ! and a zero range as in the fixed-k search. Its point count may exceed the capacity,
            ! hence the clamp.
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
