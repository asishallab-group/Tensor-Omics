#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The data-driven choice of the JS-Comp-Test's neighborhood parameters: a candidate is admissible
!| when its neighborhoods overlap enough and every consensus bin holds enough residuals, and the
!| search stops at the first admissible candidate whose bootstrapped JSDs have reached a plateau.
!|
!| Fixed-k neighborhoods:
!| [[tox_data_integration_js_comp_test_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| searches the `(n_points, n_neighbors)` candidates of
!| [[tox_data_integration_js_comp_test_search_impl(module):generate_js_comp_test_candidates_impl(interface)]].
!| Adaptive neighborhoods:
!| [[tox_data_integration_js_comp_test_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| searches the `(k_start, k_step, k_max)` candidates of
!| [[tox_data_integration_js_comp_test_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
!| The building blocks are published too: the gate
!| [[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]
!| and the plateau checks
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds` and `calc_adaptive_js_comp_test_bounds` size the work
!| arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_search_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_math_impl, only: clamp
    use f42_stats_impl, only: calc_quantile_impl
    use f42_sort_impl, only: init_perm, sort_array_heapsort, sort_real_heapsort_expl_size
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, validate_in_range_int, ERR_INVALID_INPUT
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, create_mean_pmf_impl
    use tox_data_integration_preprocessing_impl, only: construct_neighborhoods_ranged_impl, pool_means_impl, &
                                                       gather_pooled_neighborhood_residuals, &
                                                       grow_adaptive_neighborhoods, materialize_pooled_neighborhood, &
                                                       check_neighborhood_overlaps_impl, ADAPTIVE_STATUS_OK, &
                                                       ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD
    use tox_data_integration_stats_impl, only: bootstrap_histogram_impl
    use tox_data_integration_js_comp_test_binning_impl, only: select_point_bin_count, determine_point_bin_count, &
                                                              build_fixed_k_study_histograms, build_adaptive_histograms
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_candidate_bounds, calc_adaptive_js_comp_test_bounds, &
             generate_js_comp_test_candidates_impl, generate_adaptive_js_comp_test_candidates_impl, &
             check_mean_pmf_min_counts_impl, check_plateau_condition_impl, check_effect_size_plateau_condition_impl, &
             run_js_comp_test_parameter_search_impl, run_js_comp_test_adaptive_parameter_search_impl, &
             METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, &
             MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, MODE_PLATEAU_BOTH, KX_FACTORS, MAX_POINTS, MIN_POINTS, &
             GAMMA, MAX_POINT_CANDIDATES, MAX_CANDIDATE_PAIRS, &
             ADAPTIVE_START_FRACTION, ADAPTIVE_GAMMA, ADAPTIVE_K_STEP_FRACTION, ADAPTIVE_K_MAX_FACTOR, &
             ADAPTIVE_K_START_MIN_ABS, ADAPTIVE_POINT_CAPACITY_FACTOR, &
             ADAPTIVE_CANDIDATE_EVALUATED, ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED, ADAPTIVE_CANDIDATE_EMPTY_STUDY, &
             ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED, ADAPTIVE_CANDIDATE_OVERLAP_FAILED, &
             ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED

    ! `join_method`'s mode table; the `METHOD_` prefix follows from the argument name `join_method`.
    integer(int32), parameter :: METHOD_JOIN_MIN = 0_int32
        !! Join method: succeeds only once every study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MAX = 1_int32
        !! Join method: succeeds once any one study's confidence-interval overlap exceeds the
        !! threshold, in [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]
    integer(int32), parameter :: METHOD_JOIN_MEDIAN = 2_int32
        !! Join method: succeeds once a majority of studies' confidence-interval overlaps exceed
        !! the threshold, in [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]

    ! Tuning constants of the fixed-k candidate grid (generate_js_comp_test_candidates_impl).
    real(real64), parameter :: KX_FACTORS(2) = [0.25_real64, 0.5_real64]
        !! Neighbor-count factors, ascending: the first gives the larger `n_neighbors` candidate
    integer(int32), parameter :: MAX_POINTS = 1500_int32
    integer(int32), parameter :: MIN_POINTS = 300_int32
    real(real64), parameter :: GAMMA = 0.8_real64
    integer(int32), parameter :: MAX_POINT_CANDIDATES = 8_int32
        !! How often n_points_high may be multiplied by GAMMA before the grid gives up
    ! CM_MAX_CANDIDATE_PAIRS is the single source of this value, used both for the parameter and
    ! for the dummy extents of wrapped routines: a parameter cannot be such an extent, because the
    ! generated C and R bindings need the number itself.
#define CM_MAX_CANDIDATE_PAIRS 16
    integer(int32), parameter :: MAX_CANDIDATE_PAIRS = CM_MAX_CANDIDATE_PAIRS
        !! = size(KX_FACTORS) * MAX_POINT_CANDIDATES = 2 * 8

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

    ! `plateau_mode`'s mode table. CM_MODE_PLATEAU_CI_OVERLAP is also its DM_DEFAULT, which cannot
    ! name a Fortran parameter.
#define CM_MODE_PLATEAU_CI_OVERLAP 0_int32
    integer(int32), parameter :: MODE_PLATEAU_CI_OVERLAP = CM_MODE_PLATEAU_CI_OVERLAP
        !! Plateau mode: CI overlap only
    integer(int32), parameter :: MODE_PLATEAU_EFFECT_SIZE = 1_int32
        !! Plateau mode: relative-effect-size stability only
    integer(int32), parameter :: MODE_PLATEAU_BOTH = 2_int32
        !! Plateau mode: either CI overlap or relative-effect-size stability

    ! Defaults of the relative-effect-size plateau criterion (provisional, to be validated
    ! empirically), shared by DM_DEFAULT and M_DEFAULT_VAL.
#define CM_DELTA_MEDIAN_THRESHOLD_DEFAULT 0.05_real64
#define CM_DELTA_MAX_THRESHOLD_DEFAULT 0.10_real64
#define CM_DELTA_EPSILON_DEFAULT 1.0e-10_real64
#define CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT 2_int32

contains

    !> M_EXPORT_C
    !| summary: Recommend upper bounds for the js-comp-test candidate-grid work arrays
    !| AUTHOR_LASZLO_LANG
    !| For the candidate grid of
    !| [[tox_data_integration_js_comp_test_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]:
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| A candidate's reference-point count is only known after its neighborhoods are grown, so this
    !| is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):generate_js_comp_test_candidates_impl(interface)]].
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
    !| The constants are [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_START_FRACTION(variable)]] (0.02),
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_GAMMA(variable)]] (0.8),
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_K_STEP_FRACTION(variable)]] (0.25),
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]] (4),
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_K_START_MIN_ABS(variable)]] (10) and
    !| [[tox_data_integration_js_comp_test_search_impl(module):MAX_CANDIDATE_PAIRS(variable)]] (CM_MAX_CANDIDATE_PAIRS).
    !| The values strictly descend without deduplication. At least one candidate is always
    !| returned, even for data smaller than `k_start_1`, whose construction then fails with a
    !| too-few-means status in the search.
    !|
    !| Not checked here: the values fit 32-bit integers only while `k_start_1 <= huge(1_int32)/4`;
    !| [[tox_data_integration_js_comp_test_search_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
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

    !> summary: Test whether every bin of a mean pmf reaches a minimum absolute count
    !| AUTHOR_LASZLO_LANG
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]].
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]].
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
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]]:
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

    !> Pass A of the fixed-k parameter search for one candidate: builds every study's
    !| neighborhoods into `neighborhood_indices` and checks the first admissibility gate
    !| ([[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]])
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
    !| [[tox_data_integration_preprocessing_impl(module):gather_pooled_neighborhood_residuals(interface)]]
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
    !| ([[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]]) must
    !| reach `min_count` in every bin, scoped to each point's own bin count
    !| ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !|    ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !|    study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !|    ([[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !|    residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !|    `min_residuals_per_bin`
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !|    A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !|    ([[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]]),
    !|    and both plateau criteria: CI overlap with the best candidate so far
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]])
    !|    and relative effect-size stability
    !|    ([[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
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
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (a majority must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_points_candidate
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_search_impl, JUST_INFO)
            !! DM_MIN(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)))
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_search_impl, JUST_INFO)
            !! DM_MIN(max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32)))
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_stats_impl, AUTO)
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
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
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
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !| Instead of a fixed `(n_points, n_neighbors)` grid, it walks the growth-knob candidates of
    !| [[tox_data_integration_js_comp_test_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
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
    !|    [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_CANDIDATE_CONSTRUCTION_FAILED(variable)]]
    !| (too few means or residuals),
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_CANDIDATE_EMPTY_STUDY(variable)]]
    !| (a study neighborhood below `min_study_neighbors`, or a study without a non-NaN residual at
    !| some point),
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_CANDIDATE_CAPACITY_EXCEEDED(variable)]],
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_CANDIDATE_OVERLAP_FAILED(variable)]],
    !| [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_CANDIDATE_MIN_COUNT_FAILED(variable)]],
    !| or else [[tox_data_integration_js_comp_test_search_impl(module):ADAPTIVE_CANDIDATE_EVALUATED(variable)]]
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
    !| [[tox_data_integration_js_comp_test_search_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]);
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
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_bounds, tox_data_integration_js_comp_test_search_impl, JUST_INFO)
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
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
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
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
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

end module tox_data_integration_js_comp_test_search_impl
