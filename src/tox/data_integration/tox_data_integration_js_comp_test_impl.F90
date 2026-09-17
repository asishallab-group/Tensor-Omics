#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
!| before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
!| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]],
!| each candidate's histogram bin count from
!| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]]), two
!| admissibility gates a candidate must pass before it is bootstrapped
!| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]),
!| and the plateau check that decides when the search has converged
!| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
!| `calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
!| allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| Later stages of the port add the K-study permutation test and the top-level orchestrator that
!| wire these building blocks together.
module tox_data_integration_js_comp_test_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_math_impl, only: clamp, is_close, LOG_2
    use f42_stats_impl, only: calc_percentile_impl
    use f42_sort_impl, only: init_perm, sort_array_heapsort, sort_real_heapsort_expl_size
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_multinomial
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, &
                                             build_residual_histograms_impl, calc_pmf_impl
    use tox_data_integration_preprocessing_impl, only: construct_neighborhoods_ranged_impl, pool_means_impl
    use tox_data_integration_stats_impl, only: gjct_permutation_test_impl
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_candidate_bounds, estimate_bin_count_impl, determine_bin_count_occupancy_impl, &
             generate_js_comp_test_candidates_impl, &
             check_neighborhood_overlaps_impl, check_mean_pmf_min_counts_impl, check_plateau_condition_impl, &
             check_effect_size_plateau_condition_impl, create_mean_pmf_impl, create_mean_pmf_only_impl, &
             calc_js_comp_test_n_top_k_jsds, bootstrap_histogram_impl, run_js_comp_test_impl, &
             run_js_comp_test_parameter_search_impl, METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, &
             MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, MODE_PLATEAU_BOTH, KX_FACTORS, MAX_POINTS, MIN_POINTS, &
             GAMMA, MAX_POINT_CANDIDATES, MAX_CANDIDATE_PAIRS, MAX_N_BINS

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
        !! bound written on candidates_n_points_n_neighbors/n_bins_candidates below -- the two
        !! must agree, and a compile-time PARAMETER expression cannot itself be documented as
        !! equal to a dimension literal the generator reads as plain text.
    integer(int32), parameter :: MAX_N_BINS = 256_int32
        !! Fixed ceiling for a candidate's estimated histogram bin count (see
        !! estimate_bin_count_impl below): later stages of this pipeline size their histogram
        !! work arrays to this fixed bound rather than to a data-dependent one.

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
    pure subroutine calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, &
                                                        max_n_neighbors_candidate)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Exact upper bound on the `n_points` candidate the grid ever produces
        integer(int32), intent(out) :: max_n_neighbors_candidate
            !! Safe upper bound on the `n_neighbors` candidate the grid ever produces

        integer(int32) :: n_points_low

        max_n_points_candidate = clamp(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), &
                                       min_val=MIN_POINTS, max_val=MAX_POINTS)
        n_points_low = max(MIN_POINTS, ceiling(0.2_real64*real(max_n_points_candidate, real64), kind=int32))
        max_n_neighbors_candidate = max(1_int32, floor(real(max_n_genes_all_studies, real64)/ &
                                                        (KX_FACTORS(1)*real(n_points_low, real64)), kind=int32))
    end subroutine calc_js_comp_test_candidate_bounds

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
    !| over `[-shared_residual_range, shared_residual_range]` has every bin at or above
    !| `min_residuals_per_bin` (the occupancy criterion), rather than the generic
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
                                                        max_n_reps_all_studies, n_neighbors, shared_residual_range, &
                                                        selected_n_bins, occupancy_failed, n_pooled_residuals, &
                                                        min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, &
                                                        sturges_bins, fd_bins, tmp_bin_counts, m_min, m_max, &
                                                        min_residuals_per_bin, gamma_occupancy)
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
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! DM_MIN(0.0_real64)
        integer(int32), intent(out) :: selected_n_bins
            !! The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
            !! every bin at or above min_residuals_per_bin; m_min when occupancy_failed
        logical(c_bool), intent(out) :: occupancy_failed
            !! `.true.` iff even m_min bins could not satisfy the occupancy criterion (including
            !! the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
            !! caller should reject this neighborhood rather than build a histogram from
            !! selected_n_bins
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
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(CM_OCCUPANCY_GAMMA_DEFAULT)

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        real(real64) :: actual_gamma_occupancy
        integer(int32) :: i_pool, trial_m, next_m, m_valid, m_invalid, best_m, min_occ, best_min_occ, best_max_occ
        integer(int32) :: combined_n_bins_unused
        logical(c_bool) :: found_valid

        M_DEFAULT_VAL(m_min, actual_m_min, CM_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, CM_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, CM_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(gamma_occupancy, actual_gamma_occupancy, CM_OCCUPANCY_GAMMA_DEFAULT)
        actual_m_max = max(actual_m_min, actual_m_max)

        ! Diagnostics only -- estimate_bin_count_impl's own all-NaN branch (n_bins=sturges_bins=
        ! fd_bins=1) already handles n_pool==0 correctly, so this call is unconditional.
        call estimate_bin_count_impl(pooled_residuals, pooled_residuals_perm, n_residuals, max_n_reps_all_studies, &
                                     n_neighbors, shared_residual_range, combined_n_bins_unused, sturges_bins, fd_bins)

        ! NaN sorts last under the ascending permutation -- same scan-back convention as
        ! estimate_bin_count_impl / determine_shared_residual_range_impl.
        n_pooled_residuals = n_residuals
        do i_pool = n_pooled_residuals, 1, -1
            if (ieee_is_nan(pooled_residuals(pooled_residuals_perm(i_pool)))) then
                n_pooled_residuals = n_pooled_residuals - 1_int32
            else
                exit
            end if
        end do

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
                                      shared_residual_range, trial_m, tmp_bin_counts(1:trial_m), min_occ)

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
                                      shared_residual_range, trial_m, tmp_bin_counts(1:trial_m), min_occ)
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
    !| `[-shared_residual_range, shared_residual_range]`, mirroring
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
                                         shared_residual_range, trial_m, bin_counts, min_occupancy)
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        integer(int32), intent(in) :: n_pool
            !! Number of non-NaN pooled residuals -- only pooled_residuals_perm(1:n_pool) is used
        real(real64), intent(in) :: pooled_residuals(n_residuals)
            !! Pooled signed residuals for one neighborhood
        integer(int32), intent(in) :: pooled_residuals_perm(n_residuals)
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
        integer(int32), intent(in) :: trial_m
            !! Candidate number of equally sized histogram bins in range [-R,R]
        integer(int32), intent(out) :: bin_counts(trial_m)
            !! Absolute count of pooled residuals per bin
        integer(int32), intent(out) :: min_occupancy
            !! minval(bin_counts) -- the occupancy criterion's own test statistic

        real(real64) :: bin_width, clamped_residual
        integer(int32) :: i_pool, bin_idx

        ! Same zero-range guard as build_residual_histograms_impl: a degenerate shared_residual_range
        ! (e.g. every pooled residual is 0) falls back to a fixed bin width instead of dividing by zero.
        if (shared_residual_range <= 0.0_real64) then
            bin_width = 1.0_real64
        else
            bin_width = 2.0_real64*shared_residual_range/real(trial_m, real64)
        end if

        bin_counts = 0_int32
        do i_pool = 1, n_pool
            clamped_residual = clamp(pooled_residuals(pooled_residuals_perm(i_pool)), min_val=-shared_residual_range, &
                                     max_val=shared_residual_range)
            bin_idx = min(trial_m, int((clamped_residual + shared_residual_range)/bin_width) + 1_int32)
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1_int32
        end do

        min_occupancy = minval(bin_counts)
    end subroutine histogram_bin_counts

    !> Pools one reference point's residuals across all studies from each study's already-known
    !| neighbor gene indices for that point -- the same gather
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] already
    !| performs per study inline (its own `do concurrent` gather from `neighborhood_indices` into
    !| `tmp_neighborhood_residuals_gathered`), generalized here to pool across every study at once
    !| into one flat array instead of overwriting one reused per-study buffer. Mirrors the nested
    !| `do concurrent`-with-`local(...)` structure of
    !| [[tox_data_integration_jsd_impl(module):determine_all_studies_shared_residual_range_impl(interface)]].
    !|
    !| Not published, and not `_impl` (so it does not itself trigger wrapper generation): it has no
    !| caller yet in production code. It exists so a later restructuring of
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]/
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] (Issue
    !| #187's own Steps 2.5/2.6) can pool a single reference point's residuals across studies before
    !| calling
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| on the result, without duplicating this gather logic when that restructuring lands.
    pure subroutine gather_pooled_neighborhood_residuals(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                         n_neighbors, n_studies, neighborhood_indices_point, &
                                                         pooled_residuals)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN explicitly allowed for missing values
        integer(int32), intent(in) :: neighborhood_indices_point(n_neighbors, n_studies)
            !! Gene indices of one reference point's neighborhood, per study -- one column of a
            !! larger neighborhood_indices_all_studies(n_neighbors, n_points, n_studies), as sliced
            !! by the caller
        real(real64), intent(out) :: pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies)
            !! The pooled residual values for this reference point, across every neighbor and every
            !! study, laid out exactly as a (max_n_reps_all_studies, n_neighbors, n_studies) array
            !! would be

        integer(int32) :: i_study, i_neighbor, i_rep, gene_idx, n_predecessors

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

    !> summary: Generate the GAMMA-decay (n_points, n_neighbors) candidate grid
    !| AUTHOR_LASZLO_LANG
    !| Ported from the grid-building half of 125-stabilize-jscomp's
    !| `determine_js_comp_test_n_points_n_neighbors_helper`: starting from an initial
    !| `n_points_high` (clamped between MIN_POINTS and MAX_POINTS), repeatedly multiplies by
    !| GAMMA until it would drop below `n_points_low`, and for each distinct resulting
    !| `n_points` candidate pairs it with up to `size(KX_FACTORS)` distinct `n_neighbors`
    !| candidates, calling
    !| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]] for
    !| each pair's bin count. A duplicate `n_points` or `n_neighbors` value (from clamping or
    !| integer rounding) collapses rather than repeating -- this is real, derived behavior the
    !| grid depends on to avoid redundant candidates at small `max_n_genes_all_studies`, not a
    !| bug: a small enough `max_n_genes_all_studies` collapses the whole grid down to exactly one
    !| candidate.
    pure subroutine generate_js_comp_test_candidates_impl(max_n_genes_all_studies, residuals, residuals_perm, n_residuals, &
                                                          max_n_reps_all_studies, shared_residual_range, &
                                                          candidates_n_points_n_neighbors, n_bins_candidates, n_candidates)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
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
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! DM_MIN(0.0_real64)
        integer(int32), intent(out) :: candidates_n_points_n_neighbors(2, 16)
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_bins_candidates(16)
            !! Per-candidate bin count from estimate_bin_count_impl, one per candidate pair
            !! DM_RESULT_SIZE_IS(n_candidates)
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate pairs actually filled (at most MAX_CANDIDATE_PAIRS = 16)

        real(real64) :: n_points_high, n_points_low
        integer(int32) :: i_point_candidate, point_candidate, prev_point_candidate
        integer(int32) :: i_neighbor_candidate, neighbor_candidate, prev_neighbor_candidate
        ! Checkpoint 1 of Issue #187 only adds estimate_bin_count_impl's two new diagnostic
        ! outputs; this call site itself is removed in a later checkpoint, so these are genuinely
        ! unused here rather than wired anywhere yet.
        integer(int32) :: sturges_bins_unused, fd_bins_unused

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

                        call estimate_bin_count_impl(residuals, residuals_perm, n_residuals, max_n_reps_all_studies, &
                                                     neighbor_candidate, shared_residual_range, &
                                                     n_bins_candidates(n_candidates), sturges_bins_unused, fd_bins_unused)
                    end if
                end do
            end if
            n_points_high = n_points_high*GAMMA
        end do
    end subroutine generate_js_comp_test_candidates_impl

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
    pure subroutine check_mean_pmf_min_counts_impl(mean_pmf_counts, n_bins, n_points, min_count, all_bins_have_min_count)
        integer(int32), intent(in) :: n_bins
            !! Number of histogram bins
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: mean_pmf_counts(n_bins, n_points)
            !! Absolute counts of a residual per bin for the mean pmf
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! DM_MIN(0_int32)
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` if every bin, at every reference point, reaches at least `min_count`

        integer(int32) :: i_point, i_bin
        logical(c_bool) :: all_pass

        all_pass = .true.
        do concurrent(i_point=1:n_points, i_bin=1:n_bins) shared(mean_pmf_counts, min_count) reduce(.and.:all_pass)
            all_pass = all_pass .and. (mean_pmf_counts(i_bin, i_point) >= min_count)
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

    !> summary: Run the JSD-Comp-Test pipeline for one fixed (n_points, n_neighbors, n_bins) parameter setting
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `js_comp_test_helper`: for every study, builds its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]])
    !| and residual histograms
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]), pools
    !| them into the consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), computes
    !| each study's observed JSD against that consensus
    !| ([[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]/[[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl(interface)]],
    !| called with the consensus pmf as the second argument), runs the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]), and
    !| finally re-derives each study's pmf/JSD/weights/global JSD from its own UNTOUCHED `counts`
    !| via
    !| [[tox_data_integration_jsd_impl(module):calc_pmf_impl(interface)]] -- `mean_pmf`/`mean_pmf_counts`
    !| are NOT re-derived, since they are invariant across permutations by construction (the
    !| permutation test above only resamples its own scratch copies, never `mean_pmf_counts`
    !| itself), exactly as 125 relies on.
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
    !| so this routine gathers each neighbor's actual residual values from `residuals` itself
    !| (`tmp_neighborhood_residuals_gathered`, a per-study scratch buffer) before calling
    !| `build_residual_histograms_impl`. `build_residual_histograms_impl`/`calc_pmf_impl` are
    !| POINT-major (`(n_points, n_bins)`), while `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` here
    !| are BIN-major (`(n_bins, n_points, n_studies)`) to match
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]'s own
    !| convention -- every call across that boundary bridges with an explicit `transpose`, exactly
    !| as [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] and
    !| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] already do.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is folded
    !| into `ierr` (first failure only), matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                                     n_bins, shared_residual_range, gene_means, gene_means_perms, residuals, x_star, &
                                     neighborhood_indices, neighborhood_range, pmfs, counts, included_n_reps, mean_pmf, &
                                     mean_pmf_counts, mean_pmf_included_n_reps, js_divergences, weights, &
                                     global_js_divergence, p_values, tmp_neighborhood_residuals_gathered, &
                                     tmp_counts_point_major, tmp_pmf_point_major, tmp_n_bins_per_point, &
                                     tmp_permutation_mean_pmf_counts, &
                                     tmp_permutation_counts, tmp_permutation_pmfs, tmp_permutation_js_divergences, &
                                     tmp_permutation_weights, tmp_permutation_global_js_divergence, n_permutations, &
                                     random_seed, ierr)
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
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! DM_MIN(0.0_real64)
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
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(int32), dimension(2, n_points, n_studies), intent(out) :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            !! produced by construct_neighborhoods_ranged_impl
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: pmfs
            !! `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`
        integer(int32), dimension(n_bins, n_points, n_studies), intent(out) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The consensus pmf, from create_mean_pmf_impl
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: &
            tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values, reused per study
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts from build_residual_histograms_impl
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf, reused both for
            !! build_residual_histograms_impl's output and for calc_pmf_impl's re-derived pmf
        integer(int32), dimension(n_points), intent(out) :: tmp_n_bins_per_point
            !! Working array: `n_bins` broadcast to every reference point, since
            !! build_residual_histograms_impl now takes a per-point bin count; every point uses the
            !! same `n_bins` here, so this is a pure mechanical translation with no behavior change
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to gjct_permutation_test_impl's own resampling pool
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to gjct_permutation_test_impl's own per-study resampled counts
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to gjct_permutation_test_impl's own resampled pmfs
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
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test

        integer(int32) :: i_study, i_point, i_neighbor, gene_idx, permutation_ierr, actual_n_permutations

        call set_ok(ierr)
        M_DEFAULT_VAL(n_permutations, actual_n_permutations, 1000_int32)

        do i_study = 1, n_studies
            call construct_neighborhoods_ranged_impl(n_points, x_star, max_n_genes_all_studies, gene_means(:, i_study), &
                                                     gene_means_perms(:, i_study), n_neighbors, &
                                                     neighborhood_indices(:, :, i_study), neighborhood_range(:, :, i_study))

            do concurrent(i_point=1:n_points, i_neighbor=1:n_neighbors) local(gene_idx) &
                    shared(tmp_neighborhood_residuals_gathered, residuals, neighborhood_indices, i_study)
                gene_idx = neighborhood_indices(i_neighbor, i_point, i_study)
                tmp_neighborhood_residuals_gathered(:, i_neighbor, i_point) = residuals(:, gene_idx, i_study)
            end do

            tmp_n_bins_per_point(1:n_points) = n_bins
            call build_residual_histograms_impl(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies, n_neighbors, &
                                                n_points, shared_residual_range, n_bins, tmp_n_bins_per_point, &
                                                tmp_counts_point_major, tmp_pmf_point_major, included_n_reps(:, i_study))
            counts(:, :, i_study) = transpose(tmp_counts_point_major)
            pmfs(:, :, i_study) = transpose(tmp_pmf_point_major)
        end do

        call create_mean_pmf_impl(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                                  mean_pmf_included_n_reps, mean_pmf_counts)

        do concurrent(i_study=1:n_studies) shared(pmfs, mean_pmf, n_points, n_bins, js_divergences, included_n_reps, &
                                                  mean_pmf_included_n_reps, global_js_divergence, weights)
            call compute_divergence_per_reference_point_impl(transpose(pmfs(:, :, i_study)), transpose(mean_pmf), n_points, &
                                                             n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
        end do

        call gjct_permutation_test_impl(actual_n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                        mean_pmf_included_n_reps, included_n_reps, global_js_divergence, p_values, &
                                        tmp_permutation_mean_pmf_counts, tmp_permutation_counts, tmp_permutation_pmfs, &
                                        tmp_permutation_js_divergences, tmp_permutation_weights, &
                                        tmp_permutation_global_js_divergence, tmp_pmf_point_major, &
                                        tmp_counts_point_major, random_seed, permutation_ierr)
        if (is_err(permutation_ierr)) call set_err_once(ierr, get_err_code(permutation_ierr))

        ! Re-derive each study's own pmf/JSD/weights/global JSD from its UNTOUCHED counts --
        ! mean_pmf/mean_pmf_counts are NOT re-derived, since the permutation test above only
        ! perturbs its own scratch copies (tmp_permutation_*), never mean_pmf_counts itself.
        do i_study = 1, n_studies
            call calc_pmf_impl(transpose(counts(:, :, i_study)), included_n_reps(:, i_study), n_points, n_bins, &
                               tmp_pmf_point_major)
            pmfs(:, :, i_study) = transpose(tmp_pmf_point_major)
            call compute_divergence_per_reference_point_impl(transpose(pmfs(:, :, i_study)), transpose(mean_pmf), n_points, &
                                                             n_bins, js_divergences(:, i_study))
            call compute_weighted_global_divergence_impl(js_divergences(:, i_study), n_points, included_n_reps(:, i_study), &
                                                          mean_pmf_included_n_reps, global_js_divergence(i_study), &
                                                          weights(:, i_study))
        end do
    end subroutine run_js_comp_test_impl

    !> summary: Search a GAMMA-decay (n_points, n_neighbors) candidate grid for a stable JSD parameter setting
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `determine_js_comp_test_n_points_n_neighbors_helper` and
    !| `_alloc`, merged into one implementation now that the new `_impl` rules leave no separate
    !| hand-written allocation layer. Pools all studies' residuals and gene means, sorts them once
    !| ([[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]), generates the candidate
    !| grid
    !| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]),
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
    !| Ported verbatim, including the
    !| fallback 125 relies on: if no candidate ever plateaus, the search falls back to the FIRST
    !| (finest-resolution) candidate and resets `best_candidate_pair_confidence_interval` to
    !| `-1.0`; if the grid collapsed to a single candidate (see
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]'s
    !| own small-N collapse note), that one candidate is used regardless of whether it plateaued or
    !| even passed either gate -- the plateau machinery is bypassed entirely, exactly as 125 does.
    !|
    !| Per the plan's work-array translation for this routine specifically: `max_n_bins_all_candidates`
    !| (data-dependent, not cheaply closed-form in 125) is replaced by the fixed
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] ceiling, so every
    !| bin-dimensioned work array below is sized to MAX_N_BINS and sliced `(1:n_bins, ...)` per
    !| candidate, rather than carrying a separate recommend-sized dimension argument for it.
    !| `residuals`/`gene_means` are passed to
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]/[[f42_sort_impl(module):sort_real_heapsort_expl_size(interface)]]
    !| as their own multi-dimensional selves -- both callees declare their matching dummy with an
    !| explicit shape, so standard Fortran sequence association reinterprets the contiguous actual
    !| argument as the flat 1-D array they expect, exactly as 125's own `_alloc` layer did for the
    !| same calls.
    !|
    !| Impure: calls the impure
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]. A GSL
    !| failure it reports is folded into `ierr` (first failure only) without aborting the search,
    !| matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_parameter_search_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                       gene_means, residuals, shared_residual_range, n_bootstraps, &
                                                       join_method, max_n_points_candidate, max_n_neighbors_candidate, &
                                                       n_bootstrapping_top_k_jsds, n_points, n_neighbors, n_bins, &
                                                       best_candidate_pair_confidence_interval, plateau_established, &
                                                       n_admissible_evaluated, &
                                                       trace_n_points, trace_n_neighbors, trace_global_js_divergence, &
                                                       trace_ci_lower, trace_ci_upper, trace_ci_width, &
                                                       trace_ci_width_relative, trace_delta, trace_delta_median, &
                                                       trace_delta_max, tmp_gene_means_perms, tmp_gene_means_perm_all, &
                                                       tmp_residuals_perm, tmp_x_star, tmp_neighborhood_indices, &
                                                       tmp_neighborhood_range, tmp_neighborhood_residuals_gathered, &
                                                       tmp_counts_point_major, tmp_pmf_point_major, tmp_n_bins_per_point, &
                                                       tmp_pmfs, tmp_counts, &
                                                       tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                       tmp_mean_pmf_included_n_reps, tmp_js_divergences, tmp_weights, &
                                                       tmp_global_js_divergence, tmp_confidence_interval, &
                                                       tmp_bootstrapping_top_k_jsds, tmp_prev_global_js_divergence, &
                                                       tmp_delta_perm, tmp_best_uncertainty_confidence_interval, &
                                                       min_count_per_mean_bin, min_neighbor_overlap, &
                                                       succeeding_ci_overlap, plateau_mode, delta_median_threshold, &
                                                       delta_max_threshold, delta_epsilon, &
                                                       delta_min_consecutive_transitions, &
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
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! DM_MIN(0.0_real64)
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
            !! Exact upper bound on the grid's first (largest) `n_points` candidate
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_impl, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_candidate_bounds, tox_data_integration_js_comp_test_impl, AUTO)
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
        integer(int32), intent(out) :: n_bins
            !! The finally chosen candidate's bin count
        real(real64), dimension(2, n_studies), intent(out) :: best_candidate_pair_confidence_interval
            !! The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            !! `-1.0_real64` throughout only when `plateau_established` is `.false.` and no
            !! smallest-bootstrap-uncertainty candidate could be substituted either (see
            !! `plateau_established`)
        logical(c_bool), intent(out) :: plateau_established
            !! `.true.` when a real plateau was found (by whichever criterion
            !! `plateau_mode` selected) or the candidate grid never had more than one candidate to
            !! begin with. `.false.` when the search exhausted every admissible candidate
            !! without ever finding one -- Issue #178's own "report that parameter stability could
            !! not be established". When `.false.` and `plateau_mode` is
            !! `MODE_PLATEAU_CI_OVERLAP` and at least one candidate was admissible, the routine
            !! still returns a real (non-`-1.0`) candidate and confidence interval: the admissible
            !! candidate with the smallest bootstrapped uncertainty, per the issue's own fallback
            !! recommendation -- `plateau_established` is what distinguishes that case from an
            !! actual plateau, not the confidence interval's sentinel value
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
            !! candidates_n_points_n_neighbors/n_bins_candidates are in
            !! generate_js_comp_test_candidates_impl
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
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation for the flattened, all-studies-pooled `gene_means`
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: tmp_residuals_perm
            !! Working array: sorting permutation for the flattened, all-studies-pooled `residuals`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array: reference points for the candidate whose `n_points` is current,
            !! recomputed only when `n_points` changes between candidates
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_indices
            !! Working array: one study's neighbor gene indices for the current candidate, reused
            !! per study
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
            !! Working array: the current candidate's `n_bins` broadcast to every reference point,
            !! since build_residual_histograms_impl now takes a per-point bin count; every point
            !! uses the same `n_bins` here, so this is a pure mechanical translation with no
            !! behavior change
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
        integer(int32), intent(in), optional :: min_count_per_mean_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(5_int32)
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

        integer(int32) :: candidates_n_points_n_neighbors(2, MAX_CANDIDATE_PAIRS), n_bins_candidates(MAX_CANDIDATE_PAIRS)
        integer(int32) :: i_candidate, i_study, i_point, i_neighbor, gene_idx, n_pool
        integer(int32) :: prev_n_points, best_candidate_index, best_exceeded_ci_overlap_count, n_candidates
        integer(int32) :: n_residuals, pool_size, bootstrap_ierr, actual_min_count_per_mean_bin
        integer(int32) :: actual_plateau_mode, actual_delta_min_consecutive_transitions, n_consecutive_effect_size_ok
        integer(int32) :: best_uncertainty_candidate_index
        real(real64) :: actual_min_neighbor_overlap, actual_succeeding_ci_overlap
        real(real64) :: actual_delta_median_threshold, actual_delta_max_threshold, actual_delta_epsilon
        real(real64) :: best_uncertainty_value, median_ci_width
        logical(c_bool) :: all_have_min_neighbor_overlap, all_bins_have_min_count, plateau_found
        logical(c_bool) :: ci_plateau_found, effect_size_plateau_found, has_previous_admissible

        call set_ok(ierr)
        M_DEFAULT_VAL(min_count_per_mean_bin, actual_min_count_per_mean_bin, 5_int32)
        M_DEFAULT_VAL(min_neighbor_overlap, actual_min_neighbor_overlap, 0.1_real64)
        M_DEFAULT_VAL(succeeding_ci_overlap, actual_succeeding_ci_overlap, 0.9_real64)
        M_DEFAULT_VAL(plateau_mode, actual_plateau_mode, CM_MODE_PLATEAU_CI_OVERLAP)
        M_DEFAULT_VAL(delta_median_threshold, actual_delta_median_threshold, CM_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_max_threshold, actual_delta_max_threshold, CM_DELTA_MAX_THRESHOLD_DEFAULT)
        M_DEFAULT_VAL(delta_epsilon, actual_delta_epsilon, CM_DELTA_EPSILON_DEFAULT)
        M_DEFAULT_VAL(delta_min_consecutive_transitions, actual_delta_min_consecutive_transitions, CM_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)

        n_residuals = max_n_reps_all_studies*max_n_genes_all_studies*n_studies
        pool_size = max_n_genes_all_studies*n_studies

        ! Sort the pooled residuals once, for estimate_bin_count_impl's percentile calculations
        ! inside generate_js_comp_test_candidates_impl. `residuals` sequence-associates into that
        ! routine's flat explicit-shape dummy (see this routine's own doc comment above).
        call init_perm(tmp_residuals_perm)
        call sort_real_heapsort_expl_size(residuals, tmp_residuals_perm, n_residuals)

        call generate_js_comp_test_candidates_impl(max_n_genes_all_studies, residuals, tmp_residuals_perm, n_residuals, &
                                                   max_n_reps_all_studies, shared_residual_range, &
                                                   candidates_n_points_n_neighbors, n_bins_candidates, n_candidates)

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
            n_bins = n_bins_candidates(i_candidate)

            if (prev_n_points /= n_points) then
                prev_n_points = n_points
                call pool_means_impl(gene_means, tmp_gene_means_perm_all, pool_size, n_points, n_pool, &
                                     tmp_x_star(1:n_points))
            end if

            ! Build every study's neighborhoods for this candidate and check the first
            ! admissibility gate, study by study -- exactly as 125 does, a study that fails the
            ! gate short-circuits the remaining studies for this candidate.
            all_have_min_neighbor_overlap = logical(.true., kind=c_bool)
            do i_study = 1, n_studies
                call construct_neighborhoods_ranged_impl(n_points, tmp_x_star(1:n_points), max_n_genes_all_studies, &
                                                         gene_means(:, i_study), tmp_gene_means_perms(:, i_study), &
                                                         n_neighbors, tmp_neighborhood_indices(1:n_neighbors, 1:n_points), &
                                                         tmp_neighborhood_range(1:2, 1:n_points))

                call check_neighborhood_overlaps_impl(tmp_neighborhood_range(1:2, 1:n_points), n_points, &
                                                       actual_min_neighbor_overlap, all_have_min_neighbor_overlap)
                if (.not. all_have_min_neighbor_overlap) exit

                do concurrent(i_point=1:n_points, i_neighbor=1:n_neighbors) local(gene_idx) &
                        shared(tmp_neighborhood_residuals_gathered, residuals, tmp_neighborhood_indices, i_study)
                    gene_idx = tmp_neighborhood_indices(i_neighbor, i_point)
                    tmp_neighborhood_residuals_gathered(:, i_neighbor, i_point) = residuals(:, gene_idx, i_study)
                end do

                tmp_n_bins_per_point(1:n_points) = n_bins
                call build_residual_histograms_impl(tmp_neighborhood_residuals_gathered(:, 1:n_neighbors, 1:n_points), &
                                                    max_n_reps_all_studies, n_neighbors, n_points, shared_residual_range, &
                                                    n_bins, tmp_n_bins_per_point(1:n_points), &
                                                    tmp_counts_point_major(1:n_points, 1:n_bins), &
                                                    tmp_pmf_point_major(1:n_points, 1:n_bins), &
                                                    tmp_included_n_reps(1:n_points, i_study))
                tmp_counts(1:n_bins, 1:n_points, i_study) = transpose(tmp_counts_point_major(1:n_points, 1:n_bins))
                tmp_pmfs(1:n_bins, 1:n_points, i_study) = transpose(tmp_pmf_point_major(1:n_points, 1:n_bins))
            end do
            if (.not. all_have_min_neighbor_overlap) cycle

            ! Second admissibility gate: every bin of the consensus pmf must reach the minimum count.
            call create_mean_pmf_impl(tmp_pmfs(1:n_bins, 1:n_points, 1:n_studies), tmp_counts(1:n_bins, 1:n_points, 1:n_studies), &
                                      n_bins, n_points, n_studies, tmp_included_n_reps(1:n_points, 1:n_studies), &
                                      tmp_mean_pmf(1:n_bins, 1:n_points), tmp_mean_pmf_included_n_reps(1:n_points), &
                                      tmp_mean_pmf_counts(1:n_bins, 1:n_points))
            call check_mean_pmf_min_counts_impl(tmp_mean_pmf_counts(1:n_bins, 1:n_points), n_bins, n_points, &
                                                actual_min_count_per_mean_bin, all_bins_have_min_count)
            if (.not. all_bins_have_min_count) cycle

            ! Observed JSD per study, seeding the confidence interval bootstrap_histogram_impl bootstraps in place.
            do i_study = 1, n_studies
                call compute_divergence_per_reference_point_impl(transpose(tmp_pmfs(1:n_bins, 1:n_points, i_study)), &
                                                                 transpose(tmp_mean_pmf(1:n_bins, 1:n_points)), n_points, &
                                                                 n_bins, tmp_js_divergences(1:n_points, i_study))
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
            call bootstrap_histogram_impl(n_bootstraps, n_bins, n_points, n_studies, tmp_mean_pmf_counts(1:n_bins, 1:n_points), &
                                          tmp_mean_pmf_included_n_reps(1:n_points), tmp_included_n_reps(1:n_points, 1:n_studies), &
                                          n_bootstrapping_top_k_jsds, tmp_confidence_interval, &
                                          tmp_bootstrapping_top_k_jsds, tmp_counts(1:n_bins, 1:n_points, 1), &
                                          tmp_pmfs(1:n_bins, 1:n_points, 1:n_studies), tmp_mean_pmf(1:n_bins, 1:n_points), &
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
            end if

            call check_plateau_condition_impl(tmp_confidence_interval, best_candidate_pair_confidence_interval, n_studies, &
                                              best_candidate_index, best_exceeded_ci_overlap_count, i_candidate, join_method, &
                                              actual_succeeding_ci_overlap, ci_plateau_found)

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
            end if

            if (plateau_found) exit
        end do

        ! Final candidate pair, one of three cases:
        ! 1. A plateau was reached, or the grid never had more than one candidate to begin with (the
        !    single-candidate/collapsed-grid path bypasses the plateau machinery entirely, exactly
        !    as 125 does): use the best candidate found, plateau_established = .true.
        ! 2. No plateau, plateau_mode is CI-overlap-only, and at least one candidate was ever
        !    admissible: Issue #178's own fallback -- the admissible candidate with the smallest
        !    bootstrapped uncertainty, a real (non--1.0) confidence interval, plateau_established =
        !    .false. Not extended to the effect-size/both modes yet: CI-overlap is the only mode
        !    whose search+fallback combination has been validated end-to-end against a real
        !    known-correct baseline so far.
        ! 3. No plateau, and either plateau_mode isn't CI-overlap-only or zero candidates were ever
        !    admissible (n_candidates counts the raw grid including gate-failed candidates, so this
        !    is distinct from case 2's "at least one admissible" -- nothing exists for case 2 to
        !    select from here): fall back to the finest-resolution (first) candidate and reset the
        !    confidence interval to -1.0, exactly as before this change, plateau_established =
        !    .false.
        if (plateau_found .or. n_candidates < 2_int32) then
            n_points = candidates_n_points_n_neighbors(1, best_candidate_index)
            n_neighbors = candidates_n_points_n_neighbors(2, best_candidate_index)
            n_bins = n_bins_candidates(best_candidate_index)
            plateau_established = logical(.true., kind=c_bool)
        else if (actual_plateau_mode == MODE_PLATEAU_CI_OVERLAP .and. n_admissible_evaluated >= 1_int32) then
            n_points = candidates_n_points_n_neighbors(1, best_uncertainty_candidate_index)
            n_neighbors = candidates_n_points_n_neighbors(2, best_uncertainty_candidate_index)
            n_bins = n_bins_candidates(best_uncertainty_candidate_index)
            best_candidate_pair_confidence_interval = tmp_best_uncertainty_confidence_interval
            plateau_established = logical(.false., kind=c_bool)
        else
            n_points = candidates_n_points_n_neighbors(1, 1)
            n_neighbors = candidates_n_points_n_neighbors(2, 1)
            n_bins = n_bins_candidates(1)
            best_candidate_pair_confidence_interval = -1.0_real64
            plateau_established = logical(.false., kind=c_bool)
        end if
    end subroutine run_js_comp_test_parameter_search_impl

end module tox_data_integration_js_comp_test_impl
