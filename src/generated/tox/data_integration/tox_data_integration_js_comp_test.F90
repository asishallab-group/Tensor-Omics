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
!| `calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
!| allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| Later stages of the port add the K-study permutation test and the top-level orchestrator that
!| wire these building blocks together.
!|
!| Generated from [[tox_data_integration_js_comp_test_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test
    use f42_safeguard
    use tox_data_integration_js_comp_test_impl, only: MAX_N_BINS, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN
    use tox_data_integration_js_comp_test_impl, only: MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, bootstrap_histogram_impl
    use tox_data_integration_js_comp_test_impl, only: calc_js_comp_test_n_top_k_jsds, check_effect_size_plateau_condition_impl, check_mean_pmf_min_counts_impl, check_neighborhood_overlaps_impl
    use tox_data_integration_js_comp_test_impl, only: check_plateau_condition_impl, create_mean_pmf_impl, create_mean_pmf_only_impl, determine_bin_count_occupancy_impl
    use tox_data_integration_js_comp_test_impl, only: estimate_bin_count_impl, generate_js_comp_test_candidates_impl, run_js_comp_test_impl, run_js_comp_test_parameter_search_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_int
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: estimate_bin_count
    public :: estimate_bin_count_expert
    public :: determine_bin_count_occupancy
    public :: determine_bin_count_occupancy_expert
    public :: generate_js_comp_test_candidates
    public :: check_neighborhood_overlaps
    public :: check_mean_pmf_min_counts
    public :: check_plateau_condition
    public :: check_effect_size_plateau_condition
    public :: check_effect_size_plateau_condition_expert
    public :: create_mean_pmf
    public :: create_mean_pmf_only
    public :: bootstrap_histogram
    public :: bootstrap_histogram_expert
    public :: run_js_comp_test
    public :: run_js_comp_test_expert
    public :: run_js_comp_test_parameter_search
    public :: run_js_comp_test_parameter_search_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):estimate_bin_count_expert]] to prepare it yourself.
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    !| the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    !| is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    !| range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    !| independently, each clamped on its own to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins and returned as
    !| `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    !| interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    !| skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
    pure subroutine estimate_bin_count(&
            residuals,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            shared_residual_range,&
            n_bins,&
            sturges_bins,&
            fd_bins,&
            ierr&
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), dimension(n_residuals), intent(in) :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(int32), intent(out) :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            !! to the (clamped) sturges_bins when the interquartile range is too close to zero to
            !! divide by (see the is_close guard below)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: residuals_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=2_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_all_in_range_real(residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(residuals_perm(n_residuals))
        call init_perm(residuals_perm)
        call sort_array_heapsort(residuals, residuals_perm)

        call estimate_bin_count_impl(&
            residuals = residuals,&
            residuals_perm = residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins&
        )
    end subroutine estimate_bin_count

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):estimate_bin_count]] does both.
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    !| the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    !| is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    !| range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    !| independently, each clamped on its own to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins and returned as
    !| `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    !| interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    !| skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
    pure subroutine estimate_bin_count_expert(&
            residuals,&
            residuals_perm,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            shared_residual_range,&
            n_bins,&
            sturges_bins,&
            fd_bins,&
            ierr&
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), dimension(n_residuals), intent(in) :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(int32), dimension(n_residuals), intent(in) :: residuals_perm
            !! Sorting permutation for `residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(int32), intent(out) :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        integer(int32), intent(out) :: sturges_bins
            !! Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        integer(int32), intent(out) :: fd_bins
            !! Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            !! to the (clamped) sturges_bins when the interquartile range is too close to zero to
            !! divide by (see the is_close guard below)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=3_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=6_int32, min=0.0_real64)
        call validate_all_in_range_real(residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(residuals_perm, n_residuals, ierr, arg_pos=2_int32, min=1_int32, max=n_residuals)
        if (is_err(ierr)) return
#endif

        call estimate_bin_count_impl(&
            residuals = residuals,&
            residuals_perm = residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins&
        )
    end subroutine estimate_bin_count_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_expert]] to prepare it yourself.
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
    pure subroutine determine_bin_count_occupancy(&
            pooled_residuals,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            selected_n_bins,&
            occupancy_failed,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            m_min,&
            m_max,&
            min_residuals_per_bin,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            ierr&
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), dimension(n_residuals), intent(in) :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
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
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); if a caller passes
            !! `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
            !! relying on an unconfirmed generator capability to bound one optional argument by
            !! another
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=2_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=15_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=16_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=17_int32, min=0_int32)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=18_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=19_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=20_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pooled_residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(pooled_residuals_perm(n_residuals))
        M_ALLOCATE(tmp_bin_counts(256))
        call init_perm(pooled_residuals_perm)
        call sort_array_heapsort(pooled_residuals, pooled_residuals_perm)

        call determine_bin_count_occupancy_impl(&
            pooled_residuals = pooled_residuals,&
            pooled_residuals_perm = pooled_residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            selected_n_bins = selected_n_bins,&
            occupancy_failed = occupancy_failed,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            tmp_bin_counts = tmp_bin_counts,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile&
        )
    end subroutine determine_bin_count_occupancy

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy]] does both.
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
    pure subroutine determine_bin_count_occupancy_expert(&
            pooled_residuals,&
            pooled_residuals_perm,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            selected_n_bins,&
            occupancy_failed,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            tmp_bin_counts,&
            m_min,&
            m_max,&
            min_residuals_per_bin,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            ierr&
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), dimension(n_residuals), intent(in) :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(int32), dimension(n_residuals), intent(in) :: pooled_residuals_perm
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
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
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts
            !! Working array: per-bin counts of whichever candidate bin count is currently being
            !! tested, reused throughout the search (256 = MAX_N_BINS, written as a literal since a
            !! generated wrapper's dummy dimension cannot reference a module parameter)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); if a caller passes
            !! `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
            !! relying on an unconfirmed generator capability to bound one optional argument by
            !! another
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=3_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=17_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=18_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=19_int32, min=0_int32)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=20_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=21_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=22_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pooled_residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_residuals_perm, n_residuals, ierr, arg_pos=2_int32, min=1_int32, max=n_residuals)
        if (is_err(ierr)) return
#endif

        call determine_bin_count_occupancy_impl(&
            pooled_residuals = pooled_residuals,&
            pooled_residuals_perm = pooled_residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            selected_n_bins = selected_n_bins,&
            occupancy_failed = occupancy_failed,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            tmp_bin_counts = tmp_bin_counts,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile&
        )
    end subroutine determine_bin_count_occupancy_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl]].
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
    pure subroutine generate_js_comp_test_candidates(&
            max_n_genes_all_studies,&
            candidates_n_points_n_neighbors,&
            n_candidates,&
            ierr&
        )
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(2, 16), intent(out) :: candidates_n_points_n_neighbors
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! The first `n_candidates` elements will hold the results.
        integer(int32), intent(out) :: n_candidates
            !! Number of candidate pairs actually filled (at most MAX_CANDIDATE_PAIRS = 16)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        if (is_err(ierr)) return
#endif

        call generate_js_comp_test_candidates_impl(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            candidates_n_points_n_neighbors = candidates_n_points_n_neighbors,&
            n_candidates = n_candidates&
        )
    end subroutine generate_js_comp_test_candidates

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl]].
    !| Ported from 125-stabilize-jscomp's `test_neighborhood_overlaps_helper`: the first
    !| admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, using the
    !| `[min_idx, max_idx]` neighborhood spans
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| produces. Named `check_*` rather than 125's `test_*`, a deliberate deviation from the
    !| verbatim port: the generated R binding is published under the Fortran name, and the R test
    !| harness (`r/test_helpers.R`'s `run_all_tests`) discovers every `test_`-prefixed name in the
    !| environment as a test case to run with no arguments -- a `test_`-prefixed export would be
    !| swept up and fail every R suite that sources the package, not just this module's own.
    pure subroutine check_neighborhood_overlaps(&
            neighborhood_range,&
            n_points,&
            min_neighbor_overlap,&
            all_have_min_neighbor_overlap,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), dimension(2, n_points), intent(in) :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span, as produced
            !! by construct_neighborhoods_ranged_impl
            !! The minimum valid value is `1_int32`.
        real(real64), intent(in) :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=3_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(neighborhood_range, 2 * n_points, ierr, arg_pos=1_int32, min=1_int32)
        if (is_err(ierr)) return
#endif

        call check_neighborhood_overlaps_impl(&
            neighborhood_range = neighborhood_range,&
            n_points = n_points,&
            min_neighbor_overlap = min_neighbor_overlap,&
            all_have_min_neighbor_overlap = all_have_min_neighbor_overlap&
        )
    end subroutine check_neighborhood_overlaps

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl]].
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
    pure subroutine check_mean_pmf_min_counts(&
            mean_pmf_counts,&
            n_bins,&
            n_bins_per_point,&
            n_points,&
            min_count,&
            all_bins_have_min_count,&
            ierr&
        )
        integer(int32), intent(in) :: n_bins
            !! Second extent of `mean_pmf_counts` -- the widest bin count any point uses (max_n_bins),
            !! not necessarily every point's own bin count
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: n_bins_per_point
            !! This point's own bin count -- only `mean_pmf_counts(1:n_bins_per_point(i_point), i_point)`
            !! is inspected; columns beyond it are legitimate zero-padding, not failures
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_bins`.
        integer(int32), intent(in) :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! The minimum valid value is `0_int32`.
        logical(c_bool), intent(out) :: all_bins_have_min_count
            !! `.true.` if every bin within each reference point's own `n_bins_per_point`, at every
            !! reference point, reaches at least `min_count`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_bins, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_int(min_count, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_counts, n_bins * n_points, ierr, arg_pos=1_int32, min=0_int32)
        call validate_all_in_range_int(n_bins_per_point, n_points, ierr, arg_pos=3_int32, min=1_int32, max=n_bins)
        if (is_err(ierr)) return
#endif

        call check_mean_pmf_min_counts_impl(&
            mean_pmf_counts = mean_pmf_counts,&
            n_bins = n_bins,&
            n_bins_per_point = n_bins_per_point,&
            n_points = n_points,&
            min_count = min_count,&
            all_bins_have_min_count = all_bins_have_min_count&
        )
    end subroutine check_mean_pmf_min_counts

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl]].
    !| Ported from 125-stabilize-jscomp's `check_plateau_condition_helper`. A search over
    !| candidates (finest resolution to coarsest) stops -- "plateaus" -- either when a new
    !| candidate is no better than the previous best (the short-circuit below: keep the previous
    !| best and stop searching), or once the new candidate's confidence-interval overlap with the
    !| previous best meets the condition `join_method` names. `join_method` replaces
    !| 125-stabilize-jscomp's hand-rolled join-method-range validation macro entirely: the mode
    !| table below is itself the validation, checked against exactly the values it names.
    pure subroutine check_plateau_condition(&
            confidence_interval,&
            best_candidate_pair_confidence_interval,&
            n_studies,&
            best_candidate_index,&
            best_exceeded_ci_overlap_count,&
            candidate_index,&
            join_method,&
            succeeding_ci_overlap,&
            plateau_found,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(2, n_studies), intent(in) :: confidence_interval
            !! JSD confidence interval `[lower, upper]` from bootstrapping, for the candidate
            !! pair under test
        real(real64), dimension(2, n_studies), intent(inout) :: best_candidate_pair_confidence_interval
            !! JSD confidence intervals for the current best candidate pair; overwritten with
            !! `confidence_interval` unless the new candidate is worse
        integer(int32), intent(inout) :: best_candidate_index
            !! Candidate-grid index of the current best candidate pair; overwritten with
            !! `candidate_index` unless the new candidate is worse
        integer(int32), intent(inout) :: best_exceeded_ci_overlap_count
            !! Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            !! best candidate pair; overwritten unless the new candidate is worse
            !! The minimum valid value is `0_int32`.
        integer(int32), intent(in) :: candidate_index
            !! Candidate-grid index of the candidate pair that produced `confidence_interval`
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition: METHOD_JOIN_MIN requires every study's overlap to exceed
            !! `succeeding_ci_overlap`, METHOD_JOIN_MAX requires only one study's overlap to
            !! exceed it, and METHOD_JOIN_MEDIAN requires a majority
            !! (`count > (n_studies - 1) / 2`) to exceed it
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum fractional overlap an interval in `confidence_interval` must have with its
            !! respective interval in `best_candidate_pair_confidence_interval` to count as
            !! "exceeded"
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_in_range_int(best_exceeded_ci_overlap_count, ierr, arg_pos=5_int32, min=0_int32)
        call validate_in_range_int(candidate_index, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(best_candidate_pair_confidence_interval, 2 * n_studies, ierr, arg_pos=2_int32)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return
#endif

        call check_plateau_condition_impl(&
            confidence_interval = confidence_interval,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            n_studies = n_studies,&
            best_candidate_index = best_candidate_index,&
            best_exceeded_ci_overlap_count = best_exceeded_ci_overlap_count,&
            candidate_index = candidate_index,&
            join_method = join_method,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_found = plateau_found&
        )
    end subroutine check_plateau_condition

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition_expert]] to prepare it yourself.
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
    pure subroutine check_effect_size_plateau_condition(&
            global_js_divergence,&
            prev_global_js_divergence,&
            n_studies,&
            has_previous,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            n_consecutive_ok,&
            delta,&
            delta_median,&
            delta_max,&
            plateau_found,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_studies), intent(in) :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_studies), intent(in) :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in) :: has_previous
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_studies), intent(out) :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(real64), intent(out) :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(real64), intent(out) :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_delta_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=5_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=6_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=7_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_consecutive_ok, ierr, arg_pos=9_int32, min=0_int32)
        call validate_all_in_range_real(global_js_divergence, n_studies, ierr, arg_pos=1_int32, min=0.0_real64)
        call validate_all_in_range_real(prev_global_js_divergence, n_studies, ierr, arg_pos=2_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_delta_perm(n_studies))

        call check_effect_size_plateau_condition_impl(&
            global_js_divergence = global_js_divergence,&
            prev_global_js_divergence = prev_global_js_divergence,&
            n_studies = n_studies,&
            has_previous = has_previous,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            n_consecutive_ok = n_consecutive_ok,&
            delta = delta,&
            delta_median = delta_median,&
            delta_max = delta_max,&
            plateau_found = plateau_found,&
            tmp_delta_perm = tmp_delta_perm&
        )
    end subroutine check_effect_size_plateau_condition

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition]] does both.
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
    pure subroutine check_effect_size_plateau_condition_expert(&
            global_js_divergence,&
            prev_global_js_divergence,&
            n_studies,&
            has_previous,&
            delta_median_threshold,&
            delta_max_threshold,&
            delta_epsilon,&
            delta_min_consecutive_transitions,&
            n_consecutive_ok,&
            delta,&
            delta_median,&
            delta_max,&
            plateau_found,&
            tmp_delta_perm,&
            ierr&
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_studies), intent(in) :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_studies), intent(in) :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in) :: has_previous
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_studies), intent(out) :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(real64), intent(out) :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(real64), intent(out) :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out) :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array: sorting permutation for `delta`, used to compute `delta_median`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_studies, ierr, arg_pos=3_int32)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=5_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=6_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=7_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(n_consecutive_ok, ierr, arg_pos=9_int32, min=0_int32)
        call validate_all_in_range_real(global_js_divergence, n_studies, ierr, arg_pos=1_int32, min=0.0_real64)
        call validate_all_in_range_real(prev_global_js_divergence, n_studies, ierr, arg_pos=2_int32, min=0.0_real64)
        if (is_err(ierr)) return
#endif

        call check_effect_size_plateau_condition_impl(&
            global_js_divergence = global_js_divergence,&
            prev_global_js_divergence = prev_global_js_divergence,&
            n_studies = n_studies,&
            has_previous = has_previous,&
            delta_median_threshold = delta_median_threshold,&
            delta_max_threshold = delta_max_threshold,&
            delta_epsilon = delta_epsilon,&
            delta_min_consecutive_transitions = delta_min_consecutive_transitions,&
            n_consecutive_ok = n_consecutive_ok,&
            delta = delta,&
            delta_median = delta_median,&
            delta_max = delta_max,&
            plateau_found = plateau_found,&
            tmp_delta_perm = tmp_delta_perm&
        )
    end subroutine check_effect_size_plateau_condition_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl]].
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_helper`.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    pure subroutine create_mean_pmf(&
            pmfs,&
            counts,&
            n_bins,&
            n_points,&
            n_studies,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_included_n_reps,&
            mean_pmf_counts,&
            ierr&
        )
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
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_bins, n_points, n_studies), intent(in) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for `mean_pmf`,
            !! summed across all n_studies
        integer(int32), dimension(n_bins, n_points), intent(out) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts, dim=3)`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_bins, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(pmfs, n_bins * n_points * n_studies, ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(counts, n_bins * n_points * n_studies, ierr, arg_pos=2_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=6_int32, min=0_int32)
        if (is_err(ierr)) return
#endif

        call create_mean_pmf_impl(&
            pmfs = pmfs,&
            counts = counts,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            mean_pmf_counts = mean_pmf_counts&
        )
    end subroutine create_mean_pmf

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_only_impl]].
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_only_helper`: useful where the
    !| mean pmf's own counts don't matter, e.g. for the bootstrap confidence interval a later
    !| stage of this port adds.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    pure subroutine create_mean_pmf_only(&
            pmfs,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf,&
            ierr&
        )
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
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_bins, n_points), intent(out) :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_bins, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(pmfs, n_bins * n_points * n_studies, ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call create_mean_pmf_only_impl(&
            pmfs = pmfs,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf = mean_pmf&
        )
    end subroutine create_mean_pmf_only

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):bootstrap_histogram_expert]] to prepare it yourself.
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
    subroutine bootstrap_histogram(&
            n_bootstraps,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            confidence_interval,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            ierr&
        )
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
            !! otherwise used here
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator
        integer(int32) :: n_bootstrapping_top_k_jsds
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=9_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(mean_pmf_counts, n_bins * n_points, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=6_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=7_int32, min=0_int32)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_n_top_k_jsds(&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            n_top_k = n_bootstrapping_top_k_jsds&
        )
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_counts(n_bins, n_points))
        M_ALLOCATE(tmp_pmfs(n_bins, n_points, n_studies))
        M_ALLOCATE(tmp_mean_pmf(n_bins, n_points))
        M_ALLOCATE(tmp_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_weights(n_points, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))

        call bootstrap_histogram_impl(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            confidence_interval = confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_counts = tmp_counts,&
            tmp_pmfs = tmp_pmfs,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine bootstrap_histogram

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):bootstrap_histogram]] does both.
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
    subroutine bootstrap_histogram_expert(&
            n_bootstraps,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            n_bootstrapping_top_k_jsds,&
            confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_counts,&
            tmp_pmfs,&
            tmp_mean_pmf,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            ierr&
        )
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
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
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=17_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(mean_pmf_counts, n_bins * n_points, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=6_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=7_int32, min=0_int32)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=9_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call bootstrap_histogram_impl(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            confidence_interval = confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_counts = tmp_counts,&
            tmp_pmfs = tmp_pmfs,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine bootstrap_histogram_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):run_js_comp_test_expert]] to prepare it yourself.
    !| Ported from 125-stabilize-jscomp's `js_comp_test_helper`, restructured for Issue #187's
    !| occupancy-constrained per-neighborhood histogram binning into three passes, mirroring
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]'s
    !| own Pass A/B/C split (Issue #187's own Steps 2.5/2.6):
    !|
    !| - Pass A (per study): builds every study's neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]),
    !| writing into `neighborhood_indices`/`neighborhood_range`, which already retain every
    !| study's own values simultaneously (both are real `intent(out)` arguments sized
    !| `(..., n_points, n_studies)` -- unlike `run_js_comp_test_parameter_search_impl`, no new
    !| buffer was needed for this). Unlike that routine, there is no admissibility gate here, so
    !| Pass A always runs to completion for every study.
    !| - Pass B (per point, sequential -- see the implementation body's own comment for why): pools
    !| every study's residuals for one reference point at a time (`gather_pooled_neighborhood_residuals`,
    !| a private module helper, not itself published) and runs Issue #187's occupancy-constrained
    !| bin-count search on the pooled result
    !| ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]),
    !| deciding `n_bins_per_point(i_point)` independently for every reference point, plus the
    !| `occupancy_failed`/`n_pooled_residuals`/`min_bin_occupancy`/`mean_bin_occupancy`/
    !| `max_bin_occupancy`/`sturges_bins`/`fd_bins` diagnostics. `max_n_bins_per_point`
    !| (`maxval(n_bins_per_point(1:n_points))`) is derived once after Pass B and replaces the old
    !| caller-supplied scalar `n_bins` everywhere downstream.
    !| - Pass C (per study): re-gathers this study's residual values from the neighbor indices Pass
    !| A already computed, then builds its residual histograms at the real per-point bin counts
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !|
    !| After Pass C, the pipeline continues exactly as before: pools the per-study pmfs into the
    !| consensus pmf
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
    !| (`tmp_neighborhood_residuals_gathered`, a per-study scratch buffer) before calling
    !| `build_residual_histograms_impl`. `build_residual_histograms_impl`/`calc_pmf_impl` are
    !| POINT-major (`(n_points, max_n_bins_per_point)`), while `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts`
    !| here are BIN-major (`(256, n_points, n_studies)`) to match
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]'s own
    !| convention -- every call across that boundary bridges with an explicit `transpose`, exactly
    !| as [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] and
    !| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] already do.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is folded
    !| into `ierr` (first failure only), matching that routine's own tolerant precedent.
    subroutine run_js_comp_test(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            n_neighbors,&
            gene_means,&
            gene_means_perms,&
            residuals,&
            x_star,&
            neighborhood_indices,&
            neighborhood_range,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins_per_point,&
            occupancy_failed,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            p_values,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
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
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
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
            !! leading, meaningful bins/rows in `pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`,
            !! `tmp_counts_point_major` and `tmp_pmf_point_major` below -- those are all declared
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
            !! Empirical p-value per study from gjct_permutation_test_impl
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible in Pass B's occupancy search, forwarded to
            !! determine_bin_count_occupancy_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        integer(int32), dimension(:, :), allocatable :: tmp_counts_point_major
        real(real64), dimension(:, :), allocatable :: tmp_pmf_point_major
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_permutation_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_permutation_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_permutation_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_permutation_weights
        real(real64), dimension(:), allocatable :: tmp_permutation_global_js_divergence

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=33_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=35_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=36_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=37_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=38_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=39_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=40_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(gene_means_perms, max_n_genes_all_studies * n_studies, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=8_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=9_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points))
        M_ALLOCATE(tmp_counts_point_major(n_points, 256))
        M_ALLOCATE(tmp_pmf_point_major(n_points, 256))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(256, n_points))
        M_ALLOCATE(tmp_permutation_counts(256, n_points))
        M_ALLOCATE(tmp_permutation_pmfs(256, n_points, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_weights(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_global_js_divergence(n_studies))

        call run_js_comp_test_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            gene_means = gene_means,&
            gene_means_perms = gene_means_perms,&
            residuals = residuals,&
            x_star = x_star,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins_per_point = max_n_bins_per_point,&
            occupancy_failed = occupancy_failed,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            p_values = p_values,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major = tmp_counts_point_major,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):run_js_comp_test]] does both.
    !| Ported from 125-stabilize-jscomp's `js_comp_test_helper`, restructured for Issue #187's
    !| occupancy-constrained per-neighborhood histogram binning into three passes, mirroring
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]'s
    !| own Pass A/B/C split (Issue #187's own Steps 2.5/2.6):
    !|
    !| - Pass A (per study): builds every study's neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]),
    !| writing into `neighborhood_indices`/`neighborhood_range`, which already retain every
    !| study's own values simultaneously (both are real `intent(out)` arguments sized
    !| `(..., n_points, n_studies)` -- unlike `run_js_comp_test_parameter_search_impl`, no new
    !| buffer was needed for this). Unlike that routine, there is no admissibility gate here, so
    !| Pass A always runs to completion for every study.
    !| - Pass B (per point, sequential -- see the implementation body's own comment for why): pools
    !| every study's residuals for one reference point at a time (`gather_pooled_neighborhood_residuals`,
    !| a private module helper, not itself published) and runs Issue #187's occupancy-constrained
    !| bin-count search on the pooled result
    !| ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]),
    !| deciding `n_bins_per_point(i_point)` independently for every reference point, plus the
    !| `occupancy_failed`/`n_pooled_residuals`/`min_bin_occupancy`/`mean_bin_occupancy`/
    !| `max_bin_occupancy`/`sturges_bins`/`fd_bins` diagnostics. `max_n_bins_per_point`
    !| (`maxval(n_bins_per_point(1:n_points))`) is derived once after Pass B and replaces the old
    !| caller-supplied scalar `n_bins` everywhere downstream.
    !| - Pass C (per study): re-gathers this study's residual values from the neighbor indices Pass
    !| A already computed, then builds its residual histograms at the real per-point bin counts
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !|
    !| After Pass C, the pipeline continues exactly as before: pools the per-study pmfs into the
    !| consensus pmf
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
    !| (`tmp_neighborhood_residuals_gathered`, a per-study scratch buffer) before calling
    !| `build_residual_histograms_impl`. `build_residual_histograms_impl`/`calc_pmf_impl` are
    !| POINT-major (`(n_points, max_n_bins_per_point)`), while `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts`
    !| here are BIN-major (`(256, n_points, n_studies)`) to match
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]'s own
    !| convention -- every call across that boundary bridges with an explicit `transpose`, exactly
    !| as [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] and
    !| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] already do.
    !|
    !| Impure: calls the impure `gjct_permutation_test_impl`. A GSL failure it reports is folded
    !| into `ierr` (first failure only), matching that routine's own tolerant precedent.
    subroutine run_js_comp_test_expert(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            n_neighbors,&
            gene_means,&
            gene_means_perms,&
            residuals,&
            x_star,&
            neighborhood_indices,&
            neighborhood_range,&
            n_bins_per_point,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins_per_point,&
            occupancy_failed,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            pmfs,&
            counts,&
            included_n_reps,&
            mean_pmf,&
            mean_pmf_counts,&
            mean_pmf_included_n_reps,&
            js_divergences,&
            weights,&
            global_js_divergence,&
            p_values,&
            tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major,&
            tmp_pmf_point_major,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts,&
            tmp_permutation_pmfs,&
            tmp_permutation_js_divergences,&
            tmp_permutation_weights,&
            tmp_permutation_global_js_divergence,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
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
        integer(int32), intent(in) :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
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
            !! leading, meaningful bins/rows in `pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`,
            !! `tmp_counts_point_major` and `tmp_pmf_point_major` below -- those are all declared
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
            !! Empirical p-value per study from gjct_permutation_test_impl
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values, reused per study
            !! (Pass C)
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts from
            !! build_residual_histograms_impl. `256` = MAX_N_BINS, written as a literal because a
            !! generated wrapper's dummy dimension cannot reference a module parameter
        real(real64), dimension(n_points, 256), intent(out) :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf, reused both for
            !! build_residual_histograms_impl's output and for calc_pmf_impl's re-derived pmf.
            !! `256` = MAX_N_BINS, see tmp_counts_point_major above
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
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible in Pass B's occupancy search, forwarded to
            !! determine_bin_count_occupancy_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=45_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=47_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=48_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=49_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=50_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=51_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=52_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(gene_means_perms, max_n_genes_all_studies * n_studies, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=8_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=9_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            gene_means = gene_means,&
            gene_means_perms = gene_means_perms,&
            residuals = residuals,&
            x_star = x_star,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            n_bins_per_point = n_bins_per_point,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins_per_point = max_n_bins_per_point,&
            occupancy_failed = occupancy_failed,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            pmfs = pmfs,&
            counts = counts,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            js_divergences = js_divergences,&
            weights = weights,&
            global_js_divergence = global_js_divergence,&
            p_values = p_values,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_counts_point_major = tmp_counts_point_major,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_permutation_mean_pmf_counts = tmp_permutation_mean_pmf_counts,&
            tmp_permutation_counts = tmp_permutation_counts,&
            tmp_permutation_pmfs = tmp_permutation_pmfs,&
            tmp_permutation_js_divergences = tmp_permutation_js_divergences,&
            tmp_permutation_weights = tmp_permutation_weights,&
            tmp_permutation_global_js_divergence = tmp_permutation_global_js_divergence,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search_expert]] to prepare it yourself.
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
    subroutine run_js_comp_test_parameter_search(&
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
            !! Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
            !! per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
            !! jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
            !! is no longer purely an internal sizing detail the plain wrapper can compute and
            !! hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
            !! rather than AUTO here now
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
            !! because anything returned is sized by it (nothing is), but because it comes from the
            !! same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
            !! now that that call can no longer run automatically inside this wrapper, splitting
            !! this one back into an AUTO call would just be a second, redundant call to the same
            !! routine for no benefit
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
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
            !! never ran for the returned candidate (see the final three-way branch's own comments)
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The finally chosen candidate's per-point upper residual-range bound (R_high),
            !! mirroring `shared_residual_range_low` above in every respect
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
            !! candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
            !! formula omits this floor, but the same near-zero-JSD instability that motivates
            !! `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
            !! destabilizes any ratio that divides by it, whichever candidate's `J` it is)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
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
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
            !! from determine_bin_count_occupancy_impl (never part of the occupancy search's own
            !! decision). Jagged per candidate column exactly as trace_selected_n_bins above --
            !! only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
            !! must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
            !! diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
            !! search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
            !! above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
            !! caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate. Reuses Issue #187's occupancy-search default rather than an
            !! independently-tunable threshold of its own: once
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! wires real per-neighborhood bin counts in, a separate laxer threshold here would
            !! silently let a candidate the occupancy search already marked `occupancy_failed`
            !! pass this gate anyway, defeating the FAILURE-detection mechanism
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                      | Value                                                                                 |
            !! |-------------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only       | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                          | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:, :), allocatable :: tmp_gene_means_perms
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :, :), allocatable :: tmp_neighborhood_indices_all_studies
        integer(int32), dimension(:, :), allocatable :: tmp_neighborhood_range
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        integer(int32), dimension(:, :), allocatable :: tmp_counts_point_major
        real(real64), dimension(:, :), allocatable :: tmp_pmf_point_major
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
        real(real64), dimension(:, :), allocatable :: tmp_confidence_interval
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:), allocatable :: tmp_prev_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_delta_perm
        real(real64), dimension(:, :), allocatable :: tmp_best_uncertainty_confidence_interval
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
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
        call validate_in_range_int(max_n_neighbors_candidate, ierr, arg_pos=9_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=38_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=39_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=40_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=42_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=43_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=44_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=45_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=46_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=47_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=48_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=49_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=50_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=51_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=41_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_n_top_k_jsds(&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            n_top_k = n_bootstrapping_top_k_jsds&
        )
        M_ALLOCATE(tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_x_star(max_n_points_candidate))
        M_ALLOCATE(tmp_neighborhood_indices_all_studies(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_neighborhood_range(2, max_n_points_candidate))
        M_ALLOCATE(tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
        M_ALLOCATE(tmp_counts_point_major(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_pmf_point_major(max_n_points_candidate, 256))
        M_ALLOCATE(tmp_n_bins_per_point(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_low(max_n_points_candidate))
        M_ALLOCATE(tmp_shared_residual_range_high(max_n_points_candidate))
        M_ALLOCATE(tmp_pmfs(256, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_counts(256, max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_included_n_reps(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_mean_pmf(256, max_n_points_candidate))
        M_ALLOCATE(tmp_mean_pmf_counts(256, max_n_points_candidate))
        M_ALLOCATE(tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        M_ALLOCATE(tmp_js_divergences(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_weights(max_n_points_candidate, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_delta_perm(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
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

        call run_js_comp_test_parameter_search_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
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
            tmp_counts_point_major = tmp_counts_point_major,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
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
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_parameter_search

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search]] does both.
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
    subroutine run_js_comp_test_parameter_search_expert(&
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
            tmp_counts_point_major,&
            tmp_pmf_point_major,&
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
            !! Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
            !! per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
            !! jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
            !! is no longer purely an internal sizing detail the plain wrapper can compute and
            !! hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
            !! rather than AUTO here now
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
            !! because anything returned is sized by it (nothing is), but because it comes from the
            !! same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
            !! now that that call can no longer run automatically inside this wrapper, splitting
            !! this one back into an AUTO call would just be a second, redundant call to the same
            !! routine for no benefit
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstraps to perform for a candidate pair
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: join_method
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
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
            !! never ran for the returned candidate (see the final three-way branch's own comments)
        real(real64), dimension(max_n_points_candidate), intent(out) :: shared_residual_range_high
            !! The finally chosen candidate's per-point upper residual-range bound (R_high),
            !! mirroring `shared_residual_range_low` above in every respect
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
            !! candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! (`L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! (`U_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width (`W_{i,t} = U_{i,t} -
            !! L_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
            !! formula omits this floor, but the same near-zero-JSD instability that motivates
            !! `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
            !! destabilizes any ratio that divides by it, whichever candidate's `J` it is)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
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
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
            !! from determine_bin_count_occupancy_impl (never part of the occupancy search's own
            !! decision). Jagged per candidate column exactly as trace_selected_n_bins above --
            !! only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
            !! must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
            !! diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
            !! search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
            !! above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
            !! caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation for the flattened, all-studies-pooled `gene_means`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array: reference points for the candidate whose `n_points` is current,
            !! recomputed only when `n_points` changes between candidates
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate,
            !! retained simultaneously (unlike the single-study-reused buffer this replaces) so
            !! Pass B below can pool residuals across studies for one reference point at a time,
            !! and Pass C can re-gather each study's own residuals from the already-known indices
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_range
            !! Working array: one study's `[min_idx, max_idx]` neighborhood spans for the current
            !! candidate, reused per study
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_residuals_gathered
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
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals across every neighbor and
            !! every study (Pass B), reused per point -- one small buffer, not one per point, since
            !! Pass B is a deliberate sequential loop (see the module-internal doc comment on the
            !! implementation body)
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals_perm
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
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        integer(int32), intent(in), optional :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                      | Value                                                                                 |
            !! |-------------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only       | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                          | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `0_int32`.
        real(real64), intent(in), optional :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(real64), intent(in), optional :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_int(max_n_neighbors_candidate, ierr, arg_pos=9_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=10_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=80_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=81_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=82_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=84_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=85_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=86_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=87_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=88_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=89_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=90_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=91_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=92_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=93_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=83_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_parameter_search_impl(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method,&
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
            tmp_counts_point_major = tmp_counts_point_major,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
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
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_parameter_search_expert

end module tox_data_integration_js_comp_test
