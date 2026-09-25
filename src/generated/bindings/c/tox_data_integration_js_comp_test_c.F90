#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
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
!| allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| Later stages of the port add the K-study permutation test and the top-level orchestrator that
!| wire these building blocks together.
module tox_data_integration_js_comp_test_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: estimate_bin_count_c
    public :: estimate_bin_count_expert_c
    public :: determine_bin_count_occupancy_c
    public :: determine_bin_count_occupancy_expert_c
    public :: determine_bin_count_occupancy_exhaustive_c
    public :: determine_bin_count_occupancy_exhaustive_expert_c
    public :: generate_js_comp_test_candidates_c
    public :: check_neighborhood_overlaps_c
    public :: check_mean_pmf_min_counts_c
    public :: check_plateau_condition_c
    public :: check_effect_size_plateau_condition_c
    public :: check_effect_size_plateau_condition_expert_c
    public :: create_mean_pmf_c
    public :: create_mean_pmf_only_c
    public :: bootstrap_histogram_c
    public :: bootstrap_histogram_expert_c
    public :: run_js_comp_test_c
    public :: run_js_comp_test_expert_c
    public :: run_js_comp_test_parameter_search_c
    public :: run_js_comp_test_parameter_search_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):estimate_bin_count(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    !| the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    !| is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    !| range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    !| independently, each clamped on its own to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins and returned as
    !| `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    !| interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    !| skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
    subroutine estimate_bin_count_c(&
            residuals,&
            n_residuals,&
            max_n_reps_all_studies,&
            n_neighbors,&
            shared_residual_range,&
            n_bins,&
            sturges_bins,&
            fd_bins,&
            ierr&
        ) bind(C, name="estimate_bin_count_c")
        use tox_data_integration_js_comp_test, only: estimate_bin_count

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(out), target :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        integer(c_int), intent(out), target :: sturges_bins
            !! Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            !! to the (clamped) sturges_bins when the interquartile range is too close to zero to
            !! divide by (see the is_close guard below)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(sturges_bins)
        M_CHECK_NON_NULL(fd_bins)
        M_CHECK_ARRAY_NON_NULL(residuals, n_residuals)

        call estimate_bin_count(&
            residuals = residuals,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            ierr = ierr&
        )
    end subroutine estimate_bin_count_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):estimate_bin_count_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `estimate_bin_count_helper`. Computes Sturges' rule and
    !| the Freedman-Diaconis rule (without doubling the bin width, since `shared_residual_range`
    !| is already the one-sided half of the full `[-R, R]` histogram range, so dividing the full
    !| range by the undoubled Freedman-Diaconis width already gives the doubled rule's bin count)
    !| independently, each clamped on its own to at most
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] bins and returned as
    !| `sturges_bins`/`fd_bins` for diagnostics, with their maximum returned as `n_bins`. When the
    !| interquartile range is too close to zero to safely divide by, the Freedman-Diaconis term is
    !| skipped and `fd_bins` falls back to the (clamped) Sturges estimate instead of blowing up.
    subroutine estimate_bin_count_expert_c(&
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
        ) bind(C, name="estimate_bin_count_expert_c")
        use tox_data_integration_js_comp_test, only: estimate_bin_count_expert

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: residuals
            !! Pooled signed residuals across all studies, reference points and neighbors
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: residuals_perm
            !! Sorting permutation for `residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(out), target :: n_bins
            !! Estimated number of histogram bins, at least 1 and at most MAX_N_BINS: max(sturges_bins, fd_bins)
        integer(c_int), intent(out), target :: sturges_bins
            !! Sturges' rule estimate alone, at least 1 and at most MAX_N_BINS
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate alone, at least 1 and at most MAX_N_BINS; falls back
            !! to the (clamped) sturges_bins when the interquartile range is too close to zero to
            !! divide by (see the is_close guard below)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(sturges_bins)
        M_CHECK_NON_NULL(fd_bins)
        M_CHECK_ARRAY_NON_NULL(residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(residuals_perm, n_residuals)

        call estimate_bin_count_expert(&
            residuals = residuals,&
            residuals_perm = residuals_perm,&
            n_residuals = n_residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            ierr = ierr&
        )
    end subroutine estimate_bin_count_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy(subroutine)]]
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
    subroutine determine_bin_count_occupancy_c(&
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
        ) bind(C, name="determine_bin_count_occupancy_c")
        use tox_data_integration_js_comp_test, only: determine_bin_count_occupancy

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
            !! every bin at or above min_residuals_per_bin; m_min when occupancy_failed
        logical(c_bool), intent(out), target :: occupancy_failed
            !! `.true.` iff even m_min bins could not satisfy the occupancy criterion (including
            !! the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
            !! caller should reject this neighborhood rather than build a histogram from
            !! selected_n_bins
        real(c_double), intent(out), target :: shared_residual_range_low
            !! This neighborhood's own lower residual-range bound (R_low): the
            !! lower_residual_range_quantile percentile of its own pooled signed residuals --
            !! replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            !! occupancy_failed because every pooled residual is NaN
        real(c_double), intent(out), target :: shared_residual_range_high
            !! This neighborhood's own upper residual-range bound (R_high): the
            !! upper_residual_range_quantile percentile of its own pooled signed residuals --
            !! replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            !! occupancy_failed because every pooled residual is NaN
        integer(c_int), intent(out), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(c_int), intent(out), target :: min_bin_occupancy
            !! Minimum bin count at selected_n_bins; 0 when occupancy_failed
        real(c_double), intent(out), target :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            !! occupancy_failed
        integer(c_int), intent(out), target :: max_bin_occupancy
            !! Maximum bin count at selected_n_bins; 0 when occupancy_failed
        integer(c_int), intent(out), target :: sturges_bins
            !! Sturges' rule estimate for this neighborhood's own pooled residuals, from
            !! estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
            !! estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count the search will ever test (M_max); if a caller passes
            !! `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
            !! relying on an unconfirmed generator capability to bound one optional argument by
            !! another
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(selected_n_bins)
        M_CHECK_NON_NULL(occupancy_failed)
        M_CHECK_NON_NULL(shared_residual_range_low)
        M_CHECK_NON_NULL(shared_residual_range_high)
        M_CHECK_NON_NULL(n_pooled_residuals)
        M_CHECK_NON_NULL(min_bin_occupancy)
        M_CHECK_NON_NULL(mean_bin_occupancy)
        M_CHECK_NON_NULL(max_bin_occupancy)
        M_CHECK_NON_NULL(sturges_bins)
        M_CHECK_NON_NULL(fd_bins)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, n_residuals)

        call determine_bin_count_occupancy(&
            pooled_residuals = pooled_residuals,&
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
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_bin_count_occupancy_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_expert(subroutine)]]
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
    subroutine determine_bin_count_occupancy_expert_c(&
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
        ) bind(C, name="determine_bin_count_occupancy_expert_c")
        use tox_data_integration_js_comp_test, only: determine_bin_count_occupancy_expert

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: pooled_residuals_perm
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The chosen M_j: the largest bin count in [m_min, m_max] whose pooled histogram has
            !! every bin at or above min_residuals_per_bin; m_min when occupancy_failed
        logical(c_bool), intent(out), target :: occupancy_failed
            !! `.true.` iff even m_min bins could not satisfy the occupancy criterion (including
            !! the case where every pooled residual is NaN) -- per Issue #187's FAILURE policy, the
            !! caller should reject this neighborhood rather than build a histogram from
            !! selected_n_bins
        real(c_double), intent(out), target :: shared_residual_range_low
            !! This neighborhood's own lower residual-range bound (R_low): the
            !! lower_residual_range_quantile percentile of its own pooled signed residuals --
            !! replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            !! occupancy_failed because every pooled residual is NaN
        real(c_double), intent(out), target :: shared_residual_range_high
            !! This neighborhood's own upper residual-range bound (R_high): the
            !! upper_residual_range_quantile percentile of its own pooled signed residuals --
            !! replaces the old dataset-wide symmetric shared_residual_range. 0.0 when
            !! occupancy_failed because every pooled residual is NaN
        integer(c_int), intent(out), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(c_int), intent(out), target :: min_bin_occupancy
            !! Minimum bin count at selected_n_bins; 0 when occupancy_failed
        real(c_double), intent(out), target :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            !! occupancy_failed
        integer(c_int), intent(out), target :: max_bin_occupancy
            !! Maximum bin count at selected_n_bins; 0 when occupancy_failed
        integer(c_int), intent(out), target :: sturges_bins
            !! Sturges' rule estimate for this neighborhood's own pooled residuals, from
            !! estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate for this neighborhood's own pooled residuals, from
            !! estimate_bin_count_impl -- a diagnostic only, never part of the search's own decision
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts
            !! Working array: per-bin counts of whichever candidate bin count is currently being
            !! tested, reused throughout the search (256 = MAX_N_BINS, written as a literal since a
            !! generated wrapper's dummy dimension cannot reference a module parameter)
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count the search will ever test (M_max); if a caller passes
            !! `m_max < m_min`, the implementation clamps it up to `m_min` internally rather than
            !! relying on an unconfirmed generator capability to bound one optional argument by
            !! another
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(selected_n_bins)
        M_CHECK_NON_NULL(occupancy_failed)
        M_CHECK_NON_NULL(shared_residual_range_low)
        M_CHECK_NON_NULL(shared_residual_range_high)
        M_CHECK_NON_NULL(n_pooled_residuals)
        M_CHECK_NON_NULL(min_bin_occupancy)
        M_CHECK_NON_NULL(mean_bin_occupancy)
        M_CHECK_NON_NULL(max_bin_occupancy)
        M_CHECK_NON_NULL(sturges_bins)
        M_CHECK_NON_NULL(fd_bins)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals_perm, n_residuals)
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts, 256)

        call determine_bin_count_occupancy_expert(&
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
            upper_residual_range_quantile = upper_residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_bin_count_occupancy_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_exhaustive(subroutine)]]
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
    subroutine determine_bin_count_occupancy_exhaustive_c(&
            pooled_residuals,&
            n_residuals,&
            n_pooled_residuals,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            selected_n_bins,&
            occupancy_failed,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            m_min,&
            m_max,&
            min_residuals_per_bin,&
            ierr&
        ) bind(C, name="determine_bin_count_occupancy_exhaustive_c")
        use tox_data_integration_js_comp_test, only: determine_bin_count_occupancy_exhaustive

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_residuals`.
        real(c_double), intent(in), target :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(c_double), intent(in), target :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! The minimum valid value is `shared_residual_range_low`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            !! above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        logical(c_bool), intent(out), target :: occupancy_failed
            !! `.true.` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
            !! criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
            !! is not itself 0)
        integer(c_int), intent(out), target :: min_bin_occupancy
            !! Minimum bin count at selected_n_bins; 0 when occupancy_failed
        real(c_double), intent(out), target :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins; 0 when occupancy_failed
        integer(c_int), intent(out), target :: max_bin_occupancy
            !! Maximum bin count at selected_n_bins; 0 when occupancy_failed
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
            !! implementation clamps it up to `m_min` internally rather than relying on an
            !! unconfirmed generator capability to bound one optional argument by another
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(n_pooled_residuals)
        M_CHECK_NON_NULL(shared_residual_range_low)
        M_CHECK_NON_NULL(shared_residual_range_high)
        M_CHECK_NON_NULL(selected_n_bins)
        M_CHECK_NON_NULL(occupancy_failed)
        M_CHECK_NON_NULL(min_bin_occupancy)
        M_CHECK_NON_NULL(mean_bin_occupancy)
        M_CHECK_NON_NULL(max_bin_occupancy)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, n_residuals)

        call determine_bin_count_occupancy_exhaustive(&
            pooled_residuals = pooled_residuals,&
            n_residuals = n_residuals,&
            n_pooled_residuals = n_pooled_residuals,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            selected_n_bins = selected_n_bins,&
            occupancy_failed = occupancy_failed,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            ierr = ierr&
        )
    end subroutine determine_bin_count_occupancy_exhaustive_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_exhaustive_expert(subroutine)]]
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
    subroutine determine_bin_count_occupancy_exhaustive_expert_c(&
            pooled_residuals,&
            pooled_residuals_perm,&
            n_residuals,&
            n_pooled_residuals,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            selected_n_bins,&
            occupancy_failed,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            m_min,&
            m_max,&
            min_residuals_per_bin,&
            ierr&
        ) bind(C, name="determine_bin_count_occupancy_exhaustive_expert_c")
        use tox_data_integration_js_comp_test, only: determine_bin_count_occupancy_exhaustive_expert

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: pooled_residuals_perm
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(c_int), intent(in), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_residuals`.
        real(c_double), intent(in), target :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(c_double), intent(in), target :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! The minimum valid value is `shared_residual_range_low`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            !! above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        logical(c_bool), intent(out), target :: occupancy_failed
            !! `.true.` iff no candidate bin count in [m_min, m_max] satisfies the occupancy
            !! criterion (including the case where n_pooled_residuals is 0 and min_residuals_per_bin
            !! is not itself 0)
        integer(c_int), intent(out), target :: min_bin_occupancy
            !! Minimum bin count at selected_n_bins; 0 when occupancy_failed
        real(c_double), intent(out), target :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins; 0 when occupancy_failed
        integer(c_int), intent(out), target :: max_bin_occupancy
            !! Maximum bin count at selected_n_bins; 0 when occupancy_failed
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count tested (M_max); if a caller passes `m_max < m_min`, the
            !! implementation clamps it up to `m_min` internally rather than relying on an
            !! unconfirmed generator capability to bound one optional argument by another
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(n_pooled_residuals)
        M_CHECK_NON_NULL(shared_residual_range_low)
        M_CHECK_NON_NULL(shared_residual_range_high)
        M_CHECK_NON_NULL(selected_n_bins)
        M_CHECK_NON_NULL(occupancy_failed)
        M_CHECK_NON_NULL(min_bin_occupancy)
        M_CHECK_NON_NULL(mean_bin_occupancy)
        M_CHECK_NON_NULL(max_bin_occupancy)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals_perm, n_residuals)

        call determine_bin_count_occupancy_exhaustive_expert(&
            pooled_residuals = pooled_residuals,&
            pooled_residuals_perm = pooled_residuals_perm,&
            n_residuals = n_residuals,&
            n_pooled_residuals = n_pooled_residuals,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            selected_n_bins = selected_n_bins,&
            occupancy_failed = occupancy_failed,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            ierr = ierr&
        )
    end subroutine determine_bin_count_occupancy_exhaustive_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):generate_js_comp_test_candidates(subroutine)]]
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
    subroutine generate_js_comp_test_candidates_c(&
            max_n_genes_all_studies,&
            candidates_n_points_n_neighbors,&
            n_candidates,&
            ierr&
        ) bind(C, name="generate_js_comp_test_candidates_c")
        use tox_data_integration_js_comp_test, only: generate_js_comp_test_candidates

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(2, 16), intent(out), target :: candidates_n_points_n_neighbors
            !! Candidate `[n_points, n_neighbors]` pairs, `n_points` descending
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates
            !! Number of candidate pairs actually filled (at most MAX_CANDIDATE_PAIRS = 16)
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_neighborhood_overlaps(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `test_neighborhood_overlaps_helper`: the first
    !| admissibility gate a candidate `(n_points, n_neighbors)` pair must pass, using the
    !| `[min_idx, max_idx]` neighborhood spans
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| produces. Named `check_*` rather than 125's `test_*`, a deliberate deviation from the
    !| verbatim port: the generated R binding is published under the Fortran name, and the R test
    !| harness (`r/test_helpers.R`'s `run_all_tests`) discovers every `test_`-prefixed name in the
    !| environment as a test case to run with no arguments -- a `test_`-prefixed export would be
    !| swept up and fail every R suite that sources the package, not just this module's own.
    subroutine check_neighborhood_overlaps_c(&
            neighborhood_range,&
            n_points,&
            min_neighbor_overlap,&
            all_have_min_neighbor_overlap,&
            ierr&
        ) bind(C, name="check_neighborhood_overlaps_c")
        use tox_data_integration_js_comp_test, only: check_neighborhood_overlaps

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(2, n_points), intent(in), target :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span, as produced
            !! by construct_neighborhoods_ranged_impl
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out), target :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(all_have_min_neighbor_overlap)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points)

        call check_neighborhood_overlaps(&
            neighborhood_range = neighborhood_range,&
            n_points = n_points,&
            min_neighbor_overlap = min_neighbor_overlap,&
            all_have_min_neighbor_overlap = all_have_min_neighbor_overlap,&
            ierr = ierr&
        )
    end subroutine check_neighborhood_overlaps_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_mean_pmf_min_counts(subroutine)]]
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
    subroutine check_mean_pmf_min_counts_c(&
            mean_pmf_counts,&
            n_bins,&
            n_bins_per_point,&
            n_points,&
            min_count,&
            all_bins_have_min_count,&
            ierr&
        ) bind(C, name="check_mean_pmf_min_counts_c")
        use tox_data_integration_js_comp_test, only: check_mean_pmf_min_counts

        integer(c_int), intent(in), target :: n_bins
            !! Second extent of `mean_pmf_counts` -- the widest bin count any point uses (max_n_bins),
            !! not necessarily every point's own bin count
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: n_bins_per_point
            !! This point's own bin count -- only `mean_pmf_counts(1:n_bins_per_point(i_point), i_point)`
            !! is inspected; columns beyond it are legitimate zero-padding, not failures
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_bins`.
        integer(c_int), intent(in), target :: min_count
            !! Minimum count each bin of the mean pmf must reach
            !! The minimum valid value is `0_int32`.
        logical(c_bool), intent(out), target :: all_bins_have_min_count
            !! `.true.` if every bin within each reference point's own `n_bins_per_point`, at every
            !! reference point, reaches at least `min_count`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(min_count)
        M_CHECK_NON_NULL(all_bins_have_min_count)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)

        call check_mean_pmf_min_counts(&
            mean_pmf_counts = mean_pmf_counts,&
            n_bins = n_bins,&
            n_bins_per_point = n_bins_per_point,&
            n_points = n_points,&
            min_count = min_count,&
            all_bins_have_min_count = all_bins_have_min_count,&
            ierr = ierr&
        )
    end subroutine check_mean_pmf_min_counts_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_plateau_condition(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `check_plateau_condition_helper`. A search over
    !| candidates (finest resolution to coarsest) stops -- "plateaus" -- either when a new
    !| candidate is no better than the previous best (the short-circuit below: keep the previous
    !| best and stop searching), or once the new candidate's confidence-interval overlap with the
    !| previous best meets the condition `join_method` names. `join_method` replaces
    !| 125-stabilize-jscomp's hand-rolled join-method-range validation macro entirely: the mode
    !| table below is itself the validation, checked against exactly the values it names.
    subroutine check_plateau_condition_c(&
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
        ) bind(C, name="check_plateau_condition_c")
        use tox_data_integration_js_comp_test, only: check_plateau_condition
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(2, n_studies), intent(in), target :: confidence_interval
            !! JSD confidence interval `[lower, upper]` from bootstrapping, for the candidate
            !! pair under test
        real(c_double), dimension(2, n_studies), intent(inout), target :: best_candidate_pair_confidence_interval
            !! JSD confidence intervals for the current best candidate pair; overwritten with
            !! `confidence_interval` unless the new candidate is worse
        integer(c_int), intent(inout), target :: best_candidate_index
            !! Candidate-grid index of the current best candidate pair; overwritten with
            !! `candidate_index` unless the new candidate is worse
        integer(c_int), intent(inout), target :: best_exceeded_ci_overlap_count
            !! Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            !! best candidate pair; overwritten unless the new candidate is worse
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(in), target :: candidate_index
            !! Candidate-grid index of the candidate pair that produced `confidence_interval`
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
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
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap an interval in `confidence_interval` must have with its
            !! respective interval in `best_candidate_pair_confidence_interval` to count as
            !! "exceeded"
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out), target :: plateau_found
            !! `.true.` once the new candidate is no better than the previous best, or once
            !! `join_method`'s overlap condition is met by the new candidate
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32) :: join_method_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(best_candidate_index)
        M_CHECK_NON_NULL(best_exceeded_ci_overlap_count)
        M_CHECK_NON_NULL(candidate_index)
        M_CHECK_NON_NULL(succeeding_ci_overlap)
        M_CHECK_NON_NULL(plateau_found)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(best_candidate_pair_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)

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

        call check_plateau_condition(&
            confidence_interval = confidence_interval,&
            best_candidate_pair_confidence_interval = best_candidate_pair_confidence_interval,&
            n_studies = n_studies,&
            best_candidate_index = best_candidate_index,&
            best_exceeded_ci_overlap_count = best_exceeded_ci_overlap_count,&
            candidate_index = candidate_index,&
            join_method = join_method_mode_f,&
            succeeding_ci_overlap = succeeding_ci_overlap,&
            plateau_found = plateau_found,&
            ierr = ierr&
        )
    end subroutine check_plateau_condition_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition(subroutine)]]
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
    subroutine check_effect_size_plateau_condition_c(&
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
        ) bind(C, name="check_effect_size_plateau_condition_c")
        use tox_data_integration_js_comp_test, only: check_effect_size_plateau_condition

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_studies), intent(in), target :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_studies), intent(in), target :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in), target :: has_previous
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_studies), intent(out), target :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out), target :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(has_previous)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(n_consecutive_ok)
        M_CHECK_NON_NULL(delta_median)
        M_CHECK_NON_NULL(delta_max)
        M_CHECK_NON_NULL(plateau_found)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(delta, n_studies)

        call check_effect_size_plateau_condition(&
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
            ierr = ierr&
        )
    end subroutine check_effect_size_plateau_condition_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition_expert(subroutine)]]
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
    subroutine check_effect_size_plateau_condition_expert_c(&
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
        ) bind(C, name="check_effect_size_plateau_condition_expert_c")
        use tox_data_integration_js_comp_test, only: check_effect_size_plateau_condition_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_studies), intent(in), target :: global_js_divergence
            !! Current admissible candidate's observed global JSD per study
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_studies), intent(in), target :: prev_global_js_divergence
            !! Previous admissible candidate's observed global JSD per study; ignored when
            !! `has_previous` is `.false.`
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), intent(in), target :: has_previous
            !! `.false.` for the very first admissible candidate a caller has ever passed in, where
            !! no transition exists to compute a relative change from
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative change across studies must stay under for a
            !! transition to count toward a plateau
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; incremented when this
            !! transition qualifies, reset to zero otherwise (and whenever `has_previous` is
            !! `.false.`)
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_studies), intent(out), target :: delta
            !! Per-study relative JSD change from the previous admissible candidate; `-1.0_real64`
            !! throughout iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_median
            !! Median of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        real(c_double), intent(out), target :: delta_max
            !! Maximum of `delta` across studies; `-1.0_real64` iff `.not. has_previous`
        logical(c_bool), intent(out), target :: plateau_found
            !! `.true.` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_delta_perm
            !! Working array: sorting permutation for `delta`, used to compute `delta_median`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(has_previous)
        M_CHECK_NON_NULL(delta_median_threshold)
        M_CHECK_NON_NULL(delta_max_threshold)
        M_CHECK_NON_NULL(delta_epsilon)
        M_CHECK_NON_NULL(delta_min_consecutive_transitions)
        M_CHECK_NON_NULL(n_consecutive_ok)
        M_CHECK_NON_NULL(delta_median)
        M_CHECK_NON_NULL(delta_max)
        M_CHECK_NON_NULL(plateau_found)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(delta, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_delta_perm, n_studies)

        call check_effect_size_plateau_condition_expert(&
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
            tmp_delta_perm = tmp_delta_perm,&
            ierr = ierr&
        )
    end subroutine check_effect_size_plateau_condition_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):create_mean_pmf(subroutine)]]
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_helper`.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    subroutine create_mean_pmf_c(&
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
        ) bind(C, name="create_mean_pmf_c")
        use tox_data_integration_js_comp_test, only: create_mean_pmf

        integer(c_int), intent(in), target :: n_bins
            !! The array's first extent for `pmfs`/`counts`/`mean_pmf`/`mean_pmf_counts` -- the
            !! widest histogram bin count of any reference point (`max_n_bins`, for pmfs/counts
            !! built by `build_residual_histograms_impl` from its own per-point
            !! `n_bins_per_point`). Every study must share the same per-point bin count for this
            !! to be safe: their zero-padded columns beyond that count then coincide across all
            !! `n_studies`, so the averaged/summed `mean_pmf`/`mean_pmf_counts` are zero there too
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_bins, n_points, n_studies), intent(in), target :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_bins, n_points, n_studies), intent(in), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for `mean_pmf`,
            !! summed across all n_studies
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts, dim=3)`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)

        call create_mean_pmf(&
            pmfs = pmfs,&
            counts = counts,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            included_n_reps = included_n_reps,&
            mean_pmf = mean_pmf,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            mean_pmf_counts = mean_pmf_counts,&
            ierr = ierr&
        )
    end subroutine create_mean_pmf_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):create_mean_pmf_only(subroutine)]]
    !| Ported verbatim from 125-stabilize-jscomp's `create_mean_pmf_only_helper`: useful where the
    !| mean pmf's own counts don't matter, e.g. for the bootstrap confidence interval a later
    !| stage of this port adds.
    !|
    !| Known limitation: averages over all n_studies including the study being compared against it,
    !| rather than a true leave-one-out background as the manuscript specifies. Deliberately ported
    !| as-is from origin/125-stabilize-jscomp; see the project's JSD-Comp-Test follow-up issue for
    !| the fix.
    subroutine create_mean_pmf_only_c(&
            pmfs,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf,&
            ierr&
        ) bind(C, name="create_mean_pmf_only_c")
        use tox_data_integration_js_comp_test, only: create_mean_pmf_only

        integer(c_int), intent(in), target :: n_bins
            !! The array's first extent for `pmfs`/`mean_pmf` -- the widest histogram bin count of
            !! any reference point (`max_n_bins`, for pmfs built by
            !! `build_residual_histograms_impl` from its own per-point `n_bins_per_point`). Every
            !! study must share the same per-point bin count for this to be safe: their
            !! zero-padded columns beyond that count then coincide across all `n_studies`, so the
            !! averaged `mean_pmf` is zero there too
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_bins, n_points, n_studies), intent(in), target :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, built as `mean_pmf = sum(pmfs, dim=3) / n_studies` -- see the
            !! known-limitation note above
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)

        call create_mean_pmf_only(&
            pmfs = pmfs,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf = mean_pmf,&
            ierr = ierr&
        )
    end subroutine create_mean_pmf_only_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):bootstrap_histogram(subroutine)]]
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
    subroutine bootstrap_histogram_c(&
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
        ) bind(C, name="bootstrap_histogram_c")
        use tox_data_integration_js_comp_test, only: bootstrap_histogram

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
            !! otherwise used here
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)

        call bootstrap_histogram(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            confidence_interval = confidence_interval,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):bootstrap_histogram_expert(subroutine)]]
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
    subroutine bootstrap_histogram_expert_c(&
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
        ) bind(C, name="bootstrap_histogram_expert_c")
        use tox_data_integration_js_comp_test, only: bootstrap_histogram_expert

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the pooled/consensus pmf, from
            !! create_mean_pmf_impl -- resampled with replacement each bootstrap
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the pooled pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study --
            !! how many elements are drawn (with replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! Confidence interval to be bootstrapped -- incoming values are the reference values
            !! that seed the top-k/bottom-k heaps (see above); overwritten with the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array used as top-k/bottom-k heaps for efficient percentile detection in
            !! the bootstrapped values -- `(:, 1, :)` the bottom-k (lower bound), `(:, 2, :)` the
            !! top-k (upper bound)
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: tmp_counts
            !! Working array that holds one bootstrap's resampled histogram counts, one study at a time
        real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: tmp_pmfs
            !! Working array that holds one bootstrap's resampled pmfs, all studies
        real(c_double), dimension(n_bins, n_points), intent(out), target :: tmp_mean_pmf
            !! Working array for one bootstrap's own (resampled) mean pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array for one bootstrap's per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_weights
            !! Working array for one bootstrap's per-point global-JSD weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array for one bootstrap's global weighted JSD values
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds to size n_bootstrapping_top_k_jsds; not
            !! otherwise used here
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)

        call bootstrap_histogram_expert(&
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
    end subroutine bootstrap_histogram_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test(subroutine)]]
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
    subroutine run_js_comp_test_c(&
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
        ) bind(C, name="run_js_comp_test_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study (Pass A)
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            !! produced by construct_neighborhoods_ranged_impl (Pass A)
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! This reference point's own selected histogram bin count (Issue #187's `M_j`), from
            !! Pass B's occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood
            !! may use a different bin count
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! This reference point's own lower residual-range bound (R_low), from Pass B's
            !! occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            !! different, asymmetric range (Step 3)
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! This reference point's own upper residual-range bound (R_high), from Pass B's
            !! occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            !! different, asymmetric range (Step 3)
        integer(c_int), intent(out), target :: max_n_bins_per_point
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
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! This reference point's `occupancy_failed` flag from Pass B
            !! (determine_bin_count_occupancy_impl) -- `.true.` iff even `m_min` bins could not
            !! satisfy the occupancy criterion for it. See this routine's own doc block above for
            !! the behavioral asymmetry this implies vs. run_js_comp_test_parameter_search_impl: a
            !! `.true.` point here still gets a real histogram and still contributes to
            !! `global_js_divergence`, it is never rejected
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! This reference point's pooled residual count (N_j) from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! This reference point's minimum bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! This reference point's mean bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! This reference point's maximum bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! This reference point's Sturges' rule bin-count diagnostic from Pass B
            !! (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            !! decision
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! This reference point's Freedman-Diaconis rule bin-count diagnostic from Pass B
            !! (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            !! decision
        real(c_double), dimension(256, n_points, n_studies), intent(out), target :: pmfs
            !! `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`. `256`
            !! = MAX_N_BINS, a fixed ceiling (see `max_n_bins_per_point` above) -- only rows `1:max_n_bins_per_point`
            !! are meaningful; a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves
        integer(c_int), dimension(256, n_points, n_studies), intent(out), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
        real(c_double), dimension(256, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, from create_mean_pmf_impl. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(c_int), dimension(256, n_points), intent(out), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS;
            !! only rows `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible in Pass B's occupancy search, forwarded to
            !! determine_bin_count_occupancy_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(max_n_bins_per_point)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, n_points)
        M_CHECK_ARRAY_NON_NULL(occupancy_failed, n_points)
        M_CHECK_ARRAY_NON_NULL(n_pooled_residuals, n_points)
        M_CHECK_ARRAY_NON_NULL(min_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(max_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(sturges_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(fd_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(pmfs, 256 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, 256 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, 256 * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, 256 * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)

        call run_js_comp_test(&
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
    end subroutine run_js_comp_test_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_expert(subroutine)]]
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
    subroutine run_js_comp_test_expert_c(&
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
        ) bind(C, name="run_js_comp_test_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points for neighborhoods
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Per-gene mean expression values for all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means_perms
            !! Per-study sorting permutation for `gene_means` (ascending, NaN last)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study
            !! NaN is permitted for this value.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study (Pass A)
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span, as
            !! produced by construct_neighborhoods_ranged_impl (Pass A)
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! This reference point's own selected histogram bin count (Issue #187's `M_j`), from
            !! Pass B's occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood
            !! may use a different bin count
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! This reference point's own lower residual-range bound (R_low), from Pass B's
            !! occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            !! different, asymmetric range (Step 3)
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! This reference point's own upper residual-range bound (R_high), from Pass B's
            !! occupancy search (determine_bin_count_occupancy_impl) -- every neighborhood may use a
            !! different, asymmetric range (Step 3)
        integer(c_int), intent(out), target :: max_n_bins_per_point
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
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! This reference point's `occupancy_failed` flag from Pass B
            !! (determine_bin_count_occupancy_impl) -- `.true.` iff even `m_min` bins could not
            !! satisfy the occupancy criterion for it. See this routine's own doc block above for
            !! the behavioral asymmetry this implies vs. run_js_comp_test_parameter_search_impl: a
            !! `.true.` point here still gets a real histogram and still contributes to
            !! `global_js_divergence`, it is never rejected
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! This reference point's pooled residual count (N_j) from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! This reference point's minimum bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! This reference point's mean bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! This reference point's maximum bin occupancy at `n_bins_per_point`, from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! This reference point's Sturges' rule bin-count diagnostic from Pass B
            !! (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            !! decision
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! This reference point's Freedman-Diaconis rule bin-count diagnostic from Pass B
            !! (determine_bin_count_occupancy_impl) -- never part of the occupancy search's own
            !! decision
        real(c_double), dimension(256, n_points, n_studies), intent(out), target :: pmfs
            !! `counts` normalized to `0 <= pmfs(:, :, i) <= 1` and `sum(pmfs(:, j, i)) == 1`. `256`
            !! = MAX_N_BINS, a fixed ceiling (see `max_n_bins_per_point` above) -- only rows `1:max_n_bins_per_point`
            !! are meaningful; a Python/R caller must slice `[:max_n_bins_per_point, ...]` themselves
        integer(c_int), dimension(256, n_points, n_studies), intent(out), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study
        real(c_double), dimension(256, n_points), intent(out), target :: mean_pmf
            !! The consensus pmf, from create_mean_pmf_impl. `256` = MAX_N_BINS; only rows
            !! `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(c_int), dimension(256, n_points), intent(out), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf. `256` = MAX_N_BINS;
            !! only rows `1:max_n_bins_per_point` are meaningful -- see `pmfs` above
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study from gjct_permutation_test_impl
        real(c_double), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values, reused per study
            !! (Pass C)
        integer(c_int), dimension(n_points, 256), intent(out), target :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts from
            !! build_residual_histograms_impl. `256` = MAX_N_BINS, written as a literal because a
            !! generated wrapper's dummy dimension cannot reference a module parameter
        real(c_double), dimension(n_points, 256), intent(out), target :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf, reused both for
            !! build_residual_histograms_impl's output and for calc_pmf_impl's re-derived pmf.
            !! `256` = MAX_N_BINS, see tmp_counts_point_major above
        real(c_double), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals across every neighbor and
            !! every study (Pass B), reused per point -- one small buffer, not one per point, since
            !! Pass B is a deliberate sequential loop (see the implementation body's own comment)
        integer(c_int), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array forwarded to determine_bin_count_occupancy_impl's own per-bin-count
            !! search scratch, reused per point. `256` = MAX_N_BINS, matching
            !! determine_bin_count_occupancy_impl's own tmp_bin_counts dummy -- written as a literal
            !! because a generated wrapper's dummy dimension cannot reference a module parameter
        integer(c_int), dimension(256, n_points), intent(out), target :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to gjct_permutation_test_impl's own resampling pool. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(256, n_points), intent(out), target :: tmp_permutation_counts
            !! Working array forwarded to gjct_permutation_test_impl's own per-study resampled
            !! counts. `256` = MAX_N_BINS, see tmp_counts_point_major above
        real(c_double), dimension(256, n_points, n_studies), intent(out), target :: tmp_permutation_pmfs
            !! Working array forwarded to gjct_permutation_test_impl's own resampled pmfs. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_js_divergences
            !! Working array forwarded to gjct_permutation_test_impl's own per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_weights
            !! Working array forwarded to gjct_permutation_test_impl's own per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_permutation_global_js_divergence
            !! Working array forwarded to gjct_permutation_test_impl's own resampled global JSD values
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations, forwarded to gjct_permutation_test_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible in Pass B's occupancy search, forwarded to
            !! determine_bin_count_occupancy_impl
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator for
            !! the permutation test

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(max_n_bins_per_point)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_NON_NULL(min_residuals_per_bin)
        M_CHECK_NON_NULL(m_min)
        M_CHECK_NON_NULL(m_max)
        M_CHECK_NON_NULL(gamma_occupancy)
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(gene_means_perms, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, n_points)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, n_points)
        M_CHECK_ARRAY_NON_NULL(occupancy_failed, n_points)
        M_CHECK_ARRAY_NON_NULL(n_pooled_residuals, n_points)
        M_CHECK_ARRAY_NON_NULL(min_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(max_bin_occupancy, n_points)
        M_CHECK_ARRAY_NON_NULL(sturges_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(fd_bins, n_points)
        M_CHECK_ARRAY_NON_NULL(pmfs, 256 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, 256 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, 256 * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, 256 * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_counts_point_major, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_point_major, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_mean_pmf_counts, 256 * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_counts, 256 * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_pmfs, 256 * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_global_js_divergence, n_studies)

        call run_js_comp_test_expert(&
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
    end subroutine run_js_comp_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search(subroutine)]]
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
        ) bind(C, name="run_js_comp_test_parameter_search_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_parameter_search
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

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
            !! Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
            !! per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
            !! jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
            !! is no longer purely an internal sizing detail the plain wrapper can compute and
            !! hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
            !! rather than AUTO here now
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
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
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
            !! because anything returned is sized by it (nothing is), but because it comes from the
            !! same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
            !! now that that call can no longer run automatically inside this wrapper, splitting
            !! this one back into an AUTO call would just be a second, redundant call to the same
            !! routine for no benefit
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: n_points
            !! The finally chosen candidate's `n_points`
        integer(c_int), intent(out), target :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: n_bins_per_point
            !! The finally chosen candidate's per-point histogram bin count, one per reference
            !! point (Issue #187: every neighborhood may use a different bin count). Only the
            !! leading `n_points` entries are meaningful, mirroring how `n_points`/`n_neighbors`
            !! above are the finally chosen candidate's own values
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_low
            !! The finally chosen candidate's per-point lower residual-range bound (R_low), one per
            !! reference point (Step 3: every neighborhood may use a different, asymmetric range).
            !! Only the leading `n_points` entries are meaningful, mirroring `n_bins_per_point`
            !! above; `0.0_real64` throughout in the two genuinely-degenerate cases where Pass B
            !! never ran for the returned candidate (see the final three-way branch's own comments)
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_high
            !! The finally chosen candidate's per-point upper residual-range bound (R_high),
            !! mirroring `shared_residual_range_low` above in every respect
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
            !! The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            !! `-1.0_real64` throughout only when `plateau_established` is `.false.` and no
            !! smallest-bootstrap-uncertainty candidate could be substituted either (see
            !! `plateau_established`)
        logical(c_bool), intent(out), target :: plateau_established
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
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of candidates that passed both admissibility gates and got a JSD/confidence
            !! interval computed before the search stopped (by plateau or grid exhaustion) -- the
            !! number of leading, valid columns/elements in every `trace_*` array below. This is
            !! Issue #178's own index `t` domain: "position in the ordered sequence of ADMISSIBLE
            !! parameter pairs" -- a candidate that failed either gate has no `trace_*` entry at
            !! all, rather than a zero-filled one
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
            !! arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
            !! candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
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
            !! (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
            !! formula omits this floor, but the same near-zero-JSD instability that motivates
            !! `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
            !! destabilizes any ratio that divides by it, whichever candidate's `J` it is)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_selected_n_bins
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
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out), target :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
            !! from determine_bin_count_occupancy_impl (never part of the occupancy search's own
            !! decision). Jagged per candidate column exactly as trace_selected_n_bins above --
            !! only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
            !! must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
            !! diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
            !! search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
            !! above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
            !! caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate. Reuses Issue #187's occupancy-search default rather than an
            !! independently-tunable threshold of its own: once
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! wires real per-neighborhood bin counts in, a separate laxer threshold here would
            !! silently let a candidate the occupancy search already marked `occupancy_failed`
            !! pass this gate anyway, defeating the FAILURE-detection mechanism
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                      | Value                                                                                 |
            !! |-------------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only       | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                          | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_ci_overlap'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f

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
            ierr = ierr&
        )
    end subroutine run_js_comp_test_parameter_search_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search_expert(subroutine)]]
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
        ) bind(C, name="run_js_comp_test_parameter_search_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_parameter_search_expert
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

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
            !! Exact upper bound on the grid's first (largest) `n_points` candidate. Issue #187's
            !! per-point outputs below (`n_bins_per_point`, `trace_selected_n_bins`, and the other
            !! jagged `trace_*` arrays) are sized by this argument, so unlike before Issue #187 it
            !! is no longer purely an internal sizing detail the plain wrapper can compute and
            !! hide -- the caller must know it up front to receive those arrays, hence JUST_INFO
            !! rather than AUTO here now
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Safe upper bound on the grid's largest `n_neighbors` candidate. Also JUST_INFO, not
            !! because anything returned is sized by it (nothing is), but because it comes from the
            !! same `calc_js_comp_test_candidate_bounds` call as `max_n_points_candidate` above --
            !! now that that call can no longer run automatically inside this wrapper, splitting
            !! this one back into an AUTO call would just be a second, redundant call to the same
            !! routine for no benefit
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of elements kept at each end of the bootstrap distribution (top-k/bottom-k
            !! heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
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
            !! The way to evaluate all studies' confidence-interval overlaps for the plateau
            !! condition, forwarded to check_plateau_condition_impl
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(out), target :: n_points
            !! The finally chosen candidate's `n_points`
        integer(c_int), intent(out), target :: n_neighbors
            !! The finally chosen candidate's `n_neighbors`
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: n_bins_per_point
            !! The finally chosen candidate's per-point histogram bin count, one per reference
            !! point (Issue #187: every neighborhood may use a different bin count). Only the
            !! leading `n_points` entries are meaningful, mirroring how `n_points`/`n_neighbors`
            !! above are the finally chosen candidate's own values
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_low
            !! The finally chosen candidate's per-point lower residual-range bound (R_low), one per
            !! reference point (Step 3: every neighborhood may use a different, asymmetric range).
            !! Only the leading `n_points` entries are meaningful, mirroring `n_bins_per_point`
            !! above; `0.0_real64` throughout in the two genuinely-degenerate cases where Pass B
            !! never ran for the returned candidate (see the final three-way branch's own comments)
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_high
            !! The finally chosen candidate's per-point upper residual-range bound (R_high),
            !! mirroring `shared_residual_range_low` above in every respect
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
            !! The bootstrapped JSD confidence interval for the finally chosen candidate pair;
            !! `-1.0_real64` throughout only when `plateau_established` is `.false.` and no
            !! smallest-bootstrap-uncertainty candidate could be substituted either (see
            !! `plateau_established`)
        logical(c_bool), intent(out), target :: plateau_established
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
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of candidates that passed both admissibility gates and got a JSD/confidence
            !! interval computed before the search stopped (by plateau or grid exhaustion) -- the
            !! number of leading, valid columns/elements in every `trace_*` array below. This is
            !! Issue #178's own index `t` domain: "position in the ordered sequence of ADMISSIBLE
            !! parameter pairs" -- a candidate that failed either gate has no `trace_*` entry at
            !! all, rather than a zero-filled one
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate `n_points`, one entry per column of the other `trace_*`
            !! arrays. `16` = MAX_CANDIDATE_PAIRS, written as a literal for the same reason
            !! candidates_n_points_n_neighbors is in generate_js_comp_test_candidates_impl
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`, paired with trace_n_points above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}` in Issue #178)
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
            !! (`W_{i,t} / J_{i,t}`), denominator floored at `delta_epsilon` -- the issue's own
            !! formula omits this floor, but the same near-zero-JSD instability that motivates
            !! `delta_epsilon` in the `Delta_{i,t}` formula applies here too (a near-zero `J`
            !! destabilizes any ratio that divides by it, whichever candidate's `J` it is)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`), from check_effect_size_plateau_condition_impl;
            !! `-1.0_real64` throughout at the first admissible candidate specifically (no
            !! predecessor to diff against) -- every other column within `1:n_admissible_evaluated`
            !! holds a real value
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies (Delta-tilde_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(16), intent(out), target :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies (Delta^max_t);
            !! `-1.0_real64` at the first admissible candidate, see trace_delta above
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_selected_n_bins
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
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out), target :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point `occupancy_failed` flag from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`) from
            !! determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy at
            !! trace_selected_n_bins, from determine_bin_count_occupancy_impl. Jagged per candidate
            !! column exactly as trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are
            !! meaningful for column `t`; a Python/R caller must slice `[:trace_n_points[t], t]`
            !! themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule bin-count diagnostic
            !! from determine_bin_count_occupancy_impl (never part of the occupancy search's own
            !! decision). Jagged per candidate column exactly as trace_selected_n_bins above --
            !! only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R caller
            !! must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule bin-count
            !! diagnostic from determine_bin_count_occupancy_impl (never part of the occupancy
            !! search's own decision). Jagged per candidate column exactly as trace_selected_n_bins
            !! above -- only rows `1:trace_n_points(t)` are meaningful for column `t`; a Python/R
            !! caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high)
            !! from determine_bin_count_occupancy_impl. Jagged per candidate column exactly as
            !! trace_selected_n_bins above -- only rows `1:trace_n_points(t)` are meaningful for
            !! column `t`; a Python/R caller must slice `[:trace_n_points[t], t]` themselves
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation for the flattened, all-studies-pooled `gene_means`
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_x_star
            !! Working array: reference points for the candidate whose `n_points` is current,
            !! recomputed only when `n_points` changes between candidates
        integer(c_int), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out), target :: tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate,
            !! retained simultaneously (unlike the single-study-reused buffer this replaces) so
            !! Pass B below can pool residuals across studies for one reference point at a time,
            !! and Pass C can re-gather each study's own residuals from the already-known indices
        integer(c_int), dimension(2, max_n_points_candidate), intent(out), target :: tmp_neighborhood_range
            !! Working array: one study's `[min_idx, max_idx]` neighborhood spans for the current
            !! candidate, reused per study
        real(c_double), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residual values for the current
            !! candidate, reused per study
        integer(c_int), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_counts_point_major
            !! Working array: one study's point-major histogram counts for the current candidate,
            !! reused per study. The `256` is
            !! [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]], written as a
            !! literal because a dimension naming a module parameter the generated wrapper never
            !! `use`s would not compile there -- the same reason MAX_CANDIDATE_PAIRS is written as
            !! `16` on generate_js_comp_test_candidates_impl's own dummy arguments
        real(c_double), dimension(max_n_points_candidate, 256), intent(out), target :: tmp_pmf_point_major
            !! Working array: one study's point-major pmf for the current candidate, reused per
            !! study. `256` = MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_bins_per_point
            !! Working array: this candidate's per-point histogram bin count, decided by Pass B's
            !! occupancy search
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]])
            !! for each reference point independently
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_low
            !! Working array: this candidate's per-point lower residual-range bound (R_low), from
            !! Pass B's occupancy search, mirroring tmp_n_bins_per_point above
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_high
            !! Working array: this candidate's per-point upper residual-range bound (R_high),
            !! mirroring tmp_shared_residual_range_low above
        real(c_double), dimension(256, max_n_points_candidate, n_studies), intent(out), target :: tmp_pmfs
            !! Working array: every study's bin-major pmf for the current candidate. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(256, max_n_points_candidate, n_studies), intent(out), target :: tmp_counts
            !! Working array: every study's bin-major histogram counts for the current candidate.
            !! `256` = MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_included_n_reps
            !! Working array: every study's included-replicate counts for the current candidate
        real(c_double), dimension(256, max_n_points_candidate), intent(out), target :: tmp_mean_pmf
            !! Working array: the current candidate's consensus pmf. `256` = MAX_N_BINS, see
            !! tmp_counts_point_major above
        integer(c_int), dimension(256, max_n_points_candidate), intent(out), target :: tmp_mean_pmf_counts
            !! Working array: the current candidate's consensus histogram counts. `256` =
            !! MAX_N_BINS, see tmp_counts_point_major above
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_pmf_included_n_reps
            !! Working array: the current candidate's consensus included-replicate counts
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array: the current candidate's per-reference-point JSD values, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(c_double), dimension(max_n_points_candidate, n_studies), intent(out), target :: tmp_weights
            !! Working array: the current candidate's per-reference-point weights, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array: the current candidate's observed global JSD per study, reused as
            !! bootstrap_histogram_impl's own scratch once consumed
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_confidence_interval
            !! Working array: the current candidate's confidence interval, seeded with the observed
            !! global JSD and then bootstrapped in place
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array forwarded to bootstrap_histogram_impl's own top-k/bottom-k heaps
        real(c_double), dimension(n_studies), intent(out), target :: tmp_prev_global_js_divergence
            !! Working array: the previous admissible candidate's observed global JSD per study,
            !! forwarded to check_effect_size_plateau_condition_impl
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_delta_perm
            !! Working array forwarded to check_effect_size_plateau_condition_impl's own median
            !! computation
        real(c_double), dimension(2, n_studies), intent(out), target :: tmp_best_uncertainty_confidence_interval
            !! Working array: the confidence interval of the admissible candidate with the smallest
            !! bootstrapped uncertainty seen so far, used only internally by the no-plateau fallback
            !! (see plateau_established)
        real(c_double), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals across every neighbor and
            !! every study (Pass B), reused per point -- one small buffer, not one per point, since
            !! Pass B is a deliberate sequential loop (see the module-internal doc comment on the
            !! implementation body)
        integer(c_int), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array forwarded to determine_bin_count_occupancy_impl's own per-bin-count
            !! search scratch, reused per point. `256` = MAX_N_BINS, matching
            !! determine_bin_count_occupancy_impl's own tmp_bin_counts dummy -- written as a
            !! literal because a generated wrapper's dummy dimension cannot reference a module
            !! parameter
        logical(c_bool), dimension(max_n_points_candidate), intent(out), target :: tmp_occupancy_failed
            !! Working array: this candidate's per-point occupancy_failed flag from Pass B
            !! (determine_bin_count_occupancy_impl)
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_pooled_residuals
            !! Working array: this candidate's per-point pooled residual count (N_j) from Pass B
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_min_bin_occupancy
            !! Working array: this candidate's per-point minimum bin occupancy from Pass B
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_bin_occupancy
            !! Working array: this candidate's per-point mean bin occupancy from Pass B
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_max_bin_occupancy
            !! Working array: this candidate's per-point maximum bin occupancy from Pass B
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_sturges_bins
            !! Working array: this candidate's per-point Sturges' rule diagnostic from Pass B
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_fd_bins
            !! Working array: this candidate's per-point Freedman-Diaconis rule diagnostic from Pass B
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_n_bins_per_point
            !! Working array: a snapshot of tmp_n_bins_per_point for whichever candidate is
            !! currently best_candidate_index, updated in lockstep with
            !! check_plateau_condition_impl's own best-candidate bookkeeping (and with the
            !! effect-size override below) -- never read from tmp_n_bins_per_point at loop exit,
            !! which would be stale once the loop has moved on to a later, non-best candidate
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: a snapshot of tmp_n_bins_per_point for whichever admissible candidate
            !! currently has the smallest bootstrapped uncertainty, mirroring how
            !! tmp_best_uncertainty_confidence_interval already works
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_low
            !! Working array: a snapshot of tmp_shared_residual_range_low for whichever candidate is
            !! currently best_candidate_index, updated at the exact same sites as
            !! tmp_best_n_bins_per_point above
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_high
            !! Working array: a snapshot of tmp_shared_residual_range_high, mirroring
            !! tmp_best_shared_residual_range_low above
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: a snapshot of tmp_shared_residual_range_low for whichever admissible
            !! candidate currently has the smallest bootstrapped uncertainty, updated at the exact
            !! same site as tmp_best_uncertainty_n_bins_per_point above
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: a snapshot of tmp_shared_residual_range_high, mirroring
            !! tmp_best_uncertainty_shared_residual_range_low above
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum count each bin of the consensus pmf must reach to pass the second
            !! admissibility gate. Reuses Issue #187's occupancy-search default rather than an
            !! independently-tunable threshold of its own: once
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! wires real per-neighborhood bin counts in, a separate laxer threshold here would
            !! silently let a candidate the occupancy search already marked `occupancy_failed`
            !! pass this gate anyway, defeating the FAILURE-detection mechanism
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have to pass the first
            !! admissibility gate
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap a candidate's confidence interval must have with the
            !! running best, per `join_method`, to plateau
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.9_real64`.
        character(len=1, kind=c_char), dimension(19), intent(in), target :: plateau_mode
            !! Which plateau criterion decides when the search stops
            !!
            !! | Mode                                      | Value                                                                                 |
            !! |-------------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only (pre-Issue-#178 behavior) | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only       | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                          | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
            !! The default value is `'plateau_ci_overlap'`.
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound the median relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound the largest relative JSD change across studies must stay under for a
            !! transition to count toward an effect-size plateau, forwarded to
            !! check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `0.10_real64`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous admissible
            !! candidate's JSD was zero, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `1.0e-10_real64`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare an effect-size
            !! plateau, forwarded to check_effect_size_plateau_condition_impl
            !! The minimum valid value is `1_int32`.
            !! The default value is `2_int32`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count Pass B's occupancy search will ever test (M_min),
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count Pass B's occupancy search will ever test (M_max),
            !! forwarded to determine_bin_count_occupancy_impl; if a caller passes `m_max < m_min`,
            !! determine_bin_count_occupancy_impl clamps it up to `m_min` internally
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        real(c_double), intent(in), target :: gamma_occupancy
            !! Geometric growth factor for Pass B's occupancy search's coarse search stage,
            !! forwarded to determine_bin_count_occupancy_impl; must exceed 1 or the search never
            !! advances
            !! The minimum valid value is `above(1.0_real64)`.
            !! The default value is `1.25_real64`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own lower residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for each reference point's own upper residual-range bound,
            !! forwarded to determine_bin_count_occupancy_impl
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Forwarded to calc_js_comp_test_n_top_k_jsds (sizing n_bootstrapping_top_k_jsds) and
            !! to bootstrap_histogram_impl itself
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; folds any GSL allocation failure bootstrap_histogram_impl reports
        integer(int32) :: join_method_mode_f
        integer(int32) :: plateau_mode_mode_f

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
        M_CHECK_ARRAY_NON_NULL(tmp_counts_point_major, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_point_major, max_n_points_candidate * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, 256 * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, 256 * max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps, max_n_points_candidate * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, 256 * max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, 256 * max_n_points_candidate)
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
            ierr = ierr&
        )
    end subroutine run_js_comp_test_parameter_search_expert_c

end module tox_data_integration_js_comp_test_c
#endif
