#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_binning(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bin Counts
!|
!| How many histogram bins each reference point of the JS-Comp-Test gets, chosen from the pooled
!| residuals of all studies in its neighborhood. The occupancy search
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
!| (the default) takes the largest bin count whose every bin holds at least a minimum number of
!| residuals; Stone's cross-validation
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]
!| takes the bin count with the smallest estimated risk instead (`bin_count_mode`).
!| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
!| is the brute-force reference of the occupancy search, and
!| [[tox_data_integration_js_comp_test_binning_impl(module):estimate_bin_count_impl(interface)]]
!| the Sturges/Freedman-Diaconis rule-of-thumb estimate. Every bin count is at most
!| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]].
module tox_data_integration_js_comp_test_binning_c
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
    public :: determine_bin_count_stone_c
    public :: determine_bin_count_stone_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):estimate_bin_count(subroutine)]]
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]] and returned as
    !| `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    !| `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    !| from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    !| bin count.
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
        use tox_data_integration_js_comp_test_binning, only: estimate_bin_count

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
            !! to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):estimate_bin_count_expert(subroutine)]]
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]] and returned as
    !| `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    !| `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    !| from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    !| bin count.
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
        use tox_data_integration_js_comp_test_binning, only: estimate_bin_count_expert

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
            !! to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):determine_bin_count_occupancy(subroutine)]]
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):estimate_bin_count_impl(interface)]] and are
    !| diagnostics only.
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
        use tox_data_integration_js_comp_test_binning, only: determine_bin_count_occupancy

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
            !! `.true.` iff even m_min bins fail the occupancy criterion (or every pooled residual is
            !! NaN); the caller should then reject this neighborhood
        real(c_double), intent(out), target :: shared_residual_range_low
            !! This neighborhood's lower residual-range bound (R_low): the
            !! lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        real(c_double), intent(out), target :: shared_residual_range_high
            !! This neighborhood's upper residual-range bound (R_high): the
            !! upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
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
            !! Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
            !! smaller
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):determine_bin_count_occupancy_expert(subroutine)]]
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):estimate_bin_count_impl(interface)]] and are
    !| diagnostics only.
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
        use tox_data_integration_js_comp_test_binning, only: determine_bin_count_occupancy_expert

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
            !! `.true.` iff even m_min bins fail the occupancy criterion (or every pooled residual is
            !! NaN); the caller should then reject this neighborhood
        real(c_double), intent(out), target :: shared_residual_range_low
            !! This neighborhood's lower residual-range bound (R_low): the
            !! lower_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
        real(c_double), intent(out), target :: shared_residual_range_high
            !! This neighborhood's upper residual-range bound (R_high): the
            !! upper_residual_range_quantile quantile of its pooled residuals; 0.0 if all are NaN
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
            !! Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts
            !! Working array: per-bin counts of the candidate bin count under test
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
            !! smaller
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
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):determine_bin_count_occupancy_exhaustive(subroutine)]]
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]],
    !| whose faster search could in principle miss a larger admissible `M`. Takes the range and
    !| `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
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
        use tox_data_integration_js_comp_test_binning, only: determine_bin_count_occupancy_exhaustive

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from a prior
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_residuals`.
        real(c_double), intent(in), target :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(c_double), intent(in), target :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! The minimum valid value is `shared_residual_range_low`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            !! above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        logical(c_bool), intent(out), target :: occupancy_failed
            !! `.true.` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion
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
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):determine_bin_count_occupancy_exhaustive_expert(subroutine)]]
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]],
    !| whose faster search could in principle miss a larger admissible `M`. Takes the range and
    !| `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
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
        use tox_data_integration_js_comp_test_binning, only: determine_bin_count_occupancy_exhaustive_expert

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
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_residuals`.
        real(c_double), intent(in), target :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(c_double), intent(in), target :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! The minimum valid value is `shared_residual_range_low`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The largest bin count in [m_min, m_max] whose pooled histogram has every bin at or
            !! above min_residuals_per_bin, found by exhaustive search; m_min when occupancy_failed
        logical(c_bool), intent(out), target :: occupancy_failed
            !! `.true.` iff no candidate bin count in [m_min, m_max] satisfies the occupancy criterion
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
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):determine_bin_count_stone(subroutine)]]
    !| Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
    !| Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
    !| `n` values over `m` bins of width `h`, the bin count minimising
    !| `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
    !| chosen. An alternative to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]],
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] (the
    !| default) minimises the sum of every study's own score `J_k(m)`,
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]]
    !| the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] the
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
    subroutine determine_bin_count_stone_c(&
            pooled_residuals,&
            n_residuals,&
            residual_to_study,&
            n_studies,&
            max_n_reps_all_studies,&
            n_neighbors,&
            selected_n_bins,&
            stone_failed,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            bin_count_mode,&
            m_min,&
            m_max,&
            min_residuals_per_bin,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_bin_count_stone_c")
        use tox_data_integration_js_comp_test_binning, only: determine_bin_count_stone
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: residual_to_study
            !! Study each pooled residual comes from
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_studies`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
            !! Stone score; m_min when stone_failed
        logical(c_bool), intent(out), target :: stone_failed
            !! `.true.` iff no bin count could be chosen (see above); the caller should then reject
            !! this neighborhood rather than build a histogram from selected_n_bins
        real(c_double), intent(out), target :: shared_residual_range_low
            !! This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN
        real(c_double), intent(out), target :: shared_residual_range_high
            !! This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN
        integer(c_int), intent(out), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(c_int), intent(out), target :: min_bin_occupancy
            !! Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        real(c_double), intent(out), target :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            !! stone_failed
        integer(c_int), intent(out), target :: max_bin_occupancy
            !! Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        integer(c_int), intent(out), target :: sturges_bins
            !! Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How the studies' scores are combined into one bin count
            !!
            !! | Mode                                                                                  | Value                                                                                                   |
            !! |---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Minimise the sum of the per-study scores (default)                                    | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled                                   | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_stone_equal'`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin of the pooled clamped histogram must
            !! reach for a candidate bin count to be admissible (0 turns the floor off)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32) :: bin_count_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(selected_n_bins)
        M_CHECK_NON_NULL(stone_failed)
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
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(residual_to_study, n_residuals)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
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

        call determine_bin_count_stone(&
            pooled_residuals = pooled_residuals,&
            n_residuals = n_residuals,&
            residual_to_study = residual_to_study,&
            n_studies = n_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            selected_n_bins = selected_n_bins,&
            stone_failed = stone_failed,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            bin_count_mode = bin_count_mode_mode_f,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_bin_count_stone_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_binning(module):determine_bin_count_stone_expert(subroutine)]]
    !| Chooses the number of equal-width histogram bins for one neighborhood's pooled residuals by
    !| Stone's leave-one-out cross-validation (Duke STA 293.3 lecture notes, ch. 6, Theorem 2): for
    !| `n` values over `m` bins of width `h`, the bin count minimising
    !| `J(m) = 2/(h(n-1)) - (n+1)/(h(n-1)) * sum_j p_j^2` (`p_j` the share of values in bin `j`) is
    !| chosen. An alternative to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]],
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] (the
    !| default) minimises the sum of every study's own score `J_k(m)`,
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]]
    !| the sum of `n_k * J_k(m)` with `n_k` the study's inside count, and
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] the
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
    subroutine determine_bin_count_stone_expert_c(&
            pooled_residuals,&
            pooled_residuals_perm,&
            n_residuals,&
            residual_to_study,&
            n_studies,&
            max_n_reps_all_studies,&
            n_neighbors,&
            selected_n_bins,&
            stone_failed,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            n_pooled_residuals,&
            min_bin_occupancy,&
            mean_bin_occupancy,&
            max_bin_occupancy,&
            sturges_bins,&
            fd_bins,&
            tmp_bin_counts,&
            tmp_residual_to_study_perm,&
            tmp_study_block_sizes,&
            bin_count_mode,&
            m_min,&
            m_max,&
            min_residuals_per_bin,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_bin_count_stone_expert_c")
        use tox_data_integration_js_comp_test_binning, only: determine_bin_count_stone_expert
        use tox_data_integration_js_comp_test_binning_impl, only: MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

        integer(c_int), intent(in), target :: n_residuals
            !! Number of pooled residuals
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_residuals), intent(in), target :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_residuals), intent(in), target :: pooled_residuals_perm
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(c_int), dimension(n_residuals), intent(in), target :: residual_to_study
            !! Study each pooled residual comes from
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_studies`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: selected_n_bins
            !! The chosen bin count: the admissible bin count in [m_min, m_max] with the smallest
            !! Stone score; m_min when stone_failed
        logical(c_bool), intent(out), target :: stone_failed
            !! `.true.` iff no bin count could be chosen (see above); the caller should then reject
            !! this neighborhood rather than build a histogram from selected_n_bins
        real(c_double), intent(out), target :: shared_residual_range_low
            !! This neighborhood's lower residual-range bound (R_low); 0.0 if every residual is NaN
        real(c_double), intent(out), target :: shared_residual_range_high
            !! This neighborhood's upper residual-range bound (R_high); 0.0 if every residual is NaN
        integer(c_int), intent(out), target :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j)
        integer(c_int), intent(out), target :: min_bin_occupancy
            !! Minimum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        real(c_double), intent(out), target :: mean_bin_occupancy
            !! Mean bin count at selected_n_bins (== n_pooled_residuals / selected_n_bins); 0 when
            !! stone_failed
        integer(c_int), intent(out), target :: max_bin_occupancy
            !! Maximum bin count of the pooled clamped histogram at selected_n_bins; 0 when stone_failed
        integer(c_int), intent(out), target :: sturges_bins
            !! Sturges' rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), intent(out), target :: fd_bins
            !! Freedman-Diaconis rule estimate from estimate_bin_count_impl (diagnostic only)
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts
            !! Working array: one study's bin counts at the bin count under test
        integer(c_int), dimension(n_residuals), intent(out), target :: tmp_residual_to_study_perm
            !! Working array: indices of `pooled_residuals` grouped by study, each study sorted
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array: number of pooled residuals of each study
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How the studies' scores are combined into one bin count
            !!
            !! | Mode                                                                                  | Value                                                                                                   |
            !! |---------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------|
            !! | Minimise the sum of the per-study scores (default)                                    | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled                                   | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_stone_equal'`.
        integer(c_int), intent(in), target :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(c_int), intent(in), target :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin of the pooled clamped histogram must
            !! reach for a candidate bin count to be admissible (0 turns the floor off)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile; otherwise the call fails with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(c_double), intent(in), target :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32) :: bin_count_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_residuals)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(selected_n_bins)
        M_CHECK_NON_NULL(stone_failed)
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
        M_CHECK_NON_NULL(lower_residual_range_quantile)
        M_CHECK_NON_NULL(upper_residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, n_residuals)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals_perm, n_residuals)
        M_CHECK_ARRAY_NON_NULL(residual_to_study, n_residuals)
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, n_residuals)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

        block
            character(len=:), pointer :: bin_count_mode_f
            bin_count_mode_f => c_char_as_view(bin_count_mode)

            select case (bin_count_mode_f)
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

        call determine_bin_count_stone_expert(&
            pooled_residuals = pooled_residuals,&
            pooled_residuals_perm = pooled_residuals_perm,&
            n_residuals = n_residuals,&
            residual_to_study = residual_to_study,&
            n_studies = n_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            selected_n_bins = selected_n_bins,&
            stone_failed = stone_failed,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            n_pooled_residuals = n_pooled_residuals,&
            min_bin_occupancy = min_bin_occupancy,&
            mean_bin_occupancy = mean_bin_occupancy,&
            max_bin_occupancy = max_bin_occupancy,&
            sturges_bins = sturges_bins,&
            fd_bins = fd_bins,&
            tmp_bin_counts = tmp_bin_counts,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            bin_count_mode = bin_count_mode_mode_f,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_bin_count_stone_expert_c

end module tox_data_integration_js_comp_test_binning_c
#endif
