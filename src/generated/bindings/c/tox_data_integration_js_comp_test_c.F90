#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
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
    public :: determine_bin_count_stone_c
    public :: determine_bin_count_stone_expert_c
    public :: generate_js_comp_test_candidates_c
    public :: generate_adaptive_js_comp_test_candidates_c
    public :: construct_adaptive_neighborhoods_c
    public :: construct_adaptive_neighborhoods_expert_c
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
    public :: run_js_comp_test_adaptive_c
    public :: run_js_comp_test_adaptive_expert_c
    public :: run_js_comp_test_parameter_search_c
    public :: run_js_comp_test_parameter_search_expert_c
    public :: run_js_comp_test_adaptive_parameter_search_c
    public :: run_js_comp_test_adaptive_parameter_search_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):estimate_bin_count(subroutine)]]
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] and returned as
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):estimate_bin_count_expert(subroutine)]]
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] and returned as
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy(subroutine)]]
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_expert(subroutine)]]
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_exhaustive(subroutine)]]
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]],
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_exhaustive_expert(subroutine)]]
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]],
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_stone(subroutine)]]
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
        use tox_data_integration_js_comp_test, only: determine_bin_count_stone
        use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

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
            !! | Mode                                                                                  | Value                                                                                           |
            !! |---------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Minimise the sum of the per-study scores (default)                                    | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled                                   | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):determine_bin_count_stone_expert(subroutine)]]
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
        use tox_data_integration_js_comp_test, only: determine_bin_count_stone_expert
        use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

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
            !! | Mode                                                                                  | Value                                                                                           |
            !! |---------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Minimise the sum of the per-study scores (default)                                    | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled                                   | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):generate_js_comp_test_candidates(subroutine)]]
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
        use tox_data_integration_js_comp_test, only: generate_js_comp_test_candidates

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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):generate_adaptive_js_comp_test_candidates(subroutine)]]
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
    !| The values strictly descend without deduplication. At least one candidate is always
    !| returned, even for data smaller than `k_start_1`, whose construction then fails with a
    !| too-few-means status in the search.
    !|
    !| Not checked here: the values fit 32-bit integers only while `k_start_1 <= huge(1_int32)/4`;
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]
    !| rejects inputs beyond that.
    subroutine generate_adaptive_js_comp_test_candidates_c(&
            max_n_genes_all_studies,&
            n_studies,&
            candidates_k_start_k_step_k_max,&
            n_candidates,&
            ierr&
        ) bind(C, name="generate_adaptive_js_comp_test_candidates_c")
        use tox_data_integration_js_comp_test, only: generate_adaptive_js_comp_test_candidates

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(3, 16), intent(out), target :: candidates_k_start_k_step_k_max
            !! Candidate `[k_start, k_step, k_max]` triples, `k_start` strictly descending
            !! The first `n_candidates` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates
            !! Number of candidate triples actually filled (at least 1, at most 16)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(n_candidates)
        M_CHECK_ARRAY_NON_NULL(candidates_k_start_k_step_k_max, 3 * 16)

        call generate_adaptive_js_comp_test_candidates(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            candidates_k_start_k_step_k_max = candidates_k_start_k_step_k_max,&
            n_candidates = n_candidates,&
            ierr = ierr&
        )
    end subroutine generate_adaptive_js_comp_test_candidates_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):construct_adaptive_neighborhoods(subroutine)]]
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
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !| (`10_int32`) non-NaN residuals the neighborhood stops with
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !| (`6_int32`); with all residuals zero, with
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !| (`5_int32`).
    !| 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    !| `(S_new - S_old) / S_old > tau` is discarded and growth stops
    !| ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TAU(variable)]], `1_int32`);
    !| otherwise it is committed. Growth also stops at `k_max` entries
    !| ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `2_int32`,
    !| checked first) or when no entry is left
    !| ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `4_int32`).
    !| 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    !| non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    !| is still evaluated, and growth then stops with
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`3_int32`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    !| sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    !| `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    !| absolute deviation of the neighborhood's pooled means. The construction ends with the
    !| neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`1_int32`) when fewer than `k_start` pooled means are non-NaN: nothing
    !| is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    !| and reports, in this priority,
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`3_int32`) if any neighborhood stopped with too few residuals,
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`2_int32`) if any has fewer than `min_study_neighbors`
    !| entries of some study, else
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`0_int32`).
    !|
    !| Outputs sized by `N` are filled only for their first `n_points` reference points.
    subroutine construct_adaptive_neighborhoods_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            k_start,&
            k_step,&
            k_max,&
            n_points,&
            x_star,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
            stop_reason,&
            neighborhood_dispersion,&
            neighborhood_mad,&
            max_n_neighbors,&
            construction_status,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            ierr&
        ) bind(C, name="construct_adaptive_neighborhoods_c")
        use tox_data_integration_js_comp_test, only: construct_adaptive_neighborhoods

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        integer(c_int), intent(out), target :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(2, max_n_genes_all_studies*n_studies), intent(out), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: stop_reason
            !! Why each neighborhood stopped growing: `1_int32` tau exceeded,
            !! `2_int32` `k_max` reached, `3_int32` residual cap
            !! reached, `4_int32` pool exhausted, `5_int32`
            !! zero dispersion, `6_int32` too few residuals
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), intent(out), target :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(c_int), intent(out), target :: construction_status
            !! Overall outcome: `0_int32` ok, `1_int32` too few
            !! means, `2_int32` a study neighborhood below
            !! `min_study_neighbors`, `3_int32` a neighborhood with too
            !! few residuals
        real(c_double), intent(in), target :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(c_int), intent(in), target :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(c_int), intent(in), target :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(k_step)
        M_CHECK_NON_NULL(k_max)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(max_n_neighbors)
        M_CHECK_NON_NULL(construction_status)
        M_CHECK_NON_NULL(tau)
        M_CHECK_NON_NULL(mad_distance_factor)
        M_CHECK_NON_NULL(max_pooled_residuals)
        M_CHECK_NON_NULL(min_study_neighbors)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(stop_reason, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_dispersion, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_mad, (max_n_genes_all_studies*n_studies))

        call construct_adaptive_neighborhoods(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            n_points = n_points,&
            x_star = x_star,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
            stop_reason = stop_reason,&
            neighborhood_dispersion = neighborhood_dispersion,&
            neighborhood_mad = neighborhood_mad,&
            max_n_neighbors = max_n_neighbors,&
            construction_status = construction_status,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            ierr = ierr&
        )
    end subroutine construct_adaptive_neighborhoods_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):construct_adaptive_neighborhoods_expert(subroutine)]]
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
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !| (`10_int32`) non-NaN residuals the neighborhood stops with
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !| (`6_int32`); with all residuals zero, with
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !| (`5_int32`).
    !| 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    !| `(S_new - S_old) / S_old > tau` is discarded and growth stops
    !| ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_TAU(variable)]], `1_int32`);
    !| otherwise it is committed. Growth also stops at `k_max` entries
    !| ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `2_int32`,
    !| checked first) or when no entry is left
    !| ([[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `4_int32`).
    !| 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    !| non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    !| is still evaluated, and growth then stops with
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`3_int32`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    !| sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    !| `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    !| absolute deviation of the neighborhood's pooled means. The construction ends with the
    !| neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`1_int32`) when fewer than `k_start` pooled means are non-NaN: nothing
    !| is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    !| and reports, in this priority,
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`3_int32`) if any neighborhood stopped with too few residuals,
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`2_int32`) if any has fewer than `min_study_neighbors`
    !| entries of some study, else
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`0_int32`).
    !|
    !| Outputs sized by `N` are filled only for their first `n_points` reference points.
    subroutine construct_adaptive_neighborhoods_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            k_start,&
            k_step,&
            k_max,&
            n_points,&
            x_star,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
            stop_reason,&
            neighborhood_dispersion,&
            neighborhood_mad,&
            max_n_neighbors,&
            construction_status,&
            tmp_gene_means_perm_all,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            ierr&
        ) bind(C, name="construct_adaptive_neighborhoods_expert_c")
        use tox_data_integration_js_comp_test, only: construct_adaptive_neighborhoods_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        integer(c_int), intent(out), target :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(2, max_n_genes_all_studies*n_studies), intent(out), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: stop_reason
            !! Why each neighborhood stopped growing: `1_int32` tau exceeded,
            !! `2_int32` `k_max` reached, `3_int32` residual cap
            !! reached, `4_int32` pool exhausted, `5_int32`
            !! zero dispersion, `6_int32` too few residuals
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), intent(out), target :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(c_int), intent(out), target :: construction_status
            !! Overall outcome: `0_int32` ok, `1_int32` too few
            !! means, `2_int32` a study neighborhood below
            !! `min_study_neighbors`, `3_int32` a neighborhood with too
            !! few residuals
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, seeded and sorted here
        real(c_double), intent(in), target :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(c_int), intent(in), target :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(c_int), intent(in), target :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(k_step)
        M_CHECK_NON_NULL(k_max)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(max_n_neighbors)
        M_CHECK_NON_NULL(construction_status)
        M_CHECK_NON_NULL(tau)
        M_CHECK_NON_NULL(mad_distance_factor)
        M_CHECK_NON_NULL(max_pooled_residuals)
        M_CHECK_NON_NULL(min_study_neighbors)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(stop_reason, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_dispersion, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_mad, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))

        call construct_adaptive_neighborhoods_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            n_points = n_points,&
            x_star = x_star,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
            stop_reason = stop_reason,&
            neighborhood_dispersion = neighborhood_dispersion,&
            neighborhood_mad = neighborhood_mad,&
            max_n_neighbors = max_n_neighbors,&
            construction_status = construction_status,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            ierr = ierr&
        )
    end subroutine construct_adaptive_neighborhoods_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_neighborhood_overlaps(subroutine)]]
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
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
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
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
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]].
    !| Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    !| them does not count as a failure.
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
            !! Second extent of `mean_pmf_counts`: the widest bin count of any point
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: n_bins_per_point
            !! Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
            !! is inspected
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
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
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
    !| For each study, the overlap of the new candidate's interval with the best one so far is
    !| measured as a fraction of the new interval's width, and counted when it reaches
    !| `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    !| best (the best is then kept), or, after the new candidate became the best, when the count
    !| meets the condition `join_method` names.
    subroutine check_plateau_condition_c(&
            confidence_interval,&
            best_candidate_pair_confidence_interval,&
            n_studies,&
            best_candidate_index,&
            best_exceeded_ci_overlap_count,&
            candidate_index,&
            join_method,&
            succeeding_ci_overlap,&
            ci_overlap,&
            plateau_found,&
            ierr&
        ) bind(C, name="check_plateau_condition_c")
        use tox_data_integration_js_comp_test, only: check_plateau_condition
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(2, n_studies), intent(in), target :: confidence_interval
            !! Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
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
            !! How the studies' overlaps are joined into the plateau condition; the median method
            !! requires `count > (n_studies - 1) / 2`
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(c_double), intent(in), target :: succeeding_ci_overlap
            !! Minimum fractional overlap with the best candidate's interval for a study to pass
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_studies), intent(out), target :: ci_overlap
            !! Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            !! it was on entry, relative to the new interval's width; always filled. Meaningless
            !! for the very first candidate, whose best interval is only a seed
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
        M_CHECK_ARRAY_NON_NULL(ci_overlap, n_studies)

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
            ci_overlap = ci_overlap,&
            plateau_found = plateau_found,&
            ierr = ierr&
        )
    end subroutine check_plateau_condition_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition(subroutine)]]
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
            !! `.false.` for the first admissible candidate, which has no transition
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
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
            !! `.false.` for the first admissible candidate, which has no transition
        real(c_double), intent(in), target :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(c_double), intent(in), target :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(c_int), intent(in), target :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(inout), target :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
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
    !| `mean_pmf` is the equal-weight mean of all studies' pmfs, `mean_pmf_counts` and
    !| `mean_pmf_included_n_reps` the sums of their counts. Known limitation: every study is
    !| compared against a consensus that includes itself, not a leave-one-out background.
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
            !! Second extent of the bin-sized arrays: the widest bin count of any reference point.
            !! All studies must share each point's bin count, so their zero padding coincides
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_points, n_bins, n_studies), intent(in), target :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_points, n_bins, n_studies), intent(in), target :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study; every study needs at
            !! least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_points, n_bins), intent(out), target :: mean_pmf
            !! The consensus pmf, `sum(pmfs, dim=3) / n_studies`
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point, summed over the studies
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus, `sum(counts, dim=3)`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)

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
    !| As [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]], for
    !| callers that need no counts (the bootstrap). Same known limitation: the consensus includes
    !| the study compared against it.
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
            !! Second extent of the bin-sized arrays: the widest bin count of any reference point.
            !! All studies must share each point's bin count, so their zero padding coincides
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(n_points, n_bins, n_studies), intent(in), target :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_points, n_bins), intent(out), target :: mean_pmf
            !! The consensus pmf, `sum(pmfs, dim=3) / n_studies`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * n_bins)

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
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
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
            tmp_point_pool,&
            tmp_point_draw,&
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
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array: per study the bottom-k heap `(:, 1, :)` and the top-k heap `(:, 2, :)`
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Working array: one study's resampled histogram counts
        real(c_double), dimension(n_points, n_bins, n_studies), intent(out), target :: tmp_pmfs
            !! Working array: one bootstrap's resampled pmfs, all studies
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_mean_pmf
            !! Working array: one bootstrap's consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array: one bootstrap's per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_weights
            !! Working array: one bootstrap's per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array: one bootstrap's weighted global JSD values
        integer(c_int), dimension(n_bins), intent(out), target :: tmp_point_pool
            !! Working array: one reference point's pooled counts, the pool of one draw
        integer(c_int), dimension(n_bins), intent(out), target :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, n_bins)

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
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test(subroutine)]]
    !| The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
    !|
    !| 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    !| the caller, e.g. with
    !| [[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) form its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| 2. Per reference point, all studies' residuals are pooled and the bin-count search
    !| `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| by default, or Stone's rule
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]).
    !| 3. Each study's histogram and pmf is built per point
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), each
    !| study's per-point JSD against it and its weighted global JSD, and the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]) give
    !| `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
    !| `n_permutations = 0` both p-values are NaN.
    !|
    !| There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
    !| still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
    !| caller that needs such neighborhoods rejected should use
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !|
    !| The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
    !| second extent [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]]
    !| (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
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
            p_values_observed_consensus,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test
        use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

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
            !! The maximum valid value is `max_n_genes_all_studies`.
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
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's histogram bin count `M_j`, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts` (see above)
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes (see above)
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(c_double), dimension(n_studies), intent(out), target :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_occupancy'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32) :: bin_count_mode_mode_f

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
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values_observed_consensus, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

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
            p_values_observed_consensus = p_values_observed_consensus,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_expert(subroutine)]]
    !| The final JS-Comp-Test of K studies for one fixed `(n_points, n_neighbors)` setting:
    !|
    !| 1. Each study's `n_neighbors` genes nearest every reference point `x_star` (computed once by
    !| the caller, e.g. with
    !| [[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) form its
    !| neighborhoods
    !| ([[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| 2. Per reference point, all studies' residuals are pooled and the bin-count search
    !| `bin_count_mode` names chooses the point's bin count and residual range (occupancy search
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| by default, or Stone's rule
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]).
    !| 3. Each study's histogram and pmf is built per point
    !| ([[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]).
    !| 4. The consensus pmf
    !| ([[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]]), each
    !| study's per-point JSD against it and its weighted global JSD, and the permutation test
    !| ([[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]]) give
    !| `p_values` (primary) and `p_values_observed_consensus` (comparison only). With
    !| `n_permutations = 0` both p-values are NaN.
    !|
    !| There is no admissibility gate: a point whose bin-count search fails (`occupancy_failed`)
    !| still gets a histogram at `m_min` bins and still contributes to `global_js_divergence`. A
    !| caller that needs such neighborhoods rejected should use
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]].
    !|
    !| The bin-sized arrays (`pmfs`, `counts`, `mean_pmf`, `mean_pmf_counts`) have the fixed
    !| second extent [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]]
    !| (`256`); only their first `max_n_bins_per_point` bins are meaningful, the rest is
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
            p_values_observed_consensus,&
            tmp_neighborhood_residuals_gathered,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool,&
            tmp_permutation_point_draw,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
            tmp_point_n_neighbors,&
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
            bin_count_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_expert
        use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

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
            !! The maximum valid value is `max_n_genes_all_studies`.
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
        integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_indices
            !! Gene indices of the selected neighborhood, per reference point, per study
        integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_range
            !! For each reference point and study, the `[min_idx, max_idx]` neighborhood span
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's histogram bin count `M_j`, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts` (see above)
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes (see above)
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(c_double), dimension(n_studies), intent(out), target :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        real(c_double), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out), target :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        real(c_double), dimension(n_points, 256), intent(out), target :: tmp_permutation_mean_pmf
            !! Working array of the permutation test: consensus pmf of the permuted studies
        integer(c_int), dimension(256), intent(out), target :: tmp_permutation_point_pool
            !! Working array of the permutation test: one point's pool
        integer(c_int), dimension(256), intent(out), target :: tmp_permutation_point_draw
            !! Working array of the permutation test: one study's draw at one point
        real(c_double), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array: the last index of each study's block in tmp_pooled_residuals
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in tmp_pooled_residuals
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_point_n_neighbors
            !! Working array: every study's neighbor count at a reference point (all `n_neighbors`)
        integer(c_int), dimension(n_points, 256), intent(out), target :: tmp_permutation_mean_pmf_counts
            !! Working array of the permutation test: the remaining pool
        integer(c_int), dimension(n_points, 256), intent(out), target :: tmp_permutation_counts
            !! Working array of the permutation test: one study's resampled counts
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: tmp_permutation_pmfs
            !! Working array of the permutation test: resampled pmfs
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_js_divergences
            !! Working array of the permutation test: per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_weights
            !! Working array of the permutation test: per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_permutation_global_js_divergence
            !! Working array of the permutation test: global JSD values
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_occupancy'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32) :: bin_count_mode_mode_f

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
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values_observed_consensus, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_residuals_gathered, max_n_reps_all_studies * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_point_pool, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_point_draw, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*n_neighbors*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_n_neighbors, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

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
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_neighborhood_residuals_gathered = tmp_neighborhood_residuals_gathered,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool = tmp_permutation_point_pool,&
            tmp_permutation_point_draw = tmp_permutation_point_draw,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
            tmp_point_n_neighbors = tmp_point_n_neighbors,&
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
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive(subroutine)]]
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
    subroutine run_js_comp_test_adaptive_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            gene_means,&
            residuals,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
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
            p_values_observed_consensus,&
            n_permutations,&
            random_seed,&
            min_residuals_per_bin,&
            m_min,&
            m_max,&
            gamma_occupancy,&
            lower_residual_range_quantile,&
            upper_residual_range_quantile,&
            bin_count_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_adaptive_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_adaptive
        use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

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
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), dimension(2, n_points), intent(in), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
        integer(c_int), dimension(n_studies, n_points), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(c_double), dimension(n_studies), intent(out), target :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_occupancy'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32) :: bin_count_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
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
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * n_points)
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * n_points)
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
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values_observed_consensus, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

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

        call run_js_comp_test_adaptive(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            gene_means = gene_means,&
            residuals = residuals,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
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
            p_values_observed_consensus = p_values_observed_consensus,&
            n_permutations = n_permutations,&
            random_seed = random_seed,&
            min_residuals_per_bin = min_residuals_per_bin,&
            m_min = m_min,&
            m_max = m_max,&
            gamma_occupancy = gamma_occupancy,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile,&
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_adaptive_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive_expert(subroutine)]]
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
    subroutine run_js_comp_test_adaptive_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            n_points,&
            gene_means,&
            residuals,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
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
            p_values_observed_consensus,&
            tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals,&
            tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool,&
            tmp_permutation_point_draw,&
            tmp_pooled_residuals,&
            tmp_pooled_residuals_perm,&
            tmp_bin_counts_search,&
            tmp_residual_to_study_perm,&
            tmp_study_block_ends,&
            tmp_study_block_sizes,&
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
            bin_count_mode,&
            ierr&
        ) bind(C, name="run_js_comp_test_adaptive_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_adaptive_expert
        use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED

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
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), dimension(2, n_points), intent(in), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
        integer(c_int), dimension(n_studies, n_points), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its genes belong to each study (at least 1 each)
        integer(c_int), dimension(n_points), intent(out), target :: n_bins_per_point
            !! Each reference point's selected histogram bin count, from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_low
            !! Each reference point's lower residual-range bound (R_low), from the bin-count search
        real(c_double), dimension(n_points), intent(out), target :: shared_residual_range_high
            !! Each reference point's upper residual-range bound (R_high), from the bin-count search
        integer(c_int), intent(out), target :: max_n_bins_per_point
            !! `maxval(n_bins_per_point)`: the number of meaningful leading bins of `pmfs`, `counts`,
            !! `mean_pmf` and `mean_pmf_counts`
        logical(c_bool), dimension(n_points), intent(out), target :: occupancy_failed
            !! `.true.` iff the bin-count search failed for this point; it then uses `m_min` bins and
            !! still contributes
        integer(c_int), dimension(n_points), intent(out), target :: n_pooled_residuals
            !! Each reference point's non-NaN pooled residual count (N_j)
        integer(c_int), dimension(n_points), intent(out), target :: min_bin_occupancy
            !! Each reference point's minimum bin occupancy at `n_bins_per_point`
        real(c_double), dimension(n_points), intent(out), target :: mean_bin_occupancy
            !! Each reference point's mean bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: max_bin_occupancy
            !! Each reference point's maximum bin occupancy at `n_bins_per_point`
        integer(c_int), dimension(n_points), intent(out), target :: sturges_bins
            !! Each reference point's Sturges' rule diagnostic (not used for the decision)
        integer(c_int), dimension(n_points), intent(out), target :: fd_bins
            !! Each reference point's Freedman-Diaconis rule diagnostic (not used for the decision)
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(c_int), dimension(n_points, 256, n_studies), intent(out), target :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(c_double), dimension(n_points, 256), intent(out), target :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(c_int), dimension(n_points, 256), intent(out), target :: mean_pmf_counts
            !! Residual counts per bin of the consensus pmf; only bins `1:max_n_bins_per_point` are
            !! meaningful
        integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
            !! Count of non-NaN residuals per reference point for the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
            !! Per-reference-point JSD of each study against the consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
            !! Per-reference-point weights for `global_js_divergence`
        real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
            !! Weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, the primary result; NaN when `n_permutations` is 0
        real(c_double), dimension(n_studies), intent(out), target :: p_values_observed_consensus
            !! Empirical p-value per study against the fixed observed consensus `mean_pmf`, kept for
            !! comparison only; NaN when `n_permutations` is 0
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out), target :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        real(c_double), dimension(n_points, 256), intent(out), target :: tmp_permutation_mean_pmf
            !! Working array forwarded to the permutation test
        integer(c_int), dimension(256), intent(out), target :: tmp_permutation_point_pool
            !! Working array forwarded to the permutation test
        integer(c_int), dimension(256), intent(out), target :: tmp_permutation_point_draw
            !! Working array forwarded to the permutation test
        real(c_double), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        integer(c_int), dimension(n_points, 256), intent(out), target :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to the permutation test
        integer(c_int), dimension(n_points, 256), intent(out), target :: tmp_permutation_counts
            !! Working array forwarded to the permutation test
        real(c_double), dimension(n_points, 256, n_studies), intent(out), target :: tmp_permutation_pmfs
            !! Working array forwarded to the permutation test
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_js_divergences
            !! Working array forwarded to the permutation test
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_permutation_weights
            !! Working array forwarded to the permutation test
        real(c_double), dimension(n_studies), intent(out), target :: tmp_permutation_global_js_divergence
            !! Working array forwarded to the permutation test
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `'bin_count_occupancy'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32) :: bin_count_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_points)
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
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * n_points)
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * n_points)
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
        M_CHECK_ARRAY_NON_NULL(pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values_observed_consensus, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_point_neighborhood_indices, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbor_residuals, max_n_reps_all_studies * max_n_genes_all_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_mean_pmf, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_point_pool, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_point_draw, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_ends, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_study_block_sizes, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_mean_pmf_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_counts, n_points * 256)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_pmfs, n_points * 256 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(bin_count_mode, 29)

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

        call run_js_comp_test_adaptive_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_points = n_points,&
            gene_means = gene_means,&
            residuals = residuals,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
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
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tmp_point_neighborhood_indices = tmp_point_neighborhood_indices,&
            tmp_neighbor_residuals = tmp_neighbor_residuals,&
            tmp_permutation_mean_pmf = tmp_permutation_mean_pmf,&
            tmp_permutation_point_pool = tmp_permutation_point_pool,&
            tmp_permutation_point_draw = tmp_permutation_point_draw,&
            tmp_pooled_residuals = tmp_pooled_residuals,&
            tmp_pooled_residuals_perm = tmp_pooled_residuals_perm,&
            tmp_bin_counts_search = tmp_bin_counts_search,&
            tmp_residual_to_study_perm = tmp_residual_to_study_perm,&
            tmp_study_block_ends = tmp_study_block_ends,&
            tmp_study_block_sizes = tmp_study_block_sizes,&
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
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_adaptive_expert_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search(subroutine)]]
    !| Walks the candidates of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !| study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !| residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !| `min_residuals_per_bin`
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !| ([[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]),
    !| and both plateau criteria: CI overlap with the best candidate so far
    !| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]])
    !| and relative effect-size stability
    !| ([[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
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
        use tox_data_integration_js_comp_test, only: run_js_comp_test_parameter_search
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

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
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
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
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search_expert(subroutine)]]
    !| Walks the candidates of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| from finest to coarsest resolution. For each candidate:
    !|
    !| 1. Reference points from the pooled gene means
    !| ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]) and every
    !| study's neighborhoods; gate 1: consecutive neighborhoods overlap enough
    !| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]).
    !| 2. Per reference point, the bin-count search `bin_count_mode` names on all studies' pooled
    !| residuals, then every study's histograms.
    !| 3. Gate 2: every study has a non-NaN residual at every point, and every consensus bin reaches
    !| `min_residuals_per_bin`
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| A candidate that passes both gates is *admissible*.
    !| 4. The observed global JSD per study, its bootstrap confidence interval
    !| ([[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]),
    !| and both plateau criteria: CI overlap with the best candidate so far
    !| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]])
    !| and relative effect-size stability
    !| ([[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]]).
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
        use tox_data_integration_js_comp_test, only: run_js_comp_test_parameter_search_expert
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

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
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
        integer(c_int), intent(in), target :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
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
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive_parameter_search(subroutine)]]
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
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    !| `min_study_neighbors`.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !| reference points emerged.
    !| 3. Gate 1:
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !| on the pooled ranges.
    !| 4. The pipeline of
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !| up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    !| (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    !| count
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    !| first plateau `plateau_mode` selects.
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
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]);
    !| the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
    !| `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
    !| `max_n_reps_all_studies * N` values each, the size of `residuals`.
    !|
    !| A failure to create the bootstrap's random number generator is returned in `ierr` and
    !| stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_adaptive_parameter_search_c(&
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
        ) bind(C, name="run_js_comp_test_adaptive_parameter_search_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_adaptive_parameter_search
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

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
            !! Most reference points a candidate may emerge with; a candidate with more is rejected
            !! with the capacity status. The per-point outputs and traces below are sized by it, so
            !! the caller has to know it up front
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds]].
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps to perform for a candidate
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(out), target :: k_start
            !! The finally chosen candidate's `k_start`
        integer(c_int), intent(out), target :: k_step
            !! The finally chosen candidate's `k_step`
        integer(c_int), intent(out), target :: k_max
            !! The finally chosen candidate's `k_max`
        integer(c_int), intent(out), target :: n_points
            !! Number of reference points the finally chosen candidate's construction produced
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful. Without an admissible candidate, the first
            !! `min(n_points, max_n_points_candidate)` entries are `m_min`
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` without an admissible candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_confidence_interval
            !! The bootstrapped JSD confidence interval of the finally chosen candidate;
            !! `-1.0_real64` throughout where no admissible candidate backs the result
        logical(c_bool), intent(out), target :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array
        integer(c_int), dimension(16), intent(out), target :: trace_k_start
            !! Per-admissible-candidate `k_start`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width; `-1.0` for the
            !! first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0` for the first
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
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            !! as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(c_int), dimension(16), intent(out), target :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! The first `n_candidates_tried` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: candidate_n_points
            !! Number of reference points each tried candidate's construction produced (0 for too
            !! few means)
            !! The first `n_candidates_tried` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: candidate_status
            !! What happened to every tried candidate: `0_int32` evaluated,
            !! `1_int32` construction failed,
            !! `2_int32` a study neighborhood below `min_study_neighbors`
            !! or a study without a non-NaN residual at some point,
            !! `3_int32` capacity exceeded,
            !! `4_int32` overlap gate failed,
            !! `5_int32` consensus-pmf count gate failed
            !! The first `n_candidates_tried` elements will hold the results.
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum count every bin must reach, both in the bin-count search (the floor in every
            !! `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
        real(c_double), intent(in), target :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(c_int), intent(in), target :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(c_int), intent(in), target :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(k_step)
        M_CHECK_NON_NULL(k_max)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(plateau_established)
        M_CHECK_NON_NULL(n_admissible_evaluated)
        M_CHECK_NON_NULL(n_candidates_tried)
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
        M_CHECK_NON_NULL(tau)
        M_CHECK_NON_NULL(mad_distance_factor)
        M_CHECK_NON_NULL(max_pooled_residuals)
        M_CHECK_NON_NULL(min_study_neighbors)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(best_candidate_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(trace_k_start, 16)
        M_CHECK_ARRAY_NON_NULL(trace_k_step, 16)
        M_CHECK_ARRAY_NON_NULL(trace_k_max, 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_points, 16)
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
        M_CHECK_ARRAY_NON_NULL(candidate_k_start, 16)
        M_CHECK_ARRAY_NON_NULL(candidate_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(candidate_status, 16)
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

        call run_js_comp_test_adaptive_parameter_search(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
            max_n_points_candidate = max_n_points_candidate,&
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
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_adaptive_parameter_search_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive_parameter_search_expert(subroutine)]]
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
    !| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
    !| does with the same `tau`, `mad_distance_factor`, `max_pooled_residuals` and
    !| `min_study_neighbors`.
    !| 2. Reject it if the construction status is not ok, or if more than `max_n_points_candidate`
    !| reference points emerged.
    !| 3. Gate 1:
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]]
    !| on the pooled ranges.
    !| 4. The pipeline of
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]]
    !| up to the consensus pmf, then gate 2: every study has a non-NaN residual at every point
    !| (otherwise logged as an empty-study candidate), and every consensus bin reaches the minimum
    !| count
    !| ([[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]).
    !| 5. As in the fixed-k search: traces, bootstrap, both plateau criteria, and the stop at the
    !| first plateau `plateau_mode` selects.
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
    !| [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(interface)]]);
    !| the construction buffers by `N = max_n_genes_all_studies * n_studies` (the largest,
    !| `tmp_n_neighbors_per_point`, `n_studies * N` integers); the pooling buffers hold
    !| `max_n_reps_all_studies * N` values each, the size of `residuals`.
    !|
    !| A failure to create the bootstrap's random number generator is returned in `ierr` and
    !| stops the search; the outputs are then undefined.
    subroutine run_js_comp_test_adaptive_parameter_search_expert_c(&
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
        ) bind(C, name="run_js_comp_test_adaptive_parameter_search_expert_c")
        use tox_data_integration_js_comp_test, only: run_js_comp_test_adaptive_parameter_search_expert
        use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, MODE_PLATEAU_BOTH, MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE

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
            !! Most reference points a candidate may emerge with; a candidate with more is rejected
            !! with the capacity status. The per-point outputs and traces below are sized by it, so
            !! the caller has to know it up front
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps to perform for a candidate
            !! The minimum valid value is `1_int32`.
        character(len=1, kind=c_char), dimension(11), intent(in), target :: join_method
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(c_int), intent(out), target :: k_start
            !! The finally chosen candidate's `k_start`
        integer(c_int), intent(out), target :: k_step
            !! The finally chosen candidate's `k_step`
        integer(c_int), intent(out), target :: k_max
            !! The finally chosen candidate's `k_max`
        integer(c_int), intent(out), target :: n_points
            !! Number of reference points the finally chosen candidate's construction produced
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: n_bins_per_point
            !! The chosen candidate's bin count per reference point; the first `n_points` entries
            !! are meaningful. Without an admissible candidate, the first
            !! `min(n_points, max_n_points_candidate)` entries are `m_min`
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_low
            !! The chosen candidate's lower residual-range bound (R_low) per reference point; the
            !! first `n_points` entries are meaningful, `0.0` without an admissible candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: shared_residual_range_high
            !! The chosen candidate's upper residual-range bound (R_high), like `shared_residual_range_low`
        real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_confidence_interval
            !! The bootstrapped JSD confidence interval of the finally chosen candidate;
            !! `-1.0_real64` throughout where no admissible candidate backs the result
        logical(c_bool), intent(out), target :: plateau_established
            !! `.true.` when a plateau was found (by the criterion `plateau_mode` selects) or the
            !! only candidate was admissible; `.false.` otherwise
        integer(c_int), intent(out), target :: n_admissible_evaluated
            !! Number of admissible candidates evaluated before the search stopped: the number of
            !! valid columns of every `trace_*` array
        integer(c_int), dimension(16), intent(out), target :: trace_k_start
            !! Per-admissible-candidate `k_start`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_k_step
            !! Per-admissible-candidate `k_step`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_k_max
            !! Per-admissible-candidate `k_max`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: trace_n_points
            !! Per-admissible-candidate number of reference points
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_lower
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval lower bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_upper
            !! Per-admissible-candidate, per-study bootstrapped confidence-interval upper bound
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width
            !! Per-admissible-candidate, per-study confidence-interval width
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_width_relative
            !! Per-admissible-candidate, per-study relative confidence-interval width
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width; `-1.0` for the
            !! first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(n_studies, 16), intent(out), target :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate; `-1.0` for the first
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
            !! Per-admissible-candidate, per-reference-point pooled non-NaN residual count (`N_j`);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic; jagged as
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), dimension(max_n_points_candidate, 16), intent(out), target :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic; jagged
            !! as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(c_double), dimension(max_n_points_candidate, 16), intent(out), target :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high);
            !! jagged as trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(c_int), intent(out), target :: n_candidates_tried
            !! Number of candidates the search tried before it stopped: the number of leading,
            !! valid entries of `candidate_k_start`, `candidate_n_points` and `candidate_status`
        integer(c_int), dimension(16), intent(out), target :: candidate_k_start
            !! `k_start` of every candidate tried, in sequence order
            !! The first `n_candidates_tried` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: candidate_n_points
            !! Number of reference points each tried candidate's construction produced (0 for too
            !! few means)
            !! The first `n_candidates_tried` elements will hold the results.
        integer(c_int), dimension(16), intent(out), target :: candidate_status
            !! What happened to every tried candidate: `0_int32` evaluated,
            !! `1_int32` construction failed,
            !! `2_int32` a study neighborhood below `min_study_neighbors`
            !! or a study without a non-NaN residual at some point,
            !! `3_int32` capacity exceeded,
            !! `4_int32` overlap gate failed,
            !! `5_int32` consensus-pmf count gate failed
            !! The first `n_candidates_tried` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, as in the construction
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_x_star
            !! Working array: the current candidate's reference points
        integer(c_int), dimension(2, max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_neighborhood_range
            !! Working array: the current candidate's `[first, last]` pooled position per neighborhood
        integer(c_int), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out), target :: tmp_n_neighbors_per_point
            !! Working array: the current candidate's per-study entry count per neighborhood
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_stop_reason
            !! Working array: the current candidate's stop reason per neighborhood
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_neighborhood_dispersion
            !! Working array: the current candidate's dispersion per neighborhood
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_neighborhood_mad
            !! Working array: the current candidate's MAD of the pooled means per neighborhood
        integer(c_int), dimension(max_n_genes_all_studies, n_studies), intent(out), target :: tmp_point_neighborhood_indices
            !! Working array: one reference point's gene indices per study, reused per point
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out), target :: tmp_neighbor_residuals
            !! Working array: one (reference point, study) pair's gathered residuals, reused per pair
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_bins_per_point
            !! Working array: the current candidate's per-point bin count
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_low
            !! Working array: the current candidate's per-point R_low
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_shared_residual_range_high
            !! Working array: the current candidate's per-point R_high
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
        integer(c_int), dimension(256), intent(out), target :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(c_int), dimension(256), intent(out), target :: tmp_point_draw
            !! Working array of the bootstrap: one study's draw at one point
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
        real(c_double), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for `tmp_pooled_residuals`, reused per point
        integer(c_int), dimension(256), intent(out), target :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(c_int), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out), target :: tmp_residual_to_study_perm
            !! Working array (Stone modes): pool indices grouped by study, each study sorted, reused
            !! per point
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_ends
            !! Working array: the last index of each study's block in `tmp_pooled_residuals`
        integer(c_int), dimension(n_studies), intent(out), target :: tmp_study_block_sizes
            !! Working array (Stone modes): each study's block size in `tmp_pooled_residuals`
        logical(c_bool), dimension(max_n_points_candidate), intent(out), target :: tmp_occupancy_failed
            !! Working array: the current candidate's per-point occupancy_failed flag
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_n_pooled_residuals
            !! Working array: the current candidate's per-point pooled residual count
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_min_bin_occupancy
            !! Working array: the current candidate's per-point minimum bin occupancy
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_mean_bin_occupancy
            !! Working array: the current candidate's per-point mean bin occupancy
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_max_bin_occupancy
            !! Working array: the current candidate's per-point maximum bin occupancy
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_sturges_bins
            !! Working array: the current candidate's per-point Sturges' rule diagnostic
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_fd_bins
            !! Working array: the current candidate's per-point Freedman-Diaconis rule diagnostic
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_n_bins_per_point
            !! Working array: snapshot of tmp_n_bins_per_point for the current best candidate
        integer(c_int), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_n_bins_per_point
            !! Working array: snapshot of tmp_n_bins_per_point for the smallest-uncertainty candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_low
            !! Working array: snapshot of tmp_shared_residual_range_low for the current best candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_shared_residual_range_high
            !! Working array: snapshot of tmp_shared_residual_range_high for the current best candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array: R_low per point of the smallest-uncertainty admissible candidate
        real(c_double), dimension(max_n_points_candidate), intent(out), target :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array: R_high per point of the smallest-uncertainty admissible candidate
        integer(c_int), intent(in), target :: min_residuals_per_bin
            !! Minimum count every bin must reach, both in the bin-count search (the floor in every
            !! `bin_count_mode`) and in the second admissibility gate (the consensus pmf)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods' pooled ranges (gate 1)
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
        real(c_double), intent(in), target :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(c_int), intent(in), target :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(c_int), intent(in), target :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for a candidate to be
            !! admissible
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        character(len=1, kind=c_char), dimension(29), intent(in), target :: bin_count_mode
            !! How each reference point's histogram bin count is chosen, for every candidate; the
            !! Stone modes test every bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(k_step)
        M_CHECK_NON_NULL(k_max)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(plateau_established)
        M_CHECK_NON_NULL(n_admissible_evaluated)
        M_CHECK_NON_NULL(n_candidates_tried)
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
        M_CHECK_NON_NULL(tau)
        M_CHECK_NON_NULL(mad_distance_factor)
        M_CHECK_NON_NULL(max_pooled_residuals)
        M_CHECK_NON_NULL(min_study_neighbors)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(join_method, 11)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_low, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(shared_residual_range_high, max_n_points_candidate)
        M_CHECK_ARRAY_NON_NULL(best_candidate_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(trace_k_start, 16)
        M_CHECK_ARRAY_NON_NULL(trace_k_step, 16)
        M_CHECK_ARRAY_NON_NULL(trace_k_max, 16)
        M_CHECK_ARRAY_NON_NULL(trace_n_points, 16)
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
        M_CHECK_ARRAY_NON_NULL(candidate_k_start, 16)
        M_CHECK_ARRAY_NON_NULL(candidate_n_points, 16)
        M_CHECK_ARRAY_NON_NULL(candidate_status, 16)
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_x_star, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_neighborhood_range, 2 * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_n_neighbors_per_point, n_studies * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_stop_reason, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_dispersion, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_neighborhood_mad, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_point_neighborhood_indices, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_neighbor_residuals, max_n_reps_all_studies * max_n_genes_all_studies)
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
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_prev_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_delta_perm, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_best_uncertainty_confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_residuals_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_bin_counts_search, 256)
        M_CHECK_ARRAY_NON_NULL(tmp_residual_to_study_perm, (max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
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

        call run_js_comp_test_adaptive_parameter_search_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            n_bootstraps = n_bootstraps,&
            join_method = join_method_mode_f,&
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
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            bin_count_mode = bin_count_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine run_js_comp_test_adaptive_parameter_search_expert_c

end module tox_data_integration_js_comp_test_c
#endif
