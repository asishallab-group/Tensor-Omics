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
!|
!| Generated from [[tox_data_integration_js_comp_test_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test
    use f42_safeguard
    use tox_data_integration_js_comp_test_impl, only: KX_FACTORS, MAX_N_BINS, MAX_POINTS, METHOD_JOIN_MAX
    use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MEDIAN, METHOD_JOIN_MIN, MIN_POINTS, MODE_BIN_COUNT_OCCUPANCY
    use tox_data_integration_js_comp_test_impl, only: MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_POOLED, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, MODE_PLATEAU_BOTH
    use tox_data_integration_js_comp_test_impl, only: MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, bootstrap_histogram_impl, calc_js_comp_test_n_top_k_jsds
    use tox_data_integration_js_comp_test_impl, only: check_effect_size_plateau_condition_impl, check_mean_pmf_min_counts_impl, check_neighborhood_overlaps_impl, check_plateau_condition_impl
    use tox_data_integration_js_comp_test_impl, only: construct_adaptive_neighborhoods_impl, create_mean_pmf_impl, create_mean_pmf_only_impl, determine_bin_count_occupancy_exhaustive_impl
    use tox_data_integration_js_comp_test_impl, only: determine_bin_count_occupancy_impl, determine_bin_count_stone_impl, estimate_bin_count_impl, generate_adaptive_js_comp_test_candidates_impl
    use tox_data_integration_js_comp_test_impl, only: generate_js_comp_test_candidates_impl, run_js_comp_test_adaptive_impl, run_js_comp_test_adaptive_parameter_search_impl, run_js_comp_test_impl
    use tox_data_integration_js_comp_test_impl, only: run_js_comp_test_parameter_search_impl
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
    public :: determine_bin_count_occupancy_exhaustive
    public :: determine_bin_count_occupancy_exhaustive_expert
    public :: determine_bin_count_stone
    public :: determine_bin_count_stone_expert
    public :: generate_js_comp_test_candidates
    public :: generate_adaptive_js_comp_test_candidates
    public :: construct_adaptive_neighborhoods
    public :: construct_adaptive_neighborhoods_expert
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
    public :: run_js_comp_test_adaptive
    public :: run_js_comp_test_adaptive_expert
    public :: run_js_comp_test_parameter_search
    public :: run_js_comp_test_parameter_search_expert
    public :: run_js_comp_test_adaptive_parameter_search
    public :: run_js_comp_test_adaptive_parameter_search_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):estimate_bin_count_expert]] to prepare it yourself.
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] and returned as
    !| `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    !| `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    !| from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    !| bin count.
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
            !! to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN
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
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_impl(module):MAX_N_BINS(variable)]] and returned as
    !| `sturges_bins`/`fd_bins`; `n_bins` is their maximum. Both take the sample size to be
    !| `max_n_reps_all_studies*n_neighbors` (one study's padded slot count), so they can differ
    !| from the textbook rules. The JS-Comp-Test uses them as diagnostics only, never to decide a
    !| bin count.
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
            !! to sturges_bins when the interquartile range is (close to) zero or the ratio is NaN
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
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
            !! smaller
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
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
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=20_int32, min=lower_residual_range_quantile, max=1.0_real64)
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
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts
            !! Working array: per-bin counts of the candidate bin count under test
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
            !! smaller
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
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
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=22_int32, min=lower_residual_range_quantile, max=1.0_real64)
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

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_exhaustive_expert]] to prepare it yourself.
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]],
    !| whose faster search could in principle miss a larger admissible `M`. Takes the range and
    !| `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
    pure subroutine determine_bin_count_occupancy_exhaustive(&
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
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), dimension(n_residuals), intent(in) :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_residuals`.
        real(real64), intent(in) :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(real64), intent(in) :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! The minimum valid value is `shared_residual_range_low`.
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
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: pooled_residuals_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_pooled_residuals, ierr, arg_pos=3_int32, min=0_int32, max=n_residuals)
        call validate_in_range_real(shared_residual_range_low, ierr, arg_pos=4_int32)
        call validate_in_range_real(shared_residual_range_high, ierr, arg_pos=5_int32, min=shared_residual_range_low)
        call validate_in_range_int(m_min, ierr, arg_pos=11_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=12_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=13_int32, min=0_int32)
        call validate_all_in_range_real(pooled_residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(pooled_residuals_perm(n_residuals))
        call init_perm(pooled_residuals_perm)
        call sort_array_heapsort(pooled_residuals, pooled_residuals_perm)

        call determine_bin_count_occupancy_exhaustive_impl(&
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
            min_residuals_per_bin = min_residuals_per_bin&
        )
    end subroutine determine_bin_count_occupancy_exhaustive

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):determine_bin_count_occupancy_exhaustive]] does both.
    !| Tests every candidate bin count `M` in `[m_min, m_max]` and keeps the largest one whose
    !| pooled histogram satisfies the occupancy criterion. A reference for validating
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]],
    !| whose faster search could in principle miss a larger admissible `M`. Takes the range and
    !| `n_pooled_residuals` from a prior call of that routine, so only the `M` selection is compared.
    pure subroutine determine_bin_count_occupancy_exhaustive_expert(&
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
        integer(int32), intent(in) :: n_pooled_residuals
            !! Count of non-NaN pooled residuals (N_j), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_residuals`.
        real(real64), intent(in) :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
        real(real64), intent(in) :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high), from the same prior call as
            !! `shared_residual_range_low`
            !! The minimum valid value is `shared_residual_range_low`.
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
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_pooled_residuals, ierr, arg_pos=4_int32, min=0_int32, max=n_residuals)
        call validate_in_range_real(shared_residual_range_low, ierr, arg_pos=5_int32)
        call validate_in_range_real(shared_residual_range_high, ierr, arg_pos=6_int32, min=shared_residual_range_low)
        call validate_in_range_int(m_min, ierr, arg_pos=12_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=13_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=14_int32, min=0_int32)
        call validate_all_in_range_real(pooled_residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_residuals_perm, n_residuals, ierr, arg_pos=2_int32, min=1_int32, max=n_residuals)
        if (is_err(ierr)) return
#endif

        call determine_bin_count_occupancy_exhaustive_impl(&
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
            min_residuals_per_bin = min_residuals_per_bin&
        )
    end subroutine determine_bin_count_occupancy_exhaustive_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):determine_bin_count_stone_expert]] to prepare it yourself.
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
    pure subroutine determine_bin_count_stone(&
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
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        real(real64), dimension(n_residuals), intent(in) :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(int32), dimension(n_residuals), intent(in) :: residual_to_study
            !! Study each pooled residual comes from
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_studies`.
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
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
        integer(int32), intent(in), optional :: bin_count_mode
            !! How the studies' scores are combined into one bin count
            !!
            !! | Mode                                                                                  | Value                                                                                           |
            !! |---------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Minimise the sum of the per-study scores (default)                                    | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled                                   | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `1_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin of the pooled clamped histogram must
            !! reach for a candidate bin count to be admissible (0 turns the floor off)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=18_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=19_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=20_int32, min=0_int32)
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=21_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=22_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(pooled_residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(residual_to_study, n_residuals, ierr, arg_pos=3_int32, min=1_int32, max=n_studies)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=17_int32); end if
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(pooled_residuals_perm(n_residuals))
        M_ALLOCATE(tmp_bin_counts(256))
        M_ALLOCATE(tmp_residual_to_study_perm(n_residuals))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        call init_perm(pooled_residuals_perm)
        call sort_array_heapsort(pooled_residuals, pooled_residuals_perm)

        call determine_bin_count_stone_impl(&
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
            bin_count_mode = bin_count_mode,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile&
        )
    end subroutine determine_bin_count_stone

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):determine_bin_count_stone]] does both.
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
    pure subroutine determine_bin_count_stone_expert(&
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
        )
        integer(int32), intent(in) :: n_residuals
            !! Number of pooled residuals
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_residuals), intent(in) :: pooled_residuals
            !! Pooled signed residuals for one neighborhood, across all its neighbors and all studies
            !! NaN is permitted for this value.
        integer(int32), dimension(n_residuals), intent(in) :: pooled_residuals_perm
            !! Sorting permutation for `pooled_residuals`, ascending, NaN last
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_residuals`.
        integer(int32), dimension(n_residuals), intent(in) :: residual_to_study
            !! Study each pooled residual comes from
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_studies`.
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Neighborhood size of the candidate under test
            !! The minimum valid value is `1_int32`.
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
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts
            !! Working array: one study's bin counts at the bin count under test
        integer(int32), dimension(n_residuals), intent(out) :: tmp_residual_to_study_perm
            !! Working array: indices of `pooled_residuals` grouped by study, each study sorted
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array: number of pooled residuals of each study
        integer(int32), intent(in), optional :: bin_count_mode
            !! How the studies' scores are combined into one bin count
            !!
            !! | Mode                                                                                  | Value                                                                                           |
            !! |---------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Minimise the sum of the per-study scores (default)                                    | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled                                   | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
            !! The default value is `1_int32`.
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `3_int32`.
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_N_BINS`.
            !! The default value is `120_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin of the pooled clamped histogram must
            !! reach for a candidate bin count to be admissible (0 turns the floor off)
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_residuals, ierr, arg_pos=3_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=5_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=22_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=23_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=24_int32, min=0_int32)
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=25_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=26_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(pooled_residuals, n_residuals, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_residuals_perm, n_residuals, ierr, arg_pos=2_int32, min=1_int32, max=n_residuals)
        call validate_all_in_range_int(residual_to_study, n_residuals, ierr, arg_pos=4_int32, min=1_int32, max=n_studies)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=21_int32); end if
        if (is_err(ierr)) return
#endif

        call determine_bin_count_stone_impl(&
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
            bin_count_mode = bin_count_mode,&
            m_min = m_min,&
            m_max = m_max,&
            min_residuals_per_bin = min_residuals_per_bin,&
            lower_residual_range_quantile = lower_residual_range_quantile,&
            upper_residual_range_quantile = upper_residual_range_quantile&
        )
    end subroutine determine_bin_count_stone_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl]].
    !| Starting from `n_points_high = clamp(ceiling(4*sqrt(G)), MIN_POINTS, MAX_POINTS)` (`G` =
    !| `max_n_genes_all_studies`), multiplies it by GAMMA until it drops below
    !| `n_points_low = max(MIN_POINTS, ceiling(0.2*n_points_high))`, and pairs each distinct
    !| `n_points` with up to `size(KX_FACTORS)` distinct `n_neighbors = max(1, floor(G/(f*n_points_high)))`.
    !| Duplicate values from clamping or rounding collapse, so a small enough `G` yields a single
    !| candidate. Bin counts are not part of a candidate; they are chosen per reference point later.
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
            !! Number of candidate pairs actually filled (at most 16)
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

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl]].
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

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):construct_adaptive_neighborhoods_expert]] to prepare it yourself.
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
    pure subroutine construct_adaptive_neighborhoods(&
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
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        integer(int32), intent(out) :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! The first `n_points` elements will hold the results.
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! The first `n_points` elements will hold the results.
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! The first `n_points` elements will hold the results.
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: stop_reason
            !! Why each neighborhood stopped growing: `1_int32` tau exceeded,
            !! `2_int32` `k_max` reached, `3_int32` residual cap
            !! reached, `4_int32` pool exhausted, `5_int32`
            !! zero dispersion, `6_int32` too few residuals
            !! The first `n_points` elements will hold the results.
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! The first `n_points` elements will hold the results.
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! The first `n_points` elements will hold the results.
        integer(int32), intent(out) :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(int32), intent(out) :: construction_status
            !! Overall outcome: `0_int32` ok, `1_int32` too few
            !! means, `2_int32` a study neighborhood below
            !! `min_study_neighbors`, `3_int32` a neighborhood with too
            !! few residuals
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(k_step, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(k_max, ierr, arg_pos=8_int32, min=k_start)
        call validate_in_range_real(tau, ierr, arg_pos=18_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=19_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=20_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=21_int32, min=1_int32)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))

        call construct_adaptive_neighborhoods_impl(&
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
            min_study_neighbors = min_study_neighbors&
        )
    end subroutine construct_adaptive_neighborhoods

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):construct_adaptive_neighborhoods]] does both.
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
    pure subroutine construct_adaptive_neighborhoods_expert(&
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
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        integer(int32), intent(out) :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! The first `n_points` elements will hold the results.
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! The first `n_points` elements will hold the results.
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! The first `n_points` elements will hold the results.
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: stop_reason
            !! Why each neighborhood stopped growing: `1_int32` tau exceeded,
            !! `2_int32` `k_max` reached, `3_int32` residual cap
            !! reached, `4_int32` pool exhausted, `5_int32`
            !! zero dispersion, `6_int32` too few residuals
            !! The first `n_points` elements will hold the results.
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! The first `n_points` elements will hold the results.
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! The first `n_points` elements will hold the results.
        integer(int32), intent(out) :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(int32), intent(out) :: construction_status
            !! Overall outcome: `0_int32` ok, `1_int32` too few
            !! means, `2_int32` a study neighborhood below
            !! `min_study_neighbors`, `3_int32` a neighborhood with too
            !! few residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, seeded and sorted here
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(k_step, ierr, arg_pos=7_int32, min=1_int32)
        call validate_in_range_int(k_max, ierr, arg_pos=8_int32, min=k_start)
        call validate_in_range_real(tau, ierr, arg_pos=19_int32, min=0.0_real64)
        call validate_in_range_real(mad_distance_factor, ierr, arg_pos=20_int32, min=0.0_real64)
        call validate_in_range_int(max_pooled_residuals, ierr, arg_pos=21_int32, min=0_int32)
        call validate_in_range_int(min_study_neighbors, ierr, arg_pos=22_int32, min=1_int32)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        call construct_adaptive_neighborhoods_impl(&
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
            min_study_neighbors = min_study_neighbors&
        )
    end subroutine construct_adaptive_neighborhoods_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl]].
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
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
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
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
    !| The second admissibility gate of a parameter-search candidate, checked after
    !| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]].
    !| Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    !| them does not count as a failure.
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
            !! Second extent of `mean_pmf_counts`: the widest bin count of any point
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the mean pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: n_bins_per_point
            !! Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
            !! is inspected
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
        call validate_all_in_range_int(mean_pmf_counts, n_points * n_bins, ierr, arg_pos=1_int32, min=0_int32)
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
    !| For each study, the overlap of the new candidate's interval with the best one so far is
    !| measured as a fraction of the new interval's width, and counted when it reaches
    !| `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    !| best (the best is then kept), or, after the new candidate became the best, when the count
    !| meets the condition `join_method` names.
    pure subroutine check_plateau_condition(&
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
        )
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(2, n_studies), intent(in) :: confidence_interval
            !! Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
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
            !! How the studies' overlaps are joined into the plateau condition; the median method
            !! requires `count > (n_studies - 1) / 2`
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        real(real64), intent(in) :: succeeding_ci_overlap
            !! Minimum fractional overlap with the best candidate's interval for a study to pass
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_studies), intent(out) :: ci_overlap
            !! Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            !! it was on entry, relative to the new interval's width; always filled. Meaningless
            !! for the very first candidate, whose best interval is only a seed
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
            ci_overlap = ci_overlap,&
            plateau_found = plateau_found&
        )
    end subroutine check_plateau_condition

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):check_effect_size_plateau_condition_expert]] to prepare it yourself.
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
            !! `.false.` for the first admissible candidate, which has no transition
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
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
            !! `.false.` for the first admissible candidate, which has no transition
        real(real64), intent(in) :: delta_median_threshold
            !! Upper bound on the median relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_max_threshold
            !! Upper bound on the largest relative change for a transition to qualify
            !! The minimum valid value is `above(0.0_real64)`.
        real(real64), intent(in) :: delta_epsilon
            !! Small constant preventing division by zero when a study's previous JSD was zero
            !! The minimum valid value is `above(0.0_real64)`.
        integer(int32), intent(in) :: delta_min_consecutive_transitions
            !! Number of consecutive qualifying transitions required to declare a plateau
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(inout) :: n_consecutive_ok
            !! Running count of consecutive qualifying transitions; reset to zero when this one does
            !! not qualify or `has_previous` is `.false.`
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
    !| `mean_pmf` is the equal-weight mean of all studies' pmfs, `mean_pmf_counts` and
    !| `mean_pmf_included_n_reps` the sums of their counts. Known limitation: every study is
    !| compared against a consensus that includes itself, not a leave-one-out background.
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
            !! Second extent of the bin-sized arrays: the widest bin count of any reference point.
            !! All studies must share each point's bin count, so their zero padding coincides
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), dimension(n_points, n_bins, n_studies), intent(in) :: pmfs
            !! Per-study probabilities of each bin per reference point, from
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_points, n_bins, n_studies), intent(in) :: counts
            !! Absolute counts of a residual per bin for `pmfs`
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study; every study needs at
            !! least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_points, n_bins), intent(out) :: mean_pmf
            !! The consensus pmf, `sum(pmfs, dim=3) / n_studies`
        integer(int32), dimension(n_points), intent(out) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point, summed over the studies
        integer(int32), dimension(n_points, n_bins), intent(out) :: mean_pmf_counts
            !! Residual counts per bin of the consensus, `sum(counts, dim=3)`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_bins, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(pmfs, n_points * n_bins * n_studies, ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(counts, n_points * n_bins * n_studies, ierr, arg_pos=2_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=6_int32, min=1_int32)
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
    !| As [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]], for
    !| callers that need no counts (the bootstrap). Same known limitation: the consensus includes
    !| the study compared against it.
    pure subroutine create_mean_pmf_only(&
            pmfs,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf,&
            ierr&
        )
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
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_points, n_bins), intent(out) :: mean_pmf
            !! The consensus pmf, `sum(pmfs, dim=3) / n_studies`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_bins, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(pmfs, n_points * n_bins * n_studies, ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
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
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)
        integer(int32) :: n_bootstrapping_top_k_jsds
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=9_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(mean_pmf_counts, n_points * n_bins, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=6_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=7_int32, min=1_int32)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_n_top_k_jsds(&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            n_top_k = n_bootstrapping_top_k_jsds&
        )
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_counts(n_points, n_bins))
        M_ALLOCATE(tmp_pmfs(n_points, n_bins, n_studies))
        M_ALLOCATE(tmp_mean_pmf(n_points, n_bins))
        M_ALLOCATE(tmp_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_weights(n_points, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_point_pool(n_bins))
        M_ALLOCATE(tmp_point_draw(n_bins))

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
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine bootstrap_histogram

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):bootstrap_histogram]] does both.
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
            tmp_point_pool,&
            tmp_point_draw,&
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
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
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
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=19_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(mean_pmf_counts, n_points * n_bins, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=6_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=7_int32, min=1_int32)
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
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine bootstrap_histogram_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):run_js_comp_test_expert]] to prepare it yourself.
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
            !! The maximum valid value is `max_n_genes_all_studies`.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
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
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_pool
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_draw
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
        integer(int32), dimension(:), allocatable :: tmp_point_n_neighbors
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
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_permutations, ierr, arg_pos=34_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=36_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=37_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=38_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=39_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=40_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=41_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(gene_means_perms, max_n_genes_all_studies * n_studies, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=8_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=9_int32)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=42_int32); end if
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, n_neighbors, n_points))
        M_ALLOCATE(tmp_permutation_mean_pmf(n_points, 256))
        M_ALLOCATE(tmp_permutation_point_pool(256))
        M_ALLOCATE(tmp_permutation_point_draw(256))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*n_neighbors*n_studies))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_point_n_neighbors(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_pmfs(n_points, 256, n_studies))
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
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):run_js_comp_test]] does both.
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
            !! The maximum valid value is `max_n_genes_all_studies`.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
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
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points), intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        real(real64), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array of the permutation test: consensus pmf of the permuted studies
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_pool
            !! Working array of the permutation test: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_draw
            !! Working array of the permutation test: one study's draw at one point
        real(real64), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
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
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array of the permutation test: the remaining pool
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_counts
            !! Working array of the permutation test: one study's resampled counts
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array of the permutation test: resampled pmfs
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array of the permutation test: per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array of the permutation test: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array of the permutation test: global JSD values
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=5_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_permutations, ierr, arg_pos=51_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=53_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=54_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=55_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=56_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=57_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=58_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(gene_means_perms, max_n_genes_all_studies * n_studies, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=8_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=9_int32)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=59_int32); end if
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
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive_expert]] to prepare it yourself.
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
    subroutine run_js_comp_test_adaptive(&
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
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
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
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        integer(int32), dimension(:, :), allocatable :: tmp_point_neighborhood_indices
        real(real64), dimension(:, :), allocatable :: tmp_neighbor_residuals
        real(real64), dimension(:, :), allocatable :: tmp_permutation_mean_pmf
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_pool
        integer(int32), dimension(:), allocatable :: tmp_permutation_point_draw
        real(real64), dimension(:), allocatable :: tmp_pooled_residuals
        integer(int32), dimension(:), allocatable :: tmp_pooled_residuals_perm
        integer(int32), dimension(:), allocatable :: tmp_bin_counts_search
        integer(int32), dimension(:), allocatable :: tmp_residual_to_study_perm
        integer(int32), dimension(:), allocatable :: tmp_study_block_ends
        integer(int32), dimension(:), allocatable :: tmp_study_block_sizes
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
        call validate_in_range_int(n_permutations, ierr, arg_pos=31_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=33_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=34_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=35_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=36_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=37_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=38_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_neighborhood_range, 2 * n_points, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies*n_studies)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=39_int32); end if
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        M_ALLOCATE(tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf(n_points, 256))
        M_ALLOCATE(tmp_permutation_point_pool(256))
        M_ALLOCATE(tmp_permutation_point_draw(256))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        M_ALLOCATE(tmp_study_block_ends(n_studies))
        M_ALLOCATE(tmp_study_block_sizes(n_studies))
        M_ALLOCATE(tmp_permutation_mean_pmf_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_counts(n_points, 256))
        M_ALLOCATE(tmp_permutation_pmfs(n_points, 256, n_studies))
        M_ALLOCATE(tmp_permutation_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_weights(n_points, n_studies))
        M_ALLOCATE(tmp_permutation_global_js_divergence(n_studies))

        call run_js_comp_test_adaptive_impl(&
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
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive]] does both.
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
    subroutine run_js_comp_test_adaptive_expert(&
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
            !! Number of reference points (neighborhoods), the construction's own `n_points`
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene; the array the
            !! neighborhoods were constructed from
            !! NaN is permitted for this value.
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(int32), dimension(2, n_points), intent(in) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last sorted pooled position, from
            !! construct_adaptive_neighborhoods
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies*n_studies`.
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
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: pmfs
            !! `counts` normalized per reference point and study; only bins `1:max_n_bins_per_point`
            !! are meaningful
        integer(int32), dimension(n_points, 256, n_studies), intent(out) :: counts
            !! Residual counts per bin for `pmfs`; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, n_studies), intent(out) :: included_n_reps
            !! Count of non-NaN residuals binned per reference point, per study
        real(real64), dimension(n_points, 256), intent(out) :: mean_pmf
            !! The consensus pmf; only bins `1:max_n_bins_per_point` are meaningful
        integer(int32), dimension(n_points, 256), intent(out) :: mean_pmf_counts
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
        real(real64), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf
            !! Working array forwarded to the permutation test
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_pool
            !! Working array forwarded to the permutation test
        integer(int32), dimension(256), intent(out) :: tmp_permutation_point_draw
            !! Working array forwarded to the permutation test
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
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_mean_pmf_counts
            !! Working array forwarded to the permutation test
        integer(int32), dimension(n_points, 256), intent(out) :: tmp_permutation_counts
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, 256, n_studies), intent(out) :: tmp_permutation_pmfs
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_js_divergences
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_permutation_weights
            !! Working array forwarded to the permutation test
        real(real64), dimension(n_studies), intent(out) :: tmp_permutation_global_js_divergence
            !! Working array forwarded to the permutation test
        integer(int32), intent(in), optional :: n_permutations
            !! Number of permutations of the permutation test; `0` skips it (p-values NaN)
            !! The minimum valid value is `0_int32`.
            !! The default value is `1000_int32`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a bin count to be
            !! admissible, in every `bin_count_mode`
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.05_real64`.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! The minimum valid value is `lower_residual_range_quantile`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT or ERR_ALLOC_FAIL (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, arg_pos=49_int32, min=0_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=51_int32, min=0_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=52_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=53_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=54_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=55_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=56_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=6_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(pooled_neighborhood_range, 2 * n_points, ierr, arg_pos=7_int32, min=1_int32, max=max_n_genes_all_studies*n_studies)
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=57_int32); end if
        if (is_err(ierr)) return
#endif

        call run_js_comp_test_adaptive_impl(&
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
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_adaptive_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search_expert]] to prepare it yourself.
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
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
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
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
            !!
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
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
        integer(int32), dimension(16), intent(out) :: trace_n_points
            !! Per-admissible-candidate `n_points`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
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
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width (the values
            !! `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            !! like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum count of every consensus bin (gate 2) and the floor of the bin-count search,
            !! one threshold for both, so a point the search marks as failed also fails the gate
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods (gate 1)
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created
        integer(int32) :: n_bootstrapping_top_k_jsds
        integer(int32), dimension(:, :), allocatable :: tmp_gene_means_perms
        integer(int32), dimension(:), allocatable :: tmp_gene_means_perm_all
        real(real64), dimension(:), allocatable :: tmp_x_star
        integer(int32), dimension(:, :, :), allocatable :: tmp_neighborhood_indices_all_studies
        integer(int32), dimension(:, :), allocatable :: tmp_neighborhood_range
        real(real64), dimension(:, :, :), allocatable :: tmp_neighborhood_residuals_gathered
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw
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
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=8_int32, min=max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)))
        call validate_in_range_int(max_n_neighbors_candidate, ierr, arg_pos=9_int32, min=max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32)))
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=39_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=40_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=41_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=43_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=44_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=45_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=46_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=47_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=48_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=49_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=50_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=51_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=52_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=42_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=54_int32); end if
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
        M_ALLOCATE(tmp_point_pool(256))
        M_ALLOCATE(tmp_point_draw(256))
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
        M_ALLOCATE(tmp_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_prev_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_delta_perm(n_studies))
        M_ALLOCATE(tmp_best_uncertainty_confidence_interval(2, n_studies))
        M_ALLOCATE(tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        M_ALLOCATE(tmp_bin_counts_search(256))
        M_ALLOCATE(tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
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
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_parameter_search

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):run_js_comp_test_parameter_search]] does both.
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
            !! Upper bound on the grid's `n_points` candidates, which sizes the per-point outputs;
            !! at least what `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS))`.
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Upper bound on the grid's `n_neighbors` candidates; at least what
            !! `calc_js_comp_test_candidate_bounds` returns, else invalid input
            !! It is recommended to compute this argument from the `max_n_neighbors_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds]].
            !! The minimum valid value is `max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32))`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
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
            !! How the studies' confidence-interval overlaps are joined into the plateau condition
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
        integer(int32), dimension(16), intent(out) :: trace_n_points
            !! Per-admissible-candidate `n_points`
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(16), intent(out) :: trace_n_neighbors
            !! Per-admissible-candidate `n_neighbors`
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_global_js_divergence
            !! Per-admissible-candidate, per-study observed global JSD (`J_{i,t}`)
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
            !! (`W_{i,t} / max(J_{i,t}, delta_epsilon)`)
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_ci_overlap
            !! Per-admissible-candidate, per-study overlap in `[0, 1]` of the candidate's confidence
            !! interval with the best candidate's so far, relative to its own width (the values
            !! `succeeding_ci_overlap` is tested against); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(n_studies, 16), intent(out) :: trace_delta
            !! Per-admissible-candidate, per-study relative JSD change from the previous admissible
            !! candidate (`Delta_{i,t}`); `-1.0` for the first admissible candidate
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_median
            !! Per-admissible-candidate median of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(16), intent(out) :: trace_delta_max
            !! Per-admissible-candidate maximum of trace_delta across studies; `-1.0` for the first
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_selected_n_bins
            !! Per-admissible-candidate, per-reference-point bin count `M_j`. Jagged: only rows
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_n_pooled_residuals
            !! Per-admissible-candidate, per-reference-point pooled residual count (`N_j`). Jagged
            !! like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_min_bin_occupancy
            !! Per-admissible-candidate, per-reference-point minimum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_mean_bin_occupancy
            !! Per-admissible-candidate, per-reference-point mean bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_max_bin_occupancy
            !! Per-admissible-candidate, per-reference-point maximum bin occupancy. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_sturges_bins
            !! Per-admissible-candidate, per-reference-point Sturges' rule diagnostic. Jagged like
            !! trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_points_candidate, 16), intent(out) :: trace_fd_bins
            !! Per-admissible-candidate, per-reference-point Freedman-Diaconis rule diagnostic.
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_low
            !! Per-admissible-candidate, per-reference-point lower residual-range bound (R_low).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        real(real64), dimension(max_n_points_candidate, 16), intent(out) :: trace_shared_residual_range_high
            !! Per-admissible-candidate, per-reference-point upper residual-range bound (R_high).
            !! Jagged like trace_selected_n_bins
            !! The first `n_admissible_evaluated` elements will hold the results.
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array: each study's own sorting permutation for `gene_means`
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array: the current candidate's reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_indices_all_studies
            !! Working array: every study's neighbor gene indices for the current candidate
        integer(int32), dimension(2, max_n_points_candidate), intent(out) :: tmp_neighborhood_range
            !! Working array: one study's neighborhood spans, reused per study
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array: one study's gathered neighborhood residuals, reused per study
        integer(int32), dimension(256), intent(out) :: tmp_point_pool
            !! Working array of the bootstrap: one point's pool
        integer(int32), dimension(256), intent(out) :: tmp_point_draw
            !! Working array of the bootstrap: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array: the current candidate's bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array: the current candidate's R_low per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array: the current candidate's R_high per reference point
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
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals
            !! Working array: one reference point's pooled residuals, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_pooled_residuals_perm
            !! Working array: sorting permutation for tmp_pooled_residuals, reused per point
        integer(int32), dimension(256), intent(out) :: tmp_bin_counts_search
            !! Working array: bin counts of the bin-count search, reused per point
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: tmp_residual_to_study_perm
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
            !! The minimum valid value is `0_int32`.
            !! The default value is `10_int32`.
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Minimum fractional overlap of two consecutive neighborhoods (gate 1)
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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
        integer(int32), intent(in), optional :: bin_count_mode
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
            !! The default value is `0_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if the bootstrap's random number generator could not be created

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(max_n_reps_all_studies, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=6_int32, min=1_int32)
        call validate_in_range_int(max_n_points_candidate, ierr, arg_pos=8_int32, min=max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)))
        call validate_in_range_int(max_n_neighbors_candidate, ierr, arg_pos=9_int32, min=max(1_int32, floor(real(max_n_genes_all_studies, real64)/(KX_FACTORS(1)*real(max(MIN_POINTS, ceiling(0.2_real64*real(max(MIN_POINTS, min(ceiling(4.0_real64*sqrt(real(max_n_genes_all_studies, real64)), kind=int32), MAX_POINTS)), real64), kind=int32)), real64)), kind=int32)))
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=10_int32, min=1_int32)
        call validate_in_range_int(min_residuals_per_bin, ierr, arg_pos=84_int32, min=0_int32)
        call validate_in_range_real(min_neighbor_overlap, ierr, arg_pos=85_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(succeeding_ci_overlap, ierr, arg_pos=86_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(delta_median_threshold, ierr, arg_pos=88_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_max_threshold, ierr, arg_pos=89_int32, min=above(0.0_real64))
        call validate_in_range_real(delta_epsilon, ierr, arg_pos=90_int32, min=above(0.0_real64))
        call validate_in_range_int(delta_min_consecutive_transitions, ierr, arg_pos=91_int32, min=1_int32)
        call validate_in_range_int(m_min, ierr, arg_pos=92_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_int(m_max, ierr, arg_pos=93_int32, min=1_int32, max=MAX_N_BINS)
        call validate_in_range_real(gamma_occupancy, ierr, arg_pos=94_int32, min=above(1.0_real64))
        call validate_in_range_real(lower_residual_range_quantile, ierr, arg_pos=95_int32, min=0.0_real64, max=1.0_real64)
        call validate_in_range_real(upper_residual_range_quantile, ierr, arg_pos=96_int32, min=lower_residual_range_quantile, max=1.0_real64)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=97_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_real(gene_means, max_n_genes_all_studies * n_studies, ierr, arg_pos=4_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies, ierr, arg_pos=5_int32, allow_nan=.true._c_bool)
        if (join_method /= METHOD_JOIN_MIN .and. join_method /= METHOD_JOIN_MAX .and. join_method /= METHOD_JOIN_MEDIAN) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (present(plateau_mode)) then; if (plateau_mode /= MODE_PLATEAU_CI_OVERLAP .and. plateau_mode /= MODE_PLATEAU_EFFECT_SIZE .and. plateau_mode /= MODE_PLATEAU_BOTH) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=87_int32); end if
        if (present(bin_count_mode)) then; if (bin_count_mode /= MODE_BIN_COUNT_OCCUPANCY .and. bin_count_mode /= MODE_BIN_COUNT_STONE_EQUAL .and. bin_count_mode /= MODE_BIN_COUNT_STONE_SIZE_WEIGHTED .and. bin_count_mode /= MODE_BIN_COUNT_STONE_POOLED) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=99_int32); end if
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
            bin_count_mode = bin_count_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine run_js_comp_test_parameter_search_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_parameter_search_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive_parameter_search_expert]] to prepare it yourself.
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
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds]].
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
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
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
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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
            n_top_k = n_bootstrapping_top_k_jsds&
        )
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

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_parameter_search_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_js_comp_test(module):run_js_comp_test_adaptive_parameter_search]] does both.
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
            !! It is recommended to compute this argument from the `max_n_points_candidate` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds]].
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
            !! | Method                                  | Value                                                                           |
            !! |-----------------------------------------|---------------------------------------------------------------------------------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MIN(variable)]]    |
            !! | Maximum overlap (any one study passes)  | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MAX(variable)]]    |
            !! | Median overlap (a majority must pass)   | [[tox_data_integration_js_comp_test_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
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
            !! `1:trace_n_points(t)` of column `t` are meaningful, so a Python/R caller slices
            !! `[:trace_n_points[t], t]`
            !! The first `n_admissible_evaluated` elements will hold the results.
        logical(c_bool), dimension(max_n_points_candidate, 16), intent(out) :: trace_occupancy_failed
            !! Per-admissible-candidate, per-reference-point flag: `.true.` iff the bin-count search
            !! failed and the point used `m_min` bins. Jagged like trace_selected_n_bins
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
            !! | Mode                                | Value                                                                                 |
            !! |-------------------------------------|---------------------------------------------------------------------------------------|
            !! | CI overlap only                     | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]]  |
            !! | Relative-effect-size stability only | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion                    | [[tox_data_integration_js_comp_test_impl(module):MODE_PLATEAU_BOTH(variable)]]        |
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
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
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
            !! ([[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode                                                                                                                       | Value                                                                                           |
            !! |----------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` (default) | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]           |
            !! | Stone's cross-validation: minimise the sum of every study's own score                                                      | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]]         |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range        | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled                                              | [[tox_data_integration_js_comp_test_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]]        |
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

end module tox_data_integration_js_comp_test
