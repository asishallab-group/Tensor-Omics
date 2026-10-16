#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bin Counts
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
module tox_data_integration_js_comp_test_binning_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_math_impl, only: is_close, LOG_2
    use f42_stats_impl, only: calc_quantile_impl, histogram_bin_counts, sorted_bin_counts
    use f42_sort_impl, only: init_perm, sort_array_heapsort, block_containing_index
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code
    use tox_data_integration_jsd_impl, only: build_residual_histograms_impl
    use tox_data_integration_preprocessing_impl, only: gather_pooled_neighborhood_residuals, &
                                                       materialize_pooled_neighborhood
    M_IMPLICIT_NONE
    private
    public :: estimate_bin_count_impl, determine_bin_count_occupancy_impl, &
             determine_bin_count_occupancy_exhaustive_impl, determine_bin_count_stone_impl, &
             group_sorted_perm_by_study, determine_point_bin_count, &
             build_point_study_histogram, build_fixed_k_study_histograms, build_adaptive_histograms, &
             MAX_N_BINS, MODE_BIN_COUNT_OCCUPANCY, MODE_BIN_COUNT_STONE_EQUAL, MODE_BIN_COUNT_STONE_SIZE_WEIGHTED, &
             MODE_BIN_COUNT_STONE_POOLED

    ! M_MAX_N_BINS and the M_MODE_BIN_COUNT_* values are defined in src/tox/data_integration/macros.h:
    ! the search and final-run modules need the same literals as dummy extents and DM_DEFAULTs.
    integer(int32), parameter :: MAX_N_BINS = M_MAX_N_BINS
        !! Fixed ceiling on any histogram bin count; histogram work arrays are sized by it

    ! `bin_count_mode`'s values: how a reference point's histogram bin count is chosen. Occupancy
    ! is the default of all four pipeline entry points; determine_bin_count_stone_impl rejects it.
    integer(int32), parameter :: MODE_BIN_COUNT_OCCUPANCY = M_MODE_BIN_COUNT_OCCUPANCY
        !! Bin-count mode: the largest bin count whose pooled histogram has every bin at or above
        !! the occupancy floor
    integer(int32), parameter :: MODE_BIN_COUNT_STONE_EQUAL = M_MODE_BIN_COUNT_STONE_EQUAL
        !! Bin-count mode: Stone's cross-validation score of every study, summed with equal weights
    integer(int32), parameter :: MODE_BIN_COUNT_STONE_SIZE_WEIGHTED = M_MODE_BIN_COUNT_STONE_SIZE_WEIGHTED
        !! Bin-count mode: Stone's cross-validation score of every study, weighted by the study's
        !! residual count inside the range
    integer(int32), parameter :: MODE_BIN_COUNT_STONE_POOLED = M_MODE_BIN_COUNT_STONE_POOLED
        !! Bin-count mode: Stone's cross-validation score of the pooled residuals of all studies

contains

    !> summary: Estimate the histogram bin count for one (n_points, n_neighbors) candidate
    !| AUTHOR_LASZLO_LANG
    !| Sturges' rule and the Freedman-Diaconis rule for a histogram over `[-R, R]`
    !| (`R = shared_residual_range`), each clamped to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):MAX_N_BINS(variable)]] and returned as
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):estimate_bin_count_impl(interface)]] and are
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
        integer(int32), intent(out) :: tmp_bin_counts(M_MAX_N_BINS)
            !! Working array: per-bin counts of the candidate bin count under test
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count the search will ever test (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count the search will ever test (M_max); raised to `m_min` if
            !! smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor for the coarse search stage; must exceed 1 or the search
            !! never advances
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(M_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        real(real64) :: actual_gamma_occupancy, actual_lower_residual_range_quantile, actual_upper_residual_range_quantile
        integer(int32) :: i_pool, trial_m, next_m, m_valid, m_invalid, best_m, min_occ, best_min_occ, best_max_occ
        integer(int32) :: combined_n_bins_unused
        real(real64) :: half_span
        logical(c_bool) :: found_valid

        M_DEFAULT_VAL(m_min, actual_m_min, M_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, M_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(gamma_occupancy, actual_gamma_occupancy, M_OCCUPANCY_GAMMA_DEFAULT)
        M_DEFAULT_VAL(lower_residual_range_quantile, actual_lower_residual_range_quantile, M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        M_DEFAULT_VAL(upper_residual_range_quantile, actual_upper_residual_range_quantile, M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]],
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
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
            !! call
            !! DM_MIN(0_int32)
            !! DM_MAX(n_residuals)
        real(real64), intent(in) :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low), from a prior
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
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
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach for a candidate bin count to
            !! be admissible (n_min)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        integer(int32) :: trial_m, best_m, min_occ
        integer(int32) :: candidate_bin_counts(MAX_N_BINS)

        M_DEFAULT_VAL(m_min, actual_m_min, M_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, M_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]],
    !| also used by callers whose pool already holds every study as one contiguous block. Inputs:
    !|
    !| - `pooled_residuals_perm` orders all of `pooled_residuals` ascending, NaN last.
    !| - `residual_to_study_perm` holds indices into `pooled_residuals`, grouped by study (study `k`
    !|   fills positions `sum(study_block_sizes(1:k-1)) + 1` to `sum(study_block_sizes(1:k))`), each
    !|   block ascending with NaN last; `group_sorted_perm_by_study` builds it.
    !|
    !| The range is as in
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]].
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
            !! Smallest bin count tested; default M_OCCUPANCY_M_MIN_DEFAULT
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count tested, raised to m_min if smaller; default M_OCCUPANCY_M_MAX_DEFAULT
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Floor on every bin of the pooled clamped histogram; default
            !! M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile of R_low; default M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile of R_high; default M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT

        integer(int32) :: actual_m_min, actual_m_max, actual_min_residuals_per_bin
        real(real64) :: actual_lower_residual_range_quantile, actual_upper_residual_range_quantile
        integer(int32) :: i_pool, i_study, block_start, block_end, trial_m, best_m, best_min_occ, best_max_occ
        integer(int32) :: n_eligible_studies, n_total_inside, n_study_below, n_study_inside
        integer(int32) :: pooled_n_below, pooled_n_inside, pooled_n_above, min_occ, max_occ, combined_n_bins_unused
        integer(int32) :: pooled_counts(MAX_N_BINS)
        real(real64) :: half_span, range_width, score, best_score
        logical(c_bool) :: pooled_scoring, size_weighted, found_admissible

        M_DEFAULT_VAL(m_min, actual_m_min, M_OCCUPANCY_M_MIN_DEFAULT)
        M_DEFAULT_VAL(m_max, actual_m_max, M_OCCUPANCY_M_MAX_DEFAULT)
        M_DEFAULT_VAL(min_residuals_per_bin, actual_min_residuals_per_bin, M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        M_DEFAULT_VAL(lower_residual_range_quantile, actual_lower_residual_range_quantile, M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        M_DEFAULT_VAL(upper_residual_range_quantile, actual_upper_residual_range_quantile, M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
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
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]].
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
        integer(int32), intent(out) :: tmp_bin_counts(M_MAX_N_BINS)
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
            !! | Minimise the sum of the per-study scores (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Minimise the sum of the per-study scores, each multiplied by the study's inside count | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_STONE_EQUAL)
        integer(int32), intent(in), optional :: m_min
            !! Smallest candidate bin count tested (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest candidate bin count tested (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin of the pooled clamped histogram must
            !! reach for a candidate bin count to be admissible (0 turns the floor off)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own lower residual-range bound
            !! (shared_residual_range_low)
            !! Must not exceed upper_residual_range_quantile (default 0.95); otherwise the call fails
            !! with invalid input reported on upper_residual_range_quantile.
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] for this neighborhood's own upper residual-range bound
            !! (shared_residual_range_high)
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)

        integer(int32) :: actual_bin_count_mode

        M_DEFAULT_VAL(bin_count_mode, actual_bin_count_mode, M_MODE_BIN_COUNT_STONE_EQUAL)

        call build_residual_to_study_perm(pooled_residuals_perm, n_residuals, residual_to_study, n_studies, &
                                          tmp_residual_to_study_perm, tmp_study_block_sizes)

        call stone_bin_count_core(pooled_residuals, pooled_residuals_perm, n_residuals, n_studies, &
                                  tmp_residual_to_study_perm, tmp_study_block_sizes, max_n_reps_all_studies, n_neighbors, &
                                  actual_bin_count_mode, selected_n_bins, stone_failed, shared_residual_range_low, &
                                  shared_residual_range_high, n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                  max_bin_occupancy, sturges_bins, fd_bins, tmp_bin_counts, m_min, m_max, &
                                  min_residuals_per_bin, lower_residual_range_quantile, upper_residual_range_quantile)
    end subroutine determine_bin_count_stone_impl

    !> summary: Choose one reference point's histogram bin count with the rule `bin_count_mode` names
    !| AUTHOR_LASZLO_LANG
    !| Chooses one reference point's bin count from its gathered, sorted pool with the rule
    !| `bin_count_mode` names.
    !|
    !| The pool `pooled_residuals(1:n_residuals)` holds every study's residuals as one contiguous
    !| block, study after study: study `s` fills pool indices `study_block_ends(s-1) + 1` to
    !| `study_block_ends(s)` (from 1 for study 1). `pooled_residuals_perm` orders the whole pool
    !| ascending, NaN last.
    !|
    !| - [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]]:
    !|   calls
    !|   [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !|   with exactly these arguments, so the result is the occupancy search's, bit for bit. The
    !|   two per-study work arrays are not touched.
    !| - The three Stone modes: `group_sorted_perm_by_study` groups the sorted pool by block into
    !|   `tmp_residual_to_study_perm`, the block sizes go into `tmp_study_block_sizes`, and
    !|   `stone_bin_count_core` chooses the bin count. Its `stone_failed` is returned as
    !|   `bin_count_failed`. `gamma_occupancy` is not used.
    !|
    !| Both rules get `n_neighbors` only for the Sturges/Freedman-Diaconis diagnostics. Private
    !| helper of `determine_point_bin_count`.
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
            !! M_MODE_BIN_COUNT_OCCUPANCY

        integer(int32) :: i_study, n_study_neighbors, pool_offset, gather_ierr, diagnostic_n_neighbors
        integer(int32) :: actual_bin_count_mode

        call set_ok(ierr)
        M_DEFAULT_VAL(bin_count_mode, actual_bin_count_mode, M_MODE_BIN_COUNT_OCCUPANCY)

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

        ! point_counts/point_pmf passed whole: with one row, the callee's (1, max_n_bins) dummies are
        ! exactly their first max_n_bins elements (max_n_bins <= MAX_N_BINS by the precondition).
        call build_residual_histograms_impl(tmp_neighbor_residuals, max_n_reps_all_studies, n_neighbors, 1_int32, &
                                            point_range_low, point_range_high, max_n_bins, point_n_bins, &
                                            point_counts, point_pmf, &
                                            point_included_n_reps)

        counts(i_point, :) = point_counts(1, 1:max_n_bins)
        pmf(i_point, :) = point_pmf(1, 1:max_n_bins)
        included_n_reps = point_included_n_reps(1)
    end subroutine build_point_study_histogram

    !> summary: Build one study's fixed-k histograms for every reference point
    !| AUTHOR_LASZLO_LANG
    !| Pass C of the fixed-k pipeline for ONE study: gathers every reference point's neighbor
    !| residuals into `tmp_neighborhood_residuals_gathered`, then builds all of the study's
    !| histograms with one
    !| [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] call at
    !| the per-point bin counts. Internal helper, not published; public because
    !| `run_js_comp_test_impl` and the fixed-k parameter search call it. Precondition, not checked here: every index is in
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

    !> summary: Build every (reference point, study) histogram on adaptive neighborhoods
    !| AUTHOR_LASZLO_LANG
    !| Pass C on adaptive neighborhoods: decodes every reference point's pooled range into
    !| per-study gene lists
    !| ([[tox_data_integration_preprocessing_impl(module):materialize_pooled_neighborhood(subroutine)]])
    !| and builds each (point, study) histogram with `build_point_study_histogram` into bins
    !| `1:max_n_bins` of `counts`/`pmfs`. Sequential, so one gather buffer serves every pair.
    !| Internal helper, not published; public because `run_js_comp_test_adaptive_impl`
    !| (`n_bins_extent = M_MAX_N_BINS`) and the adaptive parameter search (`n_bins_extent =
    !| max_n_bins`, its buffers viewed packed) call it.
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

end module tox_data_integration_js_comp_test_binning_impl
