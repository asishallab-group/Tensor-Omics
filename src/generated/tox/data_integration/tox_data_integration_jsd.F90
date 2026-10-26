#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation
!|
!| The step that turns neighborhood residuals -- as produced by the preprocessing module -- into
!| a JSD value: residuals are binned into histograms per study, the two distributions are
!| compared per reference point, and the per-point divergences are weighted into one global
!| figure for the pair of studies. For K studies, each study is compared against the consensus pmf
!| of all of them, built by
!| [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]] (with its pooled
!| counts) or [[tox_data_integration_jsd_impl(module):create_mean_pmf_only_impl(interface)]].
!|
!| Generated from [[tox_data_integration_jsd_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_jsd
    use f42_safeguard
    use tox_data_integration_jsd_impl, only: build_residual_histograms_impl, calc_pmf_impl, compute_divergence_per_reference_point_impl, compute_weighted_global_divergence_impl
    use tox_data_integration_jsd_impl, only: create_mean_pmf_impl, create_mean_pmf_only_impl, determine_shared_residual_range_impl, determine_study_shared_residual_range_impl
    use tox_data_integration_jsd_impl, only: jct_compute_jsd_pipeline_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: determine_shared_residual_range
    public :: determine_shared_residual_range_expert
    public :: determine_study_shared_residual_range
    public :: determine_study_shared_residual_range_expert
    public :: build_residual_histograms
    public :: calc_pmf
    public :: compute_divergence_per_reference_point
    public :: compute_weighted_global_divergence
    public :: jct_compute_jsd_pipeline
    public :: jct_compute_jsd_pipeline_expert
    public :: create_mean_pmf
    public :: create_mean_pmf_only

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_jsd_impl(module):determine_shared_residual_range_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_jsd(module):determine_shared_residual_range_expert]] to prepare it yourself.
    !| This takes the pool already built; `determine_study_shared_residual_range` builds it from
    !| the neighborhood residuals of two studies first, if that is what is at hand.
    pure subroutine determine_shared_residual_range(&
            abs_residual_pool,&
            pool_size,&
            shared_residual_range,&
            residual_range_quantile_level,&
            ierr&
        )
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
            !! NaN is permitted for this value.
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile_level
            !! Quantile level in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: abs_residual_pool_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(pool_size, ierr, arg_pos=2_int32)
        call validate_in_range_real(residual_range_quantile_level, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(abs_residual_pool, pool_size, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(abs_residual_pool_perm(pool_size))
        call init_perm(abs_residual_pool_perm)
        call sort_array_heapsort(abs_residual_pool, abs_residual_pool_perm)

        call determine_shared_residual_range_impl(&
            abs_residual_pool = abs_residual_pool,&
            abs_residual_pool_perm = abs_residual_pool_perm,&
            pool_size = pool_size,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile_level = residual_range_quantile_level&
        )
    end subroutine determine_shared_residual_range

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):determine_shared_residual_range_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_jsd(module):determine_shared_residual_range]] does both.
    !| This takes the pool already built; `determine_study_shared_residual_range` builds it from
    !| the neighborhood residuals of two studies first, if that is what is at hand.
    pure subroutine determine_shared_residual_range_expert(&
            abs_residual_pool,&
            abs_residual_pool_perm,&
            pool_size,&
            shared_residual_range,&
            residual_range_quantile_level,&
            ierr&
        )
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
            !! NaN is permitted for this value.
        integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `pool_size`.
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile_level
            !! Quantile level in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(pool_size, ierr, arg_pos=3_int32)
        call validate_in_range_real(residual_range_quantile_level, ierr, arg_pos=5_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(abs_residual_pool, pool_size, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_int(abs_residual_pool_perm, pool_size, ierr, arg_pos=2_int32, min=1_int32, max=pool_size)
        if (is_err(ierr)) return
#endif

        call determine_shared_residual_range_impl(&
            abs_residual_pool = abs_residual_pool,&
            abs_residual_pool_perm = abs_residual_pool_perm,&
            pool_size = pool_size,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile_level = residual_range_quantile_level&
        )
    end subroutine determine_shared_residual_range_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_jsd_impl(module):determine_study_shared_residual_range_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_jsd(module):determine_study_shared_residual_range_expert]] to prepare it yourself.
    !| Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
    !| as `determine_shared_residual_range` does.
    pure subroutine determine_study_shared_residual_range(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            residual_range_quantile_level,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile_level
            !! Quantile level in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_abs_residual_pool
        integer(int32), dimension(:), allocatable :: tmp_abs_residual_pool_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_in_range_real(residual_range_quantile_level, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_abs_residual_pool((n_reps_S1 + n_reps_S2)*n_neighbors*n_points))
        M_ALLOCATE(tmp_abs_residual_pool_perm((n_reps_S1 + n_reps_S2)*n_neighbors*n_points))

        call determine_study_shared_residual_range_impl(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            tmp_abs_residual_pool = tmp_abs_residual_pool,&
            tmp_abs_residual_pool_perm = tmp_abs_residual_pool_perm,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile_level = residual_range_quantile_level&
        )
    end subroutine determine_study_shared_residual_range

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):determine_study_shared_residual_range_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_jsd(module):determine_study_shared_residual_range]] does both.
    !| Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
    !| as `determine_shared_residual_range` does.
    pure subroutine determine_study_shared_residual_range_expert(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            tmp_abs_residual_pool,&
            tmp_abs_residual_pool_perm,&
            shared_residual_range,&
            residual_range_quantile_level,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension((n_reps_S1 + n_reps_S2)*n_neighbors*n_points), intent(out) :: tmp_abs_residual_pool
            !! Work array holding the pooled absolute residuals of both studies
        integer(int32), dimension((n_reps_S1 + n_reps_S2)*n_neighbors*n_points), intent(out) :: tmp_abs_residual_pool_perm
            !! Work array for the permutation that sorts `tmp_abs_residual_pool`
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile_level
            !! Quantile level in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_in_range_real(residual_range_quantile_level, ierr, arg_pos=10_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        call determine_study_shared_residual_range_impl(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            tmp_abs_residual_pool = tmp_abs_residual_pool,&
            tmp_abs_residual_pool_perm = tmp_abs_residual_pool_perm,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile_level = residual_range_quantile_level&
        )
    end subroutine determine_study_shared_residual_range_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl]].
    !| The probability mass function `pmf(residual, bin)` is actually a matrix.
    pure subroutine build_residual_histograms(&
            neighborhood_residuals,&
            n_reps,&
            n_neighbors,&
            n_points,&
            shared_residual_range_low,&
            shared_residual_range_high,&
            max_n_bins,&
            n_bins_per_point,&
            counts,&
            pmf,&
            included_n_reps,&
            neighbor_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps
            !! Number of replicates of the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the study
        integer(int32), intent(in) :: max_n_bins
            !! Widest histogram bin count used by any reference point in this call -- the array
            !! extent `counts`/`pmf` are declared with. A reference point whose own
            !! `n_bins_per_point` is smaller has its remaining columns zero-padded.
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_points), intent(in) :: shared_residual_range_low
            !! Lower bound of the histogram range (R_low) for reference point i_point -- e.g. from
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
            !! own `shared_residual_range_low` output
        real(real64), dimension(n_points), intent(in) :: shared_residual_range_high
            !! Upper bound of the histogram range (R_high) for reference point i_point -- e.g. from
            !! [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
            !! own `shared_residual_range_high` output
        integer(int32), dimension(n_points), intent(in) :: n_bins_per_point
            !! Number of equally sized histogram bins in range [shared_residual_range_low(i_point),
            !! shared_residual_range_high(i_point)] to use for this reference point
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_bins`.
        integer(int32), dimension(n_points, max_n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin, per reference point. Zero-padded beyond
            !! `n_bins_per_point(i_point)` for each reference point `i_point`
        real(real64), dimension(n_points, max_n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`.
            !! Zero-padded beyond `n_bins_per_point(i_point)` for each reference point `i_point`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_dimension_size(max_n_bins, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(neighborhood_residuals, n_reps * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(shared_residual_range_low, n_points, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(shared_residual_range_high, n_points, ierr, arg_pos=6_int32)
        call validate_all_in_range_int(n_bins_per_point, n_points, ierr, arg_pos=8_int32, min=1_int32, max=max_n_bins)
        if (is_err(ierr)) return
#endif

        call build_residual_histograms_impl(&
            neighborhood_residuals = neighborhood_residuals,&
            n_reps = n_reps,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range_low = shared_residual_range_low,&
            shared_residual_range_high = shared_residual_range_high,&
            max_n_bins = max_n_bins,&
            n_bins_per_point = n_bins_per_point,&
            counts = counts,&
            pmf = pmf,&
            included_n_reps = included_n_reps,&
            neighbor_mask = neighbor_mask&
        )
    end subroutine build_residual_histograms

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):calc_pmf_impl]].
    !| The counts-to-pmf half of `build_residual_histograms_impl`, factored out so a caller that
    !| already holds a study's `counts` -- untouched by anything that perturbed scratch copies
    !| downstream of it -- can re-derive `pmf` without re-binning residuals.
    pure subroutine calc_pmf(&
            counts,&
            included_n_reps,&
            n_points,&
            n_bins,&
            pmf,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the study
        integer(int32), intent(in) :: n_bins
            !! The array's second extent for `counts`/`pmf` -- the widest histogram bin count of
            !! any reference point (`max_n_bins`, for a `counts` built by
            !! `build_residual_histograms_impl` from its own per-point `n_bins_per_point`).
            !! Columns beyond a given point's own bin count are zero-padded by that caller and are
            !! read/written here as ordinary zero entries
        integer(int32), dimension(n_points, n_bins), intent(in) :: counts
            !! Absolute counts of a residual per bin
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= pmf(:, i) <= 1` and `sum(pmf(:, i)) == 1`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=4_int32)
        call validate_all_in_range_int(counts, n_points * n_bins, ierr, arg_pos=1_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points, ierr, arg_pos=2_int32, min=0_int32)
        if (is_err(ierr)) return
#endif

        call calc_pmf_impl(&
            counts = counts,&
            included_n_reps = included_n_reps,&
            n_points = n_points,&
            n_bins = n_bins,&
            pmf = pmf&
        )
    end subroutine calc_pmf

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl]].
    !| Takes the probabilities `pmf` produced by `build_residual_histograms`. The natural JSD
    !| computed from the KL divergences lies in `[0, ln 2]`; the final step below rescales it by
    !| `LOG_2` onto `[0, 1]`.
    pure subroutine compute_divergence_per_reference_point(&
            pmf_S1,&
            pmf_S2,&
            n_points,&
            n_bins,&
            js_divergences,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        integer(int32), intent(in) :: n_bins
            !! The array's second extent for `pmf_S1`/`pmf_S2` -- the widest histogram bin count of
            !! any reference point (`max_n_bins`, for pmfs built by
            !! `build_residual_histograms_impl` from its own per-point `n_bins_per_point`). Both
            !! operands must have been built from the same per-point bin count for this to be
            !! safe: their zero-padded columns beyond that count then coincide, `S_mean` is zero
            !! there, and neither term contributes
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized histogram counts for study 1
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
            !! Computed normalized histogram counts for study 2
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, rescaled by `LOG_2` onto `[0, 1]`
            !! (rather than its natural `[0, ln 2]` range)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(pmf_S1, n_points * n_bins, ierr, arg_pos=1_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(pmf_S2, n_points * n_bins, ierr, arg_pos=2_int32, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return
#endif

        call compute_divergence_per_reference_point_impl(&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            n_points = n_points,&
            n_bins = n_bins,&
            js_divergences = js_divergences&
        )
    end subroutine compute_divergence_per_reference_point

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl]].
    !| Takes the divergences produced by `compute_divergence_per_reference_point`.
    pure subroutine compute_weighted_global_divergence(&
            js_divergences,&
            n_points,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            global_js_divergence,&
            weights,&
            ierr&
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points (k)
        real(real64), dimension(n_points), intent(in) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
            !! The minimum valid value is `0_int32`.
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(js_divergences, n_points, ierr, arg_pos=1_int32, min=0.0_real64)
        call validate_all_in_range_int(included_n_reps_S1, n_points, ierr, arg_pos=3_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps_S2, n_points, ierr, arg_pos=4_int32, min=0_int32)
        if (is_err(ierr)) return
#endif

        call compute_weighted_global_divergence_impl(&
            js_divergences = js_divergences,&
            n_points = n_points,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            global_js_divergence = global_js_divergence,&
            weights = weights&
        )
    end subroutine compute_weighted_global_divergence

    !> summary: Validates its inputs, prepares what [[tox_data_integration_jsd_impl(module):jct_compute_jsd_pipeline_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_jsd(module):jct_compute_jsd_pipeline_expert]] to prepare it yourself.
    !| Runs [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] for
    !| each study, [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]
    !| and [[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl(interface)]],
    !| with one bin count `n_bins` and one symmetric range `[-shared_residual_range, shared_residual_range]`
    !| for every reference point.
    pure subroutine jct_compute_jsd_pipeline(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            global_js_divergence,&
            weights,&
            pmf_S1,&
            pmf_S2,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Normalized histogram counts for study 1
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Normalized histogram counts for study 2
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        integer(int32), dimension(:), allocatable :: tmp_n_bins_per_point
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_low
        real(real64), dimension(:), allocatable :: tmp_shared_residual_range_high

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=7_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=8_int32, min=0.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_counts(n_points, n_bins))
        M_ALLOCATE(tmp_n_bins_per_point(n_points))
        M_ALLOCATE(tmp_shared_residual_range_low(n_points))
        M_ALLOCATE(tmp_shared_residual_range_high(n_points))

        call jct_compute_jsd_pipeline_impl(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            tmp_counts = tmp_counts,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            neighbor_mask_S1 = neighbor_mask_S1,&
            neighbor_mask_S2 = neighbor_mask_S2&
        )
    end subroutine jct_compute_jsd_pipeline

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):jct_compute_jsd_pipeline_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_jsd(module):jct_compute_jsd_pipeline]] does both.
    !| Runs [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]] for
    !| each study, [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl(interface)]]
    !| and [[tox_data_integration_jsd_impl(module):compute_weighted_global_divergence_impl(interface)]],
    !| with one bin count `n_bins` and one symmetric range `[-shared_residual_range, shared_residual_range]`
    !| for every reference point.
    pure subroutine jct_compute_jsd_pipeline_expert(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            global_js_divergence,&
            weights,&
            pmf_S1,&
            pmf_S2,&
            tmp_counts,&
            tmp_n_bins_per_point,&
            tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            ierr&
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
        integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
        real(real64), intent(out) :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(real64), dimension(n_points), intent(out) :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
            !! Normalized histogram counts for study 1
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
            !! Normalized histogram counts for study 2
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for the histogram counts
        integer(int32), dimension(n_points), intent(out) :: tmp_n_bins_per_point
            !! Working array: `n_bins` broadcast to every reference point, since this routine uses
            !! one uniform bin count -- passed on to
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]'s
            !! own `n_bins_per_point`
        real(real64), dimension(n_points), intent(out) :: tmp_shared_residual_range_low
            !! Working array: `-shared_residual_range` broadcast to every reference point, since
            !! this routine uses one uniform symmetric range
            !! `[-shared_residual_range, shared_residual_range]` -- passed on to
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]'s
            !! own `shared_residual_range_low`
        real(real64), dimension(n_points), intent(out) :: tmp_shared_residual_range_high
            !! Working array: `shared_residual_range` broadcast to every reference point, mirroring
            !! `tmp_shared_residual_range_low` above -- passed on to
            !! [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl(interface)]]'s
            !! own `shared_residual_range_high`
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=7_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=8_int32, min=0.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true._c_bool)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true._c_bool)
        if (is_err(ierr)) return
#endif

        call jct_compute_jsd_pipeline_impl(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            tmp_counts = tmp_counts,&
            tmp_n_bins_per_point = tmp_n_bins_per_point,&
            tmp_shared_residual_range_low = tmp_shared_residual_range_low,&
            tmp_shared_residual_range_high = tmp_shared_residual_range_high,&
            neighbor_mask_S1 = neighbor_mask_S1,&
            neighbor_mask_S2 = neighbor_mask_S2&
        )
    end subroutine jct_compute_jsd_pipeline_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl]].
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

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):create_mean_pmf_only_impl]].
    !| As [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]], for
    !| callers that need no counts (the bootstrap and the permutation test). Same known limitation: the consensus includes
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

end module tox_data_integration_jsd
