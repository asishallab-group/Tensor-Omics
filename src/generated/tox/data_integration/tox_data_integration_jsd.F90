#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation
!|
!| The step that turns neighborhood residuals -- as produced by the preprocessing module -- into
!| a JSD value: residuals are binned into histograms per study, the two distributions are
!| compared per reference point, and the per-point divergences are weighted into one global
!| figure for the pair of studies.
!|
!| Generated from [[tox_data_integration_jsd_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_jsd
    use tox_data_integration_jsd_impl, only: build_residual_histograms_impl, compute_divergence_per_reference_point_impl, compute_weighted_global_divergence_impl, determine_shared_residual_range_impl
    use tox_data_integration_jsd_impl, only: determine_study_shared_residual_range_impl
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
    public :: compute_divergence_per_reference_point
    public :: compute_weighted_global_divergence

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_jsd_impl(module):determine_shared_residual_range_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_jsd(module):determine_shared_residual_range_expert]] to prepare it yourself.
    !| This takes the pool already built; `determine_study_shared_residual_range` builds it from
    !| the neighborhood residuals of two studies first, if that is what is at hand.
    pure subroutine determine_shared_residual_range(&
            abs_residual_pool,&
            pool_size,&
            shared_residual_range,&
            residual_range_quantile,&
            ierr&
        )
        integer(int32), intent(in) :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
            !! NaN is permitted for this value.
        real(real64), intent(out) :: shared_residual_range
            !! Computed residual range (R)
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: abs_residual_pool_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(pool_size, ierr, arg_pos=2_int32)
        call validate_in_range_real(residual_range_quantile, ierr, arg_pos=4_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(abs_residual_pool, pool_size, ierr, arg_pos=1_int32, allow_nan=.true.)
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
            residual_range_quantile = residual_range_quantile&
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
            residual_range_quantile,&
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
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(pool_size, ierr, arg_pos=3_int32)
        call validate_in_range_real(residual_range_quantile, ierr, arg_pos=5_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(abs_residual_pool, pool_size, ierr, arg_pos=1_int32, allow_nan=.true.)
        call validate_all_in_range_int(abs_residual_pool_perm, pool_size, ierr, arg_pos=2_int32, min=1_int32, max=pool_size)
        if (is_err(ierr)) return
#endif

        call determine_shared_residual_range_impl(&
            abs_residual_pool = abs_residual_pool,&
            abs_residual_pool_perm = abs_residual_pool_perm,&
            pool_size = pool_size,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile&
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
            residual_range_quantile,&
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
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95`.
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
        call validate_in_range_real(residual_range_quantile, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true.)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true.)
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
            residual_range_quantile = residual_range_quantile&
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
            residual_range_quantile,&
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
        real(real64), intent(in), optional :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_in_range_real(residual_range_quantile, ierr, arg_pos=10_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true.)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true.)
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
            residual_range_quantile = residual_range_quantile&
        )
    end subroutine determine_study_shared_residual_range_expert

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):build_residual_histograms_impl]].
    !| The probability mass function `pmf(residual, bin)` is actually a matrix.
    pure subroutine build_residual_histograms(&
            neighborhood_residuals,&
            n_reps,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            n_bins,&
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
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
            !! Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_points, n_bins), intent(out) :: counts
            !! Absolute counts of a residual per bin
        real(real64), dimension(n_points, n_bins), intent(out) :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(int32), dimension(n_points), intent(out) :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_dimension_size(n_bins, ierr, arg_pos=6_int32)
        call validate_all_in_range_real(neighborhood_residuals, n_reps * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        call build_residual_histograms_impl(&
            neighborhood_residuals = neighborhood_residuals,&
            n_reps = n_reps,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            counts = counts,&
            pmf = pmf,&
            included_n_reps = included_n_reps,&
            neighbor_mask = neighbor_mask&
        )
    end subroutine build_residual_histograms

    !> summary: Validates its inputs, then calls [[tox_data_integration_jsd_impl(module):compute_divergence_per_reference_point_impl]].
    !| Takes the probabilities `pmf` produced by `build_residual_histograms`.
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
            !! Number of equally sized histogram bins in range [-R,R]
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
            !! Computed normalized histogram counts for study 1
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
            !! Computed normalized histogram counts for study 2
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), dimension(n_points), intent(out) :: js_divergences
            !! Jensen-Shannon divergence per reference point
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

end module tox_data_integration_jsd
