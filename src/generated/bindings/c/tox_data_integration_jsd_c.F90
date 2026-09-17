#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_jsd(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation
!|
!| The step that turns neighborhood residuals -- as produced by the preprocessing module -- into
!| a JSD value: residuals are binned into histograms per study, the two distributions are
!| compared per reference point, and the per-point divergences are weighted into one global
!| figure for the pair of studies.
module tox_data_integration_jsd_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: determine_shared_residual_range_c
    public :: determine_shared_residual_range_expert_c
    public :: determine_study_shared_residual_range_c
    public :: determine_study_shared_residual_range_expert_c
    public :: determine_all_studies_shared_residual_range_c
    public :: determine_all_studies_shared_residual_range_expert_c
    public :: build_residual_histograms_c
    public :: calc_pmf_c
    public :: compute_divergence_per_reference_point_c
    public :: compute_weighted_global_divergence_c

contains

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):determine_shared_residual_range(subroutine)]]
    !| This takes the pool already built; `determine_study_shared_residual_range` builds it from
    !| the neighborhood residuals of two studies first, if that is what is at hand.
    subroutine determine_shared_residual_range_c(&
            abs_residual_pool,&
            pool_size,&
            shared_residual_range,&
            residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_shared_residual_range_c")
        use tox_data_integration_jsd, only: determine_shared_residual_range

        integer(c_int), intent(in), target :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(c_double), dimension(pool_size), intent(in), target :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
            !! NaN is permitted for this value.
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(pool_size)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(abs_residual_pool, pool_size)

        call determine_shared_residual_range(&
            abs_residual_pool = abs_residual_pool,&
            pool_size = pool_size,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_shared_residual_range_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):determine_shared_residual_range_expert(subroutine)]]
    !| This takes the pool already built; `determine_study_shared_residual_range` builds it from
    !| the neighborhood residuals of two studies first, if that is what is at hand.
    subroutine determine_shared_residual_range_expert_c(&
            abs_residual_pool,&
            abs_residual_pool_perm,&
            pool_size,&
            shared_residual_range,&
            residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_shared_residual_range_expert_c")
        use tox_data_integration_jsd, only: determine_shared_residual_range_expert

        integer(c_int), intent(in), target :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(c_double), dimension(pool_size), intent(in), target :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
            !! NaN is permitted for this value.
        integer(c_int), dimension(pool_size), intent(in), target :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `pool_size`.
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(pool_size)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(abs_residual_pool, pool_size)
        M_CHECK_ARRAY_NON_NULL(abs_residual_pool_perm, pool_size)

        call determine_shared_residual_range_expert(&
            abs_residual_pool = abs_residual_pool,&
            abs_residual_pool_perm = abs_residual_pool_perm,&
            pool_size = pool_size,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_shared_residual_range_expert_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):determine_study_shared_residual_range(subroutine)]]
    !| Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
    !| as `determine_shared_residual_range` does.
    subroutine determine_study_shared_residual_range_c(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_study_shared_residual_range_c")
        use tox_data_integration_jsd, only: determine_study_shared_residual_range

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)

        call determine_study_shared_residual_range(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_study_shared_residual_range_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):determine_study_shared_residual_range_expert(subroutine)]]
    !| Pools the absolute residuals of both studies, sorts them, and takes the quantile exactly
    !| as `determine_shared_residual_range` does.
    subroutine determine_study_shared_residual_range_expert_c(&
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
        ) bind(C, name="determine_study_shared_residual_range_expert_c")
        use tox_data_integration_jsd, only: determine_study_shared_residual_range_expert

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension((n_reps_S1 + n_reps_S2)*n_neighbors*n_points), intent(out), target :: tmp_abs_residual_pool
            !! Work array holding the pooled absolute residuals of both studies
        integer(c_int), dimension((n_reps_S1 + n_reps_S2)*n_neighbors*n_points), intent(out), target :: tmp_abs_residual_pool_perm
            !! Work array for the permutation that sorts `tmp_abs_residual_pool`
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_abs_residual_pool, ((n_reps_S1 + n_reps_S2)*n_neighbors*n_points))
        M_CHECK_ARRAY_NON_NULL(tmp_abs_residual_pool_perm, ((n_reps_S1 + n_reps_S2)*n_neighbors*n_points))

        call determine_study_shared_residual_range_expert(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            tmp_abs_residual_pool = tmp_abs_residual_pool,&
            tmp_abs_residual_pool_perm = tmp_abs_residual_pool_perm,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_study_shared_residual_range_expert_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):determine_all_studies_shared_residual_range(subroutine)]]
    !| N-study generalization of `determine_study_shared_residual_range_impl`: pools the absolute
    !| residuals of every study, sorts them, and takes the quantile exactly as
    !| `determine_shared_residual_range` does.
    subroutine determine_all_studies_shared_residual_range_c(&
            neighborhood_residuals,&
            n_studies,&
            max_n_reps_all_studies,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_all_studies_shared_residual_range_c")
        use tox_data_integration_jsd, only: determine_all_studies_shared_residual_range

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        real(c_double), dimension(max_n_reps_all_studies, n_neighbors, n_points, n_studies), intent(in), target :: neighborhood_residuals
            !! Computed neighborhood residuals for every study, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals, max_n_reps_all_studies * n_neighbors * n_points * n_studies)

        call determine_all_studies_shared_residual_range(&
            neighborhood_residuals = neighborhood_residuals,&
            n_studies = n_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_all_studies_shared_residual_range_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):determine_all_studies_shared_residual_range_expert(subroutine)]]
    !| N-study generalization of `determine_study_shared_residual_range_impl`: pools the absolute
    !| residuals of every study, sorts them, and takes the quantile exactly as
    !| `determine_shared_residual_range` does.
    subroutine determine_all_studies_shared_residual_range_expert_c(&
            neighborhood_residuals,&
            n_studies,&
            max_n_reps_all_studies,&
            n_neighbors,&
            n_points,&
            tmp_abs_residual_pool,&
            tmp_abs_residual_pool_perm,&
            shared_residual_range,&
            residual_range_quantile,&
            ierr&
        ) bind(C, name="determine_all_studies_shared_residual_range_expert_c")
        use tox_data_integration_jsd, only: determine_all_studies_shared_residual_range_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        real(c_double), dimension(max_n_reps_all_studies, n_neighbors, n_points, n_studies), intent(in), target :: neighborhood_residuals
            !! Computed neighborhood residuals for every study, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies*n_neighbors*n_points*n_studies), intent(out), target :: tmp_abs_residual_pool
            !! Work array holding the pooled absolute residuals of every study
        integer(c_int), dimension(max_n_reps_all_studies*n_neighbors*n_points*n_studies), intent(out), target :: tmp_abs_residual_pool_perm
            !! Work array for the permutation that sorts `tmp_abs_residual_pool`
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile in [0,1] for determining the residual range
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
            !! The default value is `0.95_real64`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals, max_n_reps_all_studies * n_neighbors * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_abs_residual_pool, (max_n_reps_all_studies*n_neighbors*n_points*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_abs_residual_pool_perm, (max_n_reps_all_studies*n_neighbors*n_points*n_studies))

        call determine_all_studies_shared_residual_range_expert(&
            neighborhood_residuals = neighborhood_residuals,&
            n_studies = n_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            tmp_abs_residual_pool = tmp_abs_residual_pool,&
            tmp_abs_residual_pool_perm = tmp_abs_residual_pool_perm,&
            shared_residual_range = shared_residual_range,&
            residual_range_quantile = residual_range_quantile,&
            ierr = ierr&
        )
    end subroutine determine_all_studies_shared_residual_range_expert_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):build_residual_histograms(subroutine)]]
    !| The probability mass function `pmf(residual, bin)` is actually a matrix.
    subroutine build_residual_histograms_c(&
            neighborhood_residuals,&
            n_reps,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            max_n_bins,&
            n_bins_per_point,&
            counts,&
            pmf,&
            included_n_reps,&
            neighbor_mask,&
            ierr&
        ) bind(C, name="build_residual_histograms_c")
        use tox_data_integration_jsd, only: build_residual_histograms

        integer(c_int), intent(in), target :: n_reps
            !! Number of replicates of the study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of reference points (k)
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the study
        integer(c_int), intent(in), target :: max_n_bins
            !! Widest histogram bin count used by any reference point in this call -- the array
            !! extent `counts`/`pmf` are declared with. A reference point whose own
            !! `n_bins_per_point` is smaller has its remaining columns zero-padded.
        real(c_double), dimension(n_reps, n_neighbors, n_points), intent(in), target :: neighborhood_residuals
            !! Computed neighborhood residuals for a study, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R)
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_points), intent(in), target :: n_bins_per_point
            !! Number of equally sized histogram bins in range [-R,R] to use for this reference
            !! point
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_bins`.
        integer(c_int), dimension(n_points, max_n_bins), intent(out), target :: counts
            !! Absolute counts of a residual per bin, per reference point. Zero-padded beyond
            !! `n_bins_per_point(i_point)` for each reference point `i_point`
        real(c_double), dimension(n_points, max_n_bins), intent(out), target :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`.
            !! Zero-padded beyond `n_bins_per_point(i_point)` for each reference point `i_point`
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(max_n_bins)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals, n_reps * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(n_bins_per_point, n_points)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * max_n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf, n_points * max_n_bins)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points)

        call build_residual_histograms(&
            neighborhood_residuals = neighborhood_residuals,&
            n_reps = n_reps,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range = shared_residual_range,&
            max_n_bins = max_n_bins,&
            n_bins_per_point = n_bins_per_point,&
            counts = counts,&
            pmf = pmf,&
            included_n_reps = included_n_reps,&
            neighbor_mask = neighbor_mask,&
            ierr = ierr&
        )
    end subroutine build_residual_histograms_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):calc_pmf(subroutine)]]
    !| The counts-to-pmf half of `build_residual_histograms_impl`, factored out so a caller that
    !| already holds a study's `counts` -- untouched by anything that perturbed scratch copies
    !| downstream of it -- can re-derive `pmf` without re-binning residuals.
    subroutine calc_pmf_c(&
            counts,&
            included_n_reps,&
            n_points,&
            n_bins,&
            pmf,&
            ierr&
        ) bind(C, name="calc_pmf_c")
        use tox_data_integration_jsd, only: calc_pmf

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the study
        integer(c_int), intent(in), target :: n_bins
            !! The array's second extent for `counts`/`pmf` -- the widest histogram bin count of
            !! any reference point (`max_n_bins`, for a `counts` built by
            !! `build_residual_histograms_impl` from its own per-point `n_bins_per_point`).
            !! Columns beyond a given point's own bin count are zero-padded by that caller and are
            !! read/written here as ordinary zero entries
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: counts
            !! Absolute counts of a residual per bin
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf
            !! `counts` normalized to `0 <= pmf(:, i) <= 1` and `sum(pmf(:, i)) == 1`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(pmf, n_points * n_bins)

        call calc_pmf(&
            counts = counts,&
            included_n_reps = included_n_reps,&
            n_points = n_points,&
            n_bins = n_bins,&
            pmf = pmf,&
            ierr = ierr&
        )
    end subroutine calc_pmf_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):compute_divergence_per_reference_point(subroutine)]]
    !| Takes the probabilities `pmf` produced by `build_residual_histograms`. The natural JSD
    !| computed from the KL divergences lies in `[0, ln 2]`; the final step below rescales it by
    !| `LOG_2` onto `[0, 1]`.
    subroutine compute_divergence_per_reference_point_c(&
            pmf_S1,&
            pmf_S2,&
            n_points,&
            n_bins,&
            js_divergences,&
            ierr&
        ) bind(C, name="compute_divergence_per_reference_point_c")
        use tox_data_integration_jsd, only: compute_divergence_per_reference_point

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points (k)
        integer(c_int), intent(in), target :: n_bins
            !! The array's second extent for `pmf_S1`/`pmf_S2` -- the widest histogram bin count of
            !! any reference point (`max_n_bins`, for pmfs built by
            !! `build_residual_histograms_impl` from its own per-point `n_bins_per_point`). Both
            !! operands must have been built from the same per-point bin count for this to be
            !! safe: their zero-padded columns beyond that count then coincide, `S_mean` is zero
            !! there, and neither term contributes
        real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S1
            !! Computed normalized histogram counts for study 1
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S2
            !! Computed normalized histogram counts for study 2
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, rescaled by `LOG_2` onto `[0, 1]`
            !! (rather than its natural `[0, ln 2]` range)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)

        call compute_divergence_per_reference_point(&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            n_points = n_points,&
            n_bins = n_bins,&
            js_divergences = js_divergences,&
            ierr = ierr&
        )
    end subroutine compute_divergence_per_reference_point_c

    !> summary: C-wrapper for [[tox_data_integration_jsd(module):compute_weighted_global_divergence(subroutine)]]
    !| Takes the divergences produced by `compute_divergence_per_reference_point`.
    subroutine compute_weighted_global_divergence_c(&
            js_divergences,&
            n_points,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            global_js_divergence,&
            weights,&
            ierr&
        ) bind(C, name="compute_weighted_global_divergence_c")
        use tox_data_integration_jsd, only: compute_weighted_global_divergence

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points (k)
        real(c_double), dimension(n_points), intent(in), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2
            !! The minimum valid value is `0_int32`.
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)

        call compute_weighted_global_divergence(&
            js_divergences = js_divergences,&
            n_points = n_points,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            ierr = ierr&
        )
    end subroutine compute_weighted_global_divergence_c

end module tox_data_integration_jsd_c
#endif
