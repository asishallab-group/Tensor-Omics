#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| This module implements the pipeline to obtain neighborhood residuals from expression vectors, to be used for JCT based data integration.
submodule (tox_data_integration) tox_data_integration_stats
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: validate_dimension_size, validate_in_range_int, validate_in_range_real, set_ok, set_err, is_err, ERR_ALLOC_FAIL
    use f42_utils, only: init_random, shuffle_vector
contains
    
    !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
    module subroutine gjct_permutation_test_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        real(real64), dimension(:, :, :), allocatable :: S1
        real(real64), dimension(:, :, :), allocatable :: S2
        real(real64), dimension(:, :), allocatable :: tmp_pool
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S1
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S2
        integer(int32), dimension(:), allocatable :: tmp_included_n_reps_S1
        integer(int32), dimension(:), allocatable :: tmp_included_n_reps_S2
        real(real64), dimension(:), allocatable :: tmp_js_divergences
        real(real64), dimension(:), allocatable :: tmp_weights

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)

        if (is_err(ierr)) return

        ! Allocate working arrays
        M_ALLOCATE(tmp_pool(n_reps_S1 + n_reps_S2, n_neighbors))
        M_ALLOCATE(tmp_counts(n_points, n_bins))
        M_ALLOCATE(tmp_pmf_S1(n_points, n_bins))
        M_ALLOCATE(tmp_pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_included_n_reps_S1(n_points))
        M_ALLOCATE(tmp_included_n_reps_S2(n_points))
        M_ALLOCATE(tmp_js_divergences(n_points))
        M_ALLOCATE(tmp_weights(n_points))

        ! Allocate the residual copies
        M_ALLOCATE(S1(n_reps_S1, n_neighbors, n_points))
        M_ALLOCATE(S2(n_reps_S2, n_neighbors, n_points))

        S1 = neighborhood_residuals_S1
        S2 = neighborhood_residuals_S2

        call gjct_permutation_test(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    end subroutine gjct_permutation_test_alloc

    !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
    module subroutine gjct_permutation_test( &
            neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
            tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
            ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2 &
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S1_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S2_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
        real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out) :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
            !! Working array for [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
        real(real64), dimension(n_points), intent(out) :: tmp_weights
            !! Working array for [[tox_data_integration(module):compute_weighted_global_divergence(interface)]]
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        call set_ok(ierr)

        call validate_dimension_size(n_reps_S1, ierr)
        call validate_dimension_size(n_reps_S2, ierr)
        call validate_dimension_size(n_neighbors, ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_in_range_int(n_bins, ierr, min=1_int32)
        call validate_in_range_int(n_permutations, ierr, min=1_int32)
        call validate_in_range_real(shared_residual_range, ierr, min=0.0_real64)

        if (is_err(ierr)) return

        call gjct_permutation_test_helper(neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, random_seed, neighbor_mask_S1, neighbor_mask_S2)
    end subroutine gjct_permutation_test

    !> (no input validation) Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
    module subroutine gjct_permutation_test_helper( &
            neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
            tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
            random_seed, neighbor_mask_S1, neighbor_mask_S2 &
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S1_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S2_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
        real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out), target :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
            !! Working array for [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
        real(real64), dimension(n_points), intent(out) :: tmp_weights
            !! Working array for [[tox_data_integration(module):compute_weighted_global_divergence(interface)]]
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

        integer(int32) :: n_jsd_exceeding_observed, i_permutation, n_residuals_S1, n_residuals_S2, i_point
        real(real64), dimension(:), pointer :: reference_point_S1, reference_point_S2, pool_flat

        if (present(random_seed)) then
            call init_random(random_seed)
        end if

        pool_flat(1:size(tmp_pool, kind=int32)) => tmp_pool
        n_residuals_S1 = n_reps_S1 * n_neighbors
        n_residuals_S2 = n_reps_S2 * n_neighbors

        n_jsd_exceeding_observed = 0_int32
        do i_permutation = 1, n_permutations
            ! 1. shuffle residuals
            do i_point = 1, n_points
                reference_point_S1(1:n_residuals_S1) => neighborhood_residuals_S1_copy(:, :, i_point)
                reference_point_S2(1:n_residuals_S2) => neighborhood_residuals_S2_copy(:, :, i_point)
                call shuffle_reference_point_helper( &
                    reference_point_S1, reference_point_S2, &
                    n_reps_S1, n_reps_S2, n_neighbors, pool_flat &
                )
            end do

            ! 2. Pipeline to determine the global jsd for current permutation
            call jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, tmp_js_divergences, tmp_included_n_reps_S1, tmp_included_n_reps_S2, jsd_null(i_permutation), tmp_weights, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)

            if (jsd_null(i_permutation) >= global_jsd_observed) then
                n_jsd_exceeding_observed = n_jsd_exceeding_observed + 1
            end if
        end do

        p_value = real(n_jsd_exceeding_observed + 1, real64) / real(n_permutations + 1, real64)
    end subroutine gjct_permutation_test_helper

    !> Helper for [[tox_data_integration(module):gjct_permutation_test_helper(interface)]] to shuffle reference points
    module subroutine shuffle_reference_point_helper(reference_point_S1, reference_point_S2, n_reps_S1, n_reps_S2, n_neighbors, pool_flat)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        real(real64), dimension(n_reps_S1 * n_neighbors), intent(inout) :: reference_point_S1
            !! Residuals for one reference point in study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2 * n_neighbors), intent(inout) :: reference_point_S2
            !! Residuals for one reference point in study 2, will be shuffled in-place
        real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors), intent(out) :: pool_flat
            !! Working array for shuffling the concatenated residuals from both studies per reference point

        integer(int32) :: pool_size, n_residuals_S1

        pool_size = size(pool_flat, kind=int32)
        n_residuals_S1 = size(reference_point_S1, kind=int32)

        pool_flat(1:n_residuals_S1) = reference_point_S1
        pool_flat(n_residuals_S1+1:pool_size) = reference_point_S2

        call shuffle_vector(pool_flat)

        reference_point_S1 = pool_flat(1:n_residuals_S1)
        reference_point_S2 = pool_flat(n_residuals_S1+1:pool_size)
    end subroutine shuffle_reference_point_helper

end submodule tox_data_integration_stats

!> C-compatible wrapper for [[tox_data_integration(module):gjct_permutation_test_alloc(interface)]]
subroutine gjct_permutation_test_c( &
    neighborhood_residuals_S1, neighborhood_residuals_S2, &
    n_reps_S1, n_reps_S2, n_neighbors, n_points, &
    global_jsd_observed, n_bins, shared_residual_range, n_permutations, &
    jsd_null, p_value, ierr, random_seed) &
    bind(C, name="gjct_permutation_test_c")

    use tox_data_integration, only: gjct_permutation_test_alloc
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_reps_S1
        !! Number of replicates in study 1
    integer(c_int), intent(in), target :: n_reps_S2
        !! Number of replicates in study 2
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors in the studies
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points in the studies
    real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
        !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
        !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
    real(c_double), intent(in), target :: global_jsd_observed
        !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
    integer(c_int), intent(in), target :: n_bins
        !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
    real(c_double), intent(in), target :: shared_residual_range
        !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
    integer(c_int), intent(in), target :: n_permutations
        !! Number of permutations to perform
    real(c_double), dimension(n_permutations), intent(out), target :: jsd_null
        !! Vector of global divergence values obtained under the null hypothesis
    real(c_double), intent(out), target :: p_value
        !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
    integer(c_int), intent(out), target :: ierr
        !! Error code
    integer(c_int), intent(in), target :: random_seed
        !! Seed to use for shuffling

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_reps_S1)
    M_CHECK_NON_NULL(n_reps_S2)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(neighborhood_residuals_S1)
    M_CHECK_NON_NULL(neighborhood_residuals_S2)
    M_CHECK_NON_NULL(global_jsd_observed)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_permutations)
    M_CHECK_NON_NULL(jsd_null)
    M_CHECK_NON_NULL(p_value)
    M_CHECK_NON_NULL(random_seed)

    call gjct_permutation_test_alloc( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, &
        n_reps_S1, n_reps_S2, n_neighbors, n_points, &
        global_jsd_observed, n_bins, shared_residual_range, n_permutations, &
        jsd_null, p_value, ierr, random_seed)

end subroutine gjct_permutation_test_c
