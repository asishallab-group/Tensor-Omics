#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test
!|
!| Kernel for a permutation test that estimates an empirical p-value for the weighted global JSD
!| computed by the JSD kernels. Under the null hypothesis that both studies are exchangeable,
!| S1/S2 residuals are repeatedly shuffled within each reference point's pooled neighborhood and
!| the JSD pipeline is recomputed, giving a null distribution against which the observed JSD is
!| compared. The generator turns `gjct_permutation_test_kernel` into the validating and
!| allocating wrappers in module tox_data_integration_stats.
module tox_data_integration_stats_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_utils, only: init_random, shuffle_vector
    use tox_data_integration_jsd_kernel, only: jct_compute_jsd_pipeline_helper
    M_IMPLICIT_NONE
contains

    !> summary: Estimate how likely the observed divergence is to occur by chance
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
    !| the work copies, so the caller's own arrays are left untouched.
    subroutine gjct_permutation_test_kernel( &
        neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
        tmp_residuals_S1, tmp_residuals_S2, tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
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
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! DM_ALLOW_NAN
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(out), target :: tmp_residuals_S1
            !! Work copy of study 1's residuals, shuffled in place
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(out), target :: tmp_residuals_S2
            !! Work copy of study 2's residuals, shuffled in place
        real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out), target :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
            !! Working array for study 1's normalized histogram counts
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
            !! Working array for study 2's normalized histogram counts
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for the histogram counts
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
            !! Working array for study 1's included replicate counts
        integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
            !! Working array for study 2's included replicate counts
        real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
            !! Working array for the per-reference-point divergences
        real(real64), dimension(n_points), intent(out) :: tmp_weights
            !! Working array for the divergence weights
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

        ! The residuals are shuffled in place, so work on copies and leave the caller's arrays alone
        tmp_residuals_S1 = neighborhood_residuals_S1
        tmp_residuals_S2 = neighborhood_residuals_S2

        pool_flat(1:size(tmp_pool, kind=int32)) => tmp_pool
        n_residuals_S1 = n_reps_S1*n_neighbors
        n_residuals_S2 = n_reps_S2*n_neighbors

        n_jsd_exceeding_observed = 0_int32
        do i_permutation = 1, n_permutations
            ! 1. shuffle residuals
            do i_point = 1, n_points
                reference_point_S1(1:n_residuals_S1) => tmp_residuals_S1(:, :, i_point)
                reference_point_S2(1:n_residuals_S2) => tmp_residuals_S2(:, :, i_point)
                call shuffle_reference_point_helper( &
                    reference_point_S1, reference_point_S2, &
                    n_reps_S1, n_reps_S2, n_neighbors, pool_flat &
                    )
            end do

            ! 2. Pipeline to determine the global jsd for current permutation
            call jct_compute_jsd_pipeline_helper(tmp_residuals_S1, tmp_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, tmp_js_divergences, tmp_included_n_reps_S1, tmp_included_n_reps_S2, jsd_null(i_permutation), tmp_weights, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)

            if (jsd_null(i_permutation) >= global_jsd_observed) then
                n_jsd_exceeding_observed = n_jsd_exceeding_observed + 1
            end if
        end do

        ! Add-one (Laplace) smoothing on both the numerator and denominator: this treats the observed
        ! statistic itself as one additional permutation draw, so the p-value can never be exactly zero
        ! (which would otherwise happen whenever no null draw reaches the observed JSD).
        p_value = real(n_jsd_exceeding_observed + 1, real64)/real(n_permutations + 1, real64)
    end subroutine gjct_permutation_test_kernel

    !> summary: Shuffle the residuals of one reference point between the two studies
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Internal helper of `gjct_permutation_test_kernel`.
    subroutine shuffle_reference_point_helper(reference_point_S1, reference_point_S2, n_reps_S1, n_reps_S2, n_neighbors, pool_flat)
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        real(real64), dimension(n_reps_S1*n_neighbors), intent(inout) :: reference_point_S1
            !! Residuals for one reference point in study 1, will be shuffled in-place
        real(real64), dimension(n_reps_S2*n_neighbors), intent(inout) :: reference_point_S2
            !! Residuals for one reference point in study 2, will be shuffled in-place
        real(real64), dimension((n_reps_S1 + n_reps_S2)*n_neighbors), intent(out) :: pool_flat
            !! Working array for shuffling the concatenated residuals from both studies per reference point

        integer(int32) :: pool_size, n_residuals_S1

        pool_size = size(pool_flat, kind=int32)
        n_residuals_S1 = size(reference_point_S1, kind=int32)

        pool_flat(1:n_residuals_S1) = reference_point_S1
        pool_flat(n_residuals_S1 + 1:pool_size) = reference_point_S2

        call shuffle_vector(pool_flat)

        reference_point_S1 = pool_flat(1:n_residuals_S1)
        reference_point_S2 = pool_flat(n_residuals_S1 + 1:pool_size)
    end subroutine shuffle_reference_point_helper

end module tox_data_integration_stats_kernel
