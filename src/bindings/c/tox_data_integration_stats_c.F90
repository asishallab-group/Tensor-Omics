#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_stats(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_data_integration_stats_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: gjct_permutation_test_expert_c
    public :: gjct_permutation_test_c

contains

    !> summary: C-wrapper for [[tox_data_integration_stats(module):gjct_permutation_test(subroutine)]]
    !| Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
    !| the work copies, so the caller's own arrays are left untouched.
    subroutine gjct_permutation_test_expert_c(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            global_jsd_observed,&
            n_bins,&
            shared_residual_range,&
            n_permutations,&
            jsd_null,&
            p_value,&
            tmp_residuals_S1,&
            tmp_residuals_S2,&
            tmp_pool,&
            tmp_pmf_S1,&
            tmp_pmf_S2,&
            tmp_counts,&
            tmp_included_n_reps_S1,&
            tmp_included_n_reps_S2,&
            tmp_js_divergences,&
            tmp_weights,&
            random_seed,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            ierr&
        ) bind(C, name="gjct_permutation_test_expert_c")
        use tox_data_integration_stats, only: gjct_permutation_test

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations to perform
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), intent(in), target :: global_jsd_observed
            !! Observed global JSD value for both studies
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_permutations), intent(out), target :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(c_double), intent(out), target :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(out), target :: tmp_residuals_S1
            !! Work copy of study 1's residuals, shuffled in place
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(out), target :: tmp_residuals_S2
            !! Work copy of study 2's residuals, shuffled in place
        real(c_double), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out), target :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_S1
            !! Working array for study 1's normalized histogram counts
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_S2
            !! Working array for study 2's normalized histogram counts
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Working array for the histogram counts
        integer(c_int), dimension(n_points), intent(out), target :: tmp_included_n_reps_S1
            !! Working array for study 1's included replicate counts
        integer(c_int), dimension(n_points), intent(out), target :: tmp_included_n_reps_S2
            !! Working array for study 2's included replicate counts
        real(c_double), dimension(n_points), intent(out), target :: tmp_js_divergences
            !! Working array for the per-reference-point divergences
        real(c_double), dimension(n_points), intent(out), target :: tmp_weights
            !! Working array for the divergence weights
        integer(c_int), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(:, :), allocatable :: neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(global_jsd_observed)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(p_value)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(jsd_null, n_permutations)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_pool, (n_reps_S1 + n_reps_S2) * n_neighbors)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points)

        if (present(neighbor_mask_S1)) neighbor_mask_S1_f = neighbor_mask_S1
        if (present(neighbor_mask_S2)) neighbor_mask_S2_f = neighbor_mask_S2

        call gjct_permutation_test(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            global_jsd_observed = global_jsd_observed,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            n_permutations = n_permutations,&
            jsd_null = jsd_null,&
            p_value = p_value,&
            tmp_residuals_S1 = tmp_residuals_S1,&
            tmp_residuals_S2 = tmp_residuals_S2,&
            tmp_pool = tmp_pool,&
            tmp_pmf_S1 = tmp_pmf_S1,&
            tmp_pmf_S2 = tmp_pmf_S2,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps_S1 = tmp_included_n_reps_S1,&
            tmp_included_n_reps_S2 = tmp_included_n_reps_S2,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            random_seed = random_seed,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f,&
            ierr = ierr&
        )
    end subroutine gjct_permutation_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration_stats(module):gjct_permutation_test_alloc(subroutine)]]
    !| Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
    !| the work copies, so the caller's own arrays are left untouched.
    subroutine gjct_permutation_test_c(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            global_jsd_observed,&
            n_bins,&
            shared_residual_range,&
            n_permutations,&
            jsd_null,&
            p_value,&
            random_seed,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            ierr&
        ) bind(C, name="gjct_permutation_test_c")
        use tox_data_integration_stats, only: gjct_permutation_test_alloc

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations to perform
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(c_double), intent(in), target :: global_jsd_observed
            !! Observed global JSD value for both studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_permutations), intent(out), target :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(c_double), intent(out), target :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        integer(c_int), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(:, :), allocatable :: neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(global_jsd_observed)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(p_value)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(jsd_null, n_permutations)

        if (present(neighbor_mask_S1)) neighbor_mask_S1_f = neighbor_mask_S1
        if (present(neighbor_mask_S2)) neighbor_mask_S2_f = neighbor_mask_S2

        call gjct_permutation_test_alloc(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            global_jsd_observed = global_jsd_observed,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            n_permutations = n_permutations,&
            jsd_null = jsd_null,&
            p_value = p_value,&
            random_seed = random_seed,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f,&
            ierr = ierr&
        )
    end subroutine gjct_permutation_test_c

end module tox_data_integration_stats_c
#endif
