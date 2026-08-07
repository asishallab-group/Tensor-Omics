#include <src/macros.h>

!> summary: Wrappers for [[tox_data_integration_stats_impl(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_data_integration_stats
    use tox_data_integration_stats_impl, only: gjct_permutation_test_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: gjct_permutation_test
    public :: gjct_permutation_test_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_stats(module):gjct_permutation_test_expert]] to prepare it yourself.
    !| Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
    !| the work copies, so the caller's own arrays are left untouched.
    subroutine gjct_permutation_test(&
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
        )
        integer(int32), intent(in) :: n_reps_S1
            !! Number of replicates in study 1
        integer(int32), intent(in) :: n_reps_S2
            !! Number of replicates in study 2
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors in the studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points in the studies
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins used for the studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        integer(int32), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:, :, :), allocatable :: tmp_residuals_S1
        real(real64), dimension(:, :, :), allocatable :: tmp_residuals_S2
        real(real64), dimension(:, :), allocatable :: tmp_pool
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S1
        real(real64), dimension(:, :), allocatable :: tmp_pmf_S2
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        integer(int32), dimension(:), allocatable :: tmp_included_n_reps_S1
        integer(int32), dimension(:), allocatable :: tmp_included_n_reps_S2
        real(real64), dimension(:), allocatable :: tmp_js_divergences
        real(real64), dimension(:), allocatable :: tmp_weights

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_in_range_real(global_jsd_observed, ierr, arg_pos=7_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_dimension_size(n_permutations, ierr, arg_pos=10_int32)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true.)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_residuals_S1(n_reps_S1, n_neighbors, n_points))
        M_ALLOCATE(tmp_residuals_S2(n_reps_S2, n_neighbors, n_points))
        M_ALLOCATE(tmp_pool(n_reps_S1 + n_reps_S2, n_neighbors))
        M_ALLOCATE(tmp_pmf_S1(n_points, n_bins))
        M_ALLOCATE(tmp_pmf_S2(n_points, n_bins))
        M_ALLOCATE(tmp_counts(n_points, n_bins))
        M_ALLOCATE(tmp_included_n_reps_S1(n_points))
        M_ALLOCATE(tmp_included_n_reps_S2(n_points))
        M_ALLOCATE(tmp_js_divergences(n_points))
        M_ALLOCATE(tmp_weights(n_points))

        call gjct_permutation_test_impl(&
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
            neighbor_mask_S1 = neighbor_mask_S1,&
            neighbor_mask_S2 = neighbor_mask_S2&
        )
    end subroutine gjct_permutation_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_stats(module):gjct_permutation_test]] does both.
    !| Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
    !| the work copies, so the caller's own arrays are left untouched.
    subroutine gjct_permutation_test_expert(&
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
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
            !! NaN is permitted for this value.
        real(real64), intent(in) :: global_jsd_observed
            !! Observed global JSD value for both studies
        real(real64), intent(in) :: shared_residual_range
            !! Computed residual range for both studies
            !! The minimum valid value is `0.0_real64`.
        real(real64), dimension(n_permutations), intent(out) :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(real64), intent(out) :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(out) :: tmp_residuals_S1
            !! Work copy of study 1's residuals, shuffled in place
        real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(out) :: tmp_residuals_S2
            !! Work copy of study 2's residuals, shuffled in place
        real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out) :: tmp_pool
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
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_reps_S1, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S2, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=6_int32)
        call validate_in_range_real(global_jsd_observed, ierr, arg_pos=7_int32)
        call validate_dimension_size(n_bins, ierr, arg_pos=8_int32)
        call validate_in_range_real(shared_residual_range, ierr, arg_pos=9_int32, min=0.0_real64)
        call validate_dimension_size(n_permutations, ierr, arg_pos=10_int32)
        call validate_all_in_range_real(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points, ierr, arg_pos=1_int32, allow_nan=.true.)
        call validate_all_in_range_real(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points, ierr, arg_pos=2_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        call gjct_permutation_test_impl(&
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
            neighbor_mask_S1 = neighbor_mask_S1,&
            neighbor_mask_S2 = neighbor_mask_S2&
        )
    end subroutine gjct_permutation_test_expert

end module tox_data_integration_stats
