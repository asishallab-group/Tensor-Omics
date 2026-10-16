#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test
!|
!| The permutation test of the JS-Comp-Test:
!| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] estimates an
!| empirical p-value for each of K studies' weighted global JSD against the consensus pmf, by
!| resampling every reference point's pooled histogram counts into the studies without
!| replacement. `p_values` is the primary result; `p_values_observed_consensus` uses the older
!| null (a fixed observed consensus) and is kept for comparison.
!|
!| Generated from [[tox_data_integration_stats_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_stats
    use tox_data_integration_stats_impl, only: gjct_permutation_test_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, clear_err_arg_pos
    use tox_errors, only: set_err, validate_all_in_range_int, validate_all_in_range_real, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: gjct_permutation_test
    public :: gjct_permutation_test_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_stats(module):gjct_permutation_test_expert]] to prepare it yourself.
    !| In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
    !| drawn without replacement (multivariate hypergeometric) into the studies, each study taking
    !| as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
    !| consensus (built exactly as
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] builds the
    !| observed one) and each study's weighted global JSD against that consensus are recomputed, and
    !| compared with the observed `global_jsd_observed`. `p_values_observed_consensus` comes from the
    !| same draws, taking each permuted study's JSD against the fixed observed `mean_pmf` instead.
    !|
    !| Both p-values are `(count + 1)/(n_permutations + 1)`, so never exactly 0. With
    !| `n_permutations = 0` no permutation runs and both are NaN ("not computed").
    !|
    !| Errors (`ierr`): `ERR_INVALID_INPUT`, before any draw, when at some reference point the
    !| studies' counts do not add up to the pool's total (`mean_pmf_included_n_reps`) or that total
    !| is not the sum of the pool's bin counts; `ERR_ALLOC_FAIL` when the random number generator
    !| cannot be created. The p-values are then NaN. A failed draw later on is recorded (first
    !| failure only) without stopping the remaining permutations.
    subroutine gjct_permutation_test(&
            n_permutations,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf_counts,&
            mean_pmf,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            global_jsd_observed,&
            p_values,&
            p_values_observed_consensus,&
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
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations; `0` skips the test (both p-value outputs are then NaN)
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_points, n_bins), intent(in) :: mean_pmf
            !! Observed consensus pmf; used only for `p_values_observed_consensus`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
            !! whose global JSD reached the observed one; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Like `p_values`, from the same permutations, but against the fixed observed consensus
            !! `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_permuted_mean_pmf
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_permutations, ierr, arg_pos=1_int32, min=0_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_all_in_range_int(mean_pmf_counts, n_points * n_bins, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_real(mean_pmf, n_points * n_bins, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=7_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=8_int32, min=1_int32)
        call validate_all_in_range_real(global_jsd_observed, n_studies, ierr, arg_pos=9_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_mean_pmf_counts(n_points, n_bins))
        M_ALLOCATE(tmp_counts(n_points, n_bins))
        M_ALLOCATE(tmp_pmfs(n_points, n_bins, n_studies))
        M_ALLOCATE(tmp_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_weights(n_points, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_permuted_mean_pmf(n_points, n_bins))
        M_ALLOCATE(tmp_point_pool(n_bins))
        M_ALLOCATE(tmp_point_draw(n_bins))

        call gjct_permutation_test_impl(&
            n_permutations = n_permutations,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf = mean_pmf,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            global_jsd_observed = global_jsd_observed,&
            p_values = p_values,&
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_counts = tmp_counts,&
            tmp_pmfs = tmp_pmfs,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_permuted_mean_pmf = tmp_permuted_mean_pmf,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine gjct_permutation_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_stats(module):gjct_permutation_test]] does both.
    !| In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
    !| drawn without replacement (multivariate hypergeometric) into the studies, each study taking
    !| as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
    !| consensus (built exactly as
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] builds the
    !| observed one) and each study's weighted global JSD against that consensus are recomputed, and
    !| compared with the observed `global_jsd_observed`. `p_values_observed_consensus` comes from the
    !| same draws, taking each permuted study's JSD against the fixed observed `mean_pmf` instead.
    !|
    !| Both p-values are `(count + 1)/(n_permutations + 1)`, so never exactly 0. With
    !| `n_permutations = 0` no permutation runs and both are NaN ("not computed").
    !|
    !| Errors (`ierr`): `ERR_INVALID_INPUT`, before any draw, when at some reference point the
    !| studies' counts do not add up to the pool's total (`mean_pmf_included_n_reps`) or that total
    !| is not the sum of the pool's bin counts; `ERR_ALLOC_FAIL` when the random number generator
    !| cannot be created. The p-values are then NaN. A failed draw later on is recorded (first
    !| failure only) without stopping the remaining permutations.
    subroutine gjct_permutation_test_expert(&
            n_permutations,&
            n_bins,&
            n_points,&
            n_studies,&
            mean_pmf_counts,&
            mean_pmf,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            global_jsd_observed,&
            p_values,&
            p_values_observed_consensus,&
            tmp_mean_pmf_counts,&
            tmp_counts,&
            tmp_pmfs,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_permuted_mean_pmf,&
            tmp_point_pool,&
            tmp_point_draw,&
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
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations; `0` skips the test (both p-value outputs are then NaN)
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_points, n_bins), intent(in) :: mean_pmf
            !! Observed consensus pmf; used only for `p_values_observed_consensus`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
            !! whose global JSD reached the observed one; NaN when `n_permutations` is 0
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Like `p_values`, from the same permutations, but against the fixed observed consensus
            !! `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_mean_pmf_counts
            !! Working array: the remaining pool, reset to `mean_pmf_counts` every permutation
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array for one study's resampled histogram counts, one permutation at a time
        real(real64), dimension(n_points, n_bins, n_studies), intent(out) :: tmp_pmfs
            !! Working array that holds one permutation's resampled pmfs, all studies
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Working array for one permutation's per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Working array for one permutation's per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array for one permutation's global weighted JSD values
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_permuted_mean_pmf
            !! Working array for one permutation's equal-weight consensus pmf of the permuted studies
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_pool
            !! Working array: one reference point's remaining pool, drawn from and then written back
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_permutations, ierr, arg_pos=1_int32, min=0_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_all_in_range_int(mean_pmf_counts, n_points * n_bins, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_real(mean_pmf, n_points * n_bins, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=7_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=8_int32, min=1_int32)
        call validate_all_in_range_real(global_jsd_observed, n_studies, ierr, arg_pos=9_int32)
        if (is_err(ierr)) return
#endif

        call gjct_permutation_test_impl(&
            n_permutations = n_permutations,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf = mean_pmf,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            global_jsd_observed = global_jsd_observed,&
            p_values = p_values,&
            p_values_observed_consensus = p_values_observed_consensus,&
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_counts = tmp_counts,&
            tmp_pmfs = tmp_pmfs,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_permuted_mean_pmf = tmp_permuted_mean_pmf,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            random_seed = random_seed,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine gjct_permutation_test_expert

end module tox_data_integration_stats
