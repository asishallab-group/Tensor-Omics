#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_stats(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test and Bootstrap
!|
!| The permutation test of the JS-Comp-Test:
!| [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl(interface)]] estimates an
!| empirical p-value for each of K studies' weighted global JSD against the consensus pmf, by
!| resampling every reference point's pooled histogram counts into the studies without
!| replacement. `p_values` is the primary result; `p_values_observed_consensus` uses the older
!| null (a fixed observed consensus) and is kept for comparison.
!|
!| The bootstrap
!| [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]] gives each
!| study's weighted global JSD a confidence interval, by drawing each study's histograms with
!| replacement from the pooled consensus counts; `calc_js_comp_test_n_top_k_jsds` sizes its top-/bottom-k
!| heaps for a caller that allocates its own.
module tox_data_integration_stats_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: gjct_permutation_test_c
    public :: gjct_permutation_test_expert_c
    public :: bootstrap_histogram_c
    public :: bootstrap_histogram_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_stats(module):gjct_permutation_test(subroutine)]]
    !| In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
    !| drawn without replacement (multivariate hypergeometric) into the studies, each study taking
    !| as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
    !| consensus (built exactly as
    !| [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]] builds the
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
    subroutine gjct_permutation_test_c(&
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
        ) bind(C, name="gjct_permutation_test_c")
        use tox_data_integration_stats, only: gjct_permutation_test

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations; `0` skips the test (both p-value outputs are then NaN)
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_points, n_bins), intent(in), target :: mean_pmf
            !! Observed consensus pmf; used only for `p_values_observed_consensus`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_studies), intent(in), target :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
            !! whose global JSD reached the observed one; NaN when `n_permutations` is 0
        real(c_double), dimension(n_studies), intent(out), target :: p_values_observed_consensus
            !! Like `p_values`, from the same permutations, but against the fixed observed consensus
            !! `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_jsd_observed, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values_observed_consensus, n_studies)

        call gjct_permutation_test(&
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
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine gjct_permutation_test_c

    !> summary: C-wrapper for [[tox_data_integration_stats(module):gjct_permutation_test_expert(subroutine)]]
    !| In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
    !| drawn without replacement (multivariate hypergeometric) into the studies, each study taking
    !| as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
    !| consensus (built exactly as
    !| [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]] builds the
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
    subroutine gjct_permutation_test_expert_c(&
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
        ) bind(C, name="gjct_permutation_test_expert_c")
        use tox_data_integration_stats, only: gjct_permutation_test_expert

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations; `0` skips the test (both p-value outputs are then NaN)
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_points, n_bins), intent(in), target :: mean_pmf
            !! Observed consensus pmf; used only for `p_values_observed_consensus`
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_studies), intent(in), target :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
            !! whose global JSD reached the observed one; NaN when `n_permutations` is 0
        real(c_double), dimension(n_studies), intent(out), target :: p_values_observed_consensus
            !! Like `p_values`, from the same permutations, but against the fixed observed consensus
            !! `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_mean_pmf_counts
            !! Working array: the remaining pool, reset to `mean_pmf_counts` every permutation
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Working array for one study's resampled histogram counts, one permutation at a time
        real(c_double), dimension(n_points, n_bins, n_studies), intent(out), target :: tmp_pmfs
            !! Working array that holds one permutation's resampled pmfs, all studies
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array for one permutation's per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_weights
            !! Working array for one permutation's per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array for one permutation's global weighted JSD values
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_permuted_mean_pmf
            !! Working array for one permutation's equal-weight consensus pmf of the permuted studies
        integer(c_int), dimension(n_bins), intent(out), target :: tmp_point_pool
            !! Working array: one reference point's remaining pool, drawn from and then written back
        integer(c_int), dimension(n_bins), intent(out), target :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_jsd_observed, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values_observed_consensus, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_permuted_mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, n_bins)

        call gjct_permutation_test_expert(&
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
    end subroutine gjct_permutation_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration_stats(module):bootstrap_histogram(subroutine)]]
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
    !| [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]] produces
    !| them); `ERR_ALLOC_FAIL` when the random number generator cannot be created. A failed draw
    !| later on is recorded (first failure only) without stopping the bootstrap.
    subroutine bootstrap_histogram_c(&
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
        ) bind(C, name="bootstrap_histogram_c")
        use tox_data_integration_stats, only: bootstrap_histogram

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)

        call bootstrap_histogram(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            mean_pmf_counts = mean_pmf_counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            confidence_interval = confidence_interval,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_c

    !> summary: C-wrapper for [[tox_data_integration_stats(module):bootstrap_histogram_expert(subroutine)]]
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
    !| [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]] produces
    !| them); `ERR_ALLOC_FAIL` when the random number generator cannot be created. A failed draw
    !| later on is recorded (first failure only) without stopping the bootstrap.
    subroutine bootstrap_histogram_expert_c(&
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
        ) bind(C, name="bootstrap_histogram_expert_c")
        use tox_data_integration_stats, only: bootstrap_histogram_expert

        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_points, n_bins), intent(in), target :: mean_pmf_counts
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(2, n_studies), intent(inout), target :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out), target :: tmp_bootstrapping_top_k_jsds
            !! Working array: per study the bottom-k heap `(:, 1, :)` and the top-k heap `(:, 2, :)`
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Working array: one study's resampled histogram counts
        real(c_double), dimension(n_points, n_bins, n_studies), intent(out), target :: tmp_pmfs
            !! Working array: one bootstrap's resampled pmfs, all studies
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_mean_pmf
            !! Working array: one bootstrap's consensus pmf
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array: one bootstrap's per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_weights
            !! Working array: one bootstrap's per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array: one bootstrap's weighted global JSD values
        integer(c_int), dimension(n_bins), intent(out), target :: tmp_point_pool
            !! Working array: one reference point's pooled counts, the pool of one draw
        integer(c_int), dimension(n_bins), intent(out), target :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, n_bins)

        call bootstrap_histogram_expert(&
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
    end subroutine bootstrap_histogram_expert_c

end module tox_data_integration_stats_c
#endif
