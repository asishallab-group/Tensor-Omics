#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_stats(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test
!|
!| A permutation test estimating an empirical p-value for each study's weighted global JSD
!| against the consensus pmf. Under the null hypothesis that a study is exchangeable with the
!| pool, every reference point's pooled consensus histogram counts are repeatedly resampled
!| without replacement (GSL's multivariate hypergeometric distribution) into each study's own
!| draw size, the pmf/JSD/weighted global JSD are recomputed from that resample, and the observed
!| value is compared against the resulting null distribution.
!|
!| Generalized to K studies; the pooled resampling pool (`tmp_mean_pmf_counts`) is a work copy
!| reset every permutation, so the caller's own `mean_pmf_counts` is left untouched.
module tox_data_integration_stats_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: gjct_permutation_test_c
    public :: gjct_permutation_test_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_stats(module):gjct_permutation_test(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `gjct_permutation_test_helper`, generalized from its
    !| hardcoded 2-study shuffle to a K-study GSL `random_multiv_hypergeom` resample from the
    !| pooled `mean_pmf_counts` (without replacement, per reference point, per study). Impure:
    !| draws random numbers, so it carries `ierr` for the one genuine runtime failure a validated
    !| caller cannot foresee -- `create_rng` failing to allocate the GSL generator -- folded in as
    !| `ERR_ALLOC_FAIL`; a later `random_multiv_hypergeom` failure (which validated,
    !| internally-consistent inputs should never trigger) is likewise folded in, first failure
    !| only, without stopping the resampling already in flight, matching
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]'s own
    !| precedent.
    !|
    !| Known limitation: ported verbatim from 125-stabilize-jscomp's p-value formula, WITHOUT the
    !| `(1+count)/(n+1)` Laplace add-one correction -- `p_values(i) = anint(count(i))/n_permutations`,
    !| so a study whose observed JSD is never reached by any permutation gets `p = 0.0` exactly,
    !| not the `1/(n_permutations+1)` a Laplace-corrected formula would give. Deliberately ported
    !| as-is; see the project's JSD-Comp-Test follow-up issue for the fix.
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
            !! Number of permutations to perform
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf -- the pool each
            !! permutation resamples from without replacement, per reference point
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_bins, n_points), intent(in), target :: mean_pmf
            !! The consensus pmf built from all studies' pmfs
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study -- how
            !! many elements are drawn (without replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_studies), intent(in), target :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study: the fraction of permutations whose resampled global
            !! JSD reached or exceeded the observed value -- see the known-limitation note above
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_jsd_observed, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)

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
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine gjct_permutation_test_c

    !> summary: C-wrapper for [[tox_data_integration_stats(module):gjct_permutation_test_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `gjct_permutation_test_helper`, generalized from its
    !| hardcoded 2-study shuffle to a K-study GSL `random_multiv_hypergeom` resample from the
    !| pooled `mean_pmf_counts` (without replacement, per reference point, per study). Impure:
    !| draws random numbers, so it carries `ierr` for the one genuine runtime failure a validated
    !| caller cannot foresee -- `create_rng` failing to allocate the GSL generator -- folded in as
    !| `ERR_ALLOC_FAIL`; a later `random_multiv_hypergeom` failure (which validated,
    !| internally-consistent inputs should never trigger) is likewise folded in, first failure
    !| only, without stopping the resampling already in flight, matching
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]'s own
    !| precedent.
    !|
    !| Known limitation: ported verbatim from 125-stabilize-jscomp's p-value formula, WITHOUT the
    !| `(1+count)/(n+1)` Laplace add-one correction -- `p_values(i) = anint(count(i))/n_permutations`,
    !| so a study whose observed JSD is never reached by any permutation gets `p = 0.0` exactly,
    !| not the `1/(n_permutations+1)` a Laplace-corrected formula would give. Deliberately ported
    !| as-is; see the project's JSD-Comp-Test follow-up issue for the fix.
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
            tmp_mean_pmf_counts,&
            tmp_counts,&
            tmp_pmfs,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_pmf_point_major,&
            tmp_counts_point_major,&
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
            !! Number of permutations to perform
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_bins, n_points), intent(in), target :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf -- the pool each
            !! permutation resamples from without replacement, per reference point
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_bins, n_points), intent(in), target :: mean_pmf
            !! The consensus pmf built from all studies' pmfs
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points, n_studies), intent(in), target :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study -- how
            !! many elements are drawn (without replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(c_double), dimension(n_studies), intent(in), target :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(c_double), dimension(n_studies), intent(out), target :: p_values
            !! Empirical p-value per study: the fraction of permutations whose resampled global
            !! JSD reached or exceeded the observed value -- see the known-limitation note above
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: tmp_mean_pmf_counts
            !! Working array for proper resampling per permutation: the remaining pool to draw
            !! from, reset to `mean_pmf_counts` at the start of every permutation
        integer(c_int), dimension(n_bins, n_points), intent(out), target :: tmp_counts
            !! Working array for one study's resampled histogram counts, one permutation at a time
        real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: tmp_pmfs
            !! Working array that holds one permutation's resampled pmfs, all studies
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_js_divergences
            !! Working array for one permutation's per-point JSD values
        real(c_double), dimension(n_points, n_studies), intent(out), target :: tmp_weights
            !! Working array for one permutation's per-point weights
        real(c_double), dimension(n_studies), intent(out), target :: tmp_global_js_divergence
            !! Working array for one permutation's global weighted JSD values
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_point_major
            !! Working array: one study's resampled pmf, point-major, for calc_pmf_impl's
            !! point-major convention before transposing into tmp_pmfs
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts_point_major
            !! Working array: one study's resampled counts, point-major, transposed from tmp_counts
            !! before calc_pmf_impl
        integer(c_int), intent(in), target :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(global_jsd_observed, n_studies)
        M_CHECK_ARRAY_NON_NULL(p_values, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_bins * n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_bins * n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_point_major, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts_point_major, n_points * n_bins)

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
            tmp_mean_pmf_counts = tmp_mean_pmf_counts,&
            tmp_counts = tmp_counts,&
            tmp_pmfs = tmp_pmfs,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_pmf_point_major = tmp_pmf_point_major,&
            tmp_counts_point_major = tmp_counts_point_major,&
            random_seed = random_seed,&
            ierr = ierr&
        )
    end subroutine gjct_permutation_test_expert_c

end module tox_data_integration_stats_c
#endif
