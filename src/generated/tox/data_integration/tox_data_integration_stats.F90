#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test
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
    !| The p-value applies the `(1+count)/(n+1)` Laplace add-one correction --
    !| `p_values(i) = anint(count(i)+1)/(n_permutations+1)` -- so a study whose observed JSD is
    !| never reached by any permutation gets `p = 1/(n_permutations+1)`, never `p = 0.0` exactly.
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
            !! Number of permutations to perform
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf -- the pool each
            !! permutation resamples from without replacement, per reference point
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_bins, n_points), intent(in) :: mean_pmf
            !! The consensus pmf built from all studies' pmfs
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study -- how
            !! many elements are drawn (without replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study: the Laplace-corrected fraction of permutations whose
            !! resampled global JSD reached or exceeded the observed value -- see the correction
            !! note above
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator
        integer(int32), dimension(:, :), allocatable :: tmp_mean_pmf_counts
        integer(int32), dimension(:, :), allocatable :: tmp_counts
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        real(real64), dimension(:, :), allocatable :: tmp_pmf_point_major
        integer(int32), dimension(:, :), allocatable :: tmp_counts_point_major

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_permutations, ierr, arg_pos=1_int32, min=0_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_all_in_range_int(mean_pmf_counts, n_bins * n_points, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_real(mean_pmf, n_bins * n_points, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=7_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=8_int32, min=0_int32)
        call validate_all_in_range_real(global_jsd_observed, n_studies, ierr, arg_pos=9_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_mean_pmf_counts(n_bins, n_points))
        M_ALLOCATE(tmp_counts(n_bins, n_points))
        M_ALLOCATE(tmp_pmfs(n_bins, n_points, n_studies))
        M_ALLOCATE(tmp_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_weights(n_points, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_pmf_point_major(n_points, n_bins))
        M_ALLOCATE(tmp_counts_point_major(n_points, n_bins))

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
        call clear_err_arg_pos(ierr)
    end subroutine gjct_permutation_test

    !> summary: Validates its inputs, then calls [[tox_data_integration_stats_impl(module):gjct_permutation_test_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_stats(module):gjct_permutation_test]] does both.
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
    !| The p-value applies the `(1+count)/(n+1)` Laplace add-one correction --
    !| `p_values(i) = anint(count(i)+1)/(n_permutations+1)` -- so a study whose observed JSD is
    !| never reached by any permutation gets `p = 1/(n_permutations+1)`, never `p = 0.0` exactly.
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
            !! Number of permutations to perform
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf -- the pool each
            !! permutation resamples from without replacement, per reference point
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_bins, n_points), intent(in) :: mean_pmf
            !! The consensus pmf built from all studies' pmfs
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study -- how
            !! many elements are drawn (without replacement) from the pooled pool per study
            !! The minimum valid value is `0_int32`.
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study: the Laplace-corrected fraction of permutations whose
            !! resampled global JSD reached or exceeded the observed value -- see the correction
            !! note above
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_mean_pmf_counts
            !! Working array for proper resampling per permutation: the remaining pool to draw
            !! from, reset to `mean_pmf_counts` at the start of every permutation
        integer(int32), dimension(n_bins, n_points), intent(out) :: tmp_counts
            !! Working array for one study's resampled histogram counts, one permutation at a time
        real(real64), dimension(n_bins, n_points, n_studies), intent(out) :: tmp_pmfs
            !! Working array that holds one permutation's resampled pmfs, all studies
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Working array for one permutation's per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Working array for one permutation's per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array for one permutation's global weighted JSD values
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_point_major
            !! Working array: one study's resampled pmf, point-major, for calc_pmf_impl's
            !! point-major convention before transposing into tmp_pmfs
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts_point_major
            !! Working array: one study's resampled counts, point-major, transposed from tmp_counts
            !! before calc_pmf_impl
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_permutations, ierr, arg_pos=1_int32, min=0_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_all_in_range_int(mean_pmf_counts, n_bins * n_points, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_real(mean_pmf, n_bins * n_points, ierr, arg_pos=6_int32, min=0.0_real64, max=1.0_real64)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=7_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=8_int32, min=0_int32)
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
        call clear_err_arg_pos(ierr)
    end subroutine gjct_permutation_test_expert

end module tox_data_integration_stats
