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
module tox_data_integration_stats_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_multiv_hypergeom
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, calc_pmf_impl
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code
    M_IMPLICIT_NONE
contains

    !> summary: Estimate how likely each study's observed weighted global JSD is to occur by chance
    !| AUTHOR_LASZLO_LANG
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
    subroutine gjct_permutation_test_impl(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                          mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values, &
                                          tmp_mean_pmf_counts, tmp_counts, tmp_pmfs, tmp_js_divergences, tmp_weights, &
                                          tmp_global_js_divergence, tmp_pmf_point_major, tmp_counts_point_major, &
                                          random_seed, ierr)
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform
            !! DM_MIN(0_int32)
        integer(int32), intent(in) :: n_bins
            !! Number of equally sized histogram bins
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), dimension(n_bins, n_points), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf -- the pool each
            !! permutation resamples from without replacement, per reference point
            !! DM_MIN(0_int32)
        real(real64), dimension(n_bins, n_points), intent(in) :: mean_pmf
            !! The consensus pmf built from all studies' pmfs
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study -- how
            !! many elements are drawn (without replacement) from the pooled pool per study
            !! DM_MIN(0_int32)
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study: the fraction of permutations whose resampled global
            !! JSD reached or exceeded the observed value -- see the known-limitation note above
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
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator

        integer(int32) :: i_permutation, i_point, i_study, draw_ierr
        type(rng_t) :: rng

        call set_ok(ierr)
        rng = create_rng(ierr, random_seed)
        if (is_err(ierr)) then
            ! A genuine runtime failure -- leave every work array in a defined (zeroed) state
            ! rather than resample with an uninitialized generator.
            p_values = 0.0_real64
            tmp_mean_pmf_counts = 0_int32
            tmp_counts = 0_int32
            tmp_pmfs = 0.0_real64
            tmp_js_divergences = 0.0_real64
            tmp_weights = 0.0_real64
            tmp_global_js_divergence = 0.0_real64
            tmp_pmf_point_major = 0.0_real64
            tmp_counts_point_major = 0_int32
            return
        end if

        p_values = 0.0_real64
        do i_permutation = 1, n_permutations
            ! Resample the histogram: each reference point's pooled consensus counts become a
            ! pool to draw from without replacement, per study. tmp_mean_pmf_counts is refreshed
            ! every permutation, so previous studies' draws never leak into the next permutation.
            tmp_mean_pmf_counts = mean_pmf_counts
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    call random_multiv_hypergeom(rng, n_bins, tmp_mean_pmf_counts(:, i_point), &
                                                 sum(tmp_mean_pmf_counts(:, i_point)), included_n_reps(i_point, i_study), &
                                                 tmp_counts(:, i_point), draw_ierr)
                    if (is_err(draw_ierr)) call set_err_once(ierr, get_err_code(draw_ierr))
                end do
                tmp_counts_point_major = transpose(tmp_counts)
                call calc_pmf_impl(tmp_counts_point_major, included_n_reps(:, i_study), n_points, n_bins, tmp_pmf_point_major)
                tmp_pmfs(:, :, i_study) = transpose(tmp_pmf_point_major)
            end do

            do concurrent(i_study=1:n_studies) shared(tmp_pmfs, mean_pmf, n_points, n_bins, tmp_js_divergences, &
                                                      included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, &
                                                      tmp_weights, global_jsd_observed, p_values)
                call compute_divergence_per_reference_point_impl(transpose(tmp_pmfs(:, :, i_study)), transpose(mean_pmf), &
                                                                 n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_impl(tmp_js_divergences(:, i_study), n_points, &
                                                              included_n_reps(:, i_study), mean_pmf_included_n_reps, &
                                                              tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                if (tmp_global_js_divergence(i_study) >= global_jsd_observed(i_study)) then
                    p_values(i_study) = p_values(i_study) + 1.0_real64
                end if
            end do
        end do

        if (n_permutations /= 0) then
            do concurrent(i_study=1:n_studies) shared(p_values, n_permutations)
                p_values(i_study) = anint(p_values(i_study))/real(n_permutations, real64)
            end do
        end if

        call destroy_rng(rng)
    end subroutine gjct_permutation_test_impl

end module tox_data_integration_stats_impl
