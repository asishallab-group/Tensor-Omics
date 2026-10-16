#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test
!|
!| A permutation test estimating an empirical p-value for each study's weighted global JSD
!| against the consensus pmf. Under the null hypothesis that a study is exchangeable with the
!| pool, every reference point's pooled consensus histogram counts are repeatedly resampled
!| without replacement (GSL's multivariate hypergeometric distribution) into each study's own
!| draw size, the pmfs, their equal-weight consensus and each study's JSD/weighted global JSD
!| against it are recomputed from that resample, and the observed value is compared against the
!| resulting null distribution (`p_values`). From the same resamples, a second p-value
!| (`p_values_observed_consensus`) compares each resampled study against the fixed observed
!| consensus pmf instead, the null this test used before; it is kept for comparison, and
!| `p_values` is the primary result.
!|
!| Generalized to K studies; the pooled resampling pool (`tmp_mean_pmf_counts`) is a work copy
!| reset every permutation, so the caller's own `mean_pmf_counts` is left untouched.
module tox_data_integration_stats_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_multiv_hypergeom
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, calc_pmf_impl
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private
    public :: gjct_permutation_test_impl
contains

    !> summary: Estimate how likely each study's observed weighted global JSD is to occur by chance
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `gjct_permutation_test_helper`, generalized from its
    !| hardcoded 2-study shuffle to a K-study GSL `random_multiv_hypergeom` resample from the
    !| pooled `mean_pmf_counts` (without replacement, per reference point, per study). Impure:
    !| draws random numbers, so it carries `ierr` for the genuine runtime failure a validated
    !| caller cannot foresee -- `create_rng` failing to allocate the GSL generator -- folded in as
    !| `ERR_ALLOC_FAIL` -- and for inconsistent inputs (see below); a later `random_multiv_hypergeom` failure (which validated,
    !| internally-consistent inputs should never trigger) is likewise folded in, first failure
    !| only, without stopping the resampling already in flight, matching
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]]'s own
    !| precedent.
    !|
    !| Inconsistent inputs are rejected before any draw. At every reference point, the studies'
    !| counts (`included_n_reps`) must add up to the pool's total (`mean_pmf_included_n_reps`), and
    !| that total must equal the sum of the pool's bin counts (`mean_pmf_counts`), as
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] produces
    !| them: a permutation splits the pool exactly between the studies. Otherwise the call fails
    !| with `ERR_INVALID_INPUT`, not tied to any one argument, and no permutation runs. No p-values
    !| are computed then (a Fortran caller gets them as NaN).
    !|
    !| The consensus pmf is recomputed from the permuted pmfs in every permutation, as the
    !| equal-weight mean of all studies' pmfs, exactly as
    !| [[tox_data_integration_js_comp_test_impl(module):create_mean_pmf_impl(interface)]] builds
    !| the observed one. So each permuted study's JSD is taken against the consensus of the
    !| permuted studies, and the null statistic is computed exactly like the observed one. This is
    !| required because the equal-weight mixture is not invariant under the permutation when the
    !| studies differ in size: the permuted studies follow the pooled, replicate-weighted
    !| distribution, whose equal-weight mean moves with every draw, so a fixed observed consensus
    !| would shift the whole null distribution.
    !|
    !| A second p-value, `p_values_observed_consensus`, is computed from the same permutations (no
    !| further random draws) with the null this test used before: each permuted study's JSD
    !| against the fixed observed consensus `mean_pmf`. It is kept to compare the two nulls on data
    !| with known batch effects; `p_values` stays the primary result.
    !|
    !| Both p-values apply the `(1+count)/(n+1)` Laplace add-one correction --
    !| `p_values(i) = anint(count(i)+1)/(n_permutations+1)` -- so a study whose observed JSD is
    !| never reached by any permutation gets `p = 1/(n_permutations+1)`, never `p = 0.0` exactly.
    !| With `n_permutations = 0` there is nothing to test: no permutation runs, and both p-values
    !| are NaN ("not computed"), not `0.0`.
    subroutine gjct_permutation_test_impl(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                          mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values, &
                                          p_values_observed_consensus, tmp_mean_pmf_counts, tmp_counts, tmp_pmfs, &
                                          tmp_js_divergences, tmp_weights, tmp_global_js_divergence, tmp_permuted_mean_pmf, &
                                          tmp_point_pool, tmp_point_draw, random_seed, ierr)
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations to perform. `0` skips the test: no permutation runs, and
            !! `p_values` and `p_values_observed_consensus` are then NaN ("not computed")
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
        integer(int32), dimension(n_points, n_bins), intent(in) :: mean_pmf_counts
            !! Absolute counts of a residual per bin for the consensus pmf -- the pool each
            !! permutation resamples from without replacement, per reference point
            !! DM_MIN(0_int32)
        real(real64), dimension(n_points, n_bins), intent(in) :: mean_pmf
            !! The consensus pmf built from all studies' pmfs; each permuted study is compared
            !! against it for `p_values_observed_consensus` only
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point, per study -- how
            !! many elements are drawn (without replacement) from the pooled pool per study. Every
            !! study needs at least one at every reference point
            !! DM_MIN(1_int32)
        real(real64), dimension(n_studies), intent(in) :: global_jsd_observed
            !! Observed weighted global JSD of each study against the consensus pmf
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study: the Laplace-corrected fraction of permutations whose
            !! resampled global JSD reached or exceeded the observed value -- see the correction
            !! note above. NaN when `n_permutations` is 0: no permutation was run, so no p-value
            !! was computed
        real(real64), dimension(n_studies), intent(out) :: p_values_observed_consensus
            !! Empirical p-value per study from the same permutations as `p_values`, but with each
            !! permuted study's JSD taken against the fixed observed consensus `mean_pmf` instead of
            !! the consensus of the permuted studies -- the null this test used before it recomputed
            !! the consensus per permutation. Kept for comparison only; `p_values` is the primary
            !! result. The same Laplace correction applies. NaN when `n_permutations` is 0: no
            !! permutation was run, so no p-value was computed
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_mean_pmf_counts
            !! Working array for proper resampling per permutation: the remaining pool to draw
            !! from, reset to `mean_pmf_counts` at the start of every permutation
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
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_ALLOC_FAIL if GSL could not allocate the random number generator;
            !! ERR_INVALID_INPUT, before any draw, if at some point the studies' counts do not add
            !! up to the pool's total, or that total is not the sum of the pool's bin counts (see above)

        integer(int32) :: i_permutation, i_point, i_study, i_bin, draw_ierr
        type(rng_t) :: rng

        call set_ok(ierr)
        ! Inconsistent inputs, checked before the generator is created: at every point the studies'
        ! draws must add up to the pool's total, and that total must be the sum of its bin counts,
        ! so that a permutation splits the pool exactly. Not one argument's fault, so no position
        ! (codegen_guide.md Sec 5.14).
        do i_point = 1, n_points
            if (sum(included_n_reps(i_point, :)) /= mean_pmf_included_n_reps(i_point) .or. &
                mean_pmf_included_n_reps(i_point) /= sum(mean_pmf_counts(i_point, :))) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                exit
            end if
        end do
        if (.not. is_err(ierr)) rng = create_rng(ierr, random_seed)
        if (is_err(ierr)) then
            ! Inconsistent inputs or a genuine runtime failure -- nothing was tested, so the
            ! p-values are NaN ("not computed"), not 0.0 (which would read as maximally
            ! significant); every work array is left in a defined (zeroed) state.
            p_values = M_NAN
            p_values_observed_consensus = M_NAN
            tmp_mean_pmf_counts = 0_int32
            tmp_counts = 0_int32
            tmp_pmfs = 0.0_real64
            tmp_js_divergences = 0.0_real64
            tmp_weights = 0.0_real64
            tmp_global_js_divergence = 0.0_real64
            tmp_permuted_mean_pmf = 0.0_real64
            tmp_point_pool = 0_int32
            tmp_point_draw = 0_int32
            return
        end if

        p_values = 0.0_real64
        p_values_observed_consensus = 0.0_real64
        do i_permutation = 1, n_permutations
            ! Resample the histogram: each reference point's pooled consensus counts become a
            ! pool to draw from without replacement, per study. tmp_mean_pmf_counts is refreshed
            ! every permutation, so previous studies' draws never leak into the next permutation.
            ! A point's bins are a strided row here, so each draw goes through the contiguous
            ! tmp_point_pool/tmp_point_draw, and the reduced pool is written back.
            tmp_mean_pmf_counts = mean_pmf_counts
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    tmp_point_pool = tmp_mean_pmf_counts(i_point, :)
                    call random_multiv_hypergeom(rng, n_bins, tmp_point_pool, sum(tmp_point_pool), &
                                                 included_n_reps(i_point, i_study), tmp_point_draw, draw_ierr)
                    if (is_err(draw_ierr)) call set_err_once(ierr, get_err_code(draw_ierr))
                    tmp_mean_pmf_counts(i_point, :) = tmp_point_pool
                    tmp_counts(i_point, :) = tmp_point_draw
                end do
                call calc_pmf_impl(tmp_counts, included_n_reps(:, i_study), n_points, n_bins, tmp_pmfs(:, :, i_study))
            end do

            ! The observed-consensus null, as before the consensus was recomputed per permutation:
            ! each permuted study against the fixed observed consensus mean_pmf, from the same draws.
            ! The block below overwrites the tmp_* arrays with the recomputed-consensus values.
            do concurrent(i_study=1:n_studies) shared(tmp_pmfs, mean_pmf, n_points, n_bins, tmp_js_divergences, &
                                                      included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, &
                                                      tmp_weights, global_jsd_observed, p_values_observed_consensus)
                call compute_divergence_per_reference_point_impl(tmp_pmfs(:, :, i_study), mean_pmf, &
                                                                 n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_impl(tmp_js_divergences(:, i_study), n_points, &
                                                              included_n_reps(:, i_study), mean_pmf_included_n_reps, &
                                                              tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                if (tmp_global_js_divergence(i_study) >= global_jsd_observed(i_study)) then
                    p_values_observed_consensus(i_study) = p_values_observed_consensus(i_study) + 1.0_real64
                end if
            end do

            ! The permuted studies' equal-weight consensus: the same formula and accumulation order
            ! as create_mean_pmf_impl (start at 0, then add pmf/n_studies study by study), so it
            ! matches that routine bit for bit. Inlined rather than called:
            ! tox_data_integration_js_comp_test_impl already uses this module.
            ! For some study counts (e.g. 9, 11, 18, 20) this sum exceeds 1.0 by a few ulps (one or
            ! two up to 41 studies) where every study has all its mass in the same bin. Harmless here,
            ! as it is never returned, but such a consensus, if fed to the public gjct_permutation_test,
            ! would trip its DM_MAX(1.0) on mean_pmf. Deliberately not clamped (bit-identical outputs).
            tmp_permuted_mean_pmf = 0.0_real64
            do i_study = 1, n_studies
                do i_bin = 1, n_bins
                    do concurrent(i_point=1:n_points) shared(tmp_permuted_mean_pmf, tmp_pmfs, i_bin, i_study, n_studies)
                        tmp_permuted_mean_pmf(i_point, i_bin) = tmp_permuted_mean_pmf(i_point, i_bin) &
                                                                + tmp_pmfs(i_point, i_bin, i_study)/real(n_studies, real64)
                    end do
                end do
            end do

            do concurrent(i_study=1:n_studies) shared(tmp_pmfs, tmp_permuted_mean_pmf, n_points, n_bins, tmp_js_divergences, &
                                                      included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, &
                                                      tmp_weights, global_jsd_observed, p_values)
                call compute_divergence_per_reference_point_impl(tmp_pmfs(:, :, i_study), tmp_permuted_mean_pmf, &
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
            do concurrent(i_study=1:n_studies) shared(p_values, p_values_observed_consensus, n_permutations)
                p_values(i_study) = anint(p_values(i_study)+1.0_real64)/real(n_permutations+1.0_real64, real64)
                p_values_observed_consensus(i_study) = anint(p_values_observed_consensus(i_study)+1.0_real64) &
                                                       /real(n_permutations+1.0_real64, real64)
            end do
        else
            ! No permutation ran, so nothing was tested: "not computed", not 0.0 (which would read as
            ! maximally significant), as in the error branch above.
            p_values = M_NAN
            p_values_observed_consensus = M_NAN
        end if

        call destroy_rng(rng)
    end subroutine gjct_permutation_test_impl

end module tox_data_integration_stats_impl
