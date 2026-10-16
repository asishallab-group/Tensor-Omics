#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test and Bootstrap
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
module tox_data_integration_stats_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_multiv_hypergeom, random_multinomial
    use f42_sort_impl, only: top_k_heap_push, bottom_k_heap_push
    use tox_data_integration_jsd_impl, only: compute_divergence_per_reference_point_impl, &
                                             compute_weighted_global_divergence_impl, calc_pmf_impl, &
                                             create_mean_pmf_only_impl
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private
    public :: gjct_permutation_test_impl, calc_js_comp_test_n_top_k_jsds, bootstrap_histogram_impl
contains

    !> summary: Estimate how likely each study's observed weighted global JSD is to occur by chance
    !| AUTHOR_LASZLO_LANG
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
    subroutine gjct_permutation_test_impl(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                          mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values, &
                                          p_values_observed_consensus, tmp_mean_pmf_counts, tmp_counts, tmp_pmfs, &
                                          tmp_js_divergences, tmp_weights, tmp_global_js_divergence, tmp_permuted_mean_pmf, &
                                          tmp_point_pool, tmp_point_draw, random_seed, ierr)
        integer(int32), intent(in) :: n_permutations
            !! Number of permutations; `0` skips the test (both p-value outputs are then NaN)
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
            !! Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
            !! DM_MIN(0_int32)
        real(real64), dimension(n_points, n_bins), intent(in) :: mean_pmf
            !! Observed consensus pmf; used only for `p_values_observed_consensus`
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates (included ones) per reference point for the consensus pmf
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! DM_MIN(1_int32)
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
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        integer(int32) :: i_permutation, i_point, i_study, draw_ierr
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

            ! Observed-consensus null first: the block below overwrites the tmp_* arrays.
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

            ! The observed consensus's arithmetic, so the two match bit for bit (a few-ulp excess
            ! over 1.0 included, deliberately not clamped).
            call create_mean_pmf_only_impl(tmp_pmfs, n_bins, n_points, n_studies, tmp_permuted_mean_pmf)

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

    !> M_EXPORT_C
    !| summary: Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
    !| AUTHOR_LASZLO_LANG
    !| Sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    !| the default 2.5 and `n_bootstraps=1000`.
    pure subroutine calc_js_comp_test_n_top_k_jsds(n_bootstraps, two_sided_bootstrapping_significance_level, n_top_k)
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples that will be performed
            !! DM_MIN(1_int32)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level, as a percentage reserved at each end of the bootstrap
            !! distribution (e.g. 2.5 for a 95% two-sided interval)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(out) :: n_top_k
            !! Number of elements to keep at each end of the bootstrap distribution, at least 1

        real(real64) :: sig_level

        M_DEFAULT_VAL(two_sided_bootstrapping_significance_level, sig_level, 2.5_real64)

        n_top_k = max(1_int32, floor(sig_level/100.0_real64*real(n_bootstraps, real64), kind=int32))
    end subroutine calc_js_comp_test_n_top_k_jsds

    !> summary: Bootstrap a confidence interval for each study's global JSD by resampling the pooled consensus histogram
    !| AUTHOR_LASZLO_LANG
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
    subroutine bootstrap_histogram_impl(n_bootstraps, n_bins, n_points, n_studies, mean_pmf_counts, &
                                        mean_pmf_included_n_reps, included_n_reps, n_bootstrapping_top_k_jsds, &
                                        confidence_interval, tmp_bootstrapping_top_k_jsds, tmp_counts, tmp_pmfs, &
                                        tmp_mean_pmf, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                        tmp_point_pool, tmp_point_draw, two_sided_bootstrapping_significance_level, &
                                        random_seed, ierr)
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! DM_MIN(1_int32)
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
            !! Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of the pool
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! DM_OUTPUT_FROM(n_top_k, calc_js_comp_test_n_top_k_jsds, tox_data_integration_stats_impl, AUTO)
            !! DM_MIN(1_int32)
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array: per study the bottom-k heap `(:, 1, :)` and the top-k heap `(:, 2, :)`
        integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
            !! Working array: one study's resampled histogram counts
        real(real64), dimension(n_points, n_bins, n_studies), intent(out) :: tmp_pmfs
            !! Working array: one bootstrap's resampled pmfs, all studies
        real(real64), dimension(n_points, n_bins), intent(out) :: tmp_mean_pmf
            !! Working array: one bootstrap's consensus pmf
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_js_divergences
            !! Working array: one bootstrap's per-point JSD values
        real(real64), dimension(n_points, n_studies), intent(out) :: tmp_weights
            !! Working array: one bootstrap's per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array: one bootstrap's weighted global JSD values
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_pool
            !! Working array: one reference point's pooled counts, the pool of one draw
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(2.5_real64)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(42_int32)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        integer(int32) :: i_bootstrap, i_study, i_point, draw_ierr
        type(rng_t) :: rng

        call set_ok(ierr)
        ! The pool's total is the multinomial total and the weights' denominator, so it must match
        ! the studies' draws and its own bin counts (also covers drawing from an empty pool). Not
        ! one argument's fault, so no position.
        do i_point = 1, n_points
            if (sum(included_n_reps(i_point, :)) /= mean_pmf_included_n_reps(i_point) .or. &
                mean_pmf_included_n_reps(i_point) /= sum(mean_pmf_counts(i_point, :))) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                exit
            end if
        end do
        if (.not. is_err(ierr)) rng = create_rng(ierr, random_seed)
        if (is_err(ierr)) then
            ! Leave the work arrays defined (zeroed) and confidence_interval untouched.
            tmp_bootstrapping_top_k_jsds = 0.0_real64
            tmp_counts = 0_int32
            tmp_pmfs = 0.0_real64
            tmp_mean_pmf = 0.0_real64
            tmp_js_divergences = 0.0_real64
            tmp_weights = 0.0_real64
            tmp_global_js_divergence = 0.0_real64
            tmp_point_pool = 0_int32
            tmp_point_draw = 0_int32
            return
        end if

        do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, confidence_interval)
            tmp_bootstrapping_top_k_jsds(:, 1, i_study) = confidence_interval(1, i_study)
            tmp_bootstrapping_top_k_jsds(:, 2, i_study) = confidence_interval(2, i_study)
        end do

        do i_bootstrap = 1, n_bootstraps
            ! Sequential: every draw advances the rng state. A point's bins are a strided row here,
            ! so each draw goes through the contiguous tmp_point_pool/tmp_point_draw.
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    tmp_point_pool = mean_pmf_counts(i_point, :)
                    call random_multinomial(rng, n_bins, tmp_point_pool, mean_pmf_included_n_reps(i_point), &
                                            included_n_reps(i_point, i_study), tmp_point_draw, draw_ierr)
                    if (is_err(draw_ierr)) call set_err_once(ierr, get_err_code(draw_ierr))
                    tmp_counts(i_point, :) = tmp_point_draw
                end do

                ! A study without residuals at a point is rejected by the wrappers (DM_MIN(1));
                ! calc_pmf_impl gives such a point an all-zero pmf instead of dividing by 0.
                call calc_pmf_impl(tmp_counts, included_n_reps(:, i_study), n_points, n_bins, tmp_pmfs(:, :, i_study))
            end do

            ! Drawn with replacement, so the consensus changes too.
            call create_mean_pmf_only_impl(tmp_pmfs, n_bins, n_points, n_studies, tmp_mean_pmf)

            do concurrent(i_study=1:n_studies) shared(tmp_pmfs, tmp_mean_pmf, n_points, n_bins, tmp_js_divergences, &
                                                      included_n_reps, mean_pmf_included_n_reps, tmp_global_js_divergence, &
                                                      tmp_weights, tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds)
                call compute_divergence_per_reference_point_impl(tmp_pmfs(:, :, i_study), tmp_mean_pmf, &
                                                                 n_points, n_bins, tmp_js_divergences(:, i_study))
                call compute_weighted_global_divergence_impl(tmp_js_divergences(:, i_study), n_points, &
                                                              included_n_reps(:, i_study), mean_pmf_included_n_reps, &
                                                              tmp_global_js_divergence(i_study), tmp_weights(:, i_study))

                ! Keep the highest and lowest values seen so far.
                call bottom_k_heap_push(tmp_bootstrapping_top_k_jsds(:, 1, i_study), n_bootstrapping_top_k_jsds, &
                                        tmp_global_js_divergence(i_study))
                call top_k_heap_push(tmp_bootstrapping_top_k_jsds(:, 2, i_study), n_bootstrapping_top_k_jsds, &
                                     tmp_global_js_divergence(i_study))
            end do
        end do

        ! The heap roots: [largest of the small values, smallest of the large values].
        do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, confidence_interval)
            confidence_interval(1, i_study) = tmp_bootstrapping_top_k_jsds(1, 1, i_study)
            confidence_interval(2, i_study) = tmp_bootstrapping_top_k_jsds(1, 2, i_study)
        end do

        call destroy_rng(rng)
    end subroutine bootstrap_histogram_impl

end module tox_data_integration_stats_impl
