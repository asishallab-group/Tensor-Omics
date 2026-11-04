#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bootstrap
!|
!| The bootstrap
!| [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]] gives each
!| study's weighted global JSD an interval, by drawing each study's histograms again with
!| replacement from the study's own counts: a robustness interval, how much the observed value
!| would move if each study's residuals were drawn again. `calc_js_comp_test_n_top_k_jsds` sizes
!| its top-/bottom-k heaps for a caller that allocates its own. The significance test of the
!| JS-Comp-Test is the study-label Permutation-Test,
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]]
!| and its adaptive counterpart.
module tox_data_integration_stats_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_multinomial
    use f42_sort_impl, only: top_k_heap_push, bottom_k_heap_push
    use tox_data_integration_jsd_impl, only: compute_consensus_divergences, create_mean_pmf_only_impl, pmf_entry
    use tox_errors, only: set_ok, set_err_once, is_err, validate_in_range_int, validate_in_range_real, &
                          ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_n_top_k_jsds, bootstrap_heap_size, bootstrap_histogram_impl
contains

    !> M_EXPORT_C
    !| summary: Recommend the bootstrap top-/bottom-k heap size for a two-sided confidence interval
    !| AUTHOR_LASZLO_LANG
    !| Sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    !| the default 2.5 and `n_bootstraps=1000`.
    !|
    !| Rejects an out-of-range argument with invalid input and a non-finite significance level with
    !| `ERR_NAN_INF`; `n_top_k` is then left undefined.
    pure subroutine calc_js_comp_test_n_top_k_jsds(n_bootstraps, two_sided_bootstrapping_significance_level, n_top_k, &
                                                   ierr)
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples that will be performed
            !! DM_MIN(1_int32)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level, as a percentage reserved at each end of the bootstrap
            !! distribution (e.g. 2.5 for a 95% two-sided interval)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(out) :: n_top_k
            !! Number of elements to keep at each end of the bootstrap distribution, at least 1
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        real(real64) :: sig_level

        call set_ok(ierr)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=2_int32, &
                                    min=0.0_real64, max=100.0_real64)
        if (is_err(ierr)) return

        M_DEFAULT_VAL(two_sided_bootstrapping_significance_level, sig_level, M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)

        n_top_k = bootstrap_heap_size(n_bootstraps, sig_level)
    end subroutine calc_js_comp_test_n_top_k_jsds

    !> summary: The bootstrap top-/bottom-k heap size, without validation
    !| AUTHOR_LASZLO_LANG
    !| `max(1, floor(level/100 * n_bootstraps))`, the formula of
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]], for a
    !| caller whose inputs are already validated (`n_bootstraps >= 1`, the level in `[0, 100]`).
    !| Internal helper, not published; public because the study-label Permutation-Test's sizing
    !| routines call it.
    pure integer(int32) function bootstrap_heap_size(n_bootstraps, two_sided_bootstrapping_significance_level) &
        result(n_top_k)
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples, at least 1
        real(real64), intent(in) :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level in percent per end, in `[0, 100]`

        n_top_k = max(1_int32, floor(two_sided_bootstrapping_significance_level/100.0_real64*real(n_bootstraps, real64), &
                                     kind=int32))
    end function bootstrap_heap_size

    !> summary: Whether every reference point's pool total is the sum of the studies' draw sizes
    !| AUTHOR_LASZLO_LANG
    !| Private. `.true.` iff at every point the studies' `included_n_reps` add up to
    !| `mean_pmf_included_n_reps`. The scan stops at the first point that fails; a point's studies
    !| are a strided row, so the sum is a plain loop rather than a reduction over the section.
    pure logical(c_bool) function pool_totals_consistent(n_points, n_studies, mean_pmf_included_n_reps, included_n_reps) &
        result(consistent)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: mean_pmf_included_n_reps(n_points)
            !! Pool total per reference point
        integer(int32), intent(in) :: included_n_reps(n_points, n_studies)
            !! Each study's draw size per reference point

        integer(int32) :: i_point, i_study, total

        consistent = .false._c_bool
        ! Sequential with an early exit: the first inconsistent point decides.
        do i_point = 1, n_points
            total = 0_int32
            do i_study = 1, n_studies
                total = total + included_n_reps(i_point, i_study)
            end do
            if (total /= mean_pmf_included_n_reps(i_point)) return
        end do
        consistent = .true._c_bool
    end function pool_totals_consistent

    !> summary: Whether every study's own counts add up to its draw size at every reference point
    !| AUTHOR_LASZLO_LANG
    !| Private. `.true.` iff `counts(i_point, :, i_study)` sums to `included_n_reps(i_point, i_study)`
    !| for every point and study. The scan stops at the first mismatch; the bins are a strided row,
    !| so each sum is a plain loop rather than a reduction over the section.
    pure logical(c_bool) function study_totals_consistent(n_points, n_bins, n_studies, counts, included_n_reps) &
        result(consistent)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_bins
            !! Number of histogram bins
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: counts(n_points, n_bins, n_studies)
            !! Every study's counts per bin
        integer(int32), intent(in) :: included_n_reps(n_points, n_studies)
            !! Each study's draw size per reference point

        integer(int32) :: i_point, i_study, i_bin, total

        consistent = .false._c_bool
        ! Sequential with an early exit: the first mismatch decides.
        do i_study = 1, n_studies
            do i_point = 1, n_points
                total = 0_int32
                do i_bin = 1, n_bins
                    total = total + counts(i_point, i_bin, i_study)
                end do
                if (total /= included_n_reps(i_point, i_study)) return
            end do
        end do
        consistent = .true._c_bool
    end function study_totals_consistent

    !> summary: Bootstrap an interval for each study's global JSD by resampling the studies' histograms
    !| AUTHOR_LASZLO_LANG
    !| In every bootstrap, each study's histogram at each reference point is drawn again with
    !| replacement (multinomial), with as many residuals as the study observed
    !| (`included_n_reps`); the consensus pmf, each study's per-point JSD and its weighted global
    !| JSD are recomputed from the draw, with the divergence `divergence_mode` selects (by default
    !| for exactly two studies their pairwise JSD, otherwise against the consensus). Each study is
    !| drawn from its own counts (`counts`): the interval shows how much the observed JSD would move
    !| if each study's own residuals were drawn again -- a robustness interval of the observed
    !| value. Drawing again adds sampling noise, so the bootstrapped JSDs are larger than the
    !| observed one on average: when the studies hardly differ, the interval can lie entirely above
    !| the observed JSD.
    !|
    !| `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
    !| slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
    !| values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
    !| smallest of the n_bootstrapping_top_k_jsds largest values]`.
    !|
    !| Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
    !| draw, when at some reference point the studies' `included_n_reps` do not add up to
    !| `mean_pmf_included_n_reps` (as
    !| [[tox_data_integration_jsd_impl(module):create_mean_pmf_impl(interface)]] produces
    !| them), or a study's `counts` at some reference point do not add up to its
    !| `included_n_reps`; `ERR_ALLOC_FAIL` when the random number generator cannot be
    !| created. The count checks, together with the wrappers' bounds on the counts, establish every
    !| draw's preconditions, so no draw can fail later on.
    subroutine bootstrap_histogram_impl(n_bootstraps, n_bins, n_points, n_studies, counts, &
                                        mean_pmf_included_n_reps, included_n_reps, n_bootstrapping_top_k_jsds, &
                                        confidence_interval, tmp_bootstrapping_top_k_jsds, tmp_pmfs, &
                                        tmp_mean_pmf, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                        tmp_point_pool, tmp_point_draw, two_sided_bootstrapping_significance_level, &
                                        random_seed, divergence_mode, ierr)
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
        integer(int32), dimension(n_points, n_bins, n_studies), intent(in) :: counts
            !! Every study's own residual counts per bin: what each study is drawn from; at each
            !! reference point they must add up to the study's `included_n_reps` there
            !! DM_MIN(0_int32)
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of all studies together (the
            !! consensus's count, the denominator of the weights); must be the sum of the studies'
            !! `included_n_reps` there
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
            !! Working array: one reference point's counts, the pool of one draw
        integer(int32), dimension(n_bins), intent(out) :: tmp_point_draw
            !! Working array: one study's draw at one reference point
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence of every bootstrap draw is computed, as in the
            !! routine whose interval this is
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Each study's JSD against the consensus of all studies | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]] |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! DM_DEFAULT(M_MODE_DIVERGENCE_PAIRWISE_K2)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        integer(int32) :: i_bootstrap, i_study, i_point, i_bin, actual_random_seed, actual_divergence_mode
        type(rng_t) :: rng

        call set_ok(ierr)
        M_DEFAULT_VAL(random_seed, actual_random_seed, M_RANDOM_SEED_DEFAULT)
        M_DEFAULT_VAL(divergence_mode, actual_divergence_mode, M_MODE_DIVERGENCE_PAIRWISE_K2)
        ! The pool's total is the weights' denominator, so it must be the sum of the studies' draws.
        ! Not one argument's fault, so no position.
        if (.not. pool_totals_consistent(n_points, n_studies, mean_pmf_included_n_reps, included_n_reps)) &
            call set_err_once(ierr, ERR_INVALID_INPUT)
        ! Per study, its own counts are the pool and its own total the multinomial total: they must
        ! agree in the same way (also covers drawing from an empty pool).
        if (.not. is_err(ierr)) then
            if (.not. study_totals_consistent(n_points, n_bins, n_studies, counts, included_n_reps)) &
                call set_err_once(ierr, ERR_INVALID_INPUT)
        end if
        if (.not. is_err(ierr)) rng = create_rng(ierr, actual_random_seed)
        if (is_err(ierr)) then
            ! Leave the work arrays defined (zeroed) and confidence_interval untouched.
            tmp_bootstrapping_top_k_jsds = 0.0_real64
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
            ! Sequential: every draw advances the rng state (study outer, point inner). A point's
            ! bins are a strided row here, so each draw goes through the contiguous
            ! tmp_point_pool/tmp_point_draw.
            do i_study = 1, n_studies
                do i_point = 1, n_points
                    tmp_point_pool = counts(i_point, :, i_study)
                    call random_multinomial(rng, n_bins, tmp_point_pool, included_n_reps(i_point, i_study), &
                                            included_n_reps(i_point, i_study), tmp_point_draw)

                    ! The draw's pmf, entry by entry with calc_pmf_impl's own formula (pmf_entry).
                    do i_bin = 1, n_bins
                        tmp_pmfs(i_point, i_bin, i_study) = pmf_entry(tmp_point_draw(i_bin), &
                                                                      included_n_reps(i_point, i_study))
                    end do
                end do
            end do

            ! Drawn with replacement, so the consensus changes too.
            ! (Built but not read in the pairwise mode with two studies, see compute_consensus_divergences.)
            call create_mean_pmf_only_impl(tmp_pmfs, n_bins, n_points, n_studies, tmp_mean_pmf)

            call compute_consensus_divergences(tmp_pmfs, tmp_mean_pmf, n_points, n_bins, n_bins, n_studies, &
                                               included_n_reps, mean_pmf_included_n_reps, actual_divergence_mode, &
                                               tmp_js_divergences, tmp_global_js_divergence, tmp_weights)
            do concurrent(i_study=1:n_studies) shared(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, &
                                                      tmp_global_js_divergence)
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
