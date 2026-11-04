#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_stats(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bootstrap
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
module tox_data_integration_stats_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: bootstrap_histogram_c
    public :: bootstrap_histogram_expert_c

contains

    !> summary: C-wrapper for [[tox_data_integration_stats(module):bootstrap_histogram(subroutine)]]
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
    subroutine bootstrap_histogram_c(&
            n_bootstraps,&
            n_bins,&
            n_points,&
            n_studies,&
            counts,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            confidence_interval,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            divergence_mode,&
            ierr&
        ) bind(C, name="bootstrap_histogram_c")
        use tox_data_integration_stats, only: bootstrap_histogram
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

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
        integer(c_int), dimension(n_points, n_bins, n_studies), intent(in), target :: counts
            !! Every study's own residual counts per bin: what each study is drawn from; at each
            !! reference point they must add up to the study's `included_n_reps` there
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of all studies together (the
            !! consensus's count, the denominator of the weights); must be the sum of the studies'
            !! `included_n_reps` there
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
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence of every bootstrap draw is computed, as in the
            !! routine whose interval this is
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

        block
            character(len=:), pointer :: divergence_mode_f
            divergence_mode_f => c_char_as_view(divergence_mode)

            select case (divergence_mode_f)
                case ("divergence_symmetric")
                    divergence_mode_mode_f = MODE_DIVERGENCE_SYMMETRIC
                case ("divergence_pairwise_k2")
                    divergence_mode_mode_f = MODE_DIVERGENCE_PAIRWISE_K2
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call bootstrap_histogram(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            counts = counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            confidence_interval = confidence_interval,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            divergence_mode = divergence_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_c

    !> summary: C-wrapper for [[tox_data_integration_stats(module):bootstrap_histogram_expert(subroutine)]]
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
    subroutine bootstrap_histogram_expert_c(&
            n_bootstraps,&
            n_bins,&
            n_points,&
            n_studies,&
            counts,&
            mean_pmf_included_n_reps,&
            included_n_reps,&
            n_bootstrapping_top_k_jsds,&
            confidence_interval,&
            tmp_bootstrapping_top_k_jsds,&
            tmp_pmfs,&
            tmp_mean_pmf,&
            tmp_js_divergences,&
            tmp_weights,&
            tmp_global_js_divergence,&
            tmp_point_pool,&
            tmp_point_draw,&
            two_sided_bootstrapping_significance_level,&
            random_seed,&
            divergence_mode,&
            ierr&
        ) bind(C, name="bootstrap_histogram_expert_c")
        use tox_data_integration_stats, only: bootstrap_histogram_expert
        use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC

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
        integer(c_int), dimension(n_points, n_bins, n_studies), intent(in), target :: counts
            !! Every study's own residual counts per bin: what each study is drawn from; at each
            !! reference point they must add up to the study's `included_n_reps` there
            !! The minimum valid value is `0_int32`.
        integer(c_int), dimension(n_points), intent(in), target :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of all studies together (the
            !! consensus's count, the denominator of the weights); must be the sum of the studies'
            !! `included_n_reps` there
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
            !! Working array: one reference point's counts, the pool of one draw
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
        character(len=1, kind=c_char), dimension(22), intent(in), target :: divergence_mode
            !! How each study's per-point divergence of every bootstrap draw is computed, as in the
            !! routine whose interval this is
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `'divergence_pairwise_k2'`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)
        integer(int32) :: divergence_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(random_seed)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(mean_pmf_included_n_reps, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(confidence_interval, 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds * 2 * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_pmfs, n_points * n_bins * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_mean_pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points * n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_global_js_divergence, n_studies)
        M_CHECK_ARRAY_NON_NULL(tmp_point_pool, n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_point_draw, n_bins)
        M_CHECK_ARRAY_NON_NULL(divergence_mode, 22)

        block
            character(len=:), pointer :: divergence_mode_f
            divergence_mode_f => c_char_as_view(divergence_mode)

            select case (divergence_mode_f)
                case ("divergence_symmetric")
                    divergence_mode_mode_f = MODE_DIVERGENCE_SYMMETRIC
                case ("divergence_pairwise_k2")
                    divergence_mode_mode_f = MODE_DIVERGENCE_PAIRWISE_K2
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call bootstrap_histogram_expert(&
            n_bootstraps = n_bootstraps,&
            n_bins = n_bins,&
            n_points = n_points,&
            n_studies = n_studies,&
            counts = counts,&
            mean_pmf_included_n_reps = mean_pmf_included_n_reps,&
            included_n_reps = included_n_reps,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            confidence_interval = confidence_interval,&
            tmp_bootstrapping_top_k_jsds = tmp_bootstrapping_top_k_jsds,&
            tmp_pmfs = tmp_pmfs,&
            tmp_mean_pmf = tmp_mean_pmf,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            tmp_global_js_divergence = tmp_global_js_divergence,&
            tmp_point_pool = tmp_point_pool,&
            tmp_point_draw = tmp_point_draw,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            random_seed = random_seed,&
            divergence_mode = divergence_mode_mode_f,&
            ierr = ierr&
        )
    end subroutine bootstrap_histogram_expert_c

end module tox_data_integration_stats_c
#endif
