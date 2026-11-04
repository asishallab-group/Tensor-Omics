#include <src/macros.h>

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
!|
!| Generated from [[tox_data_integration_stats_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_stats
    use tox_data_integration_stats_impl, only: bootstrap_histogram_impl, calc_js_comp_test_n_top_k_jsds
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_PAIRWISE_K2, MODE_DIVERGENCE_SYMMETRIC
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    use tox_errors, only: clear_err_arg_pos, set_err, set_err_once, validate_all_in_range_int
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: bootstrap_histogram
    public :: bootstrap_histogram_expert

contains

    !> summary: Validates its inputs, prepares what [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_stats(module):bootstrap_histogram_expert]] to prepare it yourself.
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
    subroutine bootstrap_histogram(&
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
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_points, n_bins, n_studies), intent(in) :: counts
            !! Every study's own residual counts per bin: what each study is drawn from; at each
            !! reference point they must add up to the study's `included_n_reps` there
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of all studies together (the
            !! consensus's count, the denominator of the weights); must be the sum of the studies'
            !! `included_n_reps` there
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence of every bootstrap draw is computed, as in the
            !! routine whose interval this is
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `1_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)
        integer(int32) :: n_bootstrapping_top_k_jsds
        real(real64), dimension(:, :, :), allocatable :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(:, :, :), allocatable :: tmp_pmfs
        real(real64), dimension(:, :), allocatable :: tmp_mean_pmf
        real(real64), dimension(:, :), allocatable :: tmp_js_divergences
        real(real64), dimension(:, :), allocatable :: tmp_weights
        real(real64), dimension(:), allocatable :: tmp_global_js_divergence
        integer(int32), dimension(:), allocatable :: tmp_point_pool
        integer(int32), dimension(:), allocatable :: tmp_point_draw

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=9_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(counts, n_points * n_bins * n_studies, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=6_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=7_int32, min=1_int32)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=8_int32, min=0.0_real64, max=1.0_real64)
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=11_int32); end if
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_n_top_k_jsds(&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            n_top_k = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        if (is_err(ierr)) return
        M_ALLOCATE(tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        M_ALLOCATE(tmp_pmfs(n_points, n_bins, n_studies))
        M_ALLOCATE(tmp_mean_pmf(n_points, n_bins))
        M_ALLOCATE(tmp_js_divergences(n_points, n_studies))
        M_ALLOCATE(tmp_weights(n_points, n_studies))
        M_ALLOCATE(tmp_global_js_divergence(n_studies))
        M_ALLOCATE(tmp_point_pool(n_bins))
        M_ALLOCATE(tmp_point_draw(n_bins))

        call bootstrap_histogram_impl(&
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
            divergence_mode = divergence_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine bootstrap_histogram

    !> summary: Validates its inputs, then calls [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_stats(module):bootstrap_histogram]] does both.
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
    subroutine bootstrap_histogram_expert(&
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
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Number of values kept at each end of the bootstrap distribution (heap size)
            !! It is *VERY IMPORTANT* to compute this argument from the `n_top_k` output produced by [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds]].
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_bootstraps
            !! Number of bootstrap resamples to perform
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_points, n_bins, n_studies), intent(in) :: counts
            !! Every study's own residual counts per bin: what each study is drawn from; at each
            !! reference point they must add up to the study's `included_n_reps` there
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points), intent(in) :: mean_pmf_included_n_reps
            !! Count of non-NaN replicates per reference point of all studies together (the
            !! consensus's count, the denominator of the weights); must be the sum of the studies'
            !! `included_n_reps` there
            !! The minimum valid value is `0_int32`.
        integer(int32), dimension(n_points, n_studies), intent(in) :: included_n_reps
            !! Count of non-NaN replicates per reference point, per study: the size of each study's
            !! draw. Every study needs at least one at every reference point
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(2, n_studies), intent(inout) :: confidence_interval
            !! On entry the seed values of the heaps (see above); on return the bootstrapped
            !! `[lower, upper]` interval per study
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
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
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(in), optional :: random_seed
            !! Seed for the GSL random number generator
            !! The default value is `42_int32`.
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence of every bootstrap draw is computed, as in the
            !! routine whose interval this is
            !!
            !! | Mode                                                                                                             | Value                                                                           |
            !! |------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------|
            !! | Each study's JSD against the consensus of all studies                                                            | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]]   |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! The default value is `1_int32`.
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT for inconsistent counts, ERR_ALLOC_FAIL if the random
            !! number generator could not be created (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_bins, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_int(n_bootstrapping_top_k_jsds, ierr, arg_pos=8_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=18_int32, min=0.0_real64, max=100.0_real64)
        call validate_all_in_range_int(counts, n_points * n_bins * n_studies, ierr, arg_pos=5_int32, min=0_int32)
        call validate_all_in_range_int(mean_pmf_included_n_reps, n_points, ierr, arg_pos=6_int32, min=0_int32)
        call validate_all_in_range_int(included_n_reps, n_points * n_studies, ierr, arg_pos=7_int32, min=1_int32)
        call validate_all_in_range_real(confidence_interval, 2 * n_studies, ierr, arg_pos=9_int32, min=0.0_real64, max=1.0_real64)
        if (present(divergence_mode)) then; if (divergence_mode /= MODE_DIVERGENCE_SYMMETRIC .and. divergence_mode /= MODE_DIVERGENCE_PAIRWISE_K2) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=20_int32); end if
        if (is_err(ierr)) return
#endif

        call bootstrap_histogram_impl(&
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
            divergence_mode = divergence_mode,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine bootstrap_histogram_expert

end module tox_data_integration_stats
