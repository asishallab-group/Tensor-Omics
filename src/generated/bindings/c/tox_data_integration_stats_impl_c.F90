#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_stats_impl(module)]]
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
module tox_data_integration_stats_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_n_top_k_jsds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]]
    !| Sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_stats_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    !| the default 2.5 and `n_bootstraps=1000`.
    subroutine calc_js_comp_test_n_top_k_jsds_c(&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            n_top_k,&
            ierr&
        ) bind(C, name="calc_js_comp_test_n_top_k_jsds_c")
        use tox_data_integration_stats_impl, only: calc_js_comp_test_n_top_k_jsds

        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstrap resamples that will be performed
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level, as a percentage reserved at each end of the bootstrap
            !! distribution (e.g. 2.5 for a 95% two-sided interval)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(out), target :: n_top_k
            !! Number of elements to keep at each end of the bootstrap distribution, at least 1
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(n_top_k)

        call calc_js_comp_test_n_top_k_jsds(&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            n_top_k = n_top_k&
        )
    end subroutine calc_js_comp_test_n_top_k_jsds_c

end module tox_data_integration_stats_impl_c
#endif
