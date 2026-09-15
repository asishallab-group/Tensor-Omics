#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
!| before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
!| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]],
!| each candidate's histogram bin count from
!| [[tox_data_integration_js_comp_test_impl(module):estimate_bin_count_impl(interface)]]), two
!| admissibility gates a candidate must pass before it is bootstrapped
!| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]),
!| and the plateau check that decides when the search has converged
!| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
!| `calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
!| allocates its own. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| Later stages of the port add the K-study permutation test and the top-level orchestrator that
!| wire these building blocks together.
module tox_data_integration_js_comp_test_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_candidate_bounds_c
    public :: calc_js_comp_test_n_top_k_jsds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds(subroutine)]]
    !| Closed-form from the GAMMA-decay candidate-grid formula
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| uses: `max_n_points_candidate` is exactly the grid's first (largest) `n_points` candidate.
    !| `max_n_neighbors_candidate` is a safe, not necessarily tight, upper bound on the grid's
    !| largest `n_neighbors` candidate -- reached with the smallest `KX_FACTORS` entry at the
    !| smallest `n_points_high` the grid loop ever uses, which by construction never drops below
    !| `n_points_low`.
    subroutine calc_js_comp_test_candidate_bounds_c(&
            max_n_genes_all_studies,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            ierr&
        ) bind(C, name="calc_js_comp_test_candidate_bounds_c")
        use tox_data_integration_js_comp_test_impl, only: calc_js_comp_test_candidate_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Exact upper bound on the `n_points` candidate the grid ever produces
        integer(c_int), intent(out), target :: max_n_neighbors_candidate
            !! Safe upper bound on the `n_neighbors` candidate the grid ever produces
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)

        call calc_js_comp_test_candidate_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate&
        )
    end subroutine calc_js_comp_test_candidate_bounds_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]]
    !| Ported from 125-stabilize-jscomp's inline `n_bootstrapping_top_k_jsds` computation in
    !| `determine_js_comp_test_n_points_n_neighbors_alloc`: sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| to track a two-sided bootstrap confidence interval's endpoints. E.g. for a 95% two-sided
    !| interval (2.5% reserved at each end, the default) with `n_bootstraps=1000`,
    !| `n_top_k = max(1, floor(0.025*1000)) = 25`.
    subroutine calc_js_comp_test_n_top_k_jsds_c(&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            n_top_k,&
            ierr&
        ) bind(C, name="calc_js_comp_test_n_top_k_jsds_c")
        use tox_data_integration_js_comp_test_impl, only: calc_js_comp_test_n_top_k_jsds

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

end module tox_data_integration_js_comp_test_impl_c
#endif
