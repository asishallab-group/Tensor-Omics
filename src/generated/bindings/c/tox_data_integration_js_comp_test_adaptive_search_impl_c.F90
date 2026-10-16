#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_adaptive_search_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Adaptive Parameter Search
!|
!| The data-driven choice of the JS-Comp-Test's adaptive neighborhoods, each grown from `k_start`
!| pooled entries in rounds of `k_step` up to `k_max` while the residual dispersion stays stable.
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| walks the `(k_start, k_step, k_max)` candidates of
!| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
!| and stops at the first admissible candidate whose bootstrapped JSDs have reached a plateau, by
!| the criteria of
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]];
!| every candidate it tries is logged with an `ADAPTIVE_CANDIDATE_*` status.
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| builds the neighborhoods for the chosen knobs. `calc_adaptive_js_comp_test_bounds` sizes the
!| work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_adaptive_search_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_adaptive_js_comp_test_bounds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds(subroutine)]]
    !| Sizes the per-point work arrays of the adaptive parameter search over the candidates of
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| A candidate's reference-point count is only known after its neighborhoods are grown, so this
    !| is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
    !| `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
    !| candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.
    !|
    !| Rejects, with invalid input naming neither argument when the bound is joint, a
    !| `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
    !| first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
    !| `floor(huge(1_int32)/ADAPTIVE_K_MAX_FACTOR)` ([[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_K_MAX_FACTOR(variable)]]), so that its
    !| `k_max = ceiling(ADAPTIVE_K_MAX_FACTOR*k_start_1)` fits; the output is then left undefined.
    subroutine calc_adaptive_js_comp_test_bounds_c(&
            max_n_genes_all_studies,&
            n_studies,&
            max_n_points_candidate,&
            ierr&
        ) bind(C, name="calc_adaptive_js_comp_test_bounds_c")
        use tox_data_integration_js_comp_test_adaptive_search_impl, only: calc_adaptive_js_comp_test_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Practical upper bound on the number of reference points any adaptive candidate is
            !! allowed to emerge with
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_points_candidate)

        call calc_adaptive_js_comp_test_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            max_n_points_candidate = max_n_points_candidate,&
            ierr = ierr&
        )
    end subroutine calc_adaptive_js_comp_test_bounds_c

end module tox_data_integration_js_comp_test_adaptive_search_impl_c
#endif
