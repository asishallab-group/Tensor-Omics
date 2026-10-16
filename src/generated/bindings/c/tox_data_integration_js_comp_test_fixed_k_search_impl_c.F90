#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_fixed_k_search_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Fixed-k Parameter Search
!|
!| The data-driven choice of the JS-Comp-Test's fixed-k neighborhoods, `n_points` reference
!| points of `n_neighbors` genes each.
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| walks the `(n_points, n_neighbors)` candidate grid of
!| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| from finest to coarsest and stops at the first admissible candidate whose bootstrapped JSDs
!| have reached a plateau, by the criteria of
!| [[tox_data_integration_js_comp_test_search_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_search_impl(module):check_effect_size_plateau_condition_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_fixed_k_search_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_candidate_bounds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds(subroutine)]]
    !| For the candidate grid of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]:
    !| `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    !| `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`.
    !|
    !| Rejects a non-positive `max_n_genes_all_studies` with invalid input; both outputs are then
    !| left undefined.
    subroutine calc_js_comp_test_candidate_bounds_c(&
            max_n_genes_all_studies,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            ierr&
        ) bind(C, name="calc_js_comp_test_candidate_bounds_c")
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: calc_js_comp_test_candidate_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Exact upper bound on the `n_points` candidate the grid ever produces
        integer(c_int), intent(out), target :: max_n_neighbors_candidate
            !! Safe upper bound on the `n_neighbors` candidate the grid ever produces
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)

        call calc_js_comp_test_candidate_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            ierr = ierr&
        )
    end subroutine calc_js_comp_test_candidate_bounds_c

end module tox_data_integration_js_comp_test_fixed_k_search_impl_c
#endif
