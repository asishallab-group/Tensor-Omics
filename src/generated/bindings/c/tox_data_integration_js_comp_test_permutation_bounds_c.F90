#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_permutation_bounds(module)]]
!| # Work-array capacities of the JS-Comp-Test study-label Permutation-Test
!|
!| The sizing routines of the Permutation-Test that permutes the samples' study labels
!| ([[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module)]] and [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module)]]):
!| [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_js_comp_test_permutation_bounds_impl(interface)]]
!| for the k-NN routines and
!| [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_adaptive_js_comp_test_permutation_bounds_impl(interface)]]
!| for the adaptive ones. The plain entry points of the Permutation-Test call them for you; a
!| caller who allocates the work arrays itself (the expert entry points) calls them first.
module tox_data_integration_js_comp_test_permutation_bounds_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_permutation_bounds_c
    public :: calc_adaptive_js_comp_test_permutation_bounds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_js_comp_test_permutation_bounds(subroutine)]]
    !| Sizes the per-point and per-neighbor work arrays of
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]]
    !| so that both the fixed `(n_points, n_neighbors)` run and a rerun of the parameter search fit:
    !| `max_n_points_candidate = max(n_points, P)` and `max_n_neighbors_candidate = max(n_neighbors,
    !| K)`, with `P` and `K` the bounds of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds_impl(interface)]]
    !| for the same `grid_points_factor` and `grid_min_points`.
    !| `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]] computes it.
    subroutine calc_js_comp_test_permutation_bounds_c(&
            max_n_genes_all_studies,&
            n_points,&
            n_neighbors,&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            grid_points_factor,&
            grid_min_points,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds,&
            ierr&
        ) bind(C, name="calc_js_comp_test_permutation_bounds_c")
        use tox_data_integration_js_comp_test_permutation_bounds, only: calc_js_comp_test_permutation_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points of the fixed-parameter run
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood of the fixed-parameter run
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        real(c_double), intent(in), target :: grid_points_factor
            !! The rerun parameter search's `grid_points_factor`
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(c_int), intent(in), target :: grid_min_points
            !! The rerun parameter search's `grid_min_points`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Capacity for reference points
        integer(c_int), intent(out), target :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        integer(c_int), intent(out), target :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(grid_points_factor)
        M_CHECK_NON_NULL(grid_min_points)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)

        call calc_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
    end subroutine calc_js_comp_test_permutation_bounds_c

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_permutation_bounds(module):calc_adaptive_js_comp_test_permutation_bounds(subroutine)]]
    !| Sizes the per-point work arrays of
    !| [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| and
    !| [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_permutation_impl(interface)]]:
    !| `max_n_points_candidate = max(A, min(N, ceiling(F * N / k_start) + 1))`, with `N =
    !| max_n_genes_all_studies * n_studies` (the padded pool size), `F` the factor
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]]
    !| and `A` the capacity of the adaptive parameter search,
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):calc_adaptive_js_comp_test_bounds_impl(interface)]].
    !| As there, this is a **practical capacity, not a proven bound**: how many reference points a
    !| construction builds is only known once it has run, and a shuffled data set whose
    !| construction builds more fails with invalid input (it is then counted as a failed
    !| permutation). `n_bootstrapping_top_k_jsds` is the bootstrap heap size, exactly as
    !| [[tox_data_integration_stats_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]] computes it.
    !|
    !| Reports invalid input, naming neither argument, for the joint bounds on
    !| `max_n_genes_all_studies` and `n_studies` that the adaptive search's capacity rejects (a pool
    !| size or a first candidate's `k_max` beyond the 32-bit range); the outputs are then left
    !| undefined.
    subroutine calc_adaptive_js_comp_test_permutation_bounds_c(&
            max_n_genes_all_studies,&
            n_studies,&
            k_start,&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            max_n_points_candidate,&
            n_bootstrapping_top_k_jsds,&
            ierr&
        ) bind(C, name="calc_adaptive_js_comp_test_permutation_bounds_c")
        use tox_data_integration_js_comp_test_permutation_bounds, only: calc_adaptive_js_comp_test_permutation_bounds

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_start
            !! Pooled entries every neighborhood of the fixed-knob run takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Practical capacity for reference points
        integer(c_int), intent(out), target :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, invalid input when the pool size or the adaptive search's
            !! first candidate does not fit a 32-bit integer (see above)

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(n_bootstraps)
        M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(n_bootstrapping_top_k_jsds)

        call calc_adaptive_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            k_start = k_start,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
    end subroutine calc_adaptive_js_comp_test_permutation_bounds_c

end module tox_data_integration_js_comp_test_permutation_bounds_c
#endif
