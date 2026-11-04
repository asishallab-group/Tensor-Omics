#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Work-array capacities of the JS-Comp-Test study-label Permutation-Test
!|
!| The sizing routines of the Permutation-Test that permutes the samples' study labels
!| ([[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module)]] and [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module)]]):
!| [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_js_comp_test_permutation_bounds_impl(interface)]]
!| for the k-NN routines and
!| [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_adaptive_js_comp_test_permutation_bounds_impl(interface)]]
!| for the adaptive ones. The plain entry points of the Permutation-Test call them for you; a
!| caller who allocates the work arrays itself (the expert entry points) calls them first.
module tox_data_integration_js_comp_test_permutation_bounds_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: is_err
    use tox_data_integration_stats_impl, only: bootstrap_heap_size
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: calc_js_comp_test_candidate_bounds_impl, MAX_POINTS
    use tox_data_integration_js_comp_test_adaptive_search_impl, only: calc_adaptive_js_comp_test_bounds_impl, &
                                                                      adaptive_point_capacity
    M_IMPLICIT_NONE
    private
    public :: calc_js_comp_test_permutation_bounds_impl, calc_adaptive_js_comp_test_permutation_bounds_impl

contains

    !> summary: Recommend the work-array capacities of the k-NN study-label Permutation-Test
    !| AUTHOR_LASZLO_LANG
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
    pure subroutine calc_js_comp_test_permutation_bounds_impl(max_n_genes_all_studies, n_points, n_neighbors, &
                                                              n_bootstraps, two_sided_bootstrapping_significance_level, &
                                                              grid_points_factor, grid_min_points, &
                                                              max_n_points_candidate, max_n_neighbors_candidate, &
                                                              n_bootstrapping_top_k_jsds)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the fixed-parameter run
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the fixed-parameter run
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        real(real64), intent(in), optional :: grid_points_factor
            !! The rerun parameter search's `grid_points_factor`
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_GRID_POINTS_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: grid_min_points
            !! The rerun parameter search's `grid_min_points`
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_POINTS)
            !! DM_DEFAULT(M_GRID_MIN_POINTS_DEFAULT)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Capacity for reference points
        integer(int32), intent(out) :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        integer(int32), intent(out) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)

        integer(int32) :: actual_n_bootstraps, search_n_points, search_n_neighbors
        real(real64) :: actual_significance_level

        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(two_sided_bootstrapping_significance_level, actual_significance_level, M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)

        call calc_js_comp_test_candidate_bounds_impl(max_n_genes_all_studies, search_n_points, search_n_neighbors, &
                                                     grid_points_factor, grid_min_points)
        max_n_points_candidate = max(n_points, search_n_points)
        max_n_neighbors_candidate = max(n_neighbors, search_n_neighbors)
        n_bootstrapping_top_k_jsds = bootstrap_heap_size(actual_n_bootstraps, actual_significance_level)
    end subroutine calc_js_comp_test_permutation_bounds_impl

    !> summary: Recommend the work-array capacity of the adaptive study-label Permutation-Test
    !| AUTHOR_LASZLO_LANG
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
    pure subroutine calc_adaptive_js_comp_test_permutation_bounds_impl(max_n_genes_all_studies, n_studies, k_start, &
                                                                       n_bootstraps, &
                                                                       two_sided_bootstrapping_significance_level, &
                                                                       max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                                                       ierr)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood of the fixed-knob run takes unconditionally
            !! DM_MIN(1_int32)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(out) :: max_n_points_candidate
            !! Practical capacity for reference points
        integer(int32), intent(out) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, invalid input when the pool size or the adaptive search's
            !! first candidate does not fit a 32-bit integer (see above)

        integer(int32) :: actual_n_bootstraps, search_capacity
        real(real64) :: actual_significance_level

        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(two_sided_bootstrapping_significance_level, actual_significance_level, M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)

        ! Its joint bounds name no argument, so none is named here either.
        call calc_adaptive_js_comp_test_bounds_impl(max_n_genes_all_studies, n_studies, search_capacity, ierr)
        if (is_err(ierr)) return

        ! The same capacity formula as the search's, at the fixed-knob run's k_start.
        max_n_points_candidate = max(search_capacity, &
                                     adaptive_point_capacity(max_n_genes_all_studies, n_studies, k_start))
        n_bootstrapping_top_k_jsds = bootstrap_heap_size(actual_n_bootstraps, actual_significance_level)
    end subroutine calc_adaptive_js_comp_test_permutation_bounds_impl
end module tox_data_integration_js_comp_test_permutation_bounds_impl
