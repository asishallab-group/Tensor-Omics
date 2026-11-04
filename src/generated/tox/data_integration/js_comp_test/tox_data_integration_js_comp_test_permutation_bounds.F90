#include <src/macros.h>

!> # Work-array capacities of the JS-Comp-Test study-label Permutation-Test
!|
!| The sizing routines of the Permutation-Test that permutes the samples' study labels
!| ([[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module)]] and [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module)]]):
!| [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_js_comp_test_permutation_bounds_impl(interface)]]
!| for the k-NN routines and
!| [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_adaptive_js_comp_test_permutation_bounds_impl(interface)]]
!| for the adaptive ones. The plain entry points of the Permutation-Test call them for you; a
!| caller who allocates the work arrays itself (the expert entry points) calls them first.
!|
!| Generated from [[tox_data_integration_js_comp_test_permutation_bounds_impl(module)]]; do not edit -- regenerate instead.
module tox_data_integration_js_comp_test_permutation_bounds
    use tox_data_integration_js_comp_test_permutation_bounds_impl, only: calc_adaptive_js_comp_test_permutation_bounds_impl, calc_js_comp_test_permutation_bounds_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_math_impl, only: above
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MAX_POINTS
    use tox_errors, only: set_ok, is_err, clear_err_arg_pos, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_permutation_bounds
    public :: calc_adaptive_js_comp_test_permutation_bounds

contains

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_js_comp_test_permutation_bounds_impl]].
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
    pure subroutine calc_js_comp_test_permutation_bounds(&
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
        )
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the fixed-parameter run
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the fixed-parameter run
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `max_n_genes_all_studies`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        real(real64), intent(in), optional :: grid_points_factor
            !! The rerun parameter search's `grid_points_factor`
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(int32), intent(in), optional :: grid_min_points
            !! The rerun parameter search's `grid_min_points`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(int32), intent(out) :: max_n_points_candidate
            !! Capacity for reference points
        integer(int32), intent(out) :: max_n_neighbors_candidate
            !! Capacity for neighbors per neighborhood
        integer(int32), intent(out) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_points, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=3_int32, min=1_int32, max=max_n_genes_all_studies)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=5_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_real(grid_points_factor, ierr, arg_pos=6_int32, min=above(0.0_real64))
        call validate_in_range_int(grid_min_points, ierr, arg_pos=7_int32, min=1_int32, max=MAX_POINTS)
        if (is_err(ierr)) return
#endif

        call calc_js_comp_test_permutation_bounds_impl(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_points = n_points,&
            n_neighbors = n_neighbors,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds&
        )
    end subroutine calc_js_comp_test_permutation_bounds

    !> summary: Validates its inputs, then calls [[tox_data_integration_js_comp_test_permutation_bounds_impl(module):calc_adaptive_js_comp_test_permutation_bounds_impl]].
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
    pure subroutine calc_adaptive_js_comp_test_permutation_bounds(&
            max_n_genes_all_studies,&
            n_studies,&
            k_start,&
            n_bootstraps,&
            two_sided_bootstrapping_significance_level,&
            max_n_points_candidate,&
            n_bootstrapping_top_k_jsds,&
            ierr&
        )
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood of the fixed-knob run takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps of the parameter search, when it is rerun
            !! The minimum valid value is `1_int32`.
            !! The default value is `1000_int32`.
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Two-sided significance level of the bootstrap interval, in percent per end
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
            !! The default value is `2.5_real64`.
        integer(int32), intent(out) :: max_n_points_candidate
            !! Practical capacity for reference points
        integer(int32), intent(out) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size (values kept at each end of the bootstrap distribution)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, invalid input when the pool size or the adaptive search's
            !! first candidate does not fit a 32-bit integer (see above)

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(max_n_genes_all_studies, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_studies, ierr, arg_pos=2_int32, min=1_int32)
        call validate_in_range_int(k_start, ierr, arg_pos=3_int32, min=1_int32)
        call validate_in_range_int(n_bootstraps, ierr, arg_pos=4_int32, min=1_int32)
        call validate_in_range_real(two_sided_bootstrapping_significance_level, ierr, arg_pos=5_int32, min=0.0_real64, max=100.0_real64)
        if (is_err(ierr)) return
#endif

        call calc_adaptive_js_comp_test_permutation_bounds_impl(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            k_start = k_start,&
            n_bootstraps = n_bootstraps,&
            two_sided_bootstrapping_significance_level = two_sided_bootstrapping_significance_level,&
            max_n_points_candidate = max_n_points_candidate,&
            n_bootstrapping_top_k_jsds = n_bootstrapping_top_k_jsds,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine calc_adaptive_js_comp_test_permutation_bounds

end module tox_data_integration_js_comp_test_permutation_bounds
