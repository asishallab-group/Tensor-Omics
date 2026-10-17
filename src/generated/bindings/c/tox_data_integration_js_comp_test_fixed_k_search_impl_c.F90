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
!| By default the grid holds only the large neighborhoods
!| ([[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]]);
!| `neighborhood_size_mode` selects the small ones
!| ([[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]])
!| or both
!| ([[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_BOTH(variable)]]).
!| `calc_js_comp_test_candidate_bounds` sizes the work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_fixed_k_search_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_char, c_double, c_int, c_loc
    use tox_conversions, only: c_char_as_view
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_candidate_bounds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):calc_js_comp_test_candidate_bounds(subroutine)]]
    !| For the candidate grid of
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):generate_js_comp_test_candidates_impl(interface)]]:
    !| `max_n_points_candidate` is exactly its largest `n_points` candidate, and
    !| `max_n_neighbors_candidate` a safe (not necessarily tight) bound on its largest `n_neighbors`,
    !| in every `neighborhood_size_mode` (the bound uses the large-neighborhood factor `KX_FACTORS(1)`)
    !| and every `grid_floor_mode` (the floor candidate of the clamp sits at the floor this bound
    !| already uses). Pass the same `grid_sample_size_mode`, `grid_points_factor` and
    !| `grid_min_points` as the parameter search that is sized with the result.
    !|
    !| Rejects a non-positive `max_n_genes_all_studies` or `n_studies`, an unknown
    !| `grid_sample_size_mode`, a non-positive `grid_points_factor` and a `grid_min_points` outside
    !| `[1, MAX_POINTS]` with invalid input (a non-finite factor with `ERR_NAN_INF`); both outputs
    !| are then left undefined.
    subroutine calc_js_comp_test_candidate_bounds_c(&
            max_n_genes_all_studies,&
            n_studies,&
            max_n_points_candidate,&
            max_n_neighbors_candidate,&
            grid_sample_size_mode,&
            grid_points_factor,&
            grid_min_points,&
            ierr&
        ) bind(C, name="calc_js_comp_test_candidate_bounds_c")
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: calc_js_comp_test_candidate_bounds
        use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_GRID_SAMPLE_SIZE_PER_STUDY, MODE_GRID_SAMPLE_SIZE_POOLED

        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies; used only by the pooled grid sample size
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(out), target :: max_n_points_candidate
            !! Exact upper bound on the `n_points` candidate the grid ever produces
        integer(c_int), intent(out), target :: max_n_neighbors_candidate
            !! Safe upper bound on the `n_neighbors` candidate the grid ever produces
        character(len=1, kind=c_char), dimension(26), intent(in), target :: grid_sample_size_mode
            !! Sample size `N` in the grid's first reference-point count; must match the parameter
            !! search's own. An experiment switch
            !!
            !! | Mode                                                         | Value                                                                                                       |
            !! |--------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------|
            !! | One study's gene count, `max_n_genes_all_studies` (default)  | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_GRID_SAMPLE_SIZE_PER_STUDY(variable)]] |
            !! | The pooled gene count, `max_n_genes_all_studies * n_studies` | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_GRID_SAMPLE_SIZE_POOLED(variable)]]    |
            !! The default value is `'grid_sample_size_per_study'`.
        real(c_double), intent(in), target :: grid_points_factor
            !! Factor in the grid's first reference-point count `ceiling(grid_points_factor * sqrt(N))`;
            !! must match the parameter search's own. An experiment switch
            !! The minimum valid value is `above(0.0_real64)`.
            !! The default value is `4.0_real64`.
        integer(c_int), intent(in), target :: grid_min_points
            !! The grid's floor; must match the parameter search's own. An experiment switch: 300
            !! (`MIN_POINTS`) is the floor before it existed
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `MAX_POINTS`.
            !! The default value is `200_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure
        integer(int32) :: grid_sample_size_mode_mode_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_points_candidate)
        M_CHECK_NON_NULL(max_n_neighbors_candidate)
        M_CHECK_NON_NULL(grid_points_factor)
        M_CHECK_NON_NULL(grid_min_points)
        M_CHECK_ARRAY_NON_NULL(grid_sample_size_mode, 26)

        block
            character(len=:), pointer :: grid_sample_size_mode_f
            grid_sample_size_mode_f => c_char_as_view(grid_sample_size_mode)

            select case (grid_sample_size_mode_f)
                case ("grid_sample_size_per_study")
                    grid_sample_size_mode_mode_f = MODE_GRID_SAMPLE_SIZE_PER_STUDY
                case ("grid_sample_size_pooled")
                    grid_sample_size_mode_mode_f = MODE_GRID_SAMPLE_SIZE_POOLED
                case default
                    call set_err(ierr, ERR_INVALID_INPUT)
                    return
            end select
        end block

        call calc_js_comp_test_candidate_bounds(&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_studies = n_studies,&
            max_n_points_candidate = max_n_points_candidate,&
            max_n_neighbors_candidate = max_n_neighbors_candidate,&
            grid_sample_size_mode = grid_sample_size_mode_mode_f,&
            grid_points_factor = grid_points_factor,&
            grid_min_points = grid_min_points,&
            ierr = ierr&
        )
    end subroutine calc_js_comp_test_candidate_bounds_c

end module tox_data_integration_js_comp_test_fixed_k_search_impl_c
#endif
