#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_preprocessing_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| The step that turns expression vectors into the neighborhood residuals the rest of the test
!| consumes: gene-wise means, the signed deviation of each replicate from them, and the
!| neighborhoods of reference points those residuals are grouped into so the comparison is
!| conditioned on expression level rather than pooled across it.
!|
!| Besides neighborhoods of a fixed size, the adaptive construction
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| grows each reference point's neighborhood while its residual dispersion stays stable, and
!| [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]
!| checks that consecutive neighborhoods overlap enough.
!|
!| `calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own, and
!| `gather_pooled_neighborhood_residuals` gathers a neighborhood's residuals of all studies.
module tox_data_integration_preprocessing_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_neighborhood_size_c
    public :: gather_pooled_neighborhood_residuals_c

contains

    !> summary: C-wrapper for [[tox_data_integration_preprocessing_impl(module):calc_neighborhood_size(function)]]
    !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower
    !| due to few genes with non-NaN mean.
    subroutine calc_neighborhood_size_c(&
            n_pool,&
            n_points,&
            n_genes_S,&
            mean_S,&
            desired_size,&
            n_neighbors,&
            ierr&
        ) bind(C, name="calc_neighborhood_size_c")
        use tox_data_integration_preprocessing_impl, only: calc_neighborhood_size

        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_pool
            !! Total number of pooled mean-expression values across both studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: desired_size
            !! Optional desired neighborhood size
            !! The default value is `1000_int32`.
        integer(c_int), intent(out), target :: n_neighbors
            !! Calculated neighborhood size
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_genes_S)
        M_CHECK_NON_NULL(desired_size)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(mean_S, n_genes_S)

        n_neighbors = calc_neighborhood_size(&
            n_pool = n_pool,&
            n_points = n_points,&
            n_genes_S = n_genes_S,&
            mean_S = mean_S,&
            desired_size = desired_size&
        )
    end subroutine calc_neighborhood_size_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing_impl(module):gather_pooled_neighborhood_residuals(subroutine)]]
    !| Gathers one reference point's residuals from every neighbor gene of every study into one flat
    !| array, given that point's per-study neighbor gene indices (one column of the
    !| `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
    !| `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| feed to the bin-count search, so a caller can pass the same input to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
    !| on real data.
    subroutine gather_pooled_neighborhood_residuals_c(&
            residuals,&
            max_n_reps_all_studies,&
            max_n_genes_all_studies,&
            n_neighbors,&
            n_studies,&
            neighborhood_indices_point,&
            pooled_residuals,&
            ierr&
        ) bind(C, name="gather_pooled_neighborhood_residuals_c")
        use tox_data_integration_preprocessing_impl, only: gather_pooled_neighborhood_residuals

        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors per neighborhood
        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Matrix of signed residuals per study, NaN allowed for missing values; not validated
        integer(c_int), dimension(n_neighbors, n_studies), intent(in), target :: neighborhood_indices_point
            !! Gene indices of one reference point's neighborhood, per study; each in
            !! [1, max_n_genes_all_studies]
        real(c_double), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: pooled_residuals
            !! This reference point's pooled residuals, laid out as a
            !! (max_n_reps_all_studies, n_neighbors, n_studies) array
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices_point, n_neighbors * n_studies)
        M_CHECK_ARRAY_NON_NULL(pooled_residuals, (max_n_reps_all_studies*n_neighbors*n_studies))

        call gather_pooled_neighborhood_residuals(&
            residuals = residuals,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            n_neighbors = n_neighbors,&
            n_studies = n_studies,&
            neighborhood_indices_point = neighborhood_indices_point,&
            pooled_residuals = pooled_residuals,&
            ierr = ierr&
        )
    end subroutine gather_pooled_neighborhood_residuals_c

end module tox_data_integration_preprocessing_impl_c
#endif
