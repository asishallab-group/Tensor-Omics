#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The global JS-Comp-Test of K studies, with a data-driven choice of its neighborhood parameters.
!|
!| Fixed-k neighborhoods (`n_points` reference points of `n_neighbors` genes each):
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
!| searches the `(n_points, n_neighbors)` candidates of
!| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| for a bootstrap plateau, and
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] runs the
!| final test (histograms, JSDs, weighting, permutation test) for the chosen pair.
!|
!| Adaptive neighborhoods, grown per reference point while the residual dispersion stays stable:
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| searches the `(k_start, k_step, k_max)` candidates of
!| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| builds the neighborhoods for the chosen knobs, and
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]] runs
!| the final test on them.
!|
!| Each reference point's histogram bin count is chosen by the occupancy search
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
!| (default) or by Stone's cross-validation
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_stone_impl(interface)]]
!| (`bin_count_mode`). The building blocks are published too: the admissibility gates
!| [[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]], the
!| plateau checks
!| [[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]] and
!| [[tox_data_integration_js_comp_test_impl(module):check_effect_size_plateau_condition_impl(interface)]],
!| the bootstrap
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]], and the
!| brute-force reference
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]].
!| `calc_js_comp_test_candidate_bounds`, `calc_adaptive_js_comp_test_bounds` and
!| `calc_js_comp_test_n_top_k_jsds` size the work arrays for a caller that allocates its own.
module tox_data_integration_js_comp_test_impl_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: calc_js_comp_test_candidate_bounds_c
    public :: calc_adaptive_js_comp_test_bounds_c
    public :: gather_pooled_neighborhood_residuals_c
    public :: calc_js_comp_test_n_top_k_jsds_c

contains

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_candidate_bounds(subroutine)]]
    !| For the candidate grid of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]:
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
        use tox_data_integration_js_comp_test_impl, only: calc_js_comp_test_candidate_bounds

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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_impl(module):calc_adaptive_js_comp_test_bounds(subroutine)]]
    !| Sizes the per-point work arrays of the adaptive parameter search over the candidates of
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| A candidate's reference-point count is only known after its neighborhoods are grown, so this
    !| is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the factor 2 being
    !| [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
    !| `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the smallest
    !| candidate `k_start`. The search rejects a candidate that exceeds it with a capacity status.
    !|
    !| Rejects, with invalid input naming neither argument when the bound is joint, a
    !| `max_n_genes_all_studies` or `n_studies` below 1, a product `N` beyond the 32-bit range, or a
    !| first candidate `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))` above
    !| `huge(1_int32)/4 = 536870911` (so that its `k_max = 4*k_start_1` fits); the output is then
    !| left undefined.
    subroutine calc_adaptive_js_comp_test_bounds_c(&
            max_n_genes_all_studies,&
            n_studies,&
            max_n_points_candidate,&
            ierr&
        ) bind(C, name="calc_adaptive_js_comp_test_bounds_c")
        use tox_data_integration_js_comp_test_impl, only: calc_adaptive_js_comp_test_bounds

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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_impl(module):gather_pooled_neighborhood_residuals(subroutine)]]
    !| Gathers one reference point's residuals from every neighbor gene of every study into one flat
    !| array, given that point's per-study neighbor gene indices (one column of the
    !| `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
    !| `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| feed to the bin-count search, so a caller can pass the same input to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
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
        use tox_data_integration_js_comp_test_impl, only: gather_pooled_neighborhood_residuals

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

    !> summary: C-wrapper for [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(subroutine)]]
    !| Sizes the top-k/bottom-k heaps
    !| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] uses
    !| for the interval's endpoints: `n_top_k = max(1, floor(level/100 * n_bootstraps))`, e.g. 25 for
    !| the default 2.5 and `n_bootstraps=1000`.
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
