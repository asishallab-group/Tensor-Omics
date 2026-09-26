#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_js_comp_test_impl(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter Search
!|
!| The data-driven `(n_points, n_neighbors)` parameter-stabilization search this pipeline runs
!| before the JSD-Comp-Test proper (Issue #126): a GAMMA-decay candidate grid
!| ([[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
!| generates candidate `(n_points, n_neighbors)` pairs only -- each candidate's real
!| per-neighborhood histogram bin count is decided later, per reference point, by
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]),
!| two admissibility gates a candidate must pass before it is bootstrapped
!| ([[tox_data_integration_js_comp_test_impl(module):check_neighborhood_overlaps_impl(interface)]],
!| [[tox_data_integration_js_comp_test_impl(module):check_mean_pmf_min_counts_impl(interface)]]),
!| and the plateau check that decides when the search has converged
!| ([[tox_data_integration_js_comp_test_impl(module):check_plateau_condition_impl(interface)]]).
!| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
!| is a brute-force reference implementation of the same per-point bin-count search, exhaustively
!| testing every candidate `M` instead of the fast geometric-search-then-refinement the production
!| routine uses, for validating that the fast search's own result is correct.
!| `calc_js_comp_test_candidate_bounds` sizes the candidate-grid work arrays for a caller that
!| allocates its own. For adaptive (heteroscedastic, Issue #217) neighborhood construction,
!| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]]
!| generates the ascending `(k_start, k_step, k_max)` growth-knob candidates instead, and
!| `calc_adaptive_js_comp_test_bounds` recommends the reference-point capacity their search's
!| per-point arrays are sized by;
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
!| is that search, the adaptive counterpart of the fixed-k parameter search, with the same
!| plateau criteria. For the selected knobs,
!| [[tox_data_integration_js_comp_test_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| rebuilds the reference points and their neighborhoods, pooled across all studies, and
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_adaptive_impl(interface)]]
!| runs the final test on them. Together they reproduce exactly what the search traced for that
!| candidate. Once a candidate has passed both gates, its bootstrap confidence interval
!| is resampled from the pooled consensus histogram by
!| [[tox_data_integration_js_comp_test_impl(module):bootstrap_histogram_impl(interface)]] (heap
!| size recommended by
!| [[tox_data_integration_js_comp_test_impl(module):calc_js_comp_test_n_top_k_jsds(interface)]]).
!| For fixed-k neighborhoods, the final test itself (histograms, JSDs, weighting and the K-study
!| permutation test) is
!| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]].
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
    !| Closed-form from the GAMMA-decay candidate-grid formula
    !| [[tox_data_integration_js_comp_test_impl(module):generate_js_comp_test_candidates_impl(interface)]]
    !| uses: `max_n_points_candidate` is exactly the grid's first (largest) `n_points` candidate.
    !| `max_n_neighbors_candidate` is a safe, not necessarily tight, upper bound on the grid's
    !| largest `n_neighbors` candidate -- reached with the smallest `KX_FACTORS` entry at the
    !| smallest `n_points_high` the grid loop ever uses, which by construction never drops below
    !| `n_points_low`.
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
    !| Sizes the per-point work arrays of the adaptive (Issue #217) parameter search, whose
    !| candidates come from
    !| [[tox_data_integration_js_comp_test_impl(module):generate_adaptive_js_comp_test_candidates_impl(interface)]].
    !| An adaptive candidate's reference-point count is not known before its neighborhoods are
    !| grown, so this is a **practical capacity, not a proven bound**:
    !| `max_n_points_candidate = min(N, ceiling(2 * N / k_start_last) + 1)`, with the safety factor
    !| 2 being [[tox_data_integration_js_comp_test_impl(module):ADAPTIVE_POINT_CAPACITY_FACTOR(variable)]],
    !| `N = max_n_genes_all_studies * n_studies` the padded pool size and `k_start_last` the
    !| smallest `k_start` the candidate sequence contains (the candidate expected to emerge with
    !| the most reference points). The only provable bound is `N` itself, because every new
    !| reference point's seed strictly advances through the pool; sizing every per-point array by
    !| `N` would be prohibitive, so a candidate that emerges with more points than this capacity is
    !| rejected by the adaptive search with a capacity status rather than stored.
    !|
    !| Computed in 64-bit integer / double precision throughout, so `N` itself never overflows.
    !| Enforces the candidate generator's representability bound, which the generator itself
    !| cannot check: the first candidate's `k_start_1 = max(10, 2*n_studies, ceiling(0.02 * N))`
    !| must not exceed `huge(1_int32)/4 = 536870911`, so that its `k_max = 4*k_start_1` fits a
    !| 32-bit integer; otherwise this routine reports invalid input. The bound is joint in both
    !| arguments, so the error names neither.
    !|
    !| Rejects a `max_n_genes_all_studies` or `n_studies` below 1, or a pair beyond that
    !| representability bound, with invalid input; the output is then left undefined.
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
    !| Given one reference point's own per-study neighbor gene indices (one column of a larger
    !| `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)`, as produced by
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]),
    !| gathers that point's residual values from every neighbor gene, across every study, into one
    !| flat pooled array. This is the exact same pooling
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| perform internally, per reference point, before handing the result to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
    !| own occupancy search -- published so a caller can reconstruct that exact same input directly
    !| on real data and feed it to
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
    !| (or to `determine_bin_count_occupancy` itself), to check whether the fast search and the
    !| exhaustive reference ever actually disagree in practice, not just on a synthetic fixture.
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
            !! Matrix of signed residuals per study, NaN explicitly allowed for missing values;
            !! unvalidated on this export path, so any value passes through as-is
        integer(c_int), dimension(n_neighbors, n_studies), intent(in), target :: neighborhood_indices_point
            !! Gene indices of one reference point's neighborhood, per study -- one column of a
            !! larger neighborhood_indices_all_studies(n_neighbors, n_points, n_studies), as sliced
            !! by the caller
        real(c_double), dimension(max_n_reps_all_studies*n_neighbors*n_studies), intent(out), target :: pooled_residuals
            !! The pooled residual values for this reference point, across every neighbor and every
            !! study, laid out exactly as a (max_n_reps_all_studies, n_neighbors, n_studies) array
            !! would be
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
