#include "macros.h"

!> In multi-study omics analyses, it is often unclear whether biological replicates originating from different studies can be safely treated as sampling the same biological condition.
!| Even when studies nominally target the same tissue and condition, differences in sample handling, sequencing technologies, preprocessing pipelines, or cohort, 
!| composition can introduce batch effects that are not easily detectable from mean expression levels alone.
!|
!| This ambiguity has direct consequences for downstream analyses in Tensor Omics. Integrating incompatible replicate sets can:
!|
!|  - distort expression spaces,
!|  - affect distance-based analyses,
!|  - bias machine learning models,
!|
!| while unnecessarily separating compatible datasets reduces statistical power.
!|
!| To address this, we introduce a Jensen–Shannon-Divergence based compatibility test (JSD-Comp-Test)
!| that empirically evaluates whether two sets of biological replicates exhibit comparable replicate-level variability.
!| Rather than comparing mean expression values, the method focuses on the distribution of signed residuals (replicate deviations from the gene-wise mean),
!| conditioned on mean expression levels to account for heteroscedasticity, which is a well-known property of omics data.
!|
!| The goal of this issue is to define, implement, and validate this compatibility test as a diagnostic tool that can be applied prior to data integration.
!| The test is intended to support principled decisions
!| on whether replicate sets from different studies should be merged or treated as distinct conditions within Tensor Omics workflows.
module tox_data_integration
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_data_integration_preprocessing, only: compute_gene_means_impl => compute_gene_means,&
        compute_residuals_impl => compute_residuals,&
        pool_means_alloc_impl => pool_means_alloc,&
        pool_means_impl => pool_means,&
        calc_neighborhood_size_impl => calc_neighborhood_size,&
        construct_neighborhoods_alloc_impl => construct_neighborhoods_alloc,&
        construct_neighborhoods_impl => construct_neighborhoods
    use tox_data_integration_jsd, only: build_residual_histograms_impl => build_residual_histograms,&
        compute_divergence_per_reference_point_impl => compute_divergence_per_reference_point,&
        compute_weighted_global_divergence_impl => compute_weighted_global_divergence
    use tox_data_integration_js_comp_test, only: determine_js_comp_test_n_points_n_neighbors_alloc_impl => determine_js_comp_test_n_points_n_neighbors_alloc, &
        js_comp_test_alloc_impl => js_comp_test_alloc, &
        JOIN_MIN, JOIN_MAX, JOIN_MEDIAN
    implicit none

    interface compute_gene_means
        module procedure compute_gene_means_impl
    end interface compute_gene_means

    interface compute_residuals
        module procedure compute_residuals_impl
    end interface compute_residuals

    interface pool_means_alloc
        module procedure pool_means_alloc_impl
    end interface pool_means_alloc

    interface pool_means
        module procedure pool_means_impl
    end interface pool_means

    interface calc_neighborhood_size
        module procedure calc_neighborhood_size_impl
    end interface calc_neighborhood_size

    interface construct_neighborhoods_alloc
        module procedure construct_neighborhoods_alloc_impl
    end interface construct_neighborhoods_alloc

    interface construct_neighborhoods
        module procedure construct_neighborhoods_impl
    end interface construct_neighborhoods

    interface determine_js_comp_test_n_points_n_neighbors_alloc
        module procedure determine_js_comp_test_n_points_n_neighbors_alloc_impl
    end interface determine_js_comp_test_n_points_n_neighbors_alloc

    interface js_comp_test_alloc
        module procedure js_comp_test_alloc_impl
    end interface js_comp_test_alloc

    interface build_residual_histograms
        module procedure build_residual_histograms_impl
    end interface build_residual_histograms

    interface compute_divergence_per_reference_point
        module procedure compute_divergence_per_reference_point_impl
    end interface compute_divergence_per_reference_point

    interface compute_weighted_global_divergence
        module procedure compute_weighted_global_divergence_impl
    end interface compute_weighted_global_divergence

contains

    pure subroutine get_join_method(c_method_str, join_method, ierr)
        use, intrinsic :: iso_c_binding, only: c_char
        use tox_conversions, only: c_char_1d_as_string
        use tox_errors, only: ERR_INVALID_INPUT, set_ok, is_err, set_err
        character(len=1, kind=c_char), dimension(6), intent(in) :: c_method_str
            !! join method string ("min", "max", "median")
        integer(int32), intent(out) :: join_method
            !! integer representation for the join method passed by `c_method_str`
        integer(int32), intent(out) :: ierr
            !! Error code

        character(len=:), allocatable :: mode_str_f

        call set_ok(ierr)

        call c_char_1d_as_string(c_method_str, mode_str_f, ierr)
        if (is_err(ierr)) return

        select case (trim(mode_str_f))
            case ("min")
                join_method = JOIN_MIN
            case ("max")
                join_method = JOIN_MAX
            case ("median")
                join_method = JOIN_MEDIAN
            case default
                call set_err(ierr, ERR_INVALID_INPUT)
        end select
    end subroutine get_join_method
end module tox_data_integration

!> C-compatible wrapper for [[tox_data_integration(module):compute_gene_means(interface)]]
pure subroutine compute_gene_means_c(expr, n_genes, n_reps, means, max_n_genes_all_studies, ierr) &
        bind(C, name="compute_gene_means_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: compute_gene_means
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes in the study
    integer(c_int), intent(in), target :: n_reps
        !! Number of biological replicates in the study
    integer(c_int), intent(in), target :: max_n_genes_all_studies
        !! Maximum number of genes across all studies
    real(c_double), intent(in), target :: expr(n_reps, n_genes)
        !! Expression matrix
    real(c_double), intent(out), target :: means(max_n_genes_all_studies)
        !! Per-gene mean expression values
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_reps)
    M_CHECK_NON_NULL(max_n_genes_all_studies)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(means)

    call compute_gene_means(expr, n_genes, n_reps, means, max_n_genes_all_studies, ierr)

end subroutine compute_gene_means_c

!> C-compatible wrapper for [[tox_data_integration(module):compute_residuals(interface)]]
pure subroutine compute_residuals_c(expr, n_genes, n_reps, means, max_n_genes_all_studies, &
        max_n_reps_all_studies, resid, ierr) bind(C, name="compute_residuals_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: compute_residuals
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes
        !! Number of genes in the study
    integer(c_int), intent(in), target :: n_reps
        !! Number of biological replicates in the study
    integer(c_int), intent(in), target :: max_n_genes_all_studies
        !! Maximum number of genes across all studies
    integer(c_int), intent(in), target :: max_n_reps_all_studies
        !! Maximum number of replicates across all studies
    real(c_double), intent(in), target :: expr(n_reps, n_genes)
        !! Expression matrix containing
    real(c_double), intent(in), target :: means(max_n_genes_all_studies)
        !! Per-gene mean expression values
    real(c_double), intent(out), target :: resid(max_n_reps_all_studies, max_n_genes_all_studies)
        !! Matrix of signed residuals
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_reps)
    M_CHECK_NON_NULL(max_n_genes_all_studies)
    M_CHECK_NON_NULL(max_n_reps_all_studies)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(means)
    M_CHECK_NON_NULL(resid)

    call compute_residuals(expr, n_genes, n_reps, means, max_n_genes_all_studies, &
        max_n_reps_all_studies, resid, ierr)

end subroutine compute_residuals_c

!> C-compatible wrapper for [[tox_data_integration(module):determine_js_comp_test_n_points_n_neighbors_alloc(interface)]]
subroutine determine_js_comp_test_n_points_n_neighbors_c( &
        n_points, n_neighbors, residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
        shared_residual_range, n_bins, gene_means, n_studies, n_bootstraps, &
        best_candidate_pair_confidence_interval, join_method, ierr, &
        min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, &
        two_sided_bootstrapping_significance_level, random_seed, residual_range_quantile) &
    bind(C, name="determine_js_comp_test_n_points_n_neighbors_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_integration, only: determine_js_comp_test_n_points_n_neighbors_alloc, get_join_method
    use tox_errors, only: is_err
    M_USE_NULL_VALIDATION

    integer(c_int), intent(out), target :: n_points
        !! The finally chosen candidate from `candidates_n_points`
    integer(c_int), intent(out), target :: n_neighbors
        !! The finally chosen candidate from `candidates_n_neighbors`
    integer(c_int), intent(in), target :: max_n_reps_all_studies
        !! Maximum number of replicates across all studies
    integer(c_int), intent(in), target :: max_n_genes_all_studies
        !! Maximum number of genes across all studies
    integer(c_int), intent(in), target :: n_studies
        !! Number of studies
    real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
        !! Matrix of signed residuals per study
    real(c_double), intent(out), target :: shared_residual_range
        !! Computed residual range (R)
    integer(c_int), intent(out), target :: n_bins
        !! The number of bins used for the finally chosen candidate from `candidates_n_points_n_neighbors`
    real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
        !! Per-gene mean expression values for all studies
    integer(c_int), intent(in), target :: n_bootstraps
        !! Number of bootstraps to perform for a candidate pair
    real(c_double), dimension(2, n_studies), intent(out), target :: best_candidate_pair_confidence_interval
        !! The JSD Confidence Intervals from bootstrapping for the best candidate pair. `-1.0_real64` if no candidate pair succeeded and fallback to `n_points=candidates_n_points(1)` and `n_neighbors=candidates_n_neighbors(1)`
    character(len=1, kind=c_char), dimension(6), intent(in), target :: join_method
        !! The way to evaluate all studies' confidence intervals for candidate determination
        !!
        !! 1. "min": take min overlap of all studies' CI overlaps -> succeeds only if `all(ci_overlaps > min_neighbor_overlap)`
        !! 2. "max": take max overlap of all studies' CI overlaps -> succeeds only if `any(ci_overlaps > min_neighbor_overlap)`
        !! 3. "median": take median overlap of all studies' CI overlaps -> succeeds only if `count(ci_overlaps > min_neighbor_overlap) >= (n_studies - 1) / 2 + 1`
    integer(c_int), intent(in), target :: min_count_per_mean_bin
        !! Number of minimum residuals a bin should have in the mean pmf to make a candidate pair eligible, recommended: `5`
    real(c_double), intent(in), target :: min_neighbor_overlap
        !! Minimum fractional overlap in genes a neighborhood have to its succeeding neighborhood to make a candidate pair eligible, recommended: `0.1`
    real(c_double), intent(in), target :: succeeding_ci_overlap
        !! Minimum fractional overlap the confidence intervals should have to the current best confidence intervals (respecting the `join_method`) to make a candidate pair eligible, recommended: `0.9`
    integer(c_int), intent(in), target :: random_seed
        !! Random seed to use for random number generation
    real(c_double), intent(in), target :: two_sided_bootstrapping_significance_level
        !! The significance level used for obtained values in bootstrapping, recommended: `2.5` -> `best_candidate_pair_confidence_interval` caps `95%` from all obtained values
    real(c_double), intent(in), target :: residual_range_quantile
        !! Quantile for determining the residual range, recommended: `95.0`
    integer(c_int), intent(out), target :: ierr
        !! Error code

    integer(int32) :: join_method_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(max_n_reps_all_studies)
    M_CHECK_NON_NULL(max_n_genes_all_studies)
    M_CHECK_NON_NULL(n_studies)
    M_CHECK_NON_NULL(residuals)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(gene_means)
    M_CHECK_NON_NULL(n_bootstraps)
    M_CHECK_NON_NULL(best_candidate_pair_confidence_interval)
    M_CHECK_NON_NULL(join_method)
    M_CHECK_NON_NULL(min_count_per_mean_bin)
    M_CHECK_NON_NULL(min_neighbor_overlap)
    M_CHECK_NON_NULL(succeeding_ci_overlap)
    M_CHECK_NON_NULL(random_seed)
    M_CHECK_NON_NULL(two_sided_bootstrapping_significance_level)
    M_CHECK_NON_NULL(residual_range_quantile)

    call get_join_method(join_method, join_method_f, ierr)

    if (is_err(ierr)) return

    call determine_js_comp_test_n_points_n_neighbors_alloc( &
        n_points, n_neighbors, residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
        shared_residual_range, n_bins, gene_means, n_studies, n_bootstraps, &
        best_candidate_pair_confidence_interval, join_method_f, ierr, &
        min_count_per_mean_bin, min_neighbor_overlap, succeeding_ci_overlap, &
        two_sided_bootstrapping_significance_level, random_seed, residual_range_quantile)
end subroutine determine_js_comp_test_n_points_n_neighbors_c

!> C-compatible wrapper for [[tox_data_integration(module):js_comp_test_alloc(interface)]]
subroutine js_comp_test_c(&
        gene_means, max_n_genes_all_studies, n_studies, residuals, shared_residual_range,&
        n_bins, max_n_reps_all_studies, x_star, n_pool, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
        pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
        js_divergences, weights, global_js_divergence, p_values, ierr, n_permutations, random_seed&
    ) &
        bind(C, name="js_comp_test_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: js_comp_test_alloc
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_studies
        !! Number of studies
    integer(c_int), intent(in), target :: max_n_genes_all_studies
        !! Maximum number of genes across all studies
    integer(c_int), intent(in), target :: max_n_reps_all_studies
        !! Maximum number of replicates across all studies
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points for neighborhoods
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors in neighborhoods
    real(c_double), intent(in), target :: shared_residual_range
        !! Computed residual range (R)
    integer(c_int), intent(in), target :: n_bins
        !! Appropriate number of bins to do the JSD Compatibility test for
    real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
        !! Per-gene mean expression values for all studies
    real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
        !! Matrix of signed residuals per study
    integer(c_int), dimension(2, n_points, n_studies), intent(out), target :: neighborhood_ranges
        !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
        !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
        !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
        !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
        !! If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
    real(c_double), dimension(n_points), intent(out), target :: x_star
        !! Mean-expression reference points
    integer(c_int), intent(out), target :: n_pool
        !! Total number of included (non-NaN) pooled mean-expression values
    integer(c_int), dimension(n_neighbors, n_points, n_studies), intent(out), target :: neighborhood_residuals
        !! Indices of selected neighborhood genes per reference point from [[tox_data_integration(module):construct_neighborhoods(interface)]]
    real(c_double), dimension(n_bins, n_points, n_studies), intent(out), target :: pmfs
        !! `counts` normalized to `0 <= counts(i, :) <= 1` and `sum(counts(i, :)) == 1`
    integer(c_int), dimension(n_bins, n_points, n_studies), intent(out), target :: counts
        !! Absolute counts of a residual per bin for `pmfs`
    integer(c_int), dimension(n_points, n_studies), intent(out), target :: included_n_reps
        !! The count of non-NaN replicates (included ones) per bin and point for `pmfs`
    real(c_double), dimension(n_bins, n_points), intent(out), target :: mean_pmf
        !! The mean pmf built from `pmfs` as `mean_pmf = sum(pmfs) / n_studies`
    integer(c_int), dimension(n_bins, n_points), intent(out), target :: mean_pmf_counts
        !! Absolute counts of a residual per bin for the mean pmf -> `sum(counts)`
    integer(c_int), dimension(n_points), intent(out), target :: mean_pmf_included_n_reps
        !! The count of non-NaN replicates (included ones) per bin and point for `mean_pmf`
    real(c_double), dimension(n_points, n_studies), intent(out), target :: js_divergences
        !! The per-reference-point Jensen-Shannon-Divergence -> finally `js_divergences(1:n_points, :)`
    real(c_double), dimension(n_points, n_studies), intent(out), target :: weights
        !! The per-reference-point weights for the Jensen-Shannon-Divergence -> finally `weights(1:n_points, :)`
    real(c_double), dimension(n_studies), intent(out), target :: global_js_divergence
        !! The global Jensen-Shannon-Divergence
    real(c_double), dimension(n_studies), intent(out), target :: p_values
        !! The p-values from permutation test
    integer(c_int), intent(in), target :: n_permutations
        !! Number of permutations to perform in the permutation test ([[tox_data_integration(module):gjct_permutation_test(interface)]]), recommended: `1000`
    integer(c_int), intent(in), target :: random_seed
        !! Random seed to use for random number generation
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(gene_means)
    M_CHECK_NON_NULL(max_n_genes_all_studies)
    M_CHECK_NON_NULL(n_studies)
    M_CHECK_NON_NULL(residuals)
    M_CHECK_NON_NULL(shared_residual_range)
    M_CHECK_NON_NULL(n_bins)
    M_CHECK_NON_NULL(max_n_reps_all_studies)
    M_CHECK_NON_NULL(x_star)
    M_CHECK_NON_NULL(n_pool)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_neighbors)
    M_CHECK_NON_NULL(neighborhood_ranges)
    M_CHECK_NON_NULL(neighborhood_residuals)
    M_CHECK_NON_NULL(pmfs)
    M_CHECK_NON_NULL(counts)
    M_CHECK_NON_NULL(included_n_reps)
    M_CHECK_NON_NULL(mean_pmf)
    M_CHECK_NON_NULL(mean_pmf_counts)
    M_CHECK_NON_NULL(mean_pmf_included_n_reps)
    M_CHECK_NON_NULL(js_divergences)
    M_CHECK_NON_NULL(weights)
    M_CHECK_NON_NULL(global_js_divergence)
    M_CHECK_NON_NULL(p_values)

    call js_comp_test_alloc(&
        gene_means, max_n_genes_all_studies, n_studies, residuals, shared_residual_range,&
        n_bins, max_n_reps_all_studies, x_star, n_pool, n_points, n_neighbors, neighborhood_ranges, neighborhood_residuals,&
        pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, mean_pmf_included_n_reps,&
        js_divergences, weights, global_js_divergence, p_values, ierr, n_permutations, random_seed&
    )
end subroutine js_comp_test_c
