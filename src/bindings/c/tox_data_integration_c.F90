#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration(module)]]
!| In multi-study omics analyses, it is often unclear whether biological replicates originating from different studies can be safely treated as sampling the same biological condition.
!| Even when studies nominally target the same tissue and condition, differences in sample handling, sequencing technologies, preprocessing pipelines, or cohort,
!| composition can introduce batch effects that are not easily detectable from mean expression levels alone.
!|
!| This ambiguity has direct consequences for downstream analyses in Tensor Omics. Integrating incompatible replicate sets can:
!|
!| - distort expression spaces,
!| - affect distance-based analyses,
!| - bias machine learning models,
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
module tox_data_integration_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: compute_gene_means_c
    public :: compute_residuals_c
    public :: pool_means_c
    public :: pool_means_expert_c
    public :: calc_neighborhood_size_c
    public :: construct_neighborhoods_c
    public :: gjct_permutation_test_c
    public :: gjct_permutation_test_expert_c
    public :: determine_shared_residual_range_expert_c
    public :: determine_shared_residual_range_c
    public :: build_residual_histograms_c
    public :: compute_divergence_per_reference_point_c
    public :: compute_weighted_global_divergence_c
    public :: fjct_compute_jsd_c
    public :: fjct_compute_jsd_expert_c
    public :: fjct_compute_contribution_scores_c

contains

    !> summary: C-wrapper for [[tox_data_integration(module):compute_gene_means(subroutine)]]
    subroutine compute_gene_means_c(&
            n_genes,&
            n_reps,&
            expr,&
            means,&
            ierr&
        ) bind(C, name="compute_gene_means_c")
        use tox_data_integration, only: compute_gene_means

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes in the study
        integer(c_int), intent(in), target :: n_reps
            !! Number of biological replicates in the study
        real(c_double), dimension(n_reps, n_genes), intent(in), target :: expr
            !! Expression matrix
        real(c_double), dimension(n_genes), intent(out), target :: means
            !! Per-gene mean expression values
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_reps)
        M_CHECK_ARRAY_NON_NULL(expr, n_reps * n_genes)
        M_CHECK_ARRAY_NON_NULL(means, n_genes)

        call compute_gene_means(&
            n_genes = n_genes,&
            n_reps = n_reps,&
            expr = expr,&
            means = means,&
            ierr = ierr&
        )
    end subroutine compute_gene_means_c

    !> summary: C-wrapper for [[tox_data_integration(module):compute_residuals(subroutine)]]
    subroutine compute_residuals_c(&
            n_genes,&
            n_reps,&
            expr,&
            means,&
            resid,&
            ierr&
        ) bind(C, name="compute_residuals_c")
        use tox_data_integration, only: compute_residuals

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes in the study
        integer(c_int), intent(in), target :: n_reps
            !! Number of biological replicates in the study
        real(c_double), dimension(n_reps, n_genes), intent(in), target :: expr
            !! Expression matrix containing
        real(c_double), dimension(n_genes), intent(in), target :: means
            !! Per-gene mean expression values
        real(c_double), dimension(n_reps, n_genes), intent(out), target :: resid
            !! Matrix of signed residuals
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(n_reps)
        M_CHECK_ARRAY_NON_NULL(expr, n_reps * n_genes)
        M_CHECK_ARRAY_NON_NULL(means, n_genes)
        M_CHECK_ARRAY_NON_NULL(resid, n_reps * n_genes)

        call compute_residuals(&
            n_genes = n_genes,&
            n_reps = n_reps,&
            expr = expr,&
            means = means,&
            resid = resid,&
            ierr = ierr&
        )
    end subroutine compute_residuals_c

    !> summary: C-wrapper for [[tox_data_integration(module):pool_means_alloc(subroutine)]]
    subroutine pool_means_c(&
            n_genes_S1,&
            mean_S1,&
            n_genes_S2,&
            mean_S2,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_means_c")
        use tox_data_integration, only: pool_means_alloc

        integer(c_int), intent(in), target :: n_genes_S1
            !! Number of genes in study S1
        integer(c_int), intent(in), target :: n_genes_S2
            !! Number of genes in study S2
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(n_genes_S1), intent(in), target :: mean_S1
            !! Per-gene mean expression values
        real(c_double), dimension(n_genes_S2), intent(in), target :: mean_S2
            !! Per-gene mean expression values
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes_S1)
        M_CHECK_NON_NULL(n_genes_S2)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_ARRAY_NON_NULL(mean_S1, n_genes_S1)
        M_CHECK_ARRAY_NON_NULL(mean_S2, n_genes_S2)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)

        call pool_means_alloc(&
            n_genes_S1 = n_genes_S1,&
            mean_S1 = mean_S1,&
            n_genes_S2 = n_genes_S2,&
            mean_S2 = mean_S2,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_means_c

    !> summary: C-wrapper for [[tox_data_integration(module):pool_means(subroutine)]]
    subroutine pool_means_expert_c(&
            pooled_means,&
            pooled_means_perm,&
            pool_size,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_means_expert_c")
        use tox_data_integration, only: pool_means

        integer(c_int), intent(in), target :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(pool_size), intent(in), target :: pooled_means
            !! Pooled means
        integer(c_int), dimension(pool_size), intent(in), target :: pooled_means_perm
            !! Sorting permutation for `pooled_means`
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(pool_size)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_ARRAY_NON_NULL(pooled_means, pool_size)
        M_CHECK_ARRAY_NON_NULL(pooled_means_perm, pool_size)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)

        call pool_means(&
            pooled_means = pooled_means,&
            pooled_means_perm = pooled_means_perm,&
            pool_size = pool_size,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_means_expert_c

    !> summary: C-wrapper for [[tox_data_integration(module):calc_neighborhood_size(function)]]
    !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower due to few genes with non-NaN mean.
    subroutine calc_neighborhood_size_c(&
            n_pool,&
            n_points,&
            n_genes_S,&
            mean_S,&
            desired_size,&
            n_neighbors,&
            ierr&
        ) bind(C, name="calc_neighborhood_size_c")
        use tox_data_integration, only: calc_neighborhood_size

        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_pool
            !! Total number of pooled mean-expression values across both studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
        integer(c_int), intent(in), target :: desired_size
            !! Optional desired neighborhood size.
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

    !> summary: C-wrapper for [[tox_data_integration(module):construct_neighborhoods_alloc(subroutine)]]
    subroutine construct_neighborhoods_c(&
            n_points,&
            x_star,&
            n_genes_S,&
            mean_S,&
            n_reps_S,&
            resid_S,&
            neighborhood_residuals,&
            neighborhood_indices,&
            n_neighbors,&
            ierr&
        ) bind(C, name="construct_neighborhoods_c")
        use tox_data_integration, only: construct_neighborhoods_alloc

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_reps_S
            !! Number of biological replicates in the study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
        real(c_double), dimension(n_reps_S, n_genes_S), intent(in), target :: resid_S
            !! Matrix of signed residuals
        real(c_double), dimension(n_reps_S, n_neighbors, n_points), intent(out), target :: neighborhood_residuals
            !! Collection of residual vectors for each neighborhood
        integer(c_int), dimension(n_neighbors, n_points), intent(out), target :: neighborhood_indices
            !! Indices of selected neighborhood genes
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_genes_S)
        M_CHECK_NON_NULL(n_reps_S)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_S, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(resid_S, n_reps_S * n_genes_S)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals, n_reps_S * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points)

        call construct_neighborhoods_alloc(&
            n_points = n_points,&
            x_star = x_star,&
            n_genes_S = n_genes_S,&
            mean_S = mean_S,&
            n_reps_S = n_reps_S,&
            resid_S = resid_S,&
            neighborhood_residuals = neighborhood_residuals,&
            neighborhood_indices = neighborhood_indices,&
            n_neighbors = n_neighbors,&
            ierr = ierr&
        )
    end subroutine construct_neighborhoods_c

    !> summary: C-wrapper for [[tox_data_integration(module):gjct_permutation_test_alloc(subroutine)]]
    subroutine gjct_permutation_test_c(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            global_jsd_observed,&
            n_bins,&
            shared_residual_range,&
            n_permutations,&
            jsd_null,&
            p_value,&
            ierr,&
            random_seed,&
            neighbor_mask_S1,&
            neighbor_mask_S2&
        ) bind(C, name="gjct_permutation_test_c")
        use tox_data_integration, only: gjct_permutation_test_alloc

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations to perform
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), intent(in), target :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(c_double), dimension(n_permutations), intent(out), target :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(c_double), intent(out), target :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        logical, dimension(:, :), allocatable :: neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(global_jsd_observed)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(p_value)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(jsd_null, n_permutations)

        if (present(neighbor_mask_S1)) neighbor_mask_S1_f = neighbor_mask_S1
        if (present(neighbor_mask_S2)) neighbor_mask_S2_f = neighbor_mask_S2

        call gjct_permutation_test_alloc(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            global_jsd_observed = global_jsd_observed,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            n_permutations = n_permutations,&
            jsd_null = jsd_null,&
            p_value = p_value,&
            ierr = ierr,&
            random_seed = random_seed,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f&
        )
    end subroutine gjct_permutation_test_c

    !> summary: C-wrapper for [[tox_data_integration(module):gjct_permutation_test(subroutine)]]
    subroutine gjct_permutation_test_expert_c(&
            neighborhood_residuals_S1_copy,&
            neighborhood_residuals_S2_copy,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            global_jsd_observed,&
            n_bins,&
            shared_residual_range,&
            n_permutations,&
            jsd_null,&
            p_value,&
            tmp_pool,&
            tmp_pmf_S1,&
            tmp_pmf_S2,&
            tmp_counts,&
            tmp_included_n_reps_S1,&
            tmp_included_n_reps_S2,&
            tmp_js_divergences,&
            tmp_weights,&
            ierr,&
            random_seed,&
            neighbor_mask_S1,&
            neighbor_mask_S2&
        ) bind(C, name="gjct_permutation_test_expert_c")
        use tox_data_integration, only: gjct_permutation_test

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(c_int), intent(in), target :: n_permutations
            !! Number of permutations to perform
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S1_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S2_copy
            !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
        real(c_double), intent(in), target :: global_jsd_observed
            !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(c_double), dimension(n_permutations), intent(out), target :: jsd_null
            !! Vector of global divergence values obtained under the null hypothesis
        real(c_double), intent(out), target :: p_value
            !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)
        real(c_double), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out), target :: tmp_pool
            !! Working array for shuffling the concatenated residuals from both studies per reference point
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_S1
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(c_double), dimension(n_points, n_bins), intent(out), target :: tmp_pmf_S2
            !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(c_int), dimension(n_points), intent(out), target :: tmp_included_n_reps_S1
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(c_int), dimension(n_points), intent(out), target :: tmp_included_n_reps_S2
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(c_double), dimension(n_points), intent(out), target :: tmp_js_divergences
            !! Working array for [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
        real(c_double), dimension(n_points), intent(out), target :: tmp_weights
            !! Working array for [[tox_data_integration(module):compute_weighted_global_divergence(interface)]]
        integer(c_int), intent(out), target :: ierr
            !! Error code
        integer(c_int), intent(in), optional :: random_seed
            !! Seed to use for shuffling
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        logical, dimension(:, :), allocatable :: neighbor_mask_S1_f
        logical, dimension(:, :), allocatable :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(global_jsd_observed)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_permutations)
        M_CHECK_NON_NULL(p_value)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1_copy, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2_copy, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(jsd_null, n_permutations)
        M_CHECK_ARRAY_NON_NULL(tmp_pool, (n_reps_S1 + n_reps_S2) * n_neighbors)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_weights, n_points)

        if (present(neighbor_mask_S1)) neighbor_mask_S1_f = neighbor_mask_S1
        if (present(neighbor_mask_S2)) neighbor_mask_S2_f = neighbor_mask_S2

        call gjct_permutation_test(&
            neighborhood_residuals_S1_copy = neighborhood_residuals_S1_copy,&
            neighborhood_residuals_S2_copy = neighborhood_residuals_S2_copy,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            global_jsd_observed = global_jsd_observed,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            n_permutations = n_permutations,&
            jsd_null = jsd_null,&
            p_value = p_value,&
            tmp_pool = tmp_pool,&
            tmp_pmf_S1 = tmp_pmf_S1,&
            tmp_pmf_S2 = tmp_pmf_S2,&
            tmp_counts = tmp_counts,&
            tmp_included_n_reps_S1 = tmp_included_n_reps_S1,&
            tmp_included_n_reps_S2 = tmp_included_n_reps_S2,&
            tmp_js_divergences = tmp_js_divergences,&
            tmp_weights = tmp_weights,&
            ierr = ierr,&
            random_seed = random_seed,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f&
        )
    end subroutine gjct_permutation_test_expert_c

    !> summary: C-wrapper for [[tox_data_integration(module):determine_shared_residual_range(subroutine)]]
    subroutine determine_shared_residual_range_expert_c(&
            abs_residual_pool,&
            abs_residual_pool_perm,&
            pool_size,&
            shared_residual_range,&
            ierr,&
            residual_range_quantile&
        ) bind(C, name="determine_shared_residual_range_expert_c")
        use tox_data_integration, only: determine_shared_residual_range

        integer(c_int), intent(in), target :: pool_size
            !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
        real(c_double), dimension(pool_size), intent(in), target :: abs_residual_pool
            !! The absolute residual values of the concatenated S1,S2 residuals
        integer(c_int), dimension(pool_size), intent(in), target :: abs_residual_pool_perm
            !! The permutation vector that sorts `abs_residual_pool`
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile for determining the residual range.
            !! The default value is `95.0_real64`.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(pool_size)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(abs_residual_pool, pool_size)
        M_CHECK_ARRAY_NON_NULL(abs_residual_pool_perm, pool_size)

        call determine_shared_residual_range(&
            abs_residual_pool = abs_residual_pool,&
            abs_residual_pool_perm = abs_residual_pool_perm,&
            pool_size = pool_size,&
            shared_residual_range = shared_residual_range,&
            ierr = ierr,&
            residual_range_quantile = residual_range_quantile&
        )
    end subroutine determine_shared_residual_range_expert_c

    !> summary: C-wrapper for [[tox_data_integration(module):determine_shared_residual_range_alloc(subroutine)]]
    subroutine determine_shared_residual_range_c(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            ierr,&
            residual_range_quantile&
        ) bind(C, name="determine_shared_residual_range_c")
        use tox_data_integration, only: determine_shared_residual_range_alloc

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), intent(out), target :: shared_residual_range
            !! Computed residual range (R)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        real(c_double), intent(in), target :: residual_range_quantile
            !! Quantile for determining the residual range.
            !! The default value is `95.0_real64`.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(residual_range_quantile)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)

        call determine_shared_residual_range_alloc(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range = shared_residual_range,&
            ierr = ierr,&
            residual_range_quantile = residual_range_quantile&
        )
    end subroutine determine_shared_residual_range_c

    !> summary: C-wrapper for [[tox_data_integration(module):build_residual_histograms(subroutine)]]
    subroutine build_residual_histograms_c(&
            neighborhood_residuals,&
            n_reps,&
            n_neighbors,&
            n_points,&
            shared_residual_range,&
            n_bins,&
            counts,&
            pmf,&
            included_n_reps,&
            ierr,&
            neighbor_mask&
        ) bind(C, name="build_residual_histograms_c")
        use tox_data_integration, only: build_residual_histograms

        integer(c_int), intent(in), target :: n_reps
            !! Number of replicates of the study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of reference points (k)
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(c_double), dimension(n_reps, n_neighbors, n_points), intent(in), target :: neighborhood_residuals
            !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: counts
            !! Absolute counts of a residual per bin
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf
            !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps
            !! Stores the count of non-NaN replicates (included ones)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
            !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
        logical, dimension(:, :), allocatable :: neighbor_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals, n_reps * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(counts, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(included_n_reps, n_points)

        if (present(neighbor_mask)) neighbor_mask_f = neighbor_mask

        call build_residual_histograms(&
            neighborhood_residuals = neighborhood_residuals,&
            n_reps = n_reps,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            shared_residual_range = shared_residual_range,&
            n_bins = n_bins,&
            counts = counts,&
            pmf = pmf,&
            included_n_reps = included_n_reps,&
            ierr = ierr,&
            neighbor_mask = neighbor_mask_f&
        )
    end subroutine build_residual_histograms_c

    !> summary: C-wrapper for [[tox_data_integration(module):compute_divergence_per_reference_point(subroutine)]]
    subroutine compute_divergence_per_reference_point_c(&
            pmf_S1,&
            pmf_S2,&
            n_points,&
            n_bins,&
            js_divergences,&
            ierr&
        ) bind(C, name="compute_divergence_per_reference_point_c")
        use tox_data_integration, only: compute_divergence_per_reference_point

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points (k)
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins in range [-R,R]
        real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S1
            !! Computed normalized histogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
        real(c_double), dimension(n_points, n_bins), intent(in), target :: pmf_S2
            !! Computed normalized histogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)

        call compute_divergence_per_reference_point(&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            n_points = n_points,&
            n_bins = n_bins,&
            js_divergences = js_divergences,&
            ierr = ierr&
        )
    end subroutine compute_divergence_per_reference_point_c

    !> summary: C-wrapper for [[tox_data_integration(module):compute_weighted_global_divergence(subroutine)]]
    subroutine compute_weighted_global_divergence_c(&
            js_divergences,&
            n_points,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            global_js_divergence,&
            weights,&
            ierr&
        ) bind(C, name="compute_weighted_global_divergence_c")
        use tox_data_integration, only: compute_weighted_global_divergence

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points (k)
        real(c_double), dimension(n_points), intent(in), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(c_int), dimension(n_points), intent(in), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)

        call compute_weighted_global_divergence(&
            js_divergences = js_divergences,&
            n_points = n_points,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            ierr = ierr&
        )
    end subroutine compute_weighted_global_divergence_c

    !> summary: C-wrapper for [[tox_data_integration(module):fjct_compute_jsd_alloc(subroutine)]]
    !| Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
    subroutine fjct_compute_jsd_c(&
            family_idx,&
            gene_to_family_S1,&
            gene_to_family_S2,&
            n_genes_S1,&
            n_genes_S2,&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            neighborhood_genes_S1,&
            neighborhood_genes_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            total_included_n_reps,&
            global_js_divergence,&
            weights,&
            ierr&
        ) bind(C, name="fjct_compute_jsd_c")
        use tox_data_integration, only: fjct_compute_jsd_alloc

        integer(c_int), intent(in), target :: n_genes_S1
            !! Number of genes in study 1
        integer(c_int), intent(in), target :: n_genes_S2
            !! Number of genes in study 2
        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: family_idx
            !! Index of the family that should be analyzed
        integer(c_int), dimension(n_genes_S1), intent(in), target :: gene_to_family_S1
            !! Mapping for study 1: Each index (gene) holds the index of its family
        integer(c_int), dimension(n_genes_S2), intent(in), target :: gene_to_family_S2
            !! Mapping for study 2: Each index (gene) holds the index of its family
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        integer(c_int), dimension(n_neighbors, n_points), intent(in), target :: neighborhood_genes_S1
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(c_int), dimension(n_neighbors, n_points), intent(in), target :: neighborhood_genes_S2
            !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(c_int), intent(out), target :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(family_idx)
        M_CHECK_NON_NULL(n_genes_S1)
        M_CHECK_NON_NULL(n_genes_S2)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(total_included_n_reps)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(gene_to_family_S1, n_genes_S1)
        M_CHECK_ARRAY_NON_NULL(gene_to_family_S2, n_genes_S2)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_genes_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_genes_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)

        call fjct_compute_jsd_alloc(&
            family_idx = family_idx,&
            gene_to_family_S1 = gene_to_family_S1,&
            gene_to_family_S2 = gene_to_family_S2,&
            n_genes_S1 = n_genes_S1,&
            n_genes_S2 = n_genes_S2,&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            neighborhood_genes_S1 = neighborhood_genes_S1,&
            neighborhood_genes_S2 = neighborhood_genes_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            total_included_n_reps = total_included_n_reps,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            ierr = ierr&
        )
    end subroutine fjct_compute_jsd_c

    !> summary: C-wrapper for [[tox_data_integration(module):fjct_compute_jsd(subroutine)]]
    !| Restricts residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2` (typically all neighbors belonging to one gene family; see [[tox_data_integration(module):fjct_compute_jsd_alloc(interface)]] for the family-index-based entry point that builds these masks)
    subroutine fjct_compute_jsd_expert_c(&
            neighborhood_residuals_S1,&
            neighborhood_residuals_S2,&
            n_reps_S1,&
            n_reps_S2,&
            n_neighbors,&
            n_points,&
            neighbor_mask_S1,&
            neighbor_mask_S2,&
            n_bins,&
            shared_residual_range,&
            js_divergences,&
            included_n_reps_S1,&
            included_n_reps_S2,&
            total_included_n_reps,&
            global_js_divergence,&
            weights,&
            pmf_S1,&
            pmf_S2,&
            tmp_counts,&
            ierr&
        ) bind(C, name="fjct_compute_jsd_expert_c")
        use tox_data_integration, only: fjct_compute_jsd

        integer(c_int), intent(in), target :: n_reps_S1
            !! Number of replicates in study 1
        integer(c_int), intent(in), target :: n_reps_S2
            !! Number of replicates in study 2
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors in the studies
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points in the studies
        integer(c_int), intent(in), target :: n_bins
            !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
        real(c_double), dimension(n_reps_S1, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S1
            !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        real(c_double), dimension(n_reps_S2, n_neighbors, n_points), intent(in), target :: neighborhood_residuals_S2
            !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), target :: neighbor_mask_S1
            !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
        logical(c_bool), dimension(n_neighbors, n_points), intent(in), target :: neighbor_mask_S2
            !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        real(c_double), intent(in), target :: shared_residual_range
            !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
        real(c_double), dimension(n_points), intent(out), target :: js_divergences
            !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S1
            !! Count of non-NaN residuals (included ones) in study 1 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(c_int), dimension(n_points), intent(out), target :: included_n_reps_S2
            !! Count of non-NaN residuals (included ones) in study 2 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(c_int), intent(out), target :: total_included_n_reps
            !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        real(c_double), intent(out), target :: global_js_divergence
            !! Weighted global Jensen-Shannon divergence
        real(c_double), dimension(n_points), intent(out), target :: weights
            !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf_S1
            !! Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        real(c_double), dimension(n_points, n_bins), intent(out), target :: pmf_S2
            !! Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
        integer(c_int), dimension(n_points, n_bins), intent(out), target :: tmp_counts
            !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, dimension(n_neighbors, n_points) :: neighbor_mask_S1_f
        logical, dimension(n_neighbors, n_points) :: neighbor_mask_S2_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_reps_S1)
        M_CHECK_NON_NULL(n_reps_S2)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_bins)
        M_CHECK_NON_NULL(shared_residual_range)
        M_CHECK_NON_NULL(total_included_n_reps)
        M_CHECK_NON_NULL(global_js_divergence)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S1, n_reps_S1 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals_S2, n_reps_S2 * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbor_mask_S1, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbor_mask_S2, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(js_divergences, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S1, n_points)
        M_CHECK_ARRAY_NON_NULL(included_n_reps_S2, n_points)
        M_CHECK_ARRAY_NON_NULL(weights, n_points)
        M_CHECK_ARRAY_NON_NULL(pmf_S1, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(pmf_S2, n_points * n_bins)
        M_CHECK_ARRAY_NON_NULL(tmp_counts, n_points * n_bins)

        neighbor_mask_S1_f = neighbor_mask_S1
        neighbor_mask_S2_f = neighbor_mask_S2

        call fjct_compute_jsd(&
            neighborhood_residuals_S1 = neighborhood_residuals_S1,&
            neighborhood_residuals_S2 = neighborhood_residuals_S2,&
            n_reps_S1 = n_reps_S1,&
            n_reps_S2 = n_reps_S2,&
            n_neighbors = n_neighbors,&
            n_points = n_points,&
            neighbor_mask_S1 = neighbor_mask_S1_f,&
            neighbor_mask_S2 = neighbor_mask_S2_f,&
            n_bins = n_bins,&
            shared_residual_range = shared_residual_range,&
            js_divergences = js_divergences,&
            included_n_reps_S1 = included_n_reps_S1,&
            included_n_reps_S2 = included_n_reps_S2,&
            total_included_n_reps = total_included_n_reps,&
            global_js_divergence = global_js_divergence,&
            weights = weights,&
            pmf_S1 = pmf_S1,&
            pmf_S2 = pmf_S2,&
            tmp_counts = tmp_counts,&
            ierr = ierr&
        )
    end subroutine fjct_compute_jsd_expert_c

    !> summary: C-wrapper for [[tox_data_integration(module):fjct_compute_contribution_scores(subroutine)]]
    !| Combines (1) how divergent the family is between the studies, and (2) how much residual support the family has overall, using the outputs from [[tox_data_integration(module):fjct_compute_jsd(interface)]], collected for the analyzed sub-neighborhoods.
    subroutine fjct_compute_contribution_scores_c(&
            global_js_divergences,&
            total_included_n_reps_per_f,&
            k_families,&
            support_weights,&
            contribution_scores,&
            ierr&
        ) bind(C, name="fjct_compute_contribution_scores_c")
        use tox_data_integration, only: fjct_compute_contribution_scores

        integer(c_int), intent(in), target :: k_families
            !! Number of sub-neighborhoods analyzed
        real(c_double), dimension(k_families), intent(in), target :: global_js_divergences
            !! Per-sub-neighborhood weighted global JSD
        integer(c_int), dimension(k_families), intent(in), target :: total_included_n_reps_per_f
            !! Per-sub-neighborhood `total_included_n_reps`
        real(c_double), dimension(k_families), intent(out), target :: support_weights
            !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
        real(c_double), dimension(k_families), intent(out), target :: contribution_scores
            !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(k_families)
        M_CHECK_ARRAY_NON_NULL(global_js_divergences, k_families)
        M_CHECK_ARRAY_NON_NULL(total_included_n_reps_per_f, k_families)
        M_CHECK_ARRAY_NON_NULL(support_weights, k_families)
        M_CHECK_ARRAY_NON_NULL(contribution_scores, k_families)

        call fjct_compute_contribution_scores(&
            global_js_divergences = global_js_divergences,&
            total_included_n_reps_per_f = total_included_n_reps_per_f,&
            k_families = k_families,&
            support_weights = support_weights,&
            contribution_scores = contribution_scores,&
            ierr = ierr&
        )
    end subroutine fjct_compute_contribution_scores_c

end module tox_data_integration_c
#endif
