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
    use iso_fortran_env, only: int32, real64
    implicit none

    interface compute_gene_means
        !> Compute per-gene mean expression, ignoring NaN values
        pure module subroutine compute_gene_means(n_genes, n_reps, expr, means, ierr)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix
            real(real64), intent(out) :: means(n_genes)
                !! Per-gene mean expression values
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_gene_means
    end interface compute_gene_means

    interface compute_gene_means_helper
        !> (no input validation) Compute per-gene mean expression, ignoring NaN values
        pure module subroutine compute_gene_means_helper(n_genes, n_reps, expr, means)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix
            real(real64), intent(out) :: means(n_genes)
                !! Per-gene mean expression values
        end subroutine compute_gene_means_helper
    end interface compute_gene_means_helper

    interface compute_residuals
        !> Compute signed residuals (centering by mean)
        pure module subroutine compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix containing
            real(real64), intent(in) :: means(n_genes)
                !! Per-gene mean expression values
            real(real64), intent(out) :: resid(n_reps, n_genes)
                !! Matrix of signed residuals
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_residuals
    end interface compute_residuals

    interface compute_residuals_helper
        !> (no input validation) Compute signed residuals (centering by mean)
        pure module subroutine compute_residuals_helper(n_genes, n_reps, expr, means, resid)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix containing
            real(real64), intent(in) :: means(n_genes)
                !! Per-gene mean expression values
            real(real64), intent(out) :: resid(n_reps, n_genes)
                !! Matrix of signed residuals
        end subroutine compute_residuals_helper
    end interface compute_residuals_helper

    interface pool_means_alloc
        !> Pool per-gene mean expression values across studies
        pure module subroutine pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, n_pool, x_star, ierr)
            integer(int32), intent(in) :: n_genes_S1
                !! Number of genes in study S1
            integer(int32), intent(in) :: n_genes_S2
                !! Number of genes in study S2
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            real(real64), intent(in) :: mean_S1(n_genes_S1)
                !! Per-gene mean expression values
            real(real64), intent(in) :: mean_S2(n_genes_S2)
                !! Per-gene mean expression values
            integer(int32), intent(out) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine pool_means_alloc
    end interface pool_means_alloc

    interface pool_means
        !> Pool per-gene mean expression values across studies
        pure module subroutine pool_means(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star, ierr)
            integer(int32), intent(in), target :: pool_size
                !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            integer(int32), intent(out) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(in) :: pooled_means(pool_size)
                !! Pooled means
            integer(int32), intent(in) :: pooled_means_perm(pool_size)
                !! Sorting permutation for `pooled_means`
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine pool_means
    end interface pool_means

    interface pool_means_helper
        !> (no input validation) Pool per-gene mean expression values across studies
        pure module subroutine pool_means_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
            integer(int32), intent(in), target :: pool_size
                !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            integer(int32), intent(out) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(in) :: pooled_means(pool_size)
                !! Pooled means
            integer(int32), intent(in) :: pooled_means_perm(pool_size)
                !! Sorting permutation for `pooled_means`
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
        end subroutine pool_means_helper
    end interface pool_means_helper

        !> Calculate the number of neighbors to be used for [[tox_data_integration(module):construct_neighborhoods(interface)]].
        !|
        !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower due to few genes with non-NaN mean.
    interface calc_neighborhood_size
        pure module function calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S, desired_size) result(n_neighbors)
            integer(int32), intent(in) :: n_pool
                !! Total number of pooled mean-expression values across both studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            integer(int32), intent(in), optional :: desired_size
                !! Optional desired neighborhood size, default=1000
            integer(int32) :: n_neighbors
                !! Calculated neighborhood size
        end function calc_neighborhood_size
    end interface calc_neighborhood_size

    interface construct_neighborhoods_alloc
        !> Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                                      neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_reps_S
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
                !! Matrix of signed residuals
            real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
                !! Collection of residual vectors for each neighborhood
            integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
                !! Indices of selected neighborhood genes
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine construct_neighborhoods_alloc
    end interface construct_neighborhoods_alloc

    interface construct_neighborhoods
        !> Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_reps_S
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
                !! Matrix of signed residuals
            real(real64), intent(out) :: tmp_distances(n_genes_S)
                !! Distances work array
            integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
                !! Work array for permutation vector to sort `tmp_distances_perm`
            real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
                !! Collection of residual vectors for each neighborhood
            integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
                !! Indices of selected neighborhood genes
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine construct_neighborhoods
    end interface construct_neighborhoods

    interface construct_neighborhoods_helper
        !> (no input validation) Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, &
                                                       neighborhood_residuals, neighborhood_indices, n_neighbors)
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_reps_S
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Per-gene mean expression values
            real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
                !! Matrix of signed residuals
            real(real64), intent(out) :: tmp_distances(n_genes_S)
                !! Distances work array
            integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
                !! Work array for permutation vector to sort `tmp_distances_perm`
            real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
                !! Collection of residual vectors for each neighborhood
            integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
                !! Indices of selected neighborhood genes
        end subroutine construct_neighborhoods_helper
    end interface construct_neighborhoods_helper

    interface gjct_permutation_test_alloc
        !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
        module subroutine gjct_permutation_test_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2)
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
                !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
                !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), intent(in) :: global_jsd_observed
                !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
            integer(int32), intent(in) :: n_permutations
                !! Number of permutations to perform
            real(real64), dimension(n_permutations), intent(out) :: jsd_null
                !! Vector of global divergence values obtained under the null hypothesis
            real(real64), intent(out) :: p_value
                !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
            integer(int32), intent(out) :: ierr
                !! Error code
            integer(int32), intent(in), optional :: random_seed
                !! Seed to use for shuffling
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine gjct_permutation_test_alloc
    end interface gjct_permutation_test_alloc

    interface gjct_permutation_test
        !> Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
        module subroutine gjct_permutation_test( &
                neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
                tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
                ierr, random_seed, neighbor_mask_S1, neighbor_mask_S2 &
            )
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S1_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout) :: neighborhood_residuals_S2_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
            real(real64), intent(in) :: global_jsd_observed
                !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
            integer(int32), intent(in) :: n_permutations
                !! Number of permutations to perform
            real(real64), dimension(n_permutations), intent(out) :: jsd_null
                !! Vector of global divergence values obtained under the null hypothesis
            real(real64), intent(out) :: p_value
                !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
            real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out) :: tmp_pool
                !! Working array for shuffling the concatenated residuals from both studies per reference point
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
                !! Working array for [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
            real(real64), dimension(n_points), intent(out) :: tmp_weights
                !! Working array for [[tox_data_integration(module):compute_weighted_global_divergence(interface)]]
            integer(int32), intent(out) :: ierr
                !! Error code
            integer(int32), intent(in), optional :: random_seed
                !! Seed to use for shuffling
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine gjct_permutation_test
    end interface gjct_permutation_test

    interface gjct_permutation_test_helper
        !> (no input validation) Estimates how likely the observed divergence is to occur by chance under the null hypothesis that both studies are exchangeable
        module subroutine gjct_permutation_test_helper( &
                neighborhood_residuals_S1_copy, neighborhood_residuals_S2_copy, n_reps_S1, n_reps_S2, n_neighbors, n_points, global_jsd_observed, n_bins, shared_residual_range, n_permutations, jsd_null, p_value, &
                tmp_pool, tmp_pmf_S1, tmp_pmf_S2, tmp_counts, tmp_included_n_reps_S1, tmp_included_n_reps_S2, tmp_js_divergences, tmp_weights, &
                random_seed, neighbor_mask_S1, neighbor_mask_S2 &
            )
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S1_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 1, will be shuffled in-place
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(inout), target :: neighborhood_residuals_S2_copy
                !! Copy (if wanted) of the computed neighborhood residuals for study 2, will be shuffled in-place
            real(real64), intent(in) :: global_jsd_observed
                !! Observed global JSD value for both studies (from [[tox_data_integration(module):compute_weighted_global_divergence(interface)]])
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
            integer(int32), intent(in) :: n_permutations
                !! Number of permutations to perform
            real(real64), dimension(n_permutations), intent(out) :: jsd_null
                !! Vector of global divergence values obtained under the null hypothesis
            real(real64), intent(out) :: p_value
                !! Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations} \)
            real(real64), dimension(n_reps_S1 + n_reps_S2, n_neighbors), intent(out), target :: tmp_pool
                !! Working array for shuffling the concatenated residuals from both studies per reference point
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S1
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), dimension(n_points, n_bins), intent(out) :: tmp_pmf_S2
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S1
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points), intent(out) :: tmp_included_n_reps_S2
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), dimension(n_points), intent(out) :: tmp_js_divergences
                !! Working array for [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
            real(real64), dimension(n_points), intent(out) :: tmp_weights
                !! Working array for [[tox_data_integration(module):compute_weighted_global_divergence(interface)]]
            integer(int32), intent(in), optional :: random_seed
                !! Seed to use for shuffling
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine gjct_permutation_test_helper
    end interface gjct_permutation_test_helper

    interface shuffle_reference_point_helper
        !> Helper for [[tox_data_integration(module):gjct_permutation_test_helper(interface)]] to shuffle reference points
        module subroutine shuffle_reference_point_helper(reference_point_S1, reference_point_S2, n_reps_S1, n_reps_S2, n_neighbors, pool_flat)
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            real(real64), dimension(n_reps_S1 * n_neighbors), intent(inout) :: reference_point_S1
                !! Residuals for one reference point in study 1, will be shuffled in-place
            real(real64), dimension(n_reps_S2 * n_neighbors), intent(inout) :: reference_point_S2
                !! Residuals for one reference point in study 2, will be shuffled in-place
            real(real64), dimension((n_reps_S1 + n_reps_S2) * n_neighbors), intent(out) :: pool_flat
                !! Working array for shuffling the concatenated residuals from both studies per reference point
        end subroutine shuffle_reference_point_helper
    end interface shuffle_reference_point_helper

    interface determine_shared_residual_range
        !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
        pure module subroutine determine_shared_residual_range(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, ierr, residual_range_quantile)
            integer(int32), intent(in) :: pool_size
                !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
            real(real64), intent(in), optional :: residual_range_quantile
                !! Quantile for determining the residual range, default: 95.0
            real(real64), intent(out) :: shared_residual_range
                !! Computed residual range (R)
            real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
                !! The absolute residual values of the concatenated S1,S2 residuals
            integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
                !! The permutation vector that sorts `abs_residual_pool`
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine determine_shared_residual_range
    end interface determine_shared_residual_range

    interface determine_shared_residual_range_helper
        !> (no input validation) Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
        pure module subroutine determine_shared_residual_range_helper(abs_residual_pool, abs_residual_pool_perm, pool_size, shared_residual_range, residual_range_quantile)
            integer(int32), intent(in) :: pool_size
                !! Size of pool of residuals `abs_residual_pool`, usually `(n_reps_S1 + n_reps_2)*n_neighbors*n_points`
            real(real64), intent(in), optional :: residual_range_quantile
                !! Quantile for determining the residual range, default: 95.0
            real(real64), intent(out) :: shared_residual_range
                !! Computed residual range (R)
            real(real64), dimension(pool_size), intent(in) :: abs_residual_pool
                !! The absolute residual values of the concatenated S1,S2 residuals
            integer(int32), dimension(pool_size), intent(in) :: abs_residual_pool_perm
                !! The permutation vector that sorts `abs_residual_pool`
        end subroutine determine_shared_residual_range_helper
    end interface determine_shared_residual_range_helper

    interface determine_shared_residual_range_alloc
        !> Computes the shared residual range [-R, R] for the computed residuals from studies S1 and S2
        pure module subroutine determine_shared_residual_range_alloc(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, shared_residual_range, ierr, residual_range_quantile)
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
                !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
                !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), intent(in), optional :: residual_range_quantile
                !! Quantile for determining the residual range, default: 95.0
            real(real64), intent(out) :: shared_residual_range
                !! Computed residual range (R)
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine determine_shared_residual_range_alloc
    end interface determine_shared_residual_range_alloc

    interface build_residual_histograms
        !> Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
        pure module subroutine build_residual_histograms(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, ierr, neighbor_mask)
            integer(int32), intent(in) :: n_reps
                !! Number of replicates of the study
            integer(int32), intent(in) :: n_neighbors
                !! Number of reference points (k)
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
                !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins in range [-R,R]
            integer(int32), dimension(n_points, n_bins), intent(out) :: counts
                !! Absolute counts of a residual per bin
            real(real64), dimension(n_points, n_bins), intent(out) :: pmf
                !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
            integer(int32), dimension(n_points), intent(out) :: included_n_reps
                !! Stores the count of non-NaN replicates (included ones)
            integer(int32), intent(out) :: ierr
                !! Error code
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
                !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
        end subroutine build_residual_histograms
    end interface build_residual_histograms

    interface build_residual_histograms_helper
        !> (no input validation) Summarizes the neighborhood residuals in absolute histogram counts and probability mass functions `pmf(residual, bin)` (actually a matrix)
        pure module subroutine build_residual_histograms_helper(neighborhood_residuals, n_reps, n_neighbors, n_points, shared_residual_range, n_bins, counts, pmf, included_n_reps, neighbor_mask)
            integer(int32), intent(in) :: n_reps
                !! Number of replicates of the study
            integer(int32), intent(in) :: n_neighbors
                !! Number of reference points (k)
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the study
            real(real64), dimension(n_reps, n_neighbors, n_points), intent(in) :: neighborhood_residuals
                !! Computed neighborhood residuals for a study ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range (R) from [[tox_data_integration(module):determine_shared_residual_range_alloc(interface)]]
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins in range [-R,R]
            integer(int32), dimension(n_points, n_bins), intent(out) :: counts
                !! Absolute counts of a residual per bin
            real(real64), dimension(n_points, n_bins), intent(out) :: pmf
                !! `counts` normalized to `0 <= counts(:, i) <= 1` and `sum(counts(:, i)) == 1`
            integer(int32), dimension(n_points), intent(out) :: included_n_reps
                !! Stores the count of non-NaN replicates (included ones)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask
                !! Optional mask to exclude specific neighbors (e.g. for family-wise analysis)
        end subroutine build_residual_histograms_helper
    end interface build_residual_histograms_helper

    interface compute_divergence_per_reference_point
        !> Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
        pure module subroutine compute_divergence_per_reference_point(pmf_S1, pmf_S2, n_points, n_bins, js_divergences, ierr)
            integer(int32), intent(in) :: n_points
                !! Number of reference points (k)
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins in range [-R,R]
            real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
                !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
            real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
                !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
            real(real64), dimension(n_points), intent(out) :: js_divergences
                !! Jensen-Shannon divergence per reference point
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_divergence_per_reference_point
    end interface compute_divergence_per_reference_point

    interface compute_divergence_per_reference_point_helper
        !> (no input validation) Having the probabilities `pmf` from [[tox_data_integration(module):build_residual_histograms(interface)]], this subroutine computes the Jensen-Shannon divergence per reference point/neighbor
        pure module subroutine compute_divergence_per_reference_point_helper(pmf_S1, pmf_S2, n_points, n_bins, js_divergences)
            integer(int32), intent(in) :: n_points
                !! Number of reference points (k)
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins in range [-R,R]
            real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S1
                !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 1
            real(real64), dimension(n_points, n_bins), intent(in) :: pmf_S2
                !! Computed normalized hostogram counts from [[tox_data_integration(module):build_residual_histograms(interface)]] for study 2
            real(real64), dimension(n_points), intent(out) :: js_divergences
                !! Jensen-Shannon divergence per reference point
        end subroutine compute_divergence_per_reference_point_helper
    end interface compute_divergence_per_reference_point_helper

    interface compute_weighted_global_divergence
        !> Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
        pure module subroutine compute_weighted_global_divergence(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, ierr)
            integer(int32), intent(in) :: n_points
                !! Number of reference points (k)
            real(real64), dimension(n_points), intent(in) :: js_divergences
                !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
                !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
                !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            real(real64), intent(out) :: global_js_divergence
                !! Weighted global Jensen-Shannon divergence
            real(real64), dimension(n_points), intent(out) :: weights
                !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_weighted_global_divergence
    end interface compute_weighted_global_divergence

    interface compute_weighted_global_divergence_helper
        !> (no input validation) Computes the global weighted Jensen-Shannon divergence from the per-neighbor divergences calculated by [[tox_data_integration(module):compute_divergence_per_reference_point(interface)]]
        pure module subroutine compute_weighted_global_divergence_helper(js_divergences, n_points, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights)
            integer(int32), intent(in) :: n_points
                !! Number of reference points (k)
            real(real64), dimension(n_points), intent(in) :: js_divergences
                !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            integer(int32), dimension(n_points), intent(in) :: included_n_reps_S1
                !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), dimension(n_points), intent(in) :: included_n_reps_S2
                !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            real(real64), intent(out) :: global_js_divergence
                !! Weighted global Jensen-Shannon divergence
            real(real64), dimension(n_points), intent(out) :: weights
                !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
        end subroutine compute_weighted_global_divergence_helper
    end interface compute_weighted_global_divergence_helper

    interface jct_compute_jsd_pipeline_helper
        !> Helper to run the pipeline `build_residual_histograms` \(\Rightarrow\) `compute_weighted_global_divergence`
        pure module subroutine jct_compute_jsd_pipeline_helper(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, neighbor_mask_S1, neighbor_mask_S2)
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
                !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
                !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
            real(real64), dimension(n_points), intent(out) :: js_divergences
                !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
                !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
                !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            real(real64), intent(out) :: global_js_divergence
                !! Weighted global Jensen-Shannon divergence
            real(real64), dimension(n_points), intent(out) :: weights
                !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
                !! Absolute counts of a residual per bin obtained from [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in), optional :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
        end subroutine jct_compute_jsd_pipeline_helper
    end interface jct_compute_jsd_pipeline_helper

    interface fjct_compute_jsd_alloc
        !> Computes the family-level compatibility score `global_js_divergence` between two studies for a single gene family (`family_idx`), by reusing the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
        pure module subroutine fjct_compute_jsd_alloc(family_idx, gene_to_family_S1, gene_to_family_S2, n_genes_S1, n_genes_S2, neighborhood_residuals_S1, neighborhood_residuals_S2, &
                neighborhood_genes_S1, neighborhood_genes_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, n_bins, shared_residual_range, js_divergences, &
                included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, ierr &
            )
            integer(int32), intent(in) :: n_genes_S1
                !! Number of genes in study 1
            integer(int32), intent(in) :: n_genes_S2
                !! Number of genes in study 2
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            integer(int32), intent(in) :: family_idx
                !! Index of the family that should be analyzed
            integer(int32), dimension(n_genes_S1), intent(in) :: gene_to_family_S1
                !! Mapping for study 1: Each index (gene) holds the index of its family
            integer(int32), dimension(n_genes_S2), intent(in) :: gene_to_family_S2
                !! Mapping for study 2: Each index (gene) holds the index of its family
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
                !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
                !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S1
                !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
            integer(int32), dimension(n_neighbors, n_points), intent(in) :: neighborhood_genes_S2
                !! Indices of selected neighborhood genes, obtained from `neighborhood_indices` of [[tox_data_integration(module):construct_neighborhoods(interface)]]
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
            real(real64), dimension(n_points), intent(out) :: js_divergences
                !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
                !! Count of non-NaN residuals (included ones) in study 1 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
                !! Count of non-NaN residuals (included ones) in study 2 (obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), intent(out) :: total_included_n_reps
                !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
            real(real64), intent(out) :: global_js_divergence
                !! Weighted global Jensen-Shannon divergence
            real(real64), dimension(n_points), intent(out) :: weights
                !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine fjct_compute_jsd_alloc
    end interface fjct_compute_jsd_alloc

    interface fjct_compute_jsd
        !> Computes the compatibility score `global_js_divergence` between two studies per sub-neighborhood/family for a single gene family (`family_idx`), by reusing the same conditioning-on-mean-expression pipeline as the global gJCT, but restricting residual samples to genes belonging to the specified family
        pure module subroutine fjct_compute_jsd(neighborhood_residuals_S1, neighborhood_residuals_S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, neighbor_mask_S1, neighbor_mask_S2, n_bins, shared_residual_range, js_divergences, included_n_reps_S1, included_n_reps_S2, total_included_n_reps, global_js_divergence, weights, pmf_S1, pmf_S2, tmp_counts, ierr)
            integer(int32), intent(in) :: n_reps_S1
                !! Number of replicates in study 1
            integer(int32), intent(in) :: n_reps_S2
                !! Number of replicates in study 2
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors in the studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points in the studies
            real(real64), dimension(n_reps_S1, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S1
                !! Computed neighborhood residuals for study 1 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            real(real64), dimension(n_reps_S2, n_neighbors, n_points), intent(in) :: neighborhood_residuals_S2
                !! Computed neighborhood residuals for study 2 ([[tox_data_integration(module):construct_neighborhoods(interface)]]), NaN is explicitly allowed for missing values
            logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S1
                !! Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
            logical, dimension(n_neighbors, n_points), intent(in) :: neighbor_mask_S2
                !! Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)
            integer(int32), intent(in) :: n_bins
                !! Number of equally sized histogram bins used for the studies in [[tox_data_integration(module):build_residual_histograms(interface)]]
            real(real64), intent(in) :: shared_residual_range
                !! Computed residual range for both studies, from [[tox_data_integration(module):determine_shared_residual_range(interface)]]
            real(real64), dimension(n_points), intent(out) :: js_divergences
                !! Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            integer(int32), dimension(n_points), intent(out) :: included_n_reps_S1
                !! Count of non-NaN residuals (included ones) in study 1 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), dimension(n_points), intent(out) :: included_n_reps_S2
                !! Count of non-NaN residuals (included ones) in study 2 (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), intent(out) :: total_included_n_reps
                !! Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
            real(real64), intent(out) :: global_js_divergence
                !! Weighted global Jensen-Shannon divergence
            real(real64), dimension(n_points), intent(out) :: weights
                !! Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S1
                !! Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            real(real64), dimension(n_points, n_bins), intent(out) :: pmf_S2
                !! Absolute counts of a residual per bin (will be obtained from [[tox_data_integration(module):build_residual_histograms(interface)]])
            integer(int32), dimension(n_points, n_bins), intent(out) :: tmp_counts
                !! Working array for [[tox_data_integration(module):build_residual_histograms(interface)]]
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine fjct_compute_jsd
    end interface fjct_compute_jsd

    interface fjct_compute_contribution_scores
        pure module subroutine fjct_compute_contribution_scores(global_js_divergences, total_included_n_reps_per_f, k_families, support_weights, contribution_scores, ierr)
            integer(int32), intent(in) :: k_families
                !! Number of sub-neighborhoods analyzed
            integer(int32), dimension(k_families), intent(in) :: total_included_n_reps_per_f
                !! Per-sub-neighborhood `total_included_n_reps`
            real(real64), dimension(k_families), intent(in) :: global_js_divergences
                !! Per-sub-neighborhood weighted global JSD
            real(real64), dimension(k_families), intent(out) :: support_weights
                !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
            real(real64), dimension(k_families), intent(out) :: contribution_scores
                !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
            integer(int32), intent(out), target :: ierr
                !! Error code
        end subroutine fjct_compute_contribution_scores
    end interface fjct_compute_contribution_scores

    interface fjct_compute_contribution_scores_helper
        !> (no input validation) Computes the per-family/per-sub-neighborhood contribution score that combines
        !|
        !| 1. how divergent the family is between the studies (``), and
        !| 2. how much residual support the family has overall (),
        !|
        !| using the outputs from [[tox_data_integration_per_family(module):fjct_compute_jsd(subroutine)]], collected for the analyzed sub-neighborhoods.
        pure module subroutine fjct_compute_contribution_scores_helper(global_js_divergences, total_included_n_reps_per_f, k_families, support_weights, contribution_scores)
            integer(int32), intent(in) :: k_families
                !! Number of sub-neighborhoods analyzed
            integer(int32), dimension(k_families), intent(in) :: total_included_n_reps_per_f
                !! Per-sub-neighborhood `total_included_n_reps`
            real(real64), dimension(k_families), intent(in) :: global_js_divergences
                !! Per-sub-neighborhood weighted global JSD
            real(real64), dimension(k_families), intent(out) :: support_weights
                !! Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
            real(real64), dimension(k_families), intent(out) :: contribution_scores
                !! Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
        end subroutine fjct_compute_contribution_scores_helper
    end interface fjct_compute_contribution_scores_helper
end module tox_data_integration