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
    implicit none


    interface compute_gene_means
        !> Compute per-gene mean expression, ignoring NaN values
        !|
        !| @note
        !| The means of all studies should be in contiguous memory afterwards, so for using this subroutine pass `means` as `means(:, study_idx)`
        !| @endnote
        pure module subroutine compute_gene_means(expr, n_genes, n_reps, means, max_n_genes_all_studies, ierr)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: max_n_genes_all_studies
                !! Maximum number of genes across all studies
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix
            real(real64), intent(out) :: means(max_n_genes_all_studies)
                !! Per-gene mean expression values
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_gene_means
    end interface compute_gene_means

    interface compute_gene_means_helper
        !> (no input validation) Compute per-gene mean expression, ignoring NaN values
        !|
        !| @note
        !| The means of all studies should be in contiguous memory afterwards, so for using this subroutine pass `means` as `means(:, study_idx)`
        !| @endnote
        pure module subroutine compute_gene_means_helper(expr, n_genes, n_reps, means, max_n_genes_all_studies)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: max_n_genes_all_studies
                !! Maximum number of genes across all studies
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix
            real(real64), intent(out) :: means(max_n_genes_all_studies)
                !! Per-gene mean expression values
        end subroutine compute_gene_means_helper
    end interface compute_gene_means_helper

    interface compute_residuals
        !> Compute signed residuals (centering by mean).
        !|
        !| @note
        !| The residuals of all studies should be in contiguous memory afterwards, so for using this subroutine pass `expr` as `expr(:, :, study_idx)`
        !| @endnote
        pure module subroutine compute_residuals(expr, n_genes, n_reps, means, max_n_genes_all_studies, max_n_reps_all_studies, resid, ierr)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: max_n_genes_all_studies
                !! Maximum number of genes across all studies
            integer(int32), intent(in) :: max_n_reps_all_studies
                !! Maximum number of replicates across all studies
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix containing
            real(real64), intent(in) :: means(max_n_genes_all_studies)
                !! Per-gene mean expression values
            real(real64), intent(out) :: resid(max_n_reps_all_studies, max_n_genes_all_studies)
                !! Matrix of signed residuals
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine compute_residuals
    end interface compute_residuals

    interface compute_residuals_helper
        !> (no input validation) Compute signed residuals (centering by mean)
        !|
        !| @note
        !| The residuals of all studies should be in contiguous memory afterwards, so for using this subroutine pass `expr` as `expr(:, :, study_idx)`
        !| @endnote
        pure module subroutine compute_residuals_helper(expr, n_genes, n_reps, means, max_n_genes_all_studies, max_n_reps_all_studies, resid)
            integer(int32), intent(in) :: n_genes
                !! Number of genes in the study
            integer(int32), intent(in) :: n_reps
                !! Number of biological replicates in the study
            integer(int32), intent(in) :: max_n_genes_all_studies
                !! Maximum number of genes across all studies
            integer(int32), intent(in) :: max_n_reps_all_studies
                !! Maximum number of replicates across all studies
            real(real64), intent(in) :: expr(n_reps, n_genes)
                !! Expression matrix containing
            real(real64), intent(in) :: means(n_genes)
                !! Per-gene mean expression values
            real(real64), intent(out) :: resid(max_n_reps_all_studies, max_n_genes_all_studies)
                !! Matrix of signed residuals
        end subroutine compute_residuals_helper
    end interface compute_residuals_helper

    interface find_last_non_nan
        !> Returns the last index in a sorted array that has a value unlike NaN (NaN sorted to the end)
        pure module function find_last_non_nan(arr, arr_perm, n_elements) result(idx)
            integer(int32), intent(in) :: n_elements
                !! Number of elements in `arr`
            real(real64), dimension(n_elements), intent(in) :: arr
                !! Array to find the last non-NaN index
            integer(int32), dimension(n_elements), intent(in) :: arr_perm
                !! Sorting permutation for `arr`
            integer(int32) :: idx
                !! Index in `arr` with non-NaN value
        end function find_last_non_nan
    end interface find_last_non_nan

    interface pool_means_alloc
        !> Pool per-gene mean expression values across studies
        pure module subroutine pool_means_alloc(means, n_studies, max_n_genes_all_studies, n_points, n_pool, x_star, ierr)
            integer(int32), intent(in) :: n_studies
                !! Number of studies
            integer(int32), intent(in) :: max_n_genes_all_studies
                !! Maximum number of genes across all studies
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            real(real64), intent(in) :: means(max_n_genes_all_studies * n_studies)
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
            integer(int32), intent(in) :: pool_size
                !! Number of means in the pool, usually `2 * max_n_genes_all_studies`
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
            integer(int32), intent(in) :: pool_size
                !! Number of means in the pool, usually `2 * max_n_genes_all_studies`
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

    interface pool_means_n_pool_input_helper
        !> (no input validation) Pool per-gene mean expression values across studies
        pure module subroutine pool_means_n_pool_input_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
            integer(int32), intent(in) :: pool_size
                !! Number of means in the pool, usually `2 * max_n_genes_all_studies`
            integer(int32), intent(in) :: n_points
                !! Number of reference points to define
            integer(int32), intent(in) :: n_pool
                !! Total number of included (non-NaN) pooled mean-expression values
            real(real64), intent(in) :: pooled_means(pool_size)
                !! Pooled means
            integer(int32), intent(in) :: pooled_means_perm(pool_size)
                !! Sorting permutation for `pooled_means`
            real(real64), intent(out) :: x_star(n_points)
                !! Mean-expression reference points
        end subroutine pool_means_n_pool_input_helper
    end interface pool_means_n_pool_input_helper

    interface calc_neighborhood_size
        !> Calculate the number of neighbors to be used for [[tox_data_integration(module):construct_neighborhoods(interface)]].
        !|
        !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower due to few genes with non-NaN mean.
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
        pure module subroutine construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Array of per-gene mean values
            integer(int32), intent(out) :: neighborhood_residuals(n_neighbors, n_points)
                !! Indices of selected neighborhood genes per reference point.
                !!
                !! @note 
                !! All indices in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case `n_genes_S` is lower than `n_neighbors`,
                !! remaining indices will be filled with the ones from `n_genes_S+1...n_neighbors`
                !! @endnote
            integer(int32), intent(out) :: neighborhood_range(2, n_points)
                !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
                !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
                !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
                !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
                !! If all mean values are NaN, the range is [1, 1]
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine construct_neighborhoods_alloc
    end interface construct_neighborhoods_alloc

    interface construct_neighborhoods
        !> (no input validation) Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Array of per-gene mean values
            integer(int32), intent(in) :: mean_S_perm(n_genes_S)
                !! Sorting permutation for `mean_S`
            integer(int32), intent(out) :: neighborhood_residuals(n_neighbors, n_points)
                !! Indices of selected neighborhood genes per reference point.
                !!
                !! @note 
                !! All indices in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case `n_genes_S` is lower than `n_neighbors`,
                !! remaining indices will be filled with the ones from `n_genes_S+1...n_neighbors`
                !! @endnote
            integer(int32), intent(out) :: neighborhood_range(2, n_points)
                !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
                !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
                !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
                !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
                !! If all mean values are NaN, the range is [1, 1]
            integer(int32), intent(out) :: ierr
                !! Error code
        end subroutine construct_neighborhoods
    end interface construct_neighborhoods

    interface construct_neighborhoods_helper
        !> (no input validation) Construct neighborhood-based residual sets (kNN)
        pure module subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, neighborhood_range, n_neighbors)
            integer(int32), intent(in) :: n_points
                !! Number of reference points
            integer(int32), intent(in) :: n_genes_S
                !! Number of genes in the current study
            integer(int32), intent(in) :: n_neighbors
                !! Number of neighbors, **CALCULATE IT WITH [[tox_data_integration(module):calc_neighborhood_size(interface)]]**
            real(real64), intent(in) :: x_star(n_points)
                !! Mean-expression reference points
            real(real64), intent(in) :: mean_S(n_genes_S)
                !! Array of per-gene mean values
            integer(int32), intent(in) :: mean_S_perm(n_genes_S)
                !! Sorting permutation for `mean_S`
            integer(int32), intent(out) :: neighborhood_residuals(n_neighbors, n_points)
                !! Indices of selected neighborhood genes per reference point.
                !!
                !! @note 
                !! All indices in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case `n_genes_S` is lower than `n_neighbors`,
                !! remaining indices will be filled with the ones from `n_genes_S+1...n_neighbors`
                !! @endnote
            integer(int32), intent(out) :: neighborhood_range(2, n_points)
                !! For each reference point it holds, the [min_idx, max_idx] of the included genes.
                !! The index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))` would be the min value.
                !! In case of duplicate means, `min_idx` points to the first appearance of the value and `max_idx` to the last,
                !! so even though their related mean value is the min/max in the neighborhood, the actual gene might not be included.
                !! If all mean values are NaN, the range is [1, 1]
        end subroutine construct_neighborhoods_helper
    end interface construct_neighborhoods_helper
end module tox_data_integration