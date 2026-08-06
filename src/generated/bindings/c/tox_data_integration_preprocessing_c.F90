#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_preprocessing(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_data_integration_preprocessing_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: compute_gene_means_c
    public :: compute_residuals_c
    public :: pool_means_expert_c
    public :: pool_means_c
    public :: pool_study_means_expert_c
    public :: pool_study_means_c
    public :: construct_neighborhoods_expert_c
    public :: construct_neighborhoods_c

contains

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):compute_gene_means(subroutine)]]
    subroutine compute_gene_means_c(&
            n_genes,&
            n_reps,&
            expr,&
            means,&
            ierr&
        ) bind(C, name="compute_gene_means_c")
        use tox_data_integration_preprocessing, only: compute_gene_means

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes in the study
        integer(c_int), intent(in), target :: n_reps
            !! Number of biological replicates in the study
        real(c_double), dimension(n_reps, n_genes), intent(in), target :: expr
            !! Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(out), target :: means
            !! Per-gene mean expression values
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):compute_residuals(subroutine)]]
    subroutine compute_residuals_c(&
            n_genes,&
            n_reps,&
            expr,&
            means,&
            resid,&
            ierr&
        ) bind(C, name="compute_residuals_c")
        use tox_data_integration_preprocessing, only: compute_residuals

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes in the study
        integer(c_int), intent(in), target :: n_reps
            !! Number of biological replicates in the study
        real(c_double), dimension(n_reps, n_genes), intent(in), target :: expr
            !! Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(in), target :: means
            !! Per-gene mean expression values; NaN where every replicate of a gene was NaN
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps, n_genes), intent(out), target :: resid
            !! Matrix of signed residuals
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_means(subroutine)]]
    !| This takes the pool already built; `pool_study_means` pools the means of two studies
    !| first, if that is what is at hand.
    subroutine pool_means_expert_c(&
            pooled_means,&
            pooled_means_perm,&
            pool_size,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_means_expert_c")
        use tox_data_integration_preprocessing, only: pool_means

        integer(c_int), intent(in), target :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(pool_size), intent(in), target :: pooled_means
            !! Pooled means
            !! NaN is permitted for this value.
        integer(c_int), dimension(pool_size), intent(in), target :: pooled_means_perm
            !! Sorting permutation for `pooled_means`
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_means_alloc(subroutine)]]
    !| This takes the pool already built; `pool_study_means` pools the means of two studies
    !| first, if that is what is at hand.
    subroutine pool_means_c(&
            pooled_means,&
            pool_size,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_means_c")
        use tox_data_integration_preprocessing, only: pool_means_alloc

        integer(c_int), intent(in), target :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(pool_size), intent(in), target :: pooled_means
            !! Pooled means
            !! NaN is permitted for this value.
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(pool_size)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_ARRAY_NON_NULL(pooled_means, pool_size)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)

        call pool_means_alloc(&
            pooled_means = pooled_means,&
            pool_size = pool_size,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_means_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_study_means(subroutine)]]
    !| Concatenates the two studies' means, sorts the pool, and turns it into reference
    !| points exactly as `pool_means` does.
    subroutine pool_study_means_expert_c(&
            n_genes_S1,&
            mean_S1,&
            n_genes_S2,&
            mean_S2,&
            n_points,&
            tmp_pooled_means,&
            tmp_pooled_means_perm,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_study_means_expert_c")
        use tox_data_integration_preprocessing, only: pool_study_means

        integer(c_int), intent(in), target :: n_genes_S1
            !! Number of genes in study S1
        integer(c_int), intent(in), target :: n_genes_S2
            !! Number of genes in study S2
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(n_genes_S1), intent(in), target :: mean_S1
            !! Per-gene mean expression values of study S1
            !! NaN is permitted for this value.
        real(c_double), dimension(n_genes_S2), intent(in), target :: mean_S2
            !! Per-gene mean expression values of study S2
            !! NaN is permitted for this value.
        real(c_double), dimension(n_genes_S1+n_genes_S2), intent(out), target :: tmp_pooled_means
            !! Work array holding the concatenated means of both studies
        integer(c_int), dimension(n_genes_S1+n_genes_S2), intent(out), target :: tmp_pooled_means_perm
            !! Work array for the permutation that sorts `tmp_pooled_means`
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes_S1)
        M_CHECK_NON_NULL(n_genes_S2)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_ARRAY_NON_NULL(mean_S1, n_genes_S1)
        M_CHECK_ARRAY_NON_NULL(mean_S2, n_genes_S2)
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_means, (n_genes_S1+n_genes_S2))
        M_CHECK_ARRAY_NON_NULL(tmp_pooled_means_perm, (n_genes_S1+n_genes_S2))
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)

        call pool_study_means(&
            n_genes_S1 = n_genes_S1,&
            mean_S1 = mean_S1,&
            n_genes_S2 = n_genes_S2,&
            mean_S2 = mean_S2,&
            n_points = n_points,&
            tmp_pooled_means = tmp_pooled_means,&
            tmp_pooled_means_perm = tmp_pooled_means_perm,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_study_means_expert_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_study_means_alloc(subroutine)]]
    !| Concatenates the two studies' means, sorts the pool, and turns it into reference
    !| points exactly as `pool_means` does.
    subroutine pool_study_means_c(&
            n_genes_S1,&
            mean_S1,&
            n_genes_S2,&
            mean_S2,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_study_means_c")
        use tox_data_integration_preprocessing, only: pool_study_means_alloc

        integer(c_int), intent(in), target :: n_genes_S1
            !! Number of genes in study S1
        integer(c_int), intent(in), target :: n_genes_S2
            !! Number of genes in study S2
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(n_genes_S1), intent(in), target :: mean_S1
            !! Per-gene mean expression values of study S1
            !! NaN is permitted for this value.
        real(c_double), dimension(n_genes_S2), intent(in), target :: mean_S2
            !! Per-gene mean expression values of study S2
            !! NaN is permitted for this value.
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes_S1)
        M_CHECK_NON_NULL(n_genes_S2)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_ARRAY_NON_NULL(mean_S1, n_genes_S1)
        M_CHECK_ARRAY_NON_NULL(mean_S2, n_genes_S2)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)

        call pool_study_means_alloc(&
            n_genes_S1 = n_genes_S1,&
            mean_S1 = mean_S1,&
            n_genes_S2 = n_genes_S2,&
            mean_S2 = mean_S2,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_study_means_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_neighborhoods(subroutine)]]
    subroutine construct_neighborhoods_expert_c(&
            n_points,&
            x_star,&
            n_genes_S,&
            mean_S,&
            n_reps_S,&
            resid_S,&
            tmp_distances,&
            tmp_distances_perm,&
            neighborhood_residuals,&
            neighborhood_indices,&
            n_neighbors,&
            ierr&
        ) bind(C, name="construct_neighborhoods_expert_c")
        use tox_data_integration_preprocessing, only: construct_neighborhoods

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_reps_S
            !! Number of biological replicates in the study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
            !! cannot exceed the number of genes with a defined mean
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
            !! It is recommended to compute this with
            !! [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]].
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S, n_genes_S), intent(in), target :: resid_S
            !! Matrix of signed residuals
            !! NaN is permitted for this value.
        real(c_double), dimension(n_genes_S), intent(out), target :: tmp_distances
            !! Distances work array
        integer(c_int), dimension(n_genes_S), intent(out), target :: tmp_distances_perm
            !! Work array for the permutation that sorts `tmp_distances`
        real(c_double), dimension(n_reps_S, n_neighbors, n_points), intent(out), target :: neighborhood_residuals
            !! Collection of residual vectors for each neighborhood
        integer(c_int), dimension(n_neighbors, n_points), intent(out), target :: neighborhood_indices
            !! Indices of selected neighborhood genes
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_genes_S)
        M_CHECK_NON_NULL(n_reps_S)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_S, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(resid_S, n_reps_S * n_genes_S)
        M_CHECK_ARRAY_NON_NULL(tmp_distances, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(tmp_distances_perm, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(neighborhood_residuals, n_reps_S * n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points)

        call construct_neighborhoods(&
            n_points = n_points,&
            x_star = x_star,&
            n_genes_S = n_genes_S,&
            mean_S = mean_S,&
            n_reps_S = n_reps_S,&
            resid_S = resid_S,&
            tmp_distances = tmp_distances,&
            tmp_distances_perm = tmp_distances_perm,&
            neighborhood_residuals = neighborhood_residuals,&
            neighborhood_indices = neighborhood_indices,&
            n_neighbors = n_neighbors,&
            ierr = ierr&
        )
    end subroutine construct_neighborhoods_expert_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_neighborhoods_alloc(subroutine)]]
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
        use tox_data_integration_preprocessing, only: construct_neighborhoods_alloc

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_reps_S
            !! Number of biological replicates in the study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
            !! cannot exceed the number of genes with a defined mean
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
            !! It is recommended to compute this with
            !! [[tox_data_integration_preprocessing_kernel(module):calc_neighborhood_size(function)]].
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S, n_genes_S), intent(in), target :: resid_S
            !! Matrix of signed residuals
            !! NaN is permitted for this value.
        real(c_double), dimension(n_reps_S, n_neighbors, n_points), intent(out), target :: neighborhood_residuals
            !! Collection of residual vectors for each neighborhood
        integer(c_int), dimension(n_neighbors, n_points), intent(out), target :: neighborhood_indices
            !! Indices of selected neighborhood genes
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

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

end module tox_data_integration_preprocessing_c
#endif
