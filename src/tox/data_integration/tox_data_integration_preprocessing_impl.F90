#include <src/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| The step that turns expression vectors into the neighborhood residuals the rest of the test
!| consumes: gene-wise means, the signed deviation of each replicate from them, and the
!| neighborhoods of reference points those residuals are grouped into so the comparison is
!| conditioned on expression level rather than pooled across it.
!|
!| `calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own.
module tox_data_integration_preprocessing_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
    use f42_math_impl, only: clamp
    use f42_sort_impl, only: sort_array_heapsort
    use f42_stats_impl, only: calc_percentile_impl
    M_IMPLICIT_NONE

contains

    !> summary: Compute per-gene mean expression, ignoring NaN values
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine compute_gene_means_impl(n_genes, n_reps, expr, means)
        integer(int32), intent(in) :: n_genes
            !! Number of genes in the study
        integer(int32), intent(in) :: n_reps
            !! Number of biological replicates in the study
        real(real64), intent(in) :: expr(n_reps, n_genes)
            !! Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), intent(out) :: means(n_genes)
            !! Per-gene mean expression values

        integer(int32) :: i_gene, i_rep, n_included
        real(real64) :: sum_val, expr_val

        ! Use do concurrent for parallelization across genes
        do concurrent(i_gene=1:n_genes) local(sum_val, n_included)
            sum_val = 0.0_real64
            n_included = 0

            ! Count valid (non-NaN) replicates and compute sum
            do concurrent(i_rep=1:n_reps) local(expr_val) shared(expr) reduce(+:sum_val, n_included)
                expr_val = expr(i_rep, i_gene)
                if ((.not. ieee_is_nan(expr_val)) .and. ieee_is_finite(expr_val)) then
                    sum_val = sum_val + expr(i_rep, i_gene)
                    n_included = n_included + 1
                end if
            end do

            means(i_gene) = sum_val/real(n_included, real64)
        end do
    end subroutine compute_gene_means_impl

    !> summary: Compute signed residuals (centering by mean)
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine compute_residuals_impl(n_genes, n_reps, expr, means, resid)
        integer(int32), intent(in) :: n_genes
            !! Number of genes in the study
        integer(int32), intent(in) :: n_reps
            !! Number of biological replicates in the study
        real(real64), intent(in) :: expr(n_reps, n_genes)
            !! Expression matrix
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), intent(in) :: means(n_genes)
            !! Per-gene mean expression values; NaN where every replicate of a gene was NaN
            !! DM_ALLOW_NAN
        real(real64), intent(out) :: resid(n_reps, n_genes)
            !! Matrix of signed residuals

        integer(int32) :: i_gene, i_rep

        do concurrent(i_gene=1:n_genes)
            do concurrent(i_rep=1:n_reps) shared(expr, resid, means, i_gene)
                if ((.not. ieee_is_nan(expr(i_rep, i_gene))) .and. ieee_is_finite(expr(i_rep, i_gene))) then
                    resid(i_rep, i_gene) = expr(i_rep, i_gene) - means(i_gene)
                else
                    resid(i_rep, i_gene) = M_NAN
                end if
            end do
        end do

    end subroutine compute_residuals_impl

    !> summary: Turn a sorted pool of per-gene mean expression values into reference points
    !| AUTHOR_FRANZ_ERIC_SILL
    !| This takes the pool already built; `pool_study_means` pools the means of two studies
    !| first, if that is what is at hand.
    pure subroutine pool_means_impl(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
        integer(int32), intent(in), target :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), intent(in) :: pooled_means(pool_size)
            !! Pooled means
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: pooled_means_perm(pool_size)
            !! Sorting permutation for `pooled_means`
        real(real64), intent(out) :: x_star(n_points)
            !! Mean-expression reference points

        integer(int32) :: i_gene, i_point
        real(real64) :: quantile_level

        n_pool = size(pooled_means, kind=int32)

        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_gene = n_pool, 1, -1
            if (ieee_is_nan(pooled_means(pooled_means_perm(i_gene)))) then
                n_pool = n_pool - 1
            else
                exit
            end if
        end do

        if (n_pool == 0) then
            x_star = M_NAN
        else
            ! Compute reference points as empirical quantiles using the permutation
            do concurrent(i_point=1:n_points) local(quantile_level) shared(n_points, pooled_means, pooled_means_perm, n_pool, x_star)
                quantile_level = real(i_point, real64)/real(n_points + 1, real64)*100.0_real64

                ! Use calc_percentile to compute the value
                call calc_percentile_impl(pooled_means, size(pooled_means, kind=int32), pooled_means_perm, &
                                          quantile_level, x_star(i_point), n_considered=n_pool)
            end do
        end if
    end subroutine pool_means_impl

    !> summary: Pool the per-gene mean expression values of two studies into reference points
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Concatenates the two studies' means, sorts the pool, and turns it into reference
    !| points exactly as `pool_means` does.
    pure subroutine pool_study_means_impl(n_genes_S1, mean_S1, n_genes_S2, mean_S2, n_points, &
                                            tmp_pooled_means, tmp_pooled_means_perm, n_pool, x_star)
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study S1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study S2
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        real(real64), intent(in) :: mean_S1(n_genes_S1)
            !! Per-gene mean expression values of study S1
            !! DM_ALLOW_NAN
        real(real64), intent(in) :: mean_S2(n_genes_S2)
            !! Per-gene mean expression values of study S2
            !! DM_ALLOW_NAN
        real(real64), intent(out) :: tmp_pooled_means(n_genes_S1 + n_genes_S2)
            !! Work array holding the concatenated means of both studies
        integer(int32), intent(out) :: tmp_pooled_means_perm(n_genes_S1 + n_genes_S2)
            !! Work array for the permutation that sorts `tmp_pooled_means`
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), intent(out) :: x_star(n_points)
            !! Mean-expression reference points

        integer(int32) :: i_gene, pool_idx, pool_size

        pool_size = n_genes_S1 + n_genes_S2

        do concurrent(i_gene=1:n_genes_S1) shared(tmp_pooled_means, mean_S1)
            tmp_pooled_means(i_gene) = mean_S1(i_gene)
            tmp_pooled_means_perm(i_gene) = i_gene
        end do

        do concurrent(i_gene=1:n_genes_S2) local(pool_idx) shared(tmp_pooled_means, mean_S2)
            pool_idx = n_genes_S1 + i_gene
            tmp_pooled_means(pool_idx) = mean_S2(i_gene)
            tmp_pooled_means_perm(pool_idx) = pool_idx
        end do

        call sort_array_heapsort(tmp_pooled_means, tmp_pooled_means_perm)

        call pool_means_impl(tmp_pooled_means, tmp_pooled_means_perm, pool_size, n_points, n_pool, x_star)

    end subroutine pool_study_means_impl

    !> M_EXPORT_C
    !| summary: Calculate the number of neighbors to be used for constructing neighborhoods
    !| AUTHOR_FRANZ_ERIC_SILL
    !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower
    !| due to few genes with non-NaN mean.
    pure function calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S, desired_size) result(n_neighbors)
        integer(int32), intent(in) :: n_pool
            !! Total number of pooled mean-expression values across both studies
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
            !! DM_ALLOW_NAN
        integer(int32), intent(in), optional :: desired_size
            !! Optional desired neighborhood size
            !! DM_DEFAULT(1000)
        integer(int32) :: n_neighbors
            !! Calculated neighborhood size

        integer(int32) :: max_neighbors, i_gene, min_neighbors

        M_DEFAULT_VAL(desired_size, max_neighbors, 1000_int32)

        min_neighbors = 100_int32

        ! If less neighbors than default min are desired, take it
        if (max_neighbors < min_neighbors) then
            n_neighbors = max_neighbors

            ! Take at least `min_neighbors` neighbors. It could be lower if there are too many reference points -> low steps across x_star.
        else
            n_neighbors = int(clamp(n_pool/(2*n_points), min_val=min_neighbors, max_val=max_neighbors))
        end if

        ! don't take more neighbors than genes
        n_neighbors = min(n_neighbors, n_genes_S)

        ! Exclude genes with NaN mean, so those with all replicates NaN
        do concurrent(i_gene=1:n_genes_S) shared(mean_S) reduce(+:n_neighbors)
            if (ieee_is_nan(mean_S(i_gene))) then
                n_neighbors = n_neighbors - 1
            end if
        end do

        n_neighbors = max(0_int32, n_neighbors)
    end function calc_neighborhood_size

    !> summary: Construct neighborhood-based residual sets (kNN)
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine construct_neighborhoods_impl(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, &
                                                   tmp_distances, tmp_distances_perm, &
                                                   neighborhood_residuals, neighborhood_indices, n_neighbors)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_reps_S
            !! Number of biological replicates in the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
            !! cannot exceed the number of genes with a defined mean
            !! DM_MIN(1_int32)
            !! DM_MAX(count(.not. ieee_is_nan(mean_S), kind=int32))
            !! It is recommended to compute this with
            !! [[tox_data_integration_preprocessing_impl(module):calc_neighborhood_size(function)]].
        real(real64), intent(in) :: x_star(n_points)
            !! Mean-expression reference points
            !! DM_ALLOW_NAN
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
            !! DM_ALLOW_NAN
        real(real64), intent(in) :: resid_S(n_reps_S, n_genes_S)
            !! Matrix of signed residuals
            !! DM_ALLOW_NAN
        real(real64), intent(out) :: tmp_distances(n_genes_S)
            !! Distances work array
        integer(int32), intent(out) :: tmp_distances_perm(n_genes_S)
            !! Work array for the permutation that sorts `tmp_distances`
        real(real64), intent(out) :: neighborhood_residuals(n_reps_S, n_neighbors, n_points)
            !! Collection of residual vectors for each neighborhood
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
            !! Indices of selected neighborhood genes

        integer(int32) :: i_point, i_gene, gene_idx

        ! Process each reference point
        !TODO optimize: the kNN search per reference point is embarrassingly parallel (each i_point is independent), but this is a plain sequential `do` because `tmp_distances`/`tmp_distances_perm` are single shared scratch buffers reused/sorted in-place across iterations. Sizing them `(n_genes_S, n_points)` would let this become `do concurrent` and parallelize what is likely the hottest loop in the preprocessing pipeline.
        do i_point = 1, n_points

            ! Calculate distances.
            do concurrent(i_gene=1:n_genes_S) shared(tmp_distances, tmp_distances_perm, i_point, x_star, mean_S)
                tmp_distances(i_gene) = abs(mean_S(i_gene) - x_star(i_point))

                ! Initialize perm
                tmp_distances_perm(i_gene) = i_gene
            end do

            ! Sort distances using heapsort
            ! sort_array_heapsort will reorder tmp_distances_perm so that tmp_distances(tmp_distances_perm(1:n_genes_S)) is sorted
            ! the `n_neighbors` upper bound guarantees that the NaN `mean_S` indices are not included after sorting (they come after tmp_distances_perm(:n_neighbors))
            call sort_array_heapsort(tmp_distances, tmp_distances_perm)

            ! Store the n_neighbors nearest neighbor indices
            ! tmp_distances_perm(1:n_neighbors) now contain indices of genes with smallest distances
            do concurrent(i_gene=1:n_neighbors) local(gene_idx) shared(tmp_distances_perm, neighborhood_indices, neighborhood_residuals, resid_S)
                gene_idx = tmp_distances_perm(i_gene)  ! Get the actual gene index from the permutation vector

                neighborhood_indices(i_gene, i_point) = gene_idx

                neighborhood_residuals(:, i_gene, i_point) = resid_S(:, gene_idx)
            end do
        end do
    end subroutine construct_neighborhoods_impl
end module tox_data_integration_preprocessing_impl
