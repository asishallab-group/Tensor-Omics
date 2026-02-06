#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| This module implements the pipeline to obtain neighborhood residuals from expression vectors, to be used for JCT based data integration.
submodule (tox_data_integration) tox_data_integration_preprocessing
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_is_finite
    use f42_utils, only: heapsort_real, calc_percentile_helper, clamp
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, is_err, set_ok, validate_dimension_size, ERR_ALLOC_FAIL, set_err
    implicit none

contains

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

        call set_ok(ierr)

        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_reps, ierr)
        ! expression can contain NaN
        if (is_err(ierr)) return

        call compute_gene_means_helper(n_genes, n_reps, expr, means)
    end subroutine compute_gene_means

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
    end subroutine compute_gene_means_helper

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

        call set_ok(ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_reps, ierr)
        ! family means and expr containing NaN is expected behaviour
        if (is_err(ierr)) return

        call compute_residuals_helper(n_genes, n_reps, expr, means, resid)
    end subroutine compute_residuals

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

        integer(int32) :: i_gene, i_rep

        do concurrent(i_gene=1:n_genes)
            do concurrent(i_rep=1:n_reps) shared(expr, resid, means, i_gene)
                if (.not. ieee_is_nan(expr(i_rep, i_gene))) then
                    resid(i_rep, i_gene) = expr(i_rep, i_gene) - means(i_gene)
                else
                    resid(i_rep, i_gene) = M_NAN
                end if
            end do
        end do

    end subroutine compute_residuals_helper

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

        real(real64), allocatable :: pooled_means(:)
        integer(int32), allocatable :: perm(:)
        integer(int32) :: i_gene, pool_idx, pool_size

        call set_ok(ierr)

        call validate_dimension_size(n_genes_S1, ierr)
        call validate_dimension_size(n_genes_S2, ierr)
        ! mean values can contain NaN
        if (is_err(ierr)) return

        ! Allocate arrays for pooled means
        pool_size = n_genes_S1 + n_genes_S2
        M_ALLOCATE(pooled_means(pool_size))
        M_ALLOCATE(perm(pool_size))

        do concurrent(i_gene=1:n_genes_S1) shared(pooled_means, mean_S1)
            pooled_means(i_gene) = mean_S1(i_gene)
            perm(i_gene) = i_gene
        end do

        do concurrent(i_gene=1:n_genes_S2) local(pool_idx) shared(pooled_means, mean_S2)
            pool_idx = n_genes_S1 + i_gene
            pooled_means(pool_idx) = mean_S2(i_gene)
            perm(pool_idx) = pool_idx
        end do

        call heapsort_real(pooled_means, perm)

        call pool_means(pooled_means, perm, pool_size, n_points, n_pool, x_star, ierr)

    end subroutine pool_means_alloc

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

        call set_ok(ierr)

        call validate_dimension_size(pool_size, ierr)
        call validate_dimension_size(n_points, ierr)
        ! mean values can contain NaN
        if (is_err(ierr)) return

        call pool_means_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
    end subroutine pool_means

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
                call calc_percentile_helper(pooled_means(:n_pool), pooled_means_perm(:n_pool), quantile_level, x_star(i_point))
            end do
        end if
    end subroutine pool_means_helper

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

        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:), allocatable :: tmp_distances_perm

        call set_ok(ierr)
        call validate_dimension_size(n_neighbors, ierr)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_distances(n_genes_S))
        M_ALLOCATE(tmp_distances_perm(n_genes_S))

        call construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors, ierr)

    end subroutine construct_neighborhoods_alloc

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

        integer(int32) :: i_gene, n_possible_neighbors

        call set_ok(ierr)
        call validate_dimension_size(n_points, ierr)
        call validate_dimension_size(n_genes_S, ierr)
        call validate_dimension_size(n_reps_S, ierr)
        call validate_dimension_size(n_neighbors, ierr)

        n_possible_neighbors = 0_int32
        do concurrent (i_gene = 1:n_genes_S) shared(mean_S) reduce(+:n_possible_neighbors)
            if (.not. ieee_is_nan(mean_S(i_gene))) then
                n_possible_neighbors = n_possible_neighbors + 1
            end if
        end do

        call validate_in_range_int(n_neighbors, ierr, min=1_int32, max=n_possible_neighbors)

        if (is_err(ierr)) return

        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, n_reps_S, resid_S, tmp_distances, tmp_distances_perm, neighborhood_residuals, neighborhood_indices, n_neighbors)

    end subroutine construct_neighborhoods

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

        integer(int32) :: i_point, i_gene, gene_idx

        ! Process each reference point
        do i_point = 1, n_points

            ! Calculate distances.
            do concurrent(i_gene=1:n_genes_S) shared(tmp_distances, tmp_distances_perm, i_point, x_star, mean_S)
                tmp_distances(i_gene) = abs(mean_S(i_gene) - x_star(i_point))

                ! Initialize perm
                tmp_distances_perm(i_gene) = i_gene
            end do

            ! Sort distances using heapsort
            ! heapsort_real will reorder tmp_distances_perm so that tmp_distances(tmp_distances_perm(1:n_genes_S)) is sorted
            ! `calc_neighborhood_size` guarantees that the NaN `mean_S` indices are not included after sorting (they come after tmp_distances_perm(:n_neighbors))
            call heapsort_real(tmp_distances, tmp_distances_perm)

            ! Store the n_neighbors nearest neighbor indices
            ! tmp_distances_perm(1:n_neighbors) now contain indices of genes with smallest distances
            do concurrent(i_gene=1:n_neighbors) local(gene_idx) shared(tmp_distances_perm, neighborhood_indices, neighborhood_residuals, resid_S)
                gene_idx = tmp_distances_perm(i_gene)  ! Get the actual gene index from the permutation vector

                neighborhood_indices(i_gene, i_point) = gene_idx

                neighborhood_residuals(:, i_gene, i_point) = resid_S(:, gene_idx)
            end do
        end do
    end subroutine construct_neighborhoods_helper
end submodule tox_data_integration_preprocessing

!> C-compatible wrapper for [[tox_data_integration(module):compute_gene_means(interface)]]
pure subroutine compute_gene_means_c(n_genes, n_reps, expr, means, ierr) &
    bind(C, name="compute_gene_means_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: compute_gene_means
    M_USE_NULL_VALIDATION

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
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_reps)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(means)

    call compute_gene_means(n_genes, n_reps, expr, means, ierr)
end subroutine compute_gene_means_c

!> C-compatible wrapper for [[tox_data_integration(module):compute_residuals(interface)]]
pure subroutine compute_residuals_c(n_genes, n_reps, expr, means, resid, ierr) &
    bind(C, name="compute_residuals_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: compute_residuals
    M_USE_NULL_VALIDATION

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
    M_CHECK_NON_NULL(n_genes)
    M_CHECK_NON_NULL(n_reps)
    M_CHECK_NON_NULL(expr)
    M_CHECK_NON_NULL(means)
    M_CHECK_NON_NULL(resid)

    call compute_residuals(n_genes, n_reps, expr, means, resid, ierr)
end subroutine compute_residuals_c

!> C-compatible wrapper for [[tox_data_integration(module):pool_means_alloc(interface)]]
pure subroutine pool_means_c(n_genes_S1, mean_S1, n_genes_S2, mean_S2, &
    n_points, n_pool, x_star, ierr) bind(C, name="pool_means_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: pool_means_alloc
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_genes_S1
        !! Number of genes in study S1
    real(c_double), dimension(n_genes_S1), intent(in), target :: mean_S1
        !! Per-gene mean expression values
    integer(c_int), intent(in), target :: n_genes_S2
        !! Number of genes in study S2
    real(c_double), dimension(n_genes_S2), intent(in), target :: mean_S2
        !! Per-gene mean expression values
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points to define
    integer(c_int), intent(out), target :: n_pool
        !! Total number of included (non-NaN) pooled mean-expression values
    real(c_double), dimension(n_points), intent(out), target :: x_star
        !! Mean-expression reference points
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_genes_S1)
    M_CHECK_NON_NULL(mean_S1)
    M_CHECK_NON_NULL(n_genes_S2)
    M_CHECK_NON_NULL(mean_S2)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_pool)
    M_CHECK_NON_NULL(x_star)

    call pool_means_alloc(n_genes_S1, mean_S1, n_genes_S2, mean_S2, &
                          n_points, n_pool, x_star, ierr)
end subroutine pool_means_c

!> C-compatible wrapper for [[tox_data_integration(module):pool_means(interface)]]
pure subroutine pool_means_expert_c(pooled_means, pooled_means_perm, &
    pool_size, n_points, n_pool, x_star, ierr) &
    bind(C, name="pool_means_expert_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: pool_means
    M_USE_NULL_VALIDATION

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
    M_CHECK_NON_NULL(pool_size)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(pooled_means)
    M_CHECK_NON_NULL(pooled_means_perm)
    M_CHECK_NON_NULL(n_pool)
    M_CHECK_NON_NULL(x_star)

    call pool_means(pooled_means, pooled_means_perm, &
                    pool_size, n_points, n_pool, x_star, ierr)
end subroutine pool_means_expert_c

!> C-compatible wrapper for [[tox_data_integration(module):calc_neighborhood_size(interface)]]
pure subroutine calc_neighborhood_size_c(n_pool, n_points, n_genes_S, mean_S, &
    desired_size, n_neighbors, ierr) bind(C, name="calc_neighborhood_size_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: calc_neighborhood_size
    use tox_errors, only: set_ok
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_pool
        !! Total number of pooled mean-expression values across both studies
    integer(c_int), intent(in), target :: n_points
        !! Number of reference points
    integer(c_int), intent(in), target :: n_genes_S
        !! Number of genes in the current study
    real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
        !! Per-gene mean expression values
    integer(c_int), intent(in), target :: desired_size
        !! Optional desired neighborhood size, default=1000
    integer(c_int), intent(out), target :: n_neighbors
        !! Neighborhood size used (constant for all reference points)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_pool)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(n_genes_S)
    M_CHECK_NON_NULL(mean_S)
    M_CHECK_NON_NULL(desired_size)
    M_CHECK_NON_NULL(n_neighbors)

    if (desired_size > 0_c_int) then
        n_neighbors = calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S, desired_size)
    else
        n_neighbors = calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S)
    end if

    call set_ok(ierr)
end subroutine calc_neighborhood_size_c

!> C-compatible wrapper for [[tox_data_integration(module):construct_neighborhoods_alloc(interface)]]
pure subroutine construct_neighborhoods_c(n_points, x_star, n_genes_S, mean_S, &
    n_reps_S, resid_S, neighborhood_residuals, neighborhood_indices, n_neighbors, ierr) &
    bind(C, name="construct_neighborhoods_c")

    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_data_integration, only: construct_neighborhoods_alloc
    M_USE_NULL_VALIDATION

    integer(c_int), intent(in), target :: n_points
        !! Number of reference points
    real(c_double), dimension(n_points), intent(in), target :: x_star
        !! Mean-expression reference points
    integer(c_int), intent(in), target :: n_genes_S
        !! Number of genes in the current study
    real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
        !! Per-gene mean expression values
    integer(c_int), intent(in), target :: n_reps_S
        !! Number of biological replicates in the study
    real(c_double), dimension(n_reps_S, n_genes_S), intent(in), target :: resid_S
        !! Matrix of signed residuals
    real(c_double), dimension(n_reps_S, n_neighbors, n_points), intent(out), target :: neighborhood_residuals
        !! Collection of residual vectors for each neighborhood
    integer(c_int), dimension(n_neighbors, n_points), intent(out), target :: neighborhood_indices
        !! Indices of selected neighborhood genes
    integer(c_int), intent(in), target :: n_neighbors
        !! Number of neighbors
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(x_star)
    M_CHECK_NON_NULL(n_genes_S)
    M_CHECK_NON_NULL(mean_S)
    M_CHECK_NON_NULL(n_reps_S)
    M_CHECK_NON_NULL(resid_S)
    M_CHECK_NON_NULL(neighborhood_residuals)
    M_CHECK_NON_NULL(neighborhood_indices)
    M_CHECK_NON_NULL(n_neighbors)

    call construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, &
                                       n_reps_S, resid_S, neighborhood_residuals, &
                                       neighborhood_indices, n_neighbors, ierr)
end subroutine construct_neighborhoods_c
