#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| This module implements the pipeline to obtain neighborhood residuals from expression vectors, to be used for JCT based data integration.
submodule (tox_data_integration) tox_data_integration_preprocessing
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_is_finite
    use f42_utils, only: heapsort_real, calc_percentile_helper, clamp, binary_search_insertion, LOG_2
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, is_err, set_ok, validate_dimension_size, ERR_ALLOC_FAIL, set_err, validate_all_in_range_int
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

    !> Compute signed residuals (centering by mean).
    !|
    !| @note
    !| The residuals of all studies should be in contiguous memory afterwards, so for using this subroutine pass `expr` as `expr(:, :, study_idx)`
    !| @endnote
    pure module subroutine compute_residuals(expr, n_genes, n_reps, max_n_genes_all_studies, max_n_reps_all_studies, means, resid, ierr)
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
    !|
    !| @note
    !| The residuals of all studies should be in contiguous memory afterwards, so for using this subroutine pass `expr` as `expr(:, :, study_idx)`
    !| @endnote
    pure module subroutine compute_residuals_helper(expr, n_genes, n_reps, max_n_genes_all_studies, max_n_reps_all_studies, means, resid)
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

        integer(int32) :: i_gene, i_rep, m_genes, m_reps

        m_genes = max(n_genes, max_n_genes_all_studies)
        m_reps = max(n_reps, max_n_reps_all_studies)
        resid = M_NAN

        do concurrent(i_gene=1:m_genes)
            do concurrent(i_rep=1:m_reps) shared(expr, resid, means, i_gene)
                if (.not. ieee_is_nan(expr(i_rep, i_gene))) then
                    resid(i_rep, i_gene) = expr(i_rep, i_gene) - means(i_gene)
                end if
            end do
        end do

    end subroutine compute_residuals_helper

    pure module function find_last_non_nan(arr, arr_perm, n_elements) result(idx)
        real(real64), dimension(n_elements), intent(in) :: arr
            !! Array to find the last non-NaN index
        integer(int32), dimension(n_elements), intent(in) :: arr_perm
            !! Sorting permutation for `arr`

        integer(int32) :: i_element

        idx = n_elements

        ! NaN is always last -> find last non-NaN index for percentile calculation
        do i_element = n_elements, 1, -1
            if (ieee_is_nan(arr(arr_perm(i_element)))) then
                idx = idx - 1
            else
                exit
            end if
        end do
    end function find_last_non_nan

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

        n_pool = find_last_non_nan(pooled_means, pooled_means_perm, pool_size)

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
    pure module subroutine construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, neighborhood_residuals, n_neighbors, ierr)
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
        real(real64), intent(out) :: neighborhood_residuals(n_neighbors, n_points)
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note 
            !! All indices in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case `n_genes_S` is lower than `n_neighbors`,
            !! remaining indices will be filled with the ones from `n_genes_S+1...n_neighbors`
            !! @endnote
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_gene
        integer(int32), dimension(:), allocatable :: mean_S_perm

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=6_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(mean_S_perm(n_genes_S))

        ! Initialize perm
        do concurrent (i_gene = 1:n_genes_S) shared(mean_S_perm)
            mean_S_perm(i_gene) = i_gene
        end do

        call heapsort_real(mean_S, mean_S_perm)
        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, n_neighbors)
    end subroutine construct_neighborhoods_alloc

    !> (no input validation) Construct neighborhood-based residual sets (kNN)
    pure module subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, n_neighbors, ierr)
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
        real(real64), intent(out) :: neighborhood_residuals(n_neighbors, n_points)
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note 
            !! All indices in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case `n_genes_S` is lower than `n_neighbors`,
            !! remaining indices will be filled with the ones from `n_genes_S+1...n_neighbors`
            !! @endnote
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=6_int32)
        call validate_all_in_range_int(mean_S_perm, n_genes_S, ierr, min=1_int32, max=n_genes_S)

        if (is_err(ierr)) return

        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, n_neighbors)
    end subroutine construct_neighborhoods

    !> (no input validation) Construct neighborhood-based residual sets (kNN)
    pure module subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, n_neighbors)
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
        real(real64), intent(out) :: neighborhood_residuals(n_neighbors, n_points)
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note 
            !! All indices in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case `n_genes_S` is lower than `n_neighbors`,
            !! remaining indices will be filled with the ones from `n_genes_S+1...n_neighbors`
            !! @endnote

        integer(int32) :: i_point, i_neighbor, gene_idx, left_gene, right_gene
        real(real64) :: x_star_val

        ! Process each reference point
        do concurrent (i_point = 1:n_points) local(x_star_val, x_star_idx, left_gene, right_gene) shared(x_star, neighborhood_residuals, n_neighbors, mean_S, mean_S_perm)

            x_star_val = x_star(i_point)

            ! If x_star is NaN, all gene expressions/residuals are NaN -> first element
            if (ieee_is_nan(x_star_val)) then
                x_star_idx = 1
            else
                x_star_idx = binary_search_insertion(mean_S, mean_S_perm, x_star_val)
            end if

            if (x_star_idx == 1) then
                ! The first `n_neighbors` genes are the closest to `x_star_val`
                do concurrent (i_neighbor = 1:n_neighbors) local(gene_idx) shared(n_genes_S, neighborhood_residuals, i_point)
                    if (i_neighbor <= n_genes_S) then
                        gene_idx = mean_S_perm(i_neighbor)
                    else
                        gene_idx = i_neighbor
                    end if
                    neighborhood_residuals(gene_idx, i_point) = gene_idx
                end do
            else
                ! Collect the closest values around x_star_val
                left_gene = x_star_idx - 1
                right_gene = x_star_idx
                do i_neighbor = 1, n_neighbors
                    ! if no values lower than x_star are left, fill with right side
                    if (left_gene < 1) then
                        ! In case both indices are out of range, fill up to `n_neighbors`
                        if (right_gene > n_genes_S) then
                            neighborhood_residuals(i_neighbor, i_point) = i_neighbor
                        else
                            neighborhood_residuals(i_neighbor, i_point) = mean_S_perm(right_gene)
                            right_gene = right_gene + 1
                        end if
                    ! if no values higher than x_star are left, fill with left side
                    else if (right_gene > n_genes_S) then
                        neighborhood_residuals(i_neighbor, i_point) = mean_S_perm(left_gene)
                        left_gene = left_gene - 1
                    else
                        ! if left side is closer than right side of x_star, take left side, else right side
                        if (x_star_val - mean_S(mean_S_perm(left_gene)) < mean_S(mean_S_perm(right_gene)) - x_star_val) then
                            neighborhood_residuals(i_neighbor, i_point) = mean_S_perm(left_gene)
                            left_gene = left_gene - 1
                        else
                            neighborhood_residuals(i_neighbor, i_point) = mean_S_perm(right_gene)
                            right_gene = right_gene + 1
                        end if
                    end if
                end do
            end if
        end do
    end subroutine construct_neighborhoods_helper
end submodule tox_data_integration_preprocessing