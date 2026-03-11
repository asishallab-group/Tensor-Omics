#include "macros.h"

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| This module implements the pipeline to obtain neighborhood residuals from expression vectors, to be used for JCT based data integration.
module tox_data_integration_preprocessing
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_is_finite, ieee_positive_inf
    use f42_utils, only: heapsort_real, calc_percentile_helper, clamp, binary_search_insertion
    use tox_errors, only: validate_all_in_range_real, validate_in_range_int, is_err, set_ok, validate_dimension_size, ERR_ALLOC_FAIL, set_err, validate_all_in_range_int
    implicit none

contains

    !> Compute per-gene mean expression, ignoring NaN values
    !|
    !| @note
    !| The means of all studies should be in contiguous memory afterwards, so for using this subroutine pass `means` as `means(:, study_idx)`
    !| @endnote
    pure subroutine compute_gene_means(expr, n_genes, n_reps, means, max_n_genes_all_studies, ierr)
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

        call set_ok(ierr)

        call validate_dimension_size(max_n_genes_all_studies, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_reps, ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, min=n_genes)
        ! expression can contain NaN
        if (is_err(ierr)) return

        call compute_gene_means_helper(expr, n_genes, n_reps, means, max_n_genes_all_studies)
    end subroutine compute_gene_means

    !> (no input validation) Compute per-gene mean expression, ignoring NaN values
    !|
    !| @note
    !| The means of all studies should be in contiguous memory afterwards, so for using this subroutine pass `means` as `means(:, study_idx)`
    !| @endnote
    pure subroutine compute_gene_means_helper(expr, n_genes, n_reps, means, max_n_genes_all_studies)
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

        integer(int32) :: i_gene, i_rep, n_included
        real(real64) :: sum_val, expr_val

        ! Use do concurrent for parallelization across genes
        do concurrent(i_gene=1:max_n_genes_all_studies) local(sum_val, n_included)
            if (i_gene > n_genes) then
                means(i_gene) = M_NAN
            else
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

                if (n_included > 0) then
                    means(i_gene) = sum_val/real(n_included, real64)
                else
                    means(i_gene) = M_NAN
                end if
            end if
        end do
    end subroutine compute_gene_means_helper

    !> Compute signed residuals (centering by mean).
    !|
    !| @note
    !| The residuals of all studies should be in contiguous memory afterwards, so for using this subroutine pass `expr` as `expr(:, :, study_idx)`
    !| @endnote
    pure subroutine compute_residuals(expr, n_genes, n_reps, max_n_genes_all_studies, max_n_reps_all_studies, means, resid, ierr)
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

        call set_ok(ierr)
        call validate_dimension_size(max_n_genes_all_studies, ierr)
        call validate_dimension_size(max_n_reps_all_studies, ierr)
        call validate_dimension_size(n_genes, ierr)
        call validate_dimension_size(n_reps, ierr)
        call validate_in_range_int(max_n_genes_all_studies, ierr, min=n_genes)
        call validate_in_range_int(max_n_reps_all_studies, ierr, min=n_reps)
        ! family means and expr containing NaN is expected behaviour
        if (is_err(ierr)) return

        call compute_residuals_helper(expr, n_genes, n_reps, max_n_genes_all_studies, max_n_reps_all_studies, means, resid)
    end subroutine compute_residuals

    !> (no input validation) Compute signed residuals (centering by mean)
    !|
    !| @note
    !| The residuals of all studies should be in contiguous memory afterwards, so for using this subroutine pass `expr` as `expr(:, :, study_idx)`
    !| @endnote
    pure subroutine compute_residuals_helper(expr, n_genes, n_reps, max_n_genes_all_studies, max_n_reps_all_studies, means, resid)
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

    pure integer(int32) function find_last_non_nan(arr, arr_perm, n_elements) result(idx)
        integer(int32), intent(in) :: n_elements
            !! Number of elements in `arr`
        real(real64), dimension(n_elements), intent(in) :: arr
            !! Array to find the last non-NaN index
        integer(int32), dimension(n_elements), intent(in) :: arr_perm
            !! Sorting permutation for `arr`

        ! Search for positive infinity. As NaN is always sorted to the end, the insertion index will point to the first NaN.
        idx = binary_search_insertion(arr, arr_perm, M_POS_INF) - 1
    end function find_last_non_nan

    !> Pool per-gene mean expression values across studies
    pure subroutine pool_means_alloc(means, n_studies, max_n_genes_all_studies, n_points, n_pool, x_star, ierr)
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

        integer(int32), allocatable :: perm(:)
        integer(int32) :: i_gene, pool_idx, pool_size

        call set_ok(ierr)

        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=2_int32)
        ! mean values can contain NaN
        if (is_err(ierr)) return

        ! Allocate arrays for pooled means
        pool_size = size(means, kind=int32)
        M_ALLOCATE(perm(pool_size))

        do concurrent(i_gene=1:pool_size) shared(perm)
            perm(i_gene) = i_gene
        end do

        call heapsort_real(means, perm)

        call pool_means(means, perm, pool_size, n_points, n_pool, x_star, ierr)
    end subroutine pool_means_alloc

    !> Pool per-gene mean expression values across studies
    pure subroutine pool_means(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star, ierr)
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

        call set_ok(ierr)

        call validate_dimension_size(pool_size, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_all_in_range_int(pooled_means_perm, pool_size, ierr, arg_pos=2_int32, min=1_int32, max=pool_size)
        ! mean values can contain NaN
        if (is_err(ierr)) return

        call pool_means_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
    end subroutine pool_means

    !> (no input validation) Pool per-gene mean expression values across studies
    pure subroutine pool_means_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
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

        n_pool = find_last_non_nan(pooled_means, pooled_means_perm, pool_size)

        call pool_means_n_pool_input_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
    end subroutine pool_means_helper

    !> (no input validation) Pool per-gene mean expression values across studies
    pure subroutine pool_means_n_pool_input_helper(pooled_means, pooled_means_perm, pool_size, n_points, n_pool, x_star)
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

        integer(int32) :: i_point
        real(real64) :: quantile_level

        if (n_pool <= 0) then
            x_star = M_NAN
        else
            ! Compute reference points as empirical quantiles using the permutation
            do concurrent(i_point=1:n_points) local(quantile_level) shared(n_points, pooled_means, pooled_means_perm, n_pool, x_star)
                quantile_level = real(i_point, real64)/real(n_points + 1, real64)*100.0_real64

                ! Use calc_percentile to compute the value
                call calc_percentile_helper(pooled_means(:n_pool), pooled_means_perm(:n_pool), quantile_level, x_star(i_point))
            end do
        end if
    end subroutine pool_means_n_pool_input_helper

    !> Calculate the number of neighbors to be used for [[tox_data_integration(module):construct_neighborhoods(interface)]].
    !|
    !| The `desired_size` works as upper limit, as the actual neighborhood size might be lower due to few genes with non-NaN mean.
    pure function calc_neighborhood_size(n_pool, n_points, n_genes_S, mean_S, desired_size) result(n_neighbors)
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
    pure subroutine construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)
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

        integer(int32) :: i_gene
        integer(int32), dimension(:), allocatable :: mean_S_perm

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=7_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(mean_S_perm(n_genes_S))

        ! Initialize perm
        do concurrent (i_gene = 1:n_genes_S) shared(mean_S_perm)
            mean_S_perm(i_gene) = i_gene
        end do

        call heapsort_real(mean_S, mean_S_perm)
        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, neighborhood_range, n_neighbors)
    end subroutine construct_neighborhoods_alloc

    !> (no input validation) Construct neighborhood-based residual sets (kNN)
    pure subroutine construct_neighborhoods(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)
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

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=8_int32)
        call validate_all_in_range_int(mean_S_perm, n_genes_S, ierr, arg_pos=5_int32, min=1_int32, max=n_genes_S)

        if (is_err(ierr)) return

        call construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, neighborhood_range, n_neighbors)
    end subroutine construct_neighborhoods

    !> (no input validation) Construct neighborhood-based residual sets (kNN)
    pure subroutine construct_neighborhoods_helper(n_points, x_star, n_genes_S, mean_S, mean_S_perm, neighborhood_residuals, neighborhood_range, n_neighbors)
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

        integer(int32) :: i_point, i_neighbor, gene_idx, left_gene, right_gene, x_star_idx
        real(real64) :: x_star_val

        ! Process each reference point
        ! do concurrent (i_point = 1:n_points) local(x_star_val, x_star_idx, left_gene, right_gene) shared(x_star, neighborhood_residuals, n_neighbors, mean_S, mean_S_perm)
        do i_point = 1, n_points

            x_star_val = x_star(i_point)

            ! If x_star is NaN, all gene expressions/residuals are NaN -> first element
            if (ieee_is_nan(x_star_val)) then
                x_star_idx = 1
            else
                x_star_idx = binary_search_insertion(mean_S, mean_S_perm, x_star_val)
            end if

            if (x_star_idx == 1) then
                neighborhood_range(1, i_point) = 1_int32
                if (ieee_is_nan(x_star_val)) then
                    neighborhood_range(1, i_point) = 1_int32
                else
                    neighborhood_range(2, i_point) = n_genes_S
                end if

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
                        ! if right side is closer than left side of x_star, take right side, else left side
                        ! Note: In the sorted array, NaNs are always last -> right side
                        ! -> if condition is false because of NaN it is because right value is NaN or both
                        if (mean_S(mean_S_perm(right_gene)) - x_star_val < x_star_val - mean_S(mean_S_perm(left_gene))) then
                            neighborhood_residuals(i_neighbor, i_point) = mean_S_perm(right_gene)
                            right_gene = right_gene + 1
                        else
                            neighborhood_residuals(i_neighbor, i_point) = mean_S_perm(left_gene)
                            left_gene = left_gene - 1
                        end if
                    end if
                end do
                neighborhood_range(1, i_point) = left_gene + 1
                neighborhood_range(2, i_point) = right_gene - 1
            end if
        end do
    end subroutine construct_neighborhoods_helper
end module tox_data_integration_preprocessing