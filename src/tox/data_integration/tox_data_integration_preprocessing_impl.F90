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
    use f42_sort_impl, only: sort_array_heapsort, binary_search_insertion
    use f42_stats_impl, only: calc_percentile_impl
    M_IMPLICIT_NONE
    private
    public :: compute_gene_means_impl, compute_residuals_impl, pool_means_impl, pool_study_means_impl, &
             calc_neighborhood_size, construct_neighborhoods_impl, construct_neighborhoods_ranged_impl

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
                ! Fraction in [0,1] as expected by calc_percentile_impl. Routed through a
                ! percentage round-trip (*100.0, then /100.0 inside calc_percentile_rank) rather
                ! than written as the equivalent direct fraction: 125-stabilize-jscomp's own
                ! percentile helper took a 0-100 percentage, so its quantile_level carried this
                ! same intermediate rounding. Simplifying this to a bare fraction changes nothing
                ! mathematically but rounds to a different double near many rank boundaries --
                ! confirmed to flip which gene a neighborhood search selects on real, tied data,
                ! cascading into a measurable global_js_divergence difference from 125. Keep the
                ! round-trip so this stays bit-for-bit reproducible with 125.
                quantile_level = (real(i_point, real64)/real(n_points + 1, real64)*100.0_real64)/100.0_real64

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

    !> summary: Construct neighborhood-based residual sets (kNN), with a `[min_idx, max_idx]` range per reference point
    !| AUTHOR_LASZLO_LANG
    !| Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
    !| two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl(subroutine)]]'s
    !| distance-sort algorithm above. In addition to the neighborhood gene indices, this also
    !| reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
    !| genes tied with the span's edge value are never split from it), which a candidate
    !| admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
    !| without needing the gathered residuals.
    pure subroutine construct_neighborhoods_ranged_impl(n_points, x_star, n_genes_S, mean_S, mean_S_perm, n_neighbors, &
                                                          neighborhood_indices, neighborhood_range)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors to select per reference point
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: x_star(n_points)
            !! Mean-expression reference points
            !! DM_ALLOW_NAN
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: mean_S_perm(n_genes_S)
            !! Sorting permutation for `mean_S`
        integer(int32), intent(out) :: neighborhood_indices(n_neighbors, n_points)
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note
            !! All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
            !! `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
            !! ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
            !! limitation, ported as-is from 125-stabilize-jscomp).
            !! @endnote
        integer(int32), intent(out) :: neighborhood_range(2, n_points)
            !! For each reference point, the `[min_idx, max_idx]` of the included genes. The
            !! index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
            !! would be the min value. In case of duplicate means, `min_idx` points to the first
            !! appearance of the value and `max_idx` to the last, so even though their related
            !! mean value is the min/max in the neighborhood, the actual gene might not be
            !! included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`

        integer(int32) :: i_point, i_neighbor, gene_idx, left_gene, right_gene, x_star_idx, last_gene
        real(real64) :: x_star_val

        ! Process each reference point
        do concurrent(i_point=1:n_points) local(last_gene, x_star_val, x_star_idx, left_gene, right_gene) &
                shared(x_star, neighborhood_indices, n_neighbors, mean_S, mean_S_perm)
            associate ( &
                nhood_range_min => neighborhood_range(1, i_point), &
                nhood_range_max => neighborhood_range(2, i_point) &
                )
                x_star_val = x_star(i_point)

                ! If x_star is NaN, all gene expressions/residuals are NaN -> first element
                if (ieee_is_nan(x_star_val)) then
                    x_star_idx = 1
                else
                    x_star_idx = binary_search_insertion(mean_S, mean_S_perm, x_star_val)
                end if

                if (x_star_idx == 1) then
                    nhood_range_min = 1_int32
                    nhood_range_max = min(n_genes_S, n_neighbors)

                    ! The first `n_neighbors` genes are the closest to `x_star_val`
                    do concurrent(i_neighbor=1:n_neighbors) local(gene_idx) shared(n_genes_S, neighborhood_indices, i_point)
                        if (i_neighbor <= n_genes_S) then
                            gene_idx = mean_S_perm(i_neighbor)
                        else
                            gene_idx = i_neighbor
                        end if
                        neighborhood_indices(i_neighbor, i_point) = gene_idx
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
                                neighborhood_indices(i_neighbor, i_point) = i_neighbor
                            else
                                neighborhood_indices(i_neighbor, i_point) = mean_S_perm(right_gene)
                                right_gene = right_gene + 1
                            end if
                        ! if no values higher than x_star are left, fill with left side
                        else if (right_gene > n_genes_S) then
                            neighborhood_indices(i_neighbor, i_point) = mean_S_perm(left_gene)
                            left_gene = left_gene - 1
                        else
                            ! if right side is closer than left side of x_star, take right side, else left side
                            ! Note: In the sorted array, NaNs are always last -> right side
                            ! -> if condition is false because of NaN it is because right value is NaN or both
                            if (mean_S(mean_S_perm(right_gene)) - x_star_val < x_star_val - mean_S(mean_S_perm(left_gene))) then
                                neighborhood_indices(i_neighbor, i_point) = mean_S_perm(right_gene)
                                right_gene = right_gene + 1
                            else
                                neighborhood_indices(i_neighbor, i_point) = mean_S_perm(left_gene)
                                left_gene = left_gene - 1
                            end if
                        end if
                    end do
                    nhood_range_min = left_gene + 1
                    nhood_range_max = right_gene - 1
                end if

                last_gene = nhood_range_min
                do i_neighbor = last_gene - 1, 1, -1
                    if (mean_S(mean_S_perm(i_neighbor)) == mean_S(mean_S_perm(last_gene))) then
                        last_gene = i_neighbor
                    else
                        exit
                    end if
                end do
                nhood_range_min = last_gene

                last_gene = nhood_range_max
                do i_neighbor = last_gene + 1, n_genes_S
                    if (mean_S(mean_S_perm(i_neighbor)) == mean_S(mean_S_perm(last_gene))) then
                        last_gene = i_neighbor
                    else
                        exit
                    end if
                end do
                nhood_range_max = last_gene
            end associate
        end do
    end subroutine construct_neighborhoods_ranged_impl
end module tox_data_integration_preprocessing_impl
