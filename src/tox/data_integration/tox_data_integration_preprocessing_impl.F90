#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
!|
!| The step that turns expression vectors into the neighborhood residuals the rest of the test
!| consumes: gene-wise means, the signed deviation of each replicate from them, and the
!| neighborhoods of reference points those residuals are grouped into so the comparison is
!| conditioned on expression level rather than pooled across it.
!|
!| Besides neighborhoods of a fixed size, the adaptive construction
!| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]
!| grows each reference point's neighborhood while its residual dispersion stays stable, and
!| [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]
!| checks that consecutive neighborhoods overlap enough.
!|
!| `calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own, and
!| `gather_pooled_neighborhood_residuals` gathers a neighborhood's residuals of all studies.
module tox_data_integration_preprocessing_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
    use f42_math_impl, only: clamp
    use f42_sort_impl, only: sort_array_heapsort, binary_search_insertion, init_perm, sort_real_heapsort_expl_size
    use f42_stats_impl, only: calc_quantile_impl, calc_sorted_slice_mad
    use tox_errors, only: set_ok, set_err_once, is_err, validate_dimension_size, ERR_INVALID_INPUT
    M_IMPLICIT_NONE
    private
    public :: compute_gene_means_impl, compute_residuals_impl, pool_means_impl, pool_study_means_impl, &
             calc_neighborhood_size, construct_neighborhoods_impl, construct_neighborhoods_ranged_impl, &
             gather_pooled_neighborhood_residuals, grow_adaptive_neighborhoods, &
             construct_adaptive_neighborhoods_impl, materialize_pooled_neighborhood, check_neighborhood_overlaps_impl, &
             ADAPTIVE_MIN_VALID_RESIDUALS, ADAPTIVE_STOP_TAU, ADAPTIVE_STOP_K_MAX, ADAPTIVE_STOP_RESIDUAL_CAP, &
             ADAPTIVE_STOP_EXHAUSTED, ADAPTIVE_STOP_ZERO_DISPERSION, ADAPTIVE_STOP_TOO_FEW_RESIDUALS, &
             ADAPTIVE_STATUS_OK, ADAPTIVE_STATUS_TOO_FEW_MEANS, ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD, &
             ADAPTIVE_STATUS_TOO_FEW_RESIDUALS

    ! Output codes of construct_adaptive_neighborhoods_impl (not mode arguments). Each is a CM_
    ! macro first, so the parameter and the documentation quoting it share one literal.
#define CM_ADAPTIVE_MIN_VALID_RESIDUALS 10_int32
#define CM_ADAPTIVE_STOP_TAU 1_int32
#define CM_ADAPTIVE_STOP_K_MAX 2_int32
#define CM_ADAPTIVE_STOP_RESIDUAL_CAP 3_int32
#define CM_ADAPTIVE_STOP_EXHAUSTED 4_int32
#define CM_ADAPTIVE_STOP_ZERO_DISPERSION 5_int32
#define CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS 6_int32
#define CM_ADAPTIVE_STATUS_OK 0_int32
#define CM_ADAPTIVE_STATUS_TOO_FEW_MEANS 1_int32
#define CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD 2_int32
#define CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS 3_int32
    integer(int32), parameter :: ADAPTIVE_MIN_VALID_RESIDUALS = CM_ADAPTIVE_MIN_VALID_RESIDUALS
        !! Minimum number of non-NaN residuals a neighborhood must pool after its unconditional
        !! `k_start` phase for its dispersion-controlled growth to run
    integer(int32), parameter :: ADAPTIVE_STOP_TAU = CM_ADAPTIVE_STOP_TAU
        !! Stop reason: the next round would have raised the dispersion by more than `tau`
    integer(int32), parameter :: ADAPTIVE_STOP_K_MAX = CM_ADAPTIVE_STOP_K_MAX
        !! Stop reason: the neighborhood reached `k_max` pooled entries
    integer(int32), parameter :: ADAPTIVE_STOP_RESIDUAL_CAP = CM_ADAPTIVE_STOP_RESIDUAL_CAP
        !! Stop reason: the next entry would have pushed the pooled residuals past `max_pooled_residuals`
    integer(int32), parameter :: ADAPTIVE_STOP_EXHAUSTED = CM_ADAPTIVE_STOP_EXHAUSTED
        !! Stop reason: no pooled entry is left on either side of the neighborhood
    integer(int32), parameter :: ADAPTIVE_STOP_ZERO_DISPERSION = CM_ADAPTIVE_STOP_ZERO_DISPERSION
        !! Stop reason: the `k_start` entries' residuals are all zero, so no relative change is defined
    integer(int32), parameter :: ADAPTIVE_STOP_TOO_FEW_RESIDUALS = CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS
        !! Stop reason: fewer than ADAPTIVE_MIN_VALID_RESIDUALS non-NaN residuals after the
        !! `k_start` phase, so the adaptive phase never ran
    integer(int32), parameter :: ADAPTIVE_STATUS_OK = CM_ADAPTIVE_STATUS_OK
        !! Construction status: every neighborhood was built without a failure
    integer(int32), parameter :: ADAPTIVE_STATUS_TOO_FEW_MEANS = CM_ADAPTIVE_STATUS_TOO_FEW_MEANS
        !! Construction status: fewer non-NaN pooled means than `k_start`; no neighborhood was built
    integer(int32), parameter :: ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD = CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD
        !! Construction status: some neighborhood has fewer than `min_study_neighbors` entries of some study
    integer(int32), parameter :: ADAPTIVE_STATUS_TOO_FEW_RESIDUALS = CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS
        !! Construction status: some neighborhood stopped with ADAPTIVE_STOP_TOO_FEW_RESIDUALS

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
            !! DM_MIN(1_int32)
            !! DM_MAX(pool_size)
        real(real64), intent(out) :: x_star(n_points)
            !! Mean-expression reference points. If every pooled mean is NaN, every reference
            !! point is NaN, and `construct_neighborhoods_ranged` rejects such points.

        integer(int32) :: i_gene, i_point
        real(real64) :: quantile_level

        n_pool = size(pooled_means, kind=int32)

        ! NaN is always last -> find last non-NaN index for the quantile calculation
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
                ! Fraction in [0,1] as expected by calc_quantile_impl
                quantile_level = real(i_point, real64)/real(n_points + 1, real64)

                ! Use calc_quantile to compute the value
                call calc_quantile_impl(pooled_means, size(pooled_means, kind=int32), pooled_means_perm, &
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
            !! Mean-expression reference points. If every mean of both studies is NaN, every
            !! reference point is NaN, and `construct_neighborhoods_ranged` rejects such points.

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
            !! DM_DEFAULT(1000_int32)
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
        real(real64), intent(in) :: mean_S(n_genes_S)
            !! Per-gene mean expression values
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: mean_S_perm(n_genes_S)
            !! Sorting permutation for `mean_S`
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes_S)
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
        do concurrent(i_point=1:n_points) local(last_gene, x_star_val, x_star_idx, left_gene, right_gene, i_neighbor) &
                shared(x_star, neighborhood_indices, neighborhood_range, n_neighbors, n_genes_S, mean_S, mean_S_perm)
            associate ( &
                nhood_range_min => neighborhood_range(1, i_point), &
                nhood_range_max => neighborhood_range(2, i_point) &
                )
                x_star_val = x_star(i_point)

                ! A NaN reference point takes the first `n_neighbors` genes in sorted order. The
                ! generated wrappers reject a NaN `x_star`, so only a direct caller of this
                ! implementation can pass one: the pipeline does so only when every mean is NaN,
                ! which its empty-study checks reject
                if (ieee_is_nan(x_star_val)) then
                    x_star_idx = 1
                else
                    x_star_idx = binary_search_insertion(mean_S, mean_S_perm, x_star_val)
                end if

                if (x_star_idx == 1) then
                    nhood_range_min = 1_int32
                    nhood_range_max = min(n_genes_S, n_neighbors)

                    ! The first `n_neighbors` genes are the closest to `x_star_val`
                    do concurrent(i_neighbor=1:n_neighbors) local(gene_idx) shared(n_genes_S, neighborhood_indices, mean_S_perm, i_point)
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

    !> M_EXPORT_C
    !| summary: Pool one reference point's residuals across every neighbor and every study
    !| AUTHOR_LASZLO_LANG
    !| Gathers one reference point's residuals from every neighbor gene of every study into one flat
    !| array, given that point's per-study neighbor gene indices (one column of the
    !| `neighborhood_indices_all_studies(n_neighbors, n_points, n_studies)` of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]).
    !| Study `s` fills pool block `(s - 1)*max_n_reps_all_studies*n_neighbors + 1` to
    !| `s*max_n_reps_all_studies*n_neighbors`. This is the pooling
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] and
    !| [[tox_data_integration_js_comp_test_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| feed to the bin-count search, so a caller can pass the same input to
    !| [[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_occupancy_exhaustive_impl(interface)]]
    !| on real data.
    pure subroutine gather_pooled_neighborhood_residuals(residuals, max_n_reps_all_studies, max_n_genes_all_studies, &
                                                         n_neighbors, n_studies, neighborhood_indices_point, &
                                                         pooled_residuals, ierr)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        real(real64), intent(in) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
            !! Matrix of signed residuals per study, NaN allowed for missing values; not validated
        integer(int32), intent(in) :: neighborhood_indices_point(n_neighbors, n_studies)
            !! Gene indices of one reference point's neighborhood, per study; each in
            !! [1, max_n_genes_all_studies]
        real(real64), intent(out) :: pooled_residuals(max_n_reps_all_studies*n_neighbors*n_studies)
            !! This reference point's pooled residuals, laid out as a
            !! (max_n_reps_all_studies, n_neighbors, n_studies) array
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure

        integer(int32) :: i_study, i_neighbor, i_rep, gene_idx, n_predecessors

        call set_ok(ierr)
        call validate_dimension_size(max_n_reps_all_studies, ierr, arg_pos=2_int32)
        call validate_dimension_size(max_n_genes_all_studies, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_neighbors, ierr, arg_pos=4_int32)
        call validate_dimension_size(n_studies, ierr, arg_pos=5_int32)
        if (is_err(ierr)) return

        ! Every gene index subscripts residuals directly, so a bad one would read out of bounds.
        ! Sequential: set_err_once on the shared ierr is a read-then-write, not allowed in do concurrent.
        do i_study = 1, n_studies
            do i_neighbor = 1, n_neighbors
                gene_idx = neighborhood_indices_point(i_neighbor, i_study)
                if (gene_idx < 1 .or. gene_idx > max_n_genes_all_studies) then
                    call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)
                end if
            end do
        end do
        if (is_err(ierr)) return

        do concurrent(i_study=1:n_studies)
            do concurrent(i_neighbor=1:n_neighbors) local(gene_idx, n_predecessors) &
                    shared(neighborhood_indices_point, max_n_reps_all_studies, n_neighbors, i_study)
                gene_idx = neighborhood_indices_point(i_neighbor, i_study)
                n_predecessors = ((i_study - 1)*n_neighbors + (i_neighbor - 1))*max_n_reps_all_studies
                do concurrent(i_rep=1:max_n_reps_all_studies) shared(pooled_residuals, residuals, n_predecessors, &
                                                                     gene_idx, i_study)
                    pooled_residuals(n_predecessors + i_rep) = residuals(i_rep, gene_idx, i_study)
                end do
            end do
        end do
    end subroutine gather_pooled_neighborhood_residuals

    !> summary: Grow Issue #217's adaptive neighborhoods over an already-sorted pool of gene means
    !| AUTHOR_LASZLO_LANG
    !| The construction of
    !| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]],
    !| documented there, on an already sorted permutation (ascending, NaN last), so the adaptive
    !| parameter search can grow every candidate on one shared sort.
    !| `pooled_means` and `pooled_residuals` are `gene_means(G, S)` and `residuals(R, G, S)` seen
    !| through sequence association as `G*S` flat entries, entry `(s - 1)*G + g` being gene `g` of
    !| study `s`. Internal helper, not published; public because the adaptive parameter search calls
    !| it besides `construct_adaptive_neighborhoods_impl`.
    pure subroutine grow_adaptive_neighborhoods(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                pooled_means, pooled_residuals, pooled_means_perm, k_start, k_step, &
                                                k_max, tau, mad_distance_factor, max_pooled_residuals, &
                                                min_study_neighbors, n_points, x_star, pooled_neighborhood_range, &
                                                n_neighbors_per_point, stop_reason, neighborhood_dispersion, &
                                                neighborhood_mad, max_n_neighbors, construction_status)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        real(real64), intent(in) :: pooled_means(max_n_genes_all_studies*n_studies)
            !! Every study's gene means, flattened, NaN allowed
        real(real64), intent(in) :: pooled_residuals(max_n_reps_all_studies, max_n_genes_all_studies*n_studies)
            !! Every study's residuals, one column per flat entry of `pooled_means`, NaN allowed
        integer(int32), intent(in) :: pooled_means_perm(max_n_genes_all_studies*n_studies)
            !! Permutation sorting `pooled_means` ascending, NaN last
        integer(int32), intent(in) :: k_start
            !! Entries every neighborhood takes unconditionally, at least 1
        integer(int32), intent(in) :: k_step
            !! Entries staged per adaptive round, at least 1
        integer(int32), intent(in) :: k_max
            !! Largest neighborhood, at least `k_start`
        real(real64), intent(in) :: tau
            !! Largest accepted relative dispersion increase per round, at least 0
        real(real64), intent(in) :: mad_distance_factor
            !! Multiple of a neighborhood's MAD the next target lies beyond it, at least 0
        integer(int32), intent(in) :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals, 0 for none
        integer(int32), intent(in) :: min_study_neighbors
            !! Fewest entries of any one study a neighborhood may have, at least 1
        integer(int32), intent(out) :: n_points
            !! Number of neighborhoods built
        real(real64), intent(out) :: x_star(max_n_genes_all_studies*n_studies)
            !! Reference point of each neighborhood
        integer(int32), intent(out) :: pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies)
            !! `[first, last]` sorted pooled position of each neighborhood
        integer(int32), intent(out) :: n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies)
            !! Entries of each study in each neighborhood
        integer(int32), intent(out) :: stop_reason(max_n_genes_all_studies*n_studies)
            !! Why each neighborhood stopped growing
        real(real64), intent(out) :: neighborhood_dispersion(max_n_genes_all_studies*n_studies)
            !! Mean absolute residual of each final neighborhood, NaN without a non-NaN residual
        real(real64), intent(out) :: neighborhood_mad(max_n_genes_all_studies*n_studies)
            !! Raw MAD of each final neighborhood's pooled means
        integer(int32), intent(out) :: max_n_neighbors
            !! Largest entry count of one study in one neighborhood, 0 without a neighborhood
        integer(int32), intent(out) :: construction_status
            !! Overall outcome of the construction

        integer(int32) :: pool_size, n_valid_means, i_pos, seed_pos, first_pos, last_pos, left_pos, right_pos
        integer(int32) :: n_entries, n_valid, entry_pos, entry_n_valid, n_planned, n_staged, round_n_valid
        integer(int32) :: stage_left, stage_right, stage_first, stage_last, stop_code, target_pos, study_idx
        real(real64) :: x_star_val, sum_abs, entry_sum_abs, round_sum_abs, s_old, s_new, median, mad, target
        logical(c_bool) :: capped, has_too_few_residuals, has_empty_study

        pool_size = max_n_genes_all_studies*n_studies

        ! NaN sorts last: every position after n_valid_means holds a NaN mean.
        n_valid_means = pool_size
        do i_pos = pool_size, 1_int32, -1_int32
            if (ieee_is_nan(pooled_means(pooled_means_perm(i_pos)))) then
                n_valid_means = n_valid_means - 1_int32
            else
                exit
            end if
        end do

        n_points = 0_int32
        max_n_neighbors = 0_int32
        if (n_valid_means < k_start) then
            construction_status = ADAPTIVE_STATUS_TOO_FEW_MEANS
            return
        end if

        has_too_few_residuals = .false.
        has_empty_study = .false.
        seed_pos = 1_int32
        ! Sequential by construction: every seed is placed beyond the previous neighborhood.
        do
            n_points = n_points + 1_int32
            x_star_val = pooled_means(pooled_means_perm(seed_pos))
            first_pos = seed_pos
            last_pos = seed_pos
            left_pos = seed_pos - 1_int32
            right_pos = seed_pos + 1_int32
            n_entries = 1_int32
            call pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, &
                                             pooled_means_perm(seed_pos), sum_abs, n_valid)

            ! Phase 1: the first k_start entries, nearest first, unconditionally. k_start <=
            ! n_valid_means, so a frontier is always left while this loop runs.
            do while (n_entries < k_start)
                entry_pos = nearest_pooled_position(pooled_means, pool_size, pooled_means_perm, n_valid_means, &
                                                    x_star_val, left_pos, right_pos)
                if (entry_pos == left_pos) then
                    left_pos = left_pos - 1_int32
                    first_pos = entry_pos
                else
                    right_pos = right_pos + 1_int32
                    last_pos = entry_pos
                end if
                call pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, &
                                                 pooled_means_perm(entry_pos), entry_sum_abs, entry_n_valid)
                sum_abs = sum_abs + entry_sum_abs
                n_valid = n_valid + entry_n_valid
                n_entries = n_entries + 1_int32
            end do

            if (n_valid < ADAPTIVE_MIN_VALID_RESIDUALS) then
                stop_code = ADAPTIVE_STOP_TOO_FEW_RESIDUALS
                has_too_few_residuals = .true.
            else if (sum_abs == 0.0_real64) then
                stop_code = ADAPTIVE_STOP_ZERO_DISPERSION
            else
                ! Phase 2: rounds of up to k_step entries, each judged against the previous
                ! committed round's dispersion. n_valid >= ADAPTIVE_MIN_VALID_RESIDUALS >= 1 and
                ! sum_abs > 0 here, and a committed round never lowers sum_abs, so s_old > 0.
                s_old = sum_abs/real(n_valid, real64)
                do
                    if (n_entries == k_max) then
                        stop_code = ADAPTIVE_STOP_K_MAX
                        exit
                    end if
                    if (left_pos < 1_int32 .and. right_pos > n_valid_means) then
                        stop_code = ADAPTIVE_STOP_EXHAUSTED
                        exit
                    end if

                    n_planned = min(k_step, k_max - n_entries)
                    stage_left = left_pos
                    stage_right = right_pos
                    stage_first = first_pos
                    stage_last = last_pos
                    round_sum_abs = 0.0_real64
                    round_n_valid = 0_int32
                    n_staged = 0_int32
                    capped = .false.
                    do while (n_staged < n_planned)
                        if (stage_left < 1_int32 .and. stage_right > n_valid_means) exit
                        entry_pos = nearest_pooled_position(pooled_means, pool_size, pooled_means_perm, n_valid_means, &
                                                            x_star_val, stage_left, stage_right)
                        call pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, &
                                                         pooled_means_perm(entry_pos), entry_sum_abs, entry_n_valid)
                        ! Whole entries only: one whose residuals would overflow the cap ends the
                        ! round before it, because nearest-first order may not skip it.
                        if (max_pooled_residuals > 0_int32) then
                            if (n_valid + round_n_valid + entry_n_valid > max_pooled_residuals) then
                                capped = .true.
                                exit
                            end if
                        end if
                        if (entry_pos == stage_left) then
                            stage_left = stage_left - 1_int32
                            stage_first = entry_pos
                        else
                            stage_right = stage_right + 1_int32
                            stage_last = entry_pos
                        end if
                        round_sum_abs = round_sum_abs + entry_sum_abs
                        round_n_valid = round_n_valid + entry_n_valid
                        n_staged = n_staged + 1_int32
                    end do

                    ! A frontier was left at the top of this round, so an empty round can only
                    ! come from the cap.
                    if (n_staged == 0_int32) then
                        stop_code = ADAPTIVE_STOP_RESIDUAL_CAP
                        exit
                    end if

                    ! A round without a non-NaN residual gives exactly s_old back, and commits.
                    s_new = (sum_abs + round_sum_abs)/real(n_valid + round_n_valid, real64)
                    if ((s_new - s_old)/s_old > tau) then
                        stop_code = ADAPTIVE_STOP_TAU
                        exit
                    end if

                    left_pos = stage_left
                    right_pos = stage_right
                    first_pos = stage_first
                    last_pos = stage_last
                    sum_abs = sum_abs + round_sum_abs
                    n_valid = n_valid + round_n_valid
                    n_entries = n_entries + n_staged
                    s_old = s_new
                    if (capped) then
                        stop_code = ADAPTIVE_STOP_RESIDUAL_CAP
                        exit
                    end if
                end do
            end if

            x_star(n_points) = x_star_val
            stop_reason(n_points) = stop_code
            pooled_neighborhood_range(1, n_points) = first_pos
            pooled_neighborhood_range(2, n_points) = last_pos

            ! Decode each pooled position to its study; the studies partition the neighborhood.
            n_neighbors_per_point(:, n_points) = 0_int32
            do i_pos = first_pos, last_pos
                study_idx = (pooled_means_perm(i_pos) - 1_int32)/max_n_genes_all_studies + 1_int32
                n_neighbors_per_point(study_idx, n_points) = n_neighbors_per_point(study_idx, n_points) + 1_int32
            end do
            if (any(n_neighbors_per_point(:, n_points) < min_study_neighbors)) has_empty_study = .true.
            max_n_neighbors = max(max_n_neighbors, maxval(n_neighbors_per_point(:, n_points)))

            if (n_valid > 0_int32) then
                neighborhood_dispersion(n_points) = sum_abs/real(n_valid, real64)
            else
                neighborhood_dispersion(n_points) = M_NAN
            end if
            call calc_sorted_slice_mad(pooled_means, pool_size, pooled_means_perm, first_pos, last_pos, median, mad)
            neighborhood_mad(n_points) = mad

            if (last_pos == n_valid_means) exit

            ! Next seed: the first position after last_pos whose mean is nearest the target, ties
            ! to the lower mean; every branch takes the first position of a run of equal means.
            ! Seeds strictly advance, so there are at most n_valid_means points.
            target = pooled_means(pooled_means_perm(last_pos)) + mad_distance_factor*mad
            target_pos = binary_search_insertion(pooled_means, pooled_means_perm, target, &
                                                 lower_idx=last_pos + 1_int32, upper_idx=n_valid_means)
            if (target_pos > n_valid_means) then
                ! Beyond the largest mean: the first position of the largest mean's run.
                seed_pos = binary_search_insertion(pooled_means, pooled_means_perm, &
                                                   pooled_means(pooled_means_perm(n_valid_means)), &
                                                   lower_idx=last_pos + 1_int32, upper_idx=n_valid_means)
            else if (target_pos == last_pos + 1_int32) then
                seed_pos = target_pos
            else if (pooled_means(pooled_means_perm(target_pos)) - target < &
                     target - pooled_means(pooled_means_perm(target_pos - 1_int32))) then
                seed_pos = target_pos
            else
                ! The lower candidate wins: take the first position of its run of equal means,
                ! the lowest position at that distance.
                seed_pos = binary_search_insertion(pooled_means, pooled_means_perm, &
                                                   pooled_means(pooled_means_perm(target_pos - 1_int32)), &
                                                   lower_idx=last_pos + 1_int32, upper_idx=target_pos - 1_int32)
            end if
        end do

        ! Fixed priority, independent of which neighborhood failed first: too few residuals means
        ! growth itself ran without its dispersion rule, an empty study only affects membership.
        if (has_too_few_residuals) then
            construction_status = ADAPTIVE_STATUS_TOO_FEW_RESIDUALS
        else if (has_empty_study) then
            construction_status = ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD
        else
            construction_status = ADAPTIVE_STATUS_OK
        end if
    end subroutine grow_adaptive_neighborhoods

    !> summary: The unvisited pooled position nearest a reference point, ties to the lower mean
    !| AUTHOR_LASZLO_LANG
    !| Chooses between the two frontiers `left_pos` and `right_pos` of a neighborhood growing
    !| around `x_star_val` (a frontier is gone once `left_pos < 1` or `right_pos > n_valid_means`).
    !| Does not move either frontier; the caller advances the one returned. Precondition, not
    !| checked: at least one frontier is left.
    pure integer(int32) function nearest_pooled_position(pooled_means, pool_size, pooled_means_perm, n_valid_means, &
                                                         x_star_val, left_pos, right_pos) result(entry_pos)
        integer(int32), intent(in) :: pool_size
            !! Number of pooled entries
        real(real64), intent(in) :: pooled_means(pool_size)
            !! Pooled gene means
        integer(int32), intent(in) :: pooled_means_perm(pool_size)
            !! Permutation sorting `pooled_means` ascending, NaN last
        integer(int32), intent(in) :: n_valid_means
            !! Number of non-NaN pooled means
        real(real64), intent(in) :: x_star_val
            !! Reference point the neighborhood grows around
        integer(int32), intent(in) :: left_pos
            !! Next candidate position below the neighborhood
        integer(int32), intent(in) :: right_pos
            !! Next candidate position above the neighborhood

        if (left_pos >= 1_int32 .and. right_pos <= n_valid_means) then
            if (abs(pooled_means(pooled_means_perm(left_pos)) - x_star_val) <= &
                abs(pooled_means(pooled_means_perm(right_pos)) - x_star_val)) then
                entry_pos = left_pos
            else
                entry_pos = right_pos
            end if
        else if (left_pos >= 1_int32) then
            entry_pos = left_pos
        else
            entry_pos = right_pos
        end if
    end function nearest_pooled_position

    !> summary: Sum of absolute values and count of one pooled entry's non-NaN residuals
    !| AUTHOR_LASZLO_LANG
    pure subroutine pooled_entry_residual_stats(max_n_reps_all_studies, pool_size, pooled_residuals, flat_idx, &
                                                entry_sum_abs, entry_n_valid)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
        integer(int32), intent(in) :: pool_size
            !! Number of pooled entries
        real(real64), intent(in) :: pooled_residuals(max_n_reps_all_studies, pool_size)
            !! Residuals, one column per pooled entry, NaN allowed
        integer(int32), intent(in) :: flat_idx
            !! Flat entry whose residuals are summarized
        real(real64), intent(out) :: entry_sum_abs
            !! Sum of the absolute non-NaN residuals
        integer(int32), intent(out) :: entry_n_valid
            !! Number of non-NaN residuals

        integer(int32) :: i_rep
        real(real64) :: residual

        entry_sum_abs = 0.0_real64
        entry_n_valid = 0_int32
        do concurrent(i_rep=1:max_n_reps_all_studies) local(residual) shared(pooled_residuals, flat_idx) &
            reduce(+:entry_sum_abs, entry_n_valid)
            residual = pooled_residuals(i_rep, flat_idx)
            if (.not. ieee_is_nan(residual)) then
                entry_sum_abs = entry_sum_abs + abs(residual)
                entry_n_valid = entry_n_valid + 1_int32
            end if
        end do
    end subroutine pooled_entry_residual_stats

    !> summary: Grow adaptive, dispersion-controlled neighborhoods over all studies' pooled gene means (Issue #217)
    !| AUTHOR_LASZLO_LANG
    !| Instead of one fixed neighbor count, each neighborhood grows while the mean absolute residual
    !| of its pool stays stable, and the next reference point is placed beyond it; the number of
    !| reference points emerges. `k_start`, `k_step` and `k_max` count pooled entries (one entry is
    !| one gene of one study).
    !|
    !| **Pooling.** All studies' gene means are pooled into `N = max_n_genes_all_studies *
    !| n_studies` entries (entry `(s - 1)*G + g` is gene `g` of study `s`) and sorted ascending; NaN
    !| means take no part. A neighborhood is a contiguous run `a..b` of that order
    !| (`pooled_neighborhood_range`, the input of
    !| [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]]),
    !| and its per-study counts `n_neighbors_per_point` add up to `b - a + 1`.
    !|
    !| **Growth** around the reference point `x_star` (the mean at the seed position): entries join
    !| nearest first, ties to the lower mean; the dispersion `S` is the mean absolute non-NaN
    !| residual of the pool.
    !|
    !| 1. The first `k_start` entries join unconditionally. With fewer than
    !|    [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !|    (`CM_ADAPTIVE_MIN_VALID_RESIDUALS`) non-NaN residuals the neighborhood stops with
    !|    [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !|    (`CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS`); with all residuals zero, with
    !|    [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !|    (`CM_ADAPTIVE_STOP_ZERO_DISPERSION`).
    !| 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    !|    `(S_new - S_old) / S_old > tau` is discarded and growth stops
    !|    ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_TAU(variable)]], `CM_ADAPTIVE_STOP_TAU`);
    !|    otherwise it is committed. Growth also stops at `k_max` entries
    !|    ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `CM_ADAPTIVE_STOP_K_MAX`,
    !|    checked first) or when no entry is left
    !|    ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `CM_ADAPTIVE_STOP_EXHAUSTED`).
    !| 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    !|    non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    !|    is still evaluated, and growth then stops with
    !|    [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`CM_ADAPTIVE_STOP_RESIDUAL_CAP`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    !| sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    !| `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    !| absolute deviation of the neighborhood's pooled means. The construction ends with the
    !| neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`CM_ADAPTIVE_STATUS_TOO_FEW_MEANS`) when fewer than `k_start` pooled means are non-NaN: nothing
    !| is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    !| and reports, in this priority,
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS`) if any neighborhood stopped with too few residuals,
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD`) if any has fewer than `min_study_neighbors`
    !| entries of some study, else
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`CM_ADAPTIVE_STATUS_OK`).
    !|
    !| Outputs sized by `N` are filled only for their first `n_points` reference points.
    pure subroutine construct_adaptive_neighborhoods_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                          gene_means, residuals, k_start, k_step, k_max, n_points, &
                                                          x_star, pooled_neighborhood_range, n_neighbors_per_point, &
                                                          stop_reason, neighborhood_dispersion, neighborhood_mad, &
                                                          max_n_neighbors, construction_status, &
                                                          tmp_gene_means_perm_all, tau, mad_distance_factor, &
                                                          max_pooled_residuals, min_study_neighbors)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(in) :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! DM_ALLOW_NAN
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! DM_MIN(k_start)
        integer(int32), intent(out) :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: stop_reason
            !! Why each neighborhood stopped growing: `CM_ADAPTIVE_STOP_TAU` tau exceeded,
            !! `CM_ADAPTIVE_STOP_K_MAX` `k_max` reached, `CM_ADAPTIVE_STOP_RESIDUAL_CAP` residual cap
            !! reached, `CM_ADAPTIVE_STOP_EXHAUSTED` pool exhausted, `CM_ADAPTIVE_STOP_ZERO_DISPERSION`
            !! zero dispersion, `CM_ADAPTIVE_STOP_TOO_FEW_RESIDUALS` too few residuals
            !! DM_RESULT_SIZE_IS(n_points)
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! DM_RESULT_SIZE_IS(n_points)
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! DM_RESULT_SIZE_IS(n_points)
        integer(int32), intent(out) :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(int32), intent(out) :: construction_status
            !! Overall outcome: `CM_ADAPTIVE_STATUS_OK` ok, `CM_ADAPTIVE_STATUS_TOO_FEW_MEANS` too few
            !! means, `CM_ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD` a study neighborhood below
            !! `min_study_neighbors`, `CM_ADAPTIVE_STATUS_TOO_FEW_RESIDUALS` a neighborhood with too
            !! few residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, seeded and sorted here
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)

        real(real64) :: actual_tau, actual_mad_distance_factor
        integer(int32) :: actual_max_pooled_residuals, actual_min_study_neighbors

        M_DEFAULT_VAL(tau, actual_tau, M_ADAPTIVE_TAU_DEFAULT)
        M_DEFAULT_VAL(mad_distance_factor, actual_mad_distance_factor, M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        M_DEFAULT_VAL(max_pooled_residuals, actual_max_pooled_residuals, M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        M_DEFAULT_VAL(min_study_neighbors, actual_min_study_neighbors, M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)

        ! The same pooled sort as the fixed-k search, so an adaptive search can share this order.
        call init_perm(tmp_gene_means_perm_all)
        call sort_real_heapsort_expl_size(gene_means, tmp_gene_means_perm_all, max_n_genes_all_studies*n_studies)

        call grow_adaptive_neighborhoods(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                         residuals, tmp_gene_means_perm_all, k_start, k_step, k_max, actual_tau, &
                                         actual_mad_distance_factor, actual_max_pooled_residuals, &
                                         actual_min_study_neighbors, n_points, x_star, pooled_neighborhood_range, &
                                         n_neighbors_per_point, stop_reason, neighborhood_dispersion, &
                                         neighborhood_mad, max_n_neighbors, construction_status)
    end subroutine construct_adaptive_neighborhoods_impl

    !> summary: Decode one adaptive neighborhood's pooled positions into per-study gene lists
    !| AUTHOR_LASZLO_LANG
    !| For one neighborhood built by
    !| [[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]:
    !| sorted pooled position `pos` names the flat entry `flat = perm_all(pos)`, gene
    !| `mod(flat - 1, G) + 1` of study `(flat - 1)/G + 1`. Each study's genes are appended to its
    !| column of `point_neighborhood_indices` in pooled order and counted in `point_n_neighbors`.
    !| Internal helper, not published.
    !|
    !| Precondition, not checked here: `1 <= first_pos <= last_pos <= max_n_genes_all_studies *
    !| n_studies`, and `max_n_neighbors` is at least the largest per-study count (`G` always is).
    pure subroutine materialize_pooled_neighborhood(perm_all, max_n_genes_all_studies, n_studies, first_pos, last_pos, &
                                                    max_n_neighbors, point_neighborhood_indices, point_n_neighbors)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
        integer(int32), intent(in) :: n_studies
            !! Number of studies
        integer(int32), intent(in) :: perm_all(max_n_genes_all_studies*n_studies)
            !! Permutation sorting the flattened `gene_means(G, S)` ascending, NaN last
        integer(int32), intent(in) :: first_pos
            !! First sorted pooled position of the neighborhood
        integer(int32), intent(in) :: last_pos
            !! Last sorted pooled position of the neighborhood
        integer(int32), intent(in) :: max_n_neighbors
            !! Leading extent of `point_neighborhood_indices`
        integer(int32), intent(out) :: point_neighborhood_indices(max_n_neighbors, n_studies)
            !! Gene indices of the neighborhood, per study, in pooled order; only the first
            !! `point_n_neighbors(i_study)` entries of column `i_study` are defined
        integer(int32), intent(out) :: point_n_neighbors(n_studies)
            !! Number of the neighborhood's entries that belong to each study

        integer(int32) :: i_pos, flat_idx, study_idx

        point_n_neighbors = 0_int32
        ! Sequential: each study's next free slot depends on the entries decoded before it.
        do i_pos = first_pos, last_pos
            flat_idx = perm_all(i_pos) - 1_int32
            study_idx = flat_idx/max_n_genes_all_studies + 1_int32
            point_n_neighbors(study_idx) = point_n_neighbors(study_idx) + 1_int32
            point_neighborhood_indices(point_n_neighbors(study_idx), study_idx) = mod(flat_idx, max_n_genes_all_studies) + 1_int32
        end do
    end subroutine materialize_pooled_neighborhood

    !> summary: Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction
    !| AUTHOR_LASZLO_LANG
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
    pure subroutine check_neighborhood_overlaps_impl(neighborhood_range, n_points, min_neighbor_overlap, &
                                                      all_have_min_neighbor_overlap)
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: neighborhood_range(2, n_points)
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
            !! DM_MIN(1_int32)
        real(real64), intent(in) :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        logical(c_bool), intent(out) :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`

        integer(int32) :: i_point
        real(real64) :: overlap
        logical(c_bool) :: all_pass

        all_pass = .true.
        do concurrent(i_point=1:n_points - 1) local(overlap) shared(neighborhood_range, min_neighbor_overlap) &
                reduce(.and.:all_pass)
            overlap = compute_index_span_overlap(neighborhood_range(1, i_point), neighborhood_range(2, i_point), &
                                                 neighborhood_range(1, i_point + 1), neighborhood_range(2, i_point + 1))
            all_pass = all_pass .and. (overlap >= min_neighbor_overlap)
        end do
        all_have_min_neighbor_overlap = logical(all_pass, kind=c_bool)
    end subroutine check_neighborhood_overlaps_impl

    !> Counts the positions two inclusive integer runs `[a_min, a_max]` and `[b_min, b_max]` share,
    !| relative to the number of positions in the first run:
    !| \[ \frac{\max(0, \min(\texttt{a\_max}, \texttt{b\_max}) - \max(\texttt{a\_min}, \texttt{b\_min}) + 1)}{\texttt{a\_max} - \texttt{a\_min} + 1} \]
    !| For example `[5,6]` and `[6,7]` give 0.5. A reversed run (`max < min`) is empty, so the
    !| result is 0 rather than a division by zero. Private helper of
    !| [[tox_data_integration_preprocessing_impl(module):check_neighborhood_overlaps_impl(interface)]].
    pure real(real64) function compute_index_span_overlap(a_min, a_max, b_min, b_max) result(overlap_fraction)
        integer(int32), intent(in) :: a_min
            !! First position of the first run
        integer(int32), intent(in) :: a_max
            !! Last position of the first run; the run is empty if `a_max < a_min`
        integer(int32), intent(in) :: b_min
            !! First position of the second run
        integer(int32), intent(in) :: b_max
            !! Last position of the second run; the run is empty if `b_max < b_min`

        integer(int32) :: n_shared

        ! An empty first run shares no position
        if (a_max < a_min) then
            overlap_fraction = 0.0_real64
            return
        end if

        n_shared = max(0_int32, min(a_max, b_max) - max(a_min, b_min) + 1_int32)
        overlap_fraction = real(n_shared, real64)/real(a_max - a_min + 1_int32, real64)
    end function compute_index_span_overlap

end module tox_data_integration_preprocessing_impl
