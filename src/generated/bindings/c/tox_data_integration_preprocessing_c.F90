#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_data_integration_preprocessing(module)]]
!| # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing
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
module tox_data_integration_preprocessing_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: compute_gene_means_c
    public :: compute_residuals_c
    public :: pool_means_c
    public :: pool_means_expert_c
    public :: pool_study_means_c
    public :: pool_study_means_expert_c
    public :: construct_neighborhoods_c
    public :: construct_neighborhoods_expert_c
    public :: construct_neighborhoods_ranged_c
    public :: construct_neighborhoods_ranged_expert_c
    public :: construct_adaptive_neighborhoods_c
    public :: construct_adaptive_neighborhoods_expert_c
    public :: check_neighborhood_overlaps_c

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
    subroutine pool_means_c(&
            pooled_means,&
            pool_size,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        ) bind(C, name="pool_means_c")
        use tox_data_integration_preprocessing, only: pool_means

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
            !! Mean-expression reference points. If every pooled mean is NaN, every reference
            !! point is NaN, and `construct_neighborhoods_ranged` rejects such points.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(pool_size)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_pool)
        M_CHECK_ARRAY_NON_NULL(pooled_means, pool_size)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)

        call pool_means(&
            pooled_means = pooled_means,&
            pool_size = pool_size,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_means_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_means_expert(subroutine)]]
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
        use tox_data_integration_preprocessing, only: pool_means_expert

        integer(c_int), intent(in), target :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(c_int), intent(in), target :: n_points
            !! Number of reference points to define
        real(c_double), dimension(pool_size), intent(in), target :: pooled_means
            !! Pooled means
            !! NaN is permitted for this value.
        integer(c_int), dimension(pool_size), intent(in), target :: pooled_means_perm
            !! Sorting permutation for `pooled_means`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `pool_size`.
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points. If every pooled mean is NaN, every reference
            !! point is NaN, and `construct_neighborhoods_ranged` rejects such points.
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

        call pool_means_expert(&
            pooled_means = pooled_means,&
            pooled_means_perm = pooled_means_perm,&
            pool_size = pool_size,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star,&
            ierr = ierr&
        )
    end subroutine pool_means_expert_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_study_means(subroutine)]]
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
        integer(c_int), intent(out), target :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(c_double), dimension(n_points), intent(out), target :: x_star
            !! Mean-expression reference points. If every mean of both studies is NaN, every
            !! reference point is NaN, and `construct_neighborhoods_ranged` rejects such points.
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

        call pool_study_means(&
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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):pool_study_means_expert(subroutine)]]
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
        use tox_data_integration_preprocessing, only: pool_study_means_expert

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
            !! Mean-expression reference points. If every mean of both studies is NaN, every
            !! reference point is NaN, and `construct_neighborhoods_ranged` rejects such points.
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

        call pool_study_means_expert(&
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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_neighborhoods(subroutine)]]
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
            !! [[tox_data_integration_preprocessing_impl(module):calc_neighborhood_size(function)]].
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

        call construct_neighborhoods(&
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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_neighborhoods_expert(subroutine)]]
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
        use tox_data_integration_preprocessing, only: construct_neighborhoods_expert

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
            !! [[tox_data_integration_preprocessing_impl(module):calc_neighborhood_size(function)]].
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

        call construct_neighborhoods_expert(&
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

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_neighborhoods_ranged(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
    !| two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl(subroutine)]]'s
    !| distance-sort algorithm above. In addition to the neighborhood gene indices, this also
    !| reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
    !| genes tied with the span's edge value are never split from it), which a candidate
    !| admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
    !| without needing the gathered residuals.
    subroutine construct_neighborhoods_ranged_c(&
            n_points,&
            x_star,&
            n_genes_S,&
            mean_S,&
            n_neighbors,&
            neighborhood_indices,&
            neighborhood_range,&
            ierr&
        ) bind(C, name="construct_neighborhoods_ranged_c")
        use tox_data_integration_preprocessing, only: construct_neighborhoods_ranged

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors to select per reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_neighbors, n_points), intent(out), target :: neighborhood_indices
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note
            !! All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
            !! `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
            !! ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
            !! limitation, ported as-is from 125-stabilize-jscomp).
            !! @endnote
        integer(c_int), dimension(2, n_points), intent(out), target :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` of the included genes. The
            !! index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
            !! would be the min value. In case of duplicate means, `min_idx` points to the first
            !! appearance of the value and `max_idx` to the last, so even though their related
            !! mean value is the min/max in the neighborhood, the actual gene might not be
            !! included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_genes_S)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_S, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points)

        call construct_neighborhoods_ranged(&
            n_points = n_points,&
            x_star = x_star,&
            n_genes_S = n_genes_S,&
            mean_S = mean_S,&
            n_neighbors = n_neighbors,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            ierr = ierr&
        )
    end subroutine construct_neighborhoods_ranged_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_neighborhoods_ranged_expert(subroutine)]]
    !| Ported from 125-stabilize-jscomp's `construct_neighborhoods_helper`: a binary-search +
    !| two-pointer kNN construction over a pre-sorted `mean_S`, distinct from
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl(subroutine)]]'s
    !| distance-sort algorithm above. In addition to the neighborhood gene indices, this also
    !| reports each reference point's `[min_idx, max_idx]` neighborhood span (tie-extended, so
    !| genes tied with the span's edge value are never split from it), which a candidate
    !| admissibility gate elsewhere in the JSD-Comp-Test parameter search reasons about directly,
    !| without needing the gathered residuals.
    subroutine construct_neighborhoods_ranged_expert_c(&
            n_points,&
            x_star,&
            n_genes_S,&
            mean_S,&
            mean_S_perm,&
            n_neighbors,&
            neighborhood_indices,&
            neighborhood_range,&
            ierr&
        ) bind(C, name="construct_neighborhoods_ranged_expert_c")
        use tox_data_integration_preprocessing, only: construct_neighborhoods_ranged_expert

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), intent(in), target :: n_genes_S
            !! Number of genes in the current study
        integer(c_int), intent(in), target :: n_neighbors
            !! Number of neighbors to select per reference point
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_points), intent(in), target :: x_star
            !! Mean-expression reference points
        real(c_double), dimension(n_genes_S), intent(in), target :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        integer(c_int), dimension(n_genes_S), intent(in), target :: mean_S_perm
            !! Sorting permutation for `mean_S`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes_S`.
        integer(c_int), dimension(n_neighbors, n_points), intent(out), target :: neighborhood_indices
            !! Indices of selected neighborhood genes per reference point.
            !!
            !! @note
            !! All indices are in range `1<=idx<=max(n_neighbors, n_genes_S)`. So in case
            !! `n_genes_S` is lower than `n_neighbors`, remaining indices are filled with the
            !! ones from `n_genes_S+1...n_neighbors` (a documented, deliberately-preserved
            !! limitation, ported as-is from 125-stabilize-jscomp).
            !! @endnote
        integer(c_int), dimension(2, n_points), intent(out), target :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` of the included genes. The
            !! index is related to the permutation vector, so e.g. `mean_S(mean_S_perm(min_idx))`
            !! would be the min value. In case of duplicate means, `min_idx` points to the first
            !! appearance of the value and `max_idx` to the last, so even though their related
            !! mean value is the min/max in the neighborhood, the actual gene might not be
            !! included. If all mean values are NaN, the range is `[1, min(n_genes_S, n_neighbors)]`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(n_genes_S)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(x_star, n_points)
        M_CHECK_ARRAY_NON_NULL(mean_S, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(mean_S_perm, n_genes_S)
        M_CHECK_ARRAY_NON_NULL(neighborhood_indices, n_neighbors * n_points)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points)

        call construct_neighborhoods_ranged_expert(&
            n_points = n_points,&
            x_star = x_star,&
            n_genes_S = n_genes_S,&
            mean_S = mean_S,&
            mean_S_perm = mean_S_perm,&
            n_neighbors = n_neighbors,&
            neighborhood_indices = neighborhood_indices,&
            neighborhood_range = neighborhood_range,&
            ierr = ierr&
        )
    end subroutine construct_neighborhoods_ranged_expert_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_adaptive_neighborhoods(subroutine)]]
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
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !| (`10_int32`) non-NaN residuals the neighborhood stops with
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !| (`6_int32`); with all residuals zero, with
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !| (`5_int32`).
    !| 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    !| `(S_new - S_old) / S_old > tau` is discarded and growth stops
    !| ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_TAU(variable)]], `1_int32`);
    !| otherwise it is committed. Growth also stops at `k_max` entries
    !| ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `2_int32`,
    !| checked first) or when no entry is left
    !| ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `4_int32`).
    !| 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    !| non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    !| is still evaluated, and growth then stops with
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`3_int32`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    !| sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    !| `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    !| absolute deviation of the neighborhood's pooled means. The construction ends with the
    !| neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`1_int32`) when fewer than `k_start` pooled means are non-NaN: nothing
    !| is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    !| and reports, in this priority,
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`3_int32`) if any neighborhood stopped with too few residuals,
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`2_int32`) if any has fewer than `min_study_neighbors`
    !| entries of some study, else
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`0_int32`).
    !|
    !| Outputs sized by `N` are filled only for their first `n_points` reference points.
    subroutine construct_adaptive_neighborhoods_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            k_start,&
            k_step,&
            k_max,&
            n_points,&
            x_star,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
            stop_reason,&
            neighborhood_dispersion,&
            neighborhood_mad,&
            max_n_neighbors,&
            construction_status,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            ierr&
        ) bind(C, name="construct_adaptive_neighborhoods_c")
        use tox_data_integration_preprocessing, only: construct_adaptive_neighborhoods

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        integer(c_int), intent(out), target :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(2, max_n_genes_all_studies*n_studies), intent(out), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: stop_reason
            !! Why each neighborhood stopped growing: `1_int32` tau exceeded,
            !! `2_int32` `k_max` reached, `3_int32` residual cap
            !! reached, `4_int32` pool exhausted, `5_int32`
            !! zero dispersion, `6_int32` too few residuals
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), intent(out), target :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(c_int), intent(out), target :: construction_status
            !! Overall outcome: `0_int32` ok, `1_int32` too few
            !! means, `2_int32` a study neighborhood below
            !! `min_study_neighbors`, `3_int32` a neighborhood with too
            !! few residuals
        real(c_double), intent(in), target :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(c_int), intent(in), target :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(c_int), intent(in), target :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(k_step)
        M_CHECK_NON_NULL(k_max)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(max_n_neighbors)
        M_CHECK_NON_NULL(construction_status)
        M_CHECK_NON_NULL(tau)
        M_CHECK_NON_NULL(mad_distance_factor)
        M_CHECK_NON_NULL(max_pooled_residuals)
        M_CHECK_NON_NULL(min_study_neighbors)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(stop_reason, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_dispersion, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_mad, (max_n_genes_all_studies*n_studies))

        call construct_adaptive_neighborhoods(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            n_points = n_points,&
            x_star = x_star,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
            stop_reason = stop_reason,&
            neighborhood_dispersion = neighborhood_dispersion,&
            neighborhood_mad = neighborhood_mad,&
            max_n_neighbors = max_n_neighbors,&
            construction_status = construction_status,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            ierr = ierr&
        )
    end subroutine construct_adaptive_neighborhoods_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):construct_adaptive_neighborhoods_expert(subroutine)]]
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
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_MIN_VALID_RESIDUALS(variable)]]
    !| (`10_int32`) non-NaN residuals the neighborhood stops with
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_TOO_FEW_RESIDUALS(variable)]]
    !| (`6_int32`); with all residuals zero, with
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_ZERO_DISPERSION(variable)]]
    !| (`5_int32`).
    !| 2. Then rounds of the next `min(k_step, k_max - count)` entries: a round with
    !| `(S_new - S_old) / S_old > tau` is discarded and growth stops
    !| ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_TAU(variable)]], `1_int32`);
    !| otherwise it is committed. Growth also stops at `k_max` entries
    !| ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_K_MAX(variable)]], `2_int32`,
    !| checked first) or when no entry is left
    !| ([[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_EXHAUSTED(variable)]], `4_int32`).
    !| 3. With `max_pooled_residuals > 0`, a round stages whole entries only while the pool's
    !| non-NaN residuals stay within that cap (the `k_start` phase is not capped); a capped round
    !| is still evaluated, and growth then stops with
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STOP_RESIDUAL_CAP(variable)]] (`3_int32`).
    !|
    !| **Reference points.** The first seed is the smallest pooled mean. After a neighborhood ends at
    !| sorted position `b`, the next seed is the first position after `b` whose mean is nearest
    !| `mean(b) + mad_distance_factor * MAD` (ties to the lower mean), `MAD` being the raw median
    !| absolute deviation of the neighborhood's pooled means. The construction ends with the
    !| neighborhood that reaches the largest pooled mean.
    !|
    !| **Construction status.**
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_TOO_FEW_MEANS(variable)]]
    !| (`1_int32`) when fewer than `k_start` pooled means are non-NaN: nothing
    !| is built, `n_points` and `max_n_neighbors` are 0. Otherwise the construction always completes,
    !| and reports, in this priority,
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_TOO_FEW_RESIDUALS(variable)]]
    !| (`3_int32`) if any neighborhood stopped with too few residuals,
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_EMPTY_STUDY_NEIGHBORHOOD(variable)]]
    !| (`2_int32`) if any has fewer than `min_study_neighbors`
    !| entries of some study, else
    !| [[tox_data_integration_preprocessing_impl(module):ADAPTIVE_STATUS_OK(variable)]] (`0_int32`).
    !|
    !| Outputs sized by `N` are filled only for their first `n_points` reference points.
    subroutine construct_adaptive_neighborhoods_expert_c(&
            n_studies,&
            max_n_genes_all_studies,&
            max_n_reps_all_studies,&
            gene_means,&
            residuals,&
            k_start,&
            k_step,&
            k_max,&
            n_points,&
            x_star,&
            pooled_neighborhood_range,&
            n_neighbors_per_point,&
            stop_reason,&
            neighborhood_dispersion,&
            neighborhood_mad,&
            max_n_neighbors,&
            construction_status,&
            tmp_gene_means_perm_all,&
            tau,&
            mad_distance_factor,&
            max_pooled_residuals,&
            min_study_neighbors,&
            ierr&
        ) bind(C, name="construct_adaptive_neighborhoods_expert_c")
        use tox_data_integration_preprocessing, only: construct_adaptive_neighborhoods_expert

        integer(c_int), intent(in), target :: n_studies
            !! Number of studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: max_n_reps_all_studies
            !! Maximum number of replicates across all studies
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(max_n_genes_all_studies, n_studies), intent(in), target :: gene_means
            !! Mean expression of every gene in every study, NaN for a missing gene
            !! NaN is permitted for this value.
        real(c_double), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in), target :: residuals
            !! Signed residuals of every replicate of every gene in every study, NaN for a missing value
            !! NaN is permitted for this value.
        integer(c_int), intent(in), target :: k_start
            !! Pooled entries every neighborhood takes unconditionally
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_step
            !! Pooled entries staged per adaptive growth round
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! The minimum valid value is `k_start`.
        integer(c_int), intent(out), target :: n_points
            !! Number of reference points (neighborhoods) built; 0 when the status is too few means
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: x_star
            !! Reference point of each neighborhood: the pooled mean at its seed position
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(2, max_n_genes_all_studies*n_studies), intent(out), target :: pooled_neighborhood_range
            !! For each neighborhood, its first and last position in the ascending order of the
            !! pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out), target :: n_neighbors_per_point
            !! For each neighborhood, how many of its pooled entries belong to each study
            !! The first `n_points` elements will hold the results.
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: stop_reason
            !! Why each neighborhood stopped growing: `1_int32` tau exceeded,
            !! `2_int32` `k_max` reached, `3_int32` residual cap
            !! reached, `4_int32` pool exhausted, `5_int32`
            !! zero dispersion, `6_int32` too few residuals
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_dispersion
            !! Mean absolute non-NaN residual of each final neighborhood; NaN without any
            !! The first `n_points` elements will hold the results.
        real(c_double), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: neighborhood_mad
            !! Raw median absolute deviation of each final neighborhood's pooled means
            !! The first `n_points` elements will hold the results.
        integer(c_int), intent(out), target :: max_n_neighbors
            !! Largest number of entries one study has in one neighborhood; 0 without a neighborhood
        integer(c_int), intent(out), target :: construction_status
            !! Overall outcome: `0_int32` ok, `1_int32` too few
            !! means, `2_int32` a study neighborhood below
            !! `min_study_neighbors`, `3_int32` a neighborhood with too
            !! few residuals
        integer(c_int), dimension(max_n_genes_all_studies*n_studies), intent(out), target :: tmp_gene_means_perm_all
            !! Working array: sorting permutation of the pooled `gene_means`, seeded and sorted here
        real(c_double), intent(in), target :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `0.1_real64`.
        real(c_double), intent(in), target :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond
            !! the neighborhood's largest mean
            !! The minimum valid value is `0.0_real64`.
            !! The default value is `1.0_real64`.
        integer(c_int), intent(in), target :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! The minimum valid value is `0_int32`.
            !! The default value is `0_int32`.
        integer(c_int), intent(in), target :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have for the construction
            !! status to stay ok
            !! The minimum valid value is `1_int32`.
            !! The default value is `1_int32`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_studies)
        M_CHECK_NON_NULL(max_n_genes_all_studies)
        M_CHECK_NON_NULL(max_n_reps_all_studies)
        M_CHECK_NON_NULL(k_start)
        M_CHECK_NON_NULL(k_step)
        M_CHECK_NON_NULL(k_max)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(max_n_neighbors)
        M_CHECK_NON_NULL(construction_status)
        M_CHECK_NON_NULL(tau)
        M_CHECK_NON_NULL(mad_distance_factor)
        M_CHECK_NON_NULL(max_pooled_residuals)
        M_CHECK_NON_NULL(min_study_neighbors)
        M_CHECK_ARRAY_NON_NULL(gene_means, max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(residuals, max_n_reps_all_studies * max_n_genes_all_studies * n_studies)
        M_CHECK_ARRAY_NON_NULL(x_star, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(pooled_neighborhood_range, 2 * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(n_neighbors_per_point, n_studies * (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(stop_reason, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_dispersion, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(neighborhood_mad, (max_n_genes_all_studies*n_studies))
        M_CHECK_ARRAY_NON_NULL(tmp_gene_means_perm_all, (max_n_genes_all_studies*n_studies))

        call construct_adaptive_neighborhoods_expert(&
            n_studies = n_studies,&
            max_n_genes_all_studies = max_n_genes_all_studies,&
            max_n_reps_all_studies = max_n_reps_all_studies,&
            gene_means = gene_means,&
            residuals = residuals,&
            k_start = k_start,&
            k_step = k_step,&
            k_max = k_max,&
            n_points = n_points,&
            x_star = x_star,&
            pooled_neighborhood_range = pooled_neighborhood_range,&
            n_neighbors_per_point = n_neighbors_per_point,&
            stop_reason = stop_reason,&
            neighborhood_dispersion = neighborhood_dispersion,&
            neighborhood_mad = neighborhood_mad,&
            max_n_neighbors = max_n_neighbors,&
            construction_status = construction_status,&
            tmp_gene_means_perm_all = tmp_gene_means_perm_all,&
            tau = tau,&
            mad_distance_factor = mad_distance_factor,&
            max_pooled_residuals = max_pooled_residuals,&
            min_study_neighbors = min_study_neighbors,&
            ierr = ierr&
        )
    end subroutine construct_adaptive_neighborhoods_expert_c

    !> summary: C-wrapper for [[tox_data_integration_preprocessing(module):check_neighborhood_overlaps(subroutine)]]
    !| The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    !| neighborhood spans of
    !| [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_ranged_impl(interface)]]
    !| (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    !| number of positions they share divided by the length of the earlier span; a reversed span
    !| (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    !| `test_`-prefixed export as a test case.
    subroutine check_neighborhood_overlaps_c(&
            neighborhood_range,&
            n_points,&
            min_neighbor_overlap,&
            all_have_min_neighbor_overlap,&
            ierr&
        ) bind(C, name="check_neighborhood_overlaps_c")
        use tox_data_integration_preprocessing, only: check_neighborhood_overlaps

        integer(c_int), intent(in), target :: n_points
            !! Number of reference points
        integer(c_int), dimension(2, n_points), intent(in), target :: neighborhood_range
            !! For each reference point, the `[min_idx, max_idx]` neighborhood span
            !! The minimum valid value is `1_int32`.
        real(c_double), intent(in), target :: min_neighbor_overlap
            !! Minimum fractional overlap two consecutive neighborhoods must have
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        logical(c_bool), intent(out), target :: all_have_min_neighbor_overlap
            !! `.true.` if every pair of consecutive neighborhoods overlaps by at least
            !! `min_neighbor_overlap`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(min_neighbor_overlap)
        M_CHECK_NON_NULL(all_have_min_neighbor_overlap)
        M_CHECK_ARRAY_NON_NULL(neighborhood_range, 2 * n_points)

        call check_neighborhood_overlaps(&
            neighborhood_range = neighborhood_range,&
            n_points = n_points,&
            min_neighbor_overlap = min_neighbor_overlap,&
            all_have_min_neighbor_overlap = all_have_min_neighbor_overlap,&
            ierr = ierr&
        )
    end subroutine check_neighborhood_overlaps_c

end module tox_data_integration_preprocessing_c
#endif
