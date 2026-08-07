#include <src/macros.h>

!> summary: Wrappers for [[tox_data_integration_preprocessing_impl(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module tox_data_integration_preprocessing
    use tox_data_integration_preprocessing_impl, only: compute_gene_means_impl, compute_residuals_impl, construct_neighborhoods_impl, pool_means_impl
    use tox_data_integration_preprocessing_impl, only: pool_study_means_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_sort, only: init_perm, sort_array_heapsort
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: compute_gene_means
    public :: compute_residuals
    public :: pool_means
    public :: pool_means_expert
    public :: pool_study_means
    public :: pool_study_means_expert
    public :: construct_neighborhoods
    public :: construct_neighborhoods_expert

contains

    !> summary: Validates its inputs, then calls [[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl]].
    pure subroutine compute_gene_means(&
            n_genes,&
            n_reps,&
            expr,&
            means,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes in the study
        integer(int32), intent(in) :: n_reps
            !! Number of biological replicates in the study
        real(real64), dimension(n_reps, n_genes), intent(in) :: expr
            !! Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_genes), intent(out) :: means
            !! Per-gene mean expression values
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_reps, ierr, arg_pos=2_int32)
        if (is_err(ierr)) return
#endif

        call compute_gene_means_impl(&
            n_genes = n_genes,&
            n_reps = n_reps,&
            expr = expr,&
            means = means&
        )
    end subroutine compute_gene_means

    !> summary: Validates its inputs, then calls [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl]].
    pure subroutine compute_residuals(&
            n_genes,&
            n_reps,&
            expr,&
            means,&
            resid,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes
            !! Number of genes in the study
        integer(int32), intent(in) :: n_reps
            !! Number of biological replicates in the study
        real(real64), dimension(n_reps, n_genes), intent(in) :: expr
            !! Expression matrix
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(real64), dimension(n_genes), intent(in) :: means
            !! Per-gene mean expression values; NaN where every replicate of a gene was NaN
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps, n_genes), intent(out) :: resid
            !! Matrix of signed residuals
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_reps, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(means, n_genes, ierr, arg_pos=4_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        call compute_residuals_impl(&
            n_genes = n_genes,&
            n_reps = n_reps,&
            expr = expr,&
            means = means,&
            resid = resid&
        )
    end subroutine compute_residuals

    !> summary: Validates its inputs, prepares what [[tox_data_integration_preprocessing_impl(module):pool_means_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_preprocessing(module):pool_means_expert]] to prepare it yourself.
    !| This takes the pool already built; `pool_study_means` pools the means of two studies
    !| first, if that is what is at hand.
    pure subroutine pool_means(&
            pooled_means,&
            pool_size,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        )
        integer(int32), intent(in) :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        real(real64), dimension(pool_size), intent(in) :: pooled_means
            !! Pooled means
            !! NaN is permitted for this value.
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), dimension(n_points), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: pooled_means_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(pool_size, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(pooled_means, pool_size, ierr, arg_pos=1_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(pooled_means_perm(pool_size))
        call init_perm(pooled_means_perm)
        call sort_array_heapsort(pooled_means, pooled_means_perm)

        call pool_means_impl(&
            pooled_means = pooled_means,&
            pooled_means_perm = pooled_means_perm,&
            pool_size = pool_size,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star&
        )
    end subroutine pool_means

    !> summary: Validates its inputs, then calls [[tox_data_integration_preprocessing_impl(module):pool_means_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_preprocessing(module):pool_means]] does both.
    !| This takes the pool already built; `pool_study_means` pools the means of two studies
    !| first, if that is what is at hand.
    pure subroutine pool_means_expert(&
            pooled_means,&
            pooled_means_perm,&
            pool_size,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        )
        integer(int32), intent(in) :: pool_size
            !! Number of means in the pool, usually `n_genes_S1 + n_genes_S2`
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        real(real64), dimension(pool_size), intent(in) :: pooled_means
            !! Pooled means
            !! NaN is permitted for this value.
        integer(int32), dimension(pool_size), intent(in) :: pooled_means_perm
            !! Sorting permutation for `pooled_means`
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), dimension(n_points), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(pool_size, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_all_in_range_real(pooled_means, pool_size, ierr, arg_pos=1_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        call pool_means_impl(&
            pooled_means = pooled_means,&
            pooled_means_perm = pooled_means_perm,&
            pool_size = pool_size,&
            n_points = n_points,&
            n_pool = n_pool,&
            x_star = x_star&
        )
    end subroutine pool_means_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_preprocessing_impl(module):pool_study_means_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_preprocessing(module):pool_study_means_expert]] to prepare it yourself.
    !| Concatenates the two studies' means, sorts the pool, and turns it into reference
    !| points exactly as `pool_means` does.
    pure subroutine pool_study_means(&
            n_genes_S1,&
            mean_S1,&
            n_genes_S2,&
            mean_S2,&
            n_points,&
            n_pool,&
            x_star,&
            ierr&
        )
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study S1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study S2
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        real(real64), dimension(n_genes_S1), intent(in) :: mean_S1
            !! Per-gene mean expression values of study S1
            !! NaN is permitted for this value.
        real(real64), dimension(n_genes_S2), intent(in) :: mean_S2
            !! Per-gene mean expression values of study S2
            !! NaN is permitted for this value.
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), dimension(n_points), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_pooled_means
        integer(int32), dimension(:), allocatable :: tmp_pooled_means_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes_S1, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S2, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(mean_S1, n_genes_S1, ierr, arg_pos=2_int32, allow_nan=.true.)
        call validate_all_in_range_real(mean_S2, n_genes_S2, ierr, arg_pos=4_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_pooled_means(n_genes_S1+n_genes_S2))
        M_ALLOCATE(tmp_pooled_means_perm(n_genes_S1+n_genes_S2))

        call pool_study_means_impl(&
            n_genes_S1 = n_genes_S1,&
            mean_S1 = mean_S1,&
            n_genes_S2 = n_genes_S2,&
            mean_S2 = mean_S2,&
            n_points = n_points,&
            tmp_pooled_means = tmp_pooled_means,&
            tmp_pooled_means_perm = tmp_pooled_means_perm,&
            n_pool = n_pool,&
            x_star = x_star&
        )
    end subroutine pool_study_means

    !> summary: Validates its inputs, then calls [[tox_data_integration_preprocessing_impl(module):pool_study_means_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_preprocessing(module):pool_study_means]] does both.
    !| Concatenates the two studies' means, sorts the pool, and turns it into reference
    !| points exactly as `pool_means` does.
    pure subroutine pool_study_means_expert(&
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
        )
        integer(int32), intent(in) :: n_genes_S1
            !! Number of genes in study S1
        integer(int32), intent(in) :: n_genes_S2
            !! Number of genes in study S2
        integer(int32), intent(in) :: n_points
            !! Number of reference points to define
        real(real64), dimension(n_genes_S1), intent(in) :: mean_S1
            !! Per-gene mean expression values of study S1
            !! NaN is permitted for this value.
        real(real64), dimension(n_genes_S2), intent(in) :: mean_S2
            !! Per-gene mean expression values of study S2
            !! NaN is permitted for this value.
        real(real64), dimension(n_genes_S1+n_genes_S2), intent(out) :: tmp_pooled_means
            !! Work array holding the concatenated means of both studies
        integer(int32), dimension(n_genes_S1+n_genes_S2), intent(out) :: tmp_pooled_means_perm
            !! Work array for the permutation that sorts `tmp_pooled_means`
        integer(int32), intent(out) :: n_pool
            !! Total number of included (non-NaN) pooled mean-expression values
        real(real64), dimension(n_points), intent(out) :: x_star
            !! Mean-expression reference points
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_genes_S1, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S2, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(mean_S1, n_genes_S1, ierr, arg_pos=2_int32, allow_nan=.true.)
        call validate_all_in_range_real(mean_S2, n_genes_S2, ierr, arg_pos=4_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        call pool_study_means_impl(&
            n_genes_S1 = n_genes_S1,&
            mean_S1 = mean_S1,&
            n_genes_S2 = n_genes_S2,&
            mean_S2 = mean_S2,&
            n_points = n_points,&
            tmp_pooled_means = tmp_pooled_means,&
            tmp_pooled_means_perm = tmp_pooled_means_perm,&
            n_pool = n_pool,&
            x_star = x_star&
        )
    end subroutine pool_study_means_expert

    !> summary: Validates its inputs, prepares what [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl]] needs, then calls it. The entry point to reach for first; see [[tox_data_integration_preprocessing(module):construct_neighborhoods_expert]] to prepare it yourself.
    pure subroutine construct_neighborhoods(&
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
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_reps_S
            !! Number of biological replicates in the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
            !! cannot exceed the number of genes with a defined mean
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
            !! It is recommended to compute this with
            !! [[tox_data_integration_preprocessing_impl(module):calc_neighborhood_size(function)]].
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        real(real64), dimension(n_genes_S), intent(in) :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S, n_genes_S), intent(in) :: resid_S
            !! Matrix of signed residuals
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S, n_neighbors, n_points), intent(out) :: neighborhood_residuals
            !! Collection of residual vectors for each neighborhood
        integer(int32), dimension(n_neighbors, n_points), intent(out) :: neighborhood_indices
            !! Indices of selected neighborhood genes
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        real(real64), dimension(:), allocatable :: tmp_distances
        integer(int32), dimension(:), allocatable :: tmp_distances_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S, ierr, arg_pos=5_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=9_int32, min=1_int32, max=count(.not. ieee_is_nan(mean_S), kind=int32))
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=2_int32, allow_nan=.true.)
        call validate_all_in_range_real(mean_S, n_genes_S, ierr, arg_pos=4_int32, allow_nan=.true.)
        call validate_all_in_range_real(resid_S, n_reps_S * n_genes_S, ierr, arg_pos=6_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_distances(n_genes_S))
        M_ALLOCATE(tmp_distances_perm(n_genes_S))

        call construct_neighborhoods_impl(&
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
            n_neighbors = n_neighbors&
        )
    end subroutine construct_neighborhoods

    !> summary: Validates its inputs, then calls [[tox_data_integration_preprocessing_impl(module):construct_neighborhoods_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[tox_data_integration_preprocessing(module):construct_neighborhoods]] does both.
    pure subroutine construct_neighborhoods_expert(&
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
        )
        integer(int32), intent(in) :: n_points
            !! Number of reference points
        integer(int32), intent(in) :: n_genes_S
            !! Number of genes in the current study
        integer(int32), intent(in) :: n_reps_S
            !! Number of biological replicates in the study
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors; a gene whose mean is NaN can never be a neighbor, so this
            !! cannot exceed the number of genes with a defined mean
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `count(.not. ieee_is_nan(mean_S), kind=int32)`.
            !! It is recommended to compute this with
            !! [[tox_data_integration_preprocessing_impl(module):calc_neighborhood_size(function)]].
        real(real64), dimension(n_points), intent(in) :: x_star
            !! Mean-expression reference points
            !! NaN is permitted for this value.
        real(real64), dimension(n_genes_S), intent(in) :: mean_S
            !! Per-gene mean expression values
            !! NaN is permitted for this value.
        real(real64), dimension(n_reps_S, n_genes_S), intent(in) :: resid_S
            !! Matrix of signed residuals
            !! NaN is permitted for this value.
        real(real64), dimension(n_genes_S), intent(out) :: tmp_distances
            !! Distances work array
        integer(int32), dimension(n_genes_S), intent(out) :: tmp_distances_perm
            !! Work array for the permutation that sorts `tmp_distances`
        real(real64), dimension(n_reps_S, n_neighbors, n_points), intent(out) :: neighborhood_residuals
            !! Collection of residual vectors for each neighborhood
        integer(int32), dimension(n_neighbors, n_points), intent(out) :: neighborhood_indices
            !! Indices of selected neighborhood genes
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_points, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_genes_S, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_reps_S, ierr, arg_pos=5_int32)
        call validate_in_range_int(n_neighbors, ierr, arg_pos=11_int32, min=1_int32, max=count(.not. ieee_is_nan(mean_S), kind=int32))
        call validate_all_in_range_real(x_star, n_points, ierr, arg_pos=2_int32, allow_nan=.true.)
        call validate_all_in_range_real(mean_S, n_genes_S, ierr, arg_pos=4_int32, allow_nan=.true.)
        call validate_all_in_range_real(resid_S, n_reps_S * n_genes_S, ierr, arg_pos=6_int32, allow_nan=.true.)
        if (is_err(ierr)) return
#endif

        call construct_neighborhoods_impl(&
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
            n_neighbors = n_neighbors&
        )
    end subroutine construct_neighborhoods_expert

end module tox_data_integration_preprocessing
