#include "macros.h"

#define CM_OBSERVABLE_COUNT 5
#define CM_SEEDING_COVERAGE_PERCENTILE 0.5_real64
#define CM_DEFAULT_TILE_COUNT 8_int32

#define CM_STOP_FIXED_POINT 1_int32
#define CM_STOP_REJECT_AFTER_ACCEPT 2_int32
#define CM_STOP_NEVER_ACCEPTED 3_int32
#define CM_STOP_NO_CANDIDATES 4_int32

!> Module for managing data structures and geometric calculations.
!! For shatter clustering operations within Tensor Omics.
module tox_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, set_err, set_err_once, is_err, validate_dimension_size, &
                          validate_in_range_real, validate_all_in_range_int, validate_in_range_int, &
                          ERR_ALLOC_FAIL, ERR_DIM_MISMATCH, ERR_INVALID_INPUT, validate_all_in_range_real
    use tox_gene_centroids, only: mean_vector
    use tox_euclidean_distance, only: euclidean_distance_helper
    use f42_utils, only: sort_real_heapsort, calc_percentile, calc_percentile_helper
    use f42_kd_tree, only: vicinity_vectors_helper, vicinity_vectors_count_helper, &
                           KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH
    implicit none
    private

    !> Reason a seed ensemble stopped growing, reported per seed by
    !! [[tox_shatter_cluster_data(module):obtain_ensembles_alloc(subroutine)]].
    integer(int32), parameter :: STC_STOP_FIXED_POINT = CM_STOP_FIXED_POINT
    !! Growth accepted at least one batch, then found no further candidates
    integer(int32), parameter :: STC_STOP_REJECT_AFTER_ACCEPT = CM_STOP_REJECT_AFTER_ACCEPT
    !! A candidate batch was rejected after at least one earlier batch had been accepted
    integer(int32), parameter :: STC_STOP_NEVER_ACCEPTED = CM_STOP_NEVER_ACCEPTED
    !! The first candidate batch was rejected, so the ensemble is still the bare seed
    integer(int32), parameter :: STC_STOP_NO_CANDIDATES = CM_STOP_NO_CANDIDATES
    !! No candidate was ever found, so the seed is isolated at radius `r` and never grew
    public :: calculate_density_radius_alloc, calculate_density_radius, calculate_density_radius_helper
    public :: calculate_labels_as_density_alloc, calculate_labels_as_density, calculate_labels_as_density_helper
    public :: identify_ensemble_seeds_alloc, identify_ensemble_seeds, identify_ensemble_seeds_helper
    public :: grow_ensemble_alloc, grow_ensemble, grow_ensemble_helper
    public :: compute_ambient_density_stats_helper
    public :: compute_ensemble_observable_alloc, compute_ensemble_observable, compute_ensemble_observable_helper
    public :: accept_ensemble, accept_ensemble_helper
    public :: obtain_ensembles_alloc, obtain_ensembles, obtain_ensembles_helper
    public :: grow_single_seed_helper
    public :: merge_ensembles_alloc, merge_ensembles, merge_ensembles_helper
    public :: STC_STOP_FIXED_POINT, STC_STOP_REJECT_AFTER_ACCEPT
    public :: STC_STOP_NEVER_ACCEPTED, STC_STOP_NO_CANDIDATES

contains

    !> Allocating Wrapper for calculating label density search radius.
    subroutine calculate_density_radius_alloc(vectors, n_dimensions, n_vectors, &
                                              radius, mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(out) :: radius
        !! The resulting density search radius
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
        !! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64), allocatable :: tmp_mean_vec(:)
        real(real64), allocatable :: tmp_distances(:)
        integer(int32), allocatable :: tmp_perm(:)

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)

        ! Optional parameter range validation
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))

        call calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, mean_to_other_vecs_dist_quant, ierr)

    end subroutine calculate_density_radius_alloc

    !> Validated Entry Point for calculating label density search radius.
    subroutine calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                        tmp_mean_vec, tmp_distances, tmp_perm, radius, &
                                        mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(out) :: tmp_mean_vec(n_dimensions)
        !! Preallocated workspace for the mean vector
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated workspace for calculated distances
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated workspace for permutation indices
        real(real64), intent(out) :: radius
        !! The resulting density search radius
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
        !! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)

        ! Optional parameter range validation
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        call calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                             tmp_mean_vec, tmp_distances, tmp_perm, &
                                             mean_to_other_vecs_dist_quant, radius, ierr)

    end subroutine calculate_density_radius

    !> Core Implementation for calculating label density search radius.
    pure subroutine calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                                    tmp_mean_vec, tmp_distances, tmp_perm, &
                                                    mean_to_other_vecs_dist_quant, radius, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(out) :: tmp_mean_vec(n_dimensions)
        !! Preallocated workspace for the mean vector
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated workspace for calculated distances
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated workspace for permutation indices
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
        !! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
        real(real64), intent(out) :: radius
        !! The resulting density search radius
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64)   :: actual_quant
        integer(int32) :: i_vec

        call set_ok(ierr)

        M_DEFAULT_VAL(mean_to_other_vecs_dist_quant, actual_quant, 0.15_real64)

        ! Initializing perm
        do concurrent(i_vec=1:n_vectors) shared(tmp_perm)

            tmp_perm(i_vec) = i_vec

        end do

        ! Calculating mean vector
        call mean_vector(vectors, n_dimensions, n_vectors, tmp_perm, n_vectors, tmp_mean_vec, ierr)
        if (is_err(ierr)) return

        ! Calculating each vector distance to mean vector
        do concurrent(i_vec=1:n_vectors) shared(vectors, tmp_mean_vec, n_dimensions, tmp_distances)
            call euclidean_distance_helper(tmp_mean_vec, vectors(:, i_vec), &
                                           n_dimensions, tmp_distances(i_vec))
        end do

        !Sorting perm according to distance
        call sort_real_heapsort(tmp_distances, tmp_perm)

        ! Extract the adaptive search radius corresponding to the specified distance percentile.
        call calc_percentile(tmp_distances, tmp_perm, actual_quant, radius, ierr)

    end subroutine calculate_density_radius_helper

    !> Allocating Wrapper for calculating label density distributions.
    subroutine calculate_labels_as_density_alloc(vectors, n_dimensions, n_vectors, r, &
                                                 dimension_order, kd_indices, label_densities, n_tiles, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
        !! Label-sphere search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output density tracker matching individual vector slots
        integer(int32), intent(in), optional :: n_tiles
        !! Number of concurrent vector tiles, at least 1 (defaults to CM_DEFAULT_TILE_COUNT).
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_stack(:, :, :)
        integer(int32) :: actual_n_tiles

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)

        ! Array and value range validation checks
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_in_range_int(n_tiles, ierr, min=1_int32, arg_pos=8_int32)
        if (is_err(ierr)) return

        M_DEFAULT_VAL(n_tiles, actual_n_tiles, CM_DEFAULT_TILE_COUNT)
        actual_n_tiles = min(actual_n_tiles, n_vectors)
        M_ALLOCATE(tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, actual_n_tiles))

        call calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                         dimension_order, kd_indices, actual_n_tiles, tmp_stack, &
                                         label_densities, ierr)

    end subroutine calculate_labels_as_density_alloc

    !> Validated Entry Point for calculating label density distributions.
    subroutine calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                           dimension_order, kd_indices, n_tiles, tmp_stack, &
                                           label_densities, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
        !! Label-sphere search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(in) :: n_tiles
        !! Number of concurrent vector tiles, at least 1. Excess tiles remain idle.
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, n_tiles)
        !! Preallocated traversal stack per tile, reused serially across its vectors
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output density tracker matching individual vector slots
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)

        ! Array and value range validation checks
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_in_range_int(n_tiles, ierr, min=1_int32, arg_pos=7_int32)
        if (is_err(ierr)) return

        call calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                dimension_order, kd_indices, n_tiles, tmp_stack, &
                                                label_densities, ierr)

    end subroutine calculate_labels_as_density

    !> Core Implementation for calculating label density coordinates.
    pure subroutine calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                       dimension_order, kd_indices, n_tiles, tmp_stack, &
                                                       label_densities, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
        !! Label-sphere search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking tree split axes by variance
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(in) :: n_tiles
        !! Number of concurrent vector tiles, at least 1. Excess tiles remain idle.
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, n_tiles)
        !! Preallocated traversal stack per tile, reused serially across its vectors
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output array storing generated density scalars
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i_tile, i_vec, neighbor_count, vector_start, vector_end, tile_base, tile_rem

        call set_ok(ierr)

        ! Balanced contiguous tiles avoid multiplying n_vectors by a tile index.
        tile_base = n_vectors/n_tiles
        tile_rem = mod(n_vectors, n_tiles)

        do concurrent(i_tile=1:n_tiles) &
            shared(vectors, n_dimensions, n_vectors, r, dimension_order, kd_indices, &
                   tmp_stack, label_densities, tile_base, tile_rem) &
            local(i_vec, neighbor_count, vector_start, vector_end)

            ! Avoid forming n_vectors + 1 for empty tiles at the int32 limit.
            if (tile_base == 0_int32 .and. i_tile > tile_rem) cycle
            vector_start = (i_tile - 1_int32)*tile_base + min(i_tile - 1_int32, tile_rem) + 1_int32
            vector_end = vector_start - 1_int32 + tile_base
            if (i_tile <= tile_rem) vector_end = vector_end + 1_int32

            do i_vec = vector_start, vector_end
                call vicinity_vectors_count_helper(vectors(:, i_vec), vectors, n_dimensions, n_vectors, r, &
                                                   dimension_order, kd_indices, tmp_stack(:, :, i_tile), &
                                                   neighbor_count)
                label_densities(i_vec) = real(neighbor_count, real64)
            end do
        end do

    end subroutine calculate_labels_as_density_helper

    !> Allocating Wrapper for greedy density-ranked coverage-based seed selection.
    subroutine identify_ensemble_seeds_alloc(vectors, n_dimensions, n_vectors, density_labels, &
                                             dimension_order, kd_indices, k_seeding, &
                                             sorted_perm, n_seeds, seed_mask, ierr)

        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for all ambient vectors
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical(c_bool), intent(out) :: seed_mask(n_vectors)
        !! Logical mask identifying selected seed vectors
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_perm(:)
        real(real64), allocatable :: tmp_distances(:)
        integer(int32), allocatable :: tmp_stack(:, :)
        logical(c_bool), allocatable :: tmp_visited_mask(:)
        logical(c_bool), allocatable :: tmp_newly_covered_mask(:)

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        if (is_err(ierr)) return

        call validate_in_range_int(n_vectors, ierr, min=2_int32)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, &
                                       min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, &
                                       min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        call validate_in_range_int(k_seeding, ierr, min=1_int32, &
                                   max=n_vectors - 1_int32)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH))
        M_ALLOCATE(tmp_visited_mask(n_vectors))
        M_ALLOCATE(tmp_newly_covered_mask(n_vectors))

        call identify_ensemble_seeds(vectors, n_dimensions, n_vectors, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)

    end subroutine identify_ensemble_seeds_alloc

    !> Validated Entry Point for greedy density-ranked coverage-based seed selection.
    subroutine identify_ensemble_seeds(vectors, n_dimensions, n_vectors, density_labels, &
                                       dimension_order, kd_indices, k_seeding, &
                                       tmp_perm, tmp_distances, tmp_stack, &
                                       tmp_visited_mask, tmp_newly_covered_mask, &
                                       sorted_perm, n_seeds, seed_mask, ierr)

        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for all ambient vectors
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated permutation workspace
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated candidate-to-vector distance workspace
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated K-D tree traversal stack
        logical(c_bool), intent(out) :: tmp_visited_mask(n_vectors)
        !! Preallocated workspace tracking vectors already covered during seeding
        logical(c_bool), intent(out) :: tmp_newly_covered_mask(n_vectors)
        !! Preallocated workspace for the current seed's coverage query
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical(c_bool), intent(out) :: seed_mask(n_vectors)
        !! Logical mask identifying selected seed vectors
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        if (is_err(ierr)) return

        call validate_in_range_int(n_vectors, ierr, min=2_int32)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, &
                                       min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, &
                                       min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        call validate_in_range_int(k_seeding, ierr, min=1_int32, &
                                   max=n_vectors - 1_int32)
        if (is_err(ierr)) return

        call identify_ensemble_seeds_helper(vectors, n_dimensions, n_vectors, density_labels, &
                                            dimension_order, kd_indices, k_seeding, &
                                            tmp_perm, tmp_distances, tmp_stack, &
                                            tmp_visited_mask, tmp_newly_covered_mask, &
                                            sorted_perm, n_seeds, seed_mask)

    end subroutine identify_ensemble_seeds

    !> Core Implementation for greedy density-ranked coverage-based seed selection.
    pure subroutine identify_ensemble_seeds_helper(vectors, n_dimensions, n_vectors, density_labels, &
                                                   dimension_order, kd_indices, k_seeding, &
                                                   tmp_perm, tmp_distances, tmp_stack, &
                                                   tmp_visited_mask, tmp_newly_covered_mask, &
                                                   sorted_perm, n_seeds, seed_mask)

        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for all ambient vectors
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated permutation workspace
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated candidate-to-vector distance workspace
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated K-D tree traversal stack
        logical(c_bool), intent(out) :: tmp_visited_mask(n_vectors)
        !! Preallocated workspace tracking vectors already covered during seeding
        logical(c_bool), intent(out) :: tmp_newly_covered_mask(n_vectors)
        !! Preallocated workspace for the current seed's coverage query
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical(c_bool), intent(out) :: seed_mask(n_vectors)
        !! Logical mask identifying selected seed vectors

        integer(int32) :: i_vec, i_rank, candidate_idx, self_pos, median_pos
        real(real64) :: coverage_radius

        do concurrent(i_vec=1:n_vectors) shared(tmp_perm)
            tmp_perm(i_vec) = i_vec
        end do

        call sort_real_heapsort(density_labels, tmp_perm)

        ! Heapsort provides ascending density order. Reverse the permutation.
        do concurrent(i_vec=1:n_vectors) shared(sorted_perm, tmp_perm, n_vectors)
            sorted_perm(i_vec) = tmp_perm(n_vectors - i_vec + 1_int32)
        end do

        seed_mask = .false.
        tmp_visited_mask = .false.
        n_seeds = 0_int32

        do i_rank = 1, n_vectors
            candidate_idx = sorted_perm(i_rank)

            if (tmp_visited_mask(candidate_idx)) cycle

            ! Highest-density unvisited vector becomes the next seed.
            seed_mask(candidate_idx) = .true.
            tmp_visited_mask(candidate_idx) = .true.
            n_seeds = n_seeds + 1_int32

            ! Calculate the distance from the current seed to every ambient vector.
            do concurrent(i_vec=1:n_vectors) &
                shared(vectors, n_dimensions, candidate_idx, tmp_distances)

                call euclidean_distance_helper(vectors(:, candidate_idx), vectors(:, i_vec), &
                                               n_dimensions, tmp_distances(i_vec))
            end do

            do concurrent(i_vec=1:n_vectors) shared(tmp_perm)
                tmp_perm(i_vec) = i_vec
            end do

            call sort_real_heapsort(tmp_distances, tmp_perm)

            ! Exclude the seed itself from its k_seeding nearest neighbors.
            self_pos = 0_int32

            do i_vec = 1, n_vectors
                if (tmp_perm(i_vec) == candidate_idx) then
                    self_pos = i_vec
                    exit
                end if
            end do

            if ((self_pos >= 1_int32) .and. (self_pos <= k_seeding)) then
                do i_vec = self_pos, k_seeding
                    tmp_perm(i_vec) = tmp_perm(i_vec + 1_int32)
                end do
            end if

            ! Uses the 50th percentile of the k nearest-neighbor distances
            ! as the seed coverage radius.
            if (mod(k_seeding, 2_int32) == 1_int32) then
                median_pos = (k_seeding + 1_int32)/2_int32
                coverage_radius = tmp_distances(tmp_perm(median_pos))
            else
                median_pos = k_seeding/2_int32
                coverage_radius = 0.5_real64* &
                                  (tmp_distances(tmp_perm(median_pos)) + &
                                   tmp_distances(tmp_perm(median_pos + 1_int32)))
            end if

            ! Mark the region represented by the current seed as visited.
            call vicinity_vectors_helper(vectors(:, candidate_idx), vectors, &
                                         n_dimensions, n_vectors, coverage_radius, &
                                         dimension_order, kd_indices, tmp_stack, &
                                         tmp_newly_covered_mask)

            tmp_visited_mask = tmp_visited_mask .or. tmp_newly_covered_mask
        end do

    end subroutine identify_ensemble_seeds_helper

    !> Core Implementation for the ambient density median and MAD.
    pure subroutine compute_ambient_density_stats_helper(density_labels, n_vectors, &
                                                         tmp_perm, tmp_abs_diff, &
                                                         median_ambient, mad_ambient)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array for the absolute deviations
        real(real64), intent(out) :: median_ambient
        !! Ambient median density median_{i=1..N}(rho_i)
        real(real64), intent(out) :: mad_ambient
        !! Ambient median absolute deviation median_{i=1..N}(|rho_i - median_ambient|)

        integer(int32) :: i_vec

        ! 1. Calculate ambient median density: median_{i=1..N}(rho_i)
        do concurrent(i_vec=1:n_vectors) shared(tmp_perm)
            tmp_perm(i_vec) = i_vec
        end do
        call sort_real_heapsort(density_labels, tmp_perm)
        call calc_percentile_helper(density_labels, tmp_perm, 0.5_real64, median_ambient)

        ! 2. Compute absolute deviation from ambient median for all vectors
        do concurrent(i_vec=1:n_vectors) shared(tmp_abs_diff, tmp_perm, density_labels, median_ambient)
            tmp_abs_diff(i_vec) = abs(density_labels(i_vec) - median_ambient)
            tmp_perm(i_vec) = i_vec
        end do

        ! 3. Compute MAD_ambient = median_{i=1..N}(|rho_i - median_ambient|)
        call sort_real_heapsort(tmp_abs_diff, tmp_perm)
        call calc_percentile_helper(tmp_abs_diff, tmp_perm, 0.5_real64, mad_ambient)

    end subroutine compute_ambient_density_stats_helper

    !> Allocating Wrapper for growing the active ensemble at its surface.
    subroutine grow_ensemble_alloc(vectors, n_dimensions, n_vectors, &
                                   ensemble_mask, dimension_order, kd_indices, &
                                   density_labels, r, alpha_mad, mad_ambient, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), intent(inout) :: ensemble_mask(n_vectors)
        !! Logical mask tracking active ensemble member vectors
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        real(real64), intent(in), optional :: r
        !! Optional search radius; defaults to computed density label sphere radius
        real(real64), intent(in), optional :: alpha_mad
        !! Optional multiplier factor for MAD density compatibility threshold (defaults to 0.5)
        real(real64), intent(out), optional :: mad_ambient
        !! Optional ambient MAD used to scale alpha_mad; always positive on success because zero MAD is rejected before growth.
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64) :: tmp_mad_ambient
        integer(int32), allocatable :: tmp_stack(:, :)
        logical(c_bool), allocatable :: tmp_vicinity_mask(:)
        logical(c_bool), allocatable :: tmp_surface_mask(:)
        integer(int32), allocatable :: tmp_perm(:)
        real(real64), allocatable :: tmp_abs_diff(:)
        real(real64) :: actual_r, actual_alpha_mad

        ! Workspaces for default radius calculation (if r is omitted)
        real(real64), allocatable :: tmp_mean_vec(:)
        real(real64), allocatable :: tmp_distances(:)

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)

        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)

        ! Validate optional parameter alpha_mad directly (validate_in_range_real safely checks if present)
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        ! Assign default value after validation
        M_DEFAULT_VAL(alpha_mad, actual_alpha_mad, 0.5_real64)

        ! Determine actual search radius threshold r
        if (present(r)) then
            actual_r = r
            call validate_in_range_real(actual_r, ierr, min=0.0_real64)
            if (is_err(ierr)) return
        else
            ! Dynamically compute the default density label sphere radius
            M_ALLOCATE(tmp_mean_vec(n_dimensions))
            M_ALLOCATE(tmp_distances(n_vectors))
            M_ALLOCATE(tmp_perm(n_vectors))

            call calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                          tmp_mean_vec, tmp_distances, tmp_perm, &
                                          actual_r, ierr=ierr)

            deallocate (tmp_mean_vec, tmp_distances, tmp_perm)
            if (is_err(ierr)) return
        end if

        ! Allocate preallocated workspaces
        M_ALLOCATE(tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors))
        M_ALLOCATE(tmp_surface_mask(n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_abs_diff(n_vectors))

        call grow_ensemble(vectors, n_dimensions, n_vectors, ensemble_mask, &
                           actual_r, dimension_order, kd_indices, density_labels, &
                           actual_alpha_mad, tmp_stack, tmp_vicinity_mask, &
                           tmp_surface_mask, tmp_perm, tmp_abs_diff, tmp_mad_ambient, ierr)

        if (is_err(ierr)) return

        if (present(mad_ambient)) mad_ambient = tmp_mad_ambient

    end subroutine grow_ensemble_alloc

    !> Validated Entry Point for surface-based ensemble growth.
    subroutine grow_ensemble(vectors, n_dimensions, n_vectors, ensemble_mask, &
                             r, dimension_order, kd_indices, density_labels, &
                             alpha_mad, tmp_stack, tmp_vicinity_mask, &
                             tmp_surface_mask, tmp_perm, tmp_abs_diff, mad_ambient, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), intent(inout) :: ensemble_mask(n_vectors)
        !! Logical mask tracking active ensemble member vectors
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        real(real64), intent(in) :: alpha_mad
        !! Multiplier factor for MAD density compatibility threshold
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated workspace stack for a single tree traversal
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors)
        !! Preallocated workspace mask holding one member's neighborhood query result
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors)
        !! Preallocated workspace mask accumulating the union of all member neighborhoods
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array for MAD and member density calculations
        real(real64), intent(out) :: mad_ambient
        !! Optional ambient MAD used to scale alpha_mad; always positive on success because zero MAD is rejected before growth.
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64) :: median_ambient

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)

        ! Array and value range validation checks
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        ! Standalone single-step growth still derives the ambient statistics itself.
        call compute_ambient_density_stats_helper(density_labels, n_vectors, &
                                                  tmp_perm, tmp_abs_diff, &
                                                  median_ambient, mad_ambient)

        ! Reject zero ambient MAD, which makes alpha_mad ineffective by reducing growth to exact-density matching.
        if (mad_ambient <= 0.0_real64) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call grow_ensemble_helper(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                  r, dimension_order, kd_indices, density_labels, &
                                  alpha_mad, mad_ambient, tmp_stack, tmp_vicinity_mask, &
                                  tmp_surface_mask, tmp_perm, tmp_abs_diff)

    end subroutine grow_ensemble

    !> Core Implementation for surface-based ensemble growth.
    pure subroutine grow_ensemble_helper(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                         r, dimension_order, kd_indices, density_labels, &
                                         alpha_mad, mad_ambient, tmp_stack, tmp_vicinity_mask, &
                                         tmp_surface_mask, tmp_perm, tmp_abs_diff)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), intent(inout) :: ensemble_mask(n_vectors)
        !! Logical mask tracking active ensemble member vectors
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        real(real64), intent(in) :: alpha_mad
        !! Multiplier factor for MAD density compatibility threshold
        real(real64), intent(in) :: mad_ambient
        !! Ambient MAD precalculated via [[tox_shatter_cluster_data(module):compute_ambient_density_stats_helper(subroutine)]].
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated workspace stack for a single tree traversal
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors)
        !! Workspace for one neighborhood query, then the complete trial membership
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors)
        !! Geometric union of accepted member neighborhoods before this growth step
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array holding the active member densities

        integer(int32) :: i_vec

        call propose_ensemble_growth_helper(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                            r, dimension_order, kd_indices, density_labels, &
                                            alpha_mad, mad_ambient, tmp_stack, tmp_vicinity_mask, &
                                            tmp_surface_mask, tmp_perm, tmp_abs_diff)

        ! The standalone growth API commits its single step unconditionally.
        do concurrent(i_vec=1:n_vectors) shared(ensemble_mask, tmp_vicinity_mask)
            if (tmp_vicinity_mask(i_vec) .and. .not. ensemble_mask(i_vec)) ensemble_mask(i_vec) = .true.
        end do

    end subroutine grow_ensemble_helper

    !> Discovers a standalone ensemble surface and builds trial membership in vicinity scratch.
    pure subroutine propose_ensemble_growth_helper(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                                   r, dimension_order, kd_indices, density_labels, &
                                                   alpha_mad, mad_ambient, tmp_stack, tmp_vicinity_mask, &
                                                   tmp_surface_mask, tmp_perm, tmp_abs_diff)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Accepted ensemble membership, unchanged while constructing the trial
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        real(real64), intent(in) :: alpha_mad
        !! Multiplier factor for MAD density compatibility threshold
        real(real64), intent(in) :: mad_ambient
        !! Ambient MAD precalculated via [[tox_shatter_cluster_data(module):compute_ambient_density_stats_helper(subroutine)]].
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated workspace stack for a single tree traversal
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors)
        !! Workspace for one neighborhood query, then the complete trial membership
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors)
        !! Geometric union of accepted member neighborhoods, independent of density compatibility
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array holding the active member densities

        integer(int32) :: i_vec

        tmp_surface_mask = .false.
        do i_vec = 1, n_vectors
            if (ensemble_mask(i_vec)) then
                call vicinity_vectors_helper(vectors(:, i_vec), vectors, n_dimensions, n_vectors, r, &
                                             dimension_order, kd_indices, tmp_stack, &
                                             tmp_vicinity_mask)

                tmp_surface_mask = tmp_surface_mask .or. tmp_vicinity_mask
            end if
        end do

        call build_ensemble_trial_helper(ensemble_mask, n_vectors, density_labels, alpha_mad, &
                                         mad_ambient, tmp_surface_mask, tmp_perm, tmp_abs_diff, tmp_vicinity_mask)

    end subroutine propose_ensemble_growth_helper

    !> Reevaluate density compatibility on a geometric surface using the current accepted median.
    pure subroutine build_ensemble_trial_helper(ensemble_mask, n_vectors, density_labels, alpha_mad, &
                                                mad_ambient, surface_mask, tmp_perm, tmp_abs_diff, trial_mask)
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Accepted membership, always retained in the trial
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors
        real(real64), intent(in) :: alpha_mad
        !! Multiplier for the ambient MAD density tolerance
        real(real64), intent(in) :: mad_ambient
        !! Precomputed ambient median absolute deviation
        logical(c_bool), intent(in) :: surface_mask(n_vectors)
        !! Spatial neighborhood union; never filtered by density compatibility
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Sorting workspace; any frontier indices must be consumed before this call
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Workspace for accepted member densities
        logical(c_bool), intent(out) :: trial_mask(n_vectors)
        !! Complete trial membership: accepted members and density-compatible spatial neighbors

        integer(int32) :: i_vec, k, n_active
        real(real64) :: ensemble_center_density, max_allowed_dev

        max_allowed_dev = alpha_mad*mad_ambient
        n_active = count(ensemble_mask)
        if (n_active == 0_int32) then
            trial_mask = .false.
            return
        end if

        k = 0
        do i_vec = 1, n_vectors
            if (ensemble_mask(i_vec)) then
                k = k + 1
                tmp_abs_diff(k) = density_labels(i_vec)
                tmp_perm(k) = k
            end if
        end do
        call sort_real_heapsort(tmp_abs_diff(1:k), tmp_perm(1:k))
        call calc_percentile_helper(tmp_abs_diff(1:k), tmp_perm(1:k), 0.5_real64, &
                                    ensemble_center_density)

        do concurrent(i_vec=1:n_vectors) &
            shared(ensemble_mask, surface_mask, trial_mask, density_labels, ensemble_center_density, max_allowed_dev)
            trial_mask(i_vec) = ensemble_mask(i_vec) .or. &
                                (surface_mask(i_vec) .and. &
                                 abs(density_labels(i_vec) - ensemble_center_density) <= max_allowed_dev)
        end do

    end subroutine build_ensemble_trial_helper

    !> Allocating Wrapper for calculating and storing ensemble observable trajectories.
    subroutine compute_ensemble_observable_alloc(ensemble_mask, density_labels, n_vectors, &
                                                 candidate_count, current_iter, observables, &
                                                 t_observables, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Selection mask for active ensemble members
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: candidate_count
        !! Number of candidates added in current growth step
        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(inout), allocatable :: observables(:, :)
        !! 2D observable matrix history [5 x window_capacity] allocated on demand
        integer(int32), intent(in), optional :: t_observables
        !! History tracking depth (defaults to 10; <= 0 represents infinity / all)
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: actual_t_obs, required_cols

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_vectors, ierr)
        call validate_dimension_size(current_iter, ierr)
        call validate_in_range_int(candidate_count, ierr, min=0_int32)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)
        if (is_err(ierr)) return

        M_DEFAULT_VAL(t_observables, actual_t_obs, 10_int32)

        if (actual_t_obs <= 0_int32) then
            required_cols = n_vectors
        else
            required_cols = actual_t_obs
        end if

        if (.not. allocated(observables)) then
            M_ALLOCATE(observables(CM_OBSERVABLE_COUNT, required_cols))
            observables = 0.0_real64
        end if

        ! Checked after the on-demand allocation, so only a caller-supplied matrix can fail here
        if (size(observables, dim=1, kind=int32) /= CM_OBSERVABLE_COUNT) call set_err_once(ierr, ERR_DIM_MISMATCH)
        if (is_err(ierr)) return

        call compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                candidate_count, current_iter, observables, &
                                                t_observables)

    end subroutine compute_ensemble_observable_alloc

    !> Validated Entry Point for computing and updating 2D ensemble density observables.
    subroutine compute_ensemble_observable(ensemble_mask, density_labels, n_vectors, &
                                           candidate_count, current_iter, observables, &
                                           t_observables, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Selection mask for active ensemble members
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: candidate_count
        !! Number of candidates added in current growth step
        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(inout) :: observables(:, :)
        !! 2D observable trajectory matrix where columns match iteration indices
        integer(int32), intent(in), optional :: t_observables
        !! History tracking depth (defaults to 10; <= 0 represents infinity / all)
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_vectors, ierr)
        call validate_dimension_size(current_iter, ierr)
        call validate_in_range_int(candidate_count, ierr, min=0_int32)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)
        if (size(observables, dim=1, kind=int32) /= CM_OBSERVABLE_COUNT) call set_err_once(ierr, ERR_DIM_MISMATCH)
        if (is_err(ierr)) return

        call compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                candidate_count, current_iter, observables, &
                                                t_observables)

    end subroutine compute_ensemble_observable

    !> Core Implementation for calculating 5-component observables into a 2D history matrix.
    pure subroutine compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                       candidate_count, current_iter, observables, &
                                                       t_observables)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical(c_bool), intent(in) :: ensemble_mask(n_vectors)
        !! Selection mask for active ensemble members
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: candidate_count
        !! Number of candidates added in current growth step
        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(inout) :: observables(:, :)
        !! 2D observable trajectory matrix where columns match iteration indices
        integer(int32), intent(in), optional :: t_observables
        !! History tracking depth (defaults to 10; <= 0 represents infinity / all)

        integer(int32) :: i_vec, active_count, actual_t_obs, max_cols, zero_count
        real(real64)   :: sum_arithmetic, sum_reciprocal, rho_arith, rho_harm
        real(real64)   :: obs_1, obs_2, obs_3, obs_4, obs_5

        max_cols = size(observables, dim=2, kind=int32)

        M_DEFAULT_VAL(t_observables, actual_t_obs, 10_int32)
        if (actual_t_obs <= 0_int32) then
            actual_t_obs = max_cols
        else
            actual_t_obs = min(actual_t_obs, max_cols)
        end if

        active_count = count(ensemble_mask)

        if (active_count == 0_int32) then
            obs_1 = 0.0_real64
            obs_2 = 0.0_real64
            obs_3 = 0.0_real64
            obs_4 = 0.0_real64
        else
            sum_arithmetic = 0.0_real64
            sum_reciprocal = 0.0_real64
            zero_count = 0_int32

            do concurrent(i_vec=1:n_vectors) &
                shared(ensemble_mask, density_labels) &
                reduce(+:sum_arithmetic, sum_reciprocal, zero_count)

                if (ensemble_mask(i_vec)) then
                    sum_arithmetic = sum_arithmetic + density_labels(i_vec)

                    if (density_labels(i_vec) > 0.0_real64) then
                        sum_reciprocal = sum_reciprocal + (1.0_real64/density_labels(i_vec))
                    else
                        zero_count = zero_count + 1_int32
                    end if
                end if
            end do

            ! 1. Arithmetic Mean Density
            rho_arith = sum_arithmetic/real(active_count, real64)

            ! 2. Harmonic Mean Density (strictly 0 if any active element is <= 0).
            if (zero_count == 0_int32 .and. sum_reciprocal > 0.0_real64) then
                rho_harm = real(active_count, real64)/sum_reciprocal
            else
                rho_harm = 0.0_real64
            end if

            obs_1 = rho_arith
            obs_2 = rho_harm

            ! 3. Density Heterogeneity H_rho.
            if (rho_harm > 0.0_real64) then
                obs_3 = rho_arith/rho_harm
            else
                obs_3 = 1.0_real64
            end if

            ! 4. Active Ensemble Size
            obs_4 = real(active_count, real64)
        end if

        ! 5. Candidate Count (recorded independently of active_count)
        obs_5 = real(candidate_count, real64)

        ! Store scalar observables matching iteration index and history window
        if (current_iter <= actual_t_obs) then
            observables(1, current_iter) = obs_1
            observables(2, current_iter) = obs_2
            observables(3, current_iter) = obs_3
            observables(4, current_iter) = obs_4
            observables(5, current_iter) = obs_5
        else
            observables(:, 1:actual_t_obs - 1_int32) = observables(:, 2:actual_t_obs)
            observables(1, actual_t_obs) = obs_1
            observables(2, actual_t_obs) = obs_2
            observables(3, actual_t_obs) = obs_3
            observables(4, actual_t_obs) = obs_4
            observables(5, actual_t_obs) = obs_5
        end if

    end subroutine compute_ensemble_observable_helper

    !> Validated Entry Point for decision rule evaluating step acceptance.
    subroutine accept_ensemble(observables, current_iter, alpha_accept, &
                               is_accepted, ierr)

        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(in) :: observables(:, :)
        !! 2D observables trajectory matrix [5 x window_capacity]
        real(real64), intent(in), optional :: alpha_accept
        !! Optional log2 fold change acceptance threshold (defaults to 0.5)
        logical(c_bool), intent(out) :: is_accepted
        !! Output flag indicating step acceptance decision
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate input structural dimension
        call validate_dimension_size(current_iter, ierr)

        ! Optional parameter range validation
        call validate_in_range_real(alpha_accept, ierr, min=0.0_real64)
        if (size(observables, dim=1, kind=int32) /= CM_OBSERVABLE_COUNT) call set_err_once(ierr, ERR_DIM_MISMATCH)
        if (is_err(ierr)) return

        call accept_ensemble_helper(observables, current_iter, alpha_accept, is_accepted)

    end subroutine accept_ensemble

    !> Core Implementation for decision rule evaluating step acceptance.
    pure subroutine accept_ensemble_helper(observables, current_iter, alpha_accept, &
                                           is_accepted)

        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(in) :: observables(:, :)
        !! 2D observables trajectory matrix [5 x window_capacity]
        real(real64), intent(in), optional :: alpha_accept
        !! Optional log2 fold change acceptance threshold (defaults to 0.5)
        logical(c_bool), intent(out) :: is_accepted
        !! Output flag indicating step acceptance decision

        real(real64)   :: actual_alpha, h_prev, h_curr, log2_fold_change
        integer(int32) :: max_cols, col_curr, col_prev

        M_DEFAULT_VAL(alpha_accept, actual_alpha, 0.5_real64)

        ! Initial seed iteration (iter == 1) is always accepted
        if (current_iter <= 1_int32) then
            is_accepted = .true.
            return
        end if

        max_cols = size(observables, dim=2, kind=int32)

        ! Determine column indices for current and previous iterations
        if (current_iter <= max_cols) then
            col_curr = current_iter
            col_prev = current_iter - 1_int32
        else
            ! Sliding window shifted: last column is current, second-to-last is previous
            col_curr = max_cols
            col_prev = max_cols - 1_int32
        end if

        ! Safety check for minimal matrix width
        if (col_prev < 1_int32) then
            is_accepted = .true.
            return
        end if

        ! Evaluate log2 fold change of Density Heterogeneity H_rho (Row 3)
        h_prev = observables(3, col_prev)
        h_curr = observables(3, col_curr)

        if ((h_prev > 0.0_real64) .and. (h_curr > 0.0_real64)) then
            log2_fold_change = log(h_curr/h_prev)/log(2.0_real64)

            ! Reject if fold change meets or exceeds acceptance threshold alpha_accept
            if (abs(log2_fold_change) >= actual_alpha) then
                is_accepted = .false.
            else
                is_accepted = .true.
            end if
        else
            is_accepted = .true.
        end if

    end subroutine accept_ensemble_helper

    !> Allocating Wrapper for multi-ensemble parallel extraction.
    subroutine obtain_ensembles_alloc(vectors, n_dimensions, n_vectors, &
                                      dimension_order, kd_indices, density_labels, &
                                      seed_indices, n_seeds, r, alpha_mad, &
                                      alpha_accept, t_observables, n_tiles, &
                                      ensemble_matrix, stop_reasons, mad_ambient, &
                                      n_ensembles, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for ambient vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: n_seeds
        !! Total number of precalculated seed vectors [[tox_shatter_cluster_data(module):identify_ensemble_seeds_alloc(subroutine)]].
        integer(int32), intent(in) :: seed_indices(n_seeds)
        !! Array of starting seed vector indices [[tox_shatter_cluster_data(module):identify_ensemble_seeds_alloc(subroutine)]].
        real(real64), intent(in), optional :: r
        !! Search radius threshold for surface growth
        real(real64), intent(in), optional :: alpha_mad
        !! MAD multiplier for density compatibility (defaults to 0.5)
        real(real64), intent(in), optional :: alpha_accept
        !! Acceptance log2 fold-change threshold (defaults to 0.5)
        integer(int32), intent(in), optional :: t_observables
        !! Observable history depth window (defaults to 10)
        integer(int32), intent(in), optional :: n_tiles
        !! Number of concurrently processed seed tiles, `n_tiles >= 1` (defaults to 8).
        logical(c_bool), intent(out) :: ensemble_matrix(n_vectors, n_seeds)
        !! Preallocated output matrix [n_vectors x n_seeds] tracking grown raw ensembles
        integer(int32), intent(out), optional :: stop_reasons(n_seeds)
        !! Optional per-seed growth termination reason, one of the `STC_STOP_*` constants
        real(real64), intent(out), optional :: mad_ambient
        !! Optional ambient MAD of the density labels, i.e. the value `alpha_mad` was scaled by.
        !! Always `> 0` on success: a zero MAD makes the threshold degenerate for every
        !! `alpha_mad` and is rejected before growth.
        integer(int32), intent(out) :: n_ensembles
        !! Total count of extracted raw ensembles (equals n_seeds)
        integer(int32), intent(out) :: ierr
        !! Error status flag

        real(real64) :: actual_r, actual_alpha_mad, actual_alpha_accept
        integer(int32) :: actual_t_obs, actual_n_tiles, required_cols

        ! Workspaces
        integer(int32), allocatable :: seed_perm(:)
        integer(int32), allocatable :: tmp_stack(:, :, :)
        logical(c_bool), allocatable :: tmp_vicinity_mask(:, :)
        logical(c_bool), allocatable :: tmp_surface_mask(:, :)
        integer(int32), allocatable :: tmp_perm(:, :)
        real(real64), allocatable :: tmp_abs_diff(:, :)
        real(real64), allocatable :: tmp_observables(:, :, :)
        logical(c_bool), allocatable :: tmp_current_mask(:, :)
        real(real64), allocatable :: tmp_mean_vec(:), tmp_distances(:)
        integer(int32), allocatable :: tmp_stop_reasons(:)
        real(real64) :: tmp_mad_ambient

        call set_ok(ierr)

        ! Structural dimension validation
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_seeds, ierr, min=0_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Handle zero seeds edge-case cleanly up-front
        if (n_seeds == 0_int32) then
            n_ensembles = 0_int32
            return
        end if

        ! Array range and element validation
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Optional input parameter validations
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, min=0.0_real64)
        call validate_in_range_int(n_tiles, ierr, min=1_int32)
        if (is_err(ierr)) return

        ! Resolve defaults
        M_DEFAULT_VAL(alpha_mad, actual_alpha_mad, 0.5_real64)
        M_DEFAULT_VAL(alpha_accept, actual_alpha_accept, 0.5_real64)
        M_DEFAULT_VAL(t_observables, actual_t_obs, 10_int32)
        M_DEFAULT_VAL(n_tiles, actual_n_tiles, CM_DEFAULT_TILE_COUNT)

        ! Allocation sizing only: tiles beyond n_seeds would receive an empty seed range and
        ! sit idle, so their workspace is never touched. Narrowing the extent here avoids that
        ! dead allocation without changing the result.
        actual_n_tiles = min(actual_n_tiles, n_seeds)

        ! Compute required observable tracking depth columns
        if (actual_t_obs <= 0_int32) then
            required_cols = n_vectors
        else
            required_cols = min(actual_t_obs, n_vectors)
        end if

        ! Determine search radius r
        if (present(r)) then
            actual_r = r
            call validate_in_range_real(actual_r, ierr, min=0.0_real64)
            if (is_err(ierr)) return
        else
            M_ALLOCATE(seed_perm(n_vectors))
            M_ALLOCATE(tmp_mean_vec(n_dimensions))
            M_ALLOCATE(tmp_distances(n_vectors))
            call calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                          tmp_mean_vec, tmp_distances, seed_perm, &
                                          actual_r, ierr=ierr)
            deallocate (tmp_mean_vec, tmp_distances, seed_perm)
            if (is_err(ierr)) return
        end if

        ! Allocate per-tile growth workspaces. Seeds inside a tile are processed serially,
        ! so a single set of workspaces per tile is reused across that tile's seeds.
        M_ALLOCATE(tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, actual_n_tiles))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors, actual_n_tiles))
        M_ALLOCATE(tmp_surface_mask(n_vectors, actual_n_tiles))
        M_ALLOCATE(tmp_perm(n_vectors, actual_n_tiles))
        M_ALLOCATE(tmp_abs_diff(n_vectors, actual_n_tiles))
        M_ALLOCATE(tmp_observables(CM_OBSERVABLE_COUNT, required_cols, actual_n_tiles))
        M_ALLOCATE(tmp_current_mask(n_vectors, actual_n_tiles))
        M_ALLOCATE(tmp_stop_reasons(n_seeds))

        call obtain_ensembles(vectors, n_dimensions, n_vectors, dimension_order, &
                              kd_indices, density_labels, seed_indices, n_seeds, &
                              actual_r, actual_alpha_mad, actual_alpha_accept, &
                              actual_t_obs, actual_n_tiles, tmp_stack, tmp_vicinity_mask, &
                              tmp_surface_mask, tmp_perm, tmp_abs_diff, tmp_observables, &
                              tmp_current_mask, &
                              ensemble_matrix, tmp_stop_reasons, tmp_mad_ambient, &
                              n_ensembles, ierr)

        if (is_err(ierr)) return

        if (present(stop_reasons)) stop_reasons = tmp_stop_reasons
        if (present(mad_ambient)) mad_ambient = tmp_mad_ambient

    end subroutine obtain_ensembles_alloc

    !> Validated Entry Point for parallel ensemble extraction.
    subroutine obtain_ensembles(vectors, n_dimensions, n_vectors, dimension_order, &
                                kd_indices, density_labels, seed_indices, n_seeds, &
                                r, alpha_mad, alpha_accept, t_observables, n_tiles, &
                                tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, &
                                tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                ensemble_matrix, stop_reasons, mad_ambient, n_ensembles, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of precalculated seed vectors [[tox_shatter_cluster_data(module):identify_ensemble_seeds_alloc(subroutine)]].
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for ambient vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: seed_indices(n_seeds)
        !! Array of starting seed vector indices [[tox_shatter_cluster_data(module):identify_ensemble_seeds_alloc(subroutine)]].
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        real(real64), intent(in) :: alpha_mad
        !! MAD multiplier for density compatibility threshold
        real(real64), intent(in) :: alpha_accept
        !! Acceptance log2 fold-change threshold
        integer(int32), intent(in) :: t_observables
        !! History tracking depth for observables
        integer(int32), intent(in) :: n_tiles
        !! Number of concurrently processed seed tiles, `n_tiles >= 1`. Tiles beyond `n_seeds` receive an empty seed range and stay idle.
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, n_tiles)
        !! Workspace stack for tree traversal per tile
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors, n_tiles)
        !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors, n_tiles)
        !! Cached spatial neighborhood union per tile, reset for each seed
        integer(int32), intent(inout) :: tmp_perm(n_vectors, n_tiles)
        !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors, n_tiles)
        !! Workspace array for deviation calculations per tile
        real(real64), intent(inout) :: tmp_observables(:, :, :)
        !! Workspace matrix for storing observable history per tile
        logical(c_bool), intent(inout) :: tmp_current_mask(n_vectors, n_tiles)
        !! Workspace holding accepted seed ensemble membership per tile
        logical(c_bool), intent(out) :: ensemble_matrix(n_vectors, n_seeds)
        !! Output matrix storing raw grown ensemble masks
        integer(int32), intent(out) :: stop_reasons(n_seeds)
        !! Per-seed growth termination reason, one of the `STC_STOP_*` constants
        real(real64), intent(out) :: mad_ambient
        !! Ambient median absolute deviation of the density labels, `0` when `alpha_mad` was inert
        integer(int32), intent(out) :: n_ensembles
        !! Total count of extracted raw ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        real(real64) :: median_ambient

        call set_ok(ierr)

        ! Structural dimension validation
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_dimension_size(n_seeds, ierr)
        if (is_err(ierr)) return

        ! Only the lower bound is a real constraint. Tiles beyond n_seeds receive an empty
        ! seed range and stay idle, so a larger n_tiles is wasteful but never incorrect.
        call validate_in_range_int(n_tiles, ierr, min=1_int32)
        if (size(tmp_observables, dim=1, kind=int32) /= CM_OBSERVABLE_COUNT) call set_err_once(ierr, ERR_DIM_MISMATCH)

        ! Ensure tmp_observables is not wider than t_observables, preventing acceptance from reading unwritten zero columns.
        if (t_observables > 0_int32) then
            if (size(tmp_observables, dim=2, kind=int32) > t_observables) then
                call set_err_once(ierr, ERR_DIM_MISMATCH)
            end if
        end if

        if (is_err(ierr)) return

        ! Array range and element validation
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=1.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Scalar parameter validation
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        ! Reject zero ambient MAD, which makes alpha_mad ineffective; recomputation is negligible relative to growth.
        call compute_ambient_density_stats_helper(density_labels, n_vectors, &
                                                  tmp_perm(:, 1), tmp_abs_diff(:, 1), &
                                                  median_ambient, mad_ambient)

        if (mad_ambient <= 0.0_real64) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call obtain_ensembles_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                     kd_indices, density_labels, seed_indices, n_seeds, &
                                     r, alpha_mad, alpha_accept, t_observables, n_tiles, &
                                     tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, &
                                     tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                     ensemble_matrix, stop_reasons, mad_ambient, n_ensembles)

    end subroutine obtain_ensembles

    !> Core Implementation for tiled parallel seed ensemble growth.
    pure subroutine obtain_ensembles_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                            kd_indices, density_labels, seed_indices, n_seeds, &
                                            r, alpha_mad, alpha_accept, t_observables, n_tiles, &
                                            tmp_stack, tmp_vicinity_mask, tmp_surface_mask, tmp_perm, &
                                            tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                            ensemble_matrix, stop_reasons, &
                                            mad_ambient, n_ensembles)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of precalculated seed vectors [[tox_shatter_cluster_data(module):identify_ensemble_seeds_alloc(subroutine)]].
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions for KD-tree
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for ambient vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: seed_indices(n_seeds)
        !! Array of starting seed vector indices [[tox_shatter_cluster_data(module):identify_ensemble_seeds_alloc(subroutine)]].
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        real(real64), intent(in) :: alpha_mad
        !! MAD multiplier for density compatibility threshold
        real(real64), intent(in) :: alpha_accept
        !! Acceptance log2 fold-change threshold
        integer(int32), intent(in) :: t_observables
        !! History tracking depth for observables
        integer(int32), intent(in) :: n_tiles
        !! Number of concurrently processed seed tiles, `n_tiles >= 1`. Tiles beyond `n_seeds` receive an empty seed range and stay idle.
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH, n_tiles)
        !! Workspace stack for tree traversal per tile
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors, n_tiles)
        !! Workspace mask holding one member's neighborhood query result per tile
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors, n_tiles)
        !! Cached spatial neighborhood union per tile, reset for each seed
        integer(int32), intent(inout) :: tmp_perm(n_vectors, n_tiles)
        !! Sorting workspace, reused for newly accepted member indices between growth rounds per tile
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors, n_tiles)
        !! Workspace array for deviation calculations per tile
        real(real64), intent(inout) :: tmp_observables(:, :, :)
        !! Workspace matrix for storing observable history per tile
        logical(c_bool), intent(inout) :: tmp_current_mask(n_vectors, n_tiles)
        !! Workspace holding accepted seed ensemble membership per tile
        logical(c_bool), intent(out) :: ensemble_matrix(n_vectors, n_seeds)
        !! Output matrix storing unmerged grown ensemble masks
        integer(int32), intent(out) :: stop_reasons(n_seeds)
        !! Per-seed growth termination reason, one of the `STC_STOP_*` constants
        real(real64), intent(out) :: mad_ambient
        !! Ambient MAD of density labels; zero makes alpha_mad ineffective by allowing only exactly equal densities.
        integer(int32), intent(out) :: n_ensembles
        !! Total count of extracted raw ensembles (equals n_seeds)

        integer(int32) :: i_tile, i_seed, seed_start, seed_end, tile_base, tile_rem
        real(real64) :: median_ambient

        ensemble_matrix = .false.

        ! Calculating the ambient density statistics.
        call compute_ambient_density_stats_helper(density_labels, n_vectors, &
                                                  tmp_perm(:, 1), tmp_abs_diff(:, 1), &
                                                  median_ambient, mad_ambient)

        ! Seeds are spread over n_tiles balanced tiles, the first tile_rem tiles taking one
        ! extra seed. Products stay bounded by n_seeds, so the index arithmetic cannot overflow.
        tile_base = n_seeds/n_tiles
        tile_rem = mod(n_seeds, n_tiles)

        ! Pure do concurrent outer-loop parallelization over independent seed tiles.
        ! Tiles execute concurrently while the seeds within a tile execute serially.
        do concurrent(i_tile=1:n_tiles) &
            shared(vectors, n_dimensions, n_vectors, dimension_order, kd_indices, &
                   density_labels, seed_indices, n_seeds, r, alpha_mad, mad_ambient, alpha_accept, &
                   t_observables, tile_base, tile_rem, tmp_stack, tmp_vicinity_mask, &
                   tmp_surface_mask, tmp_perm, tmp_abs_diff, tmp_observables, &
                   tmp_current_mask, ensemble_matrix, stop_reasons) &
            local(i_seed, seed_start, seed_end)

            if (i_tile <= tile_rem) then
                seed_start = (i_tile - 1_int32)*(tile_base + 1_int32) + 1_int32
                seed_end = seed_start + tile_base
            else
                seed_start = tile_rem*(tile_base + 1_int32) + &
                             (i_tile - tile_rem - 1_int32)*tile_base + 1_int32
                seed_end = seed_start + tile_base - 1_int32
            end if

            do i_seed = seed_start, seed_end

                call grow_single_seed_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                             kd_indices, density_labels, seed_indices(i_seed), &
                                             r, alpha_mad, mad_ambient, alpha_accept, t_observables, &
                                             tmp_stack(:, :, i_tile), tmp_vicinity_mask(:, i_tile), &
                                             tmp_surface_mask(:, i_tile), tmp_perm(:, i_tile), &
                                             tmp_abs_diff(:, i_tile), tmp_observables(:, :, i_tile), &
                                             tmp_current_mask(:, i_tile), &
                                             ensemble_matrix(:, i_seed), stop_reasons(i_seed))

            end do

        end do

        n_ensembles = n_seeds

    end subroutine obtain_ensembles_helper

    !> Core Implementation for growing a single seed into an ensemble.
    pure subroutine grow_single_seed_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                            kd_indices, density_labels, seed_idx, &
                                            r, alpha_mad, mad_ambient, alpha_accept, t_observables, &
                                            tmp_stack, tmp_vicinity_mask, tmp_surface_mask, &
                                            tmp_perm, tmp_abs_diff, tmp_observables, &
                                            current_mask, out_mask, stop_reason)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence array tracking split dimensions for KD-tree
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Precalculated density labels for ambient vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: seed_idx
        !! Index of starting seed vector
        real(real64), intent(in) :: r
        !! Search radius threshold for surface growth
        real(real64), intent(in) :: alpha_mad
        !! MAD multiplier for density compatibility threshold
        real(real64), intent(in) :: mad_ambient
        !! Ambient MAD precalculated via [[tox_shatter_cluster_data(module):compute_ambient_density_stats_helper(subroutine)]].
        real(real64), intent(in) :: alpha_accept
        !! Acceptance log2 fold-change threshold
        integer(int32), intent(in) :: t_observables
        !! History tracking depth for observables
        integer(int32), intent(inout) :: tmp_stack(KD_STACK_ENTRY_SIZE, KD_TRAVERSAL_STACK_DEPTH)
        !! Workspace stack for a single tree traversal
        logical(c_bool), intent(inout) :: tmp_vicinity_mask(n_vectors)
        !! Workspace mask holding one member's neighborhood query result
        logical(c_bool), intent(inout) :: tmp_surface_mask(n_vectors)
        !! Cached geometric union for this seed; reset once and never density-filtered
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Newly accepted indices between rounds; reused for median sorting after their queries
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Workspace array holding the active member densities
        real(real64), intent(inout) :: tmp_observables(:, :)
        !! Workspace matrix for storing observable history
        logical(c_bool), intent(inout) :: current_mask(n_vectors)
        !! Accepted ensemble membership; modified only after a trial passes acceptance
        logical(c_bool), intent(out) :: out_mask(n_vectors)
        !! Final accepted membership; reused as trial membership until growth terminates
        integer(int32), intent(out) :: stop_reason
        !! Reason growth terminated, one of the `STC_STOP_*` constants

        integer(int32) :: iter, prev_count, curr_count, candidate_count, i_vec, i_frontier, frontier_idx, n_frontier
        logical(c_bool) :: is_growing, is_accepted

        current_mask = .false.
        current_mask(seed_idx) = .true.
        tmp_observables = 0.0_real64
        tmp_surface_mask = .false.
        n_frontier = 1_int32
        tmp_perm(1) = seed_idx

        iter = 1_int32
        call compute_ensemble_observable_helper(current_mask, density_labels, &
                                                n_vectors, 0_int32, iter, &
                                                tmp_observables, &
                                                t_observables)

        is_growing = .true.
        do while (is_growing)
            prev_count = count(current_mask)

            ! Each accepted vector is queried once. Consume the frontier before median sorting
            ! reuses tmp_perm; candidates from this round are queried only after acceptance.
            do i_frontier = 1, n_frontier
                frontier_idx = tmp_perm(i_frontier)
                call vicinity_vectors_helper(vectors(:, frontier_idx), vectors, n_dimensions, n_vectors, r, &
                                             dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask)
                tmp_surface_mask = tmp_surface_mask .or. tmp_vicinity_mask
            end do

            call build_ensemble_trial_helper(current_mask, n_vectors, density_labels, alpha_mad, &
                                             mad_ambient, tmp_surface_mask, tmp_perm, tmp_abs_diff, out_mask)

            curr_count = count(out_mask)
            candidate_count = curr_count - prev_count

            if (candidate_count == 0_int32) then
                ! iter is still 1 only on the very first pass, so the seed never grew at all
                if (iter == 1_int32) then
                    stop_reason = CM_STOP_NO_CANDIDATES
                else
                    stop_reason = CM_STOP_FIXED_POINT
                end if
                is_growing = .false.
                exit
            end if

            iter = iter + 1_int32

            call compute_ensemble_observable_helper(out_mask, density_labels, &
                                                    n_vectors, candidate_count, iter, &
                                                    tmp_observables, &
                                                    t_observables)

            call accept_ensemble_helper(tmp_observables, iter, &
                                        alpha_accept, is_accepted)

            if (.not. is_accepted) then
                ! iter counts the growth rounds that found candidates, so iter == 2 means this
                ! was the first such round and no batch was ever accepted
                if (iter == 2_int32) then
                    stop_reason = CM_STOP_NEVER_ACCEPTED
                else
                    stop_reason = CM_STOP_REJECT_AFTER_ACCEPT
                end if

                is_growing = .false.
            else
                ! Commit and collect only new members for next round's spatial queries.
                n_frontier = 0_int32
                do i_vec = 1, n_vectors
                    if (out_mask(i_vec) .and. .not. current_mask(i_vec)) then
                        current_mask(i_vec) = .true.
                        n_frontier = n_frontier + 1_int32
                        tmp_perm(n_frontier) = i_vec
                    end if
                end do
            end if
        end do

        out_mask = current_mask

    end subroutine grow_single_seed_helper

    !> Core Implementation for a union-find root lookup with path compression.
    pure subroutine find_root_helper(parent, n_nodes, node, root)

        integer(int32), intent(in) :: n_nodes
        !! Number of nodes in the forest
        integer(int32), intent(inout) :: parent(n_nodes)
        !! Parent pointers; flattened towards the root as a side effect of the lookup
        integer(int32), intent(in) :: node
        !! Node whose component root is wanted
        integer(int32), intent(out) :: root
        !! Representative of the component containing `node`

        integer(int32) :: walker, next_node

        root = node
        do while (parent(root) /= root)
            root = parent(root)
        end do

        ! Path compression keeps later lookups near constant time and prevents long union chains from restoring quadratic behaviour.
        walker = node
        do while (parent(walker) /= root)
            next_node = parent(walker)
            parent(walker) = root
            walker = next_node
        end do

    end subroutine find_root_helper

    !> Allocating Wrapper for transitive set-union ensemble merging.
    subroutine merge_ensembles_alloc(raw_masks, n_vectors, n_seeds, &
                                     min_intersection, merged_matrix, &
                                     n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        logical(c_bool), intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Input matrix of raw boolean ensemble masks [n_vectors x n_seeds]
        integer(int32), intent(in), optional :: min_intersection
        !! Minimum points in common to trigger a merge (defaults to 1)
        logical(c_bool), allocatable, intent(out) :: merged_matrix(:, :)
        !! Output matrix storing final merged non-singleton ensemble masks [n_vectors x n_ensembles]
        integer(int32), intent(out) :: n_ensembles
        !! Final count of merged non-singleton ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: actual_min_intersect
        logical(c_bool), allocatable :: tmp_merged_masks(:, :)
        logical(c_bool), allocatable :: tmp_active_flag(:)
        integer(int32), allocatable :: tmp_parent(:)

        call set_ok(ierr)

        ! Structural dimension validation
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_seeds, ierr, min=0_int32, max=n_vectors)
        if (is_err(ierr)) return

        if (n_seeds == 0_int32) then
            M_ALLOCATE(merged_matrix(n_vectors, 0_int32))
            n_ensembles = 0_int32
            return
        end if

        call validate_in_range_int(min_intersection, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        M_DEFAULT_VAL(min_intersection, actual_min_intersect, 1_int32)

        M_ALLOCATE(tmp_merged_masks(n_vectors, n_seeds))
        M_ALLOCATE(tmp_active_flag(n_seeds))
        M_ALLOCATE(tmp_parent(n_seeds))

        call merge_ensembles(raw_masks, n_vectors, n_seeds, actual_min_intersect, &
                             tmp_merged_masks, tmp_active_flag, tmp_parent, n_ensembles, ierr)

        if (is_err(ierr)) return

        ! Slice final matrix output to actual merged count
        M_ALLOCATE(merged_matrix(n_vectors, n_ensembles))
        if (n_ensembles > 0_int32) then
            merged_matrix(:, 1:n_ensembles) = tmp_merged_masks(:, 1:n_ensembles)
        end if

    end subroutine merge_ensembles_alloc

    !> Validated Entry Point for transitive set-union ensemble merging.
    subroutine merge_ensembles(raw_masks, n_vectors, n_seeds, min_intersection, &
                               merged_masks, tmp_active_flag, tmp_parent, n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        logical(c_bool), intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Input matrix of raw boolean ensemble masks [n_vectors x n_seeds]
        integer(int32), intent(in) :: min_intersection
        !! Minimum points in common to trigger a merge
        logical(c_bool), intent(out) :: merged_masks(n_vectors, n_seeds)
        !! Workspace matrix storing final merged ensemble masks
        logical(c_bool), intent(inout) :: tmp_active_flag(n_seeds)
        !! Preallocated workspace tracking active unmerged seeds
        integer(int32), intent(inout) :: tmp_parent(n_seeds)
        !! Preallocated union-find parent pointers, used when `min_intersection` is 1
        integer(int32), intent(out) :: n_ensembles
        !! Final count of merged non-singleton ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        call set_ok(ierr)

        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_seeds, ierr, min=0_int32, max=n_vectors)
        call validate_in_range_int(min_intersection, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        call merge_ensembles_helper(raw_masks, n_vectors, n_seeds, min_intersection, &
                                    merged_masks, tmp_active_flag, tmp_parent, n_ensembles, ierr)

    end subroutine merge_ensembles

    !> Core Implementation for pairwise ensemble merging based on set intersection.
    pure subroutine merge_ensembles_helper(raw_masks, n_vectors, n_seeds, &
                                           min_intersection, merged_masks, &
                                           tmp_active_flag, tmp_parent, n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        integer(int32), intent(in) :: min_intersection
        !! Minimum overlapping vectors required to merge two ensembles
        logical(c_bool), intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Matrix of raw unmerged boolean ensemble masks
        logical(c_bool), intent(out) :: merged_masks(n_vectors, n_seeds)
        !! Output matrix storing merged boolean ensemble masks
        logical(c_bool), intent(inout) :: tmp_active_flag(n_seeds)
        !! Preallocated workspace tracking active unmerged seeds
        integer(int32), intent(inout) :: tmp_parent(n_seeds)
        !! Preallocated union-find parent pointers, used when `min_intersection` is 1
        integer(int32), intent(out) :: n_ensembles
        !! Final count of merged non-singleton ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: i, j, shared_count, root_i, root_j
        logical(c_bool) :: merged_any

        call set_ok(ierr)

        if (n_seeds == 0_int32) then
            n_ensembles = 0_int32
            merged_masks = .false.
            return
        end if

        if (min_intersection == 1_int32) then
            ! Union-find finds the overlap graph's connected components in one O(S^2*N) pass instead of repeated pairwise sweeps.
            do i = 1, n_seeds
                tmp_parent(i) = i
            end do

            do i = 1, n_seeds
                do j = i + 1, n_seeds
                    ! any() may stop at the first shared vector, count() could not
                    if (any(raw_masks(:, i) .and. raw_masks(:, j))) then
                        call find_root_helper(tmp_parent, n_seeds, i, root_i)
                        call find_root_helper(tmp_parent, n_seeds, j, root_j)

                        ! Attach to the smaller index so each component keeps the lowest column, matching the iterative branch.
                        if (root_i < root_j) then
                            tmp_parent(root_j) = root_i
                        else if (root_j < root_i) then
                            tmp_parent(root_i) = root_j
                        end if
                    end if
                end do
            end do

            merged_masks = .false.
            tmp_active_flag = .false.

            do i = 1, n_seeds
                call find_root_helper(tmp_parent, n_seeds, i, root_i)
                merged_masks(:, root_i) = merged_masks(:, root_i) .or. raw_masks(:, i)
                tmp_active_flag(root_i) = .true.
            end do
        else

            ! For larger overlap thresholds, merged overlaps accumulate, requiring iteration to a fixed point.
            merged_masks = raw_masks
            tmp_active_flag = .true.

            do
                merged_any = .false.
                do i = 1, n_seeds
                    if (.not. tmp_active_flag(i)) cycle

                    do j = i + 1, n_seeds
                        if (.not. tmp_active_flag(j)) cycle

                        ! Evaluate common members between ensemble i and ensemble j
                        shared_count = count(merged_masks(:, i) .and. merged_masks(:, j))

                        if (shared_count >= min_intersection) then
                            ! Union sets into ensemble i and deactivate ensemble j
                            merged_masks(:, i) = merged_masks(:, i) .or. merged_masks(:, j)
                            tmp_active_flag(j) = .false.
                            merged_any = .true.
                        end if
                    end do
                end do

                if (.not. merged_any) exit
            end do
        end if

        ! Compact merged ensembles to the front, omitting singletons as unassigned/noise.
        n_ensembles = 0_int32
        do i = 1, n_seeds
            if (tmp_active_flag(i) .and. count(merged_masks(:, i)) > 1_int32) then
                n_ensembles = n_ensembles + 1_int32
                if (n_ensembles /= i) then
                    merged_masks(:, n_ensembles) = merged_masks(:, i)
                end if
            end if
        end do

        ! Clear post-compaction columns to prevent stale or duplicate ensembles from being read beyond n_ensembles.
        merged_masks(:, n_ensembles + 1_int32:n_seeds) = .false.

    end subroutine merge_ensembles_helper

end module tox_shatter_cluster_data
