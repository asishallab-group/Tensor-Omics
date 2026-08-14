#include "macros.h"
#define CM_KD_STACK_ENTRY_SIZE 3
#define CM_KD_TRAVERSAL_STACK_DEPTH 64
#define CM_OBSERVABLE_COUNT 5
#define CM_SEEDING_COVERAGE_PERCENTILE 50.0_real64

!> Module for managing data structures and geometric calculations.
!! For shatter clustering operations within Tensor Omics.
module tox_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, &
                          validate_all_in_range_int, validate_in_range_int, ERR_OK, ERR_ALLOC_FAIL, &
                          ERR_DIM_MISMATCH, ERR_INVALID_INPUT, validate_all_in_range_real
    use tox_gene_centroids, only: mean_vector
    use tox_euclidean_distance, only: euclidean_distance
    use f42_utils, only: sort_real_heapsort, calc_percentile, calc_percentile_helper, init_perm
    use f42_kd_tree, only: vicinity_vectors, vicinity_vectors_count
    implicit none
    private
    public :: calculate_density_radius_alloc, calculate_density_radius, calculate_density_radius_helper
    public :: calculate_labels_as_density_alloc, calculate_labels_as_density, calculate_labels_as_density_helper
    public :: identify_ensemble_seeds_alloc, identify_ensemble_seeds, identify_ensemble_seeds_helper
    public :: grow_ensemble_alloc, grow_ensemble, grow_ensemble_helper
    public :: compute_ensemble_observable_alloc, compute_ensemble_observable, compute_ensemble_observable_helper
    public :: accept_ensemble, accept_ensemble_helper
    public :: obtain_ensembles_alloc, obtain_ensembles, obtain_ensembles_helper
    public :: grow_single_seed_helper
    public :: merge_ensembles_alloc, merge_ensembles, merge_ensembles_helper

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

        real(real64)   :: actual_quant, percentile
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
            call euclidean_distance(tmp_mean_vec, vectors(:, i_vec), n_dimensions, tmp_distances(i_vec))
        end do

        !Sorting perm according to distance
        call sort_real_heapsort(tmp_distances, tmp_perm)

        percentile = actual_quant*100.0_real64

        ! Extract the adaptive search radius corresponding to the specified distance percentile.
        call calc_percentile(tmp_distances, tmp_perm, percentile, radius, ierr)

    end subroutine calculate_density_radius_helper

    !> Allocating Wrapper for calculating label density distributions.
    subroutine calculate_labels_as_density_alloc(vectors, n_dimensions, n_vectors, r, &
                                                 dimension_order, kd_indices, label_densities, ierr)

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
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_stack(:, :, :)

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)

        ! Array and value range validation checks
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors))

        call calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                         dimension_order, kd_indices, tmp_stack, &
                                         label_densities, ierr)

    end subroutine calculate_labels_as_density_alloc

    !> Validated Entry Point for calculating label density distributions.
    subroutine calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                           dimension_order, kd_indices, tmp_stack, &
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
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Preallocated workspace stack for tree traversal
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
        if (is_err(ierr)) return

        call calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                dimension_order, kd_indices, tmp_stack, &
                                                label_densities, ierr)

    end subroutine calculate_labels_as_density

    !> Core Implementation for calculating label density coordinates.
    pure subroutine calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                       dimension_order, kd_indices, tmp_stack, &
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
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Preallocated workspace stack for tree traversal
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output array storing generated density scalars
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i_vec, neighbor_count

        call set_ok(ierr)

        do concurrent(i_vec=1:n_vectors) &
            shared(vectors, n_dimensions, n_vectors, r, dimension_order, kd_indices, tmp_stack, label_densities) &
            local(neighbor_count)

            ! Query point count directly into thread-local scalar
            call vicinity_vectors_count(vectors(:, i_vec), vectors, n_dimensions, n_vectors, r, &
                                        dimension_order, kd_indices, tmp_stack(:, :, i_vec), &
                                        neighbor_count)

            label_densities(i_vec) = real(neighbor_count, real64)
        end do

    end subroutine calculate_labels_as_density_helper

    !> Allocating wrapper for greedy, density-ranked, coverage-based seed selection.
    subroutine identify_ensemble_seeds_alloc(vectors, n_dimensions, n_vectors, density_labels, &
                                             dimension_order, kd_indices, k_seeding, &
                                             sorted_perm, n_seeds, seed_mask, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical, intent(out) :: seed_mask(n_vectors)
        !! Output logical mask; .true. marks vectors selected as seeds
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_perm(:)
        real(real64), allocatable :: tmp_distances(:)
        integer(int32), allocatable :: tmp_stack(:, :)
        logical, allocatable :: tmp_visited_mask(:)
        logical, allocatable :: tmp_newly_covered_mask(:)

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=0.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        if (n_vectors < 2_int32) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        call validate_in_range_int(k_seeding, ierr, min=1_int32, max=n_vectors - 1_int32)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH))
        M_ALLOCATE(tmp_visited_mask(n_vectors))
        M_ALLOCATE(tmp_newly_covered_mask(n_vectors))

        call identify_ensemble_seeds(vectors, n_dimensions, n_vectors, density_labels, &
                                     dimension_order, kd_indices, k_seeding, &
                                     tmp_perm, tmp_distances, tmp_stack, &
                                     tmp_visited_mask, tmp_newly_covered_mask, &
                                     sorted_perm, n_seeds, seed_mask, ierr)

    end subroutine identify_ensemble_seeds_alloc

    !> Validated entry point for greedy, density-ranked, coverage-based seed selection.
    subroutine identify_ensemble_seeds(vectors, n_dimensions, n_vectors, density_labels, &
                                       dimension_order, kd_indices, k_seeding, &
                                       tmp_perm, tmp_distances, tmp_stack, &
                                       tmp_visited_mask, tmp_newly_covered_mask, &
                                       sorted_perm, n_seeds, seed_mask, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated workspace for indirect sorting
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated workspace for candidate-to-vector distances
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated K-D tree traversal stack
        logical, intent(out) :: tmp_visited_mask(n_vectors)
        !! Preallocated workspace tracking points already covered by seeding
        logical, intent(out) :: tmp_newly_covered_mask(n_vectors)
        !! Preallocated workspace for the current seed's coverage query
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical, intent(out) :: seed_mask(n_vectors)
        !! Output logical mask; .true. marks vectors selected as seeds
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=0.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        if (n_vectors < 2_int32) call set_err(ierr, ERR_INVALID_INPUT)
        if (is_err(ierr)) return

        call validate_in_range_int(k_seeding, ierr, min=1_int32, max=n_vectors - 1_int32)
        if (is_err(ierr)) return

        call identify_ensemble_seeds_helper(vectors, n_dimensions, n_vectors, density_labels, &
                                            dimension_order, kd_indices, k_seeding, &
                                            tmp_perm, tmp_distances, tmp_stack, &
                                            tmp_visited_mask, tmp_newly_covered_mask, &
                                            sorted_perm, n_seeds, seed_mask, ierr)

    end subroutine identify_ensemble_seeds

    !> Core implementation for greedy, density-ranked, coverage-based seed selection.
    !| Starts from the highest-density unvisited vector, marks it as a seed, estimates a
    !| local coverage radius from the median distance to its `k_seeding` nearest neighbors,
    !| marks all vectors inside that radius as visited, and continues in density order.
    pure subroutine identify_ensemble_seeds_helper(vectors, n_dimensions, n_vectors, density_labels, &
                                                   dimension_order, kd_indices, k_seeding, &
                                                   tmp_perm, tmp_distances, tmp_stack, &
                                                   tmp_visited_mask, tmp_newly_covered_mask, &
                                                   sorted_perm, n_seeds, seed_mask, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: density_labels(n_vectors)
        !! Density labels for all vectors [[tox_shatter_cluster_data(module):calculate_labels_as_density_alloc(subroutine)]].
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the K-D tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! K-D tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(in) :: k_seeding
        !! Number of nearest neighbors used to estimate each seed's local coverage radius
        integer(int32), intent(out) :: tmp_perm(n_vectors)
        !! Preallocated workspace for indirect sorting
        real(real64), intent(out) :: tmp_distances(n_vectors)
        !! Preallocated workspace for candidate-to-vector distances
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated K-D tree traversal stack
        logical, intent(out) :: tmp_visited_mask(n_vectors)
        !! Preallocated workspace tracking points already covered by seeding
        logical, intent(out) :: tmp_newly_covered_mask(n_vectors)
        !! Preallocated workspace for the current seed's coverage query
        integer(int32), intent(out) :: sorted_perm(n_vectors)
        !! Vector indices sorted by density descending
        integer(int32), intent(out) :: n_seeds
        !! Number of selected seed vectors
        logical, intent(out) :: seed_mask(n_vectors)
        !! Output logical mask; .true. marks vectors selected as seeds
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i, rank, candidate, self_pos
        real(real64) :: coverage_radius

        call set_ok(ierr)

        ! Rank vectors by density, descending.
        call init_perm(tmp_perm)
        call sort_real_heapsort(density_labels, tmp_perm)

        do concurrent(i=1:n_vectors) shared(sorted_perm, tmp_perm, n_vectors)
            sorted_perm(i) = tmp_perm(n_vectors - i + 1_int32)
        end do

        seed_mask = .false.
        tmp_visited_mask = .false.
        n_seeds = 0_int32

        do rank = 1, n_vectors
            candidate = sorted_perm(rank)
            if (tmp_visited_mask(candidate)) cycle

            seed_mask(candidate) = .true.
            tmp_visited_mask(candidate) = .true.
            n_seeds = n_seeds + 1_int32

            ! Compute all candidate distances.
            do concurrent(i=1:n_vectors) shared(vectors, n_dimensions, candidate, tmp_distances)
                call euclidean_distance(vectors(:, candidate), vectors(:, i), &
                                        n_dimensions, tmp_distances(i))
            end do

            call init_perm(tmp_perm)
            call sort_real_heapsort(tmp_distances, tmp_perm)

            ! Exclude the seed itself while preserving ascending distance order.
            self_pos = 0_int32
            do i = 1, n_vectors
                if (tmp_perm(i) == candidate) then
                    self_pos = i
                    exit
                end if
            end do

            if (self_pos >= 1_int32 .and. self_pos <= k_seeding) then
                do i = self_pos, k_seeding
                    tmp_perm(i) = tmp_perm(i + 1_int32)
                end do
            end if

            call calc_percentile_helper(tmp_distances, tmp_perm(1:k_seeding), &
                                        CM_SEEDING_COVERAGE_PERCENTILE, coverage_radius)

            call vicinity_vectors(vectors(:, candidate), vectors, n_dimensions, n_vectors, &
                                  coverage_radius, dimension_order, kd_indices, tmp_stack, &
                                  tmp_newly_covered_mask)

            do concurrent(i=1:n_vectors) shared(tmp_visited_mask, tmp_newly_covered_mask)
                if (tmp_newly_covered_mask(i)) tmp_visited_mask(i) = .true.
            end do
        end do

    end subroutine identify_ensemble_seeds_helper

    !> Allocating Wrapper for growing the active ensemble at its surface.
    subroutine grow_ensemble_alloc(vectors, n_dimensions, n_vectors, &
                                   ensemble_mask, dimension_order, kd_indices, &
                                   density_labels, r, alpha_mad, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical, intent(inout) :: ensemble_mask(n_vectors)
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
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_stack(:, :, :)
        logical, allocatable :: tmp_vicinity_mask(:, :)
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
        call validate_all_in_range_real(density_labels, n_vectors, ierr)

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
        M_ALLOCATE(tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors, n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))
        M_ALLOCATE(tmp_abs_diff(n_vectors))

        call grow_ensemble(vectors, n_dimensions, n_vectors, ensemble_mask, &
                           actual_r, dimension_order, kd_indices, density_labels, &
                           actual_alpha_mad, tmp_stack, tmp_vicinity_mask, &
                           tmp_perm, tmp_abs_diff, ierr)

    end subroutine grow_ensemble_alloc

    !> Validated Entry Point for surface-based ensemble growth.
    subroutine grow_ensemble(vectors, n_dimensions, n_vectors, ensemble_mask, &
                             r, dimension_order, kd_indices, density_labels, &
                             alpha_mad, tmp_stack, tmp_vicinity_mask, &
                             tmp_perm, tmp_abs_diff, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical, intent(inout) :: ensemble_mask(n_vectors)
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
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Preallocated workspace stack for parallel tree traversal
        logical, intent(inout) :: tmp_vicinity_mask(n_vectors, n_vectors)
        !! Preallocated logical matrix tracking neighbor mappings
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array for MAD calculation
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr)

        ! Array and value range validation checks
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call grow_ensemble_helper(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                  r, dimension_order, kd_indices, density_labels, &
                                  alpha_mad, tmp_stack, tmp_vicinity_mask, &
                                  tmp_perm, tmp_abs_diff, ierr)

    end subroutine grow_ensemble

    !> Core Implementation for surface-based ensemble growth.
    pure subroutine grow_ensemble_helper(vectors, n_dimensions, n_vectors, ensemble_mask, &
                                         r, dimension_order, kd_indices, density_labels, &
                                         alpha_mad, tmp_stack, tmp_vicinity_mask, &
                                         tmp_perm, tmp_abs_diff, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        logical, intent(inout) :: ensemble_mask(n_vectors)
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
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Preallocated workspace stack for parallel tree traversal
        logical, intent(inout) :: tmp_vicinity_mask(n_vectors, n_vectors)
        !! Preallocated logical matrix tracking neighbor mappings
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Preallocated workspace array for sorting and percentiles
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Preallocated workspace array for MAD calculation
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i_vec, j, i_member, k, n_active
        logical :: is_neighbor
        real(real64) :: median_ambient, mad_ambient, ensemble_center_density, max_allowed_dev

        call set_ok(ierr)

        ! 1. Calculate ambient median density: median_{i=1..N}(rho_i)
        do concurrent(i_vec=1:n_vectors) shared(tmp_perm)
            tmp_perm(i_vec) = i_vec
        end do
        call sort_real_heapsort(density_labels, tmp_perm)
        call calc_percentile(density_labels, tmp_perm, 50.0_real64, median_ambient, ierr)
        if (is_err(ierr)) return

        ! 2. Compute absolute deviation from ambient median for all vectors
        do concurrent(i_vec=1:n_vectors) shared(tmp_abs_diff, density_labels, median_ambient)
            tmp_abs_diff(i_vec) = abs(density_labels(i_vec) - median_ambient)
            tmp_perm(i_vec) = i_vec
        end do

        ! 3. Compute MAD_ambient = median_{i=1..N}(|rho_i - median_ambient|)
        call sort_real_heapsort(tmp_abs_diff, tmp_perm)
        call calc_percentile(tmp_abs_diff, tmp_perm, 50.0_real64, mad_ambient, ierr)
        if (is_err(ierr)) return

        ! Max absolute deviation allowed from ensemble central density
        max_allowed_dev = alpha_mad*mad_ambient

        ! 4. Calculate current ensemble central density (median of active members)
        n_active = count(ensemble_mask)
        if (n_active == 0_int32) return

        k = 0
        do i_vec = 1, n_vectors
            if (ensemble_mask(i_vec)) then
                k = k + 1
                tmp_abs_diff(k) = density_labels(i_vec)
                tmp_perm(k) = k
            end if
        end do
        call sort_real_heapsort(tmp_abs_diff(1:k), tmp_perm(1:k))
        call calc_percentile(tmp_abs_diff(1:k), tmp_perm(1:k), 50.0_real64, ensemble_center_density, ierr)
        if (is_err(ierr)) return

        ! 5. Parallel Neighborhood Discovery across active members using vicinity_vectors
        do concurrent(i_vec=1:n_vectors) &
            shared(ensemble_mask, vectors, n_dimensions, n_vectors, r, dimension_order, kd_indices, tmp_stack, tmp_vicinity_mask)

            if (ensemble_mask(i_vec)) then
                call vicinity_vectors(vectors(:, i_vec), vectors, n_dimensions, n_vectors, r, &
                                      dimension_order, kd_indices, tmp_stack(:, :, i_vec), &
                                      tmp_vicinity_mask(:, i_vec))
            else
                tmp_vicinity_mask(:, i_vec) = .false.
            end if
        end do

        ! 6. Surface-growth update: Add non-members touching active surface if density-compatible
        do concurrent(j=1:n_vectors) &
            shared(ensemble_mask, tmp_vicinity_mask, density_labels, ensemble_center_density, max_allowed_dev, n_vectors) &
            local(is_neighbor, i_member)

            if (.not. ensemble_mask(j)) then
                is_neighbor = .false.
                do i_member = 1, n_vectors
                    if (tmp_vicinity_mask(j, i_member)) then
                        is_neighbor = .true.
                        exit
                    end if
                end do

                ! Check proximity and density compatibility |rho_j - rho_ensemble| <= alpha_mad * MAD_ambient
                if (is_neighbor .and. (abs(density_labels(j) - ensemble_center_density) <= max_allowed_dev)) then
                    ensemble_mask(j) = .true.
                end if
            end if
        end do

    end subroutine grow_ensemble_helper

    !> Allocating Wrapper for calculating and storing ensemble observable trajectories.
    subroutine compute_ensemble_observable_alloc(ensemble_mask, density_labels, n_vectors, &
                                                 candidate_count, current_iter, observables, &
                                                 t_observables, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical, intent(in) :: ensemble_mask(n_vectors)
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
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=0.0_real64)
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

        call compute_ensemble_observable(ensemble_mask, density_labels, n_vectors, &
                                         candidate_count, current_iter, observables, &
                                         actual_t_obs, ierr)

    end subroutine compute_ensemble_observable_alloc

    !> Validated Entry Point for computing and updating 2D ensemble density observables.
    subroutine compute_ensemble_observable(ensemble_mask, density_labels, n_vectors, &
                                           candidate_count, current_iter, observables, &
                                           t_observables, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical, intent(in) :: ensemble_mask(n_vectors)
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
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                candidate_count, current_iter, observables, &
                                                t_observables, ierr)

    end subroutine compute_ensemble_observable

    !> Core Implementation for calculating 5-component observables into a 2D history matrix.
    pure subroutine compute_ensemble_observable_helper(ensemble_mask, density_labels, n_vectors, &
                                                       candidate_count, current_iter, observables, &
                                                       t_observables, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        logical, intent(in) :: ensemble_mask(n_vectors)
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

        integer(int32) :: i_vec, active_count, actual_t_obs, max_cols, zero_count
        real(real64)   :: sum_arithmetic, sum_reciprocal, rho_arith, rho_harm
        real(real64)   :: obs_1, obs_2, obs_3, obs_4, obs_5

        call set_ok(ierr)

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

            ! 2. Harmonic Mean Density (strictly 0 if any active element is <= 0)
            if (zero_count == 0_int32 .and. sum_reciprocal > 0.0_real64) then
                rho_harm = real(active_count, real64)/sum_reciprocal
            else
                rho_harm = 0.0_real64
            end if

            obs_1 = rho_arith
            obs_2 = rho_harm

            ! 3. Density Heterogeneity H_rho
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
        logical, intent(out) :: is_accepted
        !! Output flag indicating step acceptance decision
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate input structural dimension
        call validate_dimension_size(current_iter, ierr)

        ! Optional parameter range validation
        call validate_in_range_real(alpha_accept, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call accept_ensemble_helper(observables, current_iter, alpha_accept, &
                                    is_accepted, ierr)

    end subroutine accept_ensemble

    !> Core Implementation for decision rule evaluating step acceptance.
    pure subroutine accept_ensemble_helper(observables, current_iter, alpha_accept, &
                                           is_accepted, ierr)

        integer(int32), intent(in) :: current_iter
        !! Current growth iteration index (1-based)
        real(real64), intent(in) :: observables(:, :)
        !! 2D observables trajectory matrix [5 x window_capacity]
        real(real64), intent(in), optional :: alpha_accept
        !! Optional log2 fold change acceptance threshold (defaults to 0.5)
        logical, intent(out) :: is_accepted
        !! Output flag indicating step acceptance decision
        integer(int32), intent(out) :: ierr
        !! Error code

        real(real64)   :: actual_alpha, h_prev, h_curr, log2_fold_change
        integer(int32) :: max_cols, col_curr, col_prev

        call set_ok(ierr)

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
                                      alpha_accept, t_observables, &
                                      ensemble_matrix, n_ensembles, ierr)

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
        logical, allocatable, intent(out) :: ensemble_matrix(:, :)
        !! Output 2D logical matrix [n_vectors x n_seeds] tracking grown raw ensembles
        integer(int32), intent(out) :: n_ensembles
        !! Total count of extracted raw ensembles (equals n_seeds)
        integer(int32), intent(out) :: ierr
        !! Error status flag

        real(real64) :: actual_r, actual_alpha_mad, actual_alpha_accept
        integer(int32) :: actual_t_obs, required_cols

        ! Workspaces
        integer(int32), allocatable :: seed_perm(:)
        integer(int32), allocatable :: tmp_stack(:, :, :, :)
        logical, allocatable :: tmp_vicinity_mask(:, :, :)
        integer(int32), allocatable :: tmp_perm(:, :)
        real(real64), allocatable :: tmp_abs_diff(:, :)
        real(real64), allocatable :: tmp_observables(:, :, :)
        logical, allocatable :: tmp_current_mask(:, :)
        logical, allocatable :: tmp_backup_mask(:, :)
        real(real64), allocatable :: tmp_mean_vec(:), tmp_distances(:)

        call set_ok(ierr)

        ! Structural dimension validation
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_seeds, ierr, min=0_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Handle zero seeds edge-case cleanly up-front
        if (n_seeds == 0_int32) then
            M_ALLOCATE(ensemble_matrix(n_vectors, 0_int32))
            n_ensembles = 0_int32
            return
        end if

        ! Array range and element validation
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=0.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Optional input parameter validations
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        ! Resolve defaults
        M_DEFAULT_VAL(alpha_mad, actual_alpha_mad, 0.5_real64)
        M_DEFAULT_VAL(alpha_accept, actual_alpha_accept, 0.5_real64)
        M_DEFAULT_VAL(t_observables, actual_t_obs, 10_int32)

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

        ! Allocate per-seed parallel growth workspaces
        M_ALLOCATE(tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors, n_seeds))
        M_ALLOCATE(tmp_vicinity_mask(n_vectors, n_vectors, n_seeds))
        M_ALLOCATE(tmp_perm(n_vectors, n_seeds))
        M_ALLOCATE(tmp_abs_diff(n_vectors, n_seeds))
        M_ALLOCATE(tmp_observables(CM_OBSERVABLE_COUNT, required_cols, n_seeds))
        M_ALLOCATE(tmp_current_mask(n_vectors, n_seeds))
        M_ALLOCATE(tmp_backup_mask(n_vectors, n_seeds))
        M_ALLOCATE(ensemble_matrix(n_vectors, n_seeds))

        call obtain_ensembles(vectors, n_dimensions, n_vectors, dimension_order, &
                              kd_indices, density_labels, seed_indices, n_seeds, &
                              actual_r, actual_alpha_mad, actual_alpha_accept, &
                              actual_t_obs, tmp_stack, tmp_vicinity_mask, &
                              tmp_perm, tmp_abs_diff, tmp_observables, &
                              tmp_current_mask, tmp_backup_mask, &
                              ensemble_matrix, n_ensembles, ierr)

    end subroutine obtain_ensembles_alloc

    !> Validated Entry Point for parallel ensemble extraction.
    subroutine obtain_ensembles(vectors, n_dimensions, n_vectors, dimension_order, &
                                kd_indices, density_labels, seed_indices, n_seeds, &
                                r, alpha_mad, alpha_accept, t_observables, &
                                tmp_stack, tmp_vicinity_mask, tmp_perm, &
                                tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                tmp_backup_mask, ensemble_matrix, n_ensembles, ierr)

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
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors, n_seeds)
        !! Workspace stack for tree traversal per seed thread
        logical, intent(inout) :: tmp_vicinity_mask(n_vectors, n_vectors, n_seeds)
        !! Workspace matrix for neighborhood mapping per seed thread
        integer(int32), intent(inout) :: tmp_perm(n_vectors, n_seeds)
        !! Workspace array for sorting per seed thread
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors, n_seeds)
        !! Workspace array for deviation calculations per seed thread
        real(real64), intent(inout) :: tmp_observables(:, :, :)
        !! Workspace matrix for storing observable history per seed thread
        logical, intent(inout) :: tmp_current_mask(n_vectors, n_seeds)
        !! Workspace array tracking active seed ensemble state
        logical, intent(inout) :: tmp_backup_mask(n_vectors, n_seeds)
        !! Workspace array backing up prior state
        logical, intent(out) :: ensemble_matrix(n_vectors, n_seeds)
        !! Output matrix storing raw grown ensemble masks
        integer(int32), intent(out) :: n_ensembles
        !! Total count of extracted raw ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        call set_ok(ierr)

        ! Structural dimension validation
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_dimension_size(n_seeds, ierr)
        if (is_err(ierr)) return

        ! Array range and element validation
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        call validate_all_in_range_real(density_labels, n_vectors, ierr, min=0.0_real64)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        call validate_all_in_range_int(seed_indices, n_seeds, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Scalar parameter validation
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_mad, ierr, min=0.0_real64)
        call validate_in_range_real(alpha_accept, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call obtain_ensembles_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                     kd_indices, density_labels, seed_indices, n_seeds, &
                                     r, alpha_mad, alpha_accept, t_observables, &
                                     tmp_stack, tmp_vicinity_mask, tmp_perm, &
                                     tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                     tmp_backup_mask, ensemble_matrix, n_ensembles, ierr)

    end subroutine obtain_ensembles

    !> Core Implementation for parallel seed ensemble growth.
    pure subroutine obtain_ensembles_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                            kd_indices, density_labels, seed_indices, n_seeds, &
                                            r, alpha_mad, alpha_accept, t_observables, &
                                            tmp_stack, tmp_vicinity_mask, tmp_perm, &
                                            tmp_abs_diff, tmp_observables, tmp_current_mask, &
                                            tmp_backup_mask, ensemble_matrix, n_ensembles, ierr)

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
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors, n_seeds)
        !! Workspace stack for tree traversal per seed thread
        logical, intent(inout) :: tmp_vicinity_mask(n_vectors, n_vectors, n_seeds)
        !! Workspace matrix for neighborhood mapping per seed thread
        integer(int32), intent(inout) :: tmp_perm(n_vectors, n_seeds)
        !! Workspace array for sorting per seed thread
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors, n_seeds)
        !! Workspace array for deviation calculations per seed thread
        real(real64), intent(inout) :: tmp_observables(:, :, :)
        !! Workspace matrix for storing observable history per seed thread
        logical, intent(inout) :: tmp_current_mask(n_vectors, n_seeds)
        !! Workspace array tracking active seed ensemble state per thread
        logical, intent(inout) :: tmp_backup_mask(n_vectors, n_seeds)
        !! Workspace array backing up prior state per thread
        logical, intent(out) :: ensemble_matrix(n_vectors, n_seeds)
        !! Output matrix storing unmerged grown ensemble masks
        integer(int32), intent(out) :: n_ensembles
        !! Total count of extracted raw ensembles (equals n_seeds)
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: i_seed, local_ierr

        call set_ok(ierr)

        ensemble_matrix = .false.

        ! Pure do concurrent outer-loop parallelization over independent seeds
        do concurrent(i_seed=1:n_seeds) &
            shared(vectors, n_dimensions, n_vectors, dimension_order, kd_indices, &
                   density_labels, seed_indices, r, alpha_mad, alpha_accept, &
                   t_observables, tmp_stack, tmp_vicinity_mask, tmp_perm, &
                   tmp_abs_diff, tmp_observables, tmp_current_mask, tmp_backup_mask, ensemble_matrix) &
            local(local_ierr) &
            reduce(max:ierr)

            local_ierr = ERR_OK

            call grow_single_seed_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                         kd_indices, density_labels, seed_indices(i_seed), &
                                         r, alpha_mad, alpha_accept, t_observables, &
                                         tmp_stack(:, :, :, i_seed), tmp_vicinity_mask(:, :, i_seed), &
                                         tmp_perm(:, i_seed), tmp_abs_diff(:, i_seed), &
                                         tmp_observables(:, :, i_seed), tmp_current_mask(:, i_seed), &
                                         tmp_backup_mask(:, i_seed), ensemble_matrix(:, i_seed), local_ierr)

            ierr = max(ierr, local_ierr)
        end do

        n_ensembles = n_seeds

    end subroutine obtain_ensembles_helper

    !> Core Implementation for growing a single seed into an ensemble.
    pure subroutine grow_single_seed_helper(vectors, n_dimensions, n_vectors, dimension_order, &
                                            kd_indices, density_labels, seed_idx, &
                                            r, alpha_mad, alpha_accept, t_observables, &
                                            tmp_stack, tmp_vicinity_mask, &
                                            tmp_perm, tmp_abs_diff, tmp_observables, &
                                            current_mask, backup_mask, out_mask, ierr)

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
        real(real64), intent(in) :: alpha_accept
        !! Acceptance log2 fold-change threshold
        integer(int32), intent(in) :: t_observables
        !! History tracking depth for observables
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Workspace stack for tree traversal
        logical, intent(inout) :: tmp_vicinity_mask(n_vectors, n_vectors)
        !! Workspace matrix for neighborhood mapping
        integer(int32), intent(inout) :: tmp_perm(n_vectors)
        !! Workspace array for sorting
        real(real64), intent(inout) :: tmp_abs_diff(n_vectors)
        !! Workspace array for deviation calculations
        real(real64), intent(inout) :: tmp_observables(:, :)
        !! Workspace matrix for storing observable history
        logical, intent(inout) :: current_mask(n_vectors)
        !! Workspace array tracking active seed ensemble state
        logical, intent(inout) :: backup_mask(n_vectors)
        !! Workspace array backing up prior iteration state
        logical, intent(out) :: out_mask(n_vectors)
        !! Output boolean mask for grown single-seed ensemble
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: iter, prev_count, curr_count, candidate_count
        logical :: is_growing, is_accepted

        call set_ok(ierr)

        current_mask = .false.
        current_mask(seed_idx) = .true.
        tmp_observables = 0.0_real64

        iter = 1_int32
        call compute_ensemble_observable_helper(current_mask, density_labels, &
                                                n_vectors, 0_int32, iter, &
                                                tmp_observables, &
                                                t_observables, ierr)

        is_growing = .true.
        do while (is_growing)
            backup_mask = current_mask
            prev_count = count(current_mask)

            call grow_ensemble_helper(vectors, n_dimensions, n_vectors, current_mask, &
                                      r, dimension_order, kd_indices, density_labels, &
                                      alpha_mad, tmp_stack, &
                                      tmp_vicinity_mask, &
                                      tmp_perm, tmp_abs_diff, ierr)

            curr_count = count(current_mask)
            candidate_count = curr_count - prev_count

            if (candidate_count == 0_int32) then
                is_growing = .false.
                exit
            end if

            iter = iter + 1_int32

            call compute_ensemble_observable_helper(current_mask, density_labels, &
                                                    n_vectors, candidate_count, iter, &
                                                    tmp_observables, &
                                                    t_observables, ierr)

            call accept_ensemble_helper(tmp_observables, iter, &
                                        alpha_accept, is_accepted, ierr)

            if (.not. is_accepted) then
                current_mask = backup_mask
                is_growing = .false.
            end if
        end do

        out_mask = current_mask

    end subroutine grow_single_seed_helper

    !> Allocating Wrapper for transitive set-union ensemble merging.
    subroutine merge_ensembles_alloc(raw_masks, n_vectors, n_seeds, &
                                     min_intersection, merged_matrix, &
                                     n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        logical, intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Input matrix of raw boolean ensemble masks [n_vectors x n_seeds]
        integer(int32), intent(in), optional :: min_intersection
        !! Minimum points in common to trigger a merge (defaults to 1)
        logical, allocatable, intent(out) :: merged_matrix(:, :)
        !! Output matrix storing final unique merged ensemble masks [n_vectors x n_ensembles]
        integer(int32), intent(out) :: n_ensembles
        !! Final count of unique merged ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: actual_min_intersect
        logical, allocatable :: tmp_merged_masks(:, :)
        logical, allocatable :: tmp_active_flag(:)

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

        call merge_ensembles(raw_masks, n_vectors, n_seeds, actual_min_intersect, &
                             tmp_merged_masks, tmp_active_flag, n_ensembles, ierr)

        if (is_err(ierr)) return

        ! Slice final matrix output to actual merged count
        M_ALLOCATE(merged_matrix(n_vectors, n_ensembles))
        if (n_ensembles > 0_int32) then
            merged_matrix(:, 1:n_ensembles) = tmp_merged_masks(:, 1:n_ensembles)
        end if

    end subroutine merge_ensembles_alloc

    !> Validated Entry Point for transitive set-union ensemble merging.
    subroutine merge_ensembles(raw_masks, n_vectors, n_seeds, min_intersection, &
                               merged_masks, tmp_active_flag, n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        logical, intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Input matrix of raw boolean ensemble masks [n_vectors x n_seeds]
        integer(int32), intent(in) :: min_intersection
        !! Minimum points in common to trigger a merge
        logical, intent(out) :: merged_masks(n_vectors, n_seeds)
        !! Workspace matrix storing final merged ensemble masks
        logical, intent(inout) :: tmp_active_flag(n_seeds)
        !! Preallocated workspace tracking active unmerged seeds
        integer(int32), intent(out) :: n_ensembles
        !! Final count of unique merged ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        call set_ok(ierr)

        call validate_dimension_size(n_vectors, ierr)
        call validate_in_range_int(n_seeds, ierr, min=0_int32, max=n_vectors)
        call validate_in_range_int(min_intersection, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        call merge_ensembles_helper(raw_masks, n_vectors, n_seeds, min_intersection, &
                                    merged_masks, tmp_active_flag, n_ensembles, ierr)

    end subroutine merge_ensembles

    !> Core Implementation for pairwise ensemble merging based on set intersection.
    pure subroutine merge_ensembles_helper(raw_masks, n_vectors, n_seeds, &
                                           min_intersection, merged_masks, &
                                           tmp_active_flag, n_ensembles, ierr)

        integer(int32), intent(in) :: n_vectors
        !! Total number of ambient vectors
        integer(int32), intent(in) :: n_seeds
        !! Total number of grown raw ensemble columns
        integer(int32), intent(in) :: min_intersection
        !! Minimum overlapping vectors required to merge two ensembles
        logical, intent(in) :: raw_masks(n_vectors, n_seeds)
        !! Matrix of raw unmerged boolean ensemble masks
        logical, intent(out) :: merged_masks(n_vectors, n_seeds)
        !! Output matrix storing merged boolean ensemble masks
        logical, intent(inout) :: tmp_active_flag(n_seeds)
        !! Preallocated workspace tracking active unmerged seeds
        integer(int32), intent(out) :: n_ensembles
        !! Final count of unique merged ensembles
        integer(int32), intent(out) :: ierr
        !! Error status flag

        integer(int32) :: i, j, shared_count
        logical :: merged_any

        call set_ok(ierr)

        if (n_seeds == 0_int32) then
            n_ensembles = 0_int32
            merged_masks = .false.
            return
        end if

        merged_masks = raw_masks
        tmp_active_flag = .true.

        ! Iterative transitive set-union merging
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

        ! Compact active merged ensemble columns to the front of matrix
        n_ensembles = 0_int32
        do i = 1, n_seeds
            if (tmp_active_flag(i)) then
                n_ensembles = n_ensembles + 1_int32
                if (n_ensembles /= i) then
                    merged_masks(:, n_ensembles) = merged_masks(:, i)
                end if
            end if
        end do

    end subroutine merge_ensembles_helper

end module tox_shatter_cluster_data
