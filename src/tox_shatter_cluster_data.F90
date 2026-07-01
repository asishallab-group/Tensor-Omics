#include "macros.h"

!> Module for managing data structures and geometric calculations
!! For shatter clustering operations within Tensor Omics.
module tox_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, &
                          validate_all_in_range_int, ERR_ALLOC_FAIL, ERR_DIM_MISMATCH, ERR_INVALID_INPUT
    use tox_gene_centroids, only: mean_vector
    use tox_euclidean_distance, only: euclidean_distance
    use f42_utils, only: sort_real_heapsort, calc_percentile
    implicit none
    private

    public :: calculate_density_radius_alloc, calculate_density_radius, calculate_density_radius_helper
    public :: calculate_labels_as_density

contains

    !> Convenience Allocating Wrapper for calculating label density search radius.
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
        !! Error code output (always the final argument).

        ! Local work arrays utilizing mandated tmp_ prefix
        real(real64), allocatable :: tmp_mean_vec(:)
        real(real64), allocatable :: tmp_distances(:)
        integer(int32), allocatable :: tmp_perm(:)

        call set_ok(ierr)

        ! Validate all input dimensions first
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        if (is_err(ierr)) return

        ! Perform structural allocations using the project safety macro
        M_ALLOCATE(tmp_mean_vec(n_dimensions))
        M_ALLOCATE(tmp_distances(n_vectors))
        M_ALLOCATE(tmp_perm(n_vectors))

        ! Forward execution down to Layer 2
        call calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                      tmp_mean_vec, tmp_distances, tmp_perm, &
                                      radius, mean_to_other_vecs_dist_quant, ierr)

    end subroutine calculate_density_radius_alloc

    !> Validated Entry Point for calculating label density search radius.
    subroutine calculate_density_radius(vectors, n_dimensions, n_vectors, &
                                        mean_vec, distances, perm, radius, &
                                        mean_to_other_vecs_dist_quant, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(out) :: mean_vec(n_dimensions)
        !! Preallocated workspace for the mean vector
        real(real64), intent(out) :: distances(n_vectors)
        !! Preallocated workspace for calculated distances
        integer(int32), intent(out) :: perm(n_vectors)
        !! Preallocated workspace for permutation indices
        real(real64), intent(out) :: radius
        !! The resulting density search radius
        real(real64), intent(in), optional :: mean_to_other_vecs_dist_quant
        !! Optional quantile fraction (0.0 to 1.0), defaults to 0.15
        integer(int32), intent(out) :: ierr
        !! Error code output (always the final argument).

        real(real64) :: actual_quant

        call set_ok(ierr)

        ! Validate all input dimensions first
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        if (is_err(ierr)) return

        ! Assign default value using mandated preprocessor macro
        M_DEFAULT_VAL(mean_to_other_vecs_dist_quant, actual_quant, 0.15_real64)

        ! Validate the calculated/assigned range
        call validate_in_range_real(actual_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        ! Delegate directly to pure core implementation layer
        call calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                             mean_vec, distances, perm, actual_quant, &
                                             radius, ierr)

    end subroutine calculate_density_radius

    !> Core Implementation for calculating label density search radius.
    pure subroutine calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                                    mean_vec, distances, perm, &
                                                    quantile_frac, radius, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input data matrix (n_dimensions x n_vectors)
        real(real64), intent(out) :: mean_vec(n_dimensions)
        !! Preallocated workspace for the mean vector
        real(real64), intent(out) :: distances(n_vectors)
        !! Preallocated workspace for calculated distances
        integer(int32), intent(out) :: perm(n_vectors)
        !! Preallocated workspace for permutation indices
        real(real64), intent(in) :: quantile_frac
        !! Quantile fraction verified (0.0 to 1.0)
        real(real64), intent(out) :: radius
        !! The resulting density search radius
        integer(int32), intent(out) :: ierr
        !! Error code output (always the final argument).

        real(real64) :: percentile
        integer(int32) :: i_vec

        call set_ok(ierr)

        ! Element-wise initialization executed in parallel
        do concurrent(i_vec=1:n_vectors) shared(perm)
            perm(i_vec) = i_vec
        end do

        call mean_vector(vectors, n_dimensions, n_vectors, perm, n_vectors, mean_vec, ierr)
        if (is_err(ierr)) return

        ! Multi-threaded distance execution via explicit locality specification
        do concurrent(i_vec=1:n_vectors) shared(vectors, mean_vec, n_dimensions, distances)
            call euclidean_distance(mean_vec, vectors(:, i_vec), n_dimensions, distances(i_vec))
        end do

        call sort_real_heapsort(distances, perm)

        percentile = quantile_frac*100.0_real64

        call calc_percentile(distances, perm, percentile, radius, ierr)

    end subroutine calculate_density_radius_helper

    !> Validated Entry Point for calculating label density distributions.
    !! Inspects input array properties before invoking the parallel KD-Tree search.
    subroutine calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                           dimension_order, kd_indices, label_densities, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Vector dimension (rows)
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors (columns)
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input vectors data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
        !! Label-sphere search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via build_kd_index
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output density tracker matching individual vector slots
        integer(int32), intent(out) :: ierr
        !! Error code tracker (always final parameter)

        call set_ok(ierr) ! ierr = ERR_OK = 0

        ! Validate dimension sizes using core tox modules
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        if (is_err(ierr)) return

        ! Boundary validation check for index sequence
        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        ! Validate search metrics are non-negative
        call validate_in_range_real(r, ierr, min=0.0_real64, max=huge(1.0_real64))
        if (is_err(ierr)) return

        ! 🚀 FIXED CALL ALIGNMENT: Arguments here now perfectly match the helper signature below
        call calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                dimension_order, kd_indices, label_densities, ierr)

    end subroutine calculate_labels_as_density

    !> Core Implementation for calculating label density coordinates.
    !! Traverses the flat KD-tree using an inlined thread-local stack structure.
    pure subroutine calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                       dimension_order, kd_indices, label_densities, ierr)

        integer(int32), intent(in) :: n_dimensions
        !! Vector dimension (rows)
        integer(int32), intent(in) :: n_vectors
        !! Number of vectors (columns)
        real(real64), intent(in) :: vectors(n_dimensions, n_vectors)
        !! Input vectors data matrix (n_dimensions x n_vectors)
        real(real64), intent(in) :: r
        !! Label-sphere search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Dimension split order array tracking the tree structure
        integer(int32), intent(in) :: kd_indices(n_vectors)
        !! KD-tree index sequence array computed via build_kd_index
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output vector storing generated density scalars
        integer(int32), intent(out) :: ierr
        !! Error code tracker (always final parameter)

        ! Traversal state variables declared locally for parallel tracking
        integer(int32) :: i_vec
        integer(int32) :: match_count
        integer(int32) :: stack_top
        integer(int32) :: left_idx
        integer(int32) :: right_idx
        integer(int32) :: mid_idx
        integer(int32) :: current_dim
        integer(int32) :: current_depth
        integer(int32) :: node_point_idx
        real(real64) :: dist
        real(real64) :: axis_delta

        ! Mandated tmp_ prefix applied to this loop-local work array
        integer(int32) :: tmp_stack(3, 64)

        call set_ok(ierr) ! Initialise status tracker

        ! Parallel computation over completely independent loops
        do concurrent(i_vec=1:n_vectors) &
            shared(vectors, n_dimensions, n_vectors, r, dimension_order, kd_indices, label_densities) &
            local(match_count, tmp_stack, stack_top, left_idx, right_idx, mid_idx, &
                  current_dim, current_depth, node_point_idx, dist, axis_delta) ! Strict thread isolation

            match_count = 0
            stack_top = 1

            ! Push flat implicit root parameters onto our thread-local traversal stack
            tmp_stack(1, 1) = 1
            tmp_stack(2, 1) = n_vectors
            tmp_stack(3, 1) = 0

            ! Walk the implicit tree array layout
            do while (stack_top > 0)
                left_idx = tmp_stack(1, stack_top)
                right_idx = tmp_stack(2, stack_top)
                current_depth = tmp_stack(3, stack_top)
                stack_top = stack_top - 1

                if (right_idx < left_idx) cycle

                ! Identify tracking split dimension exact to building step logic
                current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)
                mid_idx = left_idx + (right_idx - left_idx)/2
                node_point_idx = kd_indices(mid_idx)

                ! Assess full distance from target vector to current tree node point
                call euclidean_distance(vectors(:, i_vec), vectors(:, node_point_idx), n_dimensions, dist)
                if (dist <= r) then
                    match_count = match_count + 1
                end if

                ! Calculate coordinate plane offset on split axis
                axis_delta = vectors(current_dim, i_vec) - vectors(current_dim, node_point_idx)

                ! Bounding box evaluation to handle Left, Right, and Overlap plane intersections
                if (axis_delta - r <= 0.0_real64) then
                    if (left_idx <= mid_idx - 1) then
                        stack_top = stack_top + 1
                        tmp_stack(1, stack_top) = left_idx
                        tmp_stack(2, stack_top) = mid_idx - 1
                        tmp_stack(3, stack_top) = current_depth + 1
                    end if
                end if

                if (axis_delta + r >= 0.0_real64) then
                    if (mid_idx + 1 <= right_idx) then
                        stack_top = stack_top + 1
                        tmp_stack(1, stack_top) = mid_idx + 1
                        tmp_stack(2, stack_top) = right_idx
                        tmp_stack(3, stack_top) = current_depth + 1
                    end if
                end if
            end do

            ! Assign counted density cleanly back to array slot
            label_densities(i_vec) = real(match_count, real64)
        end do

    end subroutine calculate_labels_as_density_helper

end module tox_shatter_cluster_data
