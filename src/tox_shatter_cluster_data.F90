#include "macros.h"
#define CM_KD_TRAVERSAL_STACK_DEPTH 64

!> Module for managing data structures and geometric calculations
!! For shatter clustering operations within Tensor Omics.
module tox_shatter_cluster_data

    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, set_err, is_err, validate_dimension_size, validate_in_range_real, &
                          validate_all_in_range_int, ERR_ALLOC_FAIL, ERR_DIM_MISMATCH, ERR_INVALID_INPUT, &
                          validate_all_in_range_real
    use tox_gene_centroids, only: mean_vector
    use tox_euclidean_distance, only: euclidean_distance
    use f42_utils, only: sort_real_heapsort, calc_percentile
    implicit none
    private

    public :: calculate_density_radius_alloc, calculate_density_radius, calculate_density_radius_helper
    public :: calculate_labels_as_density_alloc, calculate_labels_as_density, calculate_labels_as_density_helper

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
        if (is_err(ierr)) return

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
        !! Error code

        call set_ok(ierr)

        ! Validate all input structural dimensions and data arrays
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        if (is_err(ierr)) return

        ! Optional parameter range validation
        call validate_in_range_real(mean_to_other_vecs_dist_quant, ierr, min=0.0_real64, max=1.0_real64)
        if (is_err(ierr)) return

        call calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                             mean_vec, distances, perm, &
                                             mean_to_other_vecs_dist_quant, radius, ierr)

    end subroutine calculate_density_radius

    !> Core Implementation for calculating label density search radius.
    pure subroutine calculate_density_radius_helper(vectors, n_dimensions, n_vectors, &
                                                    mean_vec, distances, perm, &
                                                    mean_to_other_vecs_dist_quant, radius, ierr)

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

        do concurrent(i_vec=1:n_vectors) shared(perm)
            perm(i_vec) = i_vec
        end do

        call mean_vector(vectors, n_dimensions, n_vectors, perm, n_vectors, mean_vec, ierr)
        if (is_err(ierr)) return

        do concurrent(i_vec=1:n_vectors) shared(vectors, mean_vec, n_dimensions, distances)
            call euclidean_distance(mean_vec, vectors(:, i_vec), n_dimensions, distances(i_vec))
        end do

        call sort_real_heapsort(distances, perm)

        percentile = actual_quant*100.0_real64

        call calc_percentile(distances, perm, percentile, radius, ierr)

    end subroutine calculate_density_radius_helper

    !> Allocating Wrapper for calculating label density distributions.
    subroutine calculate_labels_as_density_alloc(vectors, n_dimensions, n_vectors, r, &
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
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output density tracker matching individual vector slots
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_stack(:, :, :)

        call set_ok(ierr)

        ! Validate dimension sizes, metrics, and sequences across inputs
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        if (is_err(ierr)) return

        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return

        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        call validate_in_range_real(r, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        M_ALLOCATE(tmp_stack(3, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors))

        call calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                         dimension_order, kd_indices, tmp_stack, label_densities, ierr)

    end subroutine calculate_labels_as_density_alloc

    !> Validated Entry Point for calculating label density distributions.
    subroutine calculate_labels_as_density(vectors, n_dimensions, n_vectors, r, &
                                           dimension_order, kd_indices, tmp_stack, label_densities, ierr)

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
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(inout) :: tmp_stack(3, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Preallocated workspace stack for tree traversal
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output density tracker matching individual vector slots
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Validate dimension sizes, metrics, and sequences across inputs
        call validate_dimension_size(n_dimensions, ierr)
        call validate_dimension_size(n_vectors, ierr)
        call validate_all_in_range_real(vectors, size(vectors, kind=int32), ierr)
        if (is_err(ierr)) return

        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return

        call validate_all_in_range_int(kd_indices, n_vectors, ierr, min=1_int32, max=n_vectors)
        if (is_err(ierr)) return

        call validate_in_range_real(r, ierr, min=0.0_real64)
        if (is_err(ierr)) return

        call calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                dimension_order, kd_indices, tmp_stack, label_densities, ierr)

    end subroutine calculate_labels_as_density

    !> Core Implementation for calculating label density coordinates.
    pure subroutine calculate_labels_as_density_helper(vectors, n_dimensions, n_vectors, r, &
                                                       dimension_order, kd_indices, tmp_stack, label_densities, ierr)

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
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(inout) :: tmp_stack(3, CM_KD_TRAVERSAL_STACK_DEPTH, n_vectors)
        !! Preallocated workspace stack for tree traversal
        real(real64), intent(out) :: label_densities(n_vectors)
        !! Output vector storing generated density scalars
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: i_vec, match_count, stack_top, left_idx, right_idx, &
                          mid_idx, current_dim, current_depth, node_point_idx
        real(real64)   :: dist, axis_delta

        call set_ok(ierr)

        do concurrent(i_vec=1:n_vectors) &
            shared(vectors, n_dimensions, n_vectors, r, dimension_order, kd_indices, tmp_stack, label_densities) &
            local(match_count, stack_top, left_idx, right_idx, mid_idx, &
                  current_dim, current_depth, node_point_idx, dist, axis_delta)

            match_count = 0
            stack_top = 1

            ! Push implicit root parameters
            tmp_stack(1, 1, i_vec) = 1
            tmp_stack(2, 1, i_vec) = n_vectors
            tmp_stack(3, 1, i_vec) = 0

            ! Walk the implicit tree array layout
            do while (stack_top > 0)
                left_idx = tmp_stack(1, stack_top, i_vec)
                right_idx = tmp_stack(2, stack_top, i_vec)
                current_depth = tmp_stack(3, stack_top, i_vec)
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

                ! Bounding box evaluation using the coordinates
                if (axis_delta - r <= 0.0_real64) then
                    if (left_idx <= mid_idx - 1) then
                        stack_top = stack_top + 1
                        tmp_stack(1, stack_top, i_vec) = left_idx
                        tmp_stack(2, stack_top, i_vec) = mid_idx - 1
                        tmp_stack(3, stack_top, i_vec) = current_depth + 1
                    end if
                end if

                if (axis_delta + r >= 0.0_real64) then
                    if (mid_idx + 1 <= right_idx) then
                        stack_top = stack_top + 1
                        tmp_stack(1, stack_top, i_vec) = mid_idx + 1
                        tmp_stack(2, stack_top, i_vec) = right_idx
                        tmp_stack(3, stack_top, i_vec) = current_depth + 1
                    end if
                end if
            end do

            ! Assign counted density cleanly back to array slot
            label_densities(i_vec) = real(match_count, real64)
        end do

    end subroutine calculate_labels_as_density_helper

end module tox_shatter_cluster_data
