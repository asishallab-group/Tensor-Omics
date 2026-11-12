#include <src/macros.h>

#define CM_KD_STACK_ENTRY_SIZE 3
#define CM_KD_TRAVERSAL_STACK_DEPTH 64

!> k-d tree spatial index over fixed-dimensional point sets.
!| Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
!| caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
!| is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
!| index array rather than as linked nodes.
module f42_kd_tree_impl
    use f42_sort_impl, only: sort_array_heapsort, init_perm
    use tox_euclidean_distance_impl, only: euclidean_distance_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, validate_in_range_int, is_err
    M_IMPLICIT_NONE

    integer(int32), parameter :: KD_STACK_ENTRY_SIZE = CM_KD_STACK_ENTRY_SIZE
    integer(int32), parameter :: KD_TRAVERSAL_STACK_DEPTH = CM_KD_TRAVERSAL_STACK_DEPTH

contains

    !> summary: Build a k-d tree index using a stack-based, non-recursive approach
    !| AUTHOR_AARON_SCHROEDER
    pure subroutine build_kd_index_impl(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), dimension(3, n_points), intent(out) :: tmp_recursion_stack
            !! Stack for l, r, depth

        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(n_points), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(n_points), intent(out) :: tmp_value_buffer
            !! Workspace for sorting
        integer(int32), dimension(n_points), intent(out) :: tmp_permutation
            !! Workspace for sorting

        integer(int32) :: stack_top
        integer(int32) :: left_idx, right_idx, mid_idx, current_dim, current_depth

        ! Initialize kd_indices to 1:n_points (original indices)
        call init_perm(kd_indices)

        stack_top = 1
        tmp_recursion_stack(1, 1) = 1
        tmp_recursion_stack(2, 1) = n_points
        tmp_recursion_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = tmp_recursion_stack(1, stack_top)
            right_idx = tmp_recursion_stack(2, stack_top)
            current_depth = tmp_recursion_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx <= left_idx) cycle

            ! Choose split dimension by cycling through dimension_order
            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)

            ! Find median index
            mid_idx = left_idx + (right_idx - left_idx)/2

            !TODO optimize: this fully heapsorts the entire [left_idx:right_idx] subrange just to find the median split
            !               point, at every level of the recursion. That makes the overall build O(n log^2 n) instead of
            !               the O(n log n) achievable with a linear-time median-of-medians / quickselect partition (only
            !               the median element needs to be correctly placed; the two sides don't need to be fully sorted).
            ! Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
            call partial_sort_by_dimension_helper(points, n_points, n_dimensions, kd_indices, left_idx, right_idx, &
                                                  current_dim, tmp_workspace, tmp_value_buffer, tmp_permutation)

            ! Push right and left intervals onto stack
            if (mid_idx < right_idx) then
                stack_top = stack_top + 1
                tmp_recursion_stack(1, stack_top) = mid_idx + 1
                tmp_recursion_stack(2, stack_top) = right_idx
                tmp_recursion_stack(3, stack_top) = current_depth + 1
            end if
            if (left_idx < mid_idx) then
                stack_top = stack_top + 1
                tmp_recursion_stack(1, stack_top) = left_idx
                tmp_recursion_stack(2, stack_top) = mid_idx - 1
                tmp_recursion_stack(3, stack_top) = current_depth + 1
            end if
        end do
    end subroutine build_kd_index_impl

    !> AUTHOR_AARON_SCHROEDER
    !| Sorts kd_indices(left_idx:right_idx) by points(dimension, kd_indices(:)). Internal to the
    !| index build, which has already validated everything this reads.
    pure subroutine partial_sort_by_dimension_helper(points, n_points, n_dimensions, kd_indices, left_idx, right_idx, &
                                                     dim, tmp_workspace, tmp_value_buffer, tmp_permutation)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: left_idx
            !! Left index of subarray
        integer(int32), intent(in) :: right_idx
            !! Right index of subarray
        integer(int32), intent(in) :: dim
            !! Dimension to sort by
        integer(int32), intent(in) :: n_points
            !! size of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Input points array
        integer(int32), dimension(:), intent(out) :: kd_indices
            !! Index array to modify
        integer(int32), dimension(:), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(:), intent(out) :: tmp_value_buffer
            !! Buffer for dimension values
        integer(int32), dimension(:), intent(out) :: tmp_permutation
            !! Permutation array

        integer(int32) :: n_sliced_elements, i_sliced_element

        n_sliced_elements = right_idx - left_idx + 1
        if (n_sliced_elements <= 1) return

        ! Fill tmp_value_buffer with the values of points(dimension, kd_indices(left_idx:right_idx))
        do concurrent(i_sliced_element=1:n_sliced_elements) shared(tmp_value_buffer, tmp_permutation, kd_indices, left_idx, dim)
            tmp_value_buffer(i_sliced_element) = points(dim, kd_indices(left_idx + i_sliced_element - 1))
            tmp_permutation(i_sliced_element) = i_sliced_element
        end do

        call sort_array_heapsort(tmp_value_buffer(1:n_sliced_elements), tmp_permutation(1:n_sliced_elements))

        ! Reorder kd_indices(left_idx:right_idx) according to tmp_permutation
        do concurrent(i_sliced_element=1:n_sliced_elements) shared(tmp_workspace, kd_indices, left_idx, tmp_permutation)
            tmp_workspace(i_sliced_element) = kd_indices(left_idx + tmp_permutation(i_sliced_element) - 1)
        end do

        do concurrent(i_sliced_element=1:n_sliced_elements) shared(kd_indices, left_idx, tmp_workspace)
            kd_indices(left_idx + i_sliced_element - 1) = tmp_workspace(i_sliced_element)
        end do
    end subroutine partial_sort_by_dimension_helper

    !> summary: Build a k-d tree index over points on the unit sphere (unit vectors)
    !| AUTHOR_AARON_SCHROEDER
    !| This is a thin, semantically-named wrapper: partitioning is identical to
    !| [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]] (plain per-axis median splits);
    !| callers are responsible for ensuring `points` are actually unit-normalized beforehand.
    pure subroutine build_spherical_kd_impl(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                            tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(n_points), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(n_points), intent(out) :: tmp_value_buffer
            !! Workspace for sorting
        integer(int32), dimension(n_points), intent(out) :: tmp_permutation
            !! Workspace for sorting
        integer(int32), dimension(3, n_points), intent(out) :: tmp_recursion_stack
            !! Stack for l, r, depth

        call build_kd_index_impl(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                 tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack)
    end subroutine build_spherical_kd_impl

    !> AUTHOR_AARON_SCHROEDER
    !| Retrieves the coordinate vector of the point stored at a given position of a built k-d
    !| index, i.e. `point_values = points(:, kd_indices(position))`.
    pure subroutine get_kd_point(points, kd_indices, position, point_values, ierr)
        real(real64), dimension(:, :), intent(in) :: points
            !! Input points
        integer(int32), dimension(:), intent(in) :: kd_indices
            !! KD index array
            !! DM_MIN(1_int32)
            !! DM_MAX(size(points, dim=2, kind=int32))
        integer(int32), intent(in) :: position
            !! Position in index
            !! DM_MIN(1_int32)
            !! DM_MAX(size(kd_indices, kind=int32))
        real(real64), dimension(size(points, dim=1, kind=int32)), intent(out) :: point_values
            !! Output point values
        integer(int32), intent(out) :: ierr
            !! Error code

        point_values = points(:, kd_indices(position))
    end subroutine get_kd_point

    !> summary: Find reference points within a radius around a query point.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine vicinity_vectors_impl(query_point, points, n_dimensions, n_points, r, &
                                          dimension_order, kd_indices, tmp_stack, vicinity_mask)

        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Total number of points organized in the k-d tree
        real(real64), intent(in) :: query_point(n_dimensions)
            !! Coordinate vector used as the center of the search
        real(real64), intent(in) :: points(n_dimensions, n_points)
            !! Ambient point matrix
        real(real64), intent(in) :: r
            !! Search radius
            !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! Sequence of k-d tree split dimensions
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: kd_indices(n_points)
            !! K-d tree index sequence
            !! DM_MIN(1_int32)
            !! DM_MAX(n_points)
        integer(int32), intent(inout) :: tmp_stack( &
                                         CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH)
            !! Preallocated k-d tree traversal stack
        logical(c_bool), intent(out) :: vicinity_mask(n_points)
            !! Mask indicating points within the search radius

        call vicinity_vectors_helper(query_point, points, n_dimensions, n_points, r, &
                                     dimension_order, kd_indices, tmp_stack, vicinity_mask)

    end subroutine vicinity_vectors_impl

    !> summary: Count reference points within a radius around a query point.
    !| AUTHOR_SALIH_ALBAYRAK
    pure subroutine vicinity_vectors_count_impl(query_point, points, n_dimensions, n_points, r, &
                                                dimension_order, kd_indices, tmp_stack, n_neighbors)

        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Total number of points organized in the k-d tree
        real(real64), intent(in) :: query_point(n_dimensions)
            !! Coordinate vector used as the center of the search
        real(real64), intent(in) :: points(n_dimensions, n_points)
            !! Ambient point matrix
        real(real64), intent(in) :: r
            !! Search radius
            !! DM_MIN(0.0_real64)
        integer(int32), intent(in) :: dimension_order(n_dimensions)
            !! Sequence of k-d tree split dimensions
            !! DM_MIN(1_int32)
            !! DM_MAX(n_dimensions)
        integer(int32), intent(in) :: kd_indices(n_points)
            !! K-d tree index sequence
            !! DM_MIN(1_int32)
            !! DM_MAX(n_points)
        integer(int32), intent(inout) :: tmp_stack( &
                                         CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH)
            !! Preallocated k-d tree traversal stack
        integer(int32), intent(out) :: n_neighbors
            !! Number of points within the search radius

        call vicinity_vectors_count_helper(query_point, points, n_dimensions, n_points, r, &
                                           dimension_order, kd_indices, tmp_stack, n_neighbors)

    end subroutine vicinity_vectors_count_impl

    !> AUTHOR_SALIH_ALBAYRAK
    !| Find reference points within a radius around a query point.
    pure subroutine vicinity_vectors_helper(query_point, points, n_dimensions, n_points, r, &
                                            dimension_order, kd_indices, tmp_stack, vicinity_mask)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_points
        !! Total number of points organized in the k-d tree
        real(real64), intent(in) :: query_point(n_dimensions)
        !! Coordinate vector used as the center of the search
        real(real64), intent(in) :: points(n_dimensions, n_points)
        !! Ambient point matrix
        real(real64), intent(in) :: r
        !! Search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence of k-d tree split dimensions
        integer(int32), intent(in) :: kd_indices(n_points)
        !! K-d tree index sequence
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, &
                                                   CM_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated k-d tree traversal stack
        logical(c_bool), intent(out) :: vicinity_mask(n_points)
        !! Mask indicating points within the search radius

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, point_idx
        real(real64) :: distance, axis_delta

        vicinity_mask = .false.
        stack_top = 1_int32

        tmp_stack(1, 1) = 1_int32
        tmp_stack(2, 1) = n_points
        tmp_stack(3, 1) = 0_int32

        do while (stack_top > 0_int32)
            left_idx = tmp_stack(1, stack_top)
            right_idx = tmp_stack(2, stack_top)
            current_depth = tmp_stack(3, stack_top)
            stack_top = stack_top - 1_int32

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1_int32)
            mid_idx = left_idx + (right_idx - left_idx)/2_int32
            point_idx = kd_indices(mid_idx)

            call euclidean_distance_impl(query_point, points(:, point_idx), &
                                          n_dimensions, distance)

            if (distance <= r) vicinity_mask(point_idx) = .true.

            axis_delta = query_point(current_dim) - points(current_dim, point_idx)

            if (axis_delta - r <= 0.0_real64) then
                if (left_idx <= mid_idx - 1_int32) then
                    stack_top = stack_top + 1_int32
                    tmp_stack(1, stack_top) = left_idx
                    tmp_stack(2, stack_top) = mid_idx - 1_int32
                    tmp_stack(3, stack_top) = current_depth + 1_int32
                end if
            end if

            if (axis_delta + r >= 0.0_real64) then
                if (mid_idx + 1_int32 <= right_idx) then
                    stack_top = stack_top + 1_int32
                    tmp_stack(1, stack_top) = mid_idx + 1_int32
                    tmp_stack(2, stack_top) = right_idx
                    tmp_stack(3, stack_top) = current_depth + 1_int32
                end if
            end if
        end do

    end subroutine vicinity_vectors_helper

    !> AUTHOR_SALIH_ALBAYRAK
    !| Count reference points within a radius around a query point.
    pure subroutine vicinity_vectors_count_helper(query_point, points, n_dimensions, n_points, r, &
                                                  dimension_order, kd_indices, tmp_stack, n_neighbors)

        integer(int32), intent(in) :: n_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: n_points
        !! Total number of points organized in the k-d tree
        real(real64), intent(in) :: query_point(n_dimensions)
        !! Coordinate vector used as the center of the search
        real(real64), intent(in) :: points(n_dimensions, n_points)
        !! Ambient point matrix
        real(real64), intent(in) :: r
        !! Search radius
        integer(int32), intent(in) :: dimension_order(n_dimensions)
        !! Sequence of k-d tree split dimensions
        integer(int32), intent(in) :: kd_indices(n_points)
        !! K-d tree index sequence
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, &
                                                   CM_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated k-d tree traversal stack
        integer(int32), intent(out) :: n_neighbors
        !! Number of points within the search radius

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, point_idx
        real(real64) :: distance, axis_delta

        n_neighbors = 0_int32
        stack_top = 1_int32

        tmp_stack(1, 1) = 1_int32
        tmp_stack(2, 1) = n_points
        tmp_stack(3, 1) = 0_int32

        do while (stack_top > 0_int32)
            left_idx = tmp_stack(1, stack_top)
            right_idx = tmp_stack(2, stack_top)
            current_depth = tmp_stack(3, stack_top)
            stack_top = stack_top - 1_int32

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1_int32)
            mid_idx = left_idx + (right_idx - left_idx)/2_int32
            point_idx = kd_indices(mid_idx)

            call euclidean_distance_impl(query_point, points(:, point_idx), &
                                          n_dimensions, distance)

            if (distance <= r) n_neighbors = n_neighbors + 1_int32

            axis_delta = query_point(current_dim) - points(current_dim, point_idx)

            if (axis_delta - r <= 0.0_real64) then
                if (left_idx <= mid_idx - 1_int32) then
                    stack_top = stack_top + 1_int32
                    tmp_stack(1, stack_top) = left_idx
                    tmp_stack(2, stack_top) = mid_idx - 1_int32
                    tmp_stack(3, stack_top) = current_depth + 1_int32
                end if
            end if

            if (axis_delta + r >= 0.0_real64) then
                if (mid_idx + 1_int32 <= right_idx) then
                    stack_top = stack_top + 1_int32
                    tmp_stack(1, stack_top) = mid_idx + 1_int32
                    tmp_stack(2, stack_top) = right_idx
                    tmp_stack(3, stack_top) = current_depth + 1_int32
                end if
            end if
        end do

    end subroutine vicinity_vectors_count_helper

end module f42_kd_tree_impl
