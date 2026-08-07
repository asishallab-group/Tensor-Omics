#include <src/macros.h>

!> k-d tree spatial index over fixed-dimensional point sets.
!| Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
!| caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
!| is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
!| index array rather than as linked nodes.
module f42_kd_tree
    use f42_utils, only: sort_array_heapsort, init_perm
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, validate_dimension_size, validate_all_in_range_int, validate_in_range_int, &
                          validate_in_range_real, is_err, set_err, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE

contains

    !> M_EXPORT_C
    !| summary: Build a k-d tree index using a stack-based, non-recursive approach
    !| AUTHOR_AARON_SCHROEDER
    pure subroutine build_kd_index_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), dimension(:, :), allocatable :: tmp_recursion_stack
        integer(int32), dimension(:), allocatable :: tmp_workspace
        real(real64), dimension(:), allocatable :: tmp_value_buffer
        integer(int32), dimension(:), allocatable :: tmp_permutation

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_recursion_stack(3, n_points))
        M_ALLOCATE(tmp_workspace(n_points))
        M_ALLOCATE(tmp_value_buffer(n_points))
        M_ALLOCATE(tmp_permutation(n_points))

        call build_kd_index_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack)
    end subroutine build_kd_index_alloc

    !> AUTHOR_AARON_SCHROEDER
    !| Build a k-d tree index using a stack-based, non-recursive approach.
    pure subroutine build_kd_index(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(3, n_points), intent(out) :: tmp_recursion_stack
            !! Stack for l, r, depth
        integer(int32), dimension(n_points), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(n_points), intent(out) :: tmp_value_buffer
            !! Workspace for sorting
        integer(int32), dimension(n_points), intent(out) :: tmp_permutation
            !! Workspace for sorting
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)

        if (is_err(ierr)) return

        call build_kd_index_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack)
    end subroutine build_kd_index

    !> AUTHOR_AARON_SCHROEDER
    !| (no input validation) Build a k-d tree index using a stack-based, non-recursive approach.
    pure subroutine build_kd_index_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                   tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
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

        !! Initialize kd_indices to 1:n_points (original indices)
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

            !! Choose split dimension by cycling through dimension_order
            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)

            !! Find median index
            mid_idx = left_idx + (right_idx - left_idx)/2

            !TODO optimize: this fully heapsorts the entire [left_idx:right_idx] subrange just to find the median split
            !               point, at every level of the recursion. That makes the overall build O(n log^2 n) instead of
            !               the O(n log n) achievable with a linear-time median-of-medians / quickselect partition (only
            !               the median element needs to be correctly placed; the two sides don't need to be fully sorted).
            !! Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
            call partial_sort_by_dimension_helper(points, n_points, n_dimensions, kd_indices, left_idx, right_idx, &
                                           current_dim, tmp_workspace, tmp_value_buffer, tmp_permutation)

            !! Push right and left intervals onto stack
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
    end subroutine build_kd_index_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Helper: sorts kd_indices(left_idx:right_idx) by points(dimension, kd_indices(:))
    pure subroutine partial_sort_by_dimension(points, n_points, n_dimensions, kd_indices, left_idx, right_idx, &
                                              dim, tmp_workspace, tmp_value_buffer, tmp_permutation, ierr)
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
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_in_range_int(right_idx, ierr, min=1_int32, max=size(kd_indices, kind=int32), arg_pos=6_int32)
        call validate_in_range_int(left_idx, ierr, min=1_int32, max=right_idx, arg_pos=5_int32)
        call validate_in_range_int(dim, ierr, min=1_int32, max=n_dimensions, arg_pos=7_int32)

        if (is_err(ierr)) return

        call partial_sort_by_dimension_helper(points, n_points, n_dimensions, kd_indices, left_idx, right_idx, &
                                              dim, tmp_workspace, tmp_value_buffer, tmp_permutation)
    end subroutine partial_sort_by_dimension

    !> AUTHOR_AARON_SCHROEDER
    !| (no input validation) sorts kd_indices(left_idx:right_idx) by points(dimension, kd_indices(:))
    pure subroutine partial_sort_by_dimension_helper(points, n_points, n_dimensions, kd_indices, left_idx, right_idx, &
                                              dim, tmp_workspace, tmp_value_buffer, tmp_permutation)
        use f42_utils, only: sort_array
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
        do concurrent (i_sliced_element = 1:n_sliced_elements) shared(tmp_value_buffer, tmp_permutation, kd_indices, left_idx, dim)
            tmp_value_buffer(i_sliced_element) = points(dim, kd_indices(left_idx + i_sliced_element - 1))
            tmp_permutation(i_sliced_element) = i_sliced_element
        end do

        call sort_array_heapsort(tmp_value_buffer(1:n_sliced_elements), tmp_permutation(1:n_sliced_elements))

        ! Reorder kd_indices(left_idx:right_idx) according to tmp_permutation
        do concurrent (i_sliced_element = 1:n_sliced_elements) shared(tmp_workspace, kd_indices, left_idx, tmp_permutation)
            tmp_workspace(i_sliced_element) = kd_indices(left_idx + tmp_permutation(i_sliced_element) - 1)
        end do

        do concurrent (i_sliced_element = 1:n_sliced_elements) shared (kd_indices, left_idx, tmp_workspace)
            kd_indices(left_idx + i_sliced_element - 1) = tmp_workspace(i_sliced_element)
        end do
    end subroutine partial_sort_by_dimension_helper

    !> M_EXPORT_C
    !| summary: Build a k-d tree index over points on the unit sphere (unit vectors)
    !| AUTHOR_AARON_SCHROEDER
    !| This is a thin, semantically-named wrapper: partitioning is identical to
    !| [[f42_kd_tree(module):build_kd_index_alloc(subroutine)]] (plain per-axis median splits);
    !| callers are responsible for ensuring `points` are actually unit-normalized beforehand.
    pure subroutine build_spherical_kd_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), intent(out) :: ierr
            !! Error code

        call build_kd_index_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, ierr)
    end subroutine build_spherical_kd_alloc

    !> AUTHOR_AARON_SCHROEDER
    !| Build a k-d tree index over points assumed to lie on the unit sphere (unit vectors), with
    !| caller-provided workspace (no internal allocation). See
    !| [[f42_kd_tree(module):build_spherical_kd_alloc(subroutine)]] for the allocating variant and
    !| a note on why this delegates to the plain (non-spherical) k-d partitioning.
    pure subroutine build_spherical_kd(vectors, n_dimensions, n_vectors, sphere_indices, &
                                       dimension_order, tmp_workspace, tmp_value_buffer, tmp_permutation, &
                                       tmp_recursion_stack, ierr)

        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_vectors
            !! Number of vectors
        real(real64), dimension(n_dimensions, n_vectors), intent(in) :: vectors
            !! Input unit vectors
        integer(int32), dimension(3, n_vectors), intent(out) :: tmp_recursion_stack
            !! Stack for recursive calls
        integer(int32), dimension(n_vectors), intent(out) :: sphere_indices
            !! Output index array
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order
        integer(int32), dimension(n_vectors), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(n_vectors), intent(out) :: tmp_value_buffer
            !! Value buffer
        integer(int32), dimension(n_vectors), intent(out) :: tmp_permutation
            !! Permutation array
        integer(int32), intent(out) :: ierr
            !! Error code

        call build_kd_index(vectors, n_dimensions, n_vectors, sphere_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_recursion_stack, ierr)
    end subroutine build_spherical_kd

    !> AUTHOR_AARON_SCHROEDER
    !| Retrieves the coordinate vector of the point stored at a given position of a built k-d
    !| index, i.e. `point_values = points(:, kd_indices(position))`.
    pure subroutine get_kd_point(points, kd_indices, position, point_values, ierr)
        real(real64), dimension(:, :), intent(in) :: points
            !! Input points
        integer(int32), dimension(:), intent(in) :: kd_indices
            !! KD index array
        integer(int32), intent(in) :: position
            !! Position in index
        real(real64), dimension(:), intent(out) :: point_values
            !! Output point values
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_in_range_int(position, ierr, min=1_int32, max=size(kd_indices, kind=int32), arg_pos=3_int32)
        call validate_in_range_int(kd_indices(position), ierr, min=1_int32, max=size(points, dim=2, kind=int32), arg_pos=2_int32)
        call validate_in_range_int(size(point_values), ierr, min=size(points, dim=1, kind=int32), arg_pos=4_int32)

        if (is_err(ierr)) return

        point_values = points(:, kd_indices(position))
    end subroutine get_kd_point

    !> M_EXPORT_C
    !| summary: Find the k nearest neighbors of a query point in a pre-built k-d tree
    !| AUTHOR_ASIS_HALLAB
    !| Iterative, stack-based traversal with a bounded max-heap of size `k_neighbors` and
    !| splitting-plane pruning (the near side of each split is always explored, the far side
    !| only when its distance to the splitting plane no longer rules out a closer neighbor
    !| than the heap's current worst). O(log k) per heap update instead of a linear scan.
    pure subroutine kd_knn_query_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                       query_point, k_neighbors, neighbors, distances, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index, see [[f42_kd_tree(module):build_kd_index_alloc(subroutine)]]
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        integer(int32), intent(in) :: k_neighbors
            !! Number of neighbors to find
        integer(int32), dimension(k_neighbors), intent(out) :: neighbors
            !! Output: indices of the k nearest neighbors; nearest-to-farthest order is not
            !! guaranteed (max-heap order internally)
        real(real64), dimension(k_neighbors), intent(out) :: distances
            !! Output: Euclidean distances to the k nearest neighbors
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_int(k_neighbors, ierr, min=1_int32, max=n_points, arg_pos=7_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_knn_query_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                 query_point, k_neighbors, tmp_range_stack, neighbors, distances)
    end subroutine kd_knn_query_alloc

    !> AUTHOR_ASIS_HALLAB
    !| Find the k nearest neighbors of a query point in a pre-built k-d tree, with
    !| caller-provided workspace (no internal allocation). See
    !| [[f42_kd_tree(module):kd_knn_query_alloc(subroutine)]] for the allocating variant.
    pure subroutine kd_knn_query(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                 query_point, k_neighbors, tmp_range_stack, neighbors, distances, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        integer(int32), intent(in) :: k_neighbors
            !! Number of neighbors to find
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), dimension(k_neighbors), intent(out) :: neighbors
            !! Output: indices of the k nearest neighbors
        real(real64), dimension(k_neighbors), intent(out) :: distances
            !! Output: Euclidean distances to the k nearest neighbors
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_int(k_neighbors, ierr, min=1_int32, max=n_points, arg_pos=7_int32)

        if (is_err(ierr)) return

        call kd_knn_query_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                 query_point, k_neighbors, tmp_range_stack, neighbors, distances)
    end subroutine kd_knn_query

    !> AUTHOR_ASIS_HALLAB
    !| (no input validation) Find the k nearest neighbors of a query point, via a bounded
    !| max-heap kept directly in `neighbors`/`distances` and splitting-plane pruning. Requires
    !| `k_neighbors <= n_points` (checked by the validated tiers above): every point is then
    !| guaranteed visited before the heap can still have room, so no fallback for "fewer than
    !| k found" is needed.
    pure subroutine kd_knn_query_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, k_neighbors, tmp_range_stack, neighbors, distances)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        integer(int32), intent(in) :: k_neighbors
            !! Number of neighbors to find
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), dimension(k_neighbors), intent(out) :: neighbors
            !! Output: indices of the k nearest neighbors
        real(real64), dimension(k_neighbors), intent(out) :: distances
            !! Output: Euclidean distances to the k nearest neighbors

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, point_idx
        integer(int32) :: found_count, i
        real(real64)   :: dist_sq, diff_val, max_dist_in_heap, axis_dist
        logical        :: explore_left, explore_right

        found_count = 0
        max_dist_in_heap = huge(1.0_real64)

        stack_top = 1
        tmp_range_stack(1, 1) = 1
        tmp_range_stack(2, 1) = n_points
        tmp_range_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = tmp_range_stack(1, stack_top)
            right_idx = tmp_range_stack(2, stack_top)
            current_depth = tmp_range_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)
            mid_idx = left_idx + (right_idx - left_idx)/2
            point_idx = kd_indices(mid_idx)

            dist_sq = 0.0_real64
            do i = 1, n_dimensions
                diff_val = query_point(i) - points(i, point_idx)
                dist_sq = dist_sq + diff_val*diff_val
            end do

            if (found_count < k_neighbors) then
                found_count = found_count + 1
                neighbors(found_count) = point_idx
                distances(found_count) = dist_sq
                call kd_knn_heapify_up(distances, neighbors, found_count)
                max_dist_in_heap = distances(1)
            else if (dist_sq < max_dist_in_heap) then
                neighbors(1) = point_idx
                distances(1) = dist_sq
                call kd_knn_heapify_down(distances, neighbors, k_neighbors)
                max_dist_in_heap = distances(1)
            end if

            axis_dist = query_point(current_dim) - points(current_dim, point_idx)

            explore_left = .true.
            explore_right = .true.
            if (found_count >= k_neighbors) then
                if (axis_dist > 0.0_real64 .and. axis_dist*axis_dist >= max_dist_in_heap) explore_left = .false.
                if (axis_dist < 0.0_real64 .and. axis_dist*axis_dist >= max_dist_in_heap) explore_right = .false.
            end if

            ! Push the far side first, near side last (so the near side pops
            ! first): the heap then fills with genuinely close points early,
            ! which makes pruning on the far side effective sooner.
            if (axis_dist <= 0.0_real64) then
                if (explore_right .and. mid_idx < right_idx) then
                    stack_top = stack_top + 1
                    tmp_range_stack(1, stack_top) = mid_idx + 1
                    tmp_range_stack(2, stack_top) = right_idx
                    tmp_range_stack(3, stack_top) = current_depth + 1
                end if
                if (explore_left .and. left_idx < mid_idx) then
                    stack_top = stack_top + 1
                    tmp_range_stack(1, stack_top) = left_idx
                    tmp_range_stack(2, stack_top) = mid_idx - 1
                    tmp_range_stack(3, stack_top) = current_depth + 1
                end if
            else
                if (explore_left .and. left_idx < mid_idx) then
                    stack_top = stack_top + 1
                    tmp_range_stack(1, stack_top) = left_idx
                    tmp_range_stack(2, stack_top) = mid_idx - 1
                    tmp_range_stack(3, stack_top) = current_depth + 1
                end if
                if (explore_right .and. mid_idx < right_idx) then
                    stack_top = stack_top + 1
                    tmp_range_stack(1, stack_top) = mid_idx + 1
                    tmp_range_stack(2, stack_top) = right_idx
                    tmp_range_stack(3, stack_top) = current_depth + 1
                end if
            end if
        end do

        do i = 1, k_neighbors
            distances(i) = sqrt(distances(i))
        end do
    end subroutine kd_knn_query_helper

    !> AUTHOR_ASIS_HALLAB
    !| Bubble the element at `pos` up a max-heap on `distances`/`neighbors` in lockstep.
    pure subroutine kd_knn_heapify_up(distances, neighbors, pos)
        real(real64), dimension(:), intent(inout) :: distances
            !! Heap keys
        integer(int32), dimension(:), intent(inout) :: neighbors
            !! Heap payload, reordered in lockstep with `distances`
        integer(int32), intent(in) :: pos
            !! Position to bubble up from

        integer(int32) :: current_pos, parent_pos
        real(real64)   :: tmp_dist
        integer(int32) :: tmp_neighbor

        current_pos = pos
        do while (current_pos > 1)
            parent_pos = current_pos/2
            if (distances(current_pos) <= distances(parent_pos)) exit

            tmp_dist = distances(current_pos)
            tmp_neighbor = neighbors(current_pos)
            distances(current_pos) = distances(parent_pos)
            neighbors(current_pos) = neighbors(parent_pos)
            distances(parent_pos) = tmp_dist
            neighbors(parent_pos) = tmp_neighbor

            current_pos = parent_pos
        end do
    end subroutine kd_knn_heapify_up

    !> AUTHOR_ASIS_HALLAB
    !| Bubble the root down a max-heap of size `heap_size` on `distances`/`neighbors`, in
    !| lockstep, after its value has been overwritten.
    pure subroutine kd_knn_heapify_down(distances, neighbors, heap_size)
        real(real64), dimension(:), intent(inout) :: distances
            !! Heap keys
        integer(int32), dimension(:), intent(inout) :: neighbors
            !! Heap payload, reordered in lockstep with `distances`
        integer(int32), intent(in) :: heap_size
            !! Number of elements currently in the heap

        integer(int32) :: current_pos, left_child, right_child, largest_pos
        real(real64)   :: tmp_dist
        integer(int32) :: tmp_neighbor

        current_pos = 1
        do
            left_child = 2*current_pos
            right_child = 2*current_pos + 1
            largest_pos = current_pos

            if (left_child <= heap_size) then
                if (distances(left_child) > distances(largest_pos)) largest_pos = left_child
            end if
            if (right_child <= heap_size) then
                if (distances(right_child) > distances(largest_pos)) largest_pos = right_child
            end if

            if (largest_pos == current_pos) exit

            tmp_dist = distances(current_pos)
            tmp_neighbor = neighbors(current_pos)
            distances(current_pos) = distances(largest_pos)
            neighbors(current_pos) = neighbors(largest_pos)
            distances(largest_pos) = tmp_dist
            neighbors(largest_pos) = tmp_neighbor

            current_pos = largest_pos
        end do
    end subroutine kd_knn_heapify_down

    !> M_EXPORT_C
    !| summary: Mark every point within `radius` of a query point in a pre-built k-d tree
    !| AUTHOR_ASIS_HALLAB
    !| Same iterative, stack-based traversal and splitting-plane pruning as
    !| [[f42_kd_tree(module):kd_knn_query_alloc(subroutine)]] (the near side of each split is
    !| always explored, the far side only when it is still within `radius` of the splitting
    !| plane), with a fixed radius bound instead of a k-nearest-neighbor heap. Compares
    !| squared distances against a precomputed `radius**2` (no `sqrt` per node visited).
    !|
    !| Fits a caller that already does an O(n_points) pass over the result (e.g. merging it
    !| into an existing coverage mask via `.or.`). A caller issuing many independent range
    !| queries per outer step (e.g. one per candidate point in a greedy loop) should use
    !| [[f42_kd_tree(module):kd_range_query_list_alloc(subroutine)]] or
    !| [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]] instead, since the
    !| `in_radius_mask = .false.` reset here costs O(n_points) on every call regardless of how
    !| few points are actually found.
    pure subroutine kd_range_query_mask_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                              query_point, radius, in_radius_mask, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        logical, dimension(n_points), intent(out) :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_real(radius, ierr, min=0.0_real64, arg_pos=7_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_range_query_mask_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, radius, tmp_range_stack, in_radius_mask)
    end subroutine kd_range_query_mask_alloc

    !> AUTHOR_ASIS_HALLAB
    !| Mark every point within `radius` of a query point, with caller-provided workspace (no
    !| internal allocation). See
    !| [[f42_kd_tree(module):kd_range_query_mask_alloc(subroutine)]] for the allocating
    !| variant and when to prefer it over
    !| [[f42_kd_tree(module):kd_range_query_list_alloc(subroutine)]].
    pure subroutine kd_range_query_mask(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, radius, tmp_range_stack, in_radius_mask, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        logical, dimension(n_points), intent(out) :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_real(radius, ierr, min=0.0_real64, arg_pos=7_int32)

        if (is_err(ierr)) return

        call kd_range_query_mask_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, radius, tmp_range_stack, in_radius_mask)
    end subroutine kd_range_query_mask

    !> AUTHOR_ASIS_HALLAB
    !| (no input validation) Mark every point within `radius` of a query point.
    pure subroutine kd_range_query_mask_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                               query_point, radius, tmp_range_stack, in_radius_mask)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        logical, dimension(n_points), intent(out) :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, point_idx, i
        real(real64)   :: dist_sq, radius_sq, axis_dist, diff_val

        radius_sq = radius*radius
        in_radius_mask = .false.

        stack_top = 1
        tmp_range_stack(1, 1) = 1
        tmp_range_stack(2, 1) = n_points
        tmp_range_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = tmp_range_stack(1, stack_top)
            right_idx = tmp_range_stack(2, stack_top)
            current_depth = tmp_range_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)
            mid_idx = left_idx + (right_idx - left_idx)/2
            point_idx = kd_indices(mid_idx)

            dist_sq = 0.0_real64
            do i = 1, n_dimensions
                diff_val = query_point(i) - points(i, point_idx)
                dist_sq = dist_sq + diff_val*diff_val
            end do
            if (dist_sq <= radius_sq) in_radius_mask(point_idx) = .true.

            axis_dist = query_point(current_dim) - points(current_dim, point_idx)

            if (axis_dist <= radius .and. left_idx < mid_idx) then
                stack_top = stack_top + 1
                tmp_range_stack(1, stack_top) = left_idx
                tmp_range_stack(2, stack_top) = mid_idx - 1
                tmp_range_stack(3, stack_top) = current_depth + 1
            end if
            if (axis_dist >= -radius .and. mid_idx < right_idx) then
                stack_top = stack_top + 1
                tmp_range_stack(1, stack_top) = mid_idx + 1
                tmp_range_stack(2, stack_top) = right_idx
                tmp_range_stack(3, stack_top) = current_depth + 1
            end if
        end do
    end subroutine kd_range_query_mask_helper

    !> M_EXPORT_C
    !| summary: List every point within `radius` of a query point in a pre-built k-d tree
    !| AUTHOR_ASIS_HALLAB
    !| Same traversal and pruning as [[f42_kd_tree(module):kd_range_query_mask_alloc(subroutine)]],
    !| but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
    !| instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
    !| in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
    !| needs the count, not the identities, should use
    !| [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]] instead, to skip the
    !| index-buffer writes entirely.
    pure subroutine kd_range_query_list_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                              query_point, radius, neighbors, n_found, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(n_points), intent(out) :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(int32), intent(out) :: n_found
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_real(radius, ierr, min=0.0_real64, arg_pos=7_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_range_query_list_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, radius, tmp_range_stack, neighbors, n_found)
    end subroutine kd_range_query_list_alloc

    !> AUTHOR_ASIS_HALLAB
    !| List every point within `radius` of a query point, with caller-provided workspace (no
    !| internal allocation). See
    !| [[f42_kd_tree(module):kd_range_query_list_alloc(subroutine)]] for the allocating
    !| variant.
    pure subroutine kd_range_query_list(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, radius, tmp_range_stack, neighbors, n_found, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), dimension(n_points), intent(out) :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(int32), intent(out) :: n_found
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_real(radius, ierr, min=0.0_real64, arg_pos=7_int32)

        if (is_err(ierr)) return

        call kd_range_query_list_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                        query_point, radius, tmp_range_stack, neighbors, n_found)
    end subroutine kd_range_query_list

    !> AUTHOR_ASIS_HALLAB
    !| (no input validation) List every point within `radius` of a query point.
    pure subroutine kd_range_query_list_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                               query_point, radius, tmp_range_stack, neighbors, n_found)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), dimension(n_points), intent(out) :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(int32), intent(out) :: n_found
            !! Output: number of points within `radius`

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, point_idx, i
        real(real64)   :: dist_sq, radius_sq, axis_dist, diff_val

        radius_sq = radius*radius
        n_found = 0

        stack_top = 1
        tmp_range_stack(1, 1) = 1
        tmp_range_stack(2, 1) = n_points
        tmp_range_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = tmp_range_stack(1, stack_top)
            right_idx = tmp_range_stack(2, stack_top)
            current_depth = tmp_range_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)
            mid_idx = left_idx + (right_idx - left_idx)/2
            point_idx = kd_indices(mid_idx)

            dist_sq = 0.0_real64
            do i = 1, n_dimensions
                diff_val = query_point(i) - points(i, point_idx)
                dist_sq = dist_sq + diff_val*diff_val
            end do
            if (dist_sq <= radius_sq) then
                n_found = n_found + 1
                neighbors(n_found) = point_idx
            end if

            axis_dist = query_point(current_dim) - points(current_dim, point_idx)

            if (axis_dist <= radius .and. left_idx < mid_idx) then
                stack_top = stack_top + 1
                tmp_range_stack(1, stack_top) = left_idx
                tmp_range_stack(2, stack_top) = mid_idx - 1
                tmp_range_stack(3, stack_top) = current_depth + 1
            end if
            if (axis_dist >= -radius .and. mid_idx < right_idx) then
                stack_top = stack_top + 1
                tmp_range_stack(1, stack_top) = mid_idx + 1
                tmp_range_stack(2, stack_top) = right_idx
                tmp_range_stack(3, stack_top) = current_depth + 1
            end if
        end do
    end subroutine kd_range_query_list_helper

    !> M_EXPORT_C
    !| summary: Count the points within `radius` of a query point in a pre-built k-d tree
    !| AUTHOR_ASIS_HALLAB
    !| Same traversal and pruning as
    !| [[f42_kd_tree(module):kd_range_query_mask_alloc(subroutine)]], but writes no index
    !| buffer at all -- only a scalar count. The right choice when the identities of the
    !| points found are never needed, e.g. a per-point local-density label.
    pure subroutine kd_range_query_count_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                               query_point, radius, neighbor_count, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), intent(out) :: neighbor_count
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_real(radius, ierr, min=0.0_real64, arg_pos=7_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_range_query_count_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                         query_point, radius, tmp_range_stack, neighbor_count)
    end subroutine kd_range_query_count_alloc

    !> AUTHOR_ASIS_HALLAB
    !| Count the points within `radius` of a query point, with caller-provided workspace (no
    !| internal allocation). See
    !| [[f42_kd_tree(module):kd_range_query_count_alloc(subroutine)]] for the allocating
    !| variant.
    pure subroutine kd_range_query_count(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                         query_point, radius, tmp_range_stack, neighbor_count, ierr)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), intent(out) :: neighbor_count
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, min=1_int32, max=n_points, arg_pos=4_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, min=1_int32, max=n_dimensions, arg_pos=5_int32)
        call validate_in_range_real(radius, ierr, min=0.0_real64, arg_pos=7_int32)

        if (is_err(ierr)) return

        call kd_range_query_count_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                         query_point, radius, tmp_range_stack, neighbor_count)
    end subroutine kd_range_query_count

    !> AUTHOR_ASIS_HALLAB
    !| (no input validation) Count the points within `radius` of a query point.
    pure subroutine kd_range_query_count_helper(points, n_dimensions, n_points, kd_indices, dimension_order, &
                                                query_point, radius, tmp_range_stack, neighbor_count)
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), intent(out) :: neighbor_count
            !! Output: number of points within `radius`

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, point_idx, i
        real(real64)   :: dist_sq, radius_sq, axis_dist, diff_val

        radius_sq = radius*radius
        neighbor_count = 0

        stack_top = 1
        tmp_range_stack(1, 1) = 1
        tmp_range_stack(2, 1) = n_points
        tmp_range_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = tmp_range_stack(1, stack_top)
            right_idx = tmp_range_stack(2, stack_top)
            current_depth = tmp_range_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, n_dimensions) + 1)
            mid_idx = left_idx + (right_idx - left_idx)/2
            point_idx = kd_indices(mid_idx)

            dist_sq = 0.0_real64
            do i = 1, n_dimensions
                diff_val = query_point(i) - points(i, point_idx)
                dist_sq = dist_sq + diff_val*diff_val
            end do
            if (dist_sq <= radius_sq) neighbor_count = neighbor_count + 1

            axis_dist = query_point(current_dim) - points(current_dim, point_idx)

            if (axis_dist <= radius .and. left_idx < mid_idx) then
                stack_top = stack_top + 1
                tmp_range_stack(1, stack_top) = left_idx
                tmp_range_stack(2, stack_top) = mid_idx - 1
                tmp_range_stack(3, stack_top) = current_depth + 1
            end if
            if (axis_dist >= -radius .and. mid_idx < right_idx) then
                stack_top = stack_top + 1
                tmp_range_stack(1, stack_top) = mid_idx + 1
                tmp_range_stack(2, stack_top) = right_idx
                tmp_range_stack(3, stack_top) = current_depth + 1
            end if
        end do
    end subroutine kd_range_query_count_helper

end module f42_kd_tree
