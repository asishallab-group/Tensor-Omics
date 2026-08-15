#include "macros.h"
#define CM_KD_STACK_ENTRY_SIZE 3
#define CM_KD_TRAVERSAL_STACK_DEPTH 64

module f42_kd_tree
    use safeguard
    use f42_utils, only: sort_array
    use tox_euclidean_distance, only: euclidean_distance
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: ERR_INVALID_INPUT, ERR_DIM_MISMATCH, ERR_ALLOC_FAIL, &
                          set_ok, set_err, set_err_once, is_ok, is_err, &
                          validate_dimension_size, validate_all_in_range_int, &
                          validate_in_range_real
    implicit none
    private
    public :: build_kd_index, build_spherical_kd, get_kd_point
    public :: vicinity_vectors_alloc, vicinity_vectors_helper, vicinity_vectors_count_alloc, vicinity_vectors_count_helper
    integer(int32), parameter, public :: KD_STACK_ENTRY_SIZE = CM_KD_STACK_ENTRY_SIZE
    integer(int32), parameter, public :: KD_TRAVERSAL_STACK_DEPTH = CM_KD_TRAVERSAL_STACK_DEPTH

contains

    !> Build a k-d tree index using a stack-based, non-recursive approach.
    pure subroutine build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                                   workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: num_points
        !! Number of points
        real(real64), intent(in) :: points(num_dimensions, num_points)
        !! Data points
        integer(int32), intent(in) :: dimension_order(num_dimensions)
        !! Dimension order (by variance)
        integer(int32), intent(out) :: recursion_stack(3, num_points)
        !! Stack for l, r, depth

        integer(int32), intent(out) :: kd_indices(num_points)
        !! Output index array (k-d tree order)
        integer(int32), intent(out) :: workspace(num_points)
        !! Workspace array
        real(real64), intent(out) :: value_buffer(num_points)
        !! Workspace for sorting
        integer(int32), intent(out) :: permutation(num_points)
        !! Workspace for sorting
        integer(int32), intent(out) :: left_stack(num_points)
        !! Workspace for sorting
        integer(int32), intent(out) :: right_stack(num_points)
        !! Workspace for sorting
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: stack_top
        integer(int32) :: left_idx, right_idx, mid_idx, current_dim, current_depth
        integer(int32) :: i

        call set_ok(ierr)

        ! Input validation
        call validate_dimension_size(num_points, ierr)
        if (.not. is_ok(ierr)) return

        call validate_dimension_size(num_dimensions, ierr)
        if (.not. is_ok(ierr)) return

        do i = 1, size(dimension_order)
            if (dimension_order(i) < 1 .or. dimension_order(i) > num_dimensions) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                exit  ! Exit the loop as soon as an invalid value is found
            end if
        end do

        if (.not. is_ok(ierr)) return

        !! Initialize kd_indices to 1:num_points (original indices)
        do i = 1, num_points
            kd_indices(i) = i
        end do

        stack_top = 1
        recursion_stack(1, 1) = 1
        recursion_stack(2, 1) = num_points
        recursion_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = recursion_stack(1, stack_top)
            right_idx = recursion_stack(2, stack_top)
            current_depth = recursion_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx <= left_idx) cycle

            !! Choose split dimension by cycling through dimension_order
            current_dim = dimension_order(mod(current_depth, num_dimensions) + 1)

            !! Find median index
            mid_idx = left_idx + (right_idx - left_idx)/2

            !! Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
            call partial_sort_by_dimension(points, num_points, num_dimensions, kd_indices, left_idx, right_idx, &
                                           current_dim, mid_idx, workspace, value_buffer, permutation, &
                                           left_stack, right_stack, ierr)
            if (.not. is_ok(ierr)) return

            !! Push right and left intervals onto stack
            if (mid_idx < right_idx) then
                stack_top = stack_top + 1
                recursion_stack(1, stack_top) = mid_idx + 1
                recursion_stack(2, stack_top) = right_idx
                recursion_stack(3, stack_top) = current_depth + 1
            end if
            if (left_idx < mid_idx) then
                stack_top = stack_top + 1
                recursion_stack(1, stack_top) = left_idx
                recursion_stack(2, stack_top) = mid_idx - 1
                recursion_stack(3, stack_top) = current_depth + 1
            end if
        end do
    end subroutine build_kd_index

    !> Helper: sorts kd_indices(left_idx:right_idx) by points(dimension, kd_indices(:))
    pure subroutine partial_sort_by_dimension(points, n_points, num_dimensions, kd_indices, left_idx, right_idx, &
                                              dim, mid_idx, workspace, value_buffer, permutation, &
                                              left_stack, right_stack, ierr)
        use f42_utils, only: sort_array
        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: left_idx
        !! Left index of subarray
        integer(int32), intent(in) :: right_idx
        !! Right index of subarray
        integer(int32), intent(in) :: dim
        !! Dimension to sort by
        integer(int32), intent(in) :: mid_idx
        !! Target median index
        integer(int32), intent(in) :: n_points
        !! size of points
        real(real64), intent(in) :: points(num_dimensions, n_points)
        !! Input points array
        integer(int32), intent(out) :: kd_indices(:)
        !! Index array to modify
        integer(int32), intent(out) :: workspace(:)
        !! Workspace array
        real(real64), intent(out) :: value_buffer(:)
        !! Buffer for dimension values
        integer(int32), intent(out) :: permutation(:)
        !! Permutation array
        integer(int32), intent(out) :: left_stack(:)
        !! Stack for left indices
        integer(int32), intent(out) :: right_stack(:)
        !! Stack for right indices
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32) :: subarray_size, i

        call set_ok(ierr)

        ! Input validation
        if (left_idx < 1 .or. right_idx > size(kd_indices) .or. left_idx > right_idx) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (dim < 1 .or. dim > num_dimensions) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (mid_idx < left_idx .or. mid_idx > right_idx) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        subarray_size = right_idx - left_idx + 1
        if (subarray_size <= 1) return

        !! Fill value_buffer with the values of points(dimension, kd_indices(left_idx:right_idx))
        do i = 1, subarray_size
            value_buffer(i) = points(dim, kd_indices(left_idx + i - 1))
            permutation(i) = i
        end do

        call sort_array(value_buffer(1:subarray_size), permutation(1:subarray_size), left_stack, right_stack)

        !! Reorder kd_indices(left_idx:right_idx) according to permutation
        do i = 1, subarray_size
            if (permutation(i) < 1 .or. permutation(i) > subarray_size) then
                ierr = ERR_INVALID_INPUT
                return
            end if
            workspace(i) = kd_indices(left_idx + permutation(i) - 1)
        end do
        do i = 1, subarray_size
            kd_indices(left_idx + i - 1) = workspace(i)
        end do
    end subroutine partial_sort_by_dimension

    !> Build spherical k-d tree index
    pure subroutine build_spherical_kd(vectors, num_dimensions, num_vectors, sphere_indices, &
                                       dimension_order, workspace, value_buffer, permutation, &
                                       left_stack, right_stack, recursion_stack, ierr)

        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: num_vectors
        !! Number of vectors
        real(real64), intent(in) :: vectors(num_dimensions, num_vectors)
        !! Input unit vectors
        integer(int32), intent(out) :: recursion_stack(3, num_vectors)
        !! Stack for recursive calls
        integer(int32), intent(out) :: sphere_indices(num_vectors)
        !! Output index array
        integer(int32), intent(out) :: dimension_order(num_dimensions)
        !! Dimension order
        integer(int32), intent(out) :: workspace(num_vectors)
        !! Workspace array
        real(real64), intent(out) :: value_buffer(num_vectors)
        !! Value buffer
        integer(int32), intent(out) :: permutation(num_vectors)
        !! Permutation array
        integer(int32), intent(out) :: left_stack(num_vectors)
        !! Left stack
        integer(int32), intent(out) :: right_stack(num_vectors)
        !! Right stack
        integer(int32), intent(out) :: ierr
        !! Error code

        call build_kd_index(vectors, num_dimensions, num_vectors, sphere_indices, dimension_order, &
                            workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
    end subroutine build_spherical_kd

    !> Get point from KD index
    pure subroutine get_kd_point(points, kd_indices, position, point_values, ierr)
        real(real64), intent(in) :: points(:, :)
        !! Input points
        integer(int32), intent(in) :: kd_indices(:)
        !! KD index array
        integer(int32), intent(in) :: position
        !! Position in index
        real(real64), intent(out) :: point_values(:)
        !! Output point values
        integer(int32), intent(out) :: ierr
        !! Error code

        call set_ok(ierr)

        ! Input validation
        if (position < 1 .or. position > size(kd_indices)) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        if (kd_indices(position) < 1 .or. kd_indices(position) > size(points, 2)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if

        if (size(point_values) < size(points, 1)) then
            call set_err_once(ierr, ERR_DIM_MISMATCH)
            return
        end if

        point_values = points(:, kd_indices(position))
    end subroutine get_kd_point

    !> Allocating wrapper for finding reference points within a given radius.
    subroutine vicinity_vectors_alloc(query_point, points, num_dimensions, num_points, r, &
                                      dimension_order, kd_indices, vicinity_mask, ierr)

        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: num_points
        !! Total number of points organized in the K-D tree
        real(real64), intent(in) :: query_point(num_dimensions)
        !! The coordinate vector used as the center point of the search window
        real(real64), intent(in) :: points(num_dimensions, num_points)
        !! Ambient data points matrix used to construct the tree layout
        real(real64), intent(in) :: r
        !! Spatial neighborhood search radius threshold
        integer(int32), intent(in) :: dimension_order(num_dimensions)
        !! Sequence array tracking tree split axes by variance
        integer(int32), intent(in) :: kd_indices(num_points)
        !! K-D tree index sequence
        logical, intent(out) :: vicinity_mask(num_points)
        !! Output logical mask flagging points within the search radius boundary
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_stack(:, :)

        call set_ok(ierr)

        call validate_dimension_size(num_dimensions, ierr)
        call validate_dimension_size(num_points, ierr)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_all_in_range_int(dimension_order, num_dimensions, ierr, &
                                       min=1_int32, max=num_dimensions)
        call validate_all_in_range_int(kd_indices, num_points, ierr, &
                                       min=1_int32, max=num_points)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH))

        call vicinity_vectors_helper(query_point, points, num_dimensions, num_points, r, &
                                     dimension_order, kd_indices, tmp_stack, vicinity_mask)

    end subroutine vicinity_vectors_alloc

    !> Finds reference points within a given radius around a query coordinate vector.
    !! Sets elements in a logical mask to .true. for points inside the search sphere.
    pure subroutine vicinity_vectors_helper(query_point, points, num_dimensions, num_points, r, &
                                            dimension_order, kd_indices, tmp_stack, vicinity_mask)

        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: num_points
        !! Total number of points organized in the K-D tree
        real(real64), intent(in) :: query_point(num_dimensions)
        !! The coordinate vector used as the center point of the search window
        real(real64), intent(in) :: points(num_dimensions, num_points)
        !! Ambient data points matrix used to construct the tree layout
        real(real64), intent(in) :: r
        !! Spatial neighborhood search radius threshold
        integer(int32), intent(in) :: dimension_order(num_dimensions)
        !! Sequence array tracking tree split axes by variance
        integer(int32), intent(in) :: kd_indices(num_points)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated private workspace stack for tracking tree-walking frames
        logical, intent(out) :: vicinity_mask(num_points)
        !! Output logical mask flagging points within the search radius boundary

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, node_point_idx
        real(real64)   :: dist, axis_delta

        vicinity_mask = .false.
        stack_top = 1

        ! Push implicit root boundaries onto the stack frame
        tmp_stack(1, 1) = 1
        tmp_stack(2, 1) = num_points
        tmp_stack(3, 1) = 0

        ! Walk the spatial index tree structure
        do while (stack_top > 0)
            left_idx = tmp_stack(1, stack_top)
            right_idx = tmp_stack(2, stack_top)
            current_depth = tmp_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, num_dimensions) + 1)
            mid_idx = left_idx + (right_idx - left_idx)/2
            node_point_idx = kd_indices(mid_idx)

            ! Compute exact spatial distance to candidate point
            call euclidean_distance(query_point, points(:, node_point_idx), num_dimensions, dist)
            if (dist <= r) then
                vicinity_mask(node_point_idx) = .true.
            end if

            ! Calculate coordinate plane delta on the split axis for bounding box checks
            axis_delta = query_point(current_dim) - points(current_dim, node_point_idx)

            ! Assess left-hand child branch viability
            if (axis_delta - r <= 0.0_real64) then
                if (left_idx <= mid_idx - 1) then
                    stack_top = stack_top + 1
                    tmp_stack(1, stack_top) = left_idx
                    tmp_stack(2, stack_top) = mid_idx - 1
                    tmp_stack(3, stack_top) = current_depth + 1
                end if
            end if

            ! Assess right-hand child branch viability
            if (axis_delta + r >= 0.0_real64) then
                if (mid_idx + 1 <= right_idx) then
                    stack_top = stack_top + 1
                    tmp_stack(1, stack_top) = mid_idx + 1
                    tmp_stack(2, stack_top) = right_idx
                    tmp_stack(3, stack_top) = current_depth + 1
                end if
            end if
        end do

    end subroutine vicinity_vectors_helper

    !> Allocating wrapper for counting reference points within a given radius.
    subroutine vicinity_vectors_count_alloc(query_point, points, num_dimensions, num_points, r, &
                                            dimension_order, kd_indices, neighbor_count, ierr)

        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: num_points
        !! Total number of points organized in the K-D tree
        real(real64), intent(in) :: query_point(num_dimensions)
        !! The coordinate vector used as the center point of the search window
        real(real64), intent(in) :: points(num_dimensions, num_points)
        !! Ambient data points matrix used to construct the tree layout
        real(real64), intent(in) :: r
        !! Spatial neighborhood search radius threshold
        integer(int32), intent(in) :: dimension_order(num_dimensions)
        !! Sequence array tracking tree split axes by variance
        integer(int32), intent(in) :: kd_indices(num_points)
        !! K-D tree index sequence
        integer(int32), intent(out) :: neighbor_count
        !! Output scalar count of points within the search radius boundary
        integer(int32), intent(out) :: ierr
        !! Error code

        integer(int32), allocatable :: tmp_stack(:, :)

        call set_ok(ierr)

        call validate_dimension_size(num_dimensions, ierr)
        call validate_dimension_size(num_points, ierr)
        call validate_in_range_real(r, ierr, min=0.0_real64)
        call validate_all_in_range_int(dimension_order, num_dimensions, ierr, &
                                       min=1_int32, max=num_dimensions)
        call validate_all_in_range_int(kd_indices, num_points, ierr, &
                                       min=1_int32, max=num_points)

        if (is_err(ierr)) return

        M_ALLOCATE(tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH))

        call vicinity_vectors_count_helper(query_point, points, num_dimensions, num_points, r, &
                                           dimension_order, kd_indices, tmp_stack, neighbor_count)

    end subroutine vicinity_vectors_count_alloc

    !> Finds the number of reference points within a given radius around a query coordinate vector.
    !! Returns a scalar count of matching points inside the sphere.
    pure subroutine vicinity_vectors_count_helper(query_point, points, num_dimensions, num_points, r, &
                                                  dimension_order, kd_indices, tmp_stack, neighbor_count)

        integer(int32), intent(in) :: num_dimensions
        !! Number of dimensions
        integer(int32), intent(in) :: num_points
        !! Total number of points organized in the K-D tree
        real(real64), intent(in) :: query_point(num_dimensions)
        !! The coordinate vector used as the center point of the search window
        real(real64), intent(in) :: points(num_dimensions, num_points)
        !! Ambient data points matrix used to construct the tree layout
        real(real64), intent(in) :: r
        !! Spatial neighborhood search radius threshold
        integer(int32), intent(in) :: dimension_order(num_dimensions)
        !! Sequence array tracking tree split axes by variance
        integer(int32), intent(in) :: kd_indices(num_points)
        !! KD-tree index sequence array computed via [[f42_kd_tree(module):build_kd_index(subroutine)]].
        integer(int32), intent(inout) :: tmp_stack(CM_KD_STACK_ENTRY_SIZE, CM_KD_TRAVERSAL_STACK_DEPTH)
        !! Preallocated private workspace stack for tracking tree-walking frames
        integer(int32), intent(out) :: neighbor_count
        !! Output scalar count of points within the search radius boundary

        integer(int32) :: stack_top, left_idx, right_idx, mid_idx, current_dim, current_depth, node_point_idx
        real(real64)   :: dist, axis_delta

        neighbor_count = 0_int32
        stack_top = 1

        ! Push implicit root boundaries onto the stack frame
        tmp_stack(1, 1) = 1
        tmp_stack(2, 1) = num_points
        tmp_stack(3, 1) = 0

        ! Walk the spatial index tree structure
        do while (stack_top > 0)
            left_idx = tmp_stack(1, stack_top)
            right_idx = tmp_stack(2, stack_top)
            current_depth = tmp_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx < left_idx) cycle

            current_dim = dimension_order(mod(current_depth, num_dimensions) + 1)
            mid_idx = left_idx + (right_idx - left_idx)/2
            node_point_idx = kd_indices(mid_idx)

            ! Compute exact spatial distance to candidate point
            call euclidean_distance(query_point, points(:, node_point_idx), num_dimensions, dist)
            if (dist <= r) then
                neighbor_count = neighbor_count + 1_int32
            end if

            ! Calculate coordinate plane delta on the split axis for bounding box checks
            axis_delta = query_point(current_dim) - points(current_dim, node_point_idx)

            ! Assess left-hand child branch viability
            if (axis_delta - r <= 0.0_real64) then
                if (left_idx <= mid_idx - 1) then
                    stack_top = stack_top + 1
                    tmp_stack(1, stack_top) = left_idx
                    tmp_stack(2, stack_top) = mid_idx - 1
                    tmp_stack(3, stack_top) = current_depth + 1
                end if
            end if

            ! Assess right-hand child branch viability
            if (axis_delta + r >= 0.0_real64) then
                if (mid_idx + 1 <= right_idx) then
                    stack_top = stack_top + 1
                    tmp_stack(1, stack_top) = mid_idx + 1
                    tmp_stack(2, stack_top) = right_idx
                    tmp_stack(3, stack_top) = current_depth + 1
                end if
            end if
        end do

    end subroutine vicinity_vectors_count_helper

end module f42_kd_tree

!> C interface for building KD index
pure subroutine build_kd_index_C(points, num_dimensions, num_points, kd_indices, dimension_order, &
                                 workspace, value_buffer, permutation, left_stack, right_stack, ierr) &
    bind(C, name="build_kd_index_C")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32
    use f42_kd_tree, only: build_kd_index
    M_USE_NULL_VALIDATION
    implicit none
    !| Input parameters
    integer(c_int), intent(in), target :: num_dimensions
    !| number of input points
    integer(c_int), intent(in), target :: num_points
    !| input points
    real(c_double), intent(in), target :: points(num_dimensions, num_points)
    !| output kd indices
    integer(c_int), intent(in), target :: dimension_order(num_dimensions)
    !| output dimension order
    integer(c_int), intent(out), target :: kd_indices(num_points)
    !| workspace
    integer(c_int), intent(out), target :: workspace(num_points)
    !| value buffer
    real(c_double), intent(out), target :: value_buffer(num_points)
    !| permutation array
    integer(c_int), intent(out), target :: permutation(num_points)
    !| left stack
    integer(c_int), intent(out), target :: left_stack(num_points)
    !| right stack
    integer(c_int), intent(out), target :: right_stack(num_points)
    !| error code
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: recursion_stack(3, num_points)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(num_dimensions)
    M_CHECK_NON_NULL(num_points)
    M_CHECK_NON_NULL(points)
    M_CHECK_NON_NULL(dimension_order)
    M_CHECK_NON_NULL(kd_indices)
    M_CHECK_NON_NULL(workspace)
    M_CHECK_NON_NULL(value_buffer)
    M_CHECK_NON_NULL(permutation)
    M_CHECK_NON_NULL(left_stack)
    M_CHECK_NON_NULL(right_stack)

    ! Call the original implementation
    call build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                        workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
end subroutine build_kd_index_C

!> C interface for building spherical KD index
pure subroutine build_spherical_kd_C(vectors, num_dimensions, num_vectors, sphere_indices, &
                                     dimension_order, workspace, value_buffer, permutation, &
                                     left_stack, right_stack, ierr) bind(C, name="build_spherical_kd_C")
    use iso_c_binding, only: c_int, c_double
    use iso_fortran_env, only: int32
    use f42_kd_tree, only: build_spherical_kd
    M_USE_NULL_VALIDATION
    implicit none
    !| number of input dimensions
    integer(c_int), intent(in), target :: num_dimensions
    !| number of input vectors
    integer(c_int), intent(in), target :: num_vectors
    !| input unit vectors
    real(c_double), intent(in), target :: vectors(num_dimensions, num_vectors)
    !| output sphere indices
    integer(c_int), intent(out), target :: sphere_indices(num_vectors)
    !| output dimension order
    integer(c_int), intent(out), target :: dimension_order(num_dimensions)
    !| workspace
    integer(c_int), intent(out), target :: workspace(num_vectors)
    !| value buffer
    real(c_double), intent(out), target :: value_buffer(num_vectors)
    !| permutation array
    integer(c_int), intent(out), target :: permutation(num_vectors)
    !| left stack
    integer(c_int), intent(out), target :: left_stack(num_vectors)
    !| right stack
    integer(c_int), intent(out), target :: right_stack(num_vectors)
    !| error code
    integer(c_int), intent(out), target :: ierr

    integer(int32) :: recursion_stack(3, num_vectors)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(num_dimensions)
    M_CHECK_NON_NULL(num_vectors)
    M_CHECK_NON_NULL(vectors)
    M_CHECK_NON_NULL(sphere_indices)
    M_CHECK_NON_NULL(dimension_order)
    M_CHECK_NON_NULL(workspace)
    M_CHECK_NON_NULL(value_buffer)
    M_CHECK_NON_NULL(permutation)
    M_CHECK_NON_NULL(left_stack)
    M_CHECK_NON_NULL(right_stack)

    ! Call the original implementation
    call build_spherical_kd(vectors, num_dimensions, num_vectors, sphere_indices, &
                            dimension_order, workspace, value_buffer, permutation, &
                            left_stack, right_stack, recursion_stack, ierr)
end subroutine build_spherical_kd_C
