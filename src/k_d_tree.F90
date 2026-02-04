module kd_tree
    use safeguard
    use f42_utils, only: sort_array
    use iso_fortran_env, only: int32, real64
    use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, ERR_EMPTY_INPUT, ERR_DIM_MISMATCH, ERR_SIZE_MISMATCH, set_ok, set_err_once, is_ok, validate_dimension_size
    implicit none
    private
    public :: build_kd_index, build_spherical_kd, get_kd_point, kd_knn_query

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
        if(.not. is_ok(ierr)) return
        
        call validate_dimension_size(num_dimensions, ierr)
        if(.not. is_ok(ierr)) return
        
        do i = 1, size(dimension_order)
            if (dimension_order(i) < 1 .or. dimension_order(i) > num_dimensions) then
                call set_err_once(ierr, ERR_INVALID_INPUT)
                exit  ! Exit the loop as soon as an invalid value is found
            end if
        end do

        if(.not. is_ok(ierr)) return

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
            mid_idx = left_idx + (right_idx - left_idx) / 2

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

    !> Efficient k-NN query on pre-built k-d tree
    !! Finds the k nearest neighbors to a query point using the k-d tree index
    !! Uses the same range-based traversal as build_kd_index for correctness
    pure subroutine kd_knn_query(points, kd_indices, num_dimensions, num_points, dimension_order, &
                        query_point, k_neighbors, neighbors, distances, workspace, ierr)
        ! Input and output arguments 
        integer(int32), intent(in) :: num_dimensions        !! Number of dimensions
        integer(int32), intent(in) :: num_points            !! Number of points in dataset
        integer(int32), intent(in) :: k_neighbors           !! Number of neighbors to find
        real(real64), intent(in) :: points(num_dimensions, num_points)    !! Original points dataset
        integer(int32), intent(in) :: kd_indices(num_points)              !! Pre-built k-d tree indices
        integer(int32), intent(in) :: dimension_order(num_dimensions)     !! Dimension order from tree build
        real(real64), intent(in) :: query_point(num_dimensions)           !! Query point coordinates
        integer(int32), intent(out) :: neighbors(k_neighbors)             !! Output: indices of k nearest neighbors
        real(real64), intent(out) :: distances(k_neighbors)               !! Output: distances to k nearest neighbors
        real(real64), intent(inout) :: workspace(num_dimensions)          !! Workspace for distance calculations
        integer(int32), intent(out) :: ierr                               !! Error code

        ! Internal variables
        integer(int32) :: current_best_count
        integer(int32) :: point_idx, current_dim, i, j
        integer(int32) :: left_idx, right_idx, mid_idx, current_depth
        real(real64) :: current_dist, max_dist_in_heap, axis_dist
        real(real64) :: diff_val
        logical :: should_explore_left, should_explore_right

        ! Stack for tree traversal
        integer(int32) :: stack_pos
        integer(int32) :: range_stack(3, num_points)    ! [left_idx, right_idx, depth]

        call set_ok(ierr)
        ! Input validation
        if (k_neighbors < 1 .or. k_neighbors > num_points) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if
        if (num_dimensions < 1 .or. num_points < 1) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        ! Initialize neighbor list as empty
        current_best_count = 0
        max_dist_in_heap = huge(1.0_real64)
        ! print *, "Initialized max_dist_in_heap to:", max_dist_in_heap

        ! Initialize stack
        stack_pos = 1
        range_stack(1, 1) = 1             ! left_idx = 1
        range_stack(2, 1) = num_points    ! right_idx = num_points 
        range_stack(3, 1) = 0             ! depth = 0

        ! Iterative tree traversal
        do while (stack_pos > 0)
            ! Pop from stack
            left_idx = range_stack(1, stack_pos)
            right_idx = range_stack(2, stack_pos)
            current_depth = range_stack(3, stack_pos)
            stack_pos = stack_pos - 1

            ! print *, "Popped from stack: left_idx=", left_idx, "right_idx=", right_idx, "current_depth=", current_depth
            ! Skip only empty ranges (right < left). 
            if (right_idx < left_idx) cycle
            ! Calculate the splitting dimension for the current level
            current_dim = dimension_order(mod(current_depth, num_dimensions) + 1)
            ! print *, "Current split dimension:", current_dim

            if (right_idx > left_idx) then
                ! Find the median index
                mid_idx = left_idx + (right_idx - left_idx) / 2
            else
                ! leaf node case
                mid_idx = left_idx
            end if

            ! The point to evaluate is in kd_indices(mid_idx)
            point_idx = kd_indices(mid_idx)
            ! print *, "Median index (in kd_indices):", mid_idx, " Original point index:", point_idx
            ! Calculate the squared distance to the current point
            current_dist = 0.0_real64
            do i = 1, num_dimensions
                diff_val = query_point(i) - points(i, point_idx)
                current_dist = current_dist + diff_val * diff_val
            end do

            ! print *, "Distance to median point:", current_dist
            ! Update the k-NN list (Max-Heap)

            if (current_best_count < k_neighbors) then
                ! There is still space in the list - simple insertion
                current_best_count = current_best_count + 1
                neighbors(current_best_count) = point_idx
                distances(current_best_count) = current_dist
                ! print *, "Added neighbor:", point_idx, "with distance:", current_dist
                ! Maintain max-heap property
                call max_heapify_up(distances, neighbors, current_best_count)
                max_dist_in_heap = distances(1)    ! The root is the maximum
                ! print *, "Updated max_dist_in_heap to:", max_dist_in_heap

            else if (current_dist < max_dist_in_heap) then
                ! Replace the farthest element (root) with the new neighbor
                neighbors(1) = point_idx
                distances(1) = current_dist
                ! print *, "Replaced farthest neighbor with:", point_idx, "distance:", current_dist
                ! Restore max-heap property
                call max_heapify_down(distances, neighbors, k_neighbors)
                max_dist_in_heap = distances(1)    ! New root after heapify
                ! print *, "Updated max_dist_in_heap to:", max_dist_in_heap
            end if

            ! Pruning and subtree exploration only applies to Internal Nodes
            if (right_idx > left_idx) then 
                ! Decide which subtrees to explore (Pruning)
                axis_dist = query_point(current_dim) - points(current_dim, point_idx)
                ! print *, "Axis distance:", axis_dist
                should_explore_left = .true.
                should_explore_right = .true.
                ! Pruning: do not explore subtrees farther than the current worst neighbor
                if (current_best_count >= k_neighbors) then
                    ! print *, "Pruning check: axis_dist=", axis_dist, "axis_dist^2=", axis_dist * axis_dist, "max_dist_in_heap=", max_dist_in_heap
                    ! Check the left subtree
                    if (axis_dist > 0.0_real64 .and. axis_dist * axis_dist >= max_dist_in_heap) then
                        should_explore_left = .false.
                        ! print *, "Pruned left subtree"
                    end if

                    ! Check the right subtree
                    if (axis_dist < 0.0_real64 .and. axis_dist * axis_dist >= max_dist_in_heap) then
                        should_explore_right = .false.
                        ! print *, "Pruned right subtree"
                    end if
                end if

                ! Add child ranges to the stack
                ! Push right subtree first

                if (should_explore_right .and. mid_idx < right_idx) then
                    stack_pos = stack_pos + 1
                    range_stack(1, stack_pos) = mid_idx + 1      ! left_idx of right subtree
                    range_stack(2, stack_pos) = right_idx        ! right_idx of right subtree
                    range_stack(3, stack_pos) = current_depth + 1
                    ! print *, "Pushed right subtree to stack: left_idx=", mid_idx + 1, "right_idx=", right_idx
                end if

                ! Push left subtree
                if (should_explore_left .and. left_idx < mid_idx) then
                    stack_pos = stack_pos + 1
                    range_stack(1, stack_pos) = left_idx         ! left_idx of left subtree  
                    range_stack(2, stack_pos) = mid_idx - 1      ! right_idx of left subtree
                    range_stack(3, stack_pos) = current_depth + 1
                    ! print *, "Pushed left subtree to stack: left_idx=", left_idx, "right_idx=", mid_idx - 1
                end if
            end if
        end do

        ! Convert squared distances to real distances
        do i = 1, min(current_best_count, k_neighbors)
            distances(i) = sqrt(distances(i))
        end do

        ! Fill remaining slots if fewer than k neighbors were found
        do i = current_best_count + 1, k_neighbors
            neighbors(i) = 0    ! Invalid index indicates no neighbor found
            distances(i) = huge(1.0_real64)
        end do
        ! print *, "Finished kd_knn_query with current_best_count:", current_best_count

    end subroutine kd_knn_query


    !> Maintain max-heap property by bubbling element up from given position
    !! Used when inserting new elements into the heap
    pure subroutine max_heapify_up(distances, neighbors, pos)
        real(real64), intent(inout) :: distances(:)     !! Distance array (heap)
        integer(int32), intent(inout) :: neighbors(:)   !! Neighbor indices (heap)
        integer(int32), intent(in) :: pos               !! Position to bubble up from
        
        integer(int32) :: parent_pos, current_pos
        real(real64) :: temp_dist
        integer(int32) :: temp_neighbor
        
        current_pos = pos
        
        ! Bubble up while current element is larger than its parent
        do while (current_pos > 1)
            parent_pos = current_pos / 2
            
            ! If heap property is satisfied, stop
            if (distances(current_pos) <= distances(parent_pos)) exit
            
            ! Swap current with parent
            temp_dist = distances(current_pos)
            temp_neighbor = neighbors(current_pos)
            distances(current_pos) = distances(parent_pos)
            neighbors(current_pos) = neighbors(parent_pos)
            distances(parent_pos) = temp_dist
            neighbors(parent_pos) = temp_neighbor
            
            current_pos = parent_pos
        end do
    end subroutine max_heapify_up

    !> Maintain max-heap property by bubbling element down from root
    !! Used when replacing the maximum element (root) of the heap
    !> Maintain max-heap property by bubbling element down from root
    !! Used when replacing the maximum element (root) of the heap
    pure subroutine max_heapify_down(distances, neighbors, heap_size)
        real(real64), intent(inout) :: distances(:)     !! Distance array (heap)
        integer(int32), intent(inout) :: neighbors(:)   !! Neighbor indices (heap)
        integer(int32), intent(in) :: heap_size         !! Size of the heap
        
        integer(int32) :: current_pos, left_child, right_child, largest_pos
        real(real64) :: temp_dist
        integer(int32) :: temp_neighbor
        
        current_pos = 1  ! Start from root
        
        ! Bubble down while heap property is violated
        do
            left_child = 2 * current_pos
            right_child = 2 * current_pos + 1
            largest_pos = current_pos
            
            ! Find the position with largest distance among current, left child, right child
            if (left_child <= heap_size) then
                if(distances(left_child) > distances(largest_pos)) then
                    largest_pos = left_child
                end if
            end if
            
            if (right_child <= heap_size) then
                if(distances(right_child) > distances(largest_pos)) then
                    largest_pos = right_child
                end if
            end if
            
            ! If heap property is satisfied, stop
            if (largest_pos == current_pos) exit
            
            ! Swap current with the largest child
            temp_dist = distances(current_pos)
            temp_neighbor = neighbors(current_pos)
            distances(current_pos) = distances(largest_pos)
            neighbors(current_pos) = neighbors(largest_pos)
            distances(largest_pos) = temp_dist
            neighbors(largest_pos) = temp_neighbor
            
            current_pos = largest_pos
        end do
    end subroutine max_heapify_down

end module kd_tree

!> R interface for building KD index
pure subroutine build_kd_index_r(points, num_dimensions, num_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr)
    use kd_tree, only: build_kd_index
    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32), intent(in) :: num_dimensions      
    !! Number of dimensions
    integer(int32), intent(in) :: num_points          
    !! Number of points
    real(real64), intent(in) :: points(num_dimensions, num_points)  
    !! Input points
    integer(int32), intent(in) :: dimension_order(num_dimensions)   
    !! Dimension order
    integer(int32), intent(out) :: kd_indices(num_points)           
    !! Output indices
    integer(int32), intent(out) :: workspace(num_points)          
    !! Workspace
    real(real64), intent(out) :: value_buffer(num_points)         
    !! Value buffer
    integer(int32), intent(out) :: permutation(num_points)        
    !! Permutation array
    integer(int32), intent(out) :: left_stack(num_points)         
    !! Left stack
    integer(int32), intent(out) :: right_stack(num_points)        
    !! Right stack
    integer(int32), intent(out) :: ierr                
    !! Error code

    integer(int32) :: recursion_stack(3, num_points)

    call build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                      workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
end subroutine build_kd_index_r

!> R interface for building spherical KD index
pure subroutine build_spherical_kd_r(vectors, num_dimensions, num_vectors, sphere_indices, &
                              dimension_order, workspace, value_buffer, permutation, &
                              left_stack, right_stack, ierr)
    use kd_tree, only: build_spherical_kd
    use iso_fortran_env, only: int32, real64
    implicit none
    integer(int32), intent(in) :: num_dimensions      
    !! Number of dimensions
    integer(int32), intent(in) :: num_vectors         
    !! Number of vectors
    real(real64), intent(in) :: vectors(num_dimensions, num_vectors)  
    !! Input vectors
    integer(int32), intent(out) :: sphere_indices(num_vectors)  
    !! Output indices
    integer(int32), intent(out) :: dimension_order(num_dimensions)  
    !! Dimension order
    integer(int32), intent(out) :: workspace(num_vectors)     
    !! Workspace
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

    integer(int32) :: recursion_stack(3, num_vectors)

    call build_spherical_kd(vectors, num_dimensions, num_vectors, sphere_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
end subroutine build_spherical_kd_r

!> C interface for building KD index
pure subroutine build_kd_index_C(points, num_dimensions, num_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr) &
                          bind(C, name="build_kd_index_C")
    use iso_c_binding, only: c_int, c_double, c_f_pointer, c_loc
    use iso_fortran_env, only : int32
    use kd_tree, only: build_kd_index
    implicit none
    integer(c_int), value :: num_dimensions
    integer(c_int), value :: num_points
    real(c_double), intent(in) :: points(num_dimensions,num_points)
    integer(c_int), intent(in) :: dimension_order(num_dimensions)
    integer(c_int), intent(out) :: kd_indices(num_points)
    integer(c_int), intent(out) :: workspace(num_points)
    real(c_double), intent(out) :: value_buffer(num_points)
    integer(c_int), intent(out) :: permutation(num_points)
    integer(c_int), intent(out) :: left_stack(num_points)
    integer(c_int), intent(out) :: right_stack(num_points)
    integer(c_int), intent(out) :: ierr

    integer(int32) :: recursion_stack(3, num_points)

    ! Call the original implementation
    call build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                      workspace, value_buffer, permutation, left_stack, right_stack, recursion_stack, ierr)
end subroutine build_kd_index_C