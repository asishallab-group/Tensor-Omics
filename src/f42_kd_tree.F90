#include "macros.h"

module f42_kd_tree
    use safeguard
    use f42_utils, only: sort_array
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: ERR_INVALID_INPUT, ERR_DIM_MISMATCH, set_ok, set_err_once, is_ok, validate_dimension_size
    implicit none
    private
    public :: build_kd_index, build_spherical_kd, get_kd_point

contains

    !> category: C-interface
    !| Build a k-d tree index using a stack-based, non-recursive approach.
    pure subroutine build_kd_index(points, num_dimensions, num_points, kd_indices, dimension_order, &
                            tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_left_stack, tmp_right_stack, tmp_recursion_stack, ierr)
        integer(int32), intent(in) :: num_dimensions      
        !! Number of dimensions
        integer(int32), intent(in) :: num_points          
        !! Number of points
        real(real64), intent(in) :: points(num_dimensions, num_points)  
        !! Data points
        integer(int32), intent(in) :: dimension_order(num_dimensions)   
        !! Dimension order (by variance)
        integer(int32), intent(out) :: tmp_recursion_stack(3, num_points) 
        !! Stack for l, r, depth

        integer(int32), intent(out) :: kd_indices(num_points)           
        !! Output index array (k-d tree order)
        integer(int32), intent(out) :: tmp_workspace(num_points)          
        !! Workspace array
        real(real64), intent(out) :: tmp_value_buffer(num_points)         
        !! Workspace for sorting
        integer(int32), intent(out) :: tmp_permutation(num_points)        
        !! Workspace for sorting
        integer(int32), intent(out) :: tmp_left_stack(num_points)                  
        !! Workspace for sorting
        integer(int32), intent(out) :: tmp_right_stack(num_points)                 
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
        tmp_recursion_stack(1, 1) = 1
        tmp_recursion_stack(2, 1) = num_points
        tmp_recursion_stack(3, 1) = 0

        do while (stack_top > 0)
            left_idx = tmp_recursion_stack(1, stack_top)
            right_idx = tmp_recursion_stack(2, stack_top)
            current_depth = tmp_recursion_stack(3, stack_top)
            stack_top = stack_top - 1

            if (right_idx <= left_idx) cycle

            !! Choose split dimension by cycling through dimension_order
            current_dim = dimension_order(mod(current_depth, num_dimensions) + 1)

            !! Find median index
            mid_idx = left_idx + (right_idx - left_idx) / 2

            !! Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
            call partial_sort_by_dimension(points, num_points, num_dimensions, kd_indices, left_idx, right_idx, &
                                        current_dim, mid_idx, tmp_workspace, tmp_value_buffer, tmp_permutation, &
                                        tmp_left_stack, tmp_right_stack, ierr)
            if (.not. is_ok(ierr)) return

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
    end subroutine build_kd_index

    !> Helper: sorts kd_indices(left_idx:right_idx) by points(dimension, kd_indices(:))
    pure subroutine partial_sort_by_dimension(points, n_points, num_dimensions, kd_indices, left_idx, right_idx, &
                                        dim, mid_idx, tmp_workspace, tmp_value_buffer, tmp_permutation, &
                                        tmp_left_stack, tmp_right_stack, ierr)
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
        integer(int32), intent(out) :: tmp_workspace(:)          
        !! Workspace array
        real(real64), intent(out) :: tmp_value_buffer(:)         
        !! Buffer for dimension values
        integer(int32), intent(out) :: tmp_permutation(:)        
        !! Permutation array
        integer(int32), intent(out) :: tmp_left_stack(:)         
        !! Stack for left indices
        integer(int32), intent(out) :: tmp_right_stack(:)        
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
        
        !! Fill tmp_value_buffer with the values of points(dimension, kd_indices(left_idx:right_idx))
        do i = 1, subarray_size
            tmp_value_buffer(i) = points(dim, kd_indices(left_idx + i - 1))
            tmp_permutation(i) = i
        end do

        call sort_array(tmp_value_buffer(1:subarray_size), tmp_permutation(1:subarray_size), tmp_left_stack, tmp_right_stack)

        !! Reorder kd_indices(left_idx:right_idx) according to tmp_permutation
        do i = 1, subarray_size
            if (tmp_permutation(i) < 1 .or. tmp_permutation(i) > subarray_size) then
                ierr = ERR_INVALID_INPUT
                return
            end if
            tmp_workspace(i) = kd_indices(left_idx + tmp_permutation(i) - 1)
        end do
        do i = 1, subarray_size
            kd_indices(left_idx + i - 1) = tmp_workspace(i)
        end do
    end subroutine partial_sort_by_dimension

    !> category: C-interface
    !| Build spherical k-d tree index
    pure subroutine build_spherical_kd(vectors, num_dimensions, num_vectors, sphere_indices, &
                                dimension_order, tmp_workspace, tmp_value_buffer, tmp_permutation, &
                                tmp_left_stack, tmp_right_stack, tmp_recursion_stack, ierr)

        integer(int32), intent(in) :: num_dimensions      
        !! Number of dimensions
        integer(int32), intent(in) :: num_vectors         
        !! Number of vectors
        real(real64), intent(in) :: vectors(num_dimensions, num_vectors)  
        !! Input unit vectors
        integer(int32), intent(out) :: tmp_recursion_stack(3, num_vectors)
        !! Stack for recursive calls
        integer(int32), intent(out) :: sphere_indices(num_vectors)  
        !! Output index array
        integer(int32), intent(in) :: dimension_order(num_dimensions)  
        !! Dimension order
        integer(int32), intent(out) :: tmp_workspace(num_vectors)     
        !! Workspace array
        real(real64), intent(out) :: tmp_value_buffer(num_vectors)    
        !! Value buffer
        integer(int32), intent(out) :: tmp_permutation(num_vectors)   
        !! Permutation array
        integer(int32), intent(out) :: tmp_left_stack(num_vectors)              
        !! Left stack
        integer(int32), intent(out) :: tmp_right_stack(num_vectors)             
        !! Right stack
        integer(int32), intent(out) :: ierr                         
        !! Error code

        call build_kd_index(vectors, num_dimensions, num_vectors, sphere_indices, dimension_order, &
                          tmp_workspace, tmp_value_buffer, tmp_permutation, tmp_left_stack, tmp_right_stack, tmp_recursion_stack, ierr)
    end subroutine build_spherical_kd

    !> category: C-interface
    !| Get point from KD index
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

end module f42_kd_tree
