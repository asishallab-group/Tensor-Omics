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

end module f42_kd_tree



!> C interface for building KD index
pure subroutine build_kd_index_C(points, num_dimensions, num_points, kd_indices, dimension_order, &
                          workspace, value_buffer, permutation, left_stack, right_stack, ierr) &
                          bind(C, name="build_kd_index_C")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only : int32
    use f42_kd_tree, only: build_kd_index
    M_USE_NULL_VALIDATION
    implicit none
    !| Input parameters
    integer(c_int), intent(in), target :: num_dimensions
    !| number of input points
    integer(c_int), intent(in), target :: num_points
    !| input points
    real(c_double), intent(in), target :: points(num_dimensions,num_points)
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
    use iso_fortran_env, only : int32
    use f42_kd_tree, only: build_spherical_kd
    M_USE_NULL_VALIDATION
    implicit none
    !| number of input dimensions
    integer(c_int), intent(in), target :: num_dimensions
    !| number of input vectors
    integer(c_int), intent(in), target :: num_vectors
    !| input unit vectors
    real(c_double), intent(in), target :: vectors(num_dimensions,num_vectors)
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