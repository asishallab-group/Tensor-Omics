#include "src/macros.h"

module f42_kd_tree
    use safeguard
    use f42_utils, only: sort_array_heapsort, init_perm
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, validate_dimension_size, validate_all_in_range_int, validate_in_range_int, is_err, set_err, ERR_ALLOC_FAIL
    implicit none

contains

    !> AUTHOR_AARON_SCHROEDER
    !| Build a k-d tree index using a stack-based, non-recursive approach.
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

    !> AUTHOR_AARON_SCHROEDER
    !| Build spherical k-d tree index
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
    !| Build spherical k-d tree index
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
    !| Get point from KD index
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

end module f42_kd_tree

!> C interface for building KD index
pure subroutine build_kd_index_c(points, n_dimensions, n_points, kd_indices, dimension_order, ierr) bind(C, name="build_kd_index_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32
    use f42_kd_tree, only: build_kd_index_alloc
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_dimensions
        !! Number of dimensions
    integer(c_int), intent(in), target :: n_points
        !! Number of points
    real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
        !! Data points
    integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
        !! Dimension order (by variance)
    integer(c_int), dimension(n_points), intent(out), target :: kd_indices
        !! Output index array (k-d tree order)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dimensions)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(points)
    M_CHECK_NON_NULL(dimension_order)
    M_CHECK_NON_NULL(kd_indices)

    ! Call the original implementation
    call build_kd_index_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, ierr)
end subroutine build_kd_index_C

!> C interface for building spherical KD index
pure subroutine build_spherical_kd_c(points, n_dimensions, n_points, kd_indices, dimension_order, ierr) bind(C, name="build_spherical_kd_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use, intrinsic :: iso_fortran_env, only: int32
    use f42_kd_tree, only: build_spherical_kd_alloc
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_dimensions
        !! Number of dimensions
    integer(c_int), intent(in), target :: n_points
        !! Number of points
    real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
        !! Data points
    integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
        !! Dimension order (by variance)
    integer(c_int), dimension(n_points), intent(out), target :: kd_indices
        !! Output index array (k-d tree order)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dimensions)
    M_CHECK_NON_NULL(n_points)
    M_CHECK_NON_NULL(points)
    M_CHECK_NON_NULL(dimension_order)
    M_CHECK_NON_NULL(kd_indices)

    ! Call the original implementation
    call build_spherical_kd_alloc(points, n_dimensions, n_points, kd_indices, dimension_order, ierr)
end subroutine build_spherical_kd_C
