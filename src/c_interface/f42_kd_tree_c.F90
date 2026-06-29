#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[f42_kd_tree(module)]]
module f42_kd_tree_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_int, c_double, c_char, c_double_complex
    use, intrinsic :: iso_c_binding, only: c_loc, c_associated

    use tox_conversions, only: logical_as_c_int, c_int_as_logical
    use tox_conversions, only: c_char_as_char, char_as_c_char
    use tox_conversions, only: string_as_c_char_1d, c_char_1d_as_string
    use tox_conversions, only: string_as_c_char_2d, c_char_2d_as_string

    use tox_errors, only: ERR_POINTER_NULL, is_err, set_err, ERR_ALLOC_FAIL, ERR_INVALID_INPUT
    implicit none
contains

    !> summary: C-wrapper for [[f42_kd_tree(module):build_kd_index(subroutine)]]
    !| Build a k-d tree index using a stack-based, non-recursive approach.
    !| Initialize kd_indices to 1:num_points (original indices)
    !| Choose split dimension by cycling through dimension_order
    !| Find median index
    !| Partition kd_indices(left_idx:right_idx) by points(current_dim, kd_indices(:))
    !| Push right and left intervals onto stack
    subroutine build_kd_index_c(&
            points,&
            num_dimensions,&
            num_points,&
            kd_indices,&
            dimension_order,&
            tmp_workspace,&
            tmp_value_buffer,&
            tmp_permutation,&
            tmp_left_stack,&
            tmp_right_stack,&
            tmp_recursion_stack,&
            ierr&
            ) bind(C, name="build_kd_index_c")
        use f42_kd_tree, only: build_kd_index
        use f42_kd_tree
        integer(c_int), intent(in), target :: num_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: num_points
            !! Number of points
        real(c_double), intent(in), dimension(num_dimensions, num_points), target :: points
            !! Data points
        integer(c_int), intent(out), dimension(num_points), target :: kd_indices
            !! Output index array (k-d tree order)
        integer(c_int), intent(in), dimension(num_dimensions), target :: dimension_order
            !! Dimension order (by variance)
        integer(c_int), intent(out), dimension(num_points), target :: tmp_workspace
            !! Workspace array
        real(c_double), intent(out), dimension(num_points), target :: tmp_value_buffer
            !! Workspace for sorting
        integer(c_int), intent(out), dimension(num_points), target :: tmp_permutation
            !! Workspace for sorting
        integer(c_int), intent(out), dimension(num_points), target :: tmp_left_stack
            !! Workspace for sorting
        integer(c_int), intent(out), dimension(num_points), target :: tmp_right_stack
            !! Workspace for sorting
        integer(c_int), intent(out), dimension(3, num_points), target :: tmp_recursion_stack
            !! Stack for l, r, depth
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(points)
        M_CHECK_NON_NULL(num_dimensions)
        M_CHECK_NON_NULL(num_points)
        M_CHECK_NON_NULL(kd_indices)
        M_CHECK_NON_NULL(dimension_order)
        M_CHECK_NON_NULL(tmp_workspace)
        M_CHECK_NON_NULL(tmp_value_buffer)
        M_CHECK_NON_NULL(tmp_permutation)
        M_CHECK_NON_NULL(tmp_left_stack)
        M_CHECK_NON_NULL(tmp_right_stack)
        M_CHECK_NON_NULL(tmp_recursion_stack)
        call build_kd_index(&
            points = points,&
            num_dimensions = num_dimensions,&
            num_points = num_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_left_stack = tmp_left_stack,&
            tmp_right_stack = tmp_right_stack,&
            tmp_recursion_stack = tmp_recursion_stack,&
            ierr = ierr&
        )
    end subroutine build_kd_index_c

    !> summary: C-wrapper for [[f42_kd_tree(module):build_spherical_kd(subroutine)]]
    !| Build spherical k-d tree index
    subroutine build_spherical_kd_c(&
            vectors,&
            num_dimensions,&
            num_vectors,&
            sphere_indices,&
            dimension_order,&
            tmp_workspace,&
            tmp_value_buffer,&
            tmp_permutation,&
            tmp_left_stack,&
            tmp_right_stack,&
            tmp_recursion_stack,&
            ierr&
            ) bind(C, name="build_spherical_kd_c")
        use f42_kd_tree, only: build_spherical_kd
        use f42_kd_tree
        integer(c_int), intent(in), target :: num_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: num_vectors
            !! Number of vectors
        real(c_double), intent(in), dimension(num_dimensions, num_vectors), target :: vectors
            !! Input unit vectors
        integer(c_int), intent(out), dimension(num_vectors), target :: sphere_indices
            !! Output index array
        integer(c_int), intent(in), dimension(num_dimensions), target :: dimension_order
            !! Dimension order
        integer(c_int), intent(out), dimension(num_vectors), target :: tmp_workspace
            !! Workspace array
        real(c_double), intent(out), dimension(num_vectors), target :: tmp_value_buffer
            !! Value buffer
        integer(c_int), intent(out), dimension(num_vectors), target :: tmp_permutation
            !! Permutation array
        integer(c_int), intent(out), dimension(num_vectors), target :: tmp_left_stack
            !! Left stack
        integer(c_int), intent(out), dimension(num_vectors), target :: tmp_right_stack
            !! Right stack
        integer(c_int), intent(out), dimension(3, num_vectors), target :: tmp_recursion_stack
            !! Stack for recursive calls
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vectors)
        M_CHECK_NON_NULL(num_dimensions)
        M_CHECK_NON_NULL(num_vectors)
        M_CHECK_NON_NULL(sphere_indices)
        M_CHECK_NON_NULL(dimension_order)
        M_CHECK_NON_NULL(tmp_workspace)
        M_CHECK_NON_NULL(tmp_value_buffer)
        M_CHECK_NON_NULL(tmp_permutation)
        M_CHECK_NON_NULL(tmp_left_stack)
        M_CHECK_NON_NULL(tmp_right_stack)
        M_CHECK_NON_NULL(tmp_recursion_stack)
        call build_spherical_kd(&
            vectors = vectors,&
            num_dimensions = num_dimensions,&
            num_vectors = num_vectors,&
            sphere_indices = sphere_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_left_stack = tmp_left_stack,&
            tmp_right_stack = tmp_right_stack,&
            tmp_recursion_stack = tmp_recursion_stack,&
            ierr = ierr&
        )
    end subroutine build_spherical_kd_c

    !> summary: C-wrapper for [[f42_kd_tree(module):get_kd_point(subroutine)]]
    !| Get point from KD index
    subroutine get_kd_point_c(&
            points,&
            n_points_elements_dim_1,&
            n_points_elements_dim_2,&
            kd_indices,&
            n_kd_indices_elements,&
            position,&
            point_values,&
            n_point_values_elements,&
            ierr&
            ) bind(C, name="get_kd_point_c")
        use f42_kd_tree, only: get_kd_point
        use f42_kd_tree
        integer(c_int), intent(in), target :: n_points_elements_dim_1
            !! Size of the 1. dimension/extent of `points`
        integer(c_int), intent(in), target :: n_points_elements_dim_2
            !! Size of the 2. dimension/extent of `points`
        integer(c_int), intent(in), target :: n_kd_indices_elements
            !! Size of the 1. dimension/extent of `kd_indices`
        integer(c_int), intent(in), target :: n_point_values_elements
            !! Size of the 1. dimension/extent of `point_values`
        real(c_double), intent(in), dimension(n_points_elements_dim_1, n_points_elements_dim_2), target :: points
            !! Input points
        integer(c_int), intent(in), dimension(n_kd_indices_elements), target :: kd_indices
            !! KD index array
        integer(c_int), intent(in), target :: position
            !! Position in index
        real(c_double), intent(out), dimension(n_point_values_elements), target :: point_values
            !! Output point values
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(points)
        M_CHECK_NON_NULL(n_points_elements_dim_1)
        M_CHECK_NON_NULL(n_points_elements_dim_2)
        M_CHECK_NON_NULL(kd_indices)
        M_CHECK_NON_NULL(n_kd_indices_elements)
        M_CHECK_NON_NULL(position)
        M_CHECK_NON_NULL(point_values)
        M_CHECK_NON_NULL(n_point_values_elements)
        call get_kd_point(&
            points = points,&
            kd_indices = kd_indices,&
            position = position,&
            point_values = point_values,&
            ierr = ierr&
        )
    end subroutine get_kd_point_c

end module f42_kd_tree_c
#endif