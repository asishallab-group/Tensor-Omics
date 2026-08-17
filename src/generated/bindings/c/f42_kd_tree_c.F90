#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_kd_tree(module)]]
!| k-d tree spatial index over fixed-dimensional point sets.
!| Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
!| caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
!| is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
!| index array rather than as linked nodes.
module f42_kd_tree_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: build_kd_index_c
    public :: build_kd_index_expert_c
    public :: build_spherical_kd_c
    public :: build_spherical_kd_expert_c

contains

    !> summary: C-wrapper for [[f42_kd_tree(module):build_kd_index(subroutine)]]
    subroutine build_kd_index_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            ierr&
        ) bind(C, name="build_kd_index_c")
        use f42_kd_tree, only: build_kd_index

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Data points
        integer(c_int), dimension(n_points), intent(out), target :: kd_indices
            !! Output index array (k-d tree order)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)

        call build_kd_index(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            ierr = ierr&
        )
    end subroutine build_kd_index_c

    !> summary: C-wrapper for [[f42_kd_tree(module):build_kd_index_expert(subroutine)]]
    subroutine build_kd_index_expert_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            tmp_workspace,&
            tmp_value_buffer,&
            tmp_permutation,&
            tmp_recursion_stack,&
            ierr&
        ) bind(C, name="build_kd_index_expert_c")
        use f42_kd_tree, only: build_kd_index_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Data points
        integer(c_int), dimension(n_points), intent(out), target :: kd_indices
            !! Output index array (k-d tree order)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_points), intent(out), target :: tmp_workspace
            !! Workspace array
        real(c_double), dimension(n_points), intent(out), target :: tmp_value_buffer
            !! Workspace for sorting
        integer(c_int), dimension(n_points), intent(out), target :: tmp_permutation
            !! Workspace for sorting
        integer(c_int), dimension(3, n_points), intent(out), target :: tmp_recursion_stack
            !! Stack for l, r, depth
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_workspace, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_value_buffer, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_recursion_stack, 3 * n_points)

        call build_kd_index_expert(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_recursion_stack = tmp_recursion_stack,&
            ierr = ierr&
        )
    end subroutine build_kd_index_expert_c

    !> summary: C-wrapper for [[f42_kd_tree(module):build_spherical_kd(subroutine)]]
    !| This is a thin, semantically-named wrapper: partitioning is identical to
    !| [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]] (plain per-axis median splits);
    !| callers are responsible for ensuring `points` are actually unit-normalized beforehand.
    subroutine build_spherical_kd_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            ierr&
        ) bind(C, name="build_spherical_kd_c")
        use f42_kd_tree, only: build_spherical_kd

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Data points
        integer(c_int), dimension(n_points), intent(out), target :: kd_indices
            !! Output index array (k-d tree order)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)

        call build_spherical_kd(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            ierr = ierr&
        )
    end subroutine build_spherical_kd_c

    !> summary: C-wrapper for [[f42_kd_tree(module):build_spherical_kd_expert(subroutine)]]
    !| This is a thin, semantically-named wrapper: partitioning is identical to
    !| [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]] (plain per-axis median splits);
    !| callers are responsible for ensuring `points` are actually unit-normalized beforehand.
    subroutine build_spherical_kd_expert_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            tmp_workspace,&
            tmp_value_buffer,&
            tmp_permutation,&
            tmp_recursion_stack,&
            ierr&
        ) bind(C, name="build_spherical_kd_expert_c")
        use f42_kd_tree, only: build_spherical_kd_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Data points
        integer(c_int), dimension(n_points), intent(out), target :: kd_indices
            !! Output index array (k-d tree order)
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_points), intent(out), target :: tmp_workspace
            !! Workspace array
        real(c_double), dimension(n_points), intent(out), target :: tmp_value_buffer
            !! Workspace for sorting
        integer(c_int), dimension(n_points), intent(out), target :: tmp_permutation
            !! Workspace for sorting
        integer(c_int), dimension(3, n_points), intent(out), target :: tmp_recursion_stack
            !! Stack for l, r, depth
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_workspace, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_value_buffer, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_permutation, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_recursion_stack, 3 * n_points)

        call build_spherical_kd_expert(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_recursion_stack = tmp_recursion_stack,&
            ierr = ierr&
        )
    end subroutine build_spherical_kd_expert_c

end module f42_kd_tree_c
#endif
