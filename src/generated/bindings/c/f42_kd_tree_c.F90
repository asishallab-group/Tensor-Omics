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
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: build_kd_index_c
    public :: build_kd_index_expert_c
    public :: build_spherical_kd_c
    public :: build_spherical_kd_expert_c
    public :: vicinity_vectors_c
    public :: vicinity_vectors_expert_c
    public :: vicinity_vectors_count_c
    public :: vicinity_vectors_count_expert_c

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

    !> summary: C-wrapper for [[f42_kd_tree(module):vicinity_vectors(subroutine)]]
    subroutine vicinity_vectors_c(&
            query_point,&
            points,&
            n_dimensions,&
            n_points,&
            r,&
            dimension_order,&
            kd_indices,&
            vicinity_mask,&
            ierr&
        ) bind(C, name="vicinity_vectors_c")
        use f42_kd_tree, only: vicinity_vectors

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Total number of points organized in the k-d tree
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Coordinate vector used as the center of the search
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Ambient point matrix
        real(c_double), intent(in), target :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        logical(c_bool), dimension(n_points), intent(out), target :: vicinity_mask
            !! Mask indicating points within the search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(r)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(vicinity_mask, n_points)

        call vicinity_vectors(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            vicinity_mask = vicinity_mask,&
            ierr = ierr&
        )
    end subroutine vicinity_vectors_c

    !> summary: C-wrapper for [[f42_kd_tree(module):vicinity_vectors_expert(subroutine)]]
    subroutine vicinity_vectors_expert_c(&
            query_point,&
            points,&
            n_dimensions,&
            n_points,&
            r,&
            dimension_order,&
            kd_indices,&
            tmp_stack,&
            vicinity_mask,&
            ierr&
        ) bind(C, name="vicinity_vectors_expert_c")
        use f42_kd_tree, only: vicinity_vectors_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Total number of points organized in the k-d tree
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Coordinate vector used as the center of the search
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Ambient point matrix
        real(c_double), intent(in), target :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(3, 64), intent(inout), target :: tmp_stack
            !! Preallocated k-d tree traversal stack
        logical(c_bool), dimension(n_points), intent(out), target :: vicinity_mask
            !! Mask indicating points within the search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(r)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64)
        M_CHECK_ARRAY_NON_NULL(vicinity_mask, n_points)

        call vicinity_vectors_expert(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            tmp_stack = tmp_stack,&
            vicinity_mask = vicinity_mask,&
            ierr = ierr&
        )
    end subroutine vicinity_vectors_expert_c

    !> summary: C-wrapper for [[f42_kd_tree(module):vicinity_vectors_count(subroutine)]]
    subroutine vicinity_vectors_count_c(&
            query_point,&
            points,&
            n_dimensions,&
            n_points,&
            r,&
            dimension_order,&
            kd_indices,&
            n_neighbors,&
            ierr&
        ) bind(C, name="vicinity_vectors_count_c")
        use f42_kd_tree, only: vicinity_vectors_count

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Total number of points organized in the k-d tree
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Coordinate vector used as the center of the search
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Ambient point matrix
        real(c_double), intent(in), target :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), intent(out), target :: n_neighbors
            !! Number of points within the search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)

        call vicinity_vectors_count(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            n_neighbors = n_neighbors,&
            ierr = ierr&
        )
    end subroutine vicinity_vectors_count_c

    !> summary: C-wrapper for [[f42_kd_tree(module):vicinity_vectors_count_expert(subroutine)]]
    subroutine vicinity_vectors_count_expert_c(&
            query_point,&
            points,&
            n_dimensions,&
            n_points,&
            r,&
            dimension_order,&
            kd_indices,&
            tmp_stack,&
            n_neighbors,&
            ierr&
        ) bind(C, name="vicinity_vectors_count_expert_c")
        use f42_kd_tree, only: vicinity_vectors_count_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Total number of points organized in the k-d tree
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Coordinate vector used as the center of the search
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Ambient point matrix
        real(c_double), intent(in), target :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(3, 64), intent(inout), target :: tmp_stack
            !! Preallocated k-d tree traversal stack
        integer(c_int), intent(out), target :: n_neighbors
            !! Number of points within the search radius
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(r)
        M_CHECK_NON_NULL(n_neighbors)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(tmp_stack, 3 * 64)

        call vicinity_vectors_count_expert(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            tmp_stack = tmp_stack,&
            n_neighbors = n_neighbors,&
            ierr = ierr&
        )
    end subroutine vicinity_vectors_count_expert_c

end module f42_kd_tree_c
#endif
