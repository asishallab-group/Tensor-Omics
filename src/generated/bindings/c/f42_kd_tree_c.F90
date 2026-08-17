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
    public :: kd_knn_query_c
    public :: kd_knn_query_expert_c
    public :: kd_range_query_mask_c
    public :: kd_range_query_mask_expert_c
    public :: kd_range_query_list_c
    public :: kd_range_query_list_expert_c
    public :: kd_range_query_count_c
    public :: kd_range_query_count_expert_c

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

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_knn_query(subroutine)]]
    !| Via a bounded max-heap kept directly in `neighbors`/`distances` and splitting-plane
    !| pruning. Requires `k_neighbors <= n_points`: every point is then guaranteed visited
    !| before the heap can still have room, so no fallback for "fewer than k found" is needed.
    !| Does not guarantee `neighbors`/`distances` are sorted nearest-to-farthest (max-heap
    !| order internally).
    subroutine kd_knn_query_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            k_neighbors,&
            neighbors,&
            distances,&
            ierr&
        ) bind(C, name="kd_knn_query_c")
        use f42_kd_tree, only: kd_knn_query

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        integer(c_int), intent(in), target :: k_neighbors
            !! Number of neighbors to find
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index, see [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]]
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        integer(c_int), dimension(k_neighbors), intent(out), target :: neighbors
            !! Output: indices of the k nearest neighbors (nearest-to-farthest order not guaranteed, max-heap order internally)
        real(c_double), dimension(k_neighbors), intent(out), target :: distances
            !! Output: Euclidean distances to the k nearest neighbors
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(k_neighbors)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(neighbors, k_neighbors)
        M_CHECK_ARRAY_NON_NULL(distances, k_neighbors)

        call kd_knn_query(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            k_neighbors = k_neighbors,&
            neighbors = neighbors,&
            distances = distances,&
            ierr = ierr&
        )
    end subroutine kd_knn_query_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_knn_query_expert(subroutine)]]
    !| Via a bounded max-heap kept directly in `neighbors`/`distances` and splitting-plane
    !| pruning. Requires `k_neighbors <= n_points`: every point is then guaranteed visited
    !| before the heap can still have room, so no fallback for "fewer than k found" is needed.
    !| Does not guarantee `neighbors`/`distances` are sorted nearest-to-farthest (max-heap
    !| order internally).
    subroutine kd_knn_query_expert_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            k_neighbors,&
            tmp_range_stack,&
            neighbors,&
            distances,&
            ierr&
        ) bind(C, name="kd_knn_query_expert_c")
        use f42_kd_tree, only: kd_knn_query_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        integer(c_int), intent(in), target :: k_neighbors
            !! Number of neighbors to find
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index, see [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]]
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        integer(c_int), dimension(3, n_points), intent(out), target :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(c_int), dimension(k_neighbors), intent(out), target :: neighbors
            !! Output: indices of the k nearest neighbors (nearest-to-farthest order not guaranteed, max-heap order internally)
        real(c_double), dimension(k_neighbors), intent(out), target :: distances
            !! Output: Euclidean distances to the k nearest neighbors
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(k_neighbors)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbors, k_neighbors)
        M_CHECK_ARRAY_NON_NULL(distances, k_neighbors)

        call kd_knn_query_expert(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            k_neighbors = k_neighbors,&
            tmp_range_stack = tmp_range_stack,&
            neighbors = neighbors,&
            distances = distances,&
            ierr = ierr&
        )
    end subroutine kd_knn_query_expert_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_range_query_mask(subroutine)]]
    !| Same iterative, stack-based traversal and splitting-plane pruning as
    !| [[f42_kd_tree_impl(module):kd_knn_query_impl(subroutine)]] (the near side of each split is
    !| always explored, the far side only when it is still within `radius` of the splitting
    !| plane), with a fixed radius bound instead of a k-nearest-neighbor heap. Compares
    !| squared distances against a precomputed `radius**2` (no `sqrt` per node visited).
    !|
    !| Fits a caller that already does an O(n_points) pass over the result (e.g. merging it
    !| into an existing coverage mask via `.or.`). A caller issuing many independent range
    !| queries per outer step (e.g. one per candidate point in a greedy loop) should use
    !| [[f42_kd_tree_impl(module):kd_range_query_list_impl(subroutine)]] or
    !| [[f42_kd_tree_impl(module):kd_range_query_count_impl(subroutine)]] instead, since the
    !| `in_radius_mask = .false.` reset here costs O(n_points) on every call regardless of how
    !| few points are actually found.
    subroutine kd_range_query_mask_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            in_radius_mask,&
            ierr&
        ) bind(C, name="kd_range_query_mask_c")
        use f42_kd_tree, only: kd_range_query_mask

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        real(c_double), intent(in), target :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), dimension(n_points), intent(out), target :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(in_radius_mask, n_points)

        call kd_range_query_mask(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            in_radius_mask = in_radius_mask,&
            ierr = ierr&
        )
    end subroutine kd_range_query_mask_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_range_query_mask_expert(subroutine)]]
    !| Same iterative, stack-based traversal and splitting-plane pruning as
    !| [[f42_kd_tree_impl(module):kd_knn_query_impl(subroutine)]] (the near side of each split is
    !| always explored, the far side only when it is still within `radius` of the splitting
    !| plane), with a fixed radius bound instead of a k-nearest-neighbor heap. Compares
    !| squared distances against a precomputed `radius**2` (no `sqrt` per node visited).
    !|
    !| Fits a caller that already does an O(n_points) pass over the result (e.g. merging it
    !| into an existing coverage mask via `.or.`). A caller issuing many independent range
    !| queries per outer step (e.g. one per candidate point in a greedy loop) should use
    !| [[f42_kd_tree_impl(module):kd_range_query_list_impl(subroutine)]] or
    !| [[f42_kd_tree_impl(module):kd_range_query_count_impl(subroutine)]] instead, since the
    !| `in_radius_mask = .false.` reset here costs O(n_points) on every call regardless of how
    !| few points are actually found.
    subroutine kd_range_query_mask_expert_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            tmp_range_stack,&
            in_radius_mask,&
            ierr&
        ) bind(C, name="kd_range_query_mask_expert_c")
        use f42_kd_tree, only: kd_range_query_mask_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        real(c_double), intent(in), target :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(3, n_points), intent(out), target :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        logical(c_bool), dimension(n_points), intent(out), target :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(radius)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_points)
        M_CHECK_ARRAY_NON_NULL(in_radius_mask, n_points)

        call kd_range_query_mask_expert(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            in_radius_mask = in_radius_mask,&
            ierr = ierr&
        )
    end subroutine kd_range_query_mask_expert_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_range_query_list(subroutine)]]
    !| Same traversal and pruning as [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]],
    !| but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
    !| instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
    !| in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
    !| needs the count, not the identities, should use
    !| [[f42_kd_tree_impl(module):kd_range_query_count_impl(subroutine)]] instead, to skip the
    !| index-buffer writes entirely.
    subroutine kd_range_query_list_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            neighbors,&
            n_found,&
            ierr&
        ) bind(C, name="kd_range_query_list_c")
        use f42_kd_tree, only: kd_range_query_list

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        real(c_double), intent(in), target :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(n_points), intent(out), target :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(c_int), intent(out), target :: n_found
            !! Output: number of points within `radius`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(radius)
        M_CHECK_NON_NULL(n_found)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(neighbors, n_points)

        call kd_range_query_list(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            neighbors = neighbors,&
            n_found = n_found,&
            ierr = ierr&
        )
    end subroutine kd_range_query_list_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_range_query_list_expert(subroutine)]]
    !| Same traversal and pruning as [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]],
    !| but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
    !| instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
    !| in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
    !| needs the count, not the identities, should use
    !| [[f42_kd_tree_impl(module):kd_range_query_count_impl(subroutine)]] instead, to skip the
    !| index-buffer writes entirely.
    subroutine kd_range_query_list_expert_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            tmp_range_stack,&
            neighbors,&
            n_found,&
            ierr&
        ) bind(C, name="kd_range_query_list_expert_c")
        use f42_kd_tree, only: kd_range_query_list_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        real(c_double), intent(in), target :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(3, n_points), intent(out), target :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(c_int), dimension(n_points), intent(out), target :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(c_int), intent(out), target :: n_found
            !! Output: number of points within `radius`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(radius)
        M_CHECK_NON_NULL(n_found)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_points)
        M_CHECK_ARRAY_NON_NULL(neighbors, n_points)

        call kd_range_query_list_expert(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            neighbors = neighbors,&
            n_found = n_found,&
            ierr = ierr&
        )
    end subroutine kd_range_query_list_expert_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_range_query_count(subroutine)]]
    !| Same traversal and pruning as
    !| [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]], but writes no index
    !| buffer at all -- only a scalar count. The right choice when the identities of the
    !| points found are never needed, e.g. a per-point local-density label.
    subroutine kd_range_query_count_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            neighbor_count,&
            ierr&
        ) bind(C, name="kd_range_query_count_c")
        use f42_kd_tree, only: kd_range_query_count

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        real(c_double), intent(in), target :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), intent(out), target :: neighbor_count
            !! Output: number of points within `radius`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(radius)
        M_CHECK_NON_NULL(neighbor_count)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)

        call kd_range_query_count(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            neighbor_count = neighbor_count,&
            ierr = ierr&
        )
    end subroutine kd_range_query_count_c

    !> summary: C-wrapper for [[f42_kd_tree(module):kd_range_query_count_expert(subroutine)]]
    !| Same traversal and pruning as
    !| [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]], but writes no index
    !| buffer at all -- only a scalar count. The right choice when the identities of the
    !| points found are never needed, e.g. a per-point local-density label.
    subroutine kd_range_query_count_expert_c(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            tmp_range_stack,&
            neighbor_count,&
            ierr&
        ) bind(C, name="kd_range_query_count_expert_c")
        use f42_kd_tree, only: kd_range_query_count_expert

        integer(c_int), intent(in), target :: n_dimensions
            !! Number of dimensions
        integer(c_int), intent(in), target :: n_points
            !! Number of points in the pre-built index
        real(c_double), dimension(n_dimensions, n_points), intent(in), target :: points
            !! Original points dataset
        integer(c_int), dimension(n_points), intent(in), target :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(c_int), dimension(n_dimensions), intent(in), target :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(c_double), dimension(n_dimensions), intent(in), target :: query_point
            !! Query point coordinates
        real(c_double), intent(in), target :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(c_int), dimension(3, n_points), intent(out), target :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(c_int), intent(out), target :: neighbor_count
            !! Output: number of points within `radius`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dimensions)
        M_CHECK_NON_NULL(n_points)
        M_CHECK_NON_NULL(radius)
        M_CHECK_NON_NULL(neighbor_count)
        M_CHECK_ARRAY_NON_NULL(points, n_dimensions * n_points)
        M_CHECK_ARRAY_NON_NULL(kd_indices, n_points)
        M_CHECK_ARRAY_NON_NULL(dimension_order, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(query_point, n_dimensions)
        M_CHECK_ARRAY_NON_NULL(tmp_range_stack, 3 * n_points)

        call kd_range_query_count_expert(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            neighbor_count = neighbor_count,&
            ierr = ierr&
        )
    end subroutine kd_range_query_count_expert_c

end module f42_kd_tree_c
#endif
