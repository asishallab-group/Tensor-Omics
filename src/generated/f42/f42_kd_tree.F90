#include <src/macros.h>

!> k-d tree spatial index over fixed-dimensional point sets.
!| Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
!| caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
!| is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
!| index array rather than as linked nodes.
!|
!| Generated from [[f42_kd_tree_impl(module)]]; do not edit -- regenerate instead.
module f42_kd_tree
    use f42_safeguard
    use f42_kd_tree_impl, only: build_kd_index_impl, build_spherical_kd_impl, kd_knn_query_impl, kd_range_query_count_impl
    use f42_kd_tree_impl, only: kd_range_query_list_impl, kd_range_query_mask_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: build_kd_index
    public :: build_kd_index_expert
    public :: build_spherical_kd
    public :: build_spherical_kd_expert
    public :: kd_knn_query
    public :: kd_knn_query_expert
    public :: kd_range_query_mask
    public :: kd_range_query_mask_expert
    public :: kd_range_query_list
    public :: kd_range_query_list_expert
    public :: kd_range_query_count
    public :: kd_range_query_count_expert

contains

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):build_kd_index_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):build_kd_index_expert]] to prepare it yourself.
    pure subroutine build_kd_index(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_workspace
        real(real64), dimension(:), allocatable :: tmp_value_buffer
        integer(int32), dimension(:), allocatable :: tmp_permutation
        integer(int32), dimension(:, :), allocatable :: tmp_recursion_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_workspace(n_points))
        M_ALLOCATE(tmp_value_buffer(n_points))
        M_ALLOCATE(tmp_permutation(n_points))
        M_ALLOCATE(tmp_recursion_stack(3, n_points))

        call build_kd_index_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_recursion_stack = tmp_recursion_stack&
        )
    end subroutine build_kd_index

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):build_kd_index_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):build_kd_index]] does both.
    pure subroutine build_kd_index_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_points), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(n_points), intent(out) :: tmp_value_buffer
            !! Workspace for sorting
        integer(int32), dimension(n_points), intent(out) :: tmp_permutation
            !! Workspace for sorting
        integer(int32), dimension(3, n_points), intent(out) :: tmp_recursion_stack
            !! Stack for l, r, depth
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call build_kd_index_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_recursion_stack = tmp_recursion_stack&
        )
    end subroutine build_kd_index_expert

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):build_spherical_kd_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):build_spherical_kd_expert]] to prepare it yourself.
    !| This is a thin, semantically-named wrapper: partitioning is identical to
    !| [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]] (plain per-axis median splits);
    !| callers are responsible for ensuring `points` are actually unit-normalized beforehand.
    pure subroutine build_spherical_kd(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: tmp_workspace
        real(real64), dimension(:), allocatable :: tmp_value_buffer
        integer(int32), dimension(:), allocatable :: tmp_permutation
        integer(int32), dimension(:, :), allocatable :: tmp_recursion_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_workspace(n_points))
        M_ALLOCATE(tmp_value_buffer(n_points))
        M_ALLOCATE(tmp_permutation(n_points))
        M_ALLOCATE(tmp_recursion_stack(3, n_points))

        call build_spherical_kd_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_recursion_stack = tmp_recursion_stack&
        )
    end subroutine build_spherical_kd

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):build_spherical_kd_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):build_spherical_kd]] does both.
    !| This is a thin, semantically-named wrapper: partitioning is identical to
    !| [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]] (plain per-axis median splits);
    !| callers are responsible for ensuring `points` are actually unit-normalized beforehand.
    pure subroutine build_spherical_kd_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Data points
        integer(int32), dimension(n_points), intent(out) :: kd_indices
            !! Output index array (k-d tree order)
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order (by variance)
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_points), intent(out) :: tmp_workspace
            !! Workspace array
        real(real64), dimension(n_points), intent(out) :: tmp_value_buffer
            !! Workspace for sorting
        integer(int32), dimension(n_points), intent(out) :: tmp_permutation
            !! Workspace for sorting
        integer(int32), dimension(3, n_points), intent(out) :: tmp_recursion_stack
            !! Stack for l, r, depth
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        if (is_err(ierr)) return
#endif

        call build_spherical_kd_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            tmp_workspace = tmp_workspace,&
            tmp_value_buffer = tmp_value_buffer,&
            tmp_permutation = tmp_permutation,&
            tmp_recursion_stack = tmp_recursion_stack&
        )
    end subroutine build_spherical_kd_expert

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):kd_knn_query_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):kd_knn_query_expert]] to prepare it yourself.
    !| Via a bounded max-heap kept directly in `neighbors`/`distances` and splitting-plane
    !| pruning. Requires `k_neighbors <= n_points`: every point is then guaranteed visited
    !| before the heap can still have room, so no fallback for "fewer than k found" is needed.
    !| Does not guarantee `neighbors`/`distances` are sorted nearest-to-farthest (max-heap
    !| order internally).
    pure subroutine kd_knn_query(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        integer(int32), intent(in) :: k_neighbors
            !! Number of neighbors to find
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index, see [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]]
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        integer(int32), dimension(k_neighbors), intent(out) :: neighbors
            !! Output: indices of the k nearest neighbors (nearest-to-farthest order not guaranteed, max-heap order internally)
        real(real64), dimension(k_neighbors), intent(out) :: distances
            !! Output: Euclidean distances to the k nearest neighbors
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_int(k_neighbors, ierr, arg_pos=7_int32, min=1_int32, max=n_points)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_knn_query_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            k_neighbors = k_neighbors,&
            tmp_range_stack = tmp_range_stack,&
            neighbors = neighbors,&
            distances = distances&
        )
    end subroutine kd_knn_query

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):kd_knn_query_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):kd_knn_query]] does both.
    !| Via a bounded max-heap kept directly in `neighbors`/`distances` and splitting-plane
    !| pruning. Requires `k_neighbors <= n_points`: every point is then guaranteed visited
    !| before the heap can still have room, so no fallback for "fewer than k found" is needed.
    !| Does not guarantee `neighbors`/`distances` are sorted nearest-to-farthest (max-heap
    !| order internally).
    pure subroutine kd_knn_query_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        integer(int32), intent(in) :: k_neighbors
            !! Number of neighbors to find
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index, see [[f42_kd_tree_impl(module):build_kd_index_impl(subroutine)]]
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), dimension(k_neighbors), intent(out) :: neighbors
            !! Output: indices of the k nearest neighbors (nearest-to-farthest order not guaranteed, max-heap order internally)
        real(real64), dimension(k_neighbors), intent(out) :: distances
            !! Output: Euclidean distances to the k nearest neighbors
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_int(k_neighbors, ierr, arg_pos=7_int32, min=1_int32, max=n_points)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        call kd_knn_query_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            k_neighbors = k_neighbors,&
            tmp_range_stack = tmp_range_stack,&
            neighbors = neighbors,&
            distances = distances&
        )
    end subroutine kd_knn_query_expert

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):kd_range_query_mask_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):kd_range_query_mask_expert]] to prepare it yourself.
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
    pure subroutine kd_range_query_mask(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            in_radius_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        logical(c_bool), dimension(n_points), intent(out) :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_range_query_mask_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            in_radius_mask = in_radius_mask&
        )
    end subroutine kd_range_query_mask

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):kd_range_query_mask_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):kd_range_query_mask]] does both.
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
    pure subroutine kd_range_query_mask_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        logical(c_bool), dimension(n_points), intent(out) :: in_radius_mask
            !! Output: .true. for points within `radius` of `query_point`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        call kd_range_query_mask_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            in_radius_mask = in_radius_mask&
        )
    end subroutine kd_range_query_mask_expert

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):kd_range_query_list_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):kd_range_query_list_expert]] to prepare it yourself.
    !| Same traversal and pruning as [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]],
    !| but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
    !| instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
    !| in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
    !| needs the count, not the identities, should use
    !| [[f42_kd_tree_impl(module):kd_range_query_count_impl(subroutine)]] instead, to skip the
    !| index-buffer writes entirely.
    pure subroutine kd_range_query_list(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_points), intent(out) :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(int32), intent(out) :: n_found
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_range_query_list_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            neighbors = neighbors,&
            n_found = n_found&
        )
    end subroutine kd_range_query_list

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):kd_range_query_list_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):kd_range_query_list]] does both.
    !| Same traversal and pruning as [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]],
    !| but writes matches into a caller-provided compact index buffer (`neighbors(1:n_found)`)
    !| instead of a full-size logical mask, so repeated calls -- e.g. once per candidate point
    !| in an outer greedy loop -- don't each pay an O(n_points) reset. A caller that only
    !| needs the count, not the identities, should use
    !| [[f42_kd_tree_impl(module):kd_range_query_count_impl(subroutine)]] instead, to skip the
    !| index-buffer writes entirely.
    pure subroutine kd_range_query_list_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), dimension(n_points), intent(out) :: neighbors
            !! Output: indices within `radius`, valid in `neighbors(1:n_found)`
        integer(int32), intent(out) :: n_found
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        call kd_range_query_list_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            neighbors = neighbors,&
            n_found = n_found&
        )
    end subroutine kd_range_query_list_expert

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):kd_range_query_count_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):kd_range_query_count_expert]] to prepare it yourself.
    !| Same traversal and pruning as
    !| [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]], but writes no index
    !| buffer at all -- only a scalar count. The right choice when the identities of the
    !| points found are never needed, e.g. a per-point local-density label.
    pure subroutine kd_range_query_count(&
            points,&
            n_dimensions,&
            n_points,&
            kd_indices,&
            dimension_order,&
            query_point,&
            radius,&
            neighbor_count,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), intent(out) :: neighbor_count
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_range_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_range_stack(3, n_points))

        call kd_range_query_count_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            neighbor_count = neighbor_count&
        )
    end subroutine kd_range_query_count

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):kd_range_query_count_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):kd_range_query_count]] does both.
    !| Same traversal and pruning as
    !| [[f42_kd_tree_impl(module):kd_range_query_mask_impl(subroutine)]], but writes no index
    !| buffer at all -- only a scalar count. The right choice when the identities of the
    !| points found are never needed, e.g. a per-point local-density label.
    pure subroutine kd_range_query_count_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Number of points in the pre-built index
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Original points dataset
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! Pre-built k-d tree index
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Dimension order used to build `kd_indices`
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Query point coordinates
        real(real64), intent(in) :: radius
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(3, n_points), intent(out) :: tmp_range_stack
            !! Workspace: traversal stack for [left_idx, right_idx, depth] frames
        integer(int32), intent(out) :: neighbor_count
            !! Output: number of points within `radius`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=3_int32)
        call validate_in_range_real(radius, ierr, arg_pos=7_int32, min=0.0_real64)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=4_int32, min=1_int32, max=n_points)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=5_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=6_int32)
        if (is_err(ierr)) return
#endif

        call kd_range_query_count_impl(&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            kd_indices = kd_indices,&
            dimension_order = dimension_order,&
            query_point = query_point,&
            radius = radius,&
            tmp_range_stack = tmp_range_stack,&
            neighbor_count = neighbor_count&
        )
    end subroutine kd_range_query_count_expert

end module f42_kd_tree
