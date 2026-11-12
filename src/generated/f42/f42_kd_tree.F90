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
    use f42_kd_tree_impl, only: build_kd_index_impl, build_spherical_kd_impl, vicinity_vectors_count_impl, vicinity_vectors_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: build_kd_index
    public :: build_kd_index_expert
    public :: build_spherical_kd
    public :: build_spherical_kd_expert
    public :: vicinity_vectors
    public :: vicinity_vectors_expert
    public :: vicinity_vectors_count
    public :: vicinity_vectors_count_expert

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

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):vicinity_vectors_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):vicinity_vectors_expert]] to prepare it yourself.
    pure subroutine vicinity_vectors(&
            query_point,&
            points,&
            n_dimensions,&
            n_points,&
            r,&
            dimension_order,&
            kd_indices,&
            vicinity_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Total number of points organized in the k-d tree
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Coordinate vector used as the center of the search
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Ambient point matrix
        real(real64), intent(in) :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        logical(c_bool), dimension(n_points), intent(out) :: vicinity_mask
            !! Mask indicating points within the search radius
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_real(r, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=6_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=7_int32, min=1_int32, max=n_points)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack(3, 64))

        call vicinity_vectors_impl(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            tmp_stack = tmp_stack,&
            vicinity_mask = vicinity_mask&
        )
    end subroutine vicinity_vectors

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):vicinity_vectors_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):vicinity_vectors]] does both.
    pure subroutine vicinity_vectors_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Total number of points organized in the k-d tree
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Coordinate vector used as the center of the search
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Ambient point matrix
        real(real64), intent(in) :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(3, 64), intent(inout) :: tmp_stack
            !! Preallocated k-d tree traversal stack
        logical(c_bool), dimension(n_points), intent(out) :: vicinity_mask
            !! Mask indicating points within the search radius
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_real(r, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=6_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=7_int32, min=1_int32, max=n_points)
        if (is_err(ierr)) return
#endif

        call vicinity_vectors_impl(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            tmp_stack = tmp_stack,&
            vicinity_mask = vicinity_mask&
        )
    end subroutine vicinity_vectors_expert

    !> summary: Validates its inputs, prepares what [[f42_kd_tree_impl(module):vicinity_vectors_count_impl]] needs, then calls it. The entry point to reach for first; see [[f42_kd_tree(module):vicinity_vectors_count_expert]] to prepare it yourself.
    pure subroutine vicinity_vectors_count(&
            query_point,&
            points,&
            n_dimensions,&
            n_points,&
            r,&
            dimension_order,&
            kd_indices,&
            n_neighbors,&
            ierr&
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Total number of points organized in the k-d tree
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Coordinate vector used as the center of the search
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Ambient point matrix
        real(real64), intent(in) :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), intent(out) :: n_neighbors
            !! Number of points within the search radius
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:, :), allocatable :: tmp_stack

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_real(r, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=6_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=7_int32, min=1_int32, max=n_points)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(tmp_stack(3, 64))

        call vicinity_vectors_count_impl(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            tmp_stack = tmp_stack,&
            n_neighbors = n_neighbors&
        )
    end subroutine vicinity_vectors_count

    !> summary: Validates its inputs, then calls [[f42_kd_tree_impl(module):vicinity_vectors_count_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_kd_tree(module):vicinity_vectors_count]] does both.
    pure subroutine vicinity_vectors_count_expert(&
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
        )
        integer(int32), intent(in) :: n_dimensions
            !! Number of dimensions
        integer(int32), intent(in) :: n_points
            !! Total number of points organized in the k-d tree
        real(real64), dimension(n_dimensions), intent(in) :: query_point
            !! Coordinate vector used as the center of the search
        real(real64), dimension(n_dimensions, n_points), intent(in) :: points
            !! Ambient point matrix
        real(real64), intent(in) :: r
            !! Search radius
            !! The minimum valid value is `0.0_real64`.
        integer(int32), dimension(n_dimensions), intent(in) :: dimension_order
            !! Sequence of k-d tree split dimensions
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_dimensions`.
        integer(int32), dimension(n_points), intent(in) :: kd_indices
            !! K-d tree index sequence
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_points`.
        integer(int32), dimension(3, 64), intent(inout) :: tmp_stack
            !! Preallocated k-d tree traversal stack
        integer(int32), intent(out) :: n_neighbors
            !! Number of points within the search radius
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_dimensions, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_points, ierr, arg_pos=4_int32)
        call validate_in_range_real(r, ierr, arg_pos=5_int32, min=0.0_real64)
        call validate_all_in_range_real(query_point, n_dimensions, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(points, n_dimensions * n_points, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(dimension_order, n_dimensions, ierr, arg_pos=6_int32, min=1_int32, max=n_dimensions)
        call validate_all_in_range_int(kd_indices, n_points, ierr, arg_pos=7_int32, min=1_int32, max=n_points)
        if (is_err(ierr)) return
#endif

        call vicinity_vectors_count_impl(&
            query_point = query_point,&
            points = points,&
            n_dimensions = n_dimensions,&
            n_points = n_points,&
            r = r,&
            dimension_order = dimension_order,&
            kd_indices = kd_indices,&
            tmp_stack = tmp_stack,&
            n_neighbors = n_neighbors&
        )
    end subroutine vicinity_vectors_count_expert

end module f42_kd_tree
