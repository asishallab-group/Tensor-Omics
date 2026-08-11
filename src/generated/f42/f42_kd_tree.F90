#include <src/macros.h>

!> summary: Wrappers for [[f42_kd_tree_impl(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module f42_kd_tree
    use f42_kd_tree_impl, only: build_kd_index_impl, build_spherical_kd_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: build_kd_index
    public :: build_kd_index_expert
    public :: build_spherical_kd
    public :: build_spherical_kd_expert

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

end module f42_kd_tree
