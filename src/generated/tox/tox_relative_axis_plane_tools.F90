#include <src/macros.h>

!> summary: Wrappers for [[tox_relative_axis_plane_tools_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_relative_axis_plane_tools
    use tox_relative_axis_plane_tools_kernel, only: clock_hand_angle_between_vectors_kernel, clock_hand_angles_for_shift_vectors_kernel, compute_relative_axis_contributions_kernel, omics_field_RAP_projection_kernel
    use tox_relative_axis_plane_tools_kernel, only: omics_vector_RAP_projection_kernel, relative_axes_changes_from_shift_vector_kernel, relative_axes_expression_from_expression_vector_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_INVALID_INPUT, clear_err_arg_pos
    use tox_errors, only: set_err_once, validate_all_in_range_real, validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: omics_vector_RAP_projection
    public :: omics_field_RAP_projection
    public :: clock_hand_angle_between_vectors
    public :: clock_hand_angles_for_shift_vectors
    public :: compute_relative_axis_contributions
    public :: relative_axes_changes_from_shift_vector
    public :: relative_axes_expression_from_expression_vector

contains

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):omics_vector_RAP_projection_kernel]].
    subroutine omics_vector_RAP_projection(&
            vecs,&
            n_axes,&
            n_vecs,&
            vecs_selection_mask,&
            n_selected_vecs,&
            axes_selection_mask,&
            n_selected_axes,&
            projections,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_vecs
            !! number of vectors per axis
        integer(int32), intent(in) :: n_selected_vecs
            !! count of `.true.` values in `vecs_selection_mask`
        integer(int32), intent(in) :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
            !! matrix with expression vectors
        logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
            !! projected vectors
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vecs, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_vecs, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_selected_axes, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(vecs, n_axes * n_vecs, ierr, arg_pos=1_int32)
        if (count(vecs_selection_mask, kind=int32) /= n_selected_vecs) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (count(axes_selection_mask, kind=int32) /= n_selected_axes) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return

        call omics_vector_RAP_projection_kernel(&
            vecs = vecs,&
            n_axes = n_axes,&
            n_vecs = n_vecs,&
            vecs_selection_mask = vecs_selection_mask,&
            n_selected_vecs = n_selected_vecs,&
            axes_selection_mask = axes_selection_mask,&
            n_selected_axes = n_selected_axes,&
            projections = projections&
        )
    end subroutine omics_vector_RAP_projection

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):omics_field_RAP_projection_kernel]].
    subroutine omics_field_RAP_projection(&
            fields,&
            n_axes,&
            n_fields,&
            fields_selection_mask,&
            n_selected_fields,&
            axes_selection_mask,&
            n_selected_axes,&
            projections,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_fields
            !! number of vectors per axis
        integer(int32), intent(in) :: n_selected_fields
            !! count of `.true.` values in `fields_selection_mask`
        integer(int32), intent(in) :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(real64), dimension(n_axes, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(real64), dimension(n_selected_axes, n_selected_fields), intent(out) :: projections
            !! projected vectors
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_fields, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_fields, ierr, arg_pos=5_int32)
        call validate_dimension_size(n_selected_axes, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(fields, n_axes * 2 * n_fields, ierr, arg_pos=1_int32)
        if (count(fields_selection_mask, kind=int32) /= n_selected_fields) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (count(axes_selection_mask, kind=int32) /= n_selected_axes) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return

        call omics_field_RAP_projection_kernel(&
            fields = fields,&
            n_axes = n_axes,&
            n_fields = n_fields,&
            fields_selection_mask = fields_selection_mask,&
            n_selected_fields = n_selected_fields,&
            axes_selection_mask = axes_selection_mask,&
            n_selected_axes = n_selected_axes,&
            projections = projections&
        )
    end subroutine omics_field_RAP_projection

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):clock_hand_angle_between_vectors_kernel]].
    !| The unsigned angle is `acos(v1 . v2)`; `orientation_reference` supplies the sign by saying
    !| which way round the plane the two vectors span counts as positive. Reports
    !| `ERR_INVALID_INPUT` when the reference is orthogonal to the rotation and so orients nothing.
    subroutine clock_hand_angle_between_vectors(&
            v1,&
            v2,&
            n_dims,&
            orientation_reference,&
            signed_angle,&
            ierr&
        )
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(real64), intent(out) :: signed_angle
            !! Signed angle between vectors in radians [-pi, pi]
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_dims, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(v1, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(v2, n_dims, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(orientation_reference, n_dims, ierr, arg_pos=4_int32)
        if (is_err(ierr)) return

        call clock_hand_angle_between_vectors_kernel(&
            v1 = v1,&
            v2 = v2,&
            n_dims = n_dims,&
            orientation_reference = orientation_reference,&
            signed_angle = signed_angle,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine clock_hand_angle_between_vectors

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):clock_hand_angles_for_shift_vectors_kernel]].
    !| Each selected field is angled by the rule of
    !| [[tox_relative_axis_plane_tools_kernel(module):clock_hand_angle_between_vectors_kernel(subroutine)]],
    !| with one `orientation_reference` shared by the whole batch. A single field whose rotation
    !| the reference fails to orient fails the call.
    subroutine clock_hand_angles_for_shift_vectors(&
            fields,&
            n_dims,&
            n_fields,&
            fields_selection_mask,&
            n_selected_fields,&
            orientation_reference,&
            signed_angles,&
            ierr&
        )
        integer(int32), intent(in) :: n_dims
            !! Dimension of each vector in RAP space
        integer(int32), intent(in) :: n_fields
            !! Number of vector pairs
        integer(int32), intent(in) :: n_selected_fields
            !! Count of .true. values in fields_selection_mask
        real(real64), dimension(n_dims, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! .true. for vector pairs where angle should be computed
        real(real64), dimension(n_dims), intent(in) :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(real64), dimension(n_selected_fields), intent(out) :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_dims, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_fields, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_fields, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(fields, n_dims * 2 * n_fields, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(orientation_reference, n_dims, ierr, arg_pos=6_int32)
        if (count(fields_selection_mask, kind=int32) /= n_selected_fields) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (is_err(ierr)) return

        call clock_hand_angles_for_shift_vectors_kernel(&
            fields = fields,&
            n_dims = n_dims,&
            n_fields = n_fields,&
            fields_selection_mask = fields_selection_mask,&
            n_selected_fields = n_selected_fields,&
            orientation_reference = orientation_reference,&
            signed_angles = signed_angles,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine clock_hand_angles_for_shift_vectors

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):compute_relative_axis_contributions_kernel]].
    !| Shared utility: the shift-vector and expression-vector entry points below both drive it.
    subroutine compute_relative_axis_contributions(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (length of vec and contributions)
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized vector (expression or shift)
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(vec, n_axes, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return

        call compute_relative_axis_contributions_kernel(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine compute_relative_axis_contributions

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):relative_axes_changes_from_shift_vector_kernel]].
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    subroutine relative_axes_changes_from_shift_vector(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized shift vector
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(vec, n_axes, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return

        call relative_axes_changes_from_shift_vector_kernel(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine relative_axes_changes_from_shift_vector

    !> summary: Validates its inputs, then calls [[tox_relative_axis_plane_tools_kernel(module):relative_axes_expression_from_expression_vector_kernel]].
    !| Wrapper for single RAP-projected expression vectors
    subroutine relative_axes_expression_from_expression_vector(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized expression vector
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(vec, n_axes, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return

        call relative_axes_expression_from_expression_vector_kernel(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
        call clear_err_arg_pos(ierr)
    end subroutine relative_axes_expression_from_expression_vector

end module tox_relative_axis_plane_tools
