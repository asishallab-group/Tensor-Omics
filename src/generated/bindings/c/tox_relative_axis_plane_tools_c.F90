#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_relative_axis_plane_tools(module)]]
!| Relative axis planes (RAPs): planes through higher-dimensional gene expression space, and
!| what can be read off a vector once it is projected onto one.
!|
!| A RAP is picked by selecting axes (tissues) from the full expression space.
!| `vector_RAP_projection` projects a single vector onto it and `field_RAP_projection` a whole
!| field of them. Within the plane, `clock_hand_angle_between_vectors` measures the signed angle
!| between two vectors -- signed by an orientation reference, so the sign means the same thing in
!| every dimension -- and `clock_hand_angles_for_shift_vectors` does that for a whole shift
!| vector field at once.
!|
!| `relative_axes_changes_from_shift_vector` and `relative_axes_expression_from_expression_vector`
!| give the per-axis breakdown instead of the angle: how much of a change, or of an expression
!| level, falls on each selected axis.
module tox_relative_axis_plane_tools_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: omics_vector_RAP_projection_c
    public :: omics_field_RAP_projection_c
    public :: clock_hand_angle_between_vectors_c
    public :: clock_hand_angles_for_shift_vectors_c
    public :: compute_relative_axis_contributions_c
    public :: relative_axes_changes_from_shift_vector_c
    public :: relative_axes_expression_from_expression_vector_c

contains

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):omics_vector_RAP_projection(subroutine)]]
    subroutine omics_vector_RAP_projection_c(&
            vecs,&
            n_axes,&
            n_vecs,&
            vecs_selection_mask,&
            n_selected_vecs,&
            axes_selection_mask,&
            n_selected_axes,&
            projections,&
            ierr&
        ) bind(C, name="omics_vector_RAP_projection_c")
        use tox_relative_axis_plane_tools, only: omics_vector_RAP_projection

        integer(c_int), intent(in), target :: n_axes
            !! number of axes
        integer(c_int), intent(in), target :: n_vecs
            !! number of vectors per axis
        integer(c_int), intent(in), target :: n_selected_vecs
            !! count of `.true.` values in `vecs_selection_mask`
        integer(c_int), intent(in), target :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(c_double), dimension(n_axes, n_vecs), intent(in), target :: vecs
            !! matrix with expression vectors
        logical(c_bool), dimension(n_vecs), intent(in), target :: vecs_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical(c_bool), dimension(n_axes), intent(in), target :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(c_double), dimension(n_selected_axes, n_selected_vecs), intent(out), target :: projections
            !! projected vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_vecs)
        M_CHECK_NON_NULL(n_selected_vecs)
        M_CHECK_NON_NULL(n_selected_axes)
        M_CHECK_ARRAY_NON_NULL(vecs, n_axes * n_vecs)
        M_CHECK_ARRAY_NON_NULL(vecs_selection_mask, n_vecs)
        M_CHECK_ARRAY_NON_NULL(axes_selection_mask, n_axes)
        M_CHECK_ARRAY_NON_NULL(projections, n_selected_axes * n_selected_vecs)

        call omics_vector_RAP_projection(&
            vecs = vecs,&
            n_axes = n_axes,&
            n_vecs = n_vecs,&
            vecs_selection_mask = vecs_selection_mask,&
            n_selected_vecs = n_selected_vecs,&
            axes_selection_mask = axes_selection_mask,&
            n_selected_axes = n_selected_axes,&
            projections = projections,&
            ierr = ierr&
        )
    end subroutine omics_vector_RAP_projection_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):omics_field_RAP_projection(subroutine)]]
    subroutine omics_field_RAP_projection_c(&
            fields,&
            n_axes,&
            n_fields,&
            fields_selection_mask,&
            n_selected_fields,&
            axes_selection_mask,&
            n_selected_axes,&
            projections,&
            ierr&
        ) bind(C, name="omics_field_RAP_projection_c")
        use tox_relative_axis_plane_tools, only: omics_field_RAP_projection

        integer(c_int), intent(in), target :: n_axes
            !! number of axes
        integer(c_int), intent(in), target :: n_fields
            !! number of vectors per axis
        integer(c_int), intent(in), target :: n_selected_fields
            !! count of `.true.` values in `fields_selection_mask`
        integer(c_int), intent(in), target :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(c_double), dimension(n_axes, 2, n_fields), intent(in), target :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical(c_bool), dimension(n_fields), intent(in), target :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical(c_bool), dimension(n_axes), intent(in), target :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(c_double), dimension(n_selected_axes, n_selected_fields), intent(out), target :: projections
            !! projected vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_fields)
        M_CHECK_NON_NULL(n_selected_fields)
        M_CHECK_NON_NULL(n_selected_axes)
        M_CHECK_ARRAY_NON_NULL(fields, n_axes * 2 * n_fields)
        M_CHECK_ARRAY_NON_NULL(fields_selection_mask, n_fields)
        M_CHECK_ARRAY_NON_NULL(axes_selection_mask, n_axes)
        M_CHECK_ARRAY_NON_NULL(projections, n_selected_axes * n_selected_fields)

        call omics_field_RAP_projection(&
            fields = fields,&
            n_axes = n_axes,&
            n_fields = n_fields,&
            fields_selection_mask = fields_selection_mask,&
            n_selected_fields = n_selected_fields,&
            axes_selection_mask = axes_selection_mask,&
            n_selected_axes = n_selected_axes,&
            projections = projections,&
            ierr = ierr&
        )
    end subroutine omics_field_RAP_projection_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):clock_hand_angle_between_vectors(subroutine)]]
    !| The unsigned angle is `acos(v1 . v2)`; `orientation_reference` supplies the sign by saying
    !| which way round the plane the two vectors span counts as positive. Reports
    !| `ERR_INVALID_INPUT` when the reference is orthogonal to the rotation and so orients nothing.
    subroutine clock_hand_angle_between_vectors_c(&
            v1,&
            v2,&
            n_dims,&
            orientation_reference,&
            signed_angle,&
            ierr&
        ) bind(C, name="clock_hand_angle_between_vectors_c")
        use tox_relative_axis_plane_tools, only: clock_hand_angle_between_vectors

        integer(c_int), intent(in), target :: n_dims
            !! Dimension of both vectors
        real(c_double), dimension(n_dims), intent(in), target :: v1
            !! First normalized vector in RAP space
        real(c_double), dimension(n_dims), intent(in), target :: v2
            !! Second normalized vector in RAP space
        real(c_double), dimension(n_dims), intent(in), target :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(c_double), intent(out), target :: signed_angle
            !! Signed angle between vectors in radians [-pi, pi]
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(signed_angle)
        M_CHECK_ARRAY_NON_NULL(v1, n_dims)
        M_CHECK_ARRAY_NON_NULL(v2, n_dims)
        M_CHECK_ARRAY_NON_NULL(orientation_reference, n_dims)

        call clock_hand_angle_between_vectors(&
            v1 = v1,&
            v2 = v2,&
            n_dims = n_dims,&
            orientation_reference = orientation_reference,&
            signed_angle = signed_angle,&
            ierr = ierr&
        )
    end subroutine clock_hand_angle_between_vectors_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):clock_hand_angles_for_shift_vectors(subroutine)]]
    !| Each selected field is angled by the rule of
    !| [[tox_relative_axis_plane_tools_impl(module):clock_hand_angle_between_vectors_impl(subroutine)]],
    !| with one `orientation_reference` shared by the whole batch. A single field whose rotation
    !| the reference fails to orient fails the call.
    subroutine clock_hand_angles_for_shift_vectors_c(&
            fields,&
            n_dims,&
            n_fields,&
            fields_selection_mask,&
            n_selected_fields,&
            orientation_reference,&
            signed_angles,&
            ierr&
        ) bind(C, name="clock_hand_angles_for_shift_vectors_c")
        use tox_relative_axis_plane_tools, only: clock_hand_angles_for_shift_vectors

        integer(c_int), intent(in), target :: n_dims
            !! Dimension of each vector in RAP space
        integer(c_int), intent(in), target :: n_fields
            !! Number of vector pairs
        integer(c_int), intent(in), target :: n_selected_fields
            !! Count of .true. values in fields_selection_mask
        real(c_double), dimension(n_dims, 2, n_fields), intent(in), target :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical(c_bool), dimension(n_fields), intent(in), target :: fields_selection_mask
            !! .true. for vector pairs where angle should be computed
        real(c_double), dimension(n_dims), intent(in), target :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(c_double), dimension(n_selected_fields), intent(out), target :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(n_fields)
        M_CHECK_NON_NULL(n_selected_fields)
        M_CHECK_ARRAY_NON_NULL(fields, n_dims * 2 * n_fields)
        M_CHECK_ARRAY_NON_NULL(fields_selection_mask, n_fields)
        M_CHECK_ARRAY_NON_NULL(orientation_reference, n_dims)
        M_CHECK_ARRAY_NON_NULL(signed_angles, n_selected_fields)

        call clock_hand_angles_for_shift_vectors(&
            fields = fields,&
            n_dims = n_dims,&
            n_fields = n_fields,&
            fields_selection_mask = fields_selection_mask,&
            n_selected_fields = n_selected_fields,&
            orientation_reference = orientation_reference,&
            signed_angles = signed_angles,&
            ierr = ierr&
        )
    end subroutine clock_hand_angles_for_shift_vectors_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):compute_relative_axis_contributions(subroutine)]]
    !| Shared utility: the shift-vector and expression-vector entry points below both drive it.
    subroutine compute_relative_axis_contributions_c(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
        ) bind(C, name="compute_relative_axis_contributions_c")
        use tox_relative_axis_plane_tools, only: compute_relative_axis_contributions

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (length of vec and contributions)
        real(c_double), dimension(n_axes), intent(in), target :: vec
            !! RAP-projected and normalized vector (expression or shift)
        real(c_double), dimension(n_axes), intent(out), target :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_ARRAY_NON_NULL(vec, n_axes)
        M_CHECK_ARRAY_NON_NULL(contributions, n_axes)

        call compute_relative_axis_contributions(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
    end subroutine compute_relative_axis_contributions_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):relative_axes_changes_from_shift_vector(subroutine)]]
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    subroutine relative_axes_changes_from_shift_vector_c(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
        ) bind(C, name="relative_axes_changes_from_shift_vector_c")
        use tox_relative_axis_plane_tools, only: relative_axes_changes_from_shift_vector

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes
        real(c_double), dimension(n_axes), intent(in), target :: vec
            !! RAP-projected and normalized shift vector
        real(c_double), dimension(n_axes), intent(out), target :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_ARRAY_NON_NULL(vec, n_axes)
        M_CHECK_ARRAY_NON_NULL(contributions, n_axes)

        call relative_axes_changes_from_shift_vector(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
    end subroutine relative_axes_changes_from_shift_vector_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):relative_axes_expression_from_expression_vector(subroutine)]]
    !| Wrapper for single RAP-projected expression vectors
    subroutine relative_axes_expression_from_expression_vector_c(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
        ) bind(C, name="relative_axes_expression_from_expression_vector_c")
        use tox_relative_axis_plane_tools, only: relative_axes_expression_from_expression_vector

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes
        real(c_double), dimension(n_axes), intent(in), target :: vec
            !! RAP-projected and normalized expression vector
        real(c_double), dimension(n_axes), intent(out), target :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_ARRAY_NON_NULL(vec, n_axes)
        M_CHECK_ARRAY_NON_NULL(contributions, n_axes)

        call relative_axes_expression_from_expression_vector(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
    end subroutine relative_axes_expression_from_expression_vector_c

end module tox_relative_axis_plane_tools_c
#endif
