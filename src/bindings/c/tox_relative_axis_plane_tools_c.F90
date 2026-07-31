#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_relative_axis_plane_tools(module)]]
!| Module for tools related to relative axis planes (RAPs), i.e. planes in higher-dimensional gene expression space
module tox_relative_axis_plane_tools_c
    use safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, is_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: omics_vector_RAP_projection_c
    public :: omics_field_RAP_projection_c
    public :: clock_hand_angle_between_vectors_c
    public :: clock_hand_angles_for_shift_vectors_c
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
            !! Error code
        logical, dimension(n_vecs) :: vecs_selection_mask_f
        logical, dimension(n_axes) :: axes_selection_mask_f

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

        vecs_selection_mask_f = vecs_selection_mask
        axes_selection_mask_f = axes_selection_mask

        call omics_vector_RAP_projection(&
            vecs = vecs,&
            n_axes = n_axes,&
            n_vecs = n_vecs,&
            vecs_selection_mask = vecs_selection_mask_f,&
            n_selected_vecs = n_selected_vecs,&
            axes_selection_mask = axes_selection_mask_f,&
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
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        logical(c_bool), dimension(n_fields), intent(in), target :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical(c_bool), dimension(n_axes), intent(in), target :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(c_double), dimension(n_selected_axes, n_selected_fields), intent(out), target :: projections
            !! projected vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, dimension(n_fields) :: fields_selection_mask_f
        logical, dimension(n_axes) :: axes_selection_mask_f

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

        fields_selection_mask_f = fields_selection_mask
        axes_selection_mask_f = axes_selection_mask

        call omics_field_RAP_projection(&
            fields = fields,&
            n_axes = n_axes,&
            n_fields = n_fields,&
            fields_selection_mask = fields_selection_mask_f,&
            n_selected_fields = n_selected_fields,&
            axes_selection_mask = axes_selection_mask_f,&
            n_selected_axes = n_selected_axes,&
            projections = projections,&
            ierr = ierr&
        )
    end subroutine omics_field_RAP_projection_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):clock_hand_angle_between_vectors(subroutine)]]
    !| Calculates the signed rotation angle between two normalized vectors in RAP space.
    !| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
    subroutine clock_hand_angle_between_vectors_c(&
            v1,&
            v2,&
            n_dims,&
            signed_angle,&
            selected_axes_for_signed,&
            ierr&
        ) bind(C, name="clock_hand_angle_between_vectors_c")
        use tox_relative_axis_plane_tools, only: clock_hand_angle_between_vectors

        integer(c_int), intent(in), target :: n_dims
            !! Dimension of both vectors
        real(c_double), dimension(n_dims), intent(in), target :: v1
            !! First normalized vector in RAP space
        real(c_double), dimension(n_dims), intent(in), target :: v2
            !! Second normalized vector in RAP space
        real(c_double), intent(out), target :: signed_angle
            !! Signed angle between vectors in radians [-π, π]
        integer(c_int), dimension(3), intent(in), target :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
        integer(c_int), intent(out), target :: ierr
            !! Error code

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(signed_angle)
        M_CHECK_ARRAY_NON_NULL(v1, n_dims)
        M_CHECK_ARRAY_NON_NULL(v2, n_dims)
        M_CHECK_ARRAY_NON_NULL(selected_axes_for_signed, 3)

        call clock_hand_angle_between_vectors(&
            v1 = v1,&
            v2 = v2,&
            n_dims = n_dims,&
            signed_angle = signed_angle,&
            selected_axes_for_signed = selected_axes_for_signed,&
            ierr = ierr&
        )
    end subroutine clock_hand_angle_between_vectors_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):clock_hand_angles_for_shift_vectors(subroutine)]]
    subroutine clock_hand_angles_for_shift_vectors_c(&
            fields,&
            n_dims,&
            n_fields,&
            fields_selection_mask,&
            n_selected_fields,&
            selected_axes_for_signed,&
            signed_angles,&
            ierr&
        ) bind(C, name="clock_hand_angles_for_shift_vectors_c")
        use tox_relative_axis_plane_tools, only: clock_hand_angles_for_shift_vectors

        integer(c_int), intent(in), target :: n_dims
            !! Dimension of each vector in RAP space
        integer(c_int), intent(in), target :: n_fields
            !! Number of vector fields
        integer(c_int), intent(in), target :: n_selected_fields
            !! Count of .true. values in fields_selection_mask
        real(c_double), dimension(n_dims, 2, n_fields), intent(in), target :: fields
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        logical(c_bool), dimension(n_fields), intent(in), target :: fields_selection_mask
            !! .true. for vector pairs where angle should be computed
        integer(c_int), dimension(3), intent(in), target :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
        real(c_double), dimension(n_selected_fields), intent(out), target :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, dimension(n_fields) :: fields_selection_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(n_fields)
        M_CHECK_NON_NULL(n_selected_fields)
        M_CHECK_ARRAY_NON_NULL(fields, n_dims * 2 * n_fields)
        M_CHECK_ARRAY_NON_NULL(fields_selection_mask, n_fields)
        M_CHECK_ARRAY_NON_NULL(selected_axes_for_signed, 3)
        M_CHECK_ARRAY_NON_NULL(signed_angles, n_selected_fields)

        fields_selection_mask_f = fields_selection_mask

        call clock_hand_angles_for_shift_vectors(&
            fields = fields,&
            n_dims = n_dims,&
            n_fields = n_fields,&
            fields_selection_mask = fields_selection_mask_f,&
            n_selected_fields = n_selected_fields,&
            selected_axes_for_signed = selected_axes_for_signed,&
            signed_angles = signed_angles,&
            ierr = ierr&
        )
    end subroutine clock_hand_angles_for_shift_vectors_c

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
