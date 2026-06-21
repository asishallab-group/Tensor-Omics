#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_relative_axis_plane_tools(module)]]
!| Module for tools related to relative axis planes (RAPs), i.e. planes in higher-dimensional gene expression space
module tox_relative_axis_plane_tools_c
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

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):omics_vector_RAP_projection(subroutine)]]
    !| Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
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
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_axes
            !! number of axes
        integer(c_int), intent(in), target :: n_vecs
            !! number of vectors per axis
        integer(c_int), intent(in), target :: n_selected_vecs
            !! count of `.true.` values in `vecs_selection_mask`
        integer(c_int), intent(in), target :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(c_double), intent(in), dimension(n_axes, n_vecs), target :: vecs
            !! matrix with expression vectors
        integer(c_int), intent(in), dimension(n_vecs), target :: vecs_selection_mask
            !! `.true.` for vectors where projection is to be computed
        integer(c_int), intent(in), dimension(n_axes), target :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(c_double), intent(out), dimension(n_selected_axes, n_selected_vecs), target :: projections
            !! projected vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, allocatable, dimension(:) :: vecs_selection_mask_f
        logical, allocatable, dimension(:) :: axes_selection_mask_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vecs)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_vecs)
        M_CHECK_NON_NULL(vecs_selection_mask)
        M_CHECK_NON_NULL(n_selected_vecs)
        M_CHECK_NON_NULL(axes_selection_mask)
        M_CHECK_NON_NULL(n_selected_axes)
        M_CHECK_NON_NULL(projections)
        M_ALLOCATE(vecs_selection_mask_f(n_vecs))
        call c_int_as_logical(vecs_selection_mask, vecs_selection_mask_f)
        M_ALLOCATE(axes_selection_mask_f(n_axes))
        call c_int_as_logical(axes_selection_mask, axes_selection_mask_f)
        call  omics_vector_RAP_projection(&
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
    !| Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    subroutine omics_field_RAP_projection_c(&
            vecs,&
            n_axes,&
            n_vecs,&
            vecs_selection_mask,&
            n_selected_vecs,&
            axes_selection_mask,&
            n_selected_axes,&
            projections,&
            ierr&
            ) bind(C, name="omics_field_RAP_projection_c")
        use tox_relative_axis_plane_tools, only: omics_field_RAP_projection
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_axes
            !! number of axes
        integer(c_int), intent(in), target :: n_vecs
            !! number of vectors per axis
        integer(c_int), intent(in), target :: n_selected_vecs
            !! count of `.true.` values in `vecs_selection_mask`
        integer(c_int), intent(in), target :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(c_double), intent(in), dimension(2 * n_axes, n_vecs), target :: vecs
            !! matrix with vector fields, first n rows mean vector origin, last n rows vector targets
        integer(c_int), intent(in), dimension(n_vecs), target :: vecs_selection_mask
            !! `.true.` for vectors where projection is to be computed
        integer(c_int), intent(in), dimension(n_axes), target :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(c_double), intent(out), dimension(n_selected_axes, n_selected_vecs), target :: projections
            !! projected vectors
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, allocatable, dimension(:) :: vecs_selection_mask_f
        logical, allocatable, dimension(:) :: axes_selection_mask_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vecs)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_vecs)
        M_CHECK_NON_NULL(vecs_selection_mask)
        M_CHECK_NON_NULL(n_selected_vecs)
        M_CHECK_NON_NULL(axes_selection_mask)
        M_CHECK_NON_NULL(n_selected_axes)
        M_CHECK_NON_NULL(projections)
        M_ALLOCATE(vecs_selection_mask_f(n_vecs))
        call c_int_as_logical(vecs_selection_mask, vecs_selection_mask_f)
        M_ALLOCATE(axes_selection_mask_f(n_axes))
        call c_int_as_logical(axes_selection_mask, axes_selection_mask_f)
        call  omics_field_RAP_projection(&
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
    end subroutine omics_field_RAP_projection_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):clock_hand_angle_between_vectors(subroutine)]]
    !| Compute the signed clock hand angle between two RAP-projected and normalized vectors.
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
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_dims
            !! Dimension of both vectors
        real(c_double), intent(in), dimension(n_dims), target :: v1
            !! First normalized vector in RAP space
        real(c_double), intent(in), dimension(n_dims), target :: v2
            !! Second normalized vector in RAP space
        real(c_double), intent(out), target :: signed_angle
            !! Signed angle between vectors in radians [-π, π]
        integer(c_int), intent(in), dimension(3), target :: selected_axes_for_signed
            !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(v1)
        M_CHECK_NON_NULL(v2)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(signed_angle)
        M_CHECK_NON_NULL(selected_axes_for_signed)
        call  clock_hand_angle_between_vectors(&
            v1 = v1,&
            v2 = v2,&
            n_dims = n_dims,&
            signed_angle = signed_angle,&
            selected_axes_for_signed = selected_axes_for_signed,&
            ierr = ierr&
        )
    end subroutine clock_hand_angle_between_vectors_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):clock_hand_angles_for_shift_vectors(subroutine)]]
    !| Compute signed rotation angles between RAP-projected and normalized vector pairs.
    !| Takes separate arrays of RAP-projected and normalized vectors (e.g. expression centroids and paralogs) and computes the signed rotation angle between corresponding pairs.
    !| This measures both magnitude and directionality of angular separation in RAP space.
    subroutine clock_hand_angles_for_shift_vectors_c(&
            origins,&
            targets,&
            n_dims,&
            n_vecs,&
            vecs_selection_mask,&
            n_selected_vecs,&
            selected_axes_for_signed,&
            signed_angles,&
            ierr&
            ) bind(C, name="clock_hand_angles_for_shift_vectors_c")
        use tox_relative_axis_plane_tools, only: clock_hand_angles_for_shift_vectors
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_dims
            !! Dimension of each vector in RAP space
        integer(c_int), intent(in), target :: n_vecs
            !! Number of vector pairs
        integer(c_int), intent(in), target :: n_selected_vecs
            !! Count of .true. values in vecs_selection_mask
        real(c_double), intent(in), dimension(n_dims, n_vecs), target :: origins
            !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
        real(c_double), intent(in), dimension(n_dims, n_vecs), target :: targets
            !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
        integer(c_int), intent(in), dimension(n_vecs), target :: vecs_selection_mask
            !! .true. for vector pairs where angle should be computed
        integer(c_int), intent(in), dimension(3), target :: selected_axes_for_signed
            !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
        real(c_double), intent(out), dimension(n_selected_vecs), target :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]
        integer(c_int), intent(out), target :: ierr
            !! Error code
        logical, allocatable, dimension(:) :: vecs_selection_mask_f
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(origins)
        M_CHECK_NON_NULL(targets)
        M_CHECK_NON_NULL(n_dims)
        M_CHECK_NON_NULL(n_vecs)
        M_CHECK_NON_NULL(vecs_selection_mask)
        M_CHECK_NON_NULL(n_selected_vecs)
        M_CHECK_NON_NULL(selected_axes_for_signed)
        M_CHECK_NON_NULL(signed_angles)
        M_ALLOCATE(vecs_selection_mask_f(n_vecs))
        call c_int_as_logical(vecs_selection_mask, vecs_selection_mask_f)
        call  clock_hand_angles_for_shift_vectors(&
            origins = origins,&
            targets = targets,&
            n_dims = n_dims,&
            n_vecs = n_vecs,&
            vecs_selection_mask = vecs_selection_mask_f,&
            n_selected_vecs = n_selected_vecs,&
            selected_axes_for_signed = selected_axes_for_signed,&
            signed_angles = signed_angles,&
            ierr = ierr&
        )
    end subroutine clock_hand_angles_for_shift_vectors_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):compute_relative_axis_contributions(subroutine)]]
    !| Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
    subroutine compute_relative_axis_contributions_c(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
            ) bind(C, name="compute_relative_axis_contributions_c")
        use tox_relative_axis_plane_tools, only: compute_relative_axis_contributions
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (length of vec and contributions)
        real(c_double), intent(in), dimension(n_axes), target :: vec
            !! RAP-projected and normalized vector (expression or shift)
        real(c_double), intent(out), dimension(n_axes), target :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vec)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(contributions)
        call  compute_relative_axis_contributions(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
    end subroutine compute_relative_axis_contributions_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):relative_axes_changes_from_shift_vector(subroutine)]]
    !| Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    subroutine relative_axes_changes_from_shift_vector_c(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
            ) bind(C, name="relative_axes_changes_from_shift_vector_c")
        use tox_relative_axis_plane_tools, only: relative_axes_changes_from_shift_vector
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_axes
            !! Number of axes
        real(c_double), intent(in), dimension(n_axes), target :: vec
            !! RAP-projected and normalized shift vector
        real(c_double), intent(out), dimension(n_axes), target :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vec)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(contributions)
        call  relative_axes_changes_from_shift_vector(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
    end subroutine relative_axes_changes_from_shift_vector_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):relative_axes_expression_from_expression_vector(subroutine)]]
    !| Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
    !| Wrapper for single RAP-projected expression vectors
    subroutine relative_axes_expression_from_expression_vector_c(&
            vec,&
            n_axes,&
            contributions,&
            ierr&
            ) bind(C, name="relative_axes_expression_from_expression_vector_c")
        use tox_relative_axis_plane_tools, only: relative_axes_expression_from_expression_vector
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_axes
            !! Number of axes
        real(c_double), intent(in), dimension(n_axes), target :: vec
            !! RAP-projected and normalized expression vector
        real(c_double), intent(out), dimension(n_axes), target :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(c_int), intent(out), target :: ierr
            !! Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(vec)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(contributions)
        call  relative_axes_expression_from_expression_vector(&
            vec = vec,&
            n_axes = n_axes,&
            contributions = contributions,&
            ierr = ierr&
        )
    end subroutine relative_axes_expression_from_expression_vector_c

end module tox_relative_axis_plane_tools_c
#endif