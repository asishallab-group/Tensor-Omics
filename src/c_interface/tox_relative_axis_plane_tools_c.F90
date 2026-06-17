#ifndef NO_C_INTERFACE
#include <src/macros.h>

!> summary: Module for C-wrappers for [[tox_relative_axis_plane_tools(module)]]
!| Module for tools related to relative axis planes (RAPs), i.e. planes in higher-dimensional gene expression space
!| category: C-interface
!| Compute the signed clock hand angle between two RAP-projected and normalized vectors.
!| Calculates the signed rotation angle between two normalized vectors in RAP space.
!| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
!| category: C-interface
!| Compute signed rotation angles between RAP-projected and normalized vector pairs.
!| Takes separate arrays of RAP-projected and normalized vectors (e.g. expression centroids and paralogs) and computes the signed rotation angle between corresponding pairs.
!| This measures both magnitude and directionality of angular separation in RAP space.
!| category: C-interface
!| Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
!| Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
!| category: C-interface
!| Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
!| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
!| category: C-interface
!| Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
!| Wrapper for single RAP-projected expression vectors
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
    subroutine omics_vector_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr) bind(C, name="omics_vector_RAP_projection_c")
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
        call omics_vector_RAP_projection(vecs = vecs, n_axes = n_axes, n_vecs = n_vecs, vecs_selection_mask = vecs_selection_mask_f, n_selected_vecs = n_selected_vecs, axes_selection_mask = axes_selection_mask_f, n_selected_axes = n_selected_axes, projections = projections, ierr = ierr)
    end subroutine omics_vector_RAP_projection_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):omics_field_RAP_projection(subroutine)]]
    !| Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    subroutine omics_field_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr) bind(C, name="omics_field_RAP_projection_c")
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
        call omics_field_RAP_projection(vecs = vecs, n_axes = n_axes, n_vecs = n_vecs, vecs_selection_mask = vecs_selection_mask_f, n_selected_vecs = n_selected_vecs, axes_selection_mask = axes_selection_mask_f, n_selected_axes = n_selected_axes, projections = projections, ierr = ierr)
    end subroutine omics_field_RAP_projection_c

    !> summary: C-wrapper for [[tox_relative_axis_plane_tools(module):project_selected_vecs_onto_rap(subroutine)]]
    !| Projects selected vectors onto its RAP
    subroutine project_selected_vecs_onto_rap_c(selected_vecs, n_selected_axes, n_selected_vecs, ierr) bind(C, name="project_selected_vecs_onto_rap_c")
        use tox_relative_axis_plane_tools, only: project_selected_vecs_onto_rap
        use tox_relative_axis_plane_tools
        integer(c_int), intent(in), target :: n_selected_axes
            !! number of selected axes
        integer(c_int), intent(in), target :: n_selected_vecs
            !! number of selected vectors per axis
        real(c_double), intent(inout), dimension(n_selected_axes, n_selected_vecs), target :: selected_vecs
            !! matrix with vectors for selected axes
        integer(c_int), intent(out), target :: ierr
            !!  Error code
        M_CHECK_IERR_NON_NULL
        M_CHECK_NON_NULL(selected_vecs)
        M_CHECK_NON_NULL(n_selected_axes)
        M_CHECK_NON_NULL(n_selected_vecs)
        call project_selected_vecs_onto_rap(selected_vecs = selected_vecs, n_selected_axes = n_selected_axes, n_selected_vecs = n_selected_vecs)
    end subroutine project_selected_vecs_onto_rap_c

end module tox_relative_axis_plane_tools_c
#endif