#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[tox_tissue_versatility(module)]]
!| Normalized tissue (axis) versatility: how uniformly a gene is expressed across tissues.
!|
!| An angle-based metric. A gene expressed equally across every selected axis points along the
!| diagonal of that subspace and scores maximally versatile; one confined to a single tissue
!| points along that axis and scores minimally. Normalized, so scores over different numbers of
!| selected axes are comparable.
module tox_tissue_versatility_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: compute_tissue_versatility_c

contains

    !> summary: C-wrapper for [[tox_tissue_versatility(module):compute_tissue_versatility(subroutine)]]
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    !|
    !| The masks follow the `n_selected_` convention, so the generated wrapper validates that each
    !| selection count matches its mask; `n_selected_axes` (not an array extent) carries its own floor.
    subroutine compute_tissue_versatility_c(&
            n_axes,&
            n_vectors,&
            expression_vectors,&
            vectors_selection_mask,&
            n_selected_vectors,&
            axes_selection_mask,&
            n_selected_axes,&
            tissue_versatilities,&
            tissue_angles_deg,&
            ierr&
        ) bind(C, name="compute_tissue_versatility_c")
        use tox_tissue_versatility, only: compute_tissue_versatility

        integer(c_int), intent(in), target :: n_axes
            !! Number of axes (tissues/dimensions)
        integer(c_int), intent(in), target :: n_vectors
            !! Number of expression vectors (genes)
        integer(c_int), intent(in), target :: n_selected_vectors
            !! Number of selected expression vectors (count of .TRUE. in vectors_selection_mask)
        real(c_double), dimension(n_axes, n_vectors), intent(in), target :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical(c_bool), dimension(n_vectors), intent(in), target :: vectors_selection_mask
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical(c_bool), dimension(n_axes), intent(in), target :: axes_selection_mask
            !! Logical array (n_axes), .TRUE. for axes to include in calculation
        integer(c_int), intent(in), target :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection_mask)
            !! The minimum valid value is `1_int32`.
        real(c_double), dimension(n_selected_vectors), intent(out), target :: tissue_versatilities
            !! Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
        real(c_double), dimension(n_selected_vectors), intent(out), target :: tissue_angles_deg
            !! Output, real array, length = n_selected_vectors, stores the calculated angles in degrees
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.
        logical, dimension(n_vectors) :: vectors_selection_mask_f
        logical, dimension(n_axes) :: axes_selection_mask_f

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_axes)
        M_CHECK_NON_NULL(n_vectors)
        M_CHECK_NON_NULL(n_selected_vectors)
        M_CHECK_NON_NULL(n_selected_axes)
        M_CHECK_ARRAY_NON_NULL(expression_vectors, n_axes * n_vectors)
        M_CHECK_ARRAY_NON_NULL(vectors_selection_mask, n_vectors)
        M_CHECK_ARRAY_NON_NULL(axes_selection_mask, n_axes)
        M_CHECK_ARRAY_NON_NULL(tissue_versatilities, n_selected_vectors)
        M_CHECK_ARRAY_NON_NULL(tissue_angles_deg, n_selected_vectors)

        vectors_selection_mask_f = vectors_selection_mask
        axes_selection_mask_f = axes_selection_mask

        call compute_tissue_versatility(&
            n_axes = n_axes,&
            n_vectors = n_vectors,&
            expression_vectors = expression_vectors,&
            vectors_selection_mask = vectors_selection_mask_f,&
            n_selected_vectors = n_selected_vectors,&
            axes_selection_mask = axes_selection_mask_f,&
            n_selected_axes = n_selected_axes,&
            tissue_versatilities = tissue_versatilities,&
            tissue_angles_deg = tissue_angles_deg,&
            ierr = ierr&
        )
    end subroutine compute_tissue_versatility_c

end module tox_tissue_versatility_c
#endif
