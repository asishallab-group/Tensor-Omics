#include <src/macros.h>

!> summary: Wrappers for [[tox_tissue_versatility_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_tissue_versatility
    use tox_tissue_versatility_kernel, only: compute_tissue_versatility_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, ERR_INVALID_INPUT, set_err_once
    use tox_errors, only: validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: compute_tissue_versatility

contains

    !> summary: Validates its inputs, then calls [[tox_tissue_versatility_kernel(module):compute_tissue_versatility_kernel]].
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    !|
    !| The masks follow the `n_selected_` convention, so the generated wrapper validates that each
    !| selection count matches its mask; `n_selected_axes` (not an array extent) carries its own floor.
    subroutine compute_tissue_versatility(&
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
        )
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions)
        integer(int32), intent(in) :: n_vectors
            !! Number of expression vectors (genes)
        integer(int32), intent(in) :: n_selected_vectors
            !! Number of selected expression vectors (count of .TRUE. in vectors_selection_mask)
        real(real64), dimension(n_axes, n_vectors), intent(in) :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical, dimension(n_vectors), intent(in) :: vectors_selection_mask
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! Logical array (n_axes), .TRUE. for axes to include in calculation
        integer(int32), intent(in) :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection_mask)
            !! The minimum valid value is `1_int32`.
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_versatilities
            !! Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_angles_deg
            !! Output, real array, length = n_selected_vectors, stores the calculated angles in degrees
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_selected_vectors, ierr, arg_pos=5_int32)
        call validate_in_range_int(n_selected_axes, ierr, arg_pos=7_int32, min=1_int32)
        call validate_all_in_range_real(expression_vectors, n_axes * n_vectors, ierr, arg_pos=3_int32)
        if (count(vectors_selection_mask, kind=int32) /= n_selected_vectors) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (count(axes_selection_mask, kind=int32) /= n_selected_axes) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)
        if (is_err(ierr)) return

        call compute_tissue_versatility_kernel(&
            n_axes = n_axes,&
            n_vectors = n_vectors,&
            expression_vectors = expression_vectors,&
            vectors_selection_mask = vectors_selection_mask,&
            n_selected_vectors = n_selected_vectors,&
            axes_selection_mask = axes_selection_mask,&
            n_selected_axes = n_selected_axes,&
            tissue_versatilities = tissue_versatilities,&
            tissue_angles_deg = tissue_angles_deg&
        )
    end subroutine compute_tissue_versatility

end module tox_tissue_versatility
