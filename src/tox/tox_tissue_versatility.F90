#include <src/macros.h>

!> summary: Wrappers for [[tox_tissue_versatility_kernel(module)]]
!| Generated from the kernel; do not edit -- regenerate instead.
module tox_tissue_versatility
    use tox_tissue_versatility_kernel, only: compute_tissue_versatility_kernel
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use tox_errors, only: set_ok, is_err, validate_all_in_range_real, validate_dimension_size
    M_IMPLICIT_NONE
    private

    public :: compute_tissue_versatility

contains

    !> summary: Validates its inputs, then calls [[tox_tissue_versatility_kernel(module):compute_tissue_versatility_kernel]].
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    !|
    !| The selection-consistency checks (`n_selected_axes` as a dimension, and each selection count
    !| matching its claimed total) live here: they compare a `count(mask)` against a claimed size, which
    !| the generated wrapper's per-argument validators cannot express.
    subroutine compute_tissue_versatility(&
            n_axes,&
            n_vectors,&
            expression_vectors,&
            exp_vecs_selection_index,&
            n_selected_vectors,&
            axes_selection,&
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
            !! Number of selected expression vectors (count of .TRUE. in exp_vecs_selection_index)
        real(real64), dimension(n_axes, n_vectors), intent(in) :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical, dimension(n_vectors), intent(in) :: exp_vecs_selection_index
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical, dimension(n_axes), intent(in) :: axes_selection
            !! Logical array (n_axes), .TRUE. for axes to include in calculation
        integer(int32), intent(in) :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection)
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_versatilities
            !! Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_angles_deg
            !! Output, real array, length = n_selected_vectors, stores the calculated angles in degrees
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        call validate_dimension_size(n_axes, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_vectors, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_selected_vectors, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(expression_vectors, n_axes * n_vectors, ierr, arg_pos=3_int32)
        if (is_err(ierr)) return

        call compute_tissue_versatility_kernel(&
            n_axes = n_axes,&
            n_vectors = n_vectors,&
            expression_vectors = expression_vectors,&
            exp_vecs_selection_index = exp_vecs_selection_index,&
            n_selected_vectors = n_selected_vectors,&
            axes_selection = axes_selection,&
            n_selected_axes = n_selected_axes,&
            tissue_versatilities = tissue_versatilities,&
            tissue_angles_deg = tissue_angles_deg,&
            ierr = ierr&
        )
    end subroutine compute_tissue_versatility

end module tox_tissue_versatility
