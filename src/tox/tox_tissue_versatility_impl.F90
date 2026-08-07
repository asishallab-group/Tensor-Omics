#include <src/macros.h>

!> Kernel for calculating normalized tissue (axis) versatility.
!| This module implements the angle-based metric for tissue versatility, quantifying how uniformly a
!| gene is expressed across selected axes (tissues). The generator turns
!| `compute_tissue_versatility_impl` into the validating wrapper `compute_tissue_versatility` in
!| module `tox_tissue_versatility`.
module tox_tissue_versatility_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use f42_utils, only: clamp, degrees
    M_IMPLICIT_NONE
contains

    !> summary: Computes normalized tissue versatility for selected expression vectors.
    !| AUTHOR_VIVIAN_BASS
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    !|
    !| The masks follow the `n_selected_` convention, so the generated wrapper validates that each
    !| selection count matches its mask; `n_selected_axes` (not an array extent) carries its own floor.
    pure subroutine compute_tissue_versatility_impl(n_axes, n_vectors, expression_vectors, vectors_selection_mask, &
                                               n_selected_vectors, axes_selection_mask, n_selected_axes, &
                                               tissue_versatilities, tissue_angles_deg)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions)
        integer(int32), intent(in) :: n_vectors
            !! Number of expression vectors (genes)
        integer(int32), intent(in) :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection_mask)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_selected_vectors
            !! Number of selected expression vectors (count of .TRUE. in vectors_selection_mask)
        real(real64), dimension(n_axes, n_vectors), intent(in) :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical, dimension(n_vectors), intent(in) :: vectors_selection_mask
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! Logical array (n_axes), .TRUE. for axes to include in calculation
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_versatilities
            !! Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
        real(real64), dimension(n_selected_vectors), intent(out) :: tissue_angles_deg
            !! Output, real array, length = n_selected_vectors, stores the calculated angles in degrees

        ! Local variables
        integer(int32) :: i_vec, i_axis, out_idx
        real(real64) :: norm_diag, dot_prod, norm_v, cos_phi, angle_rad, norm_factor

        ! Handle edge case: when only one axis is selected, tissue versatility is always 0
        if (n_selected_axes == 1) then
            tissue_versatilities = 0.0_real64
            tissue_angles_deg = 0.0_real64
            return
        end if

        norm_diag = sqrt(real(n_selected_axes, real64))
        ! cos_phi against the space diagonal ranges from 1 (uniform expression across all selected
        ! axes) down to 1/norm_diag (all expression concentrated in a single axis), so (1 - cos_phi)
        ! ranges over [0, 1 - 1/norm_diag]. Dividing by that upper bound rescales versatility to [0, 1].
        norm_factor = 1.0_real64 - 1.0_real64/norm_diag

        ! Loop over selected expression vectors
        ! Note: If the expression vector is zero in all selected axes, tissue versatility (TV) is set to 1 (maximum specificity) and the angle is set to acos(0) = 90 degrees.
        !TODO optimize: this per-gene loop is data-parallel (each i_vec is independent) but is forced sequential here only
        ! because out_idx is an accumulating compaction counter. Since n_selected_vectors is already known, this
        ! could be parallelized with `do concurrent` by precomputing an exclusive prefix sum of
        ! vectors_selection_mask to get each vector's output slot directly, consistent with how other modules in
        ! this codebase parallelize per-gene work.
        out_idx = 0
        do i_vec = 1, n_vectors
            if (.not. vectors_selection_mask(i_vec)) cycle  ! Skip if not selected

            ! Compute dot product and norm for the vector in active axes
            dot_prod = 0.0_real64
            norm_v = 0.0_real64

            do i_axis = 1, n_axes
                if (axes_selection_mask(i_axis)) then
                    dot_prod = dot_prod + expression_vectors(i_axis, i_vec)
                    norm_v = norm_v + expression_vectors(i_axis, i_vec)**2
                end if
            end do

            out_idx = out_idx + 1
            ! If the vector is zero or numerically negligible, set TV = 1 (maximum specificity)
            ! Use sqrt(epsilon) for extra-robust threshold to avoid numerical instability in cos_phi calculation
            if (norm_v <= sqrt(epsilon(1.0_real64))) then
                tissue_versatilities(out_idx) = 1.0_real64
                tissue_angles_deg(out_idx) = 90.0_real64
                cycle
            else
                cos_phi = dot_prod/(sqrt(norm_v)*norm_diag)
            end if
            ! Clamp cos_phi for numerical safety
            cos_phi = clamp(cos_phi, min_val=-1.0_real64, max_val=1.0_real64)
            angle_rad = acos(cos_phi)
            tissue_versatilities(out_idx) = (1.0_real64 - cos_phi)/norm_factor
            tissue_angles_deg(out_idx) = degrees(angle_rad)
        end do

    end subroutine compute_tissue_versatility_impl

end module tox_tissue_versatility_impl
