#include <src/macros.h>

!> Module for calculating normalized tissue (axis) versatility.
!| This module implements the angle-based metric for tissue versatility,
!| quantifying how uniformly a gene is expressed across selected axes (tissues).
module tox_tissue_versatility
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_INVALID_INPUT, set_ok, set_err_once, validate_dimension_size, validate_all_in_range_real, is_err
    use f42_utils, only: operator(.isclose.), clamp, degrees
    implicit none
contains

    !> AUTHOR_VIVIAN_BASS
    !| Computes normalized tissue versatility for selected expression vectors.
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    pure subroutine compute_tissue_versatility(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                               n_selected_vectors, axes_selection, n_selected_axes, &
                                               tissue_versatilities, tissue_angles_deg, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions)
        integer(int32), intent(in) :: n_vectors
            !! Number of expression vectors (genes)
        integer(int32), intent(in) :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection)
        integer(int32), intent(in) :: n_selected_vectors
            !! Number of selected expression vectors (count of .TRUE. in exp_vecs_selection_index)
        real(real64), dimension(n_axes, n_vectors), intent(in) :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical, dimension(n_vectors), intent(in) :: exp_vecs_selection_index
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical, dimension(n_axes), intent(in) :: axes_selection
            !! Logical array (n_axes), .TRUE. for axes to include in calculation
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
        call validate_dimension_size(n_selected_axes, ierr, arg_pos=7_int32)
        call validate_all_in_range_real(expression_vectors, size(expression_vectors, kind=int32), ierr, arg_pos=3_int32)
        if (count(exp_vecs_selection_index, kind=int32) > n_selected_vectors) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
        if (count(axes_selection, kind=int32) > n_selected_axes) call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos=6_int32)

        if (is_err(ierr)) return

        call compute_tissue_versatility_helper(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                               n_selected_vectors, axes_selection, n_selected_axes, &
                                               tissue_versatilities, tissue_angles_deg)
    end subroutine compute_tissue_versatility

    !> AUTHOR_VIVIAN_BASS
    !| (no input validation) Computes normalized tissue versatility for selected expression vectors.
    !| The metric is based on the angle between each gene expression vector and the space diagonal.
    !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.
    pure subroutine compute_tissue_versatility_helper(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                               n_selected_vectors, axes_selection, n_selected_axes, &
                                               tissue_versatilities, tissue_angles_deg)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (tissues/dimensions)
        integer(int32), intent(in) :: n_vectors
            !! Number of expression vectors (genes)
        integer(int32), intent(in) :: n_selected_axes
            !! Number of selected axes (count of .TRUE. in axes_selection)
        integer(int32), intent(in) :: n_selected_vectors
            !! Number of selected expression vectors (count of .TRUE. in exp_vecs_selection_index)
        real(real64), dimension(n_axes, n_vectors), intent(in) :: expression_vectors
            !! 2D array (n_axes, n_vectors), each column is a gene expression vector
        logical, dimension(n_vectors), intent(in) :: exp_vecs_selection_index
            !! Logical array (n_vectors), .TRUE. for vectors to process
        logical, dimension(n_axes), intent(in) :: axes_selection
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
        ! Precompute normalization factor for tissue versatility
        norm_factor = 1.0_real64 - 1.0_real64/norm_diag

        ! Loop over selected expression vectors
        ! Note: If the expression vector is zero in all selected axes, tissue versatility (TV) is set to 1 (maximum specificity) and the angle is set to acos(0) = 90 degrees.
        out_idx = 0
        do i_vec = 1, n_vectors
            if (.not. exp_vecs_selection_index(i_vec)) cycle  ! Skip if not selected

            ! Compute dot product and norm for the vector in active axes
            dot_prod = 0.0_real64
            norm_v = 0.0_real64

            do i_axis = 1, n_axes
                if (axes_selection(i_axis)) then
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

    end subroutine compute_tissue_versatility_helper

end module tox_tissue_versatility

!> C wrapper for compute_tissue_versatility.
!| Exposes compute_tissue_versatility to C via iso_c_binding types with explicit dimensions.
pure subroutine compute_tissue_versatility_c(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                             n_selected_vectors, axes_selection, n_selected_axes, &
                                             tissue_versatilities, tissue_angles_deg, ierr) bind(C, name="compute_tissue_versatility_c")
    use, intrinsic :: iso_c_binding, only: c_int, c_double
    use tox_tissue_versatility, only: compute_tissue_versatility
    use tox_conversions, only: c_int_as_logical
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_axes
        !! Number of axes (tissues/dimensions)
    integer(c_int), intent(in), target :: n_vectors
        !! Number of expression vectors (genes)
    real(c_double), dimension(n_axes, n_vectors), intent(in), target :: expression_vectors
        !! 2D array (n_axes, n_vectors), each column is a gene expression vector (column-major)
    integer(c_int), dimension(n_vectors), intent(in), target :: exp_vecs_selection_index
        !! Integer array (n_vectors), 0/1 values. 0=not selected, 1=selected. Interpreted as logical internally.
    integer(c_int), intent(in), target :: n_selected_vectors
        !! Number of selected expression vectors (count of 1s in exp_vecs_selection_index)
    integer(c_int), dimension(n_axes), intent(in), target :: axes_selection
        !! Integer array (n_axes), 0/1 values. 0=not selected, 1=selected. Interpreted as logical internally.
    integer(c_int), intent(in), target :: n_selected_axes
        !! Number of selected axes (count of 1s in axes_selection)
    real(c_double), dimension(n_selected_vectors), intent(out), target :: tissue_versatilities
        !! Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities for selected vectors
    real(c_double), dimension(n_selected_vectors), intent(out), target :: tissue_angles_deg
        !! Output, real array, length = n_selected_vectors, stores the calculated angles in degrees for selected vectors
    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical, dimension(:), allocatable :: exp_vecs_selection_index_f, axes_selection_f

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_vectors)
    M_CHECK_NON_NULL(expression_vectors)
    M_CHECK_NON_NULL(exp_vecs_selection_index)
    M_CHECK_NON_NULL(n_selected_vectors)
    M_CHECK_NON_NULL(axes_selection)
    M_CHECK_NON_NULL(n_selected_axes)
    M_CHECK_NON_NULL(tissue_versatilities)
    M_CHECK_NON_NULL(tissue_angles_deg)

    M_ALLOCATE(exp_vecs_selection_index_f(n_vectors))
    call c_int_as_logical(exp_vecs_selection_index, exp_vecs_selection_index_f)
    M_ALLOCATE(axes_selection_f(n_axes))
    call c_int_as_logical(axes_selection, axes_selection_f)

    call compute_tissue_versatility(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index_f, n_selected_vectors, &
                                    axes_selection_f, n_selected_axes, tissue_versatilities, tissue_angles_deg, ierr)
end subroutine compute_tissue_versatility_c

