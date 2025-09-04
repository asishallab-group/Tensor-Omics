!> Module for calculating normalized tissue (axis) versatility.
!| This module implements the angle-based metric for tissue versatility,
!| quantifying how uniformly a gene is expressed across selected axes (tissues).
module avmod
  use, intrinsic :: iso_fortran_env, only: real64, int32
  use tox_errors, only: ERR_EMPTY_INPUT, set_ok, set_err_once
  implicit none
contains

  !> Computes normalized tissue versatility for selected expression vectors.
  !| The metric is based on the angle between each gene expression vector and the space diagonal.
  !| Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis. 
  
  pure subroutine compute_tissue_versatility(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                             n_selected_vectors, axes_selection, n_selected_axes, &
                                             tissue_versatilities, tissue_angles_deg, ierr)
    !| Number of axes (tissues/dimensions)
    integer(int32), intent(in) :: n_axes
    !| Number of expression vectors (genes)
    integer(int32), intent(in) :: n_vectors
    !| Number of selected axes (count of .TRUE. in axes_selection)
    integer(int32), intent(in) :: n_selected_axes
    !| Number of selected expression vectors (count of .TRUE. in exp_vecs_selection_index)
    integer(int32), intent(in) :: n_selected_vectors
    !| 2D array (n_axes, n_vectors), each column is a gene expression vector
    real(real64), intent(in) :: expression_vectors(n_axes, n_vectors)
    !| Logical array (n_vectors), .TRUE. for vectors to process
    logical, intent(in) :: exp_vecs_selection_index(n_vectors)
    !| Logical array (n_axes), .TRUE. for axes to include in calculation
    logical, intent(in) :: axes_selection(n_axes)
    !| Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
    real(real64), intent(out) :: tissue_versatilities(n_selected_vectors)
    !| Output, real array, length = n_selected_vectors, stores the calculated angles in degrees
    real(real64), intent(out) :: tissue_angles_deg(n_selected_vectors)
    !| Error code: 0 = success, non-zero = error
    integer(int32), intent(out) :: ierr
    ! Local variables
    integer(int32) :: i_vec, i_axis, out_idx
    real(real64) :: norm_diag, dot_prod, norm_v, cos_phi, angle_rad, norm_factor
    real(real64), parameter :: rad2deg = 180.0_real64 / acos(-1.0_real64)

    ! Initialize error code
    call set_ok(ierr)

    ! Compute the norm of the space diagonal (only active axes)
    if (n_selected_axes <= 0) then
      ! Set error code to 202 (Empty input arrays)
      call set_err_once(ierr, ERR_EMPTY_INPUT)
      return
    end if
    norm_diag = sqrt(real(n_selected_axes, real64))
    ! Precompute normalization factor for tissue versatility
    norm_factor = 1.0_real64 - 1.0_real64 / norm_diag

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
        cos_phi = dot_prod / (sqrt(norm_v) * norm_diag)
      end if
      ! Clamp cos_phi for numerical safety
      cos_phi = max(-1.0_real64, min(1.0_real64, cos_phi))
      if (abs(cos_phi - 1.0_real64) < 1e-12_real64) cos_phi = 1.0_real64
      angle_rad = acos(cos_phi)
      tissue_versatilities(out_idx) = (1.0_real64 - cos_phi) / norm_factor
      tissue_angles_deg(out_idx) = angle_rad * rad2deg
      if (abs(tissue_angles_deg(out_idx)) < 1e-12_real64) tissue_angles_deg(out_idx) = 0.0_real64
    end do

  end subroutine compute_tissue_versatility

end module avmod 


!> R wrapper for compute_tissue_versatility.
!| Calls compute_tissue_versatility with standard Fortran types for R interface.
pure subroutine compute_tissue_versatility_r(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                             n_selected_vectors, axes_selection, n_selected_axes, &
                                             tissue_versatilities, tissue_angles_deg, ierr)
  use avmod
  !| Number of axes (tissues/dimensions)
  integer(int32), intent(in) :: n_axes
  !| Number of expression vectors (genes)
  integer(int32), intent(in) :: n_vectors
  !| Number of selected axes (count of .TRUE. in axes_selection)
  integer(int32), intent(in) :: n_selected_axes
  !| Number of selected expression vectors (count of .TRUE. in exp_vecs_selection_index)
  integer(int32), intent(in) :: n_selected_vectors
  !| 2D array (n_axes, n_vectors), each column is a gene expression vector
  real(real64), intent(in) :: expression_vectors(n_axes, n_vectors)
  !| Logical array (n_vectors), .TRUE. for vectors to process
  logical, intent(in) :: exp_vecs_selection_index(n_vectors)
  !| Logical array (n_axes), .TRUE. for axes to include in calculation
  logical, intent(in) :: axes_selection(n_axes)
  !| Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
  real(real64), intent(out) :: tissue_versatilities(n_selected_vectors)
  !| Output, real array, length = n_selected_vectors, stores the calculated angles in degrees
  real(real64), intent(out) :: tissue_angles_deg(n_selected_vectors)
  !| Error code: 0 = success, non-zero = error
  integer(int32), intent(out) :: ierr
  call compute_tissue_versatility(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, n_selected_vectors, axes_selection, &
                                  n_selected_axes, tissue_versatilities, tissue_angles_deg, ierr)
end subroutine compute_tissue_versatility_r


!> C wrapper for compute_tissue_versatility.
!| Exposes compute_tissue_versatility to C via iso_c_binding types with explicit dimensions.
pure subroutine compute_tissue_versatility_c(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index, &
                                             n_selected_vectors, axes_selection, n_selected_axes, &
                                             tissue_versatilities, tissue_angles_deg, ierr) bind(C, name="compute_tissue_versatility_c")
  use iso_c_binding, only : c_int, c_double
  use avmod
  !| Number of axes (tissues/dimensions)
  integer(c_int), intent(in), value :: n_axes
  !| Number of expression vectors (genes)
  integer(c_int), intent(in), value :: n_vectors
  !| 2D array (n_axes, n_vectors), each column is a gene expression vector (column-major)
  real(c_double), intent(in), target :: expression_vectors(n_axes, n_vectors)
  !| Integer array (n_vectors), 0/1 values. 0=not selected, 1=selected. Interpreted as logical internally.
  integer(c_int), intent(in), target :: exp_vecs_selection_index(n_vectors)
  !| Number of selected expression vectors (count of 1s in exp_vecs_selection_index)
  integer(c_int), intent(in), value :: n_selected_vectors
  !| Integer array (n_axes), 0/1 values. 0=not selected, 1=selected. Interpreted as logical internally.
  integer(c_int), intent(in), target :: axes_selection(n_axes)
  !| Number of selected axes (count of 1s in axes_selection)
  integer(c_int), intent(in), value :: n_selected_axes
  !| Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities for selected vectors
  real(c_double), intent(out), target :: tissue_versatilities(n_selected_vectors)
  !| Output, real array, length = n_selected_vectors, stores the calculated angles in degrees for selected vectors
  real(c_double), intent(out), target :: tissue_angles_deg(n_selected_vectors)
  !| Error code: 0 = success, non-zero = error  
  integer(c_int), intent(out) :: ierr

  call compute_tissue_versatility(n_axes, n_vectors, expression_vectors, exp_vecs_selection_index /= 0, n_selected_vectors, &
                                  axes_selection /= 0, n_selected_axes, tissue_versatilities, tissue_angles_deg, ierr)
end subroutine compute_tissue_versatility_c

