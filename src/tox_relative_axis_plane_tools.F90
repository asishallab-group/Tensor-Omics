!> Module for tools related to relative axis planes (RAPs), i.e. planes in higher-dimensional gene expression space
module relative_axis_plane_tools
   use, intrinsic :: iso_fortran_env, only: real64, int32
   use, intrinsic :: ieee_arithmetic
   use tox_errors, only: ERR_OK, ERR_INVALID_INPUT, set_ok, set_err_once, is_ok
   implicit none

contains

   !> Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
   pure subroutine omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
      real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
         !! matrix with expression vectors
      integer(int32), intent(in) :: n_axes
         !! number of axes
      integer(int32), intent(in) :: n_vecs
         !! number of vectors per axis
      logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
         !! `.true.` for vectors where projection is to be computed
      integer(int32), intent(in) :: n_selected_vecs
         !! count of `.true.` values in `vecs_selection_mask`
      logical, dimension(n_axes), intent(in) :: axes_selection_mask
         !! `.true.` for axes to be included in RAP
      integer(int32), intent(in) :: n_selected_axes
         !! count of `.true.` values in `axes_selection_mask`
      real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
         !! projected vectors
      integer(int32), intent(out) :: ierr
         !! Error code

      integer(int32) :: i_vec, i_axis, i_vec_proj, i_axis_proj
      ! Error handling
      call set_ok(ierr)
      if (n_selected_vecs < 1 .or. n_selected_axes < 1) then
         projections = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if
      if (n_selected_vecs > n_vecs .or. n_selected_axes > n_axes) then
         projections = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if
      i_vec_proj = 1
      do i_vec = 1, n_vecs
         if (vecs_selection_mask(i_vec)) then
            i_axis_proj = 1
            do i_axis = 1, n_axes
               if (axes_selection_mask(i_axis)) then
                  projections(i_axis_proj, i_vec_proj) = vecs(i_axis, i_vec)

                  i_axis_proj = i_axis_proj + 1
               end if
            end do

            i_vec_proj = i_vec_proj + 1
         end if
      end do

      call project_selected_vecs_onto_rap(projections, n_selected_axes, n_selected_vecs)
   end subroutine omics_vector_RAP_projection

   !> Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
   pure subroutine omics_field_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
      real(real64), dimension(2 * n_axes, n_vecs), intent(in) :: vecs
         !! matrix with vector fields, first n rows mean vector origin, last n rows vector targets
      integer(int32), intent(in) :: n_axes
         !! number of axes
      integer(int32), intent(in) :: n_vecs
         !! number of vectors per axis
      logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
         !! `.true.` for vectors where projection is to be computed
      integer(int32), intent(in) :: n_selected_vecs
         !! count of `.true.` values in `vecs_selection_mask`
      logical, dimension(n_axes), intent(in) :: axes_selection_mask
         !! `.true.` for axes to be included in RAP
      integer(int32), intent(in) :: n_selected_axes
         !! count of `.true.` values in `axes_selection_mask`
      real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
         !! projected vectors
      integer(int32), intent(out) :: ierr
         !! Error code

      integer(int32) :: i_vec, i_axis, i_vec_proj, i_axis_proj
      ! Error handling
      call set_ok(ierr)
      if (n_selected_vecs < 1 .or. n_selected_axes < 1) then
         projections = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if
      if (n_selected_vecs > n_vecs .or. n_selected_axes > n_axes) then
         projections = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if
      i_vec_proj = 1
      do i_vec = 1, n_vecs
         if (vecs_selection_mask(i_vec)) then
            i_axis_proj = 1
            do i_axis = 1, n_axes
               if (axes_selection_mask(i_axis)) then
                  ! compute shift vector as difference between origin and target
                  projections(i_axis_proj, i_vec_proj) = vecs(i_axis, i_vec) - vecs(i_axis + n_axes, i_vec)

                  i_axis_proj = i_axis_proj + 1
               end if
            end do

            i_vec_proj = i_vec_proj + 1
         end if
      end do

      call project_selected_vecs_onto_rap(projections, n_selected_axes, n_selected_vecs)
   end subroutine omics_field_RAP_projection

   !> Projects selected vectors onto its RAP
   pure subroutine project_selected_vecs_onto_rap(selected_vecs, n_selected_axes, n_selected_vecs)
      real(real64), dimension(n_selected_axes, n_selected_vecs), intent(inout) :: selected_vecs
         !! matrix with vectors for selected axes
      integer(int32), intent(in) :: n_selected_axes
         !! number of selected axes
      integer(int32), intent(in) :: n_selected_vecs
         !! number of selected vectors per axis

      ! project selected vectors onto RAP
      integer(int32) :: i_vec, i_axis
      real(real64) :: diagonal_component

      do i_vec = 1, n_selected_vecs

         ! calculate diagonal component to be subtracted from vectors for projection
         diagonal_component = 0.0_real64
         do i_axis = 1, n_selected_axes
            diagonal_component = diagonal_component + selected_vecs(i_axis, i_vec)
         end do
         diagonal_component = diagonal_component / n_selected_axes

         ! transform vector to its projection
         do i_axis = 1, n_selected_axes
            selected_vecs(i_axis, i_vec) = selected_vecs(i_axis, i_vec) - diagonal_component
         end do
      end do
   end subroutine project_selected_vecs_onto_rap

   !> Compute the signed clock hand angle between two RAP-projected and normalized vectors.
   !! Calculates the signed rotation angle between two normalized vectors in RAP space.
   !! For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
   pure subroutine clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
      real(real64), dimension(n_dims), intent(in) :: v1
         !! First normalized vector in RAP space
      real(real64), dimension(n_dims), intent(in) :: v2
         !! Second normalized vector in RAP space  
      integer(int32), intent(in) :: n_dims
         !! Dimension of both vectors
      real(real64), intent(out) :: signed_angle
         !! Signed angle between vectors in radians [-π, π]
      integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
         !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
      integer(int32), intent(out) :: ierr
         !! Error code

      real(real64) :: dot_product, unsigned_angle, orientation_sign
      integer(int32) :: i

      ! Error handling
      call set_ok(ierr)
      if (n_dims < 1) then
         signed_angle = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if

      ! Validate that indices are in range 
      if (any(selected_axes_for_signed < 1) .or. any(selected_axes_for_signed > n_dims)) then
         signed_angle = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if
      ! Validate that indices are unique if n_dims > 3
      if (n_dims > 3) then
         if (selected_axes_for_signed(1) == selected_axes_for_signed(2) .or. &
             selected_axes_for_signed(1) == selected_axes_for_signed(3) .or. &
             selected_axes_for_signed(2) == selected_axes_for_signed(3)) then
            signed_angle = 0.0_real64
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
         end if
      end if

      ! Calculate dot product of normalized vectors
      dot_product = 0.0_real64
      do i = 1, n_dims
         dot_product = dot_product + v1(i) * v2(i)
      end do

      ! Clamp dot product to [-1, 1] to handle numerical precision issues
      dot_product = max(-1.0_real64, min(1.0_real64, dot_product))

      ! Calculate unsigned angle using arccos
      unsigned_angle = acos(dot_product)

      ! Calculate orientation sign for directionality
      select case (n_dims)
         case (2)
            ! For 2D: use determinant directly
            orientation_sign = sign(1.0_real64, v1(1) * v2(2) - v1(2) * v2(1))
         case (3)
            ! For 3D, use [1,2,3] directly
            orientation_sign = cross_product_orientation_sign(v1, v2, n_dims, [1,2,3])
         case (4:)
            ! For >3D, use selected_axes_for_signed
            orientation_sign = cross_product_orientation_sign(v1, v2, n_dims, selected_axes_for_signed)
         case default
            orientation_sign = 1.0_real64
      end select

      ! Apply sign to unsigned angle
      signed_angle = orientation_sign * unsigned_angle
   end subroutine clock_hand_angle_between_vectors

   !> Compute signed rotation angles between RAP-projected and normalized vector pairs.
   !! Takes separate arrays of RAP-projected and normalized vectors (e.g. expression centroids and paralogs) and computes the signed rotation angle between corresponding pairs.
   !! This measures both magnitude and directionality of angular separation in RAP space.
   pure subroutine clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, &
                                                      vecs_selection_mask, &
                                                      n_selected_vecs, selected_axes_for_signed, &
                                                      signed_angles, ierr)
      real(real64), dimension(n_dims, n_vecs), intent(in) :: origins
         !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
      real(real64), dimension(n_dims, n_vecs), intent(in) :: targets
         !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
      integer(int32), intent(in) :: n_dims
         !! Dimension of each vector in RAP space
      integer(int32), intent(in) :: n_vecs
         !! Number of vector pairs
      logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
         !! .true. for vector pairs where angle should be computed
      integer(int32), intent(in) :: n_selected_vecs
         !! Count of .true. values in vecs_selection_mask
      integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
         !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
      real(real64), dimension(n_selected_vecs), intent(out) :: signed_angles
         !! Signed rotation angles between vector pairs in radians [-π, π]
      integer(int32), intent(out) :: ierr
         !! Error code

      integer(int32) :: i_vec, i_result

      ! Error handling
      call set_ok(ierr)
      if (n_selected_vecs < 1 .or. n_dims < 1) then
         signed_angles = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if

      i_result = 1
      do i_vec = 1, n_vecs
         if (vecs_selection_mask(i_vec)) then
            call clock_hand_angle_between_vectors(origins(:, i_vec), targets(:, i_vec), n_dims, &
                                                 signed_angles(i_result), selected_axes_for_signed, ierr)
            if (.not. is_ok(ierr)) return
            i_result = i_result + 1
         end if
      end do
   end subroutine clock_hand_angles_for_shift_vectors

   !> Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
   !! Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
   pure subroutine compute_relative_axis_contributions(vec, n_axes, contributions, ierr)
      use, intrinsic :: ieee_arithmetic
      real(real64), dimension(n_axes), intent(in) :: vec
         !! RAP-projected and normalized vector (expression or shift)
      integer(int32), intent(in) :: n_axes
         !! Number of axes (length of vec and contributions)
      real(real64), dimension(n_axes), intent(out) :: contributions
         !! Fractional contribution of each axis (output), values in [0,1], sum to 1
      integer(int32), intent(out) :: ierr
         !! Error code

      real(real64) :: total_abs
      integer(int32) :: i_axis

      ! Error handling
      call set_ok(ierr)
      if (n_axes < 1) then
         contributions = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if

      total_abs = 0.0_real64
      do i_axis = 1, n_axes
         if (ieee_is_nan(vec(i_axis)) .or. .not. ieee_is_finite(vec(i_axis))) then
            contributions = 0.0_real64
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
         end if
         total_abs = total_abs + abs(vec(i_axis))
      end do

      if (total_abs < 1.0e-8_real64) then
         contributions = 0.0_real64
         call set_err_once(ierr, ERR_INVALID_INPUT)
         return
      end if

      do i_axis = 1, n_axes
         contributions(i_axis) = abs(vec(i_axis)) / total_abs
      end do

   end subroutine compute_relative_axis_contributions

   !> Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
   !! Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
   pure subroutine relative_axes_changes_from_shift_vector(vec, n_axes, contributions, ierr)
      real(real64), dimension(n_axes), intent(in) :: vec
         !! RAP-projected and normalized shift vector
      integer(int32), intent(in) :: n_axes
         !! Number of axes
      real(real64), dimension(n_axes), intent(out) :: contributions
         !! Fractional contribution of each axis (output), values in [0,1], sum to 1
      integer(int32), intent(out) :: ierr
         !! Error code

      call compute_relative_axis_contributions(vec, n_axes, contributions, ierr)
   end subroutine relative_axes_changes_from_shift_vector

   !> Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
   !! Wrapper for single RAP-projected expression vectors
   pure subroutine relative_axes_expression_from_expression_vector(vec, n_axes, contributions, ierr)
      real(real64), dimension(n_axes), intent(in) :: vec
         !! RAP-projected and normalized expression vector
      integer(int32), intent(in) :: n_axes
         !! Number of axes
      real(real64), dimension(n_axes), intent(out) :: contributions
         !! Fractional contribution of each axis (output), values in [0,1], sum to 1
      integer(int32), intent(out) :: ierr
         !! Error code

      call compute_relative_axis_contributions(vec, n_axes, contributions, ierr)
   end subroutine relative_axes_expression_from_expression_vector

   !> Compute orientation sign from cross product of two vectors, using selected axes 
   pure function cross_product_orientation_sign(a, b, n_dims, selected_axes) result(orientation_sign)
      real(real64), intent(in) :: a(n_dims), b(n_dims)
      integer(int32), intent(in) :: n_dims
      integer(int32), intent(in) :: selected_axes(3)
      real(real64) :: orientation_sign
      real(real64) :: cross1, cross2, cross3, dotprod
      cross1 = a(selected_axes(2))*b(selected_axes(3)) - a(selected_axes(3))*b(selected_axes(2))
      cross2 = a(selected_axes(3))*b(selected_axes(1)) - a(selected_axes(1))*b(selected_axes(3))
      cross3 = a(selected_axes(1))*b(selected_axes(2)) - a(selected_axes(2))*b(selected_axes(1))
      dotprod = cross1*a(selected_axes(1)) + cross2*a(selected_axes(2)) + cross3*a(selected_axes(3))
      orientation_sign = sign(1.0_real64, dotprod)
   end function cross_product_orientation_sign

end module relative_axis_plane_tools

! Updated wrappers to pass and return ierr
subroutine relative_axes_changes_from_shift_vector_r(vec, n_axes, contributions, ierr)
   use relative_axis_plane_tools, only: relative_axes_changes_from_shift_vector
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_axes), intent(in) :: vec
      !! RAP-projected and normalized shift vector
   integer(int32), intent(in) :: n_axes
      !! Number of axes
   real(real64), dimension(n_axes), intent(out) :: contributions
      !! Relative axis contributions (output), values in [0,1], sum to 1
   integer(int32), intent(out) :: ierr
      !! Error code

   call relative_axes_changes_from_shift_vector(vec, n_axes, contributions, ierr)
end subroutine relative_axes_changes_from_shift_vector_r

subroutine relative_axes_changes_from_shift_vector_c(vec, n_axes, contributions, ierr) bind(C, name="relative_axes_changes_from_shift_vector_c")
   use iso_c_binding, only: c_double, c_int
   use relative_axis_plane_tools, only: relative_axes_changes_from_shift_vector
   implicit none

   real(c_double), dimension(n_axes), intent(in) :: vec
      !! RAP-projected and normalized shift vector
   integer(c_int), intent(in), value :: n_axes
      !! Number of axes
   real(c_double), dimension(n_axes), intent(out) :: contributions
      !! Relative axis contributions (output), values in [0,1], sum to 1
   integer(c_int), intent(out) :: ierr
      !! Error code

   call relative_axes_changes_from_shift_vector(vec, n_axes, contributions, ierr)
end subroutine relative_axes_changes_from_shift_vector_c

subroutine relative_axes_expression_from_expression_vector_r(vec, n_axes, contributions, ierr)
   use relative_axis_plane_tools, only: relative_axes_expression_from_expression_vector
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_axes), intent(in) :: vec
      !! RAP-projected and normalized expression vector
   integer(int32), intent(in) :: n_axes
      !! Number of axes
   real(real64), dimension(n_axes), intent(out) :: contributions
      !! Relative axis contributions (output), values in [0,1], sum to 1
   integer(int32), intent(out) :: ierr
      !! Error code

   call relative_axes_expression_from_expression_vector(vec, n_axes, contributions, ierr)
end subroutine relative_axes_expression_from_expression_vector_r

subroutine relative_axes_expression_from_expression_vector_c(vec, n_axes, contributions, ierr) bind(C, name="relative_axes_expression_from_expression_vector_c")
   use iso_c_binding, only: c_double, c_int
   use relative_axis_plane_tools, only: relative_axes_expression_from_expression_vector
   implicit none

   real(c_double), dimension(n_axes), intent(in) :: vec
      !! RAP-projected and normalized expression vector
   integer(c_int), intent(in), value :: n_axes
      !! Number of axes
   real(c_double), dimension(n_axes), intent(out) :: contributions
      !! Relative axis contributions (output), values in [0,1], sum to 1
   integer(c_int), intent(out) :: ierr
      !! Error code

   call relative_axes_expression_from_expression_vector(vec, n_axes, contributions, ierr)
end subroutine relative_axes_expression_from_expression_vector_c

subroutine omics_vector_RAP_projection_r(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
   use relative_axis_plane_tools, only: omics_vector_RAP_projection
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
      !! matrix with expression vectors
   integer(int32), intent(in) :: n_axes
      !! number of axes
   integer(int32), intent(in) :: n_vecs
      !! number of vectors per axis
   logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! `.true.` for vectors where projection is to be computed
   integer(int32), intent(in) :: n_selected_vecs
      !! count of `.true.` values in `vecs_selection_mask`
   logical, dimension(n_axes), intent(in) :: axes_selection_mask
      !! `.true.` for axes to be included in RAP
   integer(int32), intent(in) :: n_selected_axes
      !! count of `.true.` values in `axes_selection_mask`
   real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
      !! projected vectors
   integer(int32), intent(out) :: ierr
      !! Error code

   call omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
end subroutine omics_vector_RAP_projection_r

subroutine omics_vector_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr) bind(C, name="omics_vector_RAP_projection_c")
   use iso_c_binding, only: c_double, c_int
   use relative_axis_plane_tools, only: omics_vector_RAP_projection
   implicit none

   real(c_double), dimension(n_axes, n_vecs), intent(in) :: vecs
      !! matrix with expression vectors
   integer(c_int), intent(in), value :: n_axes
      !! number of axes
   integer(c_int), intent(in), value :: n_vecs
      !! number of vectors per axis
   integer(c_int), dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! `.true.` for vectors where projection is to be computed
   integer(c_int), intent(in), value :: n_selected_vecs
      !! count of `.true.` values in `vecs_selection_mask`
   integer(c_int), dimension(n_axes), intent(in) :: axes_selection_mask
      !! `.true.` for axes to be included in RAP
   integer(c_int), intent(in), value :: n_selected_axes
      !! count of `.true.` values in `axes_selection_mask`
   real(c_double), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
      !! projected vectors
   integer(c_int), intent(out) :: ierr
      !! Error code

   call omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask /= 0, n_selected_vecs, axes_selection_mask /= 0, n_selected_axes, projections, ierr)
end subroutine omics_vector_RAP_projection_c

subroutine omics_field_RAP_projection_r(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
   use relative_axis_plane_tools, only: omics_field_RAP_projection
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
      !! matrix with vector fields, first n rows mean vector origin, last n rows vector targets
   integer(int32), intent(in) :: n_axes
      !! number of axes
   integer(int32), intent(in) :: n_vecs
      !! number of vectors per axis
   logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! `.true.` for vectors where projection is to be computed
   integer(int32), intent(in) :: n_selected_vecs
      !! count of `.true.` values in `vecs_selection_mask`
   logical, dimension(n_axes), intent(in) :: axes_selection_mask
      !! `.true.` for axes to be included in RAP
   integer(int32), intent(in) :: n_selected_axes
      !! count of `.true.` values in `axes_selection_mask`
   real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
      !! projected vectors
   integer(int32), intent(out) :: ierr
      !! Error code

   call omics_field_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
end subroutine omics_field_RAP_projection_r

subroutine omics_field_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr) bind(C, name="omics_field_RAP_projection_c")
   use iso_c_binding, only: c_double, c_int
   use relative_axis_plane_tools, only: omics_field_RAP_projection
   implicit none

   real(c_double), dimension(2 * n_axes, n_vecs), intent(in) :: vecs
      !! matrix with vector fields, first n rows mean vector origin, last n rows vector targets
   integer(c_int), intent(in), value :: n_axes
      !! number of axes
   integer(c_int), intent(in), value :: n_vecs
      !! number of vectors per axis
   integer(c_int), dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! `.true.` for vectors where projection is to be computed
   integer(c_int), intent(in), value :: n_selected_vecs
      !! count of `.true.` values in `vecs_selection_mask`
   integer(c_int), dimension(n_axes), intent(in) :: axes_selection_mask
      !! `.true.` for axes to be included in RAP
   integer(c_int), intent(in), value :: n_selected_axes
      !! count of `.true.` values in `axes_selection_mask`
   real(c_double), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
      !! projected vectors
   integer(c_int), intent(out) :: ierr
      !! Error code

   call omics_field_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask /= 0, n_selected_vecs, axes_selection_mask /= 0, n_selected_axes, projections, ierr)
end subroutine omics_field_RAP_projection_c

subroutine clock_hand_angle_between_vectors_r(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
   use relative_axis_plane_tools, only: clock_hand_angle_between_vectors
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_dims), intent(in) :: v1
      !! First normalized vector in RAP space
   real(real64), dimension(n_dims), intent(in) :: v2
      !! Second normalized vector in RAP space  
   integer(int32), intent(in) :: n_dims
      !! Dimension of both vectors
   real(real64), intent(out) :: signed_angle
      !! Signed angle between vectors in radians [-π, π]
   integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
   integer(int32), intent(out) :: ierr  
      !! Error code

   call clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
end subroutine clock_hand_angle_between_vectors_r

subroutine clock_hand_angle_between_vectors_c(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr) bind(C, name="clock_hand_angle_between_vectors_c")
   use iso_c_binding, only: c_double, c_int
   use relative_axis_plane_tools, only: clock_hand_angle_between_vectors
   implicit none

   real(c_double), dimension(n_dims), intent(in) :: v1
      !! First normalized vector in RAP space
   real(c_double), dimension(n_dims), intent(in) :: v2
      !! Second normalized vector in RAP space  
   integer(c_int), intent(in), value :: n_dims
      !! Dimension of both vectors
   real(c_double), intent(out) :: signed_angle
      !! Signed angle between vectors in radians [-π, π]
   integer(c_int), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
   integer(c_int), intent(out) :: ierr  
      !! Error code

   call clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
end subroutine clock_hand_angle_between_vectors_c

subroutine clock_hand_angles_for_shift_vectors_r(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles, ierr)
   use relative_axis_plane_tools, only: clock_hand_angles_for_shift_vectors
   use, intrinsic :: iso_fortran_env, only: real64, int32
   implicit none

   real(real64), dimension(n_dims, n_vecs), intent(in) :: origins
      !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
   real(real64), dimension(n_dims, n_vecs), intent(in) :: targets
      !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
   integer(int32), intent(in) :: n_dims
      !! Dimension of each vector in RAP space
   integer(int32), intent(in) :: n_vecs
      !! Number of vector pairs
   logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! .true. for vector pairs where angle should be computed
      integer(int32), intent(in) :: n_selected_vecs
         !! Count of .true. values in vecs_selection_mask
      integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
         !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
      real(real64), dimension(n_selected_vecs), intent(out) :: signed_angles
         !! Signed rotation angles between vector pairs in radians [-π, π]
      integer(int32), intent(out) :: ierr  
      !! Error code

   call clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles, ierr)
end subroutine clock_hand_angles_for_shift_vectors_r

subroutine clock_hand_angles_for_shift_vectors_c(origins, targets, n_dims, n_vecs, vecs_selection_mask, n_selected_vecs, selected_axes_for_signed, signed_angles, ierr) bind(C, name="clock_hand_angles_for_shift_vectors_c")
   use iso_c_binding, only: c_double, c_int
   use relative_axis_plane_tools, only: clock_hand_angles_for_shift_vectors
   implicit none

   real(c_double), dimension(n_dims, n_vecs), intent(in) :: origins
      !! First set of RAP-projected, normalized vectors (e.g. expression centroids)
   real(c_double), dimension(n_dims, n_vecs), intent(in) :: targets
      !! Second set of RAP-projected, normalized vectors (e.g. paralogs)
   integer(c_int), intent(in), value :: n_dims
      !! Dimension of each vector in RAP space
   integer(c_int), intent(in), value :: n_vecs
      !! Number of vector pairs
   integer(c_int), dimension(n_vecs), intent(in) :: vecs_selection_mask
      !! .true. for vector pairs where angle should be computed
   integer(c_int), intent(in), value :: n_selected_vecs
      !! Count of .true. values in vecs_selection_mask
   integer(c_int), dimension(3), intent(in) :: selected_axes_for_signed
      !! Indices of 3 axes to use for directionality calculation (ignored if n_dims <= 3)
   real(c_double), dimension(n_selected_vecs), intent(out) :: signed_angles
      !! Signed rotation angles between vector pairs in radians [-π, π]
   integer(c_int), intent(out) :: ierr  
      !! Error code

   call clock_hand_angles_for_shift_vectors(origins, targets, n_dims, n_vecs, vecs_selection_mask /= 0, n_selected_vecs, selected_axes_for_signed, signed_angles, ierr)
end subroutine clock_hand_angles_for_shift_vectors_c
