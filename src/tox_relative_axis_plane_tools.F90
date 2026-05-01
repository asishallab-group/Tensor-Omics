#include "macros.h"

!> Module for tools related to relative axis planes (RAPs), i.e. planes in higher-dimensional gene expression space
module tox_relative_axis_plane_tools
    use safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_INVALID_INPUT, set_ok, set_err, is_err, validate_dimension_size, validate_all_in_range_real, validate_all_in_range_int, set_err_once, ERR_DIVISION_BY_ZERO
    use f42_utils, only: operator(.isclose.)
    implicit none

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    pure subroutine omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr)
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_vecs
            !! number of vectors per axis
        integer(int32), intent(in) :: n_selected_vecs
            !! count of `.true.` values in `vecs_selection_mask`
        integer(int32), intent(in) :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(real64), dimension(n_axes, n_vecs), intent(in) :: vecs
            !! matrix with expression vectors
        logical, dimension(n_vecs), intent(in) :: vecs_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        real(real64), dimension(n_selected_axes, n_selected_vecs), intent(out) :: projections
            !! projected vectors
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_vecs, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_axes, ierr, arg_pos=7_int32)
        call validate_dimension_size(n_selected_vecs, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(vecs, size(vecs, kind=int32), ierr, arg_pos=1_int32)
        if (count(vecs_selection_mask, kind=int32) > n_selected_vecs) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (count(axes_selection_mask, kind=int32) > n_selected_axes) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)

        if (is_err(ierr)) return

        call omics_vector_RAP_projection_helper(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
    end subroutine omics_vector_RAP_projection

    !> AUTHOR_FRANZ_ERIC_SILL
    !| (no input validation) Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    pure subroutine omics_vector_RAP_projection_helper(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
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

        integer(int32) :: i_vec, i_axis, i_vec_proj, i_axis_proj

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

        call project_selected_vecs_onto_rap_helper(projections, n_selected_axes, n_selected_vecs)
    end subroutine omics_vector_RAP_projection_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    pure subroutine omics_field_RAP_projection(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections, ierr)
        real(real64), dimension(n_axes, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_fields
            !! number of vectors per axis
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        integer(int32), intent(in) :: n_selected_fields
            !! count of `.true.` values in `fields_selection_mask`
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        integer(int32), intent(in) :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(real64), dimension(n_selected_axes, n_selected_fields), intent(out) :: projections
            !! projected vectors
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_fields, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_axes, ierr, arg_pos=7_int32)
        call validate_dimension_size(n_selected_fields, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(fields, size(fields, kind=int32), ierr, arg_pos=1_int32)
        if (count(fields_selection_mask, kind=int32) > n_selected_fields) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=5_int32)
        if (count(axes_selection_mask, kind=int32) > n_selected_axes) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=7_int32)

        if (is_err(ierr)) return

        call omics_field_RAP_projection_helper(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections)
    end subroutine omics_field_RAP_projection

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    pure subroutine omics_field_RAP_projection_helper(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections)
        real(real64), dimension(n_axes, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_fields
            !! number of vectors per axis
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        integer(int32), intent(in) :: n_selected_fields
            !! count of `.true.` values in `fields_selection_mask`
        logical, dimension(n_axes), intent(in) :: axes_selection_mask
            !! `.true.` for axes to be included in RAP
        integer(int32), intent(in) :: n_selected_axes
            !! count of `.true.` values in `axes_selection_mask`
        real(real64), dimension(n_selected_axes, n_selected_fields), intent(out) :: projections
            !! projected vectors

        integer(int32) :: i_vec, i_axis, i_vec_proj, i_axis_proj
        i_vec_proj = 1
        do i_vec = 1, n_fields
            if (fields_selection_mask(i_vec)) then
                i_axis_proj = 1
                do i_axis = 1, n_axes
                    if (axes_selection_mask(i_axis)) then
                        ! compute shift vector as difference between origin and target
                        projections(i_axis_proj, i_vec_proj) = fields(i_axis, 1, i_vec) - fields(i_axis, 2, i_vec)

                        i_axis_proj = i_axis_proj + 1
                    end if
                end do

                i_vec_proj = i_vec_proj + 1
            end if
        end do

        call project_selected_vecs_onto_rap_helper(projections, n_selected_axes, n_selected_fields)
    end subroutine omics_field_RAP_projection_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Projects selected vectors onto its RAP
    pure subroutine project_selected_vecs_onto_rap_helper(selected_vecs, n_selected_axes, n_selected_vecs)
        real(real64), dimension(n_selected_axes, n_selected_vecs), intent(inout) :: selected_vecs
            !! matrix with vectors for selected axes
        integer(int32), intent(in) :: n_selected_axes
            !! number of selected axes
        integer(int32), intent(in) :: n_selected_vecs
            !! number of selected vectors per axis

        ! project selected vectors onto RAP
        integer(int32) :: i_vec, i_axis
        real(real64) :: diagonal_component

        do concurrent (i_vec = 1:n_selected_vecs) local(diagonal_component) shared(n_selected_axes, selected_vecs)

            ! calculate diagonal component to be subtracted from vectors for projection
            diagonal_component = 0.0_real64
            do concurrent (i_axis = 1:n_selected_axes) shared(selected_vecs, i_vec) reduce(+:diagonal_component)
                diagonal_component = diagonal_component + selected_vecs(i_axis, i_vec)
            end do
            diagonal_component = diagonal_component / n_selected_axes

            ! transform vector to its projection
            do concurrent (i_axis = 1:n_selected_axes) shared(selected_vecs, i_vec, diagonal_component)
                selected_vecs(i_axis, i_vec) = selected_vecs(i_axis, i_vec) - diagonal_component
            end do
        end do
    end subroutine project_selected_vecs_onto_rap_helper

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Does the validation of the `selected_axes_for_signed` input
    pure subroutine validate_selected_axes_for_signed_helper(selected_axes_for_signed, ierr, n_dims, arg_pos)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
        integer(int32), intent(inout) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: arg_pos
            !! Position of the validated argument that triggered the error, default: 0 -> not argument related

        ! the argument is being ignored for lower 4
        if (n_dims > 3) then
            call validate_all_in_range_int(selected_axes_for_signed, size(selected_axes_for_signed, kind=int32), ierr, min=1_int32, max=n_dims, arg_pos=5_int32)
            if (selected_axes_for_signed(1) == selected_axes_for_signed(2) .or. &
                selected_axes_for_signed(1) == selected_axes_for_signed(3) .or. &
                selected_axes_for_signed(2) == selected_axes_for_signed(3)) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
                return
            end if
        end if
    end subroutine validate_selected_axes_for_signed_helper

    !> AUTHOR_VIVIAN_BASS
    !| Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    !| Calculates the signed rotation angle between two normalized vectors in RAP space.
    !| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
    pure subroutine clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second normalized vector in RAP space
        real(real64), intent(out) :: signed_angle
            !! Signed angle between vectors in radians [-π, π]
        integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dims, ierr, arg_pos=3_int32)
        call validate_all_in_range_real(v1, n_dims, ierr, arg_pos=1_int32)
        call validate_all_in_range_real(v2, n_dims, ierr, arg_pos=2_int32)
        call validate_selected_axes_for_signed_helper(selected_axes_for_signed, ierr, n_dims, arg_pos=5_int32)

        if (is_err(ierr)) return

        call clock_hand_angle_between_vectors_helper(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
    end subroutine clock_hand_angle_between_vectors

    !> AUTHOR_VIVIAN_BASS
    !| (no input validation) Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    !| Calculates the signed rotation angle between two normalized vectors in RAP space.
    !| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
    pure subroutine clock_hand_angle_between_vectors_helper(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second normalized vector in RAP space
        real(real64), intent(out) :: signed_angle
            !! Signed angle between vectors in radians [-π, π]
        integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)

        real(real64) :: dot_product, unsigned_angle, orientation_sign
        integer(int32) :: i_dim

        ! Calculate dot product of normalized vectors
        dot_product = 0.0_real64
        do concurrent (i_dim = 1:n_dims) shared(v1, v2) reduce(+:dot_product)
            dot_product = dot_product + v1(i_dim) * v2(i_dim)
        end do

        ! Clamp dot product to [-1, 1] to handle numerical precision issues
        dot_product = max(-1.0_real64, min(1.0_real64, dot_product))

        ! Calculate unsigned angle using arccos
        unsigned_angle = acos(dot_product)

        ! Calculate orientation sign for directionality
        select case (n_dims)
        case (2)
            ! For 2D: use determinant directly
            orientation_sign = sign(1.0_real64, v1(1)*v2(2) - v1(2)*v2(1))
        case (3)
            ! For 3D, use [1,2,3] directly
            orientation_sign = cross_product_orientation_sign(v1, v2, n_dims, [1_int32, 2_int32, 3_int32])
        case (4:)
            ! For >3D, use selected_axes_for_signed
            orientation_sign = cross_product_orientation_sign(v1, v2, n_dims, selected_axes_for_signed)
        case default
            orientation_sign = 1.0_real64
        end select

        ! Apply sign to unsigned angle
        signed_angle = orientation_sign * unsigned_angle
    end subroutine clock_hand_angle_between_vectors_helper

    !> AUTHOR_VIVIAN_BASS
    !| Compute signed rotation angles between for shift vectors, so between their origin and target
    pure subroutine clock_hand_angles_for_shift_vectors(fields, n_dims, n_fields, &
                                                        fields_selection_mask, &
                                                        n_selected_fields, selected_axes_for_signed, &
                                                        signed_angles, ierr)
        integer(int32), intent(in) :: n_dims
            !! Dimension of each vector in RAP space
        integer(int32), intent(in) :: n_fields
            !! Number of vector fields
        real(real64), dimension(n_dims, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! .true. for vector pairs where angle should be computed
        integer(int32), intent(in) :: n_selected_fields
            !! Count of .true. values in fields_selection_mask
        integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
        real(real64), dimension(n_selected_fields), intent(out) :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)

        call validate_dimension_size(n_dims, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_fields, ierr, arg_pos=3_int32)
        call validate_dimension_size(n_selected_fields, ierr, arg_pos=5_int32)
        call validate_all_in_range_real(fields, size(fields, kind=int32), ierr, arg_pos=1_int32)
        if (count(fields_selection_mask, kind=int32) > n_selected_fields) call set_err(ierr, ERR_INVALID_INPUT, arg_pos=4_int32)
        call validate_selected_axes_for_signed_helper(selected_axes_for_signed, ierr, n_dims, arg_pos=6_int32)

        if (is_err(ierr)) return

        call clock_hand_angles_for_shift_vectors_helper(fields, n_dims, n_fields, &
                                                        fields_selection_mask, &
                                                        n_selected_fields, selected_axes_for_signed, &
                                                        signed_angles)
    end subroutine clock_hand_angles_for_shift_vectors

    !> AUTHOR_VIVIAN_BASS
    !| (no input validation) Compute signed rotation angles between for shift vectors, so between their origin and target
    pure subroutine clock_hand_angles_for_shift_vectors_helper(fields, n_dims, n_fields, &
                                                        fields_selection_mask, &
                                                        n_selected_fields, selected_axes_for_signed, &
                                                        signed_angles)
        real(real64), dimension(n_dims, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
        integer(int32), intent(in) :: n_dims
            !! Dimension of each vector in RAP space
        integer(int32), intent(in) :: n_fields
            !! Number of vector pairs
        logical, dimension(n_fields), intent(in) :: fields_selection_mask
            !! .true. for vector pairs where angle should be computed
        integer(int32), intent(in) :: n_selected_fields
            !! Count of .true. values in fields_selection_mask
        integer(int32), dimension(3), intent(in) :: selected_axes_for_signed
            !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
        real(real64), dimension(n_selected_fields), intent(out) :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]

        integer(int32) :: i_field, result_idx

        result_idx = 1
        do i_field = 1, n_fields
            if (fields_selection_mask(i_field)) then
                call clock_hand_angle_between_vectors_helper(fields(:, 1, i_field), fields(:, 2, i_field), n_dims, &
                                                      signed_angles(result_idx), selected_axes_for_signed)
                result_idx = result_idx + 1
            end if
        end do
    end subroutine clock_hand_angles_for_shift_vectors_helper

    !> AUTHOR_VIVIAN_BASS
    !| Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
    pure subroutine compute_relative_axis_contributions(vec, n_axes, contributions, ierr)
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized vector (expression or shift)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (length of vec and contributions)
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Error handling
        call set_ok(ierr)

        call validate_dimension_size(n_axes, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(vec, n_axes, ierr, arg_pos=1_int32)

        if (is_err(ierr)) return

        call compute_relative_axis_contributions_helper(vec, n_axes, contributions, ierr)
    end subroutine compute_relative_axis_contributions

    !> AUTHOR_VIVIAN_BASS
    !| (no input validation) Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| Shared utility: computes fractional contribution of each axis to a RAP-projected and normalized vector.
    pure subroutine compute_relative_axis_contributions_helper(vec, n_axes, contributions, ierr)
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

        total_abs = 0.0_real64
        do concurrent (i_axis = 1:n_axes) shared(vec) reduce(+:total_abs)
            total_abs = total_abs + abs(vec(i_axis))
        end do

        if (total_abs .isclose. 0.0_real64) then
            contributions = 0.0_real64
            call set_err_once(ierr, ERR_DIVISION_BY_ZERO, arg_pos=1_int32)
            return
        end if

        do concurrent (i_axis = 1:n_axes) shared(contributions, vec, total_abs)
            contributions(i_axis) = abs(vec(i_axis))/total_abs
        end do
    end subroutine compute_relative_axis_contributions_helper

    !> AUTHOR_VIVIAN_BASS
    !| Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
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

    !> AUTHOR_VIVIAN_BASS
    !| Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
    !| Wrapper for single RAP-projected expression vectors
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

    !> AUTHOR_VIVIAN_BASS
    !| Compute orientation sign from cross product of two vectors, using selected axes
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

end module tox_relative_axis_plane_tools

! Updated wrappers to pass and return ierr

pure subroutine clock_hand_angle_between_vectors_c(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr) bind(C, name="clock_hand_angle_between_vectors_c")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_relative_axis_plane_tools, only: clock_hand_angle_between_vectors
    M_USE_NULL_VALIDATION
    implicit none

    real(c_double), dimension(n_dims), intent(in), target :: v1
        !! First normalized vector in RAP space
    real(c_double), dimension(n_dims), intent(in), target :: v2
        !! Second normalized vector in RAP space
    integer(c_int), intent(in), target :: n_dims
        !! Dimension of both vectors
    real(c_double), intent(out), target :: signed_angle
        !! Signed angle between vectors in radians [-π, π]
    integer(c_int), dimension(3), intent(in), target :: selected_axes_for_signed
        !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(v1)
    M_CHECK_NON_NULL(v2)
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(signed_angle)
    M_CHECK_NON_NULL(selected_axes_for_signed)

    call clock_hand_angle_between_vectors(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
end subroutine clock_hand_angle_between_vectors_c

pure subroutine clock_hand_angles_for_shift_vectors_c(fields, n_dims, n_fields, fields_selection_mask, n_selected_fields, selected_axes_for_signed, signed_angles, ierr) bind(C, name="clock_hand_angles_for_shift_vectors_c")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_relative_axis_plane_tools, only: clock_hand_angles_for_shift_vectors
    use tox_conversions, only: c_int_as_logical
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_dims
        !! Dimension of each vector in RAP space
    integer(c_int), intent(in), target :: n_fields
        !! Number of vector pairs
    real(c_double), dimension(n_dims, 2, n_fields), intent(in), target :: fields
        !! matrix with vector fields, `fields(:, 1, i_vec)` mean vector origin, `fields(:, 2, i_vec)` mean vector targets
    integer(c_int), dimension(n_fields), intent(in), target :: fields_selection_mask
        !! .true. for vector pairs where angle should be computed
    integer(c_int), intent(in), target :: n_selected_fields
        !! Count of .true. values in fields_selection_mask
    integer(c_int), dimension(3), intent(in), target :: selected_axes_for_signed
        !! Indices of 3 different axes to use for directionality calculation (ignored if n_dims <= 3, all indices must be unique)
    real(c_double), dimension(n_selected_fields), intent(out), target :: signed_angles
        !! Signed rotation angles between vector pairs in radians [-π, π]
    integer(c_int), intent(out), target :: ierr
        !! Error code

    logical, dimension(:), allocatable :: fields_selection_mask_f(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_dims)
    M_CHECK_NON_NULL(n_fields)
    M_CHECK_NON_NULL(n_selected_fields)
    M_CHECK_NON_NULL(fields)
    M_CHECK_NON_NULL(fields_selection_mask)
    M_CHECK_NON_NULL(selected_axes_for_signed)
    M_CHECK_NON_NULL(signed_angles)

    ! Convert c_int mask to logical using tox_conversions utility
    M_ALLOCATE(fields_selection_mask_f(n_fields))
    call c_int_as_logical(fields_selection_mask, fields_selection_mask_f)

    call clock_hand_angles_for_shift_vectors(fields, n_dims, n_fields, fields_selection_mask_f, n_selected_fields, selected_axes_for_signed, signed_angles, ierr)
end subroutine clock_hand_angles_for_shift_vectors_c

!> C/Python wrapper for omics_vector_RAP_projection
pure subroutine omics_vector_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr) bind(C, name="omics_vector_RAP_projection_c")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_relative_axis_plane_tools, only: omics_vector_RAP_projection
    use tox_conversions, only: c_int_as_logical
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none
    integer(c_int), intent(in), target :: n_axes
        !! number of axes
    integer(c_int), intent(in), target :: n_vecs
        !! number of vectors
    real(c_double), dimension(n_axes, n_vecs), intent(in), target :: vecs
        !! matrix with expression vectors
    integer(c_int), dimension(n_vecs), intent(in), target :: vecs_selection_mask
        !! mask for selecting vectors
    integer(c_int), intent(in), target :: n_selected_vecs
        !! count of selected vectors
    integer(c_int), dimension(n_axes), intent(in), target :: axes_selection_mask
        !! mask for selecting axes
    integer(c_int), intent(in), target :: n_selected_axes
        !! count of selected axes
    real(c_double), dimension(n_selected_axes, n_selected_vecs), intent(out), target :: projections
        !! projected vectors
    integer(c_int), intent(out), target :: ierr
        !! error code

    logical, dimension(:), allocatable :: vecs_selection_mask_f(:), axes_selection_mask_f(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_vecs)
    M_CHECK_NON_NULL(n_selected_vecs)
    M_CHECK_NON_NULL(n_selected_axes)
    M_CHECK_NON_NULL(vecs)
    M_CHECK_NON_NULL(vecs_selection_mask)
    M_CHECK_NON_NULL(axes_selection_mask)
    M_CHECK_NON_NULL(projections)

    M_ALLOCATE(vecs_selection_mask_f(n_vecs))
    call c_int_as_logical(vecs_selection_mask, vecs_selection_mask_f)
    M_ALLOCATE(axes_selection_mask_f(n_axes))
    call c_int_as_logical(axes_selection_mask, axes_selection_mask_f)

    call omics_vector_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask_f, n_selected_vecs, axes_selection_mask_f, n_selected_axes, projections, ierr)
end subroutine omics_vector_RAP_projection_c

!> C/Python wrapper for omics_field_RAP_projection
pure subroutine omics_field_RAP_projection_c(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections, ierr) bind(C, name="omics_field_RAP_projection_c")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_relative_axis_plane_tools, only: omics_field_RAP_projection
    use tox_conversions, only: c_int_as_logical
    use tox_errors, only: is_err, set_err, ERR_ALLOC_FAIL
    M_USE_NULL_VALIDATION
    implicit none

    integer(c_int), intent(in), target :: n_axes
        !! number of axes
    integer(c_int), intent(in), target :: n_vecs
        !! number of vectors
    real(c_double), dimension(n_axes, 2, n_vecs), intent(in), target :: vecs
        !! matrix with vector fields, first n rows mean vector origin, last n rows vector targets
    integer(c_int), dimension(n_vecs), intent(in), target :: vecs_selection_mask
        !! mask for selecting vectors
    integer(c_int), intent(in), target :: n_selected_vecs
        !! count of selected vectors
    integer(c_int), dimension(n_axes), intent(in), target :: axes_selection_mask
        !! mask for selecting axes
    integer(c_int), intent(in), target :: n_selected_axes
        !! count of selected axes
    real(c_double), dimension(n_selected_axes, n_selected_vecs), intent(out), target :: projections
        !! projected vectors
    integer(c_int), intent(out), target :: ierr
        !! error code

    logical, dimension(:), allocatable :: vecs_selection_mask_f(:), axes_selection_mask_f(:)

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(n_vecs)
    M_CHECK_NON_NULL(n_selected_vecs)
    M_CHECK_NON_NULL(n_selected_axes)
    M_CHECK_NON_NULL(vecs)
    M_CHECK_NON_NULL(vecs_selection_mask)
    M_CHECK_NON_NULL(axes_selection_mask)
    M_CHECK_NON_NULL(projections)

    M_ALLOCATE(vecs_selection_mask_f(n_vecs))
    call c_int_as_logical(vecs_selection_mask, vecs_selection_mask_f)
    M_ALLOCATE(axes_selection_mask_f(n_axes))
    call c_int_as_logical(axes_selection_mask, axes_selection_mask_f)

    call omics_field_RAP_projection(vecs, n_axes, n_vecs, vecs_selection_mask_f, n_selected_vecs, axes_selection_mask_f, n_selected_axes, projections, ierr)
end subroutine omics_field_RAP_projection_c

!> C/Python wrapper for relative_axes_changes_from_shift_vector
pure subroutine relative_axes_changes_from_shift_vector_c(vec, n_axes, contributions, ierr) bind(C, name="relative_axes_changes_from_shift_vector_c")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_relative_axis_plane_tools, only: relative_axes_changes_from_shift_vector
    M_USE_NULL_VALIDATION
    implicit none

    real(c_double), dimension(n_axes), intent(in), target :: vec
        !! RAP-projected and normalized shift vector
    integer(c_int), intent(in), target :: n_axes
        !! Number of axes
    real(c_double), dimension(n_axes), intent(out), target :: contributions
        !! Fractional contribution of each axis
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(vec)
    M_CHECK_NON_NULL(contributions)

    call relative_axes_changes_from_shift_vector(vec, n_axes, contributions, ierr)
end subroutine relative_axes_changes_from_shift_vector_c

!> C/Python wrapper for relative_axes_expression_from_expression_vector
pure subroutine relative_axes_expression_from_expression_vector_c(vec, n_axes, contributions, ierr) bind(C, name="relative_axes_expression_from_expression_vector_c")
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    use tox_relative_axis_plane_tools, only: relative_axes_expression_from_expression_vector
    M_USE_NULL_VALIDATION
    implicit none

    real(c_double), dimension(n_axes), intent(in), target :: vec
        !! RAP-projected and normalized expression vector
    integer(c_int), intent(in), target :: n_axes
        !! Number of axes
    real(c_double), dimension(n_axes), intent(out), target :: contributions
        !! Fractional contribution of each axis
    integer(c_int), intent(out), target :: ierr
        !! Error code

    M_CHECK_IERR_NON_NULL
    M_CHECK_NON_NULL(n_axes)
    M_CHECK_NON_NULL(vec)
    M_CHECK_NON_NULL(contributions)

    call relative_axes_expression_from_expression_vector(vec, n_axes, contributions, ierr)
end subroutine relative_axes_expression_from_expression_vector_c
