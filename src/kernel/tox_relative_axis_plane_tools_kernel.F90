#include <src/macros.h>

!> Kernels for tools related to relative axis planes (RAPs), i.e. planes in higher-dimensional gene expression space.
!| The generator turns the `*_kernel` procedures into the validating wrappers omics_vector_RAP_projection,
!| omics_field_RAP_projection, clock_hand_angle_between_vectors, clock_hand_angles_for_shift_vectors,
!| relative_axes_changes_from_shift_vector and relative_axes_expression_from_expression_vector in module
!| tox_relative_axis_plane_tools. The kernels keep only the checks a per-argument validator cannot express
!| (a count(mask) matching its claimed size, and the unique-axes check); the mechanical dimension and
!| finiteness checks come from the generated wrappers. The compute helpers below carry no `_kernel` suffix
!| and are shared by several kernels, so they stay here untouched.
module tox_relative_axis_plane_tools_kernel
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: ERR_INVALID_INPUT, set_ok, is_err, validate_all_in_range_int, set_err_once, ERR_DIVISION_BY_ZERO
    use f42_utils, only: operator(.isclose.)
    M_IMPLICIT_NONE

contains

    !> summary: Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine omics_vector_RAP_projection_kernel(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
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
    end subroutine omics_vector_RAP_projection_kernel

    !> summary: Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine omics_field_RAP_projection_kernel(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections)
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
    end subroutine omics_field_RAP_projection_kernel

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

        ! The RAP is the hyperplane orthogonal to the space diagonal (1,1,...,1) through the origin,
        ! i.e. the set of points whose coordinates sum to zero. Orthogonal projection of a vector v
        ! onto that hyperplane is v - ((v . 1)/(1 . 1)) * 1, and since 1 . 1 = n_selected_axes and
        ! every component of the diagonal unit direction is equal, this reduces to subtracting the
        ! mean of v's coordinates ("diagonal_component") from each coordinate of v.
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
            call validate_all_in_range_int(selected_axes_for_signed, size(selected_axes_for_signed, kind=int32), ierr, min=1_int32, max=n_dims, arg_pos=arg_pos)
            if (selected_axes_for_signed(1) == selected_axes_for_signed(2) .or. &
                selected_axes_for_signed(1) == selected_axes_for_signed(3) .or. &
                selected_axes_for_signed(2) == selected_axes_for_signed(3)) then
                call set_err_once(ierr, ERR_INVALID_INPUT, arg_pos)
                return
            end if
        end if
    end subroutine validate_selected_axes_for_signed_helper

    !> summary: Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    !| AUTHOR_VIVIAN_BASS
    !| Calculates the signed rotation angle between two normalized vectors in RAP space.
    !| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
    pure subroutine clock_hand_angle_between_vectors_kernel(v1, v2, n_dims, signed_angle, selected_axes_for_signed, ierr)
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

        ! Unique-axes check the generated wrapper's per-argument validators cannot express.
        call validate_selected_axes_for_signed_helper(selected_axes_for_signed, ierr, n_dims, arg_pos=5_int32)
        if (is_err(ierr)) return

        call clock_hand_angle_between_vectors_helper(v1, v2, n_dims, signed_angle, selected_axes_for_signed)
    end subroutine clock_hand_angle_between_vectors_kernel

    !> AUTHOR_VIVIAN_BASS
    !| Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    !| Calculates the signed rotation angle between two normalized vectors in RAP space.
    !| For 2D/3D: automatic directionality calculation. For >3D: uses selected axes for directionality.
    !| Shared compute core: the single-vector angle, reused per field by
    !| [[tox_relative_axis_plane_tools_kernel(module):clock_hand_angles_for_shift_vectors_kernel(subroutine)]].
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
            block
                integer(int32), dimension(3), parameter :: default_selected_axes = [1_int32, 2_int32, 3_int32]

                ! For 3D, use [1,2,3] directly
                orientation_sign = cross_product_orientation_sign(v1, v2, n_dims, default_selected_axes)
            end block
        case (4:)

            ! For >3D, use selected_axes_for_signed
            orientation_sign = cross_product_orientation_sign(v1, v2, n_dims, selected_axes_for_signed)
        case default
            orientation_sign = 1.0_real64
        end select

        ! Apply sign to unsigned angle
        signed_angle = orientation_sign * unsigned_angle
    end subroutine clock_hand_angle_between_vectors_helper

    !> summary: Compute signed rotation angles between for shift vectors, so between their origin and target
    !| AUTHOR_VIVIAN_BASS
    pure subroutine clock_hand_angles_for_shift_vectors_kernel(fields, n_dims, n_fields, &
                                                        fields_selection_mask, &
                                                        n_selected_fields, selected_axes_for_signed, &
                                                        signed_angles, ierr)
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
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_field, result_idx

        ! The count consistency is a convention check in the wrapper; the unique-axes check is not
        ! expressible there, so it stays here (this argument's design is due for a rework).
        call validate_selected_axes_for_signed_helper(selected_axes_for_signed, ierr, n_dims, arg_pos=6_int32)
        if (is_err(ierr)) return

        result_idx = 1
        do i_field = 1, n_fields
            if (fields_selection_mask(i_field)) then
                call clock_hand_angle_between_vectors_helper(fields(:, 1, i_field), fields(:, 2, i_field), n_dims, &
                                                      signed_angles(result_idx), selected_axes_for_signed)
                result_idx = result_idx + 1
            end if
        end do
    end subroutine clock_hand_angles_for_shift_vectors_kernel

    !> AUTHOR_VIVIAN_BASS
    !| Compute fractional contribution of each axis to a RAP-projected and normalized vector.
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

    !> summary: Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| AUTHOR_VIVIAN_BASS
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    pure subroutine relative_axes_changes_from_shift_vector_kernel(vec, n_axes, contributions, ierr)
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized shift vector
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call compute_relative_axis_contributions_helper(vec, n_axes, contributions, ierr)
    end subroutine relative_axes_changes_from_shift_vector_kernel

    !> summary: Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
    !| AUTHOR_VIVIAN_BASS
    !| Wrapper for single RAP-projected expression vectors
    pure subroutine relative_axes_expression_from_expression_vector_kernel(vec, n_axes, contributions, ierr)
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized expression vector
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call compute_relative_axis_contributions_helper(vec, n_axes, contributions, ierr)
    end subroutine relative_axes_expression_from_expression_vector_kernel

    !> AUTHOR_VIVIAN_BASS
    !| Compute orientation sign from cross product of two vectors, using selected axes
    pure function cross_product_orientation_sign(a, b, n_dims, selected_axes) result(orientation_sign)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors `a` and `b`
        real(real64), intent(in) :: a(n_dims)
            !! First vector
        real(real64), intent(in) :: b(n_dims)
            !! Second vector
        integer(int32), intent(in) :: selected_axes(3)
            !! Indices of the 3 axes of `a`/`b` to project onto for the cross product (all must be unique)
        real(real64) :: orientation_sign
            !! `+1.0` or `-1.0`, the sign of the (selected-axes) cross product `a x b`

        real(real64) :: cross1, cross2, cross3

        ! Only 3 components define a cross product, so for n_dims>3 we project both vectors onto the
        ! 3 selected axes first, then compute the standard 3D cross product a x b on that subspace.
        ! Only the sign of the (summed) result is used, to get the rotation direction from a to b -
        ! the magnitude is not meaningful once the vectors have been projected down to 3 components.
        cross1 = a(selected_axes(2))*b(selected_axes(3)) - a(selected_axes(3))*b(selected_axes(2))
        cross2 = a(selected_axes(3))*b(selected_axes(1)) - a(selected_axes(1))*b(selected_axes(3))
        cross3 = a(selected_axes(1))*b(selected_axes(2)) - a(selected_axes(2))*b(selected_axes(1))

        orientation_sign = sign(1.0_real64, cross1 + cross2 + cross3)
    end function cross_product_orientation_sign

end module tox_relative_axis_plane_tools_kernel
