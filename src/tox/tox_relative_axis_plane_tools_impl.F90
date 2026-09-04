#include <src/macros.h>

!> Relative axis planes (RAPs): planes through higher-dimensional gene expression space, and
!| what can be read off a vector once it is projected onto one.
!|
!| A RAP is picked by selecting axes (tissues) from the full expression space.
!| `vector_RAP_projection` projects a single vector onto it and `field_RAP_projection` a whole
!| field of them. Within the plane, `clock_hand_angle_between_vectors` measures the signed angle
!| between two vectors -- signed by an orientation reference, so the sign means the same thing in
!| every dimension -- and `clock_hand_angles_for_shift_vectors` does that for a whole shift
!| vector field at once.
!|
!| `relative_axes_changes_from_shift_vector` and `relative_axes_expression_from_expression_vector`
!| give the per-axis breakdown instead of the angle: how much of a change, or of an expression
!| level, falls on each selected axis.
module tox_relative_axis_plane_tools_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: ERR_INVALID_INPUT, set_ok, is_err, set_err_once, ERR_DIVISION_BY_ZERO
    use f42_math_impl, only: clamp, operator(.isclose.)
    M_IMPLICIT_NONE

contains

    !> summary: Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine omics_vector_RAP_projection_impl(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
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
        logical(c_bool), dimension(n_vecs), intent(in) :: vecs_selection_mask
            !! `.true.` for vectors where projection is to be computed
        logical(c_bool), dimension(n_axes), intent(in) :: axes_selection_mask
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
    end subroutine omics_vector_RAP_projection_impl

    !> summary: Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine omics_field_RAP_projection_impl(fields, n_axes, n_fields, fields_selection_mask, n_selected_fields, axes_selection_mask, n_selected_axes, projections)
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_fields
            !! number of vectors per axis
        real(real64), dimension(n_axes, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical(c_bool), dimension(n_fields), intent(in) :: fields_selection_mask
            !! `.true.` for vectors where projection is to be computed
        integer(int32), intent(in) :: n_selected_fields
            !! count of `.true.` values in `fields_selection_mask`
        logical(c_bool), dimension(n_axes), intent(in) :: axes_selection_mask
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
    end subroutine omics_field_RAP_projection_impl

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Projects selected vectors onto its RAP
    pure subroutine project_selected_vecs_onto_rap_helper(selected_vecs, n_selected_axes, n_selected_vecs)
        integer(int32), intent(in) :: n_selected_axes
            !! number of selected axes
        integer(int32), intent(in) :: n_selected_vecs
            !! number of selected vectors per axis
        real(real64), dimension(n_selected_axes, n_selected_vecs), intent(inout) :: selected_vecs
            !! matrix with vectors for selected axes

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


    !> summary: Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    !| AUTHOR_VIVIAN_BASS
    !| The unsigned angle is `acos(v1 . v2)`; `orientation_reference` supplies the sign by saying
    !| which way round the plane the two vectors span counts as positive. Reports
    !| `ERR_INVALID_INPUT` when the reference is orthogonal to the rotation and so orients nothing.
    pure subroutine clock_hand_angle_between_vectors_impl(v1, v2, n_dims, orientation_reference, &
                                                            signed_angle, ierr)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(real64), intent(out) :: signed_angle
            !! Signed angle between vectors in radians [-pi, pi]
        integer(int32), intent(out) :: ierr
            !! Error code

        logical(c_bool) :: undefined_sign

        call set_ok(ierr)
        call clock_hand_angle_between_vectors_helper(v1, v2, n_dims, orientation_reference, &
                                                     signed_angle, undefined_sign)
        ! not a bad argument on its own -- the reference only fails to orient *this* rotation,
        ! which no check on any single argument could have foreseen
        if (undefined_sign) call set_err_once(ierr, ERR_INVALID_INPUT)
    end subroutine clock_hand_angle_between_vectors_impl

    !> AUTHOR_VIVIAN_BASS
    !| Compute the signed clock hand angle between two RAP-projected and normalized vectors.
    !|
    !| The unsigned angle is `acos(v1 . v2)`. The sign is the caller's convention: the rotation
    !| from `v1` to `v2` sweeps through the plane the two span, and `orientation_reference` says
    !| which way round that plane is positive. Concretely, the sign is that of
    !| `orientation_reference . (v2 - (v1 . v2) v1)` -- the reference measured against the part
    !| of `v2` perpendicular to `v1`, which is the direction the rotation actually moves in.
    !|
    !| In two dimensions this reduces to the familiar determinant when the reference is `v1`
    !| turned a quarter turn; there is simply no such canonical quarter turn in higher ones.
    !|
    !| `undefined_sign` reports the one case with no answer: the reference is orthogonal to the
    !| rotation, so it orients nothing. Parallel `v1`/`v2` are not that case -- the angle is
    !| then 0 or pi and the sign does not matter, so the unsigned angle is returned as is.
    !|
    !| Shared compute core: the single-vector angle, reused per field by
    !| [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]].
    pure subroutine clock_hand_angle_between_vectors_helper(v1, v2, n_dims, orientation_reference, &
                                                            signed_angle, undefined_sign)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second normalized vector in RAP space
        real(real64), dimension(n_dims), intent(in) :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(real64), intent(out) :: signed_angle
            !! Signed angle between vectors in radians [-pi, pi]
        logical(c_bool), intent(out) :: undefined_sign
            !! `.true.` when the reference orients nothing, so no sign could be given

        real(real64) :: dot_product, unsigned_angle, along_rotation, perpendicular
        integer(int32) :: i_dim

        undefined_sign = .false.

        dot_product = 0.0_real64
        do concurrent (i_dim = 1:n_dims) shared(v1, v2) reduce(+:dot_product)
            dot_product = dot_product + v1(i_dim)*v2(i_dim)
        end do

        ! numerical slack can push a normalized dot product just outside the domain of acos
        dot_product = clamp(dot_product, min_val=-1.0_real64, max_val=1.0_real64)
        unsigned_angle = acos(dot_product)

        ! the part of v2 perpendicular to v1 is the direction the rotation moves in; the
        ! reference's component along it is what gives the rotation a sign
        along_rotation = 0.0_real64
        perpendicular = 0.0_real64
        do concurrent (i_dim = 1:n_dims) shared(v1, v2, orientation_reference, dot_product) &
            reduce(+:along_rotation, perpendicular)
            along_rotation = along_rotation &
                             + orientation_reference(i_dim)*(v2(i_dim) - dot_product*v1(i_dim))
            perpendicular = perpendicular + (v2(i_dim) - dot_product*v1(i_dim))**2
        end do

        if (perpendicular .isclose. 0.0_real64) then
            ! v1 and v2 are (anti)parallel: the angle is 0 or pi and has no side to fall on
            signed_angle = unsigned_angle
            return
        end if

        if (along_rotation .isclose. 0.0_real64) then
            ! the reference is orthogonal to the rotation, so it orients nothing
            signed_angle = unsigned_angle
            undefined_sign = .true.
            return
        end if

        signed_angle = sign(1.0_real64, along_rotation)*unsigned_angle
    end subroutine clock_hand_angle_between_vectors_helper

    !> summary: Compute signed rotation angles between for shift vectors, so between their origin and target
    !| AUTHOR_VIVIAN_BASS
    !| Each selected field is angled by the rule of
    !| [[tox_relative_axis_plane_tools_impl(module):clock_hand_angle_between_vectors_impl(subroutine)]],
    !| with one `orientation_reference` shared by the whole batch. A single field whose rotation
    !| the reference fails to orient fails the call.
    pure subroutine clock_hand_angles_for_shift_vectors_impl(fields, n_dims, n_fields, &
                                                        fields_selection_mask, &
                                                        n_selected_fields, orientation_reference, &
                                                        signed_angles, ierr)
        integer(int32), intent(in) :: n_dims
            !! Dimension of each vector in RAP space
        integer(int32), intent(in) :: n_fields
            !! Number of vector pairs
        real(real64), dimension(n_dims, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical(c_bool), dimension(n_fields), intent(in) :: fields_selection_mask
            !! .true. for vector pairs where angle should be computed
        integer(int32), intent(in) :: n_selected_fields
            !! Count of .true. values in fields_selection_mask
        real(real64), dimension(n_dims), intent(in) :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(real64), dimension(n_selected_fields), intent(out) :: signed_angles
            !! Signed rotation angles between vector pairs in radians [-π, π]
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_field, result_idx
        logical(c_bool) :: undefined_sign

        call set_ok(ierr)

        result_idx = 1
        do i_field = 1, n_fields
            if (fields_selection_mask(i_field)) then
                call clock_hand_angle_between_vectors_helper(fields(:, 1, i_field), fields(:, 2, i_field), &
                                                             n_dims, orientation_reference, &
                                                             signed_angles(result_idx), undefined_sign)
                ! one reference orients every field, so a failure is the reference's and not
                ! this field's -- report it once and let the rest be computed
                if (undefined_sign) call set_err_once(ierr, ERR_INVALID_INPUT)
                result_idx = result_idx + 1
            end if
        end do
    end subroutine clock_hand_angles_for_shift_vectors_impl

    !> summary: Compute the fractional contribution of each axis to a RAP-projected and normalized vector
    !| AUTHOR_VIVIAN_BASS
    !| Shared utility: the shift-vector and expression-vector entry points below both drive it.
    pure subroutine compute_relative_axis_contributions_impl(vec, n_axes, contributions, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (length of vec and contributions)
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized vector (expression or shift)
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
    end subroutine compute_relative_axis_contributions_impl

    !> summary: Compute fractional contribution of each axis to a RAP-projected and normalized shift vector.
    !| AUTHOR_VIVIAN_BASS
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    pure subroutine relative_axes_changes_from_shift_vector_impl(vec, n_axes, contributions, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized shift vector
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call compute_relative_axis_contributions_impl(vec, n_axes, contributions, ierr)
    end subroutine relative_axes_changes_from_shift_vector_impl

    !> summary: Compute fractional contribution of each axis to a RAP-projected and normalized expression vector.
    !| AUTHOR_VIVIAN_BASS
    !| Wrapper for single RAP-projected expression vectors
    pure subroutine relative_axes_expression_from_expression_vector_impl(vec, n_axes, contributions, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected and normalized expression vector
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis (output), values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call compute_relative_axis_contributions_impl(vec, n_axes, contributions, ierr)
    end subroutine relative_axes_expression_from_expression_vector_impl


end module tox_relative_axis_plane_tools_impl
