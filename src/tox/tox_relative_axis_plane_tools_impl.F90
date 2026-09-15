#include <src/macros.h>

!> Relative axis planes (RAPs): planes through higher-dimensional gene expression space, and
!| what can be read off a vector once it is projected onto one.
!|
!| A RAP is picked by selecting axes (tissues) from the full expression space.
!| `omics_vector_RAP_projection` projects vectors onto it and `omics_field_RAP_projection` a whole
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
    use f42_math_impl, only: EPS
    M_IMPLICIT_NONE

contains

    !> summary: Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine omics_vector_RAP_projection_impl(vecs, n_axes, n_vecs, vecs_selection_mask, n_selected_vecs, axes_selection_mask, n_selected_axes, projections)
        integer(int32), intent(in) :: n_axes
            !! number of axes
        integer(int32), intent(in) :: n_vecs
            !! number of vectors
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
            !! number of fields
        real(real64), dimension(n_axes, 2, n_fields), intent(in) :: fields
            !! matrix with vector fields; each field holds two vectors, the origin first and the target second
        logical(c_bool), dimension(n_fields), intent(in) :: fields_selection_mask
            !! `.true.` for fields where projection is to be computed
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
            !! number of selected vectors
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

            ! calculate diagonal component to be subtracted from vectors for projection. Each
            ! coordinate is divided before it is added: a mean is never larger than its largest
            ! value, but the sum can overflow where the mean does not.
            diagonal_component = 0.0_real64
            do concurrent (i_axis = 1:n_selected_axes) shared(selected_vecs, i_vec, n_selected_axes) &
                reduce(+:diagonal_component)
                diagonal_component = diagonal_component + selected_vecs(i_axis, i_vec)/real(n_selected_axes, real64)
            end do

            ! transform vector to its projection
            do concurrent (i_axis = 1:n_selected_axes) shared(selected_vecs, i_vec, diagonal_component)
                selected_vecs(i_axis, i_vec) = selected_vecs(i_axis, i_vec) - diagonal_component
            end do
        end do
    end subroutine project_selected_vecs_onto_rap_helper


    !> summary: Compute the signed clock hand angle between two RAP-projected vectors.
    !| AUTHOR_VIVIAN_BASS
    !| The unsigned angle is the one between the two directions; their lengths do not matter, so
    !| the vectors need not be normalized. `orientation_reference` supplies the sign by saying
    !| which way round the plane the two vectors span counts as positive. Reports
    !| `ERR_DIVISION_BY_ZERO` when `v1` or `v2` is the zero vector, which has no direction, and
    !| `ERR_INVALID_INPUT` when the reference is orthogonal to the rotation and so orients nothing.
    pure subroutine clock_hand_angle_between_vectors_impl(v1, v2, n_dims, orientation_reference, &
                                                            signed_angle, ierr)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First vector in RAP space, of any length but zero
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second vector in RAP space, of any length but zero
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

        logical(c_bool) :: zero_vector, undefined_sign

        call set_ok(ierr)
        call clock_hand_angle_between_vectors_helper(v1, v2, n_dims, orientation_reference, &
                                                     signed_angle, zero_vector, undefined_sign)
        if (zero_vector) call set_err_once(ierr, ERR_DIVISION_BY_ZERO)
        ! not a bad argument on its own -- the reference only fails to orient *this* rotation,
        ! which no check on any single argument could have foreseen
        if (undefined_sign) call set_err_once(ierr, ERR_INVALID_INPUT)
    end subroutine clock_hand_angle_between_vectors_impl

    !> AUTHOR_VIVIAN_BASS
    !| Compute the signed clock hand angle between two RAP-projected vectors.
    !|
    !| The unsigned angle is `atan2(|v1| |p|, v1 . v2)`, where `p = v2 - ((v1 . v2)/(v1 . v1)) v1`
    !| is the part of `v2` perpendicular to `v1`: the arctangent of the angle's sine and cosine
    !| parts, both scaled by the same lengths. Unlike `acos(v1 . v2)` for unit vectors, it needs
    !| no normalized input and keeps its accuracy near 0 and PI, where acos loses half its digits.
    !| The sign is the caller's convention: the rotation from `v1` to `v2` sweeps through the plane
    !| the two span, and `orientation_reference` says which way round that plane is positive.
    !| Concretely, the sign is that of `orientation_reference . p` -- the reference measured
    !| against the direction the rotation actually moves in.
    !|
    !| In two dimensions this reduces to the familiar determinant when the reference is `v1`
    !| turned a quarter turn; there is simply no such canonical quarter turn in higher ones.
    !|
    !| `zero_vector` reports a zero `v1` or `v2`, which has no direction to turn from or to.
    !| `undefined_sign` reports the one other case with no answer: the reference is orthogonal
    !| to the rotation, so it orients nothing. Parallel `v1`/`v2` are not that case -- the angle
    !| is then 0 or pi and the sign does not matter, so the unsigned angle is returned as is.
    !|
    !| Shared compute core: the single-vector angle, reused per field by
    !| [[tox_relative_axis_plane_tools_impl(module):clock_hand_angles_for_shift_vectors_impl(subroutine)]].
    pure subroutine clock_hand_angle_between_vectors_helper(v1, v2, n_dims, orientation_reference, &
                                                            signed_angle, zero_vector, undefined_sign)
        integer(int32), intent(in) :: n_dims
            !! Dimension of both vectors
        real(real64), dimension(n_dims), intent(in) :: v1
            !! First vector in RAP space, of any length but zero
        real(real64), dimension(n_dims), intent(in) :: v2
            !! Second vector in RAP space, of any length but zero
        real(real64), dimension(n_dims), intent(in) :: orientation_reference
            !! Orients the plane the rotation happens in, so the angle can carry a sign. A
            !! rotation from one vector to another has no inherent direction above two
            !! dimensions -- and in RAP space not even in two, since the axes are tissues or
            !! factors and carry no handedness -- so the caller states which way round counts
            !! as positive. The sign is that of this vector's component along the rotation.
        real(real64), intent(out) :: signed_angle
            !! Signed angle between vectors in radians [-pi, pi]
        logical(c_bool), intent(out) :: zero_vector
            !! `.true.` when `v1` or `v2` is the zero vector, so there is no angle
        logical(c_bool), intent(out) :: undefined_sign
            !! `.true.` when the reference orients nothing, so no sign could be given

        real(real64) :: largest_1, largest_2, largest_reference, dot_12, squared_1, squared_2
        real(real64) :: squared_perpendicular, squared_reference, along_rotation, projection_factor
        real(real64) :: perpendicular_component, unsigned_angle, rounding_level
        integer(int32) :: i_dim, exponent_1, exponent_2, exponent_reference

        signed_angle = 0.0_real64
        zero_vector = .false.
        undefined_sign = .false.

        largest_1 = 0.0_real64
        largest_2 = 0.0_real64
        largest_reference = 0.0_real64
        do concurrent (i_dim = 1:n_dims) shared(v1, v2, orientation_reference) &
            reduce(max:largest_1, largest_2, largest_reference)
            largest_1 = max(largest_1, abs(v1(i_dim)))
            largest_2 = max(largest_2, abs(v2(i_dim)))
            largest_reference = max(largest_reference, abs(orientation_reference(i_dim)))
        end do

        ! Only an exactly zero vector is rejected: it alone has no direction, and any other,
        ! however short, still has one -- a tolerance would reject valid vectors by their length.
        ! A largest magnitude is never negative, so `<=` is that exact test without comparing
        ! reals for equality.
        if (largest_1 <= 0.0_real64 .or. largest_2 <= 0.0_real64) then
            zero_vector = .true.
            return
        end if

        ! The angle does not depend on the vectors' lengths, so each is scaled by a power of two
        ! near its largest component. That is exact, and it keeps the products below from
        ! overflowing or underflowing whatever the magnitudes. A zero reference stays zero.
        exponent_1 = exponent(largest_1)
        exponent_2 = exponent(largest_2)
        exponent_reference = exponent(largest_reference)

        dot_12 = 0.0_real64
        squared_1 = 0.0_real64
        squared_2 = 0.0_real64
        do concurrent (i_dim = 1:n_dims) shared(v1, v2, exponent_1, exponent_2) &
            reduce(+:dot_12, squared_1, squared_2)
            dot_12 = dot_12 + scale(v1(i_dim), -exponent_1)*scale(v2(i_dim), -exponent_2)
            squared_1 = squared_1 + scale(v1(i_dim), -exponent_1)**2
            squared_2 = squared_2 + scale(v2(i_dim), -exponent_2)**2
        end do

        ! the part of v2 perpendicular to v1 is the direction the rotation moves in; the
        ! reference's component along it is what gives the rotation a sign
        projection_factor = dot_12/squared_1
        squared_perpendicular = 0.0_real64
        squared_reference = 0.0_real64
        along_rotation = 0.0_real64
        do concurrent (i_dim = 1:n_dims) local(perpendicular_component) &
            shared(v1, v2, orientation_reference, exponent_1, exponent_2, exponent_reference, projection_factor) &
            reduce(+:squared_perpendicular, squared_reference, along_rotation)
            perpendicular_component = scale(v2(i_dim), -exponent_2) - projection_factor*scale(v1(i_dim), -exponent_1)
            squared_perpendicular = squared_perpendicular + perpendicular_component**2
            squared_reference = squared_reference + scale(orientation_reference(i_dim), -exponent_reference)**2
            along_rotation = along_rotation + scale(orientation_reference(i_dim), -exponent_reference)*perpendicular_component
        end do

        unsigned_angle = atan2(sqrt(squared_1*squared_perpendicular), dot_12)

        ! Both checks below compare with what rounding can leave behind, relative to the scaled
        ! lengths: a sum of n_dims rounded products is off by up to about n_dims*EPS of them.
        ! is_close's absolute floor would instead decide by the vectors' lengths -- every turn
        ! below 1e-6 rad counted as parallel, and a reference shorter than 1e-12 as orienting
        ! nothing.
        rounding_level = real(n_dims, real64)*EPS*sqrt(squared_2)

        if (sqrt(squared_perpendicular) <= rounding_level) then
            ! v1 and v2 are (anti)parallel, as far as rounding can tell: the angle is 0 or pi and
            ! has no side to fall on
            signed_angle = unsigned_angle
            return
        end if

        if (abs(along_rotation) <= rounding_level*sqrt(squared_reference)) then
            ! the reference is orthogonal to the rotation, so it orients nothing
            signed_angle = unsigned_angle
            undefined_sign = .true.
            return
        end if

        signed_angle = sign(1.0_real64, along_rotation)*unsigned_angle
    end subroutine clock_hand_angle_between_vectors_helper

    !> summary: Compute the signed clock hand angle of every selected field, from its origin to its target
    !| AUTHOR_VIVIAN_BASS
    !| Each selected field is angled by the rule of
    !| [[tox_relative_axis_plane_tools_impl(module):clock_hand_angle_between_vectors_impl(subroutine)]],
    !| with one `orientation_reference` shared by the whole batch. A single selected field with a
    !| zero origin or target fails the call with `ERR_DIVISION_BY_ZERO`, and one whose rotation
    !| the reference fails to orient with `ERR_INVALID_INPUT`.
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
            !! Signed rotation angles between vector pairs in radians [-pi, pi]
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_field, result_idx
        logical(c_bool) :: zero_vector, undefined_sign

        call set_ok(ierr)

        result_idx = 1
        do i_field = 1, n_fields
            if (fields_selection_mask(i_field)) then
                call clock_hand_angle_between_vectors_helper(fields(:, 1, i_field), fields(:, 2, i_field), &
                                                             n_dims, orientation_reference, &
                                                             signed_angles(result_idx), zero_vector, undefined_sign)
                if (zero_vector) call set_err_once(ierr, ERR_DIVISION_BY_ZERO)
                ! one reference orients every field, so a failure is the reference's and not
                ! this field's -- report it once and let the rest be computed
                if (undefined_sign) call set_err_once(ierr, ERR_INVALID_INPUT)
                result_idx = result_idx + 1
            end if
        end do
    end subroutine clock_hand_angles_for_shift_vectors_impl

    !> summary: Compute the fractional contribution of each axis to a RAP-projected vector
    !| AUTHOR_VIVIAN_BASS
    !| Shared utility: the shift-vector and expression-vector entry points below both drive it.
    !| The shares depend on the vector's direction alone, so it need not be normalized; the zero
    !| vector, which has no direction, is `ERR_DIVISION_BY_ZERO`.
    pure subroutine compute_relative_axis_contributions_impl(vec, n_axes, contributions, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes (length of vec and contributions)
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected vector (expression or shift), of any length but zero
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis, values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        real(real64) :: largest, total_abs
        integer(int32) :: i_axis, scale_exponent

        ! Error handling
        call set_ok(ierr)

        largest = 0.0_real64
        do concurrent (i_axis = 1:n_axes) shared(vec) reduce(max:largest)
            largest = max(largest, abs(vec(i_axis)))
        end do

        ! Only an exactly zero vector is rejected: it alone has no direction, and any other,
        ! however short, still divides into shares -- a tolerance would reject valid vectors by
        ! their length. A largest magnitude is never negative, so `<=` is that exact test without
        ! comparing reals for equality.
        if (largest <= 0.0_real64) then
            contributions = 0.0_real64
            call set_err_once(ierr, ERR_DIVISION_BY_ZERO, arg_pos=1_int32)
            return
        end if

        ! Every component is scaled by the same power of two near the largest one. That is exact,
        ! so the shares are the unscaled vector's, while the sum stays at most n_axes and cannot
        ! overflow however huge the components are.
        scale_exponent = exponent(largest)
        total_abs = 0.0_real64
        do concurrent (i_axis = 1:n_axes) shared(vec, scale_exponent) reduce(+:total_abs)
            total_abs = total_abs + scale(abs(vec(i_axis)), -scale_exponent)
        end do

        do concurrent (i_axis = 1:n_axes) shared(contributions, vec, scale_exponent, total_abs)
            contributions(i_axis) = scale(abs(vec(i_axis)), -scale_exponent)/total_abs
        end do
    end subroutine compute_relative_axis_contributions_impl

    !> summary: Compute fractional contribution of each axis to a RAP-projected shift vector.
    !| AUTHOR_VIVIAN_BASS
    !| Wrapper for shift vectors (e.g. difference between two RAP-projected vectors)
    pure subroutine relative_axes_changes_from_shift_vector_impl(vec, n_axes, contributions, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected shift vector, of any length but zero
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis, values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call compute_relative_axis_contributions_impl(vec, n_axes, contributions, ierr)
    end subroutine relative_axes_changes_from_shift_vector_impl

    !> summary: Compute fractional contribution of each axis to a RAP-projected expression vector.
    !| AUTHOR_VIVIAN_BASS
    !| Wrapper for single RAP-projected expression vectors
    pure subroutine relative_axes_expression_from_expression_vector_impl(vec, n_axes, contributions, ierr)
        integer(int32), intent(in) :: n_axes
            !! Number of axes
        real(real64), dimension(n_axes), intent(in) :: vec
            !! RAP-projected expression vector, of any length but zero
        real(real64), dimension(n_axes), intent(out) :: contributions
            !! Fractional contribution of each axis, values in [0,1], sum to 1
        integer(int32), intent(out) :: ierr
            !! Error code

        call compute_relative_axis_contributions_impl(vec, n_axes, contributions, ierr)
    end subroutine relative_axes_expression_from_expression_vector_impl


end module tox_relative_axis_plane_tools_impl
