!> The `omics_field_RAP_projection` cases: each selected field's shift, origin minus target, is
!| projected onto the RAP (subtracting the mean of its selected coordinates). The cases pin that
!| sign (an open question, see the shift_sign case), the defining property of the projection on
!| inexact data, which axes and which fields the masks pick and in what order, the diagonal and
!| in-plane special cases, each dimension on its own, the selected counts checked against their
!| masks, NaN and Inf rejected anywhere in `fields`, and magnitudes near `huge`.
!|
!| Most cases select two or four axes: the mean is then a division by a power of two, exact at
!| every optimization level, so the expected values are compared exactly. The one case on three
!| axes says why it needs a tolerance.
module mod_test_omics_field_RAP_projection
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_relative_axis_plane_tools
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_omics_field_RAP_projection() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(13))

        all_tests(1) = test_case("test_omics_field_RAP_projection_shift_sign", &
                                 test_omics_field_RAP_projection_shift_sign)
        all_tests(2) = test_case("test_omics_field_RAP_projection_all_selected", &
                                 test_omics_field_RAP_projection_all_selected)
        all_tests(3) = test_case("test_omics_field_RAP_projection_three_axes", &
                                 test_omics_field_RAP_projection_three_axes)
        all_tests(4) = test_case("test_omics_field_RAP_projection_columns_sum_to_zero", &
                                 test_omics_field_RAP_projection_columns_sum_to_zero)
        all_tests(5) = test_case("test_omics_field_RAP_projection_axis_selection", &
                                 test_omics_field_RAP_projection_axis_selection)
        all_tests(6) = test_case("test_omics_field_RAP_projection_field_selection", &
                                 test_omics_field_RAP_projection_field_selection)
        all_tests(7) = test_case("test_omics_field_RAP_projection_mixed_selection", &
                                 test_omics_field_RAP_projection_mixed_selection)
        all_tests(8) = test_case("test_omics_field_RAP_projection_single_axis", &
                                 test_omics_field_RAP_projection_single_axis)
        all_tests(9) = test_case("test_omics_field_RAP_projection_diagonal_and_in_plane", &
                                 test_omics_field_RAP_projection_diagonal_and_in_plane)
        all_tests(10) = test_case("test_omics_field_RAP_projection_dimensions", &
                                  test_omics_field_RAP_projection_dimensions)
        all_tests(11) = test_case("test_omics_field_RAP_projection_selected_count_mismatch", &
                                  test_omics_field_RAP_projection_selected_count_mismatch)
        all_tests(12) = test_case("test_omics_field_RAP_projection_rejects_nan_and_inf", &
                                  test_omics_field_RAP_projection_rejects_nan_and_inf)
        all_tests(13) = test_case("test_omics_field_RAP_projection_extreme_magnitudes", &
                                  test_omics_field_RAP_projection_extreme_magnitudes)
    end function get_all_tests_omics_field_RAP_projection

    !> The sign of the shift: the implementation projects origin minus target, as its own code
    !| comment says. The doc comment only says the origin comes first and the target second, and
    !| calls the result a projected shift, which usually means target minus origin.
    !| QUESTION (for FES): which sign is meant? This case pins the current one, origin - target.
    !|
    !| Origin [1, 2, 3, 4] and target [4, 3, 2, 1] give [-3, -1, 1, 3], whose mean is 0, so it
    !| projects to itself; target - origin would give [3, 1, -1, -3]. The second field swaps origin
    !| and target, and its projection flips sign with it.
    subroutine test_omics_field_RAP_projection_shift_sign()
        integer(int32) :: ierr
        real(real64), dimension(4, 2, 2) :: fields
        real(real64), dimension(4, 2) :: projections, expected
        logical(c_bool) :: fields_mask(2), axes_mask(4)

        fields(:, 1, 1) = [1d0, 2d0, 3d0, 4d0]
        fields(:, 2, 1) = [4d0, 3d0, 2d0, 1d0]
        fields(:, 1, 2) = [4d0, 3d0, 2d0, 1d0]
        fields(:, 2, 2) = [1d0, 2d0, 3d0, 4d0]
        fields_mask = [.true., .true.]
        axes_mask = [.true., .true., .true., .true.]
        expected(:, 1) = [-3d0, -1d0, 1d0, 3d0]
        expected(:, 2) = [3d0, 1d0, -1d0, -3d0]

        call omics_field_RAP_projection(fields, 4, 2, fields_mask, 2, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_shift_sign: ierr")
        call assert_equal_array_real(projections, expected, 8, 0d0, &
                                     "test_omics_field_RAP_projection_shift_sign: origin - target")
    end subroutine test_omics_field_RAP_projection_shift_sign

    !> Every axis and every field selected, on four axes, with shifts whose means are not 0: origin
    !| [5, 0, 0, 3] minus target [1, 2, 2, -1] is [4, -2, -2, 4], mean 1, projection
    !| [3, -3, -3, 3]; [0, 0, 0, 0] - [-3, 5, 0, 2] = [3, -5, 0, -2], mean -1, projection
    !| [4, -4, 1, -1].
    subroutine test_omics_field_RAP_projection_all_selected()
        integer(int32) :: ierr
        real(real64), dimension(4, 2, 2) :: fields
        real(real64), dimension(4, 2) :: projections, expected
        logical(c_bool) :: fields_mask(2), axes_mask(4)

        fields(:, 1, 1) = [5d0, 0d0, 0d0, 3d0]
        fields(:, 2, 1) = [1d0, 2d0, 2d0, -1d0]
        fields(:, 1, 2) = [0d0, 0d0, 0d0, 0d0]
        fields(:, 2, 2) = [-3d0, 5d0, 0d0, 2d0]
        fields_mask = [.true., .true.]
        axes_mask = [.true., .true., .true., .true.]
        expected(:, 1) = [3d0, -3d0, -3d0, 3d0]
        expected(:, 2) = [4d0, -4d0, 1d0, -1d0]

        call omics_field_RAP_projection(fields, 4, 2, fields_mask, 2, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_all_selected: ierr")
        call assert_equal_array_real(projections, expected, 8, 0d0, &
                                     "test_omics_field_RAP_projection_all_selected: projections")
    end subroutine test_omics_field_RAP_projection_all_selected

    !> The defining property on data with no exact result, ported from the Python suite's random
    !| check: on five of seven axes (a mean over a non-power of two) and three of four fields, each
    !| projected column sums to 0, and differs from its selected shift (origin - target) by one
    !| constant, the diagonal component. The two together determine the orthogonal projection
    !| onto the RAP uniquely.
    subroutine test_omics_field_RAP_projection_columns_sum_to_zero()
        integer(int32) :: ierr, i_field, i_axis
        real(real64), dimension(7, 2, 4) :: fields
        real(real64), dimension(5, 3) :: projections
        real(real64) :: diagonal_component
        logical(c_bool) :: fields_mask(4), axes_mask(7)
        integer(int32), parameter :: selected_axes(5) = [1, 3, 4, 6, 7], selected_fields(3) = [1, 2, 4]

        fields(:, 1, 1) = [0.37d0, -1.25d0, 2.6d0, 0.013d0, 7.9d0, -3.3d0, 0.71d0]
        fields(:, 2, 1) = [1.9d0, 0.4d0, -0.77d0, 5.1d0, 0.2d0, 0.9d0, -2.05d0]
        fields(:, 1, 2) = [-4.1d0, 3.7d0, 0.58d0, 1.11d0, -0.6d0, 2.2d0, 9.3d0]
        fields(:, 2, 2) = [0.3d0, 0.3d0, -6.4d0, 0.05d0, 1.7d0, -0.9d0, 4.4d0]
        fields(:, 1, 3) = 1000d0
        fields(:, 2, 3) = 0d0
        fields(:, 1, 4) = [2.5d0, -0.1d0, 0.33d0, -8.8d0, 6.6d0, 1.01d0, -0.47d0]
        fields(:, 2, 4) = [-1.2d0, 4.4d0, 0.9d0, -3.1d0, 0.07d0, 2.9d0, 1.6d0]
        fields_mask = [.true., .true., .false., .true.]
        axes_mask = [.true., .false., .true., .true., .false., .true., .true.]

        call omics_field_RAP_projection(fields, 7, 4, fields_mask, 3, axes_mask, 5, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_columns_sum_to_zero: ierr")
        ! A tolerance throughout: none of these values is exact, and the mean divides by 5.
        do i_field = 1, 3
            call assert_equal_real(sum(projections(:, i_field)), 0d0, 1d-13, &
                                   "test_omics_field_RAP_projection_columns_sum_to_zero: column sums to 0")
            diagonal_component = fields(selected_axes(1), 1, selected_fields(i_field)) &
                                 - fields(selected_axes(1), 2, selected_fields(i_field)) - projections(1, i_field)
            do i_axis = 2, 5
                call assert_equal_real(fields(selected_axes(i_axis), 1, selected_fields(i_field)) &
                                       - fields(selected_axes(i_axis), 2, selected_fields(i_field)) &
                                       - projections(i_axis, i_field), diagonal_component, 1d-13, &
                                       "test_omics_field_RAP_projection_columns_sum_to_zero: one constant per column")
            end do
        end do
    end subroutine test_omics_field_RAP_projection_columns_sum_to_zero

    !> The old suite's concrete example, on three axes: origin [1, -3, 1.1] minus target
    !| [3, 6, 2.2] is [-2, -9, -1.1], whose mean is -12.1/3 = -121/30, so the projection is
    !| [-60/30 + 121/30, -270/30 + 121/30, -33/30 + 121/30] = [61/30, -149/30, 88/30].
    subroutine test_omics_field_RAP_projection_three_axes()
        integer(int32) :: ierr
        real(real64), dimension(3, 2, 1) :: fields
        real(real64), dimension(3, 1) :: projections, expected
        logical(c_bool) :: fields_mask(1), axes_mask(3)

        fields(:, 1, 1) = [1d0, -3d0, 1.1d0]
        fields(:, 2, 1) = [3d0, 6d0, 2.2d0]
        fields_mask = [.true.]
        axes_mask = [.true., .true., .true.]
        expected(:, 1) = [61d0/30d0, -149d0/30d0, 88d0/30d0]

        call omics_field_RAP_projection(fields, 3, 1, fields_mask, 1, axes_mask, 3, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_three_axes: ierr")
        ! A few ulps: 1.1, 2.2 and 1/3 are not representable, and an optimizer may divide by 3 as a
        ! multiplication with the rounded 1/3 (ifx does).
        call assert_equal_array_real(projections, expected, 3, 1d-14, &
                                     "test_omics_field_RAP_projection_three_axes: projections")
    end subroutine test_omics_field_RAP_projection_three_axes

    !> Only the selected axes span the RAP, and they keep their order: origin
    !| [5, 7, 2, 4, 9, 8] minus target [4, -50, 0, 1, 50, 2] is [1, 57, 2, 3, -41, 6]; axes 1, 3, 4
    !| and 6 give [1, 2, 3, 6], mean 3, projection [-2, -1, 0, 3].
    subroutine test_omics_field_RAP_projection_axis_selection()
        integer(int32) :: ierr
        real(real64), dimension(6, 2, 1) :: fields
        real(real64), dimension(4, 1) :: projections, expected
        logical(c_bool) :: fields_mask(1), axes_mask(6)

        fields(:, 1, 1) = [5d0, 7d0, 2d0, 4d0, 9d0, 8d0]
        fields(:, 2, 1) = [4d0, -50d0, 0d0, 1d0, 50d0, 2d0]
        fields_mask = [.true.]
        axes_mask = [.true., .false., .true., .true., .false., .true.]
        expected(:, 1) = [-2d0, -1d0, 0d0, 3d0]

        call omics_field_RAP_projection(fields, 6, 1, fields_mask, 1, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_axis_selection: ierr")
        call assert_equal_array_real(projections, expected, 4, 0d0, &
                                     "test_omics_field_RAP_projection_axis_selection: projections")
    end subroutine test_omics_field_RAP_projection_axis_selection

    !> Only the selected fields are projected, into consecutive columns in their original order:
    !| fields 2 and 4 of four. [3, 7] - [1, 1] = [2, 6] projects to [-2, 2]; [12, 1] - [2, 1] =
    !| [10, 0] projects to [5, -5]. The unselected fields hold 1000s, which would show if used.
    subroutine test_omics_field_RAP_projection_field_selection()
        integer(int32) :: ierr
        real(real64), dimension(2, 2, 4) :: fields
        real(real64), dimension(2, 2) :: projections, expected
        logical(c_bool) :: fields_mask(4), axes_mask(2)

        fields(:, 1, 1) = [1000d0, 0d0]
        fields(:, 2, 1) = [0d0, 0d0]
        fields(:, 1, 2) = [3d0, 7d0]
        fields(:, 2, 2) = [1d0, 1d0]
        fields(:, 1, 3) = [0d0, 1000d0]
        fields(:, 2, 3) = [0d0, 0d0]
        fields(:, 1, 4) = [12d0, 1d0]
        fields(:, 2, 4) = [2d0, 1d0]
        fields_mask = [.false., .true., .false., .true.]
        axes_mask = [.true., .true.]
        expected(:, 1) = [-2d0, 2d0]
        expected(:, 2) = [5d0, -5d0]

        call omics_field_RAP_projection(fields, 2, 4, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_field_selection: ierr")
        call assert_equal_array_real(projections, expected, 4, 0d0, &
                                     "test_omics_field_RAP_projection_field_selection: projections")
    end subroutine test_omics_field_RAP_projection_field_selection

    !> Both masks partial at once: axes 1 and 3 of fields 1 and 3. [4, 9, 3] - [3, 0, 0] =
    !| [1, 9, 3] gives [1, 3], mean 2, projection [-1, 1]; [7, 8, 20] - [0, 0, 1] = [7, 8, 19]
    !| gives [7, 19], mean 13, projection [-6, 6].
    subroutine test_omics_field_RAP_projection_mixed_selection()
        integer(int32) :: ierr
        real(real64), dimension(3, 2, 3) :: fields
        real(real64), dimension(2, 2) :: projections, expected
        logical(c_bool) :: fields_mask(3), axes_mask(3)

        fields(:, 1, 1) = [4d0, 9d0, 3d0]
        fields(:, 2, 1) = [3d0, 0d0, 0d0]
        fields(:, 1, 2) = [100d0, 200d0, 300d0]
        fields(:, 2, 2) = [0d0, 0d0, 0d0]
        fields(:, 1, 3) = [7d0, 8d0, 20d0]
        fields(:, 2, 3) = [0d0, 0d0, 1d0]
        fields_mask = [.true., .false., .true.]
        axes_mask = [.true., .false., .true.]
        expected(:, 1) = [-1d0, 1d0]
        expected(:, 2) = [-6d0, 6d0]

        call omics_field_RAP_projection(fields, 3, 3, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_mixed_selection: ierr")
        call assert_equal_array_real(projections, expected, 4, 0d0, &
                                     "test_omics_field_RAP_projection_mixed_selection: projections")
    end subroutine test_omics_field_RAP_projection_mixed_selection

    !> On a single axis the RAP is the point 0: every coordinate is its own mean. First one axis of
    !| three, then the smallest valid call, every dimension 1.
    subroutine test_omics_field_RAP_projection_single_axis()
        integer(int32) :: ierr
        real(real64), dimension(3, 2, 2) :: fields
        real(real64), dimension(1, 2) :: projections, expected
        real(real64), dimension(1, 2, 1) :: single_field
        real(real64), dimension(1, 1) :: single_projection, single_expected
        logical(c_bool) :: fields_mask(2), axes_mask(3), single_mask(1)

        fields(:, 1, 1) = [1d0, 2d0, 3d0]
        fields(:, 2, 1) = [0d0, 7d0, 0d0]
        fields(:, 1, 2) = [4d0, -5d0, 6d0]
        fields(:, 2, 2) = [1d0, 1d0, 1d0]
        fields_mask = [.true., .true.]
        axes_mask = [.false., .true., .false.]
        projections = 99d0
        expected = 0d0

        call omics_field_RAP_projection(fields, 3, 2, fields_mask, 2, axes_mask, 1, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_single_axis: ierr")
        call assert_equal_array_real(projections, expected, 2, 0d0, &
                                     "test_omics_field_RAP_projection_single_axis: one axis of three")

        single_field(1, :, 1) = [7d0, 3d0]
        single_mask = [.true.]
        single_projection = 99d0
        single_expected = 0d0
        call omics_field_RAP_projection(single_field, 1, 1, single_mask, 1, single_mask, 1, single_projection, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_single_axis: ierr, 1x1")
        call assert_equal_array_real(single_projection, single_expected, 1, 0d0, &
                                     "test_omics_field_RAP_projection_single_axis: every dimension 1")
    end subroutine test_omics_field_RAP_projection_single_axis

    !> A field whose target equals its origin projects to zero, and so does one whose shift lies on
    !| the diagonal: [6, 7, 8, 9] - [1, 2, 3, 4] = [5, 5, 5, 5]. A shift already in the RAP,
    !| [2, 1, 0, 1] - [1, 1, 1, 1] = [1, 0, -1, 0], stays as it is.
    subroutine test_omics_field_RAP_projection_diagonal_and_in_plane()
        integer(int32) :: ierr
        real(real64), dimension(4, 2, 3) :: fields
        real(real64), dimension(4, 3) :: projections, expected
        logical(c_bool) :: fields_mask(3), axes_mask(4)

        fields(:, 1, 1) = [1d0, 2d0, 3d0, 4d0]
        fields(:, 2, 1) = [1d0, 2d0, 3d0, 4d0]
        fields(:, 1, 2) = [6d0, 7d0, 8d0, 9d0]
        fields(:, 2, 2) = [1d0, 2d0, 3d0, 4d0]
        fields(:, 1, 3) = [2d0, 1d0, 0d0, 1d0]
        fields(:, 2, 3) = [1d0, 1d0, 1d0, 1d0]
        fields_mask = [.true., .true., .true.]
        axes_mask = [.true., .true., .true., .true.]
        expected(:, 1) = 0d0
        expected(:, 2) = 0d0
        expected(:, 3) = [1d0, 0d0, -1d0, 0d0]

        call omics_field_RAP_projection(fields, 4, 3, fields_mask, 3, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_diagonal_and_in_plane: ierr")
        call assert_equal_array_real(projections(:, 1:2), expected(:, 1:2), 8, 0d0, &
                                     "test_omics_field_RAP_projection_diagonal_and_in_plane: no or diagonal shift to zero")
        call assert_equal_array_real(projections(:, 3), expected(:, 3), 4, 0d0, &
                                     "test_omics_field_RAP_projection_diagonal_and_in_plane: in-plane shift unchanged")
    end subroutine test_omics_field_RAP_projection_diagonal_and_in_plane

    !> Each dimension on its own: 0 is ERR_EMPTY_INPUT and -1 ERR_INVALID_INPUT, blamed on that
    !| argument. The valid call next to them projects the shift [1, 1] - [0, 0] to [0, 0].
    subroutine test_omics_field_RAP_projection_dimensions()
        integer(int32) :: ierr
        real(real64), dimension(2, 2, 2) :: fields
        real(real64), dimension(2, 2) :: projections, expected
        logical(c_bool) :: fields_mask(2), axes_mask(2)

        fields(:, 1, :) = 1d0
        fields(:, 2, :) = 0d0
        fields_mask = [.true., .true.]
        axes_mask = [.true., .true.]
        expected = 0d0

        call omics_field_RAP_projection(fields, 2, 2, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_dimensions: valid call")
        call assert_equal_array_real(projections, expected, 4, 0d0, "test_omics_field_RAP_projection_dimensions: valid call")

        call omics_field_RAP_projection(fields, 0, 2, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_field_RAP_projection_dimensions: n_axes = 0", 2)
        call omics_field_RAP_projection(fields, -1, 2, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_dimensions: n_axes = -1", 2)

        call omics_field_RAP_projection(fields, 2, 0, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_field_RAP_projection_dimensions: n_fields = 0", 3)
        call omics_field_RAP_projection(fields, 2, -1, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_dimensions: n_fields = -1", 3)

        call omics_field_RAP_projection(fields, 2, 2, fields_mask, 0, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_field_RAP_projection_dimensions: n_selected_fields = 0", 5)
        call omics_field_RAP_projection(fields, 2, 2, fields_mask, -1, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_dimensions: n_selected_fields = -1", 5)

        call omics_field_RAP_projection(fields, 2, 2, fields_mask, 2, axes_mask, 0, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_field_RAP_projection_dimensions: n_selected_axes = 0", 7)
        call omics_field_RAP_projection(fields, 2, 2, fields_mask, 2, axes_mask, -1, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_dimensions: n_selected_axes = -1", 7)
    end subroutine test_omics_field_RAP_projection_dimensions

    !> Each selected count must equal the `.true.` entries of its mask, here 2: one fewer or one
    !| more is ERR_INVALID_INPUT, blamed on the count. Too many would write past the selected
    !| columns or rows. `projections` is 3x3 so that every call has room for what it claims.
    subroutine test_omics_field_RAP_projection_selected_count_mismatch()
        integer(int32) :: ierr
        real(real64), dimension(3, 2, 3) :: fields
        real(real64), dimension(3, 3) :: projections
        logical(c_bool) :: fields_mask(3), axes_mask(3)

        fields = 1d0
        fields_mask = [.true., .false., .true.]
        axes_mask = [.true., .true., .false.]

        call omics_field_RAP_projection(fields, 3, 3, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_selected_count_mismatch: counts 2 and 2")

        call omics_field_RAP_projection(fields, 3, 3, fields_mask, 1, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_selected_count_mismatch: n_selected_fields = 1", 5)
        call omics_field_RAP_projection(fields, 3, 3, fields_mask, 3, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_selected_count_mismatch: n_selected_fields = 3", 5)

        call omics_field_RAP_projection(fields, 3, 3, fields_mask, 2, axes_mask, 1, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_selected_count_mismatch: n_selected_axes = 1", 7)
        call omics_field_RAP_projection(fields, 3, 3, fields_mask, 2, axes_mask, 3, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_field_RAP_projection_selected_count_mismatch: n_selected_axes = 3", 7)
    end subroutine test_omics_field_RAP_projection_selected_count_mismatch

    !> NaN, +Inf and -Inf are rejected with ERR_NAN_INF on `fields`, wherever they sit: in a
    !| selected origin, in a selected target, on an unselected axis, and in an unselected field.
    !| TOX takes no NaN or Inf.
    subroutine test_omics_field_RAP_projection_rejects_nan_and_inf()
        integer(int32) :: ierr, i_bad, i_position
        real(real64), dimension(3, 2, 3) :: fields
        real(real64), dimension(2, 2) :: projections
        logical(c_bool) :: fields_mask(3), axes_mask(3)
        real(real64) :: bad(3)
        integer(int32) :: positions(3, 4)
        character(len=*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        character(len=*), parameter :: position_names(4) = ["a selected origin    ", "a selected target    ", &
                                                            "an unselected axis   ", "an unselected field  "]

        fields_mask = [.true., .true., .false.]
        axes_mask = [.true., .true., .false.]
        bad(1) = ieee_value(1d0, ieee_quiet_nan)
        bad(2) = ieee_value(1d0, ieee_positive_inf)
        bad(3) = ieee_value(1d0, ieee_negative_inf)
        ! (axis, origin 1 or target 2, field) of each position
        positions(:, 1) = [1, 1, 1]
        positions(:, 2) = [1, 2, 1]
        positions(:, 3) = [3, 1, 1]
        positions(:, 4) = [1, 1, 3]

        do i_bad = 1, size(bad)
            do i_position = 1, size(positions, 2)
                fields = 1d0
                fields(positions(1, i_position), positions(2, i_position), positions(3, i_position)) = bad(i_bad)
                call omics_field_RAP_projection(fields, 3, 3, fields_mask, 2, axes_mask, 2, projections, ierr)
                call assert_err(ierr, ERR_NAN_INF, "test_omics_field_RAP_projection_rejects_nan_and_inf: " &
                                // trim(bad_names(i_bad)) // " in " // trim(position_names(i_position)), 1)
            end do
        end do
    end subroutine test_omics_field_RAP_projection_rejects_nan_and_inf

    !> Finite input near the top of the range gives a finite projection. Origin [huge, huge] minus
    !| target [0, 0] lies on the diagonal and projects to [0, 0]; [huge, -huge] - [0, 0] has mean 0
    !| and stays as it is.
    !|
    !| BUG: the mean is taken as sum/n, and huge + huge overflows to Inf, so the first shift
    !| projects to [-Inf, -Inf]: finite input, infinite output, which TOX must never write.
    !| Summing x_i/n instead cannot overflow where the mean itself does not.
    subroutine test_omics_field_RAP_projection_extreme_magnitudes()
        integer(int32) :: ierr
        real(real64), dimension(2, 2, 2) :: fields
        real(real64), dimension(2, 2) :: projections, expected
        logical(c_bool) :: fields_mask(2), axes_mask(2)

        fields(:, 1, 1) = [huge(1d0), huge(1d0)]
        fields(:, 1, 2) = [huge(1d0), -huge(1d0)]
        fields(:, 2, :) = 0d0
        fields_mask = [.true., .true.]
        axes_mask = [.true., .true.]
        expected(:, 1) = [0d0, 0d0]
        expected(:, 2) = [huge(1d0), -huge(1d0)]

        call omics_field_RAP_projection(fields, 2, 2, fields_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_field_RAP_projection_extreme_magnitudes: ierr")
        call assert_equal_array_real(projections(:, 1), expected(:, 1), 2, 0d0, &
                                     "test_omics_field_RAP_projection_extreme_magnitudes: [huge, huge] to zero")
        call assert_equal_array_real(projections(:, 2), expected(:, 2), 2, 0d0, &
                                     "test_omics_field_RAP_projection_extreme_magnitudes: [huge, -huge] unchanged")
    end subroutine test_omics_field_RAP_projection_extreme_magnitudes

end module mod_test_omics_field_RAP_projection
