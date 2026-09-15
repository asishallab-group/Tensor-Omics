!> The `omics_vector_RAP_projection` cases: hand-derived projections onto the RAP (subtracting the
!| mean of the selected coordinates), which axes and which vectors the masks pick and in what
!| order, the diagonal and in-plane special cases, each dimension on its own, the selected counts
!| checked against their masks, NaN and Inf rejected anywhere in `vecs`, and magnitudes near
!| `huge`.
!|
!| Most cases select two or four axes: the mean is then a division by a power of two, exact at
!| every optimization level, so the expected values are compared exactly. The one case on three
!| axes says why it needs a tolerance.
module mod_test_omics_vector_RAP_projection
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
    function get_all_tests_omics_vector_RAP_projection() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(12))

        all_tests(1) = test_case("test_omics_vector_RAP_projection_all_selected", &
                                 test_omics_vector_RAP_projection_all_selected)
        all_tests(2) = test_case("test_omics_vector_RAP_projection_three_axes", &
                                 test_omics_vector_RAP_projection_three_axes)
        all_tests(3) = test_case("test_omics_vector_RAP_projection_columns_sum_to_zero", &
                                 test_omics_vector_RAP_projection_columns_sum_to_zero)
        all_tests(4) = test_case("test_omics_vector_RAP_projection_axis_selection", &
                                 test_omics_vector_RAP_projection_axis_selection)
        all_tests(5) = test_case("test_omics_vector_RAP_projection_vector_selection", &
                                 test_omics_vector_RAP_projection_vector_selection)
        all_tests(6) = test_case("test_omics_vector_RAP_projection_mixed_selection", &
                                 test_omics_vector_RAP_projection_mixed_selection)
        all_tests(7) = test_case("test_omics_vector_RAP_projection_single_axis", &
                                 test_omics_vector_RAP_projection_single_axis)
        all_tests(8) = test_case("test_omics_vector_RAP_projection_diagonal_and_in_plane", &
                                 test_omics_vector_RAP_projection_diagonal_and_in_plane)
        all_tests(9) = test_case("test_omics_vector_RAP_projection_dimensions", &
                                 test_omics_vector_RAP_projection_dimensions)
        all_tests(10) = test_case("test_omics_vector_RAP_projection_selected_count_mismatch", &
                                  test_omics_vector_RAP_projection_selected_count_mismatch)
        all_tests(11) = test_case("test_omics_vector_RAP_projection_rejects_nan_and_inf", &
                                  test_omics_vector_RAP_projection_rejects_nan_and_inf)
        all_tests(12) = test_case("test_omics_vector_RAP_projection_extreme_magnitudes", &
                                  test_omics_vector_RAP_projection_extreme_magnitudes)
    end function get_all_tests_omics_vector_RAP_projection

    !> Every axis and every vector selected, on four axes. The RAP is the hyperplane orthogonal to
    !| (1, 1, 1, 1), so the projection subtracts the mean of the four coordinates from each:
    !| [1, 2, 3, 4] has mean 2.5, [4, 3, 2, 1] too, [-3, 5, 0, 2] has mean 1.
    subroutine test_omics_vector_RAP_projection_all_selected()
        integer(int32) :: ierr
        real(real64), dimension(4, 3) :: vecs, projections, expected
        logical(c_bool) :: vecs_mask(3), axes_mask(4)

        vecs(:, 1) = [1d0, 2d0, 3d0, 4d0]
        vecs(:, 2) = [4d0, 3d0, 2d0, 1d0]
        vecs(:, 3) = [-3d0, 5d0, 0d0, 2d0]
        vecs_mask = [.true., .true., .true.]
        axes_mask = [.true., .true., .true., .true.]
        expected(:, 1) = [-1.5d0, -0.5d0, 0.5d0, 1.5d0]
        expected(:, 2) = [1.5d0, 0.5d0, -0.5d0, -1.5d0]
        expected(:, 3) = [-4d0, 4d0, -1d0, 1d0]

        call omics_vector_RAP_projection(vecs, 4, 3, vecs_mask, 3, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_all_selected: ierr")
        call assert_equal_array_real(projections, expected, 12, 0d0, &
                                     "test_omics_vector_RAP_projection_all_selected: projections")
    end subroutine test_omics_vector_RAP_projection_all_selected

    !> Three axes, the case the old suite used throughout: [1, 2, 3], [4, 5, 6] and [7, 8, 9] have
    !| means 2, 5 and 8 and all project to [-1, 0, 1]; [1, 0, 0] has mean 1/3 and projects to
    !| [2/3, -1/3, -1/3].
    subroutine test_omics_vector_RAP_projection_three_axes()
        integer(int32) :: ierr
        real(real64), dimension(3, 4) :: vecs, projections, expected
        logical(c_bool) :: vecs_mask(4), axes_mask(3)

        vecs(:, 1) = [1d0, 2d0, 3d0]
        vecs(:, 2) = [4d0, 5d0, 6d0]
        vecs(:, 3) = [7d0, 8d0, 9d0]
        vecs(:, 4) = [1d0, 0d0, 0d0]
        vecs_mask = [.true., .true., .true., .true.]
        axes_mask = [.true., .true., .true.]
        expected(:, 1) = [-1d0, 0d0, 1d0]
        expected(:, 2) = [-1d0, 0d0, 1d0]
        expected(:, 3) = [-1d0, 0d0, 1d0]
        expected(:, 4) = [2d0/3d0, -1d0/3d0, -1d0/3d0]

        call omics_vector_RAP_projection(vecs, 3, 4, vecs_mask, 4, axes_mask, 3, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_three_axes: ierr")
        ! A few ulps: 1/3 is not representable, and an optimizer may divide by 3 as a multiplication
        ! with the rounded 1/3 (ifx does), so even the integer means may miss by an ulp.
        call assert_equal_array_real(projections, expected, 12, 1d-14, &
                                     "test_omics_vector_RAP_projection_three_axes: projections")
    end subroutine test_omics_vector_RAP_projection_three_axes

    !> The defining property on data with no exact result, ported from the Python suite's random
    !| check: on five of seven axes (a mean over a non-power of two) and three of four vectors,
    !| each projected column sums to 0, and differs from its selected input by one constant, the
    !| diagonal component. The two together determine the orthogonal projection onto the RAP
    !| uniquely.
    subroutine test_omics_vector_RAP_projection_columns_sum_to_zero()
        integer(int32) :: ierr, i_vec, i_axis
        real(real64), dimension(7, 4) :: vecs
        real(real64), dimension(5, 3) :: projections
        real(real64) :: diagonal_component
        logical(c_bool) :: vecs_mask(4), axes_mask(7)
        integer(int32), parameter :: selected_axes(5) = [1, 3, 4, 6, 7], selected_vecs(3) = [1, 2, 4]

        vecs(:, 1) = [0.37d0, -1.25d0, 2.6d0, 0.013d0, 7.9d0, -3.3d0, 0.71d0]
        vecs(:, 2) = [-4.1d0, 3.7d0, 0.58d0, 1.11d0, -0.6d0, 2.2d0, 9.3d0]
        vecs(:, 3) = 1000d0
        vecs(:, 4) = [2.5d0, -0.1d0, 0.33d0, -8.8d0, 6.6d0, 1.01d0, -0.47d0]
        vecs_mask = [.true., .true., .false., .true.]
        axes_mask = [.true., .false., .true., .true., .false., .true., .true.]

        call omics_vector_RAP_projection(vecs, 7, 4, vecs_mask, 3, axes_mask, 5, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_columns_sum_to_zero: ierr")
        ! A tolerance throughout: none of these values is exact, and the mean divides by 5.
        do i_vec = 1, 3
            call assert_equal_real(sum(projections(:, i_vec)), 0d0, 1d-13, &
                                   "test_omics_vector_RAP_projection_columns_sum_to_zero: column sums to 0")
            diagonal_component = vecs(selected_axes(1), selected_vecs(i_vec)) - projections(1, i_vec)
            do i_axis = 2, 5
                call assert_equal_real(vecs(selected_axes(i_axis), selected_vecs(i_vec)) - projections(i_axis, i_vec), &
                                       diagonal_component, 1d-13, &
                                       "test_omics_vector_RAP_projection_columns_sum_to_zero: one constant per column")
            end do
        end do
    end subroutine test_omics_vector_RAP_projection_columns_sum_to_zero

    !> Only the selected axes span the RAP, and they keep their order: of [1, 100, 2, 3, -100, 6],
    !| axes 1, 3, 4 and 6 are selected, [1, 2, 3, 6] has mean 3, so the projection is
    !| [-2, -1, 0, 3]. Either unselected axis leaking in would shift the mean.
    subroutine test_omics_vector_RAP_projection_axis_selection()
        integer(int32) :: ierr
        real(real64), dimension(6, 1) :: vecs
        real(real64), dimension(4, 1) :: projections, expected
        logical(c_bool) :: vecs_mask(1), axes_mask(6)

        vecs(:, 1) = [1d0, 100d0, 2d0, 3d0, -100d0, 6d0]
        vecs_mask = [.true.]
        axes_mask = [.true., .false., .true., .true., .false., .true.]
        expected(:, 1) = [-2d0, -1d0, 0d0, 3d0]

        call omics_vector_RAP_projection(vecs, 6, 1, vecs_mask, 1, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_axis_selection: ierr")
        call assert_equal_array_real(projections, expected, 4, 0d0, &
                                     "test_omics_vector_RAP_projection_axis_selection: projections")
    end subroutine test_omics_vector_RAP_projection_axis_selection

    !> Only the selected vectors are projected, into consecutive columns in their original order:
    !| vectors 2 and 4 of four. [2, 6] has mean 4 and projects to [-2, 2], [10, 0] has mean 5 and
    !| projects to [5, -5]; the unselected vectors hold 1000s, which would show if they were used.
    subroutine test_omics_vector_RAP_projection_vector_selection()
        integer(int32) :: ierr
        real(real64), dimension(2, 4) :: vecs
        real(real64), dimension(2, 2) :: projections, expected
        logical(c_bool) :: vecs_mask(4), axes_mask(2)

        vecs(:, 1) = [1000d0, 0d0]
        vecs(:, 2) = [2d0, 6d0]
        vecs(:, 3) = [0d0, 1000d0]
        vecs(:, 4) = [10d0, 0d0]
        vecs_mask = [.false., .true., .false., .true.]
        axes_mask = [.true., .true.]
        expected(:, 1) = [-2d0, 2d0]
        expected(:, 2) = [5d0, -5d0]

        call omics_vector_RAP_projection(vecs, 2, 4, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_vector_selection: ierr")
        call assert_equal_array_real(projections, expected, 4, 0d0, &
                                     "test_omics_vector_RAP_projection_vector_selection: projections")
    end subroutine test_omics_vector_RAP_projection_vector_selection

    !> Both masks partial at once: axes 1 and 3 of vectors 1 and 3. [1, 3] has mean 2 and projects
    !| to [-1, 1]; [7, 19] has mean 13 and projects to [-6, 6].
    subroutine test_omics_vector_RAP_projection_mixed_selection()
        integer(int32) :: ierr
        real(real64), dimension(3, 3) :: vecs
        real(real64), dimension(2, 2) :: projections, expected
        logical(c_bool) :: vecs_mask(3), axes_mask(3)

        vecs(:, 1) = [1d0, 2d0, 3d0]
        vecs(:, 2) = [4d0, 5d0, 6d0]
        vecs(:, 3) = [7d0, 8d0, 19d0]
        vecs_mask = [.true., .false., .true.]
        axes_mask = [.true., .false., .true.]
        expected(:, 1) = [-1d0, 1d0]
        expected(:, 2) = [-6d0, 6d0]

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_mixed_selection: ierr")
        call assert_equal_array_real(projections, expected, 4, 0d0, &
                                     "test_omics_vector_RAP_projection_mixed_selection: projections")
    end subroutine test_omics_vector_RAP_projection_mixed_selection

    !> On a single axis the RAP is the point 0: every coordinate is its own mean. First one axis of
    !| three, then the smallest valid call, every dimension 1.
    subroutine test_omics_vector_RAP_projection_single_axis()
        integer(int32) :: ierr
        real(real64), dimension(3, 2) :: vecs
        real(real64), dimension(1, 2) :: projections, expected
        real(real64), dimension(1, 1) :: single_vec, single_projection, single_expected
        logical(c_bool) :: vecs_mask(2), axes_mask(3), single_mask(1)

        vecs(:, 1) = [1d0, 2d0, 3d0]
        vecs(:, 2) = [4d0, -5d0, 6d0]
        vecs_mask = [.true., .true.]
        axes_mask = [.false., .true., .false.]
        projections = 99d0
        expected = 0d0

        call omics_vector_RAP_projection(vecs, 3, 2, vecs_mask, 2, axes_mask, 1, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_single_axis: ierr")
        call assert_equal_array_real(projections, expected, 2, 0d0, &
                                     "test_omics_vector_RAP_projection_single_axis: one axis of three")

        single_vec = 7d0
        single_mask = [.true.]
        single_projection = 99d0
        single_expected = 0d0
        call omics_vector_RAP_projection(single_vec, 1, 1, single_mask, 1, single_mask, 1, single_projection, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_single_axis: ierr, 1x1")
        call assert_equal_array_real(single_projection, single_expected, 1, 0d0, &
                                     "test_omics_vector_RAP_projection_single_axis: every dimension 1")
    end subroutine test_omics_vector_RAP_projection_single_axis

    !> A vector on the space diagonal projects to zero, and one already in the RAP (coordinates
    !| summing to zero) stays as it is, since the projection is idempotent.
    subroutine test_omics_vector_RAP_projection_diagonal_and_in_plane()
        integer(int32) :: ierr
        real(real64), dimension(4, 3) :: vecs, projections, expected
        logical(c_bool) :: vecs_mask(3), axes_mask(4)

        vecs(:, 1) = [5d0, 5d0, 5d0, 5d0]
        vecs(:, 2) = [1d0, 0d0, -1d0, 0d0]
        vecs(:, 3) = [0.5d0, -0.25d0, 0d0, -0.25d0]
        vecs_mask = [.true., .true., .true.]
        axes_mask = [.true., .true., .true., .true.]
        expected(:, 1) = 0d0
        expected(:, 2:3) = vecs(:, 2:3)

        call omics_vector_RAP_projection(vecs, 4, 3, vecs_mask, 3, axes_mask, 4, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_diagonal_and_in_plane: ierr")
        call assert_equal_array_real(projections(:, 1), expected(:, 1), 4, 0d0, &
                                     "test_omics_vector_RAP_projection_diagonal_and_in_plane: diagonal to zero")
        call assert_equal_array_real(projections(:, 2:3), expected(:, 2:3), 8, 0d0, &
                                     "test_omics_vector_RAP_projection_diagonal_and_in_plane: in-plane unchanged")
    end subroutine test_omics_vector_RAP_projection_diagonal_and_in_plane

    !> Each dimension on its own: 0 is ERR_EMPTY_INPUT and -1 ERR_INVALID_INPUT, blamed on that
    !| argument. The valid call next to them projects [1, 1] to [0, 0].
    subroutine test_omics_vector_RAP_projection_dimensions()
        integer(int32) :: ierr
        real(real64), dimension(2, 2) :: vecs, projections, expected
        logical(c_bool) :: vecs_mask(2), axes_mask(2)

        vecs = 1d0
        vecs_mask = [.true., .true.]
        axes_mask = [.true., .true.]
        expected = 0d0

        call omics_vector_RAP_projection(vecs, 2, 2, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_dimensions: valid call")
        call assert_equal_array_real(projections, expected, 4, 0d0, "test_omics_vector_RAP_projection_dimensions: valid call")

        call omics_vector_RAP_projection(vecs, 0, 2, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_vector_RAP_projection_dimensions: n_axes = 0", 2)
        call omics_vector_RAP_projection(vecs, -1, 2, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_dimensions: n_axes = -1", 2)

        call omics_vector_RAP_projection(vecs, 2, 0, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_vector_RAP_projection_dimensions: n_vecs = 0", 3)
        call omics_vector_RAP_projection(vecs, 2, -1, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_dimensions: n_vecs = -1", 3)

        call omics_vector_RAP_projection(vecs, 2, 2, vecs_mask, 0, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_vector_RAP_projection_dimensions: n_selected_vecs = 0", 5)
        call omics_vector_RAP_projection(vecs, 2, 2, vecs_mask, -1, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_dimensions: n_selected_vecs = -1", 5)

        call omics_vector_RAP_projection(vecs, 2, 2, vecs_mask, 2, axes_mask, 0, projections, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_omics_vector_RAP_projection_dimensions: n_selected_axes = 0", 7)
        call omics_vector_RAP_projection(vecs, 2, 2, vecs_mask, 2, axes_mask, -1, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_dimensions: n_selected_axes = -1", 7)
    end subroutine test_omics_vector_RAP_projection_dimensions

    !> Each selected count must equal the `.true.` entries of its mask, here 2: one fewer or one
    !| more is ERR_INVALID_INPUT, blamed on the count. Too many would write past the selected
    !| columns or rows. `projections` is 3x3 so that every call has room for what it claims.
    subroutine test_omics_vector_RAP_projection_selected_count_mismatch()
        integer(int32) :: ierr
        real(real64), dimension(3, 3) :: vecs, projections
        logical(c_bool) :: vecs_mask(3), axes_mask(3)

        vecs = 1d0
        vecs_mask = [.true., .false., .true.]
        axes_mask = [.true., .true., .false.]

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_selected_count_mismatch: counts 2 and 2")

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 1, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_selected_count_mismatch: n_selected_vecs = 1", 5)
        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 3, axes_mask, 2, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_selected_count_mismatch: n_selected_vecs = 3", 5)

        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 2, axes_mask, 1, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_selected_count_mismatch: n_selected_axes = 1", 7)
        call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 2, axes_mask, 3, projections, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_omics_vector_RAP_projection_selected_count_mismatch: n_selected_axes = 3", 7)
    end subroutine test_omics_vector_RAP_projection_selected_count_mismatch

    !> NaN, +Inf and -Inf are rejected with ERR_NAN_INF on `vecs`, wherever they sit: in a selected
    !| coordinate, on an unselected axis, and in an unselected vector. TOX takes no NaN or Inf.
    subroutine test_omics_vector_RAP_projection_rejects_nan_and_inf()
        integer(int32) :: ierr, i_bad, i_position
        real(real64), dimension(3, 3) :: vecs
        real(real64), dimension(2, 2) :: projections
        logical(c_bool) :: vecs_mask(3), axes_mask(3)
        real(real64) :: bad(3)
        integer(int32) :: positions(2, 3)
        character(len=*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        character(len=*), parameter :: position_names(3) = ["a selected coordinate", "an unselected axis   ", &
                                                            "an unselected vector "]

        vecs_mask = [.true., .true., .false.]
        axes_mask = [.true., .true., .false.]
        bad(1) = ieee_value(1d0, ieee_quiet_nan)
        bad(2) = ieee_value(1d0, ieee_positive_inf)
        bad(3) = ieee_value(1d0, ieee_negative_inf)
        ! (axis, vector) of each position
        positions(:, 1) = [1, 1]
        positions(:, 2) = [3, 1]
        positions(:, 3) = [1, 3]

        do i_bad = 1, size(bad)
            do i_position = 1, size(positions, 2)
                vecs = 1d0
                vecs(positions(1, i_position), positions(2, i_position)) = bad(i_bad)
                call omics_vector_RAP_projection(vecs, 3, 3, vecs_mask, 2, axes_mask, 2, projections, ierr)
                call assert_err(ierr, ERR_NAN_INF, "test_omics_vector_RAP_projection_rejects_nan_and_inf: " &
                                // trim(bad_names(i_bad)) // " in " // trim(position_names(i_position)), 1)
            end do
        end do
    end subroutine test_omics_vector_RAP_projection_rejects_nan_and_inf

    !> Finite input near the top of the range gives a finite projection. [huge, huge] lies on the
    !| diagonal and projects to [0, 0]; [huge, -huge] has mean 0 and stays as it is.
    !|
    !| BUG: the mean is taken as sum/n, and huge + huge overflows to Inf, so [huge, huge] projects
    !| to [-Inf, -Inf]: finite input, infinite output, which TOX must never write. Summing
    !| x_i/n instead cannot overflow where the mean itself does not (calc_tiss_avg had the same
    !| bug, see test_calc_tiss_avg_extreme_magnitudes).
    subroutine test_omics_vector_RAP_projection_extreme_magnitudes()
        integer(int32) :: ierr
        real(real64), dimension(2, 2) :: vecs, projections, expected
        logical(c_bool) :: vecs_mask(2), axes_mask(2)

        vecs(:, 1) = [huge(1d0), huge(1d0)]
        vecs(:, 2) = [huge(1d0), -huge(1d0)]
        vecs_mask = [.true., .true.]
        axes_mask = [.true., .true.]
        expected(:, 1) = [0d0, 0d0]
        expected(:, 2) = [huge(1d0), -huge(1d0)]

        call omics_vector_RAP_projection(vecs, 2, 2, vecs_mask, 2, axes_mask, 2, projections, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_omics_vector_RAP_projection_extreme_magnitudes: ierr")
        call assert_equal_array_real(projections(:, 1), expected(:, 1), 2, 0d0, &
                                     "test_omics_vector_RAP_projection_extreme_magnitudes: [huge, huge] to zero")
        call assert_equal_array_real(projections(:, 2), expected(:, 2), 2, 0d0, &
                                     "test_omics_vector_RAP_projection_extreme_magnitudes: [huge, -huge] unchanged")
    end subroutine test_omics_vector_RAP_projection_extreme_magnitudes

end module mod_test_omics_vector_RAP_projection
