!> The `compute_tissue_versatility` cases: the tissue versatility (TV) and angle to the space
!| diagonal of every gene the vectors mask selects, over the axes the axes mask selects; the
!| masks' effect on the result and on the order of the output, the one-axis special case, the
!| zero vector, the rounding, scale (subnormal to huge, also under flush-to-zero) and sign of
!| the input, the dimension checks, the checks of both selection counts against their masks,
!| and the NaN/Inf check.
!|
!| The definition, with n = n_selected_axes and v the gene restricted to the selected axes:
!| cos_phi = sum(v)/(|v|*sqrt(n)), TV = (1 - cos_phi)/(1 - 1/sqrt(n)), angle = acos(cos_phi) in
!| degrees; the expected values are derived from it. (The procedure computes the angle with
!| atan2, which is better conditioned.) Most cases select n = 4 axes: sqrt(4) = 2 is exact, so
!| the normalization 1 - 1/2 = 1/2 is exact, a uniform gene's cos_phi is exactly 1 and a
!| single-axis gene's is exactly 1/2, which makes TV exact (0 and 1). Angles other than 0 and 90
!| go through an inverse trigonometric function and a conversion to degrees, both rounded:
!| they need `ANGLE_TOLERANCE`.
!|
!| The exact expectations assume value-safe floating point, as the build uses (gfortran without
!| -ffast-math, ifx with -fp-model precise): under fast math, rewritten expressions leave some
!| single-axis TVs an ulp off 1.
module mod_test_compute_tissue_versatility
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf, &
                                            ieee_support_underflow_control, ieee_get_underflow_mode, ieee_set_underflow_mode
    use tox_tissue_versatility, only: compute_tissue_versatility
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_compute_tissue_versatility

    !> What the outputs hold before a call, so that an output the procedure never writes shows.
    real(real64), parameter :: UNWRITTEN = -1.0_real64
    !> A gene or axis the mask selects.
    logical(c_bool), parameter :: SELECTED = .true._c_bool
    !> A gene or axis the mask leaves out.
    logical(c_bool), parameter :: LEFT_OUT = .false._c_bool
    !> For a TV that is not exact: 1 - cos_phi turns cos_phi's relative rounding (a few ulps of 1)
    !| into an absolute error, which the normalization (at most 1/(1 - 1/sqrt(2)), about 3.4)
    !| enlarges. 8 epsilon covers that; a wrong formula is off by far more.
    real(real64), parameter :: TV_TOLERANCE = 8*epsilon(1.0_real64)
    !> For an angle in degrees that is not exact: the rounding of the angle in radians, enlarged
    !| by the conversion's 180/pi, leaves an angle of 30 to 60 degrees off by about 1e-14.
    !| 1e-12 degrees leaves room and still catches any wrong angle.
    real(real64), parameter :: ANGLE_TOLERANCE = 1.0e-12_real64

contains

    !> Every case of `compute_tissue_versatility`.
    function get_all_tests_compute_tissue_versatility() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(23))
        all_tests(1) = test_case("test_compute_tissue_versatility_values", test_compute_tissue_versatility_values)
        all_tests(2) = test_case("test_compute_tissue_versatility_axes_mask", test_compute_tissue_versatility_axes_mask)
        all_tests(3) = test_case("test_compute_tissue_versatility_pivot_left_out", &
                                 test_compute_tissue_versatility_pivot_left_out)
        all_tests(4) = test_case("test_compute_tissue_versatility_vectors_mask", &
                                 test_compute_tissue_versatility_vectors_mask)
        all_tests(5) = test_case("test_compute_tissue_versatility_one_selected_axis", &
                                 test_compute_tissue_versatility_one_selected_axis)
        all_tests(6) = test_case("test_compute_tissue_versatility_zero_vector", test_compute_tissue_versatility_zero_vector)
        all_tests(7) = test_case("test_compute_tissue_versatility_cos_above_one", test_compute_tissue_versatility_cos_above_one)
        all_tests(8) = test_case("test_compute_tissue_versatility_uniform_rounding", &
                                 test_compute_tissue_versatility_uniform_rounding)
        all_tests(9) = test_case("test_compute_tissue_versatility_single_axis_rounding", &
                                 test_compute_tissue_versatility_single_axis_rounding)
        all_tests(10) = test_case("test_compute_tissue_versatility_single_axis_any_n", &
                                  test_compute_tissue_versatility_single_axis_any_n)
        all_tests(11) = test_case("test_compute_tissue_versatility_near_diagonal", &
                                  test_compute_tissue_versatility_near_diagonal)
        all_tests(12) = test_case("test_compute_tissue_versatility_tiny_scale", test_compute_tissue_versatility_tiny_scale)
        all_tests(13) = test_case("test_compute_tissue_versatility_subnormal_scale", &
                                  test_compute_tissue_versatility_subnormal_scale)
        all_tests(14) = test_case("test_compute_tissue_versatility_huge_scale", test_compute_tissue_versatility_huge_scale)
        all_tests(15) = test_case("test_compute_tissue_versatility_top_of_range", &
                                  test_compute_tissue_versatility_top_of_range)
        all_tests(16) = test_case("test_compute_tissue_versatility_negative_values", &
                                  test_compute_tissue_versatility_negative_values)
        all_tests(17) = test_case("test_compute_tissue_versatility_dimensions", test_compute_tissue_versatility_dimensions)
        all_tests(18) = test_case("test_compute_tissue_versatility_no_selected_vectors", &
                                  test_compute_tissue_versatility_no_selected_vectors)
        all_tests(19) = test_case("test_compute_tissue_versatility_n_selected_vectors_mismatch", &
                                  test_compute_tissue_versatility_n_selected_vectors_mismatch)
        all_tests(20) = test_case("test_compute_tissue_versatility_n_selected_axes_below_one", &
                                  test_compute_tissue_versatility_n_selected_axes_below_one)
        all_tests(21) = test_case("test_compute_tissue_versatility_n_selected_axes_mismatch", &
                                  test_compute_tissue_versatility_n_selected_axes_mismatch)
        all_tests(22) = test_case("test_compute_tissue_versatility_rejects_nan_and_inf", &
                                  test_compute_tissue_versatility_rejects_nan_and_inf)
        all_tests(23) = test_case("test_compute_tissue_versatility_smallest_sizes", &
                                  test_compute_tissue_versatility_smallest_sizes)
    end function get_all_tests_compute_tissue_versatility

    !> Four genes over four selected axes (norm_diag = 2, normalization 1/2):
    !| - uniform [3, 3, 3, 3]: cos_phi = 12/(6*2) = 1, TV 0, angle 0, all exact;
    !| - single axis [0, 0, 5, 0]: cos_phi = 5/(5*2) = 1/2, TV (1 - 1/2)/(1/2) = 1 exact, angle 60;
    !| - two equal axes [1, 1, 0, 0]: cos_phi = 2/(sqrt(2)*2) = 1/sqrt(2),
    !|   TV (1 - 1/sqrt(2))/(1/2) = 2 - sqrt(2) = 0.58578643762690495..., angle 45;
    !| - [2, 2, 1, 0]: |v| = 3, cos_phi = 5/(3*2) = 5/6, TV (1/6)/(1/2) = 1/3,
    !|   angle acos(5/6) = 33.557309761920715... degrees (quad precision).
    subroutine test_compute_tissue_versatility_values()
        integer(int32), parameter :: n_axes = 4, n_genes = 4
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = 3.0_real64
        expression_vectors(:, 2) = [0.0_real64, 0.0_real64, 5.0_real64, 0.0_real64]
        expression_vectors(:, 3) = [1.0_real64, 1.0_real64, 0.0_real64, 0.0_real64]
        expression_vectors(:, 4) = [2.0_real64, 2.0_real64, 1.0_real64, 0.0_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = [0.0_real64, 1.0_real64, 0.58578643762690495_real64, 1.0_real64/3.0_real64]
        expected_angles = [0.0_real64, 60.0_real64, 45.0_real64, 33.557309761920715_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_values: ierr")
        ! QUESTION: the module doc says a gene expressed equally across every selected axis "scores
        ! maximally versatile" and a single-tissue gene "scores minimally"; the procedure doc and
        ! the code give the uniform gene 0 and the single-tissue gene 1. The number is a
        ! specificity, named versatility. Pinned: uniform 0, single axis 1.
        call assert_equal_real(versatilities(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_values: TV, uniform")
        call assert_equal_real(angles(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_values: angle, uniform")
        call assert_equal_real(versatilities(2), 1.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_values: TV, single axis")
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, TV_TOLERANCE, &
                                     "test_compute_tissue_versatility_values: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_values: angle")
    end subroutine test_compute_tissue_versatility_values

    !> The axes mask decides which axes count. Five axes, axis 2 left out, so n = 4:
    !| - [3, 100, 3, 3, 3] is uniform on the selected axes: TV 0, angle 0 (over all five
    !|   axes it would be far from uniform);
    !| - [0, 7, 0, 0, 0] holds all its expression on the left-out axis, so it is the zero vector
    !|   on the selected ones: TV 1, angle 90;
    !| - [0, 0, 0, 0, 5] is single-axis: cos_phi = 1/2, TV 1, angle 60 (over five axes it would
    !|   be cos_phi = 1/sqrt(5), angle 63.4).
    subroutine test_compute_tissue_versatility_axes_mask()
        integer(int32), parameter :: n_axes = 5, n_genes = 3
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = [3.0_real64, 100.0_real64, 3.0_real64, 3.0_real64, 3.0_real64]
        expression_vectors(:, 2) = [0.0_real64, 7.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
        expression_vectors(:, 3) = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 5.0_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = [SELECTED, LEFT_OUT, SELECTED, SELECTED, SELECTED]
        expected_versatilities = [0.0_real64, 1.0_real64, 1.0_real64]
        expected_angles = [0.0_real64, 90.0_real64, 60.0_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, 4, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_axes_mask: ierr")
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_axes_mask: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_axes_mask: angle")
    end subroutine test_compute_tissue_versatility_axes_mask

    !> A gene uniform on the selected axes has TV and angle exactly 0 also when the axes mask
    !| leaves out axis 1: [100, 0.1, 0.1, 0.1] with axis 1 left out, so n = 3. The rejection
    !| from the diagonal is taken of v minus its first SELECTED value, 0.1, which makes every
    !| difference, and so the rejection, exactly 0.
    subroutine test_compute_tissue_versatility_pivot_left_out()
        integer(int32), parameter :: n_axes = 4, n_genes = 1
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(n_genes), angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = [100.0_real64, 0.1_real64, 0.1_real64, 0.1_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = [LEFT_OUT, SELECTED, SELECTED, SELECTED]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, 3, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_pivot_left_out: ierr")
        ! Were the differences taken from axis 1's value, 100, or from 0, the differences and
        ! their mean would be rounded (0.1 is not a binary fraction), and TV and angle came out
        ! 2.4e-26 and 8.1e-12 degrees, or 2.3e-32 and 7.9e-15 degrees, instead of 0.
        call assert_equal_real(versatilities(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_pivot_left_out: TV")
        call assert_equal_real(angles(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_pivot_left_out: angle")
    end subroutine test_compute_tissue_versatility_pivot_left_out

    !> The vectors mask selects genes, and the output holds them compacted, in their order.
    !| Five genes over four axes: 1 uniform (TV 0, angle 0), 2 single-axis (TV 1, angle 60),
    !| 3 zero (TV 1, angle 90), 4 uniform, 5 single-axis.
    !| - mask [F, T, F, T, F]: genes 2 and 4, so TV [1, 0] and angle [60, 0];
    !| - mask [T, F, F, F, T], the selection's bounds: genes 1 and 5, so TV [0, 1], angle [0, 60].
    subroutine test_compute_tissue_versatility_vectors_mask()
        integer(int32), parameter :: n_axes = 4, n_genes = 5, n_selected = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_selected), angles(n_selected)
        real(real64) :: expected_versatilities(n_selected), expected_angles(n_selected)
        logical(c_bool) :: middle_genes(n_genes), first_and_last(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = 3.0_real64
        expression_vectors(:, 2) = [0.0_real64, 0.0_real64, 5.0_real64, 0.0_real64]
        expression_vectors(:, 3) = 0.0_real64
        expression_vectors(:, 4) = 2.0_real64
        expression_vectors(:, 5) = [4.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
        middle_genes = [LEFT_OUT, SELECTED, LEFT_OUT, SELECTED, LEFT_OUT]
        first_and_last = [SELECTED, LEFT_OUT, LEFT_OUT, LEFT_OUT, SELECTED]
        axes_selection_mask = SELECTED

        expected_versatilities = [1.0_real64, 0.0_real64]
        expected_angles = [60.0_real64, 0.0_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, middle_genes, n_selected, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_vectors_mask: ierr, genes 2 and 4")
        call assert_equal_array_real(versatilities, expected_versatilities, n_selected, 0.0_real64, &
                                     "test_compute_tissue_versatility_vectors_mask: TV, genes 2 and 4")
        call assert_equal_array_real(angles, expected_angles, n_selected, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_vectors_mask: angle, genes 2 and 4")

        expected_versatilities = [0.0_real64, 1.0_real64]
        expected_angles = [0.0_real64, 60.0_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, first_and_last, n_selected, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_vectors_mask: ierr, genes 1 and 5")
        call assert_equal_array_real(versatilities, expected_versatilities, n_selected, 0.0_real64, &
                                     "test_compute_tissue_versatility_vectors_mask: TV, genes 1 and 5")
        call assert_equal_array_real(angles, expected_angles, n_selected, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_vectors_mask: angle, genes 1 and 5")
    end subroutine test_compute_tissue_versatility_vectors_mask

    !> One selected axis (the middle of three): the normalization 1 - 1/sqrt(1) is 0, so the
    !| procedure skips the formula and returns TV 0 and angle 0 for every gene. For a positive
    !| value that agrees with the formula's limit (a one-axis gene lies on the one-axis diagonal).
    !| Genes: [9, 2, 9] (positive on the selected axis), [9, 0, 9] (zero there), [9, -2, 9]
    !| (negative there).
    subroutine test_compute_tissue_versatility_one_selected_axis()
        integer(int32), parameter :: n_axes = 3, n_genes = 3
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = [9.0_real64, 2.0_real64, 9.0_real64]
        expression_vectors(:, 2) = [9.0_real64, 0.0_real64, 9.0_real64]
        expression_vectors(:, 3) = [9.0_real64, -2.0_real64, 9.0_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = [LEFT_OUT, SELECTED, LEFT_OUT]
        expected = 0.0_real64
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, 1, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_one_selected_axis: ierr")
        ! QUESTION: gene 2 is the zero vector on the selected axis, which with two or more
        ! selected axes gets TV 1 and angle 90 (test_compute_tissue_versatility_zero_vector),
        ! here TV 0 and angle 0. Gene 3 points against the diagonal (cos_phi = -1, angle 180
        ! with the formula), here TV 0 and angle 0. Pinned: every gene 0 and 0.
        call assert_equal_array_real(versatilities, expected, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_one_selected_axis: TV")
        call assert_equal_array_real(angles, expected, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_one_selected_axis: angle")
    end subroutine test_compute_tissue_versatility_one_selected_axis

    !> The zero vector has no direction; the documented result is TV 1 and angle 90, exactly.
    !| Gene 2, single-axis, is there to compare with.
    subroutine test_compute_tissue_versatility_zero_vector()
        integer(int32), parameter :: n_axes = 4, n_genes = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = 0.0_real64
        expression_vectors(:, 2) = [0.0_real64, 6.0_real64, 0.0_real64, 0.0_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_zero_vector: ierr")
        ! QUESTION: the zero vector's TV 1 is the same as a single-tissue gene's (gene 2), so TV
        ! alone cannot tell a gene expressed nowhere from one expressed in one tissue. Its angle,
        ! 90, is out of band: a non-negative gene's angle is at most acos(1/sqrt(n)), 60 here.
        ! Pinned: TV 1 and angle 90.
        call assert_equal_real(versatilities(1), 1.0_real64, 0.0_real64, "test_compute_tissue_versatility_zero_vector: TV")
        call assert_equal_real(angles(1), 90.0_real64, 0.0_real64, "test_compute_tissue_versatility_zero_vector: angle")
        call assert_equal_real(versatilities(2), 1.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_zero_vector: TV, single axis")
        call assert_equal_real(angles(2), 60.0_real64, ANGLE_TOLERANCE, &
                               "test_compute_tissue_versatility_zero_vector: angle, single axis")
    end subroutine test_compute_tissue_versatility_zero_vector

    !> A uniform gene over three axes, [1, 1, 1]: TV and angle are exactly 0.
    !| Regression: cos_phi = 3/(sqrt(3)*sqrt(3)) used to round to 1.0000000000000002, just above
    !| 1, where acos is NaN; a clamp to [-1, 1] hid it. The angle now comes from atan2 of the
    !| rejection from the diagonal, which is exactly 0 for a uniform gene.
    subroutine test_compute_tissue_versatility_cos_above_one()
        integer(int32), parameter :: n_axes = 3, n_genes = 1
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(n_genes), angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors = 1.0_real64
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_cos_above_one: ierr")
        call assert_equal_real(versatilities(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_cos_above_one: TV")
        call assert_equal_real(angles(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_cos_above_one: angle")
    end subroutine test_compute_tissue_versatility_cos_above_one

    !> A uniform gene has angle 0 and TV 0 over any number of axes. [1, 1] over two axes and
    !| [0.3, 0.3, 0.3] over three.
    subroutine test_compute_tissue_versatility_uniform_rounding()
        real(real64) :: two_axes(2, 1), three_axes(3, 1), versatilities(1), angles(1)
        logical(c_bool) :: genes_selection_mask(1), two_axes_mask(2), three_axes_mask(3)
        integer(int32) :: ierr

        two_axes = 1.0_real64
        three_axes = 0.3_real64
        genes_selection_mask = SELECTED
        two_axes_mask = SELECTED
        three_axes_mask = SELECTED

        ! Regression: cos_phi = dot/(sqrt(norm_v)*sqrt(n)) used to round to just below 1
        ! (0.99999999999999978 for [1, 1], 0.99999999999999989 for [0.3, 0.3, 0.3]), and acos,
        ! whose slope is infinite at 1, turned one ulp into 2.1e-8 rad: the angle came out 1.2e-6
        ! and 8.5e-7 degrees instead of 0. About a third of all uniform values over 2, 5, 7 or 8
        ! axes were hit. The angle now comes from atan2(|rejection from the diagonal|, projection
        ! on it), and a uniform gene's rejection is exactly 0: TV and angle are exactly 0.
        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(2, 1, two_axes, genes_selection_mask, 1, two_axes_mask, 2, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_uniform_rounding: ierr, 2 axes")
        call assert_equal_real(versatilities(1), 0.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_uniform_rounding: TV, 2 axes")
        call assert_equal_real(angles(1), 0.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_uniform_rounding: angle, 2 axes")

        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(3, 1, three_axes, genes_selection_mask, 1, three_axes_mask, 3, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_uniform_rounding: ierr, 3 axes")
        call assert_equal_real(versatilities(1), 0.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_uniform_rounding: TV, 3 axes")
        call assert_equal_real(angles(1), 0.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_uniform_rounding: angle, 3 axes")
    end subroutine test_compute_tissue_versatility_uniform_rounding

    !> A gene expressed in a single axis has TV exactly 1 at any value, over four axes: its
    !| cos_phi = v_k/(|v_k|*2) = 1/2 exactly, since |v| = sqrt(v_k**2) = |v_k| is exact, and the
    !| normalization is 1/2. Angle 60. Five genes, each a value that is not a power of two:
    !| 0.1 on axis 1, 0.1 on axis 2, 1e-300 on axis 1, huge on axis 1, 3 on axis 4.
    subroutine test_compute_tissue_versatility_single_axis_rounding()
        integer(int32), parameter :: n_axes = 4, n_genes = 5
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors = 0.0_real64
        expression_vectors(1, 1) = 0.1_real64
        expression_vectors(2, 2) = 0.1_real64
        expression_vectors(1, 3) = 1.0e-300_real64
        expression_vectors(1, 4) = huge(1.0_real64)
        expression_vectors(4, 5) = 3.0_real64
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = 1.0_real64
        expected_angles = 60.0_real64
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_single_axis_rounding: ierr")
        ! Regression: 1 - cos_phi came from r**2/(|v|*(|v| + p)) whenever p > 0, and r**2, a
        ! rounded sum of squared differences from the mean, is not exactly 3/4 of |v|**2. TV came
        ! out 0.99999999999999989, 1.0000000000000002 (above the documented [0, 1]),
        ! 0.99999999999999967 and 1.0000000000000002 for genes 1 to 4 (gene 5 happened to come
        ! out 1); about a third of single-axis values missed 1 by an ulp or two. That formula
        ! now serves only cos_phi > 3/4; a single-axis gene's cos_phi is at most 1/sqrt(2), and it
        ! gets 1 - cos_phi, which is exact here (test_compute_tissue_versatility_single_axis_any_n).
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_single_axis_rounding: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_single_axis_rounding: angle")
    end subroutine test_compute_tissue_versatility_single_axis_rounding

    !> A single-axis gene has TV exactly 1 over any number of axes, not only four: its sum is
    !| v_k and its length sqrt(v_k**2) = |v_k| exactly (a correctly rounded square root of a
    !| rounded square), so cos_phi = (sum/|v|)/sqrt(n) is 1/sqrt(n) computed exactly as the
    !| normalization 1 - 1/sqrt(n) computes it, and TV is that normalization divided by itself.
    !| Each gene is a value that is not a power of two:
    !| - 7 on axis 2 of two: cos_phi = 1/sqrt(2), angle 45;
    !| - 0.1 on axis 1 of three: cos_phi = 1/sqrt(3), angle acos(1/sqrt(3)) = 54.735610317245346 degrees;
    !| - 3 on axis 5 of seven: cos_phi = 1/sqrt(7), angle acos(1/sqrt(7)) = 67.792345701403513 degrees
    !|   (quad precision).
    subroutine test_compute_tissue_versatility_single_axis_any_n()
        real(real64) :: two_axes(2, 1), three_axes(3, 1), seven_axes(7, 1), versatilities(1), angles(1)
        logical(c_bool) :: genes_selection_mask(1), two_axes_mask(2), three_axes_mask(3), seven_axes_mask(7)
        integer(int32) :: ierr

        two_axes = 0.0_real64
        two_axes(2, 1) = 7.0_real64
        three_axes = 0.0_real64
        three_axes(1, 1) = 0.1_real64
        seven_axes = 0.0_real64
        seven_axes(5, 1) = 3.0_real64
        genes_selection_mask = SELECTED
        two_axes_mask = SELECTED
        three_axes_mask = SELECTED
        seven_axes_mask = SELECTED

        ! Regression: over two and three axes (cos_phi above 1/2) 1 - cos_phi came from
        ! r**2/(|v|*(|v| + p)), and 7 over two axes gave TV 0.99999999999999967, 0.1 over three
        ! 1.0000000000000002 (above the documented [0, 1]). Below 1/2, 1 - p/|v|, with p =
        ! sum/sqrt(n) rounded before the division by |v|, gave 3 over seven axes 0.99999999999999978.
        ! Many single-axis values over 2 to 64 axes missed 1 by an ulp or two.
        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(2, 1, two_axes, genes_selection_mask, 1, two_axes_mask, 2, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_single_axis_any_n: ierr, 2 axes")
        call assert_equal_real(versatilities(1), 1.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_single_axis_any_n: TV, 2 axes")
        call assert_equal_real(angles(1), 45.0_real64, ANGLE_TOLERANCE, &
                               "test_compute_tissue_versatility_single_axis_any_n: angle, 2 axes")

        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(3, 1, three_axes, genes_selection_mask, 1, three_axes_mask, 3, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_single_axis_any_n: ierr, 3 axes")
        call assert_equal_real(versatilities(1), 1.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_single_axis_any_n: TV, 3 axes")
        call assert_equal_real(angles(1), 54.735610317245346_real64, ANGLE_TOLERANCE, &
                               "test_compute_tissue_versatility_single_axis_any_n: angle, 3 axes")

        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(7, 1, seven_axes, genes_selection_mask, 1, seven_axes_mask, 7, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_single_axis_any_n: ierr, 7 axes")
        call assert_equal_real(versatilities(1), 1.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_single_axis_any_n: TV, 7 axes")
        call assert_equal_real(angles(1), 67.792345701403513_real64, ANGLE_TOLERANCE, &
                               "test_compute_tissue_versatility_single_axis_any_n: angle, 7 axes")
    end subroutine test_compute_tissue_versatility_single_axis_any_n

    !> Two genes over four axes whose TV and angle are not exact, checked to within a few ulps
    !| of the value computed in quad precision:
    !| - near the diagonal, [1, 1, 1, 1 + d] with d = 2**-20: the mean is 1 + d/4, so the
    !|   rejection is d*[-1/4, -1/4, -1/4, 3/4] and r**2 = 3*d**2/4; |v| = sqrt(4 + 2*d + d**2),
    !|   p = (4 + d)/2, 1 - cos_phi = r**2/(|v|*(|v| + p)), TV = 2*(1 - cos_phi)
    !|   = 1.7053017526726838e-13, angle atan2(r, p) = 2.3660463694442089e-05 degrees;
    !| - [6, 6, 3, 1]: sum 16, |v| = sqrt(82), cos_phi = 8/sqrt(82) = 0.8835 (between the
    !|   formulas' split at 3/4 and 0.9), TV = 2*(1 - 8/sqrt(82)) = 0.23309558280245529,
    !|   angle acos(8/sqrt(82)) = 27.938352729602350 degrees.
    !| The tolerance is 4 ulps of the expected value, relative, since TV and angle of the first
    !| gene are far below `TV_TOLERANCE` and `ANGLE_TOLERANCE`.
    subroutine test_compute_tissue_versatility_near_diagonal()
        integer(int32), parameter :: n_axes = 4, n_genes = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr, i_gene

        expression_vectors(:, 1) = [1.0_real64, 1.0_real64, 1.0_real64, 1.0_real64 + scale(1.0_real64, -20)]
        expression_vectors(:, 2) = [6.0_real64, 6.0_real64, 3.0_real64, 1.0_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = [1.7053017526726838e-13_real64, 0.23309558280245529_real64]
        expected_angles = [2.3660463694442089e-05_real64, 27.938352729602350_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_near_diagonal: ierr")
        ! For the first gene, cos_phi = 1 - 8.5e-14 is rounded to a multiple of 2**-53, so
        ! 1 - cos_phi would cancel: TV would come out 1.7053026e-13 (5e-7 off, relative), and
        ! acos(cos_phi) would give 2.3660469e-05 degrees (2e-7 off). For the second gene,
        ! 1 - cos_phi (the formula below the split) is 6 ulps off, r**2/(|v|*(|v| + p)) within
        ! one: moving the split up to 0.9 fails here.
        do i_gene = 1, n_genes
            call assert_equal_real(versatilities(i_gene), expected_versatilities(i_gene), &
                                   4*spacing(expected_versatilities(i_gene)), &
                                   "test_compute_tissue_versatility_near_diagonal: TV, gene "//char(ichar('0') + i_gene))
            call assert_equal_real(angles(i_gene), expected_angles(i_gene), 4*spacing(expected_angles(i_gene)), &
                                   "test_compute_tissue_versatility_near_diagonal: angle, gene "//char(ichar('0') + i_gene))
        end do
    end subroutine test_compute_tissue_versatility_near_diagonal

    !> TV and angle depend on a gene's direction, not on its length, so a gene expressed at a
    !| tiny scale scores like the same gene at scale 1. Four axes:
    !| - uniform 2**-20 (about 9.5e-7): cos_phi = 4*2**-20/(2*2**-20*2) = 1, TV 0, angle 0;
    !| - single-axis [2**-20, 0, 0, 0]: cos_phi = 1/2, TV 1, angle 60;
    !| - uniform 2**-700: the same as the first, TV 0, angle 0.
    !| Every value is a power of two, so the right results are exact (bar the 60).
    subroutine test_compute_tissue_versatility_tiny_scale()
        integer(int32), parameter :: n_axes = 4, n_genes = 3
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = scale(1.0_real64, -20)
        expression_vectors(:, 2) = 0.0_real64
        expression_vectors(1, 2) = scale(1.0_real64, -20)
        expression_vectors(:, 3) = scale(1.0_real64, -700)
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = [0.0_real64, 1.0_real64, 0.0_real64]
        expected_angles = [0.0_real64, 60.0_real64, 0.0_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_tiny_scale: ierr")
        ! Regression: the zero-vector test compared the SQUARED norm, norm_v = sum(v**2), with
        ! sqrt(epsilon) = 1.5e-8, so every gene with |v| <= 1.2e-4 counted as the zero vector
        ! and got TV 1, angle 90: genes 1 and 2 (norm_v = 2**-38 and 2**-40). For gene 3,
        ! v**2 = 2**-1400 underflowed to 0. Now only an exact zero is the zero vector, and each
        ! gene is scaled by a power of two near its largest component before it is squared.
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_tiny_scale: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_tiny_scale: angle")
    end subroutine test_compute_tissue_versatility_tiny_scale

    !> Subnormal expression scores like the same gene at scale 1. Four axes:
    !| - [3, 0, 5, 1]*2**-1064, every nonzero value subnormal (5*2**-1064 is about 2.5e-320);
    !| - [3, 0, 5, 1], the same direction at scale 1: sum 9, |v| = sqrt(35),
    !|   cos_phi = 9/(2*sqrt(35)), TV = 2*(1 - cos_phi) = 0.47872234148867016,
    !|   angle acos(cos_phi) = 40.479451925712418 degrees (quad precision);
    !| - [0, 2**-1074, 0, 0], the smallest subnormal on one axis: TV 1, angle 60.
    !| Genes 1 and 2 differ by a power of two, so every operation on them rounds alike and
    !| their results are identical.
    subroutine test_compute_tissue_versatility_subnormal_scale()
        integer(int32), parameter :: n_axes = 4, n_genes = 3
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 2) = [3.0_real64, 0.0_real64, 5.0_real64, 1.0_real64]
        expression_vectors(:, 1) = scale(expression_vectors(:, 2), -1064)
        expression_vectors(:, 3) = 0.0_real64
        expression_vectors(2, 3) = scale(1.0_real64, -1074)
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = [0.47872234148867016_real64, 0.47872234148867016_real64, 1.0_real64]
        expected_angles = [40.479451925712418_real64, 40.479451925712418_real64, 60.0_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_subnormal_scale: ierr")
        ! Each gene is scaled by a power of two near its largest value, but 2**1061 (for gene 1)
        ! and 2**1073 (gene 3) overflow. Without the floor on that exponent, which scales such a
        ! gene by 2**1023 instead, both genes came out NaN.
        call assert_equal_real(versatilities(1), versatilities(2), 0.0_real64, &
                               "test_compute_tissue_versatility_subnormal_scale: TV, as at scale 1")
        call assert_equal_real(angles(1), angles(2), 0.0_real64, &
                               "test_compute_tissue_versatility_subnormal_scale: angle, as at scale 1")
        call assert_equal_real(versatilities(3), 1.0_real64, 0.0_real64, &
                               "test_compute_tissue_versatility_subnormal_scale: TV, single axis")
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, TV_TOLERANCE, &
                                     "test_compute_tissue_versatility_subnormal_scale: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_subnormal_scale: angle")
    end subroutine test_compute_tissue_versatility_subnormal_scale

    !> The other end of the scale: uniform 2**600 over four axes, and uniform huge. Both are
    !| uniform, so TV 0 and angle 0.
    subroutine test_compute_tissue_versatility_huge_scale()
        integer(int32), parameter :: n_axes = 4, n_genes = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = scale(1.0_real64, 600)
        expression_vectors(:, 2) = huge(1.0_real64)
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected = 0.0_real64
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_huge_scale: ierr")
        ! Regression: norm_v = sum(v**2) overflowed. Gene 1: (2**600)**2 = +Inf, so cos_phi =
        ! 2**602/(Inf*2) = 0, which gave TV 2 (outside the documented [0, 1]) and angle 90.
        ! Gene 2: the dot product overflowed as well, cos_phi = Inf/Inf = NaN, and TV and angle
        ! were NaN unoptimized (a NaN from finite input). Each gene is now scaled by a power of
        ! two near its largest component before it is squared.
        call assert_equal_array_real(versatilities, expected, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_huge_scale: TV")
        call assert_equal_array_real(angles, expected, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_huge_scale: angle")
    end subroutine test_compute_tissue_versatility_huge_scale

    !> The top of the range, where the scale factor meets the bottom: [3, 0, 5, 1]*2**1021
    !| (largest value 5*2**1021, above 2**1022) scores like [3, 0, 5, 1] at scale 1, TV
    !| 0.47872234148867016 and angle 40.479451925712418 degrees (derived in
    !| `test_compute_tissue_versatility_subnormal_scale`), both with IEEE gradual underflow and,
    !| where the processor can switch it off, with flush-to-zero, which a caller (a library
    !| built with fast math) may have left on.
    subroutine test_compute_tissue_versatility_top_of_range()
        integer(int32), parameter :: n_axes = 4, n_genes = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        logical :: gradual_on_entry
        integer(int32) :: ierr

        expression_vectors(:, 2) = [3.0_real64, 0.0_real64, 5.0_real64, 1.0_real64]
        expression_vectors(:, 1) = scale(expression_vectors(:, 2), 1021)
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = 0.47872234148867016_real64
        expected_angles = 40.479451925712418_real64

        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_top_of_range: ierr")
        call assert_equal_real(versatilities(1), versatilities(2), 0.0_real64, &
                               "test_compute_tissue_versatility_top_of_range: TV, as at scale 1")
        call assert_equal_real(angles(1), angles(2), 0.0_real64, &
                               "test_compute_tissue_versatility_top_of_range: angle, as at scale 1")
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, TV_TOLERANCE, &
                                     "test_compute_tissue_versatility_top_of_range: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_top_of_range: angle")

        if (.not. ieee_support_underflow_control(1.0_real64)) return
        ! Regression: gene 1 was scaled by 2**-1024, a subnormal, which flush-to-zero turned into
        ! 0, so every scaled value was 0 and TV came out NaN, the angle 0. The scale factor is
        ! now at least 2**-1022, the smallest normal number.
        call ieee_get_underflow_mode(gradual_on_entry)
        call ieee_set_underflow_mode(.false.)
        versatilities = UNWRITTEN
        angles = UNWRITTEN
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call ieee_set_underflow_mode(gradual_on_entry)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_top_of_range: ierr, flush-to-zero")
        call assert_equal_real(versatilities(1), versatilities(2), 0.0_real64, &
                               "test_compute_tissue_versatility_top_of_range: TV, flush-to-zero")
        call assert_equal_real(angles(1), angles(2), 0.0_real64, &
                               "test_compute_tissue_versatility_top_of_range: angle, flush-to-zero")
    end subroutine test_compute_tissue_versatility_top_of_range

    !> Negative expression over four axes:
    !| - uniform -3: cos_phi = -12/(6*2) = -1, TV (1 + 1)/(1/2) = 4, angle 180;
    !| - [5, -5, 0, 0]: sum 0, cos_phi = 0, TV 1/(1/2) = 2, angle 90.
    !| TV is exact; the angles pi and pi/2 are rounded, and their conversion lands on 180 and 90
    !| within an ulp.
    subroutine test_compute_tissue_versatility_negative_values()
        integer(int32), parameter :: n_axes = 4, n_genes = 2
        real(real64) :: expression_vectors(n_axes, n_genes)
        real(real64) :: versatilities(n_genes), angles(n_genes), expected_versatilities(n_genes), expected_angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors(:, 1) = -3.0_real64
        expression_vectors(:, 2) = [5.0_real64, -5.0_real64, 0.0_real64, 0.0_real64]
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        expected_versatilities = [4.0_real64, 2.0_real64]
        expected_angles = [180.0_real64, 90.0_real64]
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        ! QUESTION: the documented range of TV is [0, 1], which holds for non-negative expression
        ! only; the wrapper checks expression_vectors for NaN/Inf but not for sign, so a negative
        ! value is accepted and TV goes up to 2/(1 - 1/sqrt(n)) (4 here, 6.83 over two axes), the
        ! angle up to 180. Either the input must be non-negative (a DM_MIN(0) on
        ! expression_vectors) or the documented range is wrong. Pinned: accepted, TV [4, 2],
        ! angle [180, 90].
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_negative_values: ierr")
        call assert_equal_array_real(versatilities, expected_versatilities, n_genes, 0.0_real64, &
                                     "test_compute_tissue_versatility_negative_values: TV")
        call assert_equal_array_real(angles, expected_angles, n_genes, ANGLE_TOLERANCE, &
                                     "test_compute_tissue_versatility_negative_values: angle")
    end subroutine test_compute_tissue_versatility_negative_values

    !> n_axes (argument 1) and n_vectors (argument 2), each on its own: 0 is ERR_EMPTY_INPUT, -1
    !| ERR_INVALID_INPUT. An empty mask counts 0, which n_selected_axes (at least 1) and
    !| n_selected_vectors (at least 1) cannot match; the dimensions are checked first, so they
    !| are the argument blamed. 1, the smallest valid size of both, is in
    !| `test_compute_tissue_versatility_smallest_sizes`.
    subroutine test_compute_tissue_versatility_dimensions()
        real(real64) :: expression_vectors(2, 2), versatilities(1), angles(1)
        logical(c_bool) :: genes_selection_mask(2), axes_selection_mask(2)
        integer(int32) :: ierr

        expression_vectors = 1.0_real64
        genes_selection_mask = [SELECTED, LEFT_OUT]
        axes_selection_mask = [SELECTED, LEFT_OUT]

        call compute_tissue_versatility(0, 2, expression_vectors, genes_selection_mask, 1, axes_selection_mask, 1, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_tissue_versatility_dimensions: n_axes = 0", 1)
        call compute_tissue_versatility(-1, 2, expression_vectors, genes_selection_mask, 1, axes_selection_mask, 1, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_dimensions: n_axes = -1", 1)
        call compute_tissue_versatility(2, 0, expression_vectors, genes_selection_mask, 1, axes_selection_mask, 1, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_tissue_versatility_dimensions: n_vectors = 0", 2)
        call compute_tissue_versatility(2, -1, expression_vectors, genes_selection_mask, 1, axes_selection_mask, 1, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_dimensions: n_vectors = -1", 2)
    end subroutine test_compute_tissue_versatility_dimensions

    !> A vectors mask that selects nothing, n_selected_vectors (argument 5) = 0, and -1 below it.
    subroutine test_compute_tissue_versatility_no_selected_vectors()
        integer(int32), parameter :: n_axes = 2, n_genes = 2
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(1), angles(1)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors = 1.0_real64
        genes_selection_mask = LEFT_OUT
        axes_selection_mask = SELECTED

        ! QUESTION: n_selected_vectors sizes the outputs, so the generator validates it as a
        ! dimension, and an empty selection, which only asks for two empty outputs, is
        ! ERR_EMPTY_INPUT. tox_gene_centroids' mean_vector accepts an empty selection. Pinned:
        ! ERR_EMPTY_INPUT at argument 5.
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, 0, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_tissue_versatility_no_selected_vectors: 0", 5)
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, -1, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_no_selected_vectors: -1", 5)
    end subroutine test_compute_tissue_versatility_no_selected_vectors

    !> n_selected_vectors (argument 5) must be the count of .true. in vectors_selection_mask:
    !| the mask [T, F, T] counts 2, so 3 (also n_vectors, valid by itself) and 1 are
    !| ERR_INVALID_INPUT. The selected genes sit first and last, so the count covers the whole mask.
    subroutine test_compute_tissue_versatility_n_selected_vectors_mismatch()
        integer(int32), parameter :: n_axes = 2, n_genes = 3
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(3), angles(3)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors = 1.0_real64
        genes_selection_mask = [SELECTED, LEFT_OUT, SELECTED]
        axes_selection_mask = SELECTED

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, 3, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_n_selected_vectors_mismatch: count + 1", 5)
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, 1, &
                                        axes_selection_mask, n_axes, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_n_selected_vectors_mismatch: count - 1", 5)
    end subroutine test_compute_tissue_versatility_n_selected_vectors_mismatch

    !> n_selected_axes (argument 7) has a floor of 1: an axes mask that selects nothing with a
    !| count of 0, and a count of -1, are ERR_INVALID_INPUT.
    subroutine test_compute_tissue_versatility_n_selected_axes_below_one()
        integer(int32), parameter :: n_axes = 2, n_genes = 1
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(n_genes), angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors = 1.0_real64
        genes_selection_mask = SELECTED
        axes_selection_mask = LEFT_OUT

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, 0, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_n_selected_axes_below_one: 0", 7)
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, -1, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_n_selected_axes_below_one: -1", 7)
    end subroutine test_compute_tissue_versatility_n_selected_axes_below_one

    !> n_selected_axes (argument 7) must be the count of .true. in axes_selection_mask: the mask
    !| [T, F, T, T] counts 3, so 4 (also n_axes) and 2 are ERR_INVALID_INPUT.
    subroutine test_compute_tissue_versatility_n_selected_axes_mismatch()
        integer(int32), parameter :: n_axes = 4, n_genes = 1
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(n_genes), angles(n_genes)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr

        expression_vectors = 1.0_real64
        genes_selection_mask = SELECTED
        axes_selection_mask = [SELECTED, LEFT_OUT, SELECTED, SELECTED]

        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, 4, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_n_selected_axes_mismatch: count + 1", 7)
        call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, n_genes, &
                                        axes_selection_mask, 2, versatilities, angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_tissue_versatility_n_selected_axes_mismatch: count - 1", 7)
    end subroutine test_compute_tissue_versatility_n_selected_axes_mismatch

    !> NaN, +Inf and -Inf in expression_vectors (argument 3) are ERR_NAN_INF. Each sits in the
    !| last element, on an axis and of a gene both masks leave out: the whole matrix is checked.
    subroutine test_compute_tissue_versatility_rejects_nan_and_inf()
        integer(int32), parameter :: n_axes = 2, n_genes = 2
        character(*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        real(real64) :: expression_vectors(n_axes, n_genes), versatilities(1), angles(1), bad(3)
        logical(c_bool) :: genes_selection_mask(n_genes), axes_selection_mask(n_axes)
        integer(int32) :: ierr, i_bad

        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        expression_vectors = 1.0_real64
        genes_selection_mask = [SELECTED, LEFT_OUT]
        axes_selection_mask = [SELECTED, LEFT_OUT]
        do i_bad = 1, size(bad)
            expression_vectors(n_axes, n_genes) = bad(i_bad)
            call compute_tissue_versatility(n_axes, n_genes, expression_vectors, genes_selection_mask, 1, &
                                            axes_selection_mask, 1, versatilities, angles, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_tissue_versatility_rejects_nan_and_inf: " &
                            //trim(bad_names(i_bad)), 3)
        end do
    end subroutine test_compute_tissue_versatility_rejects_nan_and_inf

    !> n_axes = n_vectors = 1, the smallest valid sizes: one axis, so TV 0 and angle 0 (the
    !| one-axis rule of `test_compute_tissue_versatility_one_selected_axis`).
    subroutine test_compute_tissue_versatility_smallest_sizes()
        real(real64) :: expression_vectors(1, 1), versatilities(1), angles(1)
        logical(c_bool) :: genes_selection_mask(1), axes_selection_mask(1)
        integer(int32) :: ierr

        expression_vectors = 4.2_real64
        genes_selection_mask = SELECTED
        axes_selection_mask = SELECTED
        versatilities = UNWRITTEN
        angles = UNWRITTEN

        call compute_tissue_versatility(1, 1, expression_vectors, genes_selection_mask, 1, axes_selection_mask, 1, &
                                        versatilities, angles, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_tissue_versatility_smallest_sizes: ierr")
        call assert_equal_real(versatilities(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_smallest_sizes: TV")
        call assert_equal_real(angles(1), 0.0_real64, 0.0_real64, "test_compute_tissue_versatility_smallest_sizes: angle")
    end subroutine test_compute_tissue_versatility_smallest_sizes

end module mod_test_compute_tissue_versatility
