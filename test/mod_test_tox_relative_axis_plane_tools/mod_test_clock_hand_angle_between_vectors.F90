!> The `clock_hand_angle_between_vectors` cases: hand-derived angles for turns whose cosine is
!| exact (0, 45, 60, 90, 120 and 180 degrees, and acos(-13/14)), the sign convention -- the sign of
!| the reference's component along `v2 - (v1 . v2) v1` -- including what the reference's length and
!| its part along `v1` do not change, parallel and antiparallel vectors (unsigned), the references
!| that orient nothing (ERR_INVALID_INPUT, blaming no argument), the n_dims bound and the rejection
!| of NaN and Inf in each vector.
!|
!| The angle is acos of the dot product, so an exact angle needs an exact cosine: the inputs are
!| unit vectors along the axes or with a component 1/2, where only the component that does not
!| enter the dot product is rounded. The expected values are fractions of PI; they meet acos's
!| result within ANGLE_TOL.
!|
!| Three cases pin regressions of the old acos form and its absolute is_close floors: identical
!| vectors, a turn below 1e-6 rad and a short reference.
module mod_test_clock_hand_angle_between_vectors
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use f42_math_impl, only: PI
    use tox_relative_axis_plane_tools
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

    real(real64), parameter :: ANGLE_TOL = 4e-15_real64
        !! A few ulps of an angle up to PI: acos rounds its result, PI/3 and 2*PI/3 round the
        !! division, and inputs normalized by a square root carry its rounding into the dot
        !! product, whose sum an optimizer may reorder.

contains

    !> Get array of all available tests.
    function get_all_tests_clock_hand_angle_between_vectors() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(18))

        all_tests(1) = test_case("test_clock_hand_angle_between_vectors_axis_turns", &
                                 test_clock_hand_angle_between_vectors_axis_turns)
        all_tests(2) = test_case("test_clock_hand_angle_between_vectors_exact_cosines", &
                                 test_clock_hand_angle_between_vectors_exact_cosines)
        all_tests(3) = test_case("test_clock_hand_angle_between_vectors_reversed_reference", &
                                 test_clock_hand_angle_between_vectors_reversed_reference)
        all_tests(4) = test_case("test_clock_hand_angle_between_vectors_reference_part_along_v1", &
                                 test_clock_hand_angle_between_vectors_reference_part_along_v1)
        all_tests(5) = test_case("test_clock_hand_angle_between_vectors_anticommutative", &
                                 test_clock_hand_angle_between_vectors_anticommutative)
        all_tests(6) = test_case("test_clock_hand_angle_between_vectors_sign_taken_at_v1", &
                                 test_clock_hand_angle_between_vectors_sign_taken_at_v1)
        all_tests(7) = test_case("test_clock_hand_angle_between_vectors_rap_three_tissues", &
                                 test_clock_hand_angle_between_vectors_rap_three_tissues)
        all_tests(8) = test_case("test_clock_hand_angle_between_vectors_mixed_signs", &
                                 test_clock_hand_angle_between_vectors_mixed_signs)
        all_tests(9) = test_case("test_clock_hand_angle_between_vectors_parallel_unsigned", &
                                 test_clock_hand_angle_between_vectors_parallel_unsigned)
        all_tests(10) = test_case("test_clock_hand_angle_between_vectors_degenerate_references", &
                                  test_clock_hand_angle_between_vectors_degenerate_references)
        all_tests(11) = test_case("test_clock_hand_angle_between_vectors_dimensions", &
                                  test_clock_hand_angle_between_vectors_dimensions)
        all_tests(12) = test_case("test_clock_hand_angle_between_vectors_rejects_nan_and_inf", &
                                  test_clock_hand_angle_between_vectors_rejects_nan_and_inf)
        all_tests(13) = test_case("test_clock_hand_angle_between_vectors_identical_diagonal", &
                                  test_clock_hand_angle_between_vectors_identical_diagonal)
        all_tests(14) = test_case("test_clock_hand_angle_between_vectors_small_turn_sign", &
                                  test_clock_hand_angle_between_vectors_small_turn_sign)
        all_tests(15) = test_case("test_clock_hand_angle_between_vectors_short_reference", &
                                  test_clock_hand_angle_between_vectors_short_reference)
        all_tests(16) = test_case("test_clock_hand_angle_between_vectors_zero_vector", &
                                  test_clock_hand_angle_between_vectors_zero_vector)
        all_tests(17) = test_case("test_clock_hand_angle_between_vectors_unnormalized", &
                                  test_clock_hand_angle_between_vectors_unnormalized)
        all_tests(18) = test_case("test_clock_hand_angle_between_vectors_extreme_magnitudes", &
                                  test_clock_hand_angle_between_vectors_extreme_magnitudes)
    end function get_all_tests_clock_hand_angle_between_vectors

    !> Quarter and half turns in the plane, from e1, with the reference e2 -- e1 turned a quarter
    !| turn counter-clockwise, so the familiar determinant convention: counter-clockwise is
    !| positive. The cosines are 1, 0, -1 and 0, exact.
    subroutine test_clock_hand_angle_between_vectors_axis_turns()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle
        integer(int32) :: ierr

        v1 = [1.0_real64, 0.0_real64]
        reference = [0.0_real64, 1.0_real64]

        v2 = [1.0_real64, 0.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_axis_turns: e1 -> e1 ierr")
        ! acos(1) = 0 exactly
        call assert_equal_real(signed_angle, 0.0_real64, 0.0_real64, "test_clock_hand_angle_between_vectors_axis_turns: e1 -> e1")

        v2 = [0.0_real64, 1.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_axis_turns: e1 -> e2 ierr")
        call assert_equal_real(signed_angle, PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_axis_turns: e1 -> e2 is +PI/2")

        v2 = [0.0_real64, -1.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_axis_turns: e1 -> -e2 ierr")
        call assert_equal_real(signed_angle, -PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_axis_turns: e1 -> -e2 is -PI/2")

        v2 = [-1.0_real64, 0.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_axis_turns: e1 -> -e1 ierr")
        call assert_equal_real(signed_angle, PI, ANGLE_TOL, "test_clock_hand_angle_between_vectors_axis_turns: e1 -> -e1 is PI")
    end subroutine test_clock_hand_angle_between_vectors_axis_turns

    !> Turns of 60 and 120 degrees either way, from e1 with the reference e2. The targets are
    !| (+-1/2, +-sqrt(3)/2): their dot product with e1 is +-1/2, exact, so the angles are
    !| acos(1/2) = PI/3 and acos(-1/2) = 2*PI/3, and the sign is that of the second component.
    !| And 45 degrees, (1, 1)/sqrt(2), whose cosine carries sqrt's rounding: acos's slope there is
    !| -sqrt(2), so the angle is PI/4 within an ulp or two.
    subroutine test_clock_hand_angle_between_vectors_exact_cosines()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle, half_sqrt3, half_sqrt2
        integer(int32) :: ierr, i_case
        real(real64) :: targets(2, 5), expected(5)

        half_sqrt3 = sqrt(3.0_real64)/2
        half_sqrt2 = sqrt(2.0_real64)/2
        v1 = [1.0_real64, 0.0_real64]
        reference = [0.0_real64, 1.0_real64]
        ! (+-1/2, +-sqrt(3)/2) and (1, 1)/sqrt(2), element by element: a constructor of
        ! variables would be an array temporary
        targets(1, 1:4) = [0.5_real64, 0.5_real64, -0.5_real64, -0.5_real64]
        targets(2, 1) = half_sqrt3
        targets(2, 2) = -half_sqrt3
        targets(2, 3) = half_sqrt3
        targets(2, 4) = -half_sqrt3
        targets(:, 5) = half_sqrt2
        expected = [PI/3, -PI/3, 2*PI/3, -2*PI/3, PI/4]

        do i_case = 1, size(expected)
            v2 = targets(:, i_case)
            call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_exact_cosines: ierr")
            call assert_equal_real(signed_angle, expected(i_case), ANGLE_TOL, &
                                   "test_clock_hand_angle_between_vectors_exact_cosines: 60/120/45 degree turns")
        end do
    end subroutine test_clock_hand_angle_between_vectors_exact_cosines

    !> The reference is what picks the orientation: in seven dimensions a turn from e3 to e5,
    !| oriented by e5, is +PI/2, and the same turn oriented by -e5 is -PI/2.
    subroutine test_clock_hand_angle_between_vectors_reversed_reference()
        real(real64) :: v1(7), v2(7), reference(7), signed_angle
        integer(int32) :: ierr

        v1 = 0.0_real64
        v1(3) = 1.0_real64
        v2 = 0.0_real64
        v2(5) = 1.0_real64

        reference = 0.0_real64
        reference(5) = 1.0_real64
        call clock_hand_angle_between_vectors(v1, v2, 7, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_reversed_reference: e5 ierr")
        call assert_equal_real(signed_angle, PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_reversed_reference: e5 makes e3 -> e5 positive")

        reference(5) = -1.0_real64
        call clock_hand_angle_between_vectors(v1, v2, 7, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_reversed_reference: -e5 ierr")
        call assert_equal_real(signed_angle, -PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_reversed_reference: -e5 makes it negative")
    end subroutine test_clock_hand_angle_between_vectors_reversed_reference

    !> Only the reference's component along the rotation, v2 - (v1 . v2) v1, counts: neither its
    !| length nor its part along v1 matters. From e1 to (1/2, sqrt(3)/2) that direction is
    !| (0, sqrt(3)/2), so (5, 3) -- mostly along v1 -- orients the turn positively, (5, -3)
    !| negatively, and the target itself, whose component is 1 - 1/4 > 0, positively.
    subroutine test_clock_hand_angle_between_vectors_reference_part_along_v1()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle
        integer(int32) :: ierr

        v1 = [1.0_real64, 0.0_real64]
        v2 = [0.5_real64, sqrt(3.0_real64)/2]

        reference = [5.0_real64, 3.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_reference_part_along_v1: (5, 3)")
        call assert_equal_real(signed_angle, PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_reference_part_along_v1: (5, 3) gives +PI/3")

        reference = [5.0_real64, -3.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_clock_hand_angle_between_vectors_reference_part_along_v1: (5, -3)")
        call assert_equal_real(signed_angle, -PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_reference_part_along_v1: (5, -3) gives -PI/3")

        reference = v2
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_reference_part_along_v1: v2")
        call assert_equal_real(signed_angle, PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_reference_part_along_v1: v2 gives +PI/3")
    end subroutine test_clock_hand_angle_between_vectors_reference_part_along_v1

    !> With the reference each first vector turned a quarter turn counter-clockwise (-y, x), the
    !| sign is the determinant's, so swapping the vectors negates the angle. For v1 = e1 and
    !| v2 = (1/2, sqrt(3)/2): PI/3 one way; the other way the reference is (-sqrt(3)/2, 1/2), the
    !| direction e1 - v2/2 = (3/4, -sqrt(3)/4), their product -sqrt(3)/2 < 0, so -PI/3. Likewise
    !| for (-1/2, sqrt(3)/2): 2*PI/3 and -2*PI/3.
    subroutine test_clock_hand_angle_between_vectors_anticommutative()
        real(real64) :: v1(2), v2(2), reference(2), angle_12, angle_21
        real(real64) :: targets(2, 2), expected(2)
        integer(int32) :: ierr_12, ierr_21, i_case

        targets(:, 1) = [0.5_real64, sqrt(3.0_real64)/2]
        targets(:, 2) = [-0.5_real64, sqrt(3.0_real64)/2]
        expected = [PI/3, 2*PI/3]
        v1 = [1.0_real64, 0.0_real64]

        do i_case = 1, size(expected)
            v2 = targets(:, i_case)
            reference(1) = -v1(2)
            reference(2) = v1(1)
            call clock_hand_angle_between_vectors(v1, v2, 2, reference, angle_12, ierr_12)
            reference(1) = -v2(2)
            reference(2) = v2(1)
            call clock_hand_angle_between_vectors(v2, v1, 2, reference, angle_21, ierr_21)
            call assert_equal_int(get_err_code(ierr_12), ERR_OK, "test_clock_hand_angle_between_vectors_anticommutative: 1 -> 2")
            call assert_equal_int(get_err_code(ierr_21), ERR_OK, "test_clock_hand_angle_between_vectors_anticommutative: 2 -> 1")
            call assert_equal_real(angle_12, expected(i_case), ANGLE_TOL, &
                                   "test_clock_hand_angle_between_vectors_anticommutative: v1 -> v2")
            call assert_equal_real(angle_21, -expected(i_case), ANGLE_TOL, &
                                   "test_clock_hand_angle_between_vectors_anticommutative: v2 -> v1 is its negative")
        end do
    end subroutine test_clock_hand_angle_between_vectors_anticommutative

    !> The sign is the reference's component along the direction the turn leaves v1 in, as the
    !| contract states it, so one fixed reference does not always flip sign when the turn is
    !| reversed. With the reference e2: e1 -> (-1/2, sqrt(3)/2) leaves along (0, sqrt(3)/2), +2*PI/3;
    !| the reverse leaves along e1 + v1/2 = (3/4, sqrt(3)/4), whose e2 component is positive too,
    !| so it is +2*PI/3 as well.
    subroutine test_clock_hand_angle_between_vectors_sign_taken_at_v1()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle
        integer(int32) :: ierr

        reference = [0.0_real64, 1.0_real64]

        v1 = [1.0_real64, 0.0_real64]
        v2 = [-0.5_real64, sqrt(3.0_real64)/2]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_sign_taken_at_v1: forward")
        call assert_equal_real(signed_angle, 2*PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_sign_taken_at_v1: forward is +2*PI/3")

        v1 = [-0.5_real64, sqrt(3.0_real64)/2]
        v2 = [1.0_real64, 0.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_sign_taken_at_v1: reverse")
        call assert_equal_real(signed_angle, 2*PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_sign_taken_at_v1: reverse is +2*PI/3 too")
    end subroutine test_clock_hand_angle_between_vectors_sign_taken_at_v1

    !> What the angle is for: three tissues, projected onto the RAP. The unit projections of the
    !| axes, a = (2, -1, -1)/sqrt(6), b = (-1, 2, -1)/sqrt(6) and c = (-1, -1, 2)/sqrt(6), are
    !| 120 degrees apart: a . b = (-2 - 2 + 1)/6 = -1/2. Oriented by b, the turn a -> b is
    !| +2*PI/3. Oriented by c it is negative: the direction a -> b leaves a in is
    !| b + a/2 = (0, 3/2, -3/2)/sqrt(6), and c's component along it is (-3/2 - 3)/6 < 0.
    subroutine test_clock_hand_angle_between_vectors_rap_three_tissues()
        real(real64) :: axis_a(3), axis_b(3), axis_c(3), signed_angle, root6
        integer(int32) :: ierr

        root6 = sqrt(6.0_real64)
        axis_a = [2.0_real64, -1.0_real64, -1.0_real64]/root6
        axis_b = [-1.0_real64, 2.0_real64, -1.0_real64]/root6
        axis_c = [-1.0_real64, -1.0_real64, 2.0_real64]/root6

        call clock_hand_angle_between_vectors(axis_a, axis_b, 3, axis_b, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_rap_three_tissues: by b")
        call assert_equal_real(signed_angle, 2*PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_rap_three_tissues: a -> b oriented by b")

        call clock_hand_angle_between_vectors(axis_a, axis_b, 3, axis_c, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_rap_three_tissues: by c")
        call assert_equal_real(signed_angle, -2*PI/3, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_rap_three_tissues: a -> b oriented by c")
    end subroutine test_clock_hand_angle_between_vectors_rap_three_tissues

    !> Mixed signs, an angle that is no fraction of PI: (1, -2, 3)/sqrt(14) and (-2, 1, -3)/sqrt(14)
    !| have the cosine (-2 - 2 - 9)/14 = -13/14. Oriented by v2 the turn is positive (v2's component
    !| along it is 1 - (13/14)**2 > 0), by -v2 negative.
    subroutine test_clock_hand_angle_between_vectors_mixed_signs()
        real(real64) :: v1(3), v2(3), reference(3), signed_angle, root14
        integer(int32) :: ierr

        root14 = sqrt(14.0_real64)
        v1 = [1.0_real64, -2.0_real64, 3.0_real64]/root14
        v2 = [-2.0_real64, 1.0_real64, -3.0_real64]/root14

        reference = v2
        call clock_hand_angle_between_vectors(v1, v2, 3, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_mixed_signs: by v2")
        call assert_equal_real(signed_angle, acos(-13.0_real64/14.0_real64), ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_mixed_signs: acos(-13/14), oriented by v2")

        reference = -v2
        call clock_hand_angle_between_vectors(v1, v2, 3, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_mixed_signs: by -v2")
        call assert_equal_real(signed_angle, -acos(-13.0_real64/14.0_real64), ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_mixed_signs: -acos(-13/14), oriented by -v2")
    end subroutine test_clock_hand_angle_between_vectors_mixed_signs

    !> Parallel and antiparallel vectors turn by 0 or PI, which has no side, so the angle comes back
    !| unsigned whatever the reference -- even a reference that would orient nothing, the zero
    !| vector: an antiparallel pair is +PI under e2, -e2 and 0.
    subroutine test_clock_hand_angle_between_vectors_parallel_unsigned()
        real(real64) :: v1(3), v2(3), reference(3), references(3, 3), signed_angle
        integer(int32) :: ierr, i_reference

        references(:, 1) = [0.0_real64, 1.0_real64, 0.0_real64]
        references(:, 2) = [0.0_real64, -1.0_real64, 0.0_real64]
        references(:, 3) = 0.0_real64
        v1 = [0.0_real64, 0.0_real64, 1.0_real64]

        do i_reference = 1, size(references, 2)
            reference = references(:, i_reference)

            v2 = v1
            call clock_hand_angle_between_vectors(v1, v2, 3, reference, signed_angle, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_parallel_unsigned: parallel")
            call assert_equal_real(signed_angle, 0.0_real64, 0.0_real64, &
                                   "test_clock_hand_angle_between_vectors_parallel_unsigned: parallel is 0")

            v2 = -v1
            call clock_hand_angle_between_vectors(v1, v2, 3, reference, signed_angle, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, &
                                  "test_clock_hand_angle_between_vectors_parallel_unsigned: antiparallel")
            call assert_equal_real(signed_angle, PI, ANGLE_TOL, &
                                   "test_clock_hand_angle_between_vectors_parallel_unsigned: antiparallel is +PI")
        end do
    end subroutine test_clock_hand_angle_between_vectors_parallel_unsigned

    !> A reference with no component along the rotation orients nothing: ERR_INVALID_INPUT, and no
    !| argument blamed, since it is the pair of reference and rotation that fails, not one argument.
    !| For the turn e1 -> e2 in five dimensions: the zero vector, v1 itself (its component along
    !| the direction v2 - (v1 . v2) v1, perpendicular to v1, is 0), and e3, normal to the plane.
    subroutine test_clock_hand_angle_between_vectors_degenerate_references()
        real(real64) :: v1(5), v2(5), reference(5), references(5, 3), signed_angle
        integer(int32) :: ierr, i_reference

        v1 = 0.0_real64
        v1(1) = 1.0_real64
        v2 = 0.0_real64
        v2(2) = 1.0_real64
        references = 0.0_real64
        references(:, 2) = v1
        references(3, 3) = 1.0_real64

        do i_reference = 1, size(references, 2)
            reference = references(:, i_reference)
            call clock_hand_angle_between_vectors(v1, v2, 5, reference, signed_angle, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, &
                            "test_clock_hand_angle_between_vectors_degenerate_references: the zero vector, v1, e3", 0_int32)
        end do
    end subroutine test_clock_hand_angle_between_vectors_degenerate_references

    !> n_dims on both sides of its bound: 0 is ERR_EMPTY_INPUT and -1 ERR_INVALID_INPUT, both
    !| blaming argument 3; 1 is valid, and in one dimension every pair is parallel or antiparallel:
    !| (1) -> (1) is 0 and (1) -> (-1) is PI.
    subroutine test_clock_hand_angle_between_vectors_dimensions()
        real(real64) :: v1(1), v2(1), reference(1), signed_angle
        integer(int32) :: ierr

        v1 = [1.0_real64]
        v2 = [-1.0_real64]
        reference = [1.0_real64]

        call clock_hand_angle_between_vectors(v1, v2, 0, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_clock_hand_angle_between_vectors_dimensions: n_dims = 0", 3_int32)
        call clock_hand_angle_between_vectors(v1, v2, -1, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_clock_hand_angle_between_vectors_dimensions: n_dims = -1", 3_int32)

        call clock_hand_angle_between_vectors(v1, v2, 1, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_dimensions: n_dims = 1")
        call assert_equal_real(signed_angle, PI, ANGLE_TOL, "test_clock_hand_angle_between_vectors_dimensions: (1) -> (-1)")
        call clock_hand_angle_between_vectors(v1, v1, 1, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_dimensions: n_dims = 1, same")
        call assert_equal_real(signed_angle, 0.0_real64, 0.0_real64, "test_clock_hand_angle_between_vectors_dimensions: (1) -> (1)")
    end subroutine test_clock_hand_angle_between_vectors_dimensions

    !> NaN, +Inf and -Inf in v1, v2 or the reference are ERR_NAN_INF, blaming that argument.
    subroutine test_clock_hand_angle_between_vectors_rejects_nan_and_inf()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle, bad(3)
        integer(int32) :: ierr, i_bad

        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        do i_bad = 1, size(bad)
            v1 = [1.0_real64, 0.0_real64]
            v2 = [0.0_real64, 1.0_real64]
            reference = [0.0_real64, 1.0_real64]

            v1(2) = bad(i_bad)
            call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_clock_hand_angle_between_vectors_rejects_nan_and_inf: in v1", 1_int32)
            v1(2) = 0.0_real64

            v2(1) = bad(i_bad)
            call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_clock_hand_angle_between_vectors_rejects_nan_and_inf: in v2", 2_int32)
            v2(1) = 0.0_real64

            reference(1) = bad(i_bad)
            call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_clock_hand_angle_between_vectors_rejects_nan_and_inf: in the reference", 4_int32)
        end do
    end subroutine test_clock_hand_angle_between_vectors_rejects_nan_and_inf

    !> A vector and itself turn by 0. (1, 1)/sqrt(2) is as normalized as a double can be, but its
    !| dot product with itself rounds to 1 - 2**-52, in any order of summation, with or without
    !| a fused multiply-add.
    subroutine test_clock_hand_angle_between_vectors_identical_diagonal()
        real(real64) :: v1(2), reference(2), signed_angle
        integer(int32) :: ierr

        v1 = [1.0_real64, 1.0_real64]/sqrt(2.0_real64)
        reference(1) = -v1(2)
        reference(2) = v1(1)

        call clock_hand_angle_between_vectors(v1, v1, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_identical_diagonal: ierr")
        ! Regression: acos(v1 . v2), ill-conditioned at 1 where its slope is infinite, turned the
        ! dot product's last-bit rounding, 1 - 2**-52, into 2.1e-8 rad, where the answer is 0. The
        ! atan2 form has no such loss (the part perpendicular to v1 is 2e-16 here).
        call assert_equal_real(signed_angle, 0.0_real64, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_identical_diagonal: a vector and itself")
    end subroutine test_clock_hand_angle_between_vectors_identical_diagonal

    !> A small clockwise turn, by y = 2**-22 rad (2.4e-7), from e1 to (sqrt(1 - y**2), -y) with the
    !| reference e2, is -asin(y). The tolerance is acos's again: the cosine sqrt(1 - y**2) is
    !| rounded by up to 2**-54, which acos's slope 1/sin(y) = 2**22 turns into 3e-10 rad -- still
    !| three orders of magnitude below the 4.8e-7 between the right sign and the wrong one.
    subroutine test_clock_hand_angle_between_vectors_small_turn_sign()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle, turn
        integer(int32) :: ierr

        turn = 2.0_real64**(-22)
        v1 = [1.0_real64, 0.0_real64]
        v2(1) = sqrt(1.0_real64 - turn**2)
        v2(2) = -turn
        reference = [0.0_real64, 1.0_real64]

        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_small_turn_sign: ierr")
        ! Regression: the parallel check used to compare sin(angle)**2, the squared length of v2's
        ! part perpendicular to v1, with is_close's absolute floor 1e-12, so every turn below
        ! 1e-6 rad counted as parallel and came back unsigned: this one as +2**-22.
        call assert_equal_real(signed_angle, -asin(turn), 1e-9_real64, &
                               "test_clock_hand_angle_between_vectors_small_turn_sign: a clockwise turn of 2**-22 rad")
    end subroutine test_clock_hand_angle_between_vectors_small_turn_sign

    !> The sign is that of the reference's component along the rotation, and a short reference has
    !| one as surely as a long one: e1 -> e2 oriented by (0, 1e-13) is +PI/2, as by e2.
    subroutine test_clock_hand_angle_between_vectors_short_reference()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle
        integer(int32) :: ierr

        v1 = [1.0_real64, 0.0_real64]
        v2 = [0.0_real64, 1.0_real64]
        reference = [0.0_real64, 1e-13_real64]

        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        ! Regression: that component, 1e-13 here, used to be compared with is_close's absolute
        ! floor 1e-12, so the reference's length decided whether it oriented anything, and this
        ! one was reported as orienting nothing (ERR_INVALID_INPUT).
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angle_between_vectors_short_reference: ierr")
        call assert_equal_real(signed_angle, PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_short_reference: oriented by (0, 1e-13)")
    end subroutine test_clock_hand_angle_between_vectors_short_reference


    !> A zero vector has no direction to turn from or to: ERR_DIVISION_BY_ZERO, whether it is v1 or
    !| v2 and whatever the reference. It used to come back as PI/2, with success.
    subroutine test_clock_hand_angle_between_vectors_zero_vector()
        real(real64) :: v1(2), v2(2), zero(2), reference(2), signed_angle
        integer(int32) :: ierr

        v1 = [1.0_real64, 0.0_real64]
        v2 = [0.0_real64, 1.0_real64]
        zero = 0.0_real64
        reference = [0.0_real64, 1.0_real64]

        call clock_hand_angle_between_vectors(zero, v2, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_clock_hand_angle_between_vectors_zero_vector: v1 = 0")
        call clock_hand_angle_between_vectors(v1, zero, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_clock_hand_angle_between_vectors_zero_vector: v2 = 0")
        call clock_hand_angle_between_vectors(zero, zero, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_clock_hand_angle_between_vectors_zero_vector: both 0")
    end subroutine test_clock_hand_angle_between_vectors_zero_vector

    !> The angle between two vectors does not depend on their lengths, so they need not be
    !| normalized: (2, 0) -> (1, 1) is PI/4 and (3, 0) -> (0, -5) is -PI/2, both oriented by e2.
    !| The first used to come back as 0, from acos of the dot product 2.
    subroutine test_clock_hand_angle_between_vectors_unnormalized()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle
        integer(int32) :: ierr

        reference = [0.0_real64, 1.0_real64]

        v1 = [2.0_real64, 0.0_real64]
        v2 = [1.0_real64, 1.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_OK, "test_clock_hand_angle_between_vectors_unnormalized: ierr for (2, 0) -> (1, 1)")
        call assert_equal_real(signed_angle, PI/4, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_unnormalized: (2, 0) -> (1, 1) is PI/4")

        v1 = [3.0_real64, 0.0_real64]
        v2 = [0.0_real64, -5.0_real64]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_OK, "test_clock_hand_angle_between_vectors_unnormalized: ierr for (3, 0) -> (0, -5)")
        call assert_equal_real(signed_angle, -PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_unnormalized: (3, 0) -> (0, -5) is -PI/2")
    end subroutine test_clock_hand_angle_between_vectors_unnormalized

    !> Neither a huge nor a tiny length may change an angle: (huge, huge) -> (-huge, huge) is PI/2,
    !| where the dot product used to overflow to Inf - Inf = NaN, and (2**-600, 0) ->
    !| (2**-600, 2**-600) is PI/4, where the dot product 2**-1200 used to underflow to 0 and give PI/2.
    subroutine test_clock_hand_angle_between_vectors_extreme_magnitudes()
        real(real64) :: v1(2), v2(2), reference(2), signed_angle
        integer(int32) :: ierr

        reference = [0.0_real64, 1.0_real64]

        v1 = [huge(1.0_real64), huge(1.0_real64)]
        v2 = [-huge(1.0_real64), huge(1.0_real64)]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_OK, "test_clock_hand_angle_between_vectors_extreme_magnitudes: ierr for huge")
        call assert_equal_real(signed_angle, PI/2, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_extreme_magnitudes: (huge, huge) -> (-huge, huge)")

        v1 = [2.0_real64**(-600), 0.0_real64]
        v2 = [2.0_real64**(-600), 2.0_real64**(-600)]
        call clock_hand_angle_between_vectors(v1, v2, 2, reference, signed_angle, ierr)
        call assert_err(ierr, ERR_OK, "test_clock_hand_angle_between_vectors_extreme_magnitudes: ierr for tiny")
        call assert_equal_real(signed_angle, PI/4, ANGLE_TOL, &
                               "test_clock_hand_angle_between_vectors_extreme_magnitudes: (2**-600, 0) -> (2**-600, 2**-600)")
    end subroutine test_clock_hand_angle_between_vectors_extreme_magnitudes

end module mod_test_clock_hand_angle_between_vectors
