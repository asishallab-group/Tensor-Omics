!> The `clock_hand_angles_for_shift_vectors` cases: hand-derived angles for a batch of fields
!| sharing one reference, the selection mask (unselected fields are neither computed nor
!| counted, and the results are packed in field order), a reference that fails to orient one
!| field failing the call, each dimension's bound, the mask count, and the rejection of NaN and
!| Inf -- in unselected fields too.
!|
!| How a single angle and its sign are derived is `clock_hand_angle_between_vectors`'s suite's
!| business (both share one compute core); these cases use turns with exact cosines only, and
!| test what the batch adds.
module mod_test_clock_hand_angles_for_shift_vectors
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use f42_math_impl, only: PI
    use tox_relative_axis_plane_tools
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

    real(real64), parameter :: ANGLE_TOL = 4e-15_real64
        !! A few ulps of an angle up to PI: acos rounds its result, and PI/2 and PI/3 round the
        !! division.

contains

    !> Get array of all available tests.
    function get_all_tests_clock_hand_angles_for_shift_vectors() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(8))

        all_tests(1) = test_case("test_clock_hand_angles_for_shift_vectors_axis_turns", &
                                 test_clock_hand_angles_for_shift_vectors_axis_turns)
        all_tests(2) = test_case("test_clock_hand_angles_for_shift_vectors_one_reference", &
                                 test_clock_hand_angles_for_shift_vectors_one_reference)
        all_tests(3) = test_case("test_clock_hand_angles_for_shift_vectors_selection_mask", &
                                 test_clock_hand_angles_for_shift_vectors_selection_mask)
        all_tests(4) = test_case("test_clock_hand_angles_for_shift_vectors_unoriented_field", &
                                 test_clock_hand_angles_for_shift_vectors_unoriented_field)
        all_tests(5) = test_case("test_clock_hand_angles_for_shift_vectors_dimensions", &
                                 test_clock_hand_angles_for_shift_vectors_dimensions)
        all_tests(6) = test_case("test_clock_hand_angles_for_shift_vectors_mask_count", &
                                 test_clock_hand_angles_for_shift_vectors_mask_count)
        all_tests(7) = test_case("test_clock_hand_angles_for_shift_vectors_rejects_nan_inf", &
                                 test_clock_hand_angles_for_shift_vectors_rejects_nan_inf)
        all_tests(8) = test_case("test_clock_hand_angles_for_shift_vectors_zero_vector", &
                                 test_clock_hand_angles_for_shift_vectors_zero_vector)
    end function get_all_tests_clock_hand_angles_for_shift_vectors

    !> Three fields from e1, with the reference e2 (counter-clockwise positive): to e2 is +PI/2,
    !| to -e1 PI (antiparallel, so unsigned), to -e2 -PI/2. Each field is angled from its origin
    !| (first vector) to its target (second).
    subroutine test_clock_hand_angles_for_shift_vectors_axis_turns()
        real(real64) :: fields(2, 2, 3), reference(2), signed_angles(3), expected(3)
        logical(c_bool) :: fields_selection_mask(3)
        integer(int32) :: ierr

        fields(:, 1, 1) = [1.0_real64, 0.0_real64]
        fields(:, 2, 1) = [0.0_real64, 1.0_real64]
        fields(:, 1, 2) = [1.0_real64, 0.0_real64]
        fields(:, 2, 2) = [-1.0_real64, 0.0_real64]
        fields(:, 1, 3) = [1.0_real64, 0.0_real64]
        fields(:, 2, 3) = [0.0_real64, -1.0_real64]
        fields_selection_mask = .true.
        reference = [0.0_real64, 1.0_real64]
        expected = [PI/2, PI, -PI/2]

        call clock_hand_angles_for_shift_vectors(fields, 2, 3, fields_selection_mask, 3, reference, signed_angles, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angles_for_shift_vectors_axis_turns: ierr")
        call assert_equal_array_real(signed_angles, expected, 3, ANGLE_TOL, &
                                     "test_clock_hand_angles_for_shift_vectors_axis_turns: +PI/2, PI, -PI/2")
    end subroutine test_clock_hand_angles_for_shift_vectors_axis_turns

    !> One reference orients fields turning in different planes. In three dimensions with the
    !| reference e3: e1 -> e3 is +PI/2 (the turn leaves e1 along e3); e2 -> -e3 is -PI/2 (along
    !| -e3); e1 -> (1/2, 0, -sqrt(3)/2) has the exact cosine 1/2 and leaves along -e3, so -PI/3;
    !| and e3 -> -e3 is antiparallel, so PI unsigned.
    subroutine test_clock_hand_angles_for_shift_vectors_one_reference()
        real(real64) :: fields(3, 2, 4), reference(3), signed_angles(4), expected(4)
        logical(c_bool) :: fields_selection_mask(4)
        integer(int32) :: ierr

        fields(:, 1, 1) = [1.0_real64, 0.0_real64, 0.0_real64]
        fields(:, 2, 1) = [0.0_real64, 0.0_real64, 1.0_real64]
        fields(:, 1, 2) = [0.0_real64, 1.0_real64, 0.0_real64]
        fields(:, 2, 2) = [0.0_real64, 0.0_real64, -1.0_real64]
        fields(:, 1, 3) = [1.0_real64, 0.0_real64, 0.0_real64]
        fields(:, 2, 3) = [0.5_real64, 0.0_real64, -sqrt(3.0_real64)/2]
        fields(:, 1, 4) = [0.0_real64, 0.0_real64, 1.0_real64]
        fields(:, 2, 4) = [0.0_real64, 0.0_real64, -1.0_real64]
        fields_selection_mask = .true.
        reference = [0.0_real64, 0.0_real64, 1.0_real64]
        expected = [PI/2, -PI/2, -PI/3, PI]

        call clock_hand_angles_for_shift_vectors(fields, 3, 4, fields_selection_mask, 4, reference, signed_angles, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angles_for_shift_vectors_one_reference: ierr")
        call assert_equal_array_real(signed_angles, expected, 4, ANGLE_TOL, &
                                     "test_clock_hand_angles_for_shift_vectors_one_reference: +PI/2, -PI/2, -PI/3, PI")
    end subroutine test_clock_hand_angles_for_shift_vectors_one_reference

    !> Only the selected fields 2 and 4 are angled, and their results fill signed_angles in field
    !| order: e1 -> -e1 is PI, e1 -> (1, 1)/sqrt(2) is PI/4 (cosine sqrt(2)/2, carrying sqrt's
    !| rounding, which acos's slope -sqrt(2) keeps within an ulp or two). Fields 1 and 3 turn e2
    !| to e1, which the reference e2 cannot orient (the turn leaves e2 along e1, and e2 . e1 = 0):
    !| had either been computed, the call would report ERR_INVALID_INPUT.
    subroutine test_clock_hand_angles_for_shift_vectors_selection_mask()
        real(real64) :: fields(2, 2, 4), reference(2), signed_angles(2), expected(2)
        logical(c_bool) :: fields_selection_mask(4)
        integer(int32) :: ierr

        fields(:, 1, 1) = [0.0_real64, 1.0_real64]
        fields(:, 2, 1) = [1.0_real64, 0.0_real64]
        fields(:, 1, 2) = [1.0_real64, 0.0_real64]
        fields(:, 2, 2) = [-1.0_real64, 0.0_real64]
        fields(:, 1, 3) = [0.0_real64, 1.0_real64]
        fields(:, 2, 3) = [1.0_real64, 0.0_real64]
        fields(:, 1, 4) = [1.0_real64, 0.0_real64]
        fields(:, 2, 4) = [1.0_real64, 1.0_real64]/sqrt(2.0_real64)
        fields_selection_mask = [.false., .true., .false., .true.]
        reference = [0.0_real64, 1.0_real64]
        expected = [PI, PI/4]

        call clock_hand_angles_for_shift_vectors(fields, 2, 4, fields_selection_mask, 2, reference, signed_angles, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_clock_hand_angles_for_shift_vectors_selection_mask: unselected fields are not angled")
        call assert_equal_array_real(signed_angles, expected, 2, ANGLE_TOL, &
                                     "test_clock_hand_angles_for_shift_vectors_selection_mask: fields 2 and 4, packed")
    end subroutine test_clock_hand_angles_for_shift_vectors_selection_mask

    !> A single selected field the reference cannot orient fails the call: fields 1 and 3 turn
    !| e1 -> e2, which the reference e2 orients, but field 2 turns e1 -> e3, normal to e2.
    !| ERR_INVALID_INPUT blames no argument: it is the reference with this rotation that fails.
    subroutine test_clock_hand_angles_for_shift_vectors_unoriented_field()
        real(real64) :: fields(3, 2, 3), reference(3), signed_angles(3)
        logical(c_bool) :: fields_selection_mask(3)
        integer(int32) :: ierr

        fields = 0.0_real64
        fields(1, 1, :) = 1.0_real64
        fields(2, 2, 1) = 1.0_real64
        fields(3, 2, 2) = 1.0_real64
        fields(2, 2, 3) = 1.0_real64
        fields_selection_mask = .true.
        reference = [0.0_real64, 1.0_real64, 0.0_real64]

        call clock_hand_angles_for_shift_vectors(fields, 3, 3, fields_selection_mask, 3, reference, signed_angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_clock_hand_angles_for_shift_vectors_unoriented_field: field 2 is normal to the reference", 0_int32)
    end subroutine test_clock_hand_angles_for_shift_vectors_unoriented_field

    !> Each dimension on its own, the others valid: 0 is ERR_EMPTY_INPUT and -1 ERR_INVALID_INPUT,
    !| blaming n_dims (argument 2), n_fields (3) or n_selected_fields (5). The last valid values,
    !| one dimension, one field, one selected: (1) -> (-1) is PI.
    subroutine test_clock_hand_angles_for_shift_vectors_dimensions()
        real(real64) :: fields(1, 2, 1), reference(1), signed_angles(1)
        logical(c_bool) :: fields_selection_mask(1)
        integer(int32) :: ierr, i_size, sizes(2), expected_codes(2)

        fields(:, 1, 1) = [1.0_real64]
        fields(:, 2, 1) = [-1.0_real64]
        fields_selection_mask = .true.
        reference = [1.0_real64]
        sizes = [0, -1]
        expected_codes = [ERR_EMPTY_INPUT, ERR_INVALID_INPUT]

        do i_size = 1, size(sizes)
            call clock_hand_angles_for_shift_vectors(fields, sizes(i_size), 1, fields_selection_mask, 1, reference, &
                                                     signed_angles, ierr)
            call assert_err(ierr, expected_codes(i_size), "test_clock_hand_angles_for_shift_vectors_dimensions: n_dims", 2_int32)
            call clock_hand_angles_for_shift_vectors(fields, 1, sizes(i_size), fields_selection_mask, 1, reference, &
                                                     signed_angles, ierr)
            call assert_err(ierr, expected_codes(i_size), "test_clock_hand_angles_for_shift_vectors_dimensions: n_fields", 3_int32)
            call clock_hand_angles_for_shift_vectors(fields, 1, 1, fields_selection_mask, sizes(i_size), reference, &
                                                     signed_angles, ierr)
            call assert_err(ierr, expected_codes(i_size), &
                            "test_clock_hand_angles_for_shift_vectors_dimensions: n_selected_fields", 5_int32)
        end do

        call clock_hand_angles_for_shift_vectors(fields, 1, 1, fields_selection_mask, 1, reference, signed_angles, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_clock_hand_angles_for_shift_vectors_dimensions: all 1")
        call assert_equal_real(signed_angles(1), PI, ANGLE_TOL, "test_clock_hand_angles_for_shift_vectors_dimensions: (1) -> (-1)")
    end subroutine test_clock_hand_angles_for_shift_vectors_dimensions

    !> n_selected_fields must be the mask's count of `.true.`: with two of two fields selected,
    !| 1 and 3 are ERR_INVALID_INPUT, blaming n_selected_fields (argument 5). Fewer would drop a
    !| selected field's result; more would leave one undefined.
    subroutine test_clock_hand_angles_for_shift_vectors_mask_count()
        real(real64) :: fields(2, 2, 2), reference(2), signed_angles(3)
        logical(c_bool) :: fields_selection_mask(2)
        integer(int32) :: ierr

        fields(:, 1, 1) = [1.0_real64, 0.0_real64]
        fields(:, 2, 1) = [0.0_real64, 1.0_real64]
        fields(:, :, 2) = fields(:, :, 1)
        fields_selection_mask = .true.
        reference = [0.0_real64, 1.0_real64]

        call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 1, reference, signed_angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_clock_hand_angles_for_shift_vectors_mask_count: 1 for 2 selected", 5_int32)
        call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 3, reference, signed_angles, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_clock_hand_angles_for_shift_vectors_mask_count: 3 for 2 selected", 5_int32)
    end subroutine test_clock_hand_angles_for_shift_vectors_mask_count

    !> NaN, +Inf and -Inf are ERR_NAN_INF: in a selected field's origin or target, blaming fields
    !| (argument 1); in an unselected field too, since TOX takes no NaN or Inf anywhere in its
    !| input; and in the reference, blaming it (argument 6).
    subroutine test_clock_hand_angles_for_shift_vectors_rejects_nan_inf()
        real(real64) :: fields(2, 2, 2), reference(2), signed_angles(1), bad(3)
        logical(c_bool) :: fields_selection_mask(2)
        integer(int32) :: ierr, i_bad

        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)
        fields_selection_mask = [.true., .false.]

        do i_bad = 1, size(bad)
            fields(:, 1, 1) = [1.0_real64, 0.0_real64]
            fields(:, 2, 1) = [0.0_real64, 1.0_real64]
            fields(:, :, 2) = fields(:, :, 1)
            reference = [0.0_real64, 1.0_real64]

            fields(2, 1, 1) = bad(i_bad)
            call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 1, reference, signed_angles, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_clock_hand_angles_for_shift_vectors_rejects_nan_inf: origin", 1_int32)
            fields(2, 1, 1) = 0.0_real64

            fields(1, 2, 1) = bad(i_bad)
            call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 1, reference, signed_angles, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_clock_hand_angles_for_shift_vectors_rejects_nan_inf: target", 1_int32)
            fields(1, 2, 1) = 0.0_real64

            fields(1, 1, 2) = bad(i_bad)
            call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 1, reference, signed_angles, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_clock_hand_angles_for_shift_vectors_rejects_nan_inf: an unselected field", 1_int32)
            fields(1, 1, 2) = 1.0_real64

            reference(2) = bad(i_bad)
            call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 1, reference, signed_angles, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_clock_hand_angles_for_shift_vectors_rejects_nan_inf: reference", 6_int32)
        end do
    end subroutine test_clock_hand_angles_for_shift_vectors_rejects_nan_inf


    !> A selected field with a zero origin or target has no turn to angle and fails the call with
    !| ERR_DIVISION_BY_ZERO; the same field unselected is not angled, so the call succeeds, and
    !| the selected field 2, e1 -> e2 under the reference e2, is +PI/2.
    subroutine test_clock_hand_angles_for_shift_vectors_zero_vector()
        real(real64) :: fields(2, 2, 2), reference(2), signed_angles(2), one_angle(1)
        logical(c_bool) :: fields_selection_mask(2)
        integer(int32) :: ierr

        fields(:, 1, 1) = 0.0_real64
        fields(:, 2, 1) = [0.0_real64, 1.0_real64]
        fields(:, 1, 2) = [1.0_real64, 0.0_real64]
        fields(:, 2, 2) = [0.0_real64, 1.0_real64]
        reference = [0.0_real64, 1.0_real64]

        fields_selection_mask = .true.
        call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 2, reference, signed_angles, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_clock_hand_angles_for_shift_vectors_zero_vector: zero origin")

        fields(:, 1, 1) = [1.0_real64, 0.0_real64]
        fields(:, 2, 1) = 0.0_real64
        call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 2, reference, signed_angles, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_clock_hand_angles_for_shift_vectors_zero_vector: zero target")

        fields_selection_mask = [.false., .true.]
        call clock_hand_angles_for_shift_vectors(fields, 2, 2, fields_selection_mask, 1, reference, one_angle, ierr)
        call assert_err(ierr, ERR_OK, "test_clock_hand_angles_for_shift_vectors_zero_vector: zero field unselected")
        call assert_equal_real(one_angle(1), PI/2, ANGLE_TOL, &
                               "test_clock_hand_angles_for_shift_vectors_zero_vector: e1 -> e2 is +PI/2")
    end subroutine test_clock_hand_angles_for_shift_vectors_zero_vector

end module mod_test_clock_hand_angles_for_shift_vectors
