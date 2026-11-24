!> The `euclidean_distance` cases: the distance `sqrt(sum((vec1 - vec2)**2))` of Pythagorean
!| inputs in one to four dimensions, a single element, identical and zero vectors, many elements,
!| scales far from 1, cancellation next to 1, the largest finite magnitudes, the dimension check
!| and the NaN/Inf checks.
!|
!| Every expected value is written out by hand. The differences are integers, or integers times a
!| power of two, chosen as Pythagorean triples and quadruples, so each sum of squares is a perfect
!| square, the root is exact and every comparison has tolerance 0. Every call is made in both
!| orders, so each case checks that the distance is symmetric as well.
module mod_test_euclidean_distance
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use tox_euclidean_distance, only: euclidean_distance
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_euclidean_distance_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    implicit none
    private
    public :: get_all_tests_euclidean_distance

contains

    !> Every case of `euclidean_distance`.
    function get_all_tests_euclidean_distance() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(10))
        all_tests(1) = test_case("test_euclidean_distance_pythagorean", test_euclidean_distance_pythagorean)
        all_tests(2) = test_case("test_euclidean_distance_one_element", test_euclidean_distance_one_element)
        all_tests(3) = test_case("test_euclidean_distance_identical_vectors", test_euclidean_distance_identical_vectors)
        all_tests(4) = test_case("test_euclidean_distance_zero_vectors", test_euclidean_distance_zero_vectors)
        all_tests(5) = test_case("test_euclidean_distance_many_elements", test_euclidean_distance_many_elements)
        all_tests(6) = test_case("test_euclidean_distance_power_of_two_scales", test_euclidean_distance_power_of_two_scales)
        all_tests(7) = test_case("test_euclidean_distance_cancellation_next_to_one", &
                                 test_euclidean_distance_cancellation_next_to_one)
        all_tests(8) = test_case("test_euclidean_distance_largest_magnitudes", test_euclidean_distance_largest_magnitudes)
        all_tests(9) = test_case("test_euclidean_distance_dimension", test_euclidean_distance_dimension)
        all_tests(10) = test_case("test_euclidean_distance_rejects_nan_and_inf", test_euclidean_distance_rejects_nan_and_inf)
    end function get_all_tests_euclidean_distance

    !> Calls the procedure on (`vec1`, `vec2`) and on (`vec2`, `vec1`), and checks `ierr` OK and
    !| `result` exactly `expected` both times. `result` starts at UNWRITTEN, so a call that never
    !| writes it shows.
    subroutine check_distance(n_elements, vec1, vec2, expected, msg)
        integer(int32), intent(in) :: n_elements
        real(real64), intent(in) :: vec1(n_elements), vec2(n_elements), expected
        character(*), intent(in) :: msg
        real(real64) :: result
        integer(int32) :: ierr

        result = UNWRITTEN
        call euclidean_distance(vec1, vec2, n_elements, result, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_real(result, expected, 0.0_real64, msg//": result")

        result = UNWRITTEN
        call euclidean_distance(vec2, vec1, n_elements, result, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr, vectors swapped")
        call assert_equal_real(result, expected, 0.0_real64, msg//": result, vectors swapped")
    end subroutine check_distance

    !> Differences that are Pythagorean triples and quadruples, so the distance is an integer:
    !| - kept from the old suite (`test_euclidean_distance_2d_origin`): [3, 4] to the origin is 5;
    !| - kept from the old suite (`test_euclidean_distance_mixed_signs`): [-3, 4] to [3, -4]
    !|   differs by [-6, 8], distance 10;
    !| - in place of the old `test_euclidean_distance_3d` ([1, 2, 3] to [4, 5, 6], whose expected
    !|   value sqrt(27) is not exact): [1, 2, 3] to [3, 5, 9] differs by [2, 3, 6], distance 7;
    !| - four dimensions: [1, 2, 3, 4] to [3, 6, 8, 10] differs by [2, 4, 5, 6], distance 9.
    subroutine test_euclidean_distance_pythagorean()
        real(real64) :: plane_a(2), plane_b(2), space_a(3), space_b(3), four_a(4), four_b(4)

        plane_a = [3.0_real64, 4.0_real64]
        plane_b = 0.0_real64
        call check_distance(2, plane_a, plane_b, 5.0_real64, "test_euclidean_distance_pythagorean: [3, 4] to the origin")

        plane_a = [-3.0_real64, 4.0_real64]
        plane_b = [3.0_real64, -4.0_real64]
        call check_distance(2, plane_a, plane_b, 10.0_real64, "test_euclidean_distance_pythagorean: mixed signs")

        space_a = [1.0_real64, 2.0_real64, 3.0_real64]
        space_b = [3.0_real64, 5.0_real64, 9.0_real64]
        call check_distance(3, space_a, space_b, 7.0_real64, "test_euclidean_distance_pythagorean: three dimensions")

        four_a = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        four_b = [3.0_real64, 6.0_real64, 8.0_real64, 10.0_real64]
        call check_distance(4, four_a, four_b, 9.0_real64, "test_euclidean_distance_pythagorean: four dimensions")
    end subroutine test_euclidean_distance_pythagorean

    !> Kept from the old suite (`test_euclidean_distance_1d`): the smallest valid size, one
    !| element, where the distance is the absolute difference. [5] to [2] is 3; [-2.5] to [4] is
    !| 6.5, both exact.
    subroutine test_euclidean_distance_one_element()
        real(real64) :: vec1(1), vec2(1)

        vec1 = 5.0_real64
        vec2 = 2.0_real64
        call check_distance(1, vec1, vec2, 3.0_real64, "test_euclidean_distance_one_element: [5] to [2]")

        vec1 = -2.5_real64
        vec2 = 4.0_real64
        call check_distance(1, vec1, vec2, 6.5_real64, "test_euclidean_distance_one_element: [-2.5] to [4]")
    end subroutine test_euclidean_distance_one_element

    !> Kept from the old suite (`test_euclidean_distance_identical`, whose tolerance was 1e-15):
    !| a vector is at distance exactly 0 from itself, whether its values are exact in binary or
    !| not ([1.5, 2.7, 3.9], as the old Python and R suites had it), and at the ends of the range
    !| ([huge, -huge, tiny]).
    subroutine test_euclidean_distance_identical_vectors()
        real(real64) :: vec1(3), vec2(3)

        vec1 = [1.5_real64, 2.7_real64, 3.9_real64]
        vec2 = vec1
        call check_distance(3, vec1, vec2, 0.0_real64, "test_euclidean_distance_identical_vectors: [1.5, 2.7, 3.9]")

        vec1 = [huge(1.0_real64), -huge(1.0_real64), tiny(1.0_real64)]
        vec2 = vec1
        call check_distance(3, vec1, vec2, 0.0_real64, "test_euclidean_distance_identical_vectors: [huge, -huge, tiny]")
    end subroutine test_euclidean_distance_identical_vectors

    !> Kept from the old suite (`test_euclidean_distance_zero_vectors`): two zero vectors are at
    !| distance 0, in three dimensions and in one.
    subroutine test_euclidean_distance_zero_vectors()
        real(real64) :: vec1(3), vec2(3), single1(1), single2(1)

        vec1 = 0.0_real64
        vec2 = 0.0_real64
        call check_distance(3, vec1, vec2, 0.0_real64, "test_euclidean_distance_zero_vectors: three dimensions")

        single1 = 0.0_real64
        single2 = 0.0_real64
        call check_distance(1, single1, single2, 0.0_real64, "test_euclidean_distance_zero_vectors: one element")
    end subroutine test_euclidean_distance_zero_vectors

    !> In place of the old `test_euclidean_distance_high_dim` (10 elements, sqrt(10)): 10,000
    !| elements, `vec2 = vec1 + 1` element by element, so the sum of squares is 10,000 and the
    !| distance exactly 100. A procedure that stops early, or reads one vector short, misses it.
    subroutine test_euclidean_distance_many_elements()
        integer(int32), parameter :: n_elements = 10000
        real(real64), allocatable :: vec1(:), vec2(:)
        integer(int32) :: i_element

        allocate (vec1(n_elements), vec2(n_elements))
        do i_element = 1, n_elements
            vec1(i_element) = real(i_element, real64)
            vec2(i_element) = real(i_element + 1, real64)
        end do
        call check_distance(n_elements, vec1, vec2, 100.0_real64, "test_euclidean_distance_many_elements")
    end subroutine test_euclidean_distance_many_elements

    !> In place of the old `test_euclidean_distance_small_numbers` (1e-15, tolerance 1e-27) and
    !| `test_euclidean_distance_large_numbers` (1e15, tolerance 1e3): the triple [3, 4] times
    !| 2**-50 and times 2**50 is at distance exactly 5*2**-50 and 5*2**50 from the origin. The
    !| squares stay far inside the normal range, so these hold with or without any scaling.
    subroutine test_euclidean_distance_power_of_two_scales()
        real(real64), parameter :: SMALL_SCALE = 2.0_real64**(-50), LARGE_SCALE = 2.0_real64**50
        real(real64) :: vec1(2), vec2(2)

        vec2 = 0.0_real64
        vec1 = [3.0_real64, 4.0_real64]*SMALL_SCALE
        call check_distance(2, vec1, vec2, 5.0_real64*SMALL_SCALE, "test_euclidean_distance_power_of_two_scales: 2**-50")

        vec1 = [3.0_real64, 4.0_real64]*LARGE_SCALE
        call check_distance(2, vec1, vec2, 5.0_real64*LARGE_SCALE, "test_euclidean_distance_power_of_two_scales: 2**50")
    end subroutine test_euclidean_distance_power_of_two_scales

    !> In place of the old `test_numerical_precision_epsilon` and `test_numerical_precision_close`
    !| (decimal values a few ulps apart, tolerance 10 eps and 1e-14): [1, 1] to
    !| [1 + 3 eps, 1 + 4 eps]. Both are representable, as eps is the spacing of the numbers just
    !| above 1, and their differences 3 eps and 4 eps are exact, so the distance is exactly 5 eps.
    subroutine test_euclidean_distance_cancellation_next_to_one()
        real(real64), parameter :: EPS = epsilon(1.0_real64)
        real(real64) :: vec1(2), vec2(2)

        vec1 = 1.0_real64
        vec2 = [1.0_real64 + 3.0_real64*EPS, 1.0_real64 + 4.0_real64*EPS]
        call check_distance(2, vec1, vec2, 5.0_real64*EPS, "test_euclidean_distance_cancellation_next_to_one")
    end subroutine test_euclidean_distance_cancellation_next_to_one

    !> The largest finite magnitudes whose difference needs no square past `huge`:
    !| - [huge] to [huge] is 0;
    !| - [huge] to [-huge] is 2*huge, past `huge`: the result is +Inf with `ierr` OK, which is
    !|   what f42's `norm` documents for a true length past `huge` (the square of a difference
    !|   that is itself +Inf). Neither refused nor clamped to `huge`.
    subroutine test_euclidean_distance_largest_magnitudes()
        real(real64), parameter :: HUGE_VALUE = huge(1.0_real64)
        real(real64) :: vec1(1), vec2(1)

        vec1 = HUGE_VALUE
        vec2 = HUGE_VALUE
        call check_distance(1, vec1, vec2, 0.0_real64, "test_euclidean_distance_largest_magnitudes: huge to huge")

        vec2 = -HUGE_VALUE
        call check_distance(1, vec1, vec2, ieee_value(1.0_real64, ieee_positive_inf), &
                            "test_euclidean_distance_largest_magnitudes: huge to -huge")
    end subroutine test_euclidean_distance_largest_magnitudes

    !> `n_elements` (argument 3) is ERR_EMPTY_INPUT at 0 (kept from the old Python and R suites'
    !| empty vectors) and ERR_INVALID_INPUT at -1. 1 is in
    !| `test_euclidean_distance_one_element`.
    subroutine test_euclidean_distance_dimension()
        real(real64) :: vec1(1), vec2(1), result
        integer(int32) :: ierr

        vec1 = 1.0_real64
        vec2 = 2.0_real64
        call euclidean_distance(vec1, vec2, 0, result, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_euclidean_distance_dimension: n_elements 0", 3)
        call euclidean_distance(vec1, vec2, -1, result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_euclidean_distance_dimension: n_elements -1", 3)
    end subroutine test_euclidean_distance_dimension

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `vec1` (argument 1) and in `vec2` (argument 2).
    !| Each sits in the last element, so the check must cover every element.
    subroutine test_euclidean_distance_rejects_nan_and_inf()
        integer(int32), parameter :: n_elements = 3
        real(real64) :: vec1(n_elements), vec2(n_elements), result
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            vec1 = 1.0_real64
            vec2 = 2.0_real64
            vec1(n_elements) = non_finite_value(i_value)
            call euclidean_distance(vec1, vec2, n_elements, result, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_euclidean_distance_rejects_nan_and_inf: vec1 " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)

            vec1 = 1.0_real64
            vec2(n_elements) = non_finite_value(i_value)
            call euclidean_distance(vec1, vec2, n_elements, result, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_euclidean_distance_rejects_nan_and_inf: vec2 " &
                            //trim(NON_FINITE_NAMES(i_value)), 2)
        end do
    end subroutine test_euclidean_distance_rejects_nan_and_inf

end module mod_test_euclidean_distance
