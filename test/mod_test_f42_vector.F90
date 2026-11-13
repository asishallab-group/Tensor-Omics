!> Unit test suite for the overflow-safe length and angle of f42_vector_impl: `scaled_length`,
!| `has_direction` and `angle_to_direction`.
!|
!| Expected values are derived by hand in the comments. Most are exact: the scaling is by a
!| power of two, so a scaled vector is the vector with a shifted exponent, and a vector and a
!| power-of-two multiple of it have the same unit vector to the last bit. A tolerance appears only
!| where `atan2` or the rounded `PI` stands between input and value. Every comparison is written so
!| that a NaN fails it.
module mod_test_f42_vector
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_get_underflow_mode
    use f42_math_impl, only: PI
    use f42_vector_impl, only: scaled_length, has_direction, angle_to_direction
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_f42_vector, reads_subnormal_inputs

    !> A few units in the last place of an O(1) value: the rounding of one `atan2`, or of `PI`
    real(real64), parameter :: TOL_ULP = 1.0e-15_real64
    real(real64), parameter :: EPS = epsilon(1.0_real64)
    !> The smallest normal number, and one below it that is subnormal
    real(real64), parameter :: SMALLEST_NORMAL = tiny(1.0_real64)
    real(real64), parameter :: SUBNORMAL = tiny(1.0_real64)/4
    !> The scaling boundaries and exponents scaling_exponent documents
    real(real64), parameter :: LARGEST_UNSCALED = 2.0_real64**496, SMALLEST_UNSCALED = 2.0_real64**(-470)
    integer(int32), parameter :: SCALE_DOWN = -600_int32, SCALE_UP = 600_int32

contains

    function get_all_tests_f42_vector() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(8))
        all_tests(1) = test_case("test_scaled_length_plain", test_scaled_length_plain)
        all_tests(2) = test_case("test_scaled_length_boundaries", test_scaled_length_boundaries)
        all_tests(3) = test_case("test_has_direction", test_has_direction)
        all_tests(4) = test_case("test_angle_to_direction_exact", test_angle_to_direction_exact)
        all_tests(5) = test_case("test_angle_to_direction_near_0_and_pi", test_angle_to_direction_near_0_and_pi)
        all_tests(6) = test_case("test_angle_to_direction_extreme_magnitudes", test_angle_to_direction_extreme_magnitudes)
        all_tests(7) = test_case("test_angle_to_direction_without_direction", test_angle_to_direction_without_direction)
        all_tests(8) = test_case("test_angle_to_direction_subnormal", test_angle_to_direction_subnormal)
    end function get_all_tests_f42_vector

    !> Whether arithmetic reads a subnormal operand as itself, decided at run time: false where the
    !| floating-point mode reads subnormal inputs as zero (denormals-are-zero, as ifx's flushing mode
    !| does; gfortran's flushes results only). A vector made only of subnormal numbers has a direction
    !| exactly when this holds, since scaled_length scales it into the normal range before any other
    !| arithmetic. The subnormal 2^-1024 is loaded from its bits so that nothing can be folded.
    logical function reads_subnormal_inputs()
        integer(int64), volatile :: bits
        real(real64), volatile :: subnormal

        bits = 1125899906842624_int64  ! 2^50: the bit pattern of 2^-1024
        subnormal = transfer(bits, 1.0_real64)
        reads_subnormal_inputs = scale(subnormal, 600) > 0.0_real64
    end function reads_subnormal_inputs

    !> `abs(actual - expected) <= tol`, which a NaN fails
    subroutine assert_close(actual, expected, tol, msg)
        real(real64), intent(in) :: actual, expected, tol
        character(*), intent(in) :: msg

        call assert_true(abs(actual - expected) <= tol, msg//": got "//actual//", expected "//expected)
    end subroutine assert_close

    !> `actual` equals `expected` to the last bit (+0 and -0 alike), which a NaN fails
    subroutine assert_exact(actual, expected, msg)
        real(real64), intent(in) :: actual, expected
        character(*), intent(in) :: msg

        call assert_true(abs(actual - expected) <= 0.0_real64, msg//": got "//actual//", expected "//expected)
    end subroutine assert_exact

    !> Call scaled_length and check both outputs exactly.
    subroutine check_scaled_length(vector, expected_exponent, expected_norm, msg)
        real(real64), contiguous, intent(in) :: vector(:)
        integer(int32), intent(in) :: expected_exponent
        real(real64), intent(in) :: expected_norm
        character(*), intent(in) :: msg
        integer(int32) :: exponent
        real(real64) :: scaled_norm

        call scaled_length(size(vector, kind=int32), vector, exponent, scaled_norm)
        call assert_equal_int(exponent, expected_exponent, msg//": exponent")
        call assert_exact(scaled_norm, expected_norm, msg//": scaled norm")
    end subroutine check_scaled_length

    !> Call scaled_length on a vector of length zero: its exponent means nothing.
    subroutine check_zero_length(vector, msg)
        real(real64), contiguous, intent(in) :: vector(:)
        character(*), intent(in) :: msg
        integer(int32) :: exponent
        real(real64) :: scaled_norm

        call scaled_length(size(vector, kind=int32), vector, exponent, scaled_norm)
        call assert_exact(scaled_norm, 0.0_real64, msg//": scaled norm")
    end subroutine check_zero_length

    !> Ordinary vectors need no scaling; the zero vector has length zero.
    subroutine test_scaled_length_plain()
        ! (3, 4): 9 + 16 = 25, sqrt 5, all exact
        call check_scaled_length([3.0_real64, 4.0_real64], 0_int32, 5.0_real64, "(3, 4)")
        call check_scaled_length([-3.0_real64, 0.0_real64, 4.0_real64], 0_int32, 5.0_real64, "(-3, 0, 4)")
        call check_zero_length([0.0_real64, 0.0_real64, 0.0_real64], "zero vector")
        call check_zero_length([-0.0_real64, 0.0_real64], "-0 vector")
    end subroutine test_scaled_length_plain

    !> The three boundaries: subnormal, 2^-470 and 2^496. On each side the exponent is the documented
    !| one, and since it scales by a power of two the scaled norm is exact: (3, 4) * 2^e has the scaled
    !| norm 5 * 2^(e + exponent).
    subroutine test_scaled_length_boundaries()
        integer(int32) :: exponent
        real(real64) :: scaled_norm
        logical :: gradual_underflow

        ! largest component subnormal: below 2^-470, so scaled up, to 2^-424 and -2^-425, whose
        ! squares are normal; the norm is 2^-424 * sqrt(1 + 1/4). Where subnormal inputs read as
        ! zero, it is the zero vector instead.
        call scaled_length(3, [SUBNORMAL, -SUBNORMAL/2, 0.0_real64], exponent, scaled_norm)
        call assert_equal_int(exponent, SCALE_UP, "all subnormal: exponent")
        if (reads_subnormal_inputs()) then
            call assert_close(scaled_norm/scale(SUBNORMAL, SCALE_UP), sqrt(1.25_real64), 2*EPS, "all subnormal: norm")
        else
            call assert_exact(scaled_norm, 0.0_real64, "all subnormal, read as zero: norm")
        end if
        ! largest component the smallest normal number: below 2^-470, so scaled up; the norm of a
        ! single component is that component, tiny * 2^600 = 2^-422
        call check_scaled_length([SMALLEST_NORMAL, 0.0_real64], SCALE_UP, scale(SMALLEST_NORMAL, SCALE_UP), "tiny")
        ! a subnormal next to a normal component takes part once scaled: tiny/4 and tiny become
        ! 2^-424 and 2^-422, whose squares are normal, so the norm is 2^-422 * sqrt(1 + 1/16) --
        ! where subnormals exist at all; an environment that flushes them reads tiny/4 as zero
        call scaled_length(2, [SUBNORMAL, SMALLEST_NORMAL], exponent, scaled_norm)
        call assert_equal_int(exponent, SCALE_UP, "tiny beside a subnormal: exponent")
        call ieee_get_underflow_mode(gradual_underflow)
        if (gradual_underflow) then
            call assert_close(scaled_norm/scale(SMALLEST_NORMAL, SCALE_UP), sqrt(17.0_real64)/4, 2*EPS, &
                              "tiny beside a subnormal: norm")
        end if
        ! exactly 2^-470 needs no scaling (2^-940 is still normal); one ulp below it does
        call check_scaled_length([SMALLEST_UNSCALED, 0.0_real64], 0_int32, SMALLEST_UNSCALED, "2^-470")
        call check_scaled_length([nearest(SMALLEST_UNSCALED, -1.0_real64), 0.0_real64], SCALE_UP, &
                                 scale(nearest(SMALLEST_UNSCALED, -1.0_real64), SCALE_UP), "2^-470 - ulp")
        ! (3, 4) * 2^-1000: largest 2^-998, scaled up to (3, 4) * 2^-400, norm 5 * 2^-400
        call check_scaled_length([3.0_real64, 4.0_real64]*2.0_real64**(-1000), SCALE_UP, 5*2.0_real64**(-400), &
                                 "(3, 4) * 2^-1000")
        ! (3, 4) * 2^494: largest exactly 2^496, no scaling; 25 * 2^988 is still finite
        call check_scaled_length([3.0_real64, 4.0_real64]*2.0_real64**494, 0_int32, 5*2.0_real64**494, &
                                 "(3, 4) * 2^494")
        ! one ulp above 2^496: scaled down
        call check_scaled_length([nearest(LARGEST_UNSCALED, 1.0_real64), 0.0_real64], SCALE_DOWN, &
                                 scale(nearest(LARGEST_UNSCALED, 1.0_real64), SCALE_DOWN), "2^496 + ulp")
        ! (3, 4) * 2^1000: scaled down to (3, 4) * 2^400, norm 5 * 2^400
        call check_scaled_length([3.0_real64, 4.0_real64]*2.0_real64**1000, SCALE_DOWN, 5*2.0_real64**400, &
                                 "(3, 4) * 2^1000")
        ! (huge, huge): the plain sum of squares overflows; scaled, the norm is huge * 2^-600 * sqrt(2)
        call scaled_length(2, [huge(1.0_real64), huge(1.0_real64)], exponent, scaled_norm)
        call assert_equal_int(exponent, SCALE_DOWN, "(huge, huge): exponent")
        call assert_true(ieee_is_finite(scaled_norm), "(huge, huge): finite scaled norm")
        call assert_close(scaled_norm/scale(huge(1.0_real64), SCALE_DOWN), sqrt(2.0_real64), 2*EPS, "(huge, huge): norm")
    end subroutine test_scaled_length_boundaries

    !> Every vector but the zero vector has a direction, one made only of subnormal numbers too,
    !| unless the build reads those as zero.
    subroutine test_has_direction()
        logical :: subnormals_read

        subnormals_read = reads_subnormal_inputs()
        call assert_false(has_direction(3, [0.0_real64, 0.0_real64, 0.0_real64]), "zero vector")
        call assert_false(has_direction(2, [-0.0_real64, 0.0_real64]), "-0 vector")
        call assert_true(has_direction(2, [SUBNORMAL, -SUBNORMAL]) .eqv. subnormals_read, "all subnormal")
        call assert_true(has_direction(1, [nearest(SMALLEST_NORMAL, -1.0_real64)]) .eqv. subnormals_read, &
                         "largest subnormal")
        call assert_true(has_direction(1, [SMALLEST_NORMAL]), "smallest normal")
        call assert_true(has_direction(3, [0.0_real64, SUBNORMAL, 1.0e-300_real64]), "subnormal beside 1e-300")
        call assert_true(has_direction(2, [huge(1.0_real64), -huge(1.0_real64)]), "(huge, -huge)")
        call assert_true(has_direction(2, [3.0_real64, 4.0_real64]), "(3, 4)")
    end subroutine test_has_direction

    !> Identical, power-of-two multiple and antipodal directions give 0 and pi to the last bit.
    subroutine test_angle_to_direction_exact()
        real(real64), parameter :: v(3) = [1.0_real64, 2.0_real64, 3.0_real64]

        ! u = v: |u - v| = 0, so 2 atan2(0, 2) = 0
        call assert_exact(angle_to_direction(3, v, v), 0.0_real64, "identical")
        call assert_exact(angle_to_direction(2, [3.0_real64, 4.0_real64], [3.0_real64, 4.0_real64]), 0.0_real64, &
                          "identical (3, 4)")
        ! power-of-two multiples have the same unit vector, unscaled (2^40), across the
        ! scale-down boundary (2^600) and across the scale-up boundary (2^-700)
        call assert_exact(angle_to_direction(3, v, v*2.0_real64**40), 0.0_real64, "v and v * 2^40")
        call assert_exact(angle_to_direction(3, v*2.0_real64**600, v), 0.0_real64, "v * 2^600 and v")
        call assert_exact(angle_to_direction(3, v, v*2.0_real64**(-700)), 0.0_real64, "v and v * 2^-700")
        ! -2v has the unit vector -u exactly: u + v = 0, so 2 atan2(2, 0) = 2 fl(pi/2) = PI
        call assert_exact(angle_to_direction(3, v, -2*v), PI, "v and -2v")
        ! orthogonal: |u - v| = |u + v| = sqrt(2), 2 atan2(1, 1) = pi/2
        call assert_close(angle_to_direction(2, [1.0_real64, 0.0_real64], [0.0_real64, 5.0_real64]), PI/2, TOL_ULP, &
                          "orthogonal")
        ! 60 degrees: (1, 0) and (1, sqrt 3)/2
        call assert_close(angle_to_direction(2, [1.0_real64, 0.0_real64], [0.5_real64, sqrt(3.0_real64)/2]), PI/3, &
                          2*TOL_ULP, "60 degrees")
    end subroutine test_angle_to_direction_exact

    !> Near 0 and pi the angle keeps its full relative precision, where the acos of a dot product
    !| is off by ~1e-8: cos(1e-10) rounds to 1, so acos gives 0 there, and pi instead of pi - 1e-10.
    !| (1, 1e-10) is its own unit vector, since 1 + 1e-20 rounds to 1. Against (1, 0):
    !| |u - v| = 1e-10, |u + v| = 2, and 2 atan2(1e-10, 2) = 1e-10 to a few ulp. Against (-1, 0):
    !| |u - v| = 2, |u + v| = 1e-10, and 2 atan2(2, 1e-10) = pi - 1e-10 to an ulp of pi.
    subroutine test_angle_to_direction_near_0_and_pi()
        real(real64), parameter :: d = 1.0e-10_real64

        call assert_close(angle_to_direction(2, [1.0_real64, d], [1.0_real64, 0.0_real64]), d, 4*EPS*d, "near 0")
        call assert_close(angle_to_direction(2, [1.0_real64, 0.0_real64], [1.0_real64, d]), d, 4*EPS*d, &
                          "near 0, swapped")
        call assert_close(angle_to_direction(2, [1.0_real64, d], [-1.0_real64, 0.0_real64]), PI - d, 2*spacing(PI), &
                          "near pi")
        ! the smallest angle a direction can have from another: (1, 2^-60) against (1, 0) is 2^-60
        call assert_close(angle_to_direction(2, [1.0_real64, 2.0_real64**(-60)], [1.0_real64, 0.0_real64]), &
                          2.0_real64**(-60), 4*EPS*2.0_real64**(-60), "2^-60")
    end subroutine test_angle_to_direction_near_0_and_pi

    !> Magnitudes whose squares overflow or underflow still give the right angle: 45 degrees
    !| between (a, a) and (b, 0) for any a, b > 0 of normal size.
    subroutine test_angle_to_direction_extreme_magnitudes()
        real(real64), parameter :: magnitudes(6) = [huge(1.0_real64), 1.0e300_real64, 1.0e155_real64, &
                                                    1.0e-155_real64, 1.0e-300_real64, tiny(1.0_real64)]
        integer(int32) :: i_a, i_b
        real(real64) :: a, b, vector(2), direction(2)

        do i_a = 1, size(magnitudes)
            do i_b = 1, size(magnitudes)
                a = magnitudes(i_a)
                b = magnitudes(i_b)
                vector = a
                direction(1) = b
                direction(2) = 0.0_real64
                call assert_close(angle_to_direction(2, vector, direction), PI/4, TOL_ULP, &
                                  "45 degrees, a = "//a//", b = "//b)
            end do
        end do
    end subroutine test_angle_to_direction_extreme_magnitudes

    !> -1 where either vector has no direction: it is zero.
    subroutine test_angle_to_direction_without_direction()
        real(real64) :: no_axes(0)

        call assert_exact(angle_to_direction(2, [0.0_real64, 0.0_real64], [1.0_real64, 0.0_real64]), -1.0_real64, &
                          "zero vector")
        call assert_exact(angle_to_direction(2, [1.0_real64, 0.0_real64], [0.0_real64, 0.0_real64]), -1.0_real64, &
                          "zero direction")
        ! a zero-length vector has no direction either
        call assert_exact(angle_to_direction(0, no_axes, no_axes), -1.0_real64, "no axes")
    end subroutine test_angle_to_direction_without_direction

    !> Vectors made only of subnormal numbers have their directions: (s, s) scales to
    !| (2^-424, 2^-424), whose unit vector is that of (1, 1) to the last bit, so the angle is 0;
    !| (-s, 0) against (1, 1) is 3 pi/4. Where subnormal inputs read as zero, both are -1.
    subroutine test_angle_to_direction_subnormal()
        if (reads_subnormal_inputs()) then
            call assert_exact(angle_to_direction(2, [SUBNORMAL, SUBNORMAL], [1.0_real64, 1.0_real64]), 0.0_real64, &
                              "subnormal vector")
            call assert_close(angle_to_direction(2, [1.0_real64, 1.0_real64], [-SUBNORMAL, 0.0_real64]), 3*PI/4, &
                              2*TOL_ULP, "subnormal direction")
        else
            call assert_exact(angle_to_direction(2, [SUBNORMAL, SUBNORMAL], [1.0_real64, 1.0_real64]), -1.0_real64, &
                              "subnormal vector, read as zero")
            call assert_exact(angle_to_direction(2, [1.0_real64, 1.0_real64], [-SUBNORMAL, 0.0_real64]), -1.0_real64, &
                              "subnormal direction, read as zero")
        end if
    end subroutine test_angle_to_direction_subnormal
end module mod_test_f42_vector
