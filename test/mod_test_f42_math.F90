!> Unit test suite for the angle and logarithm helpers of f42_math_impl: `angular_distance`,
!| `one_minus_cosine` and `log_one_minus`.
!|
!| Expected values are derived by hand in the comments. Where a value is not exact, the reference
!| is the true value rounded to double precision -- for `log_one_minus` computed to 60 digits
!| (Python `decimal`, `(1 - Decimal(x)).ln()`, `x` taken exactly) -- and the tolerance is stated
!| in units in the last place. Every comparison is written so that a NaN fails it.
module mod_test_f42_math
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use f42_math_impl, only: PI, angular_distance, one_minus_cosine, log_one_minus
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_f42_math

    !> A few units in the last place of an O(1) value: the rounding of one transcendental
    !| function call, or of `PI` itself
    real(real64), parameter :: TOL_ULP = 1.0e-15_real64
    real(real64), parameter :: EPS = epsilon(1.0_real64)

contains

    function get_all_tests_f42_math() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(6))
        all_tests(1) = test_case("test_angular_distance", test_angular_distance)
        all_tests(2) = test_case("test_one_minus_cosine", test_one_minus_cosine)
        all_tests(3) = test_case("test_log_one_minus_small", test_log_one_minus_small)
        all_tests(4) = test_case("test_log_one_minus_series_terms", test_log_one_minus_series_terms)
        all_tests(5) = test_case("test_log_one_minus_at_the_switch", test_log_one_minus_at_the_switch)
        all_tests(6) = test_case("test_log_one_minus_large", test_log_one_minus_large)
    end function get_all_tests_f42_math

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

    !> The distance goes around the circle: never more than pi, the same both ways.
    subroutine test_angular_distance()
        real(real64) :: near_pi

        ! equal angles: 0 exactly
        call assert_exact(angular_distance(1.0_real64, 1.0_real64), 0.0_real64, "equal angles")
        ! pi and 0 are pi apart both ways: pi - 0 = pi stays pi; 0 - pi = -pi wraps to pi
        call assert_exact(angular_distance(PI, 0.0_real64), PI, "pi from 0")
        call assert_exact(angular_distance(0.0_real64, PI), PI, "0 from pi")
        ! across the seam: a = fl(pi - 0.1) and -a are 2a apart one way, 2 pi - 2a = 0.2 the other
        near_pi = PI - 0.1_real64
        call assert_close(angular_distance(near_pi, -near_pi), 0.2_real64, TOL_ULP, "across +pi")
        call assert_close(angular_distance(-near_pi, near_pi), 0.2_real64, TOL_ULP, "across -pi")
        ! pi and -pi/2: 3 pi/2 one way, pi/2 the other
        call assert_close(angular_distance(PI, -PI/2), PI/2, TOL_ULP, "pi from -pi/2")
        ! 3 and -3: 6 one way, 2 pi - 6 = 0.28318530717958623 the other
        call assert_close(angular_distance(3.0_real64, -3.0_real64), 0.28318530717958623_real64, TOL_ULP, &
                          "3 from -3")
        ! inside the half circle nothing wraps: 0.3 - 0.1 rounds to 0.19999999999999998
        call assert_exact(angular_distance(0.3_real64, 0.1_real64), 0.3_real64 - 0.1_real64, "0.3 from 0.1")
        call assert_exact(angular_distance(0.1_real64, 0.3_real64), 0.3_real64 - 0.1_real64, "0.1 from 0.3")
    end subroutine test_angular_distance

    !> 2 sin^2(x/2) keeps the small values that 1 - cos x rounds to zero.
    subroutine test_one_minus_cosine()
        ! 0 exactly
        call assert_exact(one_minus_cosine(0.0_real64), 0.0_real64, "x = 0")
        ! 1 - cos(1e-10) = 5e-21 (1 - 1e-20/12 + ...): cos(1e-10) rounds to 1, so the naive form
        ! gives 0. sin(5e-11) = 5e-11 to double precision; squaring and doubling add two roundings
        call assert_close(one_minus_cosine(1.0e-10_real64), 5.0e-21_real64, 4*EPS*5.0e-21_real64, "x = 1e-10")
        call assert_close(one_minus_cosine(-1.0e-10_real64), 5.0e-21_real64, 4*EPS*5.0e-21_real64, "x = -1e-10")
        ! 1 - cos(1e-4) = x^2/2 - x^4/24 + x^6/720 = 5e-9 - 4.1666...e-18 + 1.4e-27
        call assert_close(one_minus_cosine(1.0e-4_real64), 4.9999999958333333e-9_real64, 4*EPS*5.0e-9_real64, &
                          "x = 1e-4")
        ! large angles: 1 - cos(pi/2) = 1, 1 - cos(2 pi/3) = 3/2, 1 - cos(pi) = 2, 1 - cos(3 pi/2) = 1
        call assert_close(one_minus_cosine(PI/2), 1.0_real64, TOL_ULP, "x = pi/2")
        call assert_close(one_minus_cosine(2*PI/3), 1.5_real64, 2*TOL_ULP, "x = 2 pi/3")
        call assert_close(one_minus_cosine(PI), 2.0_real64, 2*TOL_ULP, "x = pi")
        call assert_close(one_minus_cosine(3*PI/2), 1.0_real64, 2*TOL_ULP, "x = 3 pi/2")
    end subroutine test_one_minus_cosine

    !> Small x: the series keeps every digit that 1 - x would round away.
    subroutine test_log_one_minus_small()
        real(real64) :: x

        ! 0 exactly, and a tiny x gives -x exactly: -x*(1 + x*(...)) with 1 + x*(...) rounding to 1
        call assert_exact(log_one_minus(0.0_real64), 0.0_real64, "x = 0")
        call assert_exact(log_one_minus(1.0e-300_real64), -1.0e-300_real64, "x = 1e-300")
        ! ln(1 - 1e-8) = -1.00000000500000003e-8; log(1 - 1e-8) could be off by up to 1e-8 relative
        x = 1.0e-8_real64
        call assert_close(log_one_minus(x), -1.0000000050000001e-8_real64, 2*spacing(1.0e-8_real64), "x = 1e-8")
        ! x = 2^-30: ln(1 - x) = -9.313225750491594e-10 (reference to 60 digits)
        x = 2.0_real64**(-30)
        call assert_close(log_one_minus(x), -9.313225750491594e-10_real64, 2*spacing(9.3e-10_real64), "x = 2^-30")
    end subroutine test_log_one_minus_small

    !> Every coefficient of the series matters, even the last. The x^6/6 term is
    !| x^6/6 = 1.2e-19 at x = 9.5e-4, about one ulp of the result; replacing x/6 by x/5 adds
    !| x^6/30, about a quarter ulp, so only a result exactly at the correctly rounded value can
    !| see it. At the three x below -- k * 2^-34, exact -- the series gives the correctly rounded
    !| ln(1 - x) whether or not the compiler fuses its multiply-adds, while x/5 in its place
    !| rounds one ulp away (searched with the Horner scheme evaluated plainly and fused).
    !| A wrong x^2 coefficient (0.5000001) is off by 1e-7 x^2, a million ulp.
    subroutine test_log_one_minus_series_terms()
        integer(int32), parameter :: n_points = 3
        real(real64), parameter :: k(n_points) = [16300343.0_real64, 16301379.0_real64, 16302485.0_real64]
        ! ln(1 - k * 2^-34), correctly rounded
        real(real64), parameter :: expected(n_points) = [-0.0009492552383447765_real64, &
                                                         -0.0009493155987535507_real64, &
                                                         -0.0009493800375723441_real64]
        integer(int32) :: i_point

        do i_point = 1, n_points
            call assert_exact(log_one_minus(k(i_point)*2.0_real64**(-34)), expected(i_point), &
                              "series correctly rounded at k = "//k(i_point))
        end do
    end subroutine test_log_one_minus_series_terms

    !> Both sides of the switch at x = 1e-3. Below, the series: within 2 ulp. From 1e-3 on,
    !| log(1 - x): 1 - x rounds by at most 2^-54 = 5.6e-17, which moves ln(1 - x) by as much
    !| again, relative 5.6e-14 at x = 1e-3; measured it is 3.9 ulp there.
    subroutine test_log_one_minus_at_the_switch()
        real(real64) :: x_below

        ! ln(1 - 1e-3) = -0.0010005003335835335; one ulp lower, -0.0010005003335835333
        x_below = nearest(1.0e-3_real64, -1.0_real64)
        call assert_close(log_one_minus(x_below), -0.0010005003335835333_real64, 2*spacing(1.0e-3_real64), &
                          "just below 1e-3")
        call assert_close(log_one_minus(1.0e-3_real64), -0.0010005003335835335_real64, 6.0e-17_real64, "x = 1e-3")
    end subroutine test_log_one_minus_at_the_switch

    !> Large x: log(1 - x), exact 1 - x from 1/2 on.
    subroutine test_log_one_minus_large()
        ! ln(1/2) = -ln 2 = -0.6931471805599453; ln(3/4) = -0.2876820724517809
        call assert_close(log_one_minus(0.5_real64), -0.6931471805599453_real64, TOL_ULP, "x = 1/2")
        call assert_close(log_one_minus(0.25_real64), -0.2876820724517809_real64, TOL_ULP, "x = 1/4")
        ! the largest x below 1: 1 - x = 2^-53 exactly, ln(2^-53) = -53 ln 2 = -36.7368005696771
        call assert_close(log_one_minus(nearest(1.0_real64, -1.0_real64)), -36.7368005696771_real64, 1.0e-13_real64, &
                          "x = 1 - 2^-53")
    end subroutine test_log_one_minus_large
end module mod_test_f42_math
