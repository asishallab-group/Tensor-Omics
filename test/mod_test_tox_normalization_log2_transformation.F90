!> The `log2_transformation` cases: log2(x + 1) where it is exact, the edge of its domain at -1,
!| the extremes of real64, and the input checks.
module mod_test_tox_normalization_log2_transformation
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tox_normalization
    use tox_errors
    use test_suite, only: test_case
    implicit none
    public

    ! log2 is computed as log(x + 1)/log(2), which is exact only to a few ulps
    real(real64), parameter :: TOL = 1.0e-12_real64

contains

    !> Get array of all available tests.
    function get_all_tests_tox_normalization_log2_transformation() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(8))
        all_tests(1) = test_case("test_log2_values", test_log2_values)
        all_tests(2) = test_case("test_log2_negative_values", test_log2_negative_values)
        all_tests(3) = test_case("test_log2_extremes", test_log2_extremes)
        all_tests(4) = test_case("test_log2_rejects_minus_one_and_below", test_log2_rejects_minus_one_and_below)
        all_tests(5) = test_case("test_log2_rejects_nan_and_inf", test_log2_rejects_nan_and_inf)
        all_tests(6) = test_case("test_log2_dimensions", test_log2_dimensions)
        all_tests(7) = test_case("test_log2_tiny_values", test_log2_tiny_values)
        all_tests(8) = test_case("test_log1p", test_log1p)
    end function get_all_tests_tox_normalization_log2_transformation

    !> f42_math's log1p, which log2_transformation is built on; it moves to the f42 suites with the
    !| f42 batch. log1p(0) = 0, and below half an ulp of 1 it is x itself; log1p(1) = log(2) and
    !| log1p(-0.5) = -log(2); for a tiny x it matches the series x - x**2/2 to a few ulps.
    subroutine test_log1p()
        use f42_math_impl, only: log1p
        real(real64), parameter :: x = 1.0e-12_real64

        call assert_equal_real(log1p(0.0_real64), 0.0_real64, 0.0_real64, "test_log1p: log1p(0)")
        call assert_equal_real(log1p(1.0e-300_real64), 1.0e-300_real64, 0.0_real64, "test_log1p: log1p(1e-300) is 1e-300")
        call assert_equal_real(log1p(1.0_real64), log(2.0_real64), 2*epsilon(1.0_real64), "test_log1p: log1p(1) = log(2)")
        call assert_equal_real(log1p(-0.5_real64), -log(2.0_real64), 2*epsilon(1.0_real64), "test_log1p: log1p(-0.5) = -log(2)")
        call assert_equal_real(log1p(x), x - x**2/2.0_real64, 4*epsilon(1.0_real64)*x, "test_log1p: log1p(1e-12)")
    end subroutine test_log1p

    !> For tiny x, forming x + 1 rounds away most of x's digits: 1 + 1e-15 is 1.00000000000000111,
    !| so log(x + 1) came out 11% high. The series log(1 + x) = x - x**2/2 + x**3/3 - ... is exact
    !| to double precision after two terms for every x here, so it gives the expected values.
    subroutine test_log2_tiny_values()
        real(real64), dimension(4, 1) :: expr, transformed, expected
        integer(int32) :: ierr, i_value

        expr(:, 1) = [1.0e-15_real64, -1.0e-15_real64, 1.0e-10_real64, 1.0e-8_real64]
        do i_value = 1, 4
            expected(i_value, 1) = (expr(i_value, 1) - expr(i_value, 1)**2/2.0_real64)/log(2.0_real64)
        end do

        call log2_transformation(1, 4, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_log2_tiny_values: ierr")
        do i_value = 1, 4
            call assert_equal_real(transformed(i_value, 1), expected(i_value, 1), 4*epsilon(1.0_real64)*abs(expected(i_value, 1)), &
                                   "test_log2_tiny_values: log2(1 + x) for tiny x, to a few ulps")
        end do
    end subroutine test_log2_tiny_values

    !> Where x + 1 is a power of two the result is its exponent: 0, 1, 3, 7, 15 and 1023 become
    !| 0, 1, 2, 3, 4 and 10, over a 2 x 3 matrix, so that both dimensions are walked.
    subroutine test_log2_values()
        integer(int32), parameter :: n_genes = 3, n_tissues = 2
        real(real64), dimension(n_tissues, n_genes) :: expr, transformed, expected
        integer(int32) :: ierr

        expr(:, 1) = [0.0_real64, 1.0_real64]
        expr(:, 2) = [3.0_real64, 7.0_real64]
        expr(:, 3) = [15.0_real64, 1023.0_real64]
        expected(:, 1) = [0.0_real64, 1.0_real64]
        expected(:, 2) = [2.0_real64, 3.0_real64]
        expected(:, 3) = [4.0_real64, 10.0_real64]

        call log2_transformation(n_genes, n_tissues, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_log2_values: ierr")
        call assert_equal_array_real(transformed, expected, n_genes*n_tissues, TOL, "test_log2_values: log2(x + 1)")
    end subroutine test_log2_values

    !> The domain is x > -1. Values in (-1, 0) map below zero: -0.5, -0.75 and -0.875 become -1, -2
    !| and -3. Just above -1 is still valid: -1 + 2**-52 gives x + 1 = 2**-52, so -52.
    subroutine test_log2_negative_values()
        real(real64), dimension(4, 1) :: expr, transformed, expected
        integer(int32) :: ierr

        expr(:, 1) = [-0.5_real64, -0.75_real64, -0.875_real64, -1.0_real64 + 2.0_real64**(-52)]
        expected(:, 1) = [-1.0_real64, -2.0_real64, -3.0_real64, -52.0_real64]

        call log2_transformation(1, 4, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_log2_negative_values: ierr")
        call assert_equal_array_real(transformed, expected, 4, TOL, "test_log2_negative_values: log2(x + 1)")
    end subroutine test_log2_negative_values

    !> The extremes of real64 stay finite: huge + 1 is huge, whose log2 is just under 1024, and tiny
    !| vanishes next to 1, so log2(tiny + 1) is exactly 0.
    subroutine test_log2_extremes()
        real(real64), dimension(2, 1) :: expr, transformed, expected
        integer(int32) :: ierr

        expr(:, 1) = [huge(1.0_real64), tiny(1.0_real64)]
        expected(:, 1) = [1024.0_real64, 0.0_real64]

        call log2_transformation(1, 2, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_log2_extremes: ierr")
        call assert_equal_array_real(transformed, expected, 2, TOL, "test_log2_extremes: log2(huge + 1) and log2(tiny + 1)")
    end subroutine test_log2_extremes

    !> x = -1 would take log2(0), and x < -1 the log2 of a negative number: both ERR_INVALID_INPUT.
    subroutine test_log2_rejects_minus_one_and_below()
        real(real64), dimension(2, 2) :: expr, transformed
        integer(int32) :: ierr

        expr(:, 1) = [0.0_real64, 1.0_real64]
        expr(:, 2) = [-1.0_real64, 3.0_real64]
        call log2_transformation(2, 2, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_log2_rejects_minus_one_and_below: x = -1")

        expr(:, 2) = [3.0_real64, -2.5_real64]
        call log2_transformation(2, 2, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_log2_rejects_minus_one_and_below: x < -1")
    end subroutine test_log2_rejects_minus_one_and_below

    !> NaN and Inf are rejected up front: TOX writes no NaN.
    subroutine test_log2_rejects_nan_and_inf()
        real(real64), dimension(2, 2) :: expr, transformed
        real(real64) :: bad(2)
        integer(int32) :: ierr, i_bad

        expr = 1.0_real64
        bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]
        do i_bad = 1, size(bad)
            expr(2, 2) = bad(i_bad)
            call log2_transformation(2, 2, expr, transformed, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                                  "test_log2_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
        end do
    end subroutine test_log2_rejects_nan_and_inf

    !> Each dimension on its own: zero is ERR_EMPTY_INPUT, negative ERR_INVALID_INPUT.
    subroutine test_log2_dimensions()
        real(real64), dimension(1, 1) :: expr, transformed
        integer(int32) :: ierr

        expr = 1.0_real64
        call log2_transformation(0, 1, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_log2_dimensions: n_genes = 0")
        call log2_transformation(1, 0, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_log2_dimensions: n_tissues = 0")
        call log2_transformation(-1, 1, expr, transformed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_log2_dimensions: n_genes = -1")
    end subroutine test_log2_dimensions

end module mod_test_tox_normalization_log2_transformation
