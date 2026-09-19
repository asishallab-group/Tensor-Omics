!> The `normalize_variable_timeseries` cases: min-max scaling of one series to [0, 1], the
!| constant series, extreme magnitudes, the dimension check and the NaN/Inf check.
!|
!| Every series has a power-of-two range (max - min), so each (v - min)/range is exact, even
!| where an optimizer multiplies by the reciprocal of the range: no case needs a tolerance.
module mod_test_normalize_variable_timeseries
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_trajectory_normalization, only: normalize_variable_timeseries
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_normalize_variable_timeseries

    !> What `v_norm` holds before a call, so that an output the procedure never writes shows.
    real(real64), parameter :: UNWRITTEN = -1.0_real64

contains

    !> Every case of `normalize_variable_timeseries`.
    function get_all_tests_normalize_variable_timeseries() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(9))
        all_tests(1) = test_case("test_normalize_variable_timeseries_ascending", &
                                 test_normalize_variable_timeseries_ascending)
        all_tests(2) = test_case("test_normalize_variable_timeseries_unordered_and_negative", &
                                 test_normalize_variable_timeseries_unordered_and_negative)
        all_tests(3) = test_case("test_normalize_variable_timeseries_large_offset", &
                                 test_normalize_variable_timeseries_large_offset)
        all_tests(4) = test_case("test_normalize_variable_timeseries_constant_series", &
                                 test_normalize_variable_timeseries_constant_series)
        all_tests(5) = test_case("test_normalize_variable_timeseries_tiny_range", &
                                 test_normalize_variable_timeseries_tiny_range)
        all_tests(6) = test_case("test_normalize_variable_timeseries_largest_range", &
                                 test_normalize_variable_timeseries_largest_range)
        all_tests(7) = test_case("test_normalize_variable_timeseries_range_overflows", &
                                 test_normalize_variable_timeseries_range_overflows)
        all_tests(8) = test_case("test_normalize_variable_timeseries_dimensions", &
                                 test_normalize_variable_timeseries_dimensions)
        all_tests(9) = test_case("test_normalize_variable_timeseries_rejects_nan_and_inf", &
                                 test_normalize_variable_timeseries_rejects_nan_and_inf)
    end function get_all_tests_normalize_variable_timeseries

    !> [1, 2, 3, 4, 5]: min 1, range 4, so (v - 1)/4 = [0, 0.25, 0.5, 0.75, 1]. And that result
    !| is a fixed point: a series already spanning exactly [0, 1] has min 0 and range 1, so it
    !| comes back unchanged.
    subroutine test_normalize_variable_timeseries_ascending()
        integer(int32), parameter :: n_points = 5
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: status, ierr

        v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        v_norm = UNWRITTEN

        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_ascending: ierr")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_ascending: status")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_ascending: v_norm")

        v = expected
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_ascending: ierr, already in [0, 1]")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_ascending: status, already in [0, 1]")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_ascending: v_norm, already in [0, 1]")
    end subroutine test_normalize_variable_timeseries_ascending

    !> Minimum and maximum inside the series, not at its ends, and negative values:
    !| [0, -8, 4, 8, -4] has min -8 (point 2), max 8 (point 4), range 16, so (v + 8)/16 =
    !| [8, 0, 12, 16, 4]/16 = [0.5, 0, 0.75, 1, 0.25].
    subroutine test_normalize_variable_timeseries_unordered_and_negative()
        integer(int32), parameter :: n_points = 5
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: status, ierr

        v = [0.0_real64, -8.0_real64, 4.0_real64, 8.0_real64, -4.0_real64]
        expected = [0.5_real64, 0.0_real64, 0.75_real64, 1.0_real64, 0.25_real64]
        v_norm = UNWRITTEN

        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_unordered_and_negative: ierr")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_unordered_and_negative: status")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_unordered_and_negative: v_norm")
    end subroutine test_normalize_variable_timeseries_unordered_and_negative

    !> A range of 1 on top of 2**50 (about 1.1e15, where the spacing of doubles is 0.25): the
    !| series is tiny relative to its values but not constant, so it is scaled like any other.
    !| v - 2**50 = [0, 0.25, 0.5, 1] exactly, and the range is 1.
    subroutine test_normalize_variable_timeseries_large_offset()
        integer(int32), parameter :: n_points = 4
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points), offset
        integer(int32) :: status, ierr

        offset = scale(1.0_real64, 50)
        v = [0.0_real64, 0.25_real64, 0.5_real64, 1.0_real64]
        v = offset + v
        expected = [0.0_real64, 0.25_real64, 0.5_real64, 1.0_real64]
        v_norm = UNWRITTEN

        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_large_offset: ierr")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_large_offset: status")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_large_offset: v_norm")
    end subroutine test_normalize_variable_timeseries_large_offset

    !> A constant series has no range to divide by.
    !| Pinned for a non-zero constant, for zeros, and for the single point (n_points = 1, the
    !| smallest valid size), which is always constant.
    subroutine test_normalize_variable_timeseries_constant_series()
        integer(int32), parameter :: n_points = 4
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: status, ierr

        ! QUESTION: the contract does not say what a constant series gives. Today it is all zeros,
        ! with `status` = ERR_DIVISION_BY_ZERO and `ierr` = ERR_OK: a per-series warning, not an
        ! error, and one that reaches neither `ierr` nor any caller that ignores `status`.
        expected = 0.0_real64

        v = 3.14_real64
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_constant_series: ierr, 3.14")
        call assert_err(status, ERR_DIVISION_BY_ZERO, "test_normalize_variable_timeseries_constant_series: status, 3.14")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_constant_series: v_norm, 3.14")

        v = 0.0_real64
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_constant_series: ierr, zeros")
        call assert_err(status, ERR_DIVISION_BY_ZERO, "test_normalize_variable_timeseries_constant_series: status, zeros")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_constant_series: v_norm, zeros")

        v = 42.0_real64
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, 1, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_constant_series: ierr, one point")
        call assert_err(status, ERR_DIVISION_BY_ZERO, &
                        "test_normalize_variable_timeseries_constant_series: status, one point")
        call assert_equal_real(v_norm(1), 0.0_real64, 0.0_real64, &
                               "test_normalize_variable_timeseries_constant_series: v_norm, one point")
        call assert_equal_real(v_norm(2), UNWRITTEN, 0.0_real64, &
                               "test_normalize_variable_timeseries_constant_series: wrote past n_points = 1")
    end subroutine test_normalize_variable_timeseries_constant_series

    !> A series with a tiny range is still a series, not a constant: [0, 2**-41, 2**-40] scales to
    !| [0, 0.5, 1], and [1, 2, 3]*tiny (range 2*tiny, differences still normal numbers) to
    !| [0, 0.5, 1] as well. Both are exact.
    subroutine test_normalize_variable_timeseries_tiny_range()
        integer(int32), parameter :: n_points = 3
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: status, ierr

        ! BUG: the constant test is `is_close(max - min, 0)`, whose absolute floor of 1e-12 calls
        ! every range below 1e-12 constant, whatever the magnitude of the data: both series come
        ! back as zeros with status ERR_DIVISION_BY_ZERO. (The old suite pinned those zeros for
        ! [1, 2, 3]*tiny on purpose; the scale-free answer is [0, 0.5, 1], as for [1, 2, 3] itself.)
        expected = [0.0_real64, 0.5_real64, 1.0_real64]

        v = [0.0_real64, scale(1.0_real64, -41), scale(1.0_real64, -40)]
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_tiny_range: ierr, range 2**-40")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_tiny_range: status, range 2**-40")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_tiny_range: v_norm, range 2**-40")

        v = [tiny(1.0_real64), 2.0_real64*tiny(1.0_real64), 3.0_real64*tiny(1.0_real64)]
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_tiny_range: ierr, range 2*tiny")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_tiny_range: status, range 2*tiny")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_tiny_range: v_norm, range 2*tiny")
    end subroutine test_normalize_variable_timeseries_tiny_range

    !> The largest ranges that still fit a double, just inside the overflow of max - min:
    !| [-2**1022, 0, 2**1022, 2**1021] has range 2**1023, so [0, 2**1022, 2**1023, 2**1022 + 2**1021]
    !| /2**1023 = [0, 0.5, 1, 0.75]; and [0, huge] has range huge, so [0, 1]. huge/huge is exact
    !| because it is a division; a reciprocal 1/huge would be subnormal, which the ifx builds'
    !| -fp-model precise rules out.
    subroutine test_normalize_variable_timeseries_largest_range()
        real(real64) :: v(4), v_norm(4), expected(4)
        integer(int32) :: status, ierr

        v = [-scale(1.0_real64, 1022), 0.0_real64, scale(1.0_real64, 1022), scale(1.0_real64, 1021)]
        expected = [0.0_real64, 0.5_real64, 1.0_real64, 0.75_real64]
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, 4, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_largest_range: ierr, +-2**1022")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_largest_range: status, +-2**1022")
        call assert_equal_array_real(v_norm, expected, 4, 0.0_real64, &
                                     "test_normalize_variable_timeseries_largest_range: v_norm, +-2**1022")

        v(1:2) = [0.0_real64, huge(1.0_real64)]
        expected(1:2) = [0.0_real64, 1.0_real64]
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, 2, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_largest_range: ierr, [0, huge]")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_largest_range: status, [0, huge]")
        call assert_equal_array_real(v_norm, expected, 2, 0.0_real64, &
                                     "test_normalize_variable_timeseries_largest_range: v_norm, [0, huge]")
    end subroutine test_normalize_variable_timeseries_largest_range

    !> Just outside: max - min overflows. The wrapper accepts every finite value, so these are
    !| valid input, and the scaled series is still [0, 0.5, 1]: [-2**1023, 0, 2**1023] has range
    !| 2**1024 (one past huge), [-huge, 0, huge] range 2*huge; 0 lies halfway in both.
    subroutine test_normalize_variable_timeseries_range_overflows()
        integer(int32), parameter :: n_points = 3
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: status, ierr

        ! Regression: the range used to overflow to Inf and v_norm came out as [0, 0, NaN] -- the
        ! middle point divided a finite value by Inf, the maximum Inf by Inf. TOX wrote a NaN, and
        ! status said OK. (assert_equal_array_real lets a NaN through, as NaN > tol is false: hence the no-NaN
        ! check.)
        expected = [0.0_real64, 0.5_real64, 1.0_real64]

        v = [-scale(1.0_real64, 1023), 0.0_real64, scale(1.0_real64, 1023)]
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_range_overflows: ierr, +-2**1023")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_range_overflows: status, +-2**1023")
        call assert_no_nan_real(v_norm, n_points, "test_normalize_variable_timeseries_range_overflows: NaN, +-2**1023")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_range_overflows: v_norm, +-2**1023")

        v = [-huge(1.0_real64), 0.0_real64, huge(1.0_real64)]
        v_norm = UNWRITTEN
        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_variable_timeseries_range_overflows: ierr, +-huge")
        call assert_err(status, ERR_OK, "test_normalize_variable_timeseries_range_overflows: status, +-huge")
        call assert_no_nan_real(v_norm, n_points, "test_normalize_variable_timeseries_range_overflows: NaN, +-huge")
        call assert_equal_array_real(v_norm, expected, n_points, 0.0_real64, &
                                     "test_normalize_variable_timeseries_range_overflows: v_norm, +-huge")
    end subroutine test_normalize_variable_timeseries_range_overflows

    !> n_points (argument 3) = 0 is ERR_EMPTY_INPUT and -1 ERR_INVALID_INPUT; 1, the smallest
    !| valid size, is in `test_normalize_variable_timeseries_constant_series`.
    subroutine test_normalize_variable_timeseries_dimensions()
        real(real64) :: v(1), v_norm(1)
        integer(int32) :: status, ierr

        v = 1.0_real64

        call normalize_variable_timeseries(v, v_norm, 0, status, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_normalize_variable_timeseries_dimensions: n_points = 0", 3)
        call normalize_variable_timeseries(v, v_norm, -1, status, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_normalize_variable_timeseries_dimensions: n_points = -1", 3)
    end subroutine test_normalize_variable_timeseries_dimensions

    !> NaN, +Inf and -Inf in `v` (argument 1) are ERR_NAN_INF. Each sits at the last point, so the
    !| check must cover the whole series.
    subroutine test_normalize_variable_timeseries_rejects_nan_and_inf()
        integer(int32), parameter :: n_points = 3
        character(*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        real(real64) :: v(n_points), v_norm(n_points), bad(3)
        integer(int32) :: status, ierr, i_bad

        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        v = [1.0_real64, 2.0_real64, 3.0_real64]
        do i_bad = 1, size(bad)
            v(n_points) = bad(i_bad)
            call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_normalize_variable_timeseries_rejects_nan_and_inf: "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_normalize_variable_timeseries_rejects_nan_and_inf

end module mod_test_normalize_variable_timeseries
