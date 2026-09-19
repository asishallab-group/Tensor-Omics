!> The `normalize_single_trajectory` cases: each factor (column) of one sample's
!| `trajectory(n_timepoints, n_factors)` scaled to [0, 1] over time on its own, a constant factor
!| among varying ones, extreme magnitudes, each dimension, and the NaN/Inf check.
!|
!| Every series but those of `test_normalize_single_trajectory_range_of_three` has a power-of-two
!| range (max - min), so each (v - min)/range is exact, even where an optimizer multiplies by the
!| reciprocal of the range, and needs no tolerance.
module mod_test_normalize_single_trajectory
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_trajectory_normalization, only: normalize_single_trajectory
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_normalize_single_trajectory

    !> What `trajectory_norm` holds before a call, so that an output the procedure never writes shows.
    real(real64), parameter :: UNWRITTEN = -1.0_real64

contains

    !> Every case of `normalize_single_trajectory`.
    function get_all_tests_normalize_single_trajectory() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(7))
        all_tests(1) = test_case("test_normalize_single_trajectory_values", &
                                 test_normalize_single_trajectory_values)
        all_tests(2) = test_case("test_normalize_single_trajectory_range_of_three", &
                                 test_normalize_single_trajectory_range_of_three)
        all_tests(3) = test_case("test_normalize_single_trajectory_constant_factor", &
                                 test_normalize_single_trajectory_constant_factor)
        all_tests(4) = test_case("test_normalize_single_trajectory_single_timepoint", &
                                 test_normalize_single_trajectory_single_timepoint)
        all_tests(5) = test_case("test_normalize_single_trajectory_extreme_magnitudes", &
                                 test_normalize_single_trajectory_extreme_magnitudes)
        all_tests(6) = test_case("test_normalize_single_trajectory_dimensions", &
                                 test_normalize_single_trajectory_dimensions)
        all_tests(7) = test_case("test_normalize_single_trajectory_rejects_nan_and_inf", &
                                 test_normalize_single_trajectory_rejects_nan_and_inf)
    end function get_all_tests_normalize_single_trajectory

    !> Two factors of four time points, each scaled over its own column:
    !| factor 1 [10, 11, 12, 14]: min 10, range 4 -> [0, 1, 2, 4]/4 = [0, 0.25, 0.5, 1];
    !| factor 2 [36, 28, 24, 20]: min 20, range 16 -> [16, 8, 4, 0]/16 = [1, 0.5, 0.25, 0].
    !| Scaling across factors instead of over time would give neither.
    subroutine test_normalize_single_trajectory_values()
        integer(int32), parameter :: n_factors = 2, n_timepoints = 4
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints, n_factors)
        integer(int32) :: status(n_factors), expected_status(n_factors), ierr

        trajectory(:, 1) = [10.0_real64, 11.0_real64, 12.0_real64, 14.0_real64]
        trajectory(:, 2) = [36.0_real64, 28.0_real64, 24.0_real64, 20.0_real64]
        expected(:, 1) = [0.0_real64, 0.25_real64, 0.5_real64, 1.0_real64]
        expected(:, 2) = [1.0_real64, 0.5_real64, 0.25_real64, 0.0_real64]
        expected_status = ERR_OK
        trajectory_norm = UNWRITTEN

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_single_trajectory_values: ierr")
        call assert_equal_array_int(status, expected_status, n_factors, "test_normalize_single_trajectory_values: status")
        call assert_equal_array_real(trajectory_norm, expected, n_timepoints*n_factors, 0.0_real64, &
                                     "test_normalize_single_trajectory_values: trajectory_norm", n_timepoints)
    end subroutine test_normalize_single_trajectory_values

    !> The one case whose range is not a power of two, kept from the old Fortran and R suites:
    !| factors [11, 12, 13, 14], [21, 22, 23, 24] and [31, 32, 33, 34] each have range 3, so all
    !| scale to [0, 1/3, 2/3, 1]. A correctly rounded (v - min)/3 is the double nearest 1/3 and
    !| 2/3, which is what 1/3.0 and 2/3.0 are; a tolerance of a few ulps anyway, because an
    !| optimizer may multiply by the rounded reciprocal 1/3 instead of dividing by 3.
    subroutine test_normalize_single_trajectory_range_of_three()
        integer(int32), parameter :: n_factors = 3, n_timepoints = 4
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints, n_factors)
        integer(int32) :: status(n_factors), expected_status(n_factors), ierr, i_factor

        trajectory(:, 1) = [11.0_real64, 12.0_real64, 13.0_real64, 14.0_real64]
        trajectory(:, 2) = [21.0_real64, 22.0_real64, 23.0_real64, 24.0_real64]
        trajectory(:, 3) = [31.0_real64, 32.0_real64, 33.0_real64, 34.0_real64]
        do i_factor = 1, n_factors
            expected(:, i_factor) = [0.0_real64, 1.0_real64/3.0_real64, 2.0_real64/3.0_real64, 1.0_real64]
        end do
        expected_status = ERR_OK
        trajectory_norm = UNWRITTEN

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_single_trajectory_range_of_three: ierr")
        call assert_equal_array_int(status, expected_status, n_factors, &
                                    "test_normalize_single_trajectory_range_of_three: status")
        call assert_equal_array_real(trajectory_norm, expected, n_timepoints*n_factors, 1.0e-15_real64, &
                                     "test_normalize_single_trajectory_range_of_three: trajectory_norm", n_timepoints)
    end subroutine test_normalize_single_trajectory_range_of_three

    !> One constant factor between two varying ones, over five time points:
    !| factor 1 [9, 1, 17, 5, 13]: min 1 (point 2), max 17 (point 3), range 16 ->
    !|   [8, 0, 16, 4, 12]/16 = [0.5, 0, 1, 0.25, 0.75];
    !| factor 2 [5, 5, 5, 5, 5]: constant;
    !| factor 3 [1, 2, 4, 8, 17]: range 16 -> [0, 1, 3, 7, 16]/16 = [0, 0.0625, 0.1875, 0.4375, 1].
    !| The constant factor must not disturb its neighbours, and only its own status is set.
    subroutine test_normalize_single_trajectory_constant_factor()
        integer(int32), parameter :: n_factors = 3, n_timepoints = 5
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints, n_factors)
        integer(int32) :: status(n_factors), expected_status(n_factors), ierr

        ! A constant series is zeros with ERR_DIVISION_BY_ZERO in its status only (see
        ! test_normalize_variable_timeseries_constant_series).
        trajectory(:, 1) = [9.0_real64, 1.0_real64, 17.0_real64, 5.0_real64, 13.0_real64]
        trajectory(:, 2) = 5.0_real64
        trajectory(:, 3) = [1.0_real64, 2.0_real64, 4.0_real64, 8.0_real64, 17.0_real64]
        expected(:, 1) = [0.5_real64, 0.0_real64, 1.0_real64, 0.25_real64, 0.75_real64]
        expected(:, 2) = 0.0_real64
        expected(:, 3) = [0.0_real64, 0.0625_real64, 0.1875_real64, 0.4375_real64, 1.0_real64]
        expected_status = [ERR_OK, ERR_DIVISION_BY_ZERO, ERR_OK]
        trajectory_norm = UNWRITTEN

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_single_trajectory_constant_factor: ierr")
        call assert_equal_array_int(status, expected_status, n_factors, &
                                    "test_normalize_single_trajectory_constant_factor: status")
        call assert_equal_array_real(trajectory_norm, expected, n_timepoints*n_factors, 0.0_real64, &
                                     "test_normalize_single_trajectory_constant_factor: trajectory_norm", n_timepoints)
    end subroutine test_normalize_single_trajectory_constant_factor

    !> n_timepoints = 1, the smallest valid size: every factor is a single point, so constant.
    subroutine test_normalize_single_trajectory_single_timepoint()
        integer(int32), parameter :: n_factors = 2, n_timepoints = 1
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints, n_factors)
        integer(int32) :: status(n_factors), expected_status(n_factors), ierr

        ! As in test_normalize_single_trajectory_constant_factor: zeros and ERR_DIVISION_BY_ZERO
        ! for every factor.
        trajectory(1, :) = [3.0_real64, -7.0_real64]
        expected = 0.0_real64
        expected_status = ERR_DIVISION_BY_ZERO
        trajectory_norm = UNWRITTEN

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_single_trajectory_single_timepoint: ierr")
        call assert_equal_array_int(status, expected_status, n_factors, &
                                    "test_normalize_single_trajectory_single_timepoint: status")
        call assert_equal_array_real(trajectory_norm, expected, n_timepoints*n_factors, 0.0_real64, &
                                     "test_normalize_single_trajectory_single_timepoint: trajectory_norm")
    end subroutine test_normalize_single_trajectory_single_timepoint

    !> Factor 1 [-2**1022, 0, 2**1022] has range 2**1023, just inside the overflow of max - min:
    !| [0, 0.5, 1]. Factor 2 [-huge, 0, huge] is valid input whose range overflows; scaled, it is
    !| [0, 0.5, 1] as well, 0 lying halfway.
    subroutine test_normalize_single_trajectory_extreme_magnitudes()
        integer(int32), parameter :: n_factors = 2, n_timepoints = 3
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints, n_factors)
        integer(int32) :: status(n_factors), expected_status(n_factors), ierr

        ! Regression: factor 2 used to come out as [0, 0, NaN] with status OK, as max - min
        ! overflowed to Inf (see test_normalize_variable_timeseries_range_overflows).
        trajectory(:, 1) = [-scale(1.0_real64, 1022), 0.0_real64, scale(1.0_real64, 1022)]
        trajectory(:, 2) = [-huge(1.0_real64), 0.0_real64, huge(1.0_real64)]
        expected(:, 1) = [0.0_real64, 0.5_real64, 1.0_real64]
        expected(:, 2) = [0.0_real64, 0.5_real64, 1.0_real64]
        expected_status = ERR_OK
        trajectory_norm = UNWRITTEN

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_single_trajectory_extreme_magnitudes: ierr")
        call assert_equal_array_int(status, expected_status, n_factors, &
                                    "test_normalize_single_trajectory_extreme_magnitudes: status")
        ! assert_equal_array_real lets a NaN through (NaN > tol is false), hence this check.
        call assert_no_nan_real(trajectory_norm, n_timepoints*n_factors, &
                                "test_normalize_single_trajectory_extreme_magnitudes: NaN in trajectory_norm")
        call assert_equal_array_real(trajectory_norm, expected, n_timepoints*n_factors, 0.0_real64, &
                                     "test_normalize_single_trajectory_extreme_magnitudes: trajectory_norm", n_timepoints)
    end subroutine test_normalize_single_trajectory_extreme_magnitudes

    !> Each dimension on its own: n_factors (argument 3) and n_timepoints (argument 4) are
    !| ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1. Just inside, a single factor
    !| [2, 4, 6] (range 4 -> [0, 0.5, 1]) is fine; n_timepoints = 1 is in
    !| `test_normalize_single_trajectory_single_timepoint`.
    subroutine test_normalize_single_trajectory_dimensions()
        real(real64) :: trajectory(3, 1), trajectory_norm(3, 1), expected(3, 1)
        integer(int32) :: status(1), ierr

        trajectory(:, 1) = [2.0_real64, 4.0_real64, 6.0_real64]

        call normalize_single_trajectory(trajectory, trajectory_norm, 0, 3, status, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_normalize_single_trajectory_dimensions: n_factors = 0", 3)
        call normalize_single_trajectory(trajectory, trajectory_norm, -1, 3, status, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_normalize_single_trajectory_dimensions: n_factors = -1", 3)
        call normalize_single_trajectory(trajectory, trajectory_norm, 1, 0, status, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_normalize_single_trajectory_dimensions: n_timepoints = 0", 4)
        call normalize_single_trajectory(trajectory, trajectory_norm, 1, -1, status, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_normalize_single_trajectory_dimensions: n_timepoints = -1", 4)

        expected(:, 1) = [0.0_real64, 0.5_real64, 1.0_real64]
        trajectory_norm = UNWRITTEN
        call normalize_single_trajectory(trajectory, trajectory_norm, 1, 3, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_single_trajectory_dimensions: n_factors = 1, ierr")
        call assert_err(status(1), ERR_OK, "test_normalize_single_trajectory_dimensions: n_factors = 1, status")
        call assert_equal_array_real(trajectory_norm, expected, 3, 0.0_real64, &
                                     "test_normalize_single_trajectory_dimensions: n_factors = 1, trajectory_norm")
    end subroutine test_normalize_single_trajectory_dimensions

    !> NaN, +Inf and -Inf in `trajectory` (argument 1) are ERR_NAN_INF. Each sits in the last
    !| element of the last factor, so the check must cover all n_timepoints*n_factors values.
    subroutine test_normalize_single_trajectory_rejects_nan_and_inf()
        integer(int32), parameter :: n_factors = 2, n_timepoints = 3
        character(*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors), bad(3)
        integer(int32) :: status(n_factors), ierr, i_bad

        trajectory(:, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectory(:, 2) = [4.0_real64, 5.0_real64, 6.0_real64]
        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        do i_bad = 1, size(bad)
            trajectory(n_timepoints, n_factors) = bad(i_bad)
            call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_normalize_single_trajectory_rejects_nan_and_inf: "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_normalize_single_trajectory_rejects_nan_and_inf

end module mod_test_normalize_single_trajectory
