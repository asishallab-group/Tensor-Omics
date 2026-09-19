!> The `normalize_all_trajectories` cases: every (factor, sample) series of
!| `trajectories(n_factors, n_samples, n_timepoints)` scaled to [0, 1] over time -- the last,
!| strided dimension -- on its own; constant series among varying ones, with `status` per
!| (factor, sample); the expert entry point; extreme magnitudes; each dimension; the NaN/Inf check.
!|
!| Every series has a power-of-two range (max - min), so each (v - min)/range is exact, even
!| where an optimizer multiplies by the reciprocal of the range: no case needs a tolerance.
module mod_test_normalize_all_trajectories
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_trajectory_normalization, only: normalize_all_trajectories, normalize_all_trajectories_expert
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_normalize_all_trajectories

    !> What `trajectories_norm` holds before a call, so that an output the procedure never writes shows.
    real(real64), parameter :: UNWRITTEN = -1.0_real64

contains

    !> Every case of `normalize_all_trajectories`.
    function get_all_tests_normalize_all_trajectories() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(7))
        all_tests(1) = test_case("test_normalize_all_trajectories_values", test_normalize_all_trajectories_values)
        all_tests(2) = test_case("test_normalize_all_trajectories_constant_series", &
                                 test_normalize_all_trajectories_constant_series)
        all_tests(3) = test_case("test_normalize_all_trajectories_expert", test_normalize_all_trajectories_expert)
        all_tests(4) = test_case("test_normalize_all_trajectories_single_of_each", &
                                 test_normalize_all_trajectories_single_of_each)
        all_tests(5) = test_case("test_normalize_all_trajectories_extreme_magnitudes", &
                                 test_normalize_all_trajectories_extreme_magnitudes)
        all_tests(6) = test_case("test_normalize_all_trajectories_dimensions", test_normalize_all_trajectories_dimensions)
        all_tests(7) = test_case("test_normalize_all_trajectories_rejects_nan_and_inf", &
                                 test_normalize_all_trajectories_rejects_nan_and_inf)
    end function get_all_tests_normalize_all_trajectories

    !> Fills the 2 factors x 2 samples x 4 time points of `test_normalize_all_trajectories_values`
    !| and `test_normalize_all_trajectories_expert`, and what they scale to:
    !| (1, 1) [1, 2, 3, 5]:     min 1, range 4   -> [0, 1, 2, 4]/4    = [0, 0.25, 0.5, 1]
    !| (2, 1) [40, 30, 20, 8]:  min 8, range 32  -> [32, 22, 12, 0]/32 = [1, 0.6875, 0.375, 0]
    !| (1, 2) [9, 7, 5, 6]:     min 5, range 4   -> [4, 2, 0, 1]/4    = [1, 0.5, 0, 0.25]
    !| (2, 2) [2, 10, 4, 6]:    min 2, range 8   -> [0, 8, 2, 4]/8    = [0, 1, 0.25, 0.5]
    !| Every series differs, so scaling along the wrong dimension, or mixing up factors and
    !| samples, changes the result.
    pure subroutine fill_two_by_two(trajectories, expected)
        real(real64), intent(out) :: trajectories(2, 2, 4)
            !! (n_factors, n_samples, n_timepoints) input
        real(real64), intent(out) :: expected(2, 2, 4)
            !! Its hand-derived scaling

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64, 3.0_real64, 5.0_real64]
        trajectories(2, 1, :) = [40.0_real64, 30.0_real64, 20.0_real64, 8.0_real64]
        trajectories(1, 2, :) = [9.0_real64, 7.0_real64, 5.0_real64, 6.0_real64]
        trajectories(2, 2, :) = [2.0_real64, 10.0_real64, 4.0_real64, 6.0_real64]
        expected(1, 1, :) = [0.0_real64, 0.25_real64, 0.5_real64, 1.0_real64]
        expected(2, 1, :) = [1.0_real64, 0.6875_real64, 0.375_real64, 0.0_real64]
        expected(1, 2, :) = [1.0_real64, 0.5_real64, 0.0_real64, 0.25_real64]
        expected(2, 2, :) = [0.0_real64, 1.0_real64, 0.25_real64, 0.5_real64]
    end subroutine fill_two_by_two

    !> Two factors, two samples, four time points; values in `fill_two_by_two`.
    subroutine test_normalize_all_trajectories_values()
        integer(int32), parameter :: n_factors = 2, n_samples = 2, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints), expected(n_factors, n_samples, n_timepoints)
        integer(int32) :: status(n_factors, n_samples), expected_status(n_factors, n_samples), ierr

        call fill_two_by_two(trajectories, expected)
        expected_status = ERR_OK
        trajectories_norm = UNWRITTEN

        call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_values: ierr")
        call assert_equal_array_int(status, expected_status, n_factors*n_samples, &
                                    "test_normalize_all_trajectories_values: status")
        call assert_equal_array_real(trajectories_norm, expected, n_factors*n_samples*n_timepoints, 0.0_real64, &
                                     "test_normalize_all_trajectories_values: trajectories_norm")
    end subroutine test_normalize_all_trajectories_values

    !> 3 factors x 4 samples x 5 time points, factors 1 and 3 constant at 7 and factor 2 the ramp
    !| [0, 1, 2, 3, 4] (range 4 -> [0, 0.25, 0.5, 0.75, 1]) in every sample -- except sample 2's
    !| factor 1, [7, 7, 7, 7, 11] (range 4 -> [0, 0, 0, 0, 1]). So `status` is set per
    !| (factor, sample), not per factor, and the constant series leave the others alone.
    subroutine test_normalize_all_trajectories_constant_series()
        integer(int32), parameter :: n_factors = 3, n_samples = 4, n_timepoints = 5
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints), ramp(n_timepoints), ramp_norm(n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints), expected(n_factors, n_samples, n_timepoints)
        integer(int32) :: status(n_factors, n_samples), expected_status(n_factors, n_samples), ierr, i_sample

        ! QUESTION: what a constant series should give is open (see
        ! test_normalize_variable_timeseries_constant_series); today zeros and ERR_DIVISION_BY_ZERO.
        ramp = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        ramp_norm = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        trajectories = 7.0_real64
        expected = 0.0_real64
        do i_sample = 1, n_samples
            trajectories(2, i_sample, :) = ramp
            expected(2, i_sample, :) = ramp_norm
        end do
        trajectories(1, 2, n_timepoints) = 11.0_real64
        expected(1, 2, n_timepoints) = 1.0_real64

        expected_status = ERR_DIVISION_BY_ZERO
        expected_status(2, :) = ERR_OK
        expected_status(1, 2) = ERR_OK
        trajectories_norm = UNWRITTEN

        call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_constant_series: ierr")
        call assert_equal_array_int(status, expected_status, n_factors*n_samples, &
                                    "test_normalize_all_trajectories_constant_series: status", n_factors)
        call assert_equal_array_real(trajectories_norm, expected, n_factors*n_samples*n_timepoints, 0.0_real64, &
                                     "test_normalize_all_trajectories_constant_series: trajectories_norm")
    end subroutine test_normalize_all_trajectories_constant_series

    !> The expert entry point takes the work arrays from the caller and gives the same scaling
    !| (`fill_two_by_two`); the work arrays hold garbage before the call.
    subroutine test_normalize_all_trajectories_expert()
        integer(int32), parameter :: n_factors = 2, n_samples = 2, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints), expected(n_factors, n_samples, n_timepoints)
        real(real64) :: tmp_series(n_timepoints), tmp_series_norm(n_timepoints)
        integer(int32) :: status(n_factors, n_samples), expected_status(n_factors, n_samples), ierr

        call fill_two_by_two(trajectories, expected)
        expected_status = ERR_OK
        trajectories_norm = UNWRITTEN
        tmp_series = 1.0e300_real64
        tmp_series_norm = -1.0e300_real64

        call normalize_all_trajectories_expert(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, &
                                               tmp_series, tmp_series_norm, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_expert: ierr")
        call assert_equal_array_int(status, expected_status, n_factors*n_samples, &
                                    "test_normalize_all_trajectories_expert: status")
        call assert_equal_array_real(trajectories_norm, expected, n_factors*n_samples*n_timepoints, 0.0_real64, &
                                     "test_normalize_all_trajectories_expert: trajectories_norm")
    end subroutine test_normalize_all_trajectories_expert

    !> Each dimension at 1, the smallest valid size, with the others larger:
    !| n_factors = 1:    3 samples [0, 2, 4], [4, 0, 2], [2, 4, 0], range 4 -> [0, 0.5, 1], [1, 0, 0.5], [0.5, 1, 0];
    !| n_samples = 1:    3 factors, the same three series;
    !| n_timepoints = 1: every series a single point, constant.
    subroutine test_normalize_all_trajectories_single_of_each()
        real(real64) :: series(3, 3), series_norm(3, 3)
        real(real64) :: one_factor(1, 3, 3), one_factor_norm(1, 3, 3), one_factor_expected(1, 3, 3)
        real(real64) :: one_sample(3, 1, 3), one_sample_norm(3, 1, 3), one_sample_expected(3, 1, 3)
        real(real64) :: one_timepoint(3, 3, 1), one_timepoint_norm(3, 3, 1), one_timepoint_expected(3, 3, 1)
        integer(int32) :: status_one_factor(1, 3), status_one_sample(3, 1), status_one_timepoint(3, 3)
        integer(int32) :: expected_ok(3, 3), expected_constant(3, 3), ierr, i_series

        ! series(:, i) is the i-th series over time, series_norm(:, i) its scaling
        series(:, 1) = [0.0_real64, 2.0_real64, 4.0_real64]
        series(:, 2) = [4.0_real64, 0.0_real64, 2.0_real64]
        series(:, 3) = [2.0_real64, 4.0_real64, 0.0_real64]
        series_norm(:, 1) = [0.0_real64, 0.5_real64, 1.0_real64]
        series_norm(:, 2) = [1.0_real64, 0.0_real64, 0.5_real64]
        series_norm(:, 3) = [0.5_real64, 1.0_real64, 0.0_real64]
        expected_ok = ERR_OK
        expected_constant = ERR_DIVISION_BY_ZERO

        do i_series = 1, 3
            one_factor(1, i_series, :) = series(:, i_series)
            one_factor_expected(1, i_series, :) = series_norm(:, i_series)
            one_sample(i_series, 1, :) = series(:, i_series)
            one_sample_expected(i_series, 1, :) = series_norm(:, i_series)
        end do
        one_factor_norm = UNWRITTEN
        call normalize_all_trajectories(one_factor, one_factor_norm, 1, 3, 3, status_one_factor, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_single_of_each: n_factors = 1, ierr")
        call assert_equal_array_int(status_one_factor, expected_ok, 3, &
                                    "test_normalize_all_trajectories_single_of_each: n_factors = 1, status")
        call assert_equal_array_real(one_factor_norm, one_factor_expected, 9, 0.0_real64, &
                                     "test_normalize_all_trajectories_single_of_each: n_factors = 1, trajectories_norm")

        one_sample_norm = UNWRITTEN
        call normalize_all_trajectories(one_sample, one_sample_norm, 3, 1, 3, status_one_sample, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_single_of_each: n_samples = 1, ierr")
        call assert_equal_array_int(status_one_sample, expected_ok, 3, &
                                    "test_normalize_all_trajectories_single_of_each: n_samples = 1, status")
        call assert_equal_array_real(one_sample_norm, one_sample_expected, 9, 0.0_real64, &
                                     "test_normalize_all_trajectories_single_of_each: n_samples = 1, trajectories_norm")

        ! QUESTION: a single time point is a constant series; today zeros and ERR_DIVISION_BY_ZERO
        ! in every status (see test_normalize_variable_timeseries_constant_series).
        one_timepoint(:, :, 1) = series
        one_timepoint_expected = 0.0_real64
        one_timepoint_norm = UNWRITTEN
        call normalize_all_trajectories(one_timepoint, one_timepoint_norm, 3, 3, 1, status_one_timepoint, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_single_of_each: n_timepoints = 1, ierr")
        call assert_equal_array_int(status_one_timepoint, expected_constant, 9, &
                                    "test_normalize_all_trajectories_single_of_each: n_timepoints = 1, status")
        call assert_equal_array_real(one_timepoint_norm, one_timepoint_expected, 9, 0.0_real64, &
                                     "test_normalize_all_trajectories_single_of_each: n_timepoints = 1, trajectories_norm")
    end subroutine test_normalize_all_trajectories_single_of_each

    !> Series (1, 1) [-2**1022, 0, 2**1022] has range 2**1023, just inside the overflow of max - min:
    !| [0, 0.5, 1]. Series (2, 1) [-huge, 0, huge] is valid input whose range overflows; scaled, it
    !| is [0, 0.5, 1] as well, 0 lying halfway.
    subroutine test_normalize_all_trajectories_extreme_magnitudes()
        integer(int32), parameter :: n_factors = 2, n_samples = 1, n_timepoints = 3
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints), expected(n_factors, n_samples, n_timepoints)
        integer(int32) :: status(n_factors, n_samples), expected_status(n_factors, n_samples), ierr

        ! Regression: series (2, 1) used to come out as [0, 0, NaN] with status OK, as max - min
        ! overflowed to Inf (see test_normalize_variable_timeseries_range_overflows).
        trajectories(1, 1, :) = [-scale(1.0_real64, 1022), 0.0_real64, scale(1.0_real64, 1022)]
        trajectories(2, 1, :) = [-huge(1.0_real64), 0.0_real64, huge(1.0_real64)]
        expected(1, 1, :) = [0.0_real64, 0.5_real64, 1.0_real64]
        expected(2, 1, :) = [0.0_real64, 0.5_real64, 1.0_real64]
        expected_status = ERR_OK
        trajectories_norm = UNWRITTEN

        call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, status, ierr)
        call assert_err(ierr, ERR_OK, "test_normalize_all_trajectories_extreme_magnitudes: ierr")
        call assert_equal_array_int(status, expected_status, n_factors*n_samples, &
                                    "test_normalize_all_trajectories_extreme_magnitudes: status")
        ! assert_equal_array_real lets a NaN through (NaN > tol is false), hence this check.
        call assert_no_nan_real(trajectories_norm, n_factors*n_samples*n_timepoints, &
                                "test_normalize_all_trajectories_extreme_magnitudes: NaN in trajectories_norm")
        call assert_equal_array_real(trajectories_norm, expected, n_factors*n_samples*n_timepoints, 0.0_real64, &
                                     "test_normalize_all_trajectories_extreme_magnitudes: trajectories_norm")
    end subroutine test_normalize_all_trajectories_extreme_magnitudes

    !> Each dimension on its own: n_factors (argument 3), n_samples (argument 4) and n_timepoints
    !| (argument 5) are ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1. The last valid size, 1,
    !| is in `test_normalize_all_trajectories_single_of_each`.
    subroutine test_normalize_all_trajectories_dimensions()
        real(real64) :: trajectories(1, 1, 2), trajectories_norm(1, 1, 2)
        integer(int32) :: status(1, 1), ierr

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64]

        call normalize_all_trajectories(trajectories, trajectories_norm, 0, 1, 2, status, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_normalize_all_trajectories_dimensions: n_factors = 0", 3)
        call normalize_all_trajectories(trajectories, trajectories_norm, -1, 1, 2, status, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_normalize_all_trajectories_dimensions: n_factors = -1", 3)
        call normalize_all_trajectories(trajectories, trajectories_norm, 1, 0, 2, status, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_normalize_all_trajectories_dimensions: n_samples = 0", 4)
        call normalize_all_trajectories(trajectories, trajectories_norm, 1, -1, 2, status, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_normalize_all_trajectories_dimensions: n_samples = -1", 4)
        call normalize_all_trajectories(trajectories, trajectories_norm, 1, 1, 0, status, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_normalize_all_trajectories_dimensions: n_timepoints = 0", 5)
        call normalize_all_trajectories(trajectories, trajectories_norm, 1, 1, -1, status, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_normalize_all_trajectories_dimensions: n_timepoints = -1", 5)
    end subroutine test_normalize_all_trajectories_dimensions

    !> NaN, +Inf and -Inf in `trajectories` (argument 1) are ERR_NAN_INF. Each sits in the last
    !| element in memory, so the check must cover all n_factors*n_samples*n_timepoints values.
    subroutine test_normalize_all_trajectories_rejects_nan_and_inf()
        integer(int32), parameter :: n_factors = 2, n_samples = 2, n_timepoints = 4
        character(*), parameter :: bad_names(3) = ["NaN ", "+Inf", "-Inf"]
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints), expected(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints), bad(3)
        integer(int32) :: status(n_factors, n_samples), ierr, i_bad

        call fill_two_by_two(trajectories, expected)
        bad(1) = ieee_value(1.0_real64, ieee_quiet_nan)
        bad(2) = ieee_value(1.0_real64, ieee_positive_inf)
        bad(3) = ieee_value(1.0_real64, ieee_negative_inf)

        do i_bad = 1, size(bad)
            trajectories(n_factors, n_samples, n_timepoints) = bad(i_bad)
            call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, status, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_normalize_all_trajectories_rejects_nan_and_inf: "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_normalize_all_trajectories_rejects_nan_and_inf

end module mod_test_normalize_all_trajectories
