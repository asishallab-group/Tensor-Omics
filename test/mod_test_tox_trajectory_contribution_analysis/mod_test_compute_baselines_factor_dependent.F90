!> The `compute_baselines_factor_dependent` cases: the baseline of each mode for a factor and a
!| dependent, each taken from its own series, a single time point, extreme magnitudes, the mode,
!| the dimension and the NaN/Inf check.
!|
!| Every mean divides by a power of two, so every expected value is exact and needs no tolerance.
module mod_test_compute_baselines_factor_dependent
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_baselines_factor_dependent
    use tox_trajectory_contribution_analysis_impl, only: MODE_BASELINE_RAW, MODE_BASELINE_MIN, MODE_BASELINE_MEAN
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_baselines_factor_dependent

    integer(int32), parameter :: N_TIMEPOINTS = 4
    !> min -1 at point 2, mean 8/4 = 2
    real(real64), parameter :: SERIES_FACTOR(N_TIMEPOINTS) = [3.0_real64, -1.0_real64, 4.0_real64, 2.0_real64]
    !> min 2 at point 3, mean 22/4 = 5.5
    real(real64), parameter :: SERIES_DEPENDENT(N_TIMEPOINTS) = [5.0_real64, 9.0_real64, 2.0_real64, 6.0_real64]

contains

    !> Every case of `compute_baselines_factor_dependent`.
    function get_all_tests_compute_baselines_factor_dependent() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(8))
        all_tests(1) = test_case("test_compute_baselines_factor_dependent_raw", &
                                 test_compute_baselines_factor_dependent_raw)
        all_tests(2) = test_case("test_compute_baselines_factor_dependent_min", &
                                 test_compute_baselines_factor_dependent_min)
        all_tests(3) = test_case("test_compute_baselines_factor_dependent_mean", &
                                 test_compute_baselines_factor_dependent_mean)
        all_tests(4) = test_case("test_compute_baselines_factor_dependent_single_timepoint", &
                                 test_compute_baselines_factor_dependent_single_timepoint)
        all_tests(5) = test_case("test_compute_baselines_factor_dependent_extreme_magnitudes", &
                                 test_compute_baselines_factor_dependent_extreme_magnitudes)
        all_tests(6) = test_case("test_compute_baselines_factor_dependent_invalid_mode", &
                                 test_compute_baselines_factor_dependent_invalid_mode)
        all_tests(7) = test_case("test_compute_baselines_factor_dependent_dimensions", &
                                 test_compute_baselines_factor_dependent_dimensions)
        all_tests(8) = test_case("test_compute_baselines_factor_dependent_rejects_nan_and_inf", &
                                 test_compute_baselines_factor_dependent_rejects_nan_and_inf)
    end function get_all_tests_compute_baselines_factor_dependent

    !> Calls the procedure on `factor` and `dependent` in `baseline_mode` and checks both
    !| baselines exactly, with `ierr` OK. The outputs start at UNWRITTEN, so a baseline the
    !| procedure never writes shows.
    subroutine check_baselines(n_points, factor, dependent, baseline_mode, expected_factor_baseline, &
                               expected_dependent_baseline, msg)
        integer(int32), intent(in) :: n_points, baseline_mode
        real(real64), intent(in) :: factor(n_points), dependent(n_points)
        real(real64), intent(in) :: expected_factor_baseline, expected_dependent_baseline
        character(*), intent(in) :: msg
        real(real64) :: factor_baseline, dependent_baseline
        integer(int32) :: ierr

        factor_baseline = UNWRITTEN
        dependent_baseline = UNWRITTEN
        call compute_baselines_factor_dependent(n_points, factor, dependent, baseline_mode, &
                                                factor_baseline, dependent_baseline, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_real(factor_baseline, expected_factor_baseline, 0.0_real64, msg//": factor_baseline")
        call assert_equal_real(dependent_baseline, expected_dependent_baseline, 0.0_real64, &
                               msg//": dependent_baseline")
    end subroutine check_baselines

    !> MODE_BASELINE_RAW measures from zero: both baselines are 0, whatever the series.
    subroutine test_compute_baselines_factor_dependent_raw()
        call check_baselines(N_TIMEPOINTS, SERIES_FACTOR, SERIES_DEPENDENT, MODE_BASELINE_RAW, 0.0_real64, 0.0_real64, &
                             "test_compute_baselines_factor_dependent_raw")
    end subroutine test_compute_baselines_factor_dependent_raw

    !> MODE_BASELINE_MIN takes each series' own minimum: factor [3, -1, 4, 2] -> -1 (point 2),
    !| dependent [5, 9, 2, 6] -> 2 (point 3). Neither sits at an end, and they differ, so a
    !| baseline taken from the other series, or from an end point, shows.
    subroutine test_compute_baselines_factor_dependent_min()
        call check_baselines(N_TIMEPOINTS, SERIES_FACTOR, SERIES_DEPENDENT, MODE_BASELINE_MIN, -1.0_real64, 2.0_real64, &
                             "test_compute_baselines_factor_dependent_min")
    end subroutine test_compute_baselines_factor_dependent_min

    !> MODE_BASELINE_MEAN takes each series' own mean: factor (3 - 1 + 4 + 2)/4 = 2,
    !| dependent (5 + 9 + 2 + 6)/4 = 5.5.
    subroutine test_compute_baselines_factor_dependent_mean()
        call check_baselines(N_TIMEPOINTS, SERIES_FACTOR, SERIES_DEPENDENT, MODE_BASELINE_MEAN, 2.0_real64, 5.5_real64, &
                             "test_compute_baselines_factor_dependent_mean")
    end subroutine test_compute_baselines_factor_dependent_mean

    !> n_timepoints = 1, the smallest valid size: the minimum and the mean of a single point are
    !| the point itself, factor [7] and dependent [-3].
    subroutine test_compute_baselines_factor_dependent_single_timepoint()
        real(real64) :: factor(1), dependent(1)

        factor = 7.0_real64
        dependent = -3.0_real64
        call check_baselines(1, factor, dependent, MODE_BASELINE_RAW, 0.0_real64, 0.0_real64, &
                             "test_compute_baselines_factor_dependent_single_timepoint: raw")
        call check_baselines(1, factor, dependent, MODE_BASELINE_MIN, 7.0_real64, -3.0_real64, &
                             "test_compute_baselines_factor_dependent_single_timepoint: min")
        call check_baselines(1, factor, dependent, MODE_BASELINE_MEAN, 7.0_real64, -3.0_real64, &
                             "test_compute_baselines_factor_dependent_single_timepoint: mean")
    end subroutine test_compute_baselines_factor_dependent_single_timepoint

    !> The largest finite magnitudes are valid input: factor [-huge, huge] has the minimum -huge
    !| and the mean 0; dependent [huge, 0] has the minimum 0 and the mean huge/2 (halving is exact).
    subroutine test_compute_baselines_factor_dependent_extreme_magnitudes()
        real(real64) :: factor(2), dependent(2)

        factor = [-huge(1.0_real64), huge(1.0_real64)]
        dependent = [huge(1.0_real64), 0.0_real64]
        call check_baselines(2, factor, dependent, MODE_BASELINE_MIN, -huge(1.0_real64), 0.0_real64, &
                             "test_compute_baselines_factor_dependent_extreme_magnitudes: min")
        call check_baselines(2, factor, dependent, MODE_BASELINE_MEAN, 0.0_real64, huge(1.0_real64)/2.0_real64, &
                             "test_compute_baselines_factor_dependent_extreme_magnitudes: mean")
    end subroutine test_compute_baselines_factor_dependent_extreme_magnitudes

    !> `baseline_mode` (argument 4) is one of the three MODE_BASELINE_* values (1 to 3); the value
    !| just outside on either side, 0 and 4, is ERR_INVALID_INPUT.
    subroutine test_compute_baselines_factor_dependent_invalid_mode()
        real(real64) :: factor_baseline, dependent_baseline
        integer(int32) :: ierr

        call compute_baselines_factor_dependent(N_TIMEPOINTS, SERIES_FACTOR, SERIES_DEPENDENT, 0_int32, &
                                                factor_baseline, dependent_baseline, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_baselines_factor_dependent_invalid_mode: 0", 4)
        call compute_baselines_factor_dependent(N_TIMEPOINTS, SERIES_FACTOR, SERIES_DEPENDENT, 4_int32, &
                                                factor_baseline, dependent_baseline, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_baselines_factor_dependent_invalid_mode: 4", 4)
    end subroutine test_compute_baselines_factor_dependent_invalid_mode

    !> `n_timepoints` (argument 1) is ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1; 1 is in
    !| `test_compute_baselines_factor_dependent_single_timepoint`.
    subroutine test_compute_baselines_factor_dependent_dimensions()
        real(real64) :: factor_baseline, dependent_baseline
        integer(int32) :: ierr

        call compute_baselines_factor_dependent(0, SERIES_FACTOR, SERIES_DEPENDENT, MODE_BASELINE_MEAN, &
                                                factor_baseline, dependent_baseline, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_baselines_factor_dependent_dimensions: 0", 1)
        call compute_baselines_factor_dependent(-1, SERIES_FACTOR, SERIES_DEPENDENT, MODE_BASELINE_MEAN, &
                                                factor_baseline, dependent_baseline, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_baselines_factor_dependent_dimensions: -1", 1)
    end subroutine test_compute_baselines_factor_dependent_dimensions

    !> NaN, +Inf and -Inf in `factor` (argument 2) or in `dependent` (argument 3) are ERR_NAN_INF.
    !| Each sits in the last point, so the check must cover all n_timepoints values. The mode is
    !| MODE_BASELINE_RAW, whose baseline does not depend on the values at all: only the input
    !| check can refuse them.
    subroutine test_compute_baselines_factor_dependent_rejects_nan_and_inf()
        real(real64) :: factor(N_TIMEPOINTS), dependent(N_TIMEPOINTS), factor_baseline, dependent_baseline
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            factor = SERIES_FACTOR
            dependent = SERIES_DEPENDENT
            factor(N_TIMEPOINTS) = non_finite_value(i_value)
            call compute_baselines_factor_dependent(N_TIMEPOINTS, factor, dependent, MODE_BASELINE_RAW, &
                                                    factor_baseline, dependent_baseline, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_baselines_factor_dependent_rejects_nan_and_inf: factor " &
                            //trim(NON_FINITE_NAMES(i_value)), 2)

            factor = SERIES_FACTOR
            dependent(N_TIMEPOINTS) = non_finite_value(i_value)
            call compute_baselines_factor_dependent(N_TIMEPOINTS, factor, dependent, MODE_BASELINE_RAW, &
                                                    factor_baseline, dependent_baseline, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_baselines_factor_dependent_rejects_nan_and_inf: dependent " &
                            //trim(NON_FINITE_NAMES(i_value)), 3)
        end do
    end subroutine test_compute_baselines_factor_dependent_rejects_nan_and_inf

end module mod_test_compute_baselines_factor_dependent
