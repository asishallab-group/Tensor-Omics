!> The `compute_velocity_acceleration_contributions` cases: the contributions of every
!| (factor, dependent) pair of factors in every sample, computed on the velocities (first
!| differences) and the accelerations (second differences) of the series. A velocity series has
!| one point fewer than the positions and its contributions start at time point 2; an
!| acceleration series has two fewer and starts at 3; the points before are 0.
!|
!| Each case derives the velocities and accelerations by hand and multiplies their deviations;
!| all values are small integers or multiples of 1/4 and every mean divides by 1 or 2, so every
!| expected value is exact. A contribution is symmetric in its two series, each measured from
!| its own baseline, so the pair (a, b) always equals (b, a): the cases cannot tell the factor
!| axis from the dependent axis, only both from the time and sample axes.
module mod_test_compute_velocity_acceleration_contributions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_velocity_acceleration_contributions, &
                                                    compute_velocity_acceleration_contributions_expert
    use tox_trajectory_contribution_analysis_impl, only: MODE_BASELINE_RAW, MODE_BASELINE_MIN, MODE_BASELINE_MEAN
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_velocity_acceleration_contributions

contains

    !> Every case of `compute_velocity_acceleration_contributions`.
    function get_all_tests_compute_velocity_acceleration_contributions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(9))
        all_tests(1) = test_case("test_compute_velocity_acceleration_contributions_raw", &
                                 test_compute_velocity_acceleration_contributions_raw)
        all_tests(2) = test_case("test_compute_velocity_acceleration_contributions_expert", &
                                 test_compute_velocity_acceleration_contributions_expert)
        all_tests(3) = test_case("test_compute_velocity_acceleration_contributions_mean", &
                                 test_compute_velocity_acceleration_contributions_mean)
        all_tests(4) = test_case("test_compute_velocity_acceleration_contributions_min", &
                                 test_compute_velocity_acceleration_contributions_min)
        all_tests(5) = test_case("test_compute_velocity_acceleration_contributions_two_points", &
                                 test_compute_velocity_acceleration_contributions_two_points)
        all_tests(6) = test_case("test_compute_velocity_acceleration_contributions_one_point", &
                                 test_compute_velocity_acceleration_contributions_one_point)
        all_tests(7) = test_case("test_compute_velocity_acceleration_contributions_dimensions", &
                                 test_compute_velocity_acceleration_contributions_dimensions)
        all_tests(8) = test_case("test_compute_velocity_acceleration_contributions_bad_mode", &
                                 test_compute_velocity_acceleration_contributions_bad_mode)
        all_tests(9) = test_case("test_compute_velocity_acceleration_contributions_nan_inf", &
                                 test_compute_velocity_acceleration_contributions_nan_inf)
    end function get_all_tests_compute_velocity_acceleration_contributions

    !> Checks `ierr` OK and the four outputs exactly.
    subroutine check_outputs(n_factors, n_samples, n_timepoints, ierr, contrib_velocity, velocity_series, &
                             contrib_acceleration, acceleration_series, expected_contrib_velocity, &
                             expected_velocity_series, expected_contrib_acceleration, expected_acceleration_series, msg)
        integer(int32), intent(in) :: n_factors, n_samples, n_timepoints, ierr
        real(real64), intent(in) :: contrib_velocity(n_factors, n_factors, n_samples)
        real(real64), intent(in) :: velocity_series(n_timepoints, n_factors, n_factors, n_samples)
        real(real64), intent(in) :: contrib_acceleration(n_factors, n_factors, n_samples)
        real(real64), intent(in) :: acceleration_series(n_timepoints, n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_contrib_velocity(n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_velocity_series(n_timepoints, n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_contrib_acceleration(n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_acceleration_series(n_timepoints, n_factors, n_factors, n_samples)
        character(*), intent(in) :: msg
        integer(int32) :: n_contribs, n_series

        n_contribs = n_factors*n_factors*n_samples
        n_series = n_timepoints*n_contribs
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(contrib_velocity, expected_contrib_velocity, n_contribs, 0.0_real64, &
                                     msg//": contrib_velocity", n_factors)
        call assert_equal_array_real(velocity_series, expected_velocity_series, n_series, 0.0_real64, &
                                     msg//": velocity_contribution_series", n_timepoints)
        call assert_equal_array_real(contrib_acceleration, expected_contrib_acceleration, n_contribs, 0.0_real64, &
                                     msg//": contrib_acceleration", n_factors)
        call assert_equal_array_real(acceleration_series, expected_acceleration_series, n_series, 0.0_real64, &
                                     msg//": acceleration_contribution_series", n_timepoints)
    end subroutine check_outputs

    !> Calls the plain tier with every output at UNWRITTEN beforehand, so an element it never
    !| writes (the leading zeros of a series included) shows, and checks the result.
    subroutine check_plain(n_factors, n_samples, n_timepoints, trajectories, baseline_mode, &
                           expected_contrib_velocity, expected_velocity_series, expected_contrib_acceleration, &
                           expected_acceleration_series, msg)
        integer(int32), intent(in) :: n_factors, n_samples, n_timepoints, baseline_mode
        real(real64), intent(in) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64), intent(in) :: expected_contrib_velocity(n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_velocity_series(n_timepoints, n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_contrib_acceleration(n_factors, n_factors, n_samples)
        real(real64), intent(in) :: expected_acceleration_series(n_timepoints, n_factors, n_factors, n_samples)
        character(*), intent(in) :: msg
        real(real64) :: contrib_velocity(n_factors, n_factors, n_samples)
        real(real64) :: velocity_series(n_timepoints, n_factors, n_factors, n_samples)
        real(real64) :: contrib_acceleration(n_factors, n_factors, n_samples)
        real(real64) :: acceleration_series(n_timepoints, n_factors, n_factors, n_samples)
        integer(int32) :: ierr

        contrib_velocity = UNWRITTEN
        velocity_series = UNWRITTEN
        contrib_acceleration = UNWRITTEN
        acceleration_series = UNWRITTEN
        call compute_velocity_acceleration_contributions(trajectories, n_factors, n_samples, n_timepoints, &
                                                         baseline_mode, contrib_velocity, velocity_series, &
                                                         contrib_acceleration, acceleration_series, ierr)
        call check_outputs(n_factors, n_samples, n_timepoints, ierr, contrib_velocity, velocity_series, &
                           contrib_acceleration, acceleration_series, expected_contrib_velocity, &
                           expected_velocity_series, expected_contrib_acceleration, expected_acceleration_series, msg)
    end subroutine check_plain

    !> Kept from the old suites: two factors in one sample over four time points, x1
    !| [1, 3, 6, 10] and x2 [1, 2, 2, 1], with MODE_BASELINE_RAW. Velocities x1 [2, 3, 4] and x2
    !| [1, 0, -1]; accelerations x1 [1, 1] and x2 [-1, -1].
    !| Velocity: (x1, x1) [4, 9, 16] = 29, (x1, x2) and (x2, x1) [2, 0, -4] = -2, (x2, x2)
    !| [1, 0, 1] = 2. Acceleration: (x1, x1) [1, 1] = 2, (x1, x2) and (x2, x1) [-1, -1] = -2,
    !| (x2, x2) [1, 1] = 2.
    subroutine fill_raw_case(trajectories, expected_contrib_velocity, expected_velocity_series, &
                             expected_contrib_acceleration, expected_acceleration_series)
        real(real64), intent(out) :: trajectories(2, 1, 4)
        real(real64), intent(out) :: expected_contrib_velocity(2, 2, 1), expected_velocity_series(4, 2, 2, 1)
        real(real64), intent(out) :: expected_contrib_acceleration(2, 2, 1), expected_acceleration_series(4, 2, 2, 1)

        trajectories(1, 1, :) = [1.0_real64, 3.0_real64, 6.0_real64, 10.0_real64]
        trajectories(2, 1, :) = [1.0_real64, 2.0_real64, 2.0_real64, 1.0_real64]

        expected_contrib_velocity(:, 1, 1) = [29.0_real64, -2.0_real64]
        expected_contrib_velocity(:, 2, 1) = [-2.0_real64, 2.0_real64]
        expected_velocity_series(:, 1, 1, 1) = [0.0_real64, 4.0_real64, 9.0_real64, 16.0_real64]
        expected_velocity_series(:, 2, 1, 1) = [0.0_real64, 2.0_real64, 0.0_real64, -4.0_real64]
        expected_velocity_series(:, 1, 2, 1) = [0.0_real64, 2.0_real64, 0.0_real64, -4.0_real64]
        expected_velocity_series(:, 2, 2, 1) = [0.0_real64, 1.0_real64, 0.0_real64, 1.0_real64]

        expected_contrib_acceleration(:, 1, 1) = [2.0_real64, -2.0_real64]
        expected_contrib_acceleration(:, 2, 1) = [-2.0_real64, 2.0_real64]
        expected_acceleration_series(:, 1, 1, 1) = [0.0_real64, 0.0_real64, 1.0_real64, 1.0_real64]
        expected_acceleration_series(:, 2, 1, 1) = [0.0_real64, 0.0_real64, -1.0_real64, -1.0_real64]
        expected_acceleration_series(:, 1, 2, 1) = [0.0_real64, 0.0_real64, -1.0_real64, -1.0_real64]
        expected_acceleration_series(:, 2, 2, 1) = [0.0_real64, 0.0_real64, 1.0_real64, 1.0_real64]
    end subroutine fill_raw_case

    !> The data of `fill_raw_case`, plain tier.
    subroutine test_compute_velocity_acceleration_contributions_raw()
        real(real64) :: trajectories(2, 1, 4)
        real(real64) :: expected_contrib_velocity(2, 2, 1), expected_velocity_series(4, 2, 2, 1)
        real(real64) :: expected_contrib_acceleration(2, 2, 1), expected_acceleration_series(4, 2, 2, 1)

        call fill_raw_case(trajectories, expected_contrib_velocity, expected_velocity_series, &
                           expected_contrib_acceleration, expected_acceleration_series)
        call check_plain(2, 1, 4, trajectories, MODE_BASELINE_RAW, expected_contrib_velocity, expected_velocity_series, &
                         expected_contrib_acceleration, expected_acceleration_series, &
                         "test_compute_velocity_acceleration_contributions_raw")
    end subroutine test_compute_velocity_acceleration_contributions_raw

    !> The data of `fill_raw_case` through the `_expert` tier, with work arrays that start out
    !| full of a stale value (it reuses them for velocities and accelerations): the same result.
    subroutine test_compute_velocity_acceleration_contributions_expert()
        real(real64) :: trajectories(2, 1, 4)
        real(real64) :: expected_contrib_velocity(2, 2, 1), expected_velocity_series(4, 2, 2, 1)
        real(real64) :: expected_contrib_acceleration(2, 2, 1), expected_acceleration_series(4, 2, 2, 1)
        real(real64) :: contrib_velocity(2, 2, 1), velocity_series(4, 2, 2, 1)
        real(real64) :: contrib_acceleration(2, 2, 1), acceleration_series(4, 2, 2, 1)
        real(real64) :: tmp_factors(3, 2), tmp_dependent(3), tmp_contributions(3)
        integer(int32) :: ierr

        call fill_raw_case(trajectories, expected_contrib_velocity, expected_velocity_series, &
                           expected_contrib_acceleration, expected_acceleration_series)
        contrib_velocity = UNWRITTEN
        velocity_series = UNWRITTEN
        contrib_acceleration = UNWRITTEN
        acceleration_series = UNWRITTEN
        tmp_factors = 1.0e300_real64
        tmp_dependent = -1.0e300_real64
        tmp_contributions = 1.0e300_real64

        call compute_velocity_acceleration_contributions_expert(trajectories, 2, 1, 4, MODE_BASELINE_RAW, &
                                                                tmp_factors, tmp_dependent, tmp_contributions, &
                                                                contrib_velocity, velocity_series, &
                                                                contrib_acceleration, acceleration_series, ierr)
        call check_outputs(2, 1, 4, ierr, contrib_velocity, velocity_series, contrib_acceleration, &
                           acceleration_series, expected_contrib_velocity, expected_velocity_series, &
                           expected_contrib_acceleration, expected_acceleration_series, &
                           "test_compute_velocity_acceleration_contributions_expert")
    end subroutine test_compute_velocity_acceleration_contributions_expert

    !> MODE_BASELINE_MEAN on three time points, the fewest with an acceleration. x1 [0, 2, 3]:
    !| velocity [2, 1], mean 1.5, deviations [0.5, -0.5]; x2 [1, 1, 4]: velocity [0, 3], mean 1.5,
    !| deviations [-1.5, 1.5]. Velocity: (x1, x1) [0.25, 0.25] = 0.5, (x1, x2) and (x2, x1)
    !| [-0.75, -0.75] = -1.5, (x2, x2) [2.25, 2.25] = 4.5. Each acceleration is a single point
    !| (x1 -1, x2 3), which is its own mean, so every acceleration contribution is 0.
    subroutine test_compute_velocity_acceleration_contributions_mean()
        real(real64) :: trajectories(2, 1, 3)
        real(real64) :: expected_contrib_velocity(2, 2, 1), expected_velocity_series(3, 2, 2, 1)
        real(real64) :: expected_contrib_acceleration(2, 2, 1), expected_acceleration_series(3, 2, 2, 1)

        trajectories(1, 1, :) = [0.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2, 1, :) = [1.0_real64, 1.0_real64, 4.0_real64]
        expected_contrib_velocity(:, 1, 1) = [0.5_real64, -1.5_real64]
        expected_contrib_velocity(:, 2, 1) = [-1.5_real64, 4.5_real64]
        expected_velocity_series(:, 1, 1, 1) = [0.0_real64, 0.25_real64, 0.25_real64]
        expected_velocity_series(:, 2, 1, 1) = [0.0_real64, -0.75_real64, -0.75_real64]
        expected_velocity_series(:, 1, 2, 1) = [0.0_real64, -0.75_real64, -0.75_real64]
        expected_velocity_series(:, 2, 2, 1) = [0.0_real64, 2.25_real64, 2.25_real64]
        expected_contrib_acceleration = 0.0_real64
        expected_acceleration_series = 0.0_real64
        call check_plain(2, 1, 3, trajectories, MODE_BASELINE_MEAN, expected_contrib_velocity, &
                         expected_velocity_series, expected_contrib_acceleration, expected_acceleration_series, &
                         "test_compute_velocity_acceleration_contributions_mean")
    end subroutine test_compute_velocity_acceleration_contributions_mean

    !> MODE_BASELINE_MIN, one factor in two samples, each sample on its own:
    !| sample 1, [0, 1, 3, 6]: velocity [1, 2, 3], minimum 1, deviations [0, 1, 2], squared
    !|   [0, 1, 4] = 5; acceleration [1, 1], minimum 1, deviations [0, 0], so 0;
    !| sample 2, [5, 1, 3, -3]: velocity [-4, 2, -6], minimum -6, deviations [2, 8, 0], squared
    !|   [4, 64, 0] = 68; acceleration [6, -8], minimum -8, deviations [14, 0], squared
    !|   [196, 0] = 196.
    subroutine test_compute_velocity_acceleration_contributions_min()
        real(real64) :: trajectories(1, 2, 4)
        real(real64) :: expected_contrib_velocity(1, 1, 2), expected_velocity_series(4, 1, 1, 2)
        real(real64) :: expected_contrib_acceleration(1, 1, 2), expected_acceleration_series(4, 1, 1, 2)

        trajectories(1, 1, :) = [0.0_real64, 1.0_real64, 3.0_real64, 6.0_real64]
        trajectories(1, 2, :) = [5.0_real64, 1.0_real64, 3.0_real64, -3.0_real64]
        expected_contrib_velocity(1, 1, :) = [5.0_real64, 68.0_real64]
        expected_velocity_series(:, 1, 1, 1) = [0.0_real64, 0.0_real64, 1.0_real64, 4.0_real64]
        expected_velocity_series(:, 1, 1, 2) = [0.0_real64, 4.0_real64, 64.0_real64, 0.0_real64]
        expected_contrib_acceleration(1, 1, :) = [0.0_real64, 196.0_real64]
        expected_acceleration_series(:, 1, 1, 1) = [0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64]
        expected_acceleration_series(:, 1, 1, 2) = [0.0_real64, 0.0_real64, 196.0_real64, 0.0_real64]
        call check_plain(1, 2, 4, trajectories, MODE_BASELINE_MIN, expected_contrib_velocity, &
                         expected_velocity_series, expected_contrib_acceleration, expected_acceleration_series, &
                         "test_compute_velocity_acceleration_contributions_min")
    end subroutine test_compute_velocity_acceleration_contributions_min

    !> Two time points: one velocity per series and no acceleration. x1 [1, 4] has the velocity
    !| 3, x2 [2, 0] has -2, so, raw, (x1, x1) is 9, (x1, x2) and (x2, x1) -6, (x2, x2) 4, each at
    !| time point 2. Every acceleration output is 0, written, not left as it was.
    subroutine test_compute_velocity_acceleration_contributions_two_points()
        real(real64) :: trajectories(2, 1, 2)
        real(real64) :: expected_contrib_velocity(2, 2, 1), expected_velocity_series(2, 2, 2, 1)
        real(real64) :: expected_contrib_acceleration(2, 2, 1), expected_acceleration_series(2, 2, 2, 1)

        trajectories(1, 1, :) = [1.0_real64, 4.0_real64]
        trajectories(2, 1, :) = [2.0_real64, 0.0_real64]
        expected_contrib_velocity(:, 1, 1) = [9.0_real64, -6.0_real64]
        expected_contrib_velocity(:, 2, 1) = [-6.0_real64, 4.0_real64]
        expected_velocity_series(:, 1, 1, 1) = [0.0_real64, 9.0_real64]
        expected_velocity_series(:, 2, 1, 1) = [0.0_real64, -6.0_real64]
        expected_velocity_series(:, 1, 2, 1) = [0.0_real64, -6.0_real64]
        expected_velocity_series(:, 2, 2, 1) = [0.0_real64, 4.0_real64]
        expected_contrib_acceleration = 0.0_real64
        expected_acceleration_series = 0.0_real64
        call check_plain(2, 1, 2, trajectories, MODE_BASELINE_RAW, expected_contrib_velocity, &
                         expected_velocity_series, expected_contrib_acceleration, expected_acceleration_series, &
                         "test_compute_velocity_acceleration_contributions_two_points")
    end subroutine test_compute_velocity_acceleration_contributions_two_points

    !> One time point, the smallest valid size: no velocity and no acceleration, so every output
    !| is 0, written, in each mode.
    subroutine test_compute_velocity_acceleration_contributions_one_point()
        real(real64) :: trajectories(2, 2, 1), expected_contribs(2, 2, 2), expected_series(1, 2, 2, 2)
        integer(int32) :: modes(3), i_mode

        trajectories(:, 1, 1) = [3.0_real64, -1.0_real64]
        trajectories(:, 2, 1) = [2.0_real64, 5.0_real64]
        expected_contribs = 0.0_real64
        expected_series = 0.0_real64
        modes = [MODE_BASELINE_RAW, MODE_BASELINE_MIN, MODE_BASELINE_MEAN]
        do i_mode = 1, 3
            call check_plain(2, 2, 1, trajectories, modes(i_mode), expected_contribs, expected_series, &
                             expected_contribs, expected_series, &
                             "test_compute_velocity_acceleration_contributions_one_point")
        end do
    end subroutine test_compute_velocity_acceleration_contributions_one_point

    !> Each dimension on its own, the others valid: `n_factors` (argument 2), `n_samples` (3) and
    !| `n_timepoints` (4) are ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1.
    subroutine test_compute_velocity_acceleration_contributions_dimensions()
        integer(int32), parameter :: SIZE_VALUES(2) = [0, -1]
        integer(int32), parameter :: EXPECTED_CODES(2) = [ERR_EMPTY_INPUT, ERR_INVALID_INPUT]
        character(*), parameter :: DIM_NAMES(3) = ["n_factors   ", "n_samples   ", "n_timepoints"]
        real(real64) :: trajectories(2, 2, 3), contrib_velocity(2, 2, 2), velocity_series(3, 2, 2, 2)
        real(real64) :: contrib_acceleration(2, 2, 2), acceleration_series(3, 2, 2, 2)
        integer(int32) :: dims(3), i_dim, i_value, ierr

        trajectories = 1.0_real64
        do i_dim = 1, 3
            do i_value = 1, 2
                dims = [2, 2, 3]
                dims(i_dim) = SIZE_VALUES(i_value)
                call compute_velocity_acceleration_contributions(trajectories, dims(1), dims(2), dims(3), &
                                                                 MODE_BASELINE_RAW, contrib_velocity, &
                                                                 velocity_series, contrib_acceleration, &
                                                                 acceleration_series, ierr)
                call assert_err(ierr, EXPECTED_CODES(i_value), &
                                "test_compute_velocity_acceleration_contributions_dimensions: " &
                                //trim(DIM_NAMES(i_dim))//" = "//merge(" 0", "-1", i_value == 1), i_dim + 1)
            end do
        end do
    end subroutine test_compute_velocity_acceleration_contributions_dimensions

    !> `baseline_mode` (argument 5) is one of the three MODE_BASELINE_* values (1 to 3); 0 and 4,
    !| just outside, are ERR_INVALID_INPUT.
    subroutine test_compute_velocity_acceleration_contributions_bad_mode()
        real(real64) :: trajectories(2, 1, 3), contrib_velocity(2, 2, 1), velocity_series(3, 2, 2, 1)
        real(real64) :: contrib_acceleration(2, 2, 1), acceleration_series(3, 2, 2, 1)
        integer(int32) :: ierr

        trajectories = 1.0_real64
        call compute_velocity_acceleration_contributions(trajectories, 2, 1, 3, 0_int32, contrib_velocity, &
                                                         velocity_series, contrib_acceleration, &
                                                         acceleration_series, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_velocity_acceleration_contributions_bad_mode: 0", 5)
        call compute_velocity_acceleration_contributions(trajectories, 2, 1, 3, 4_int32, contrib_velocity, &
                                                         velocity_series, contrib_acceleration, &
                                                         acceleration_series, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_velocity_acceleration_contributions_bad_mode: 4", 5)
    end subroutine test_compute_velocity_acceleration_contributions_bad_mode

    !> NaN, +Inf and -Inf in `trajectories` (argument 1) are ERR_NAN_INF; each sits in the last
    !| element, so the check must cover the whole array.
    subroutine test_compute_velocity_acceleration_contributions_nan_inf()
        real(real64) :: trajectories(2, 2, 3), contrib_velocity(2, 2, 2), velocity_series(3, 2, 2, 2)
        real(real64) :: contrib_acceleration(2, 2, 2), acceleration_series(3, 2, 2, 2)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            trajectories = 1.0_real64
            trajectories(2, 2, 3) = non_finite_value(i_value)
            call compute_velocity_acceleration_contributions(trajectories, 2, 2, 3, MODE_BASELINE_RAW, &
                                                             contrib_velocity, velocity_series, &
                                                             contrib_acceleration, acceleration_series, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_velocity_acceleration_contributions_nan_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_compute_velocity_acceleration_contributions_nan_inf

end module mod_test_compute_velocity_acceleration_contributions
