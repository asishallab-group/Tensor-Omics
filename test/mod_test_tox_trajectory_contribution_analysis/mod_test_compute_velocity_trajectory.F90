!> The `compute_velocity_trajectory` cases: the first differences `trajectory(t + 1) -
!| trajectory(t)` of one series, rising and falling, the two smallest sizes, extreme magnitudes,
!| and the dimension and NaN/Inf checks. Every difference is exact.
module mod_test_compute_velocity_trajectory
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_velocity_trajectory
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_velocity_trajectory

contains

    !> Every case of `compute_velocity_trajectory`.
    function get_all_tests_compute_velocity_trajectory() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(6))
        all_tests(1) = test_case("test_compute_velocity_trajectory_values", test_compute_velocity_trajectory_values)
        all_tests(2) = test_case("test_compute_velocity_trajectory_two_timepoints", &
                                 test_compute_velocity_trajectory_two_timepoints)
        all_tests(3) = test_case("test_compute_velocity_trajectory_single_timepoint", &
                                 test_compute_velocity_trajectory_single_timepoint)
        all_tests(4) = test_case("test_compute_velocity_trajectory_extreme_magnitudes", &
                                 test_compute_velocity_trajectory_extreme_magnitudes)
        all_tests(5) = test_case("test_compute_velocity_trajectory_dimensions", &
                                 test_compute_velocity_trajectory_dimensions)
        all_tests(6) = test_case("test_compute_velocity_trajectory_rejects_nan_and_inf", &
                                 test_compute_velocity_trajectory_rejects_nan_and_inf)
    end function get_all_tests_compute_velocity_trajectory

    !> Calls the procedure on `n_timepoints` points and checks `ierr` OK and the
    !| n_timepoints - 1 velocities exactly. The velocity starts at UNWRITTEN, one element longer
    !| than the result, so a step the procedure skips and a write past the end both show.
    subroutine check_velocity(n_timepoints, trajectory, expected, msg)
        integer(int32), intent(in) :: n_timepoints
        real(real64), intent(in) :: trajectory(n_timepoints), expected(n_timepoints - 1)
        character(*), intent(in) :: msg
        real(real64) :: velocity(n_timepoints)
        integer(int32) :: ierr

        velocity = UNWRITTEN
        call compute_velocity_trajectory(trajectory, velocity, n_timepoints, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(velocity, expected, n_timepoints - 1, 0.0_real64, msg//": velocity")
        call assert_equal_real(velocity(n_timepoints), UNWRITTEN, 0.0_real64, msg//": no write past the end")
    end subroutine check_velocity

    !> Kept from the old suite: [1, 2, 4, 7, 11] rises by [1, 2, 3, 4]. And a series that rises
    !| and falls, [3, 1, 4, 1, 5]: [-2, 3, -3, 4].
    subroutine test_compute_velocity_trajectory_values()
        real(real64) :: trajectory(5), expected(4)

        trajectory = [1.0_real64, 2.0_real64, 4.0_real64, 7.0_real64, 11.0_real64]
        expected = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        call check_velocity(5, trajectory, expected, "test_compute_velocity_trajectory_values: rising")

        trajectory = [3.0_real64, 1.0_real64, 4.0_real64, 1.0_real64, 5.0_real64]
        expected = [-2.0_real64, 3.0_real64, -3.0_real64, 4.0_real64]
        call check_velocity(5, trajectory, expected, "test_compute_velocity_trajectory_values: rising and falling")
    end subroutine test_compute_velocity_trajectory_values

    !> n_timepoints = 2, the smallest with a velocity: [5, -2] gives [-7].
    subroutine test_compute_velocity_trajectory_two_timepoints()
        real(real64) :: trajectory(2), expected(1)

        trajectory = [5.0_real64, -2.0_real64]
        expected = -7.0_real64
        call check_velocity(2, trajectory, expected, "test_compute_velocity_trajectory_two_timepoints")
    end subroutine test_compute_velocity_trajectory_two_timepoints

    !> n_timepoints = 1, the smallest valid size: no step, so the velocity is empty and nothing
    !| is written.
    subroutine test_compute_velocity_trajectory_single_timepoint()
        real(real64) :: trajectory(1), expected(0)

        trajectory = 4.0_real64
        call check_velocity(1, trajectory, expected, "test_compute_velocity_trajectory_single_timepoint")
    end subroutine test_compute_velocity_trajectory_single_timepoint

    !> The largest finite magnitudes, where the difference stays finite: [huge, huge, 0, -huge]
    !| steps by [0, -huge, -huge].
    subroutine test_compute_velocity_trajectory_extreme_magnitudes()
        real(real64) :: trajectory(4), expected(3)

        trajectory = [huge(1.0_real64), huge(1.0_real64), 0.0_real64, -huge(1.0_real64)]
        expected = [0.0_real64, -huge(1.0_real64), -huge(1.0_real64)]
        call check_velocity(4, trajectory, expected, "test_compute_velocity_trajectory_extreme_magnitudes")
    end subroutine test_compute_velocity_trajectory_extreme_magnitudes

    !> `n_timepoints` (argument 3) is ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1; 1 is in
    !| `test_compute_velocity_trajectory_single_timepoint`.
    subroutine test_compute_velocity_trajectory_dimensions()
        real(real64) :: trajectory(2), velocity(1)
        integer(int32) :: ierr

        trajectory = [1.0_real64, 2.0_real64]
        call compute_velocity_trajectory(trajectory, velocity, 0, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_velocity_trajectory_dimensions: 0", 3)
        call compute_velocity_trajectory(trajectory, velocity, -1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_velocity_trajectory_dimensions: -1", 3)
    end subroutine test_compute_velocity_trajectory_dimensions

    !> NaN, +Inf and -Inf in `trajectory` (argument 1) are ERR_NAN_INF; each sits in the last
    !| point, so the check must cover every point.
    subroutine test_compute_velocity_trajectory_rejects_nan_and_inf()
        real(real64) :: trajectory(3), velocity(2)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            trajectory = [1.0_real64, 2.0_real64, 3.0_real64]
            trajectory(3) = non_finite_value(i_value)
            call compute_velocity_trajectory(trajectory, velocity, 3, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_velocity_trajectory_rejects_nan_and_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_compute_velocity_trajectory_rejects_nan_and_inf

end module mod_test_compute_velocity_trajectory
