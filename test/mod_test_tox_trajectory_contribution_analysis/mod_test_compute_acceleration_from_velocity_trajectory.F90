!> The `compute_acceleration_from_velocity_trajectory` cases: the first differences of one
!| velocity series, which are the second differences of its positions, the three smallest sizes,
!| and the dimension and NaN/Inf checks. Every difference is exact.
!|
!| `n_timepoints` counts the positions, so `velocity` has n_timepoints - 1 elements and the
!| acceleration n_timepoints - 2.
module mod_test_compute_acceleration_from_velocity_trajectory
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_acceleration_from_velocity_trajectory
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_acceleration_from_velocity_trajectory

contains

    !> Every case of `compute_acceleration_from_velocity_trajectory`.
    function get_all_tests_compute_acceleration_from_velocity_trajectory() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(6))
        all_tests(1) = test_case("test_compute_acceleration_from_velocity_trajectory_values", &
                                 test_compute_acceleration_from_velocity_trajectory_values)
        all_tests(2) = test_case("test_compute_acceleration_from_velocity_trajectory_three", &
                                 test_compute_acceleration_from_velocity_trajectory_three)
        all_tests(3) = test_case("test_compute_acceleration_from_velocity_trajectory_two", &
                                 test_compute_acceleration_from_velocity_trajectory_two)
        all_tests(4) = test_case("test_compute_acceleration_from_velocity_trajectory_one", &
                                 test_compute_acceleration_from_velocity_trajectory_one)
        all_tests(5) = test_case("test_compute_acceleration_from_velocity_trajectory_nan_inf", &
                                 test_compute_acceleration_from_velocity_trajectory_nan_inf)
        all_tests(6) = test_case("test_compute_acceleration_from_velocity_trajectory_dimensions", &
                                 test_compute_acceleration_from_velocity_trajectory_dimensions)
    end function get_all_tests_compute_acceleration_from_velocity_trajectory

    !> Calls the procedure for `n_timepoints` positions and checks `ierr` OK and the
    !| n_timepoints - 2 accelerations exactly. The acceleration starts at UNWRITTEN, one element
    !| longer than the result, so a step the procedure skips and a write past the end both show.
    subroutine check_acceleration(n_timepoints, velocity, expected, msg)
        integer(int32), intent(in) :: n_timepoints
        real(real64), intent(in) :: velocity(n_timepoints - 1), expected(n_timepoints - 2)
        character(*), intent(in) :: msg
        real(real64) :: acceleration(n_timepoints - 1)
        integer(int32) :: ierr

        acceleration = UNWRITTEN
        call compute_acceleration_from_velocity_trajectory(velocity, acceleration, n_timepoints, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(acceleration, expected, n_timepoints - 2, 0.0_real64, msg//": acceleration")
        call assert_equal_real(acceleration(n_timepoints - 1), UNWRITTEN, 0.0_real64, msg//": no write past the end")
    end subroutine check_acceleration

    !> Kept from the old suite: the velocity [1, 2, 3, 4] rises by [1, 1, 1]. A velocity that
    !| changes direction, [1, 3, 0, -2], gives [2, -3, -2]. And the second difference of
    !| positions, as issue #107 defines acceleration, a(t) = p(t) - 2 p(t-1) + p(t-2): the
    !| positions [3, 1, 4, 1, 5] give 4 - 2 + 3 = 5, 1 - 8 + 1 = -6 and 5 - 2 + 4 = 7, from
    !| their velocity [-2, 3, -3, 4].
    subroutine test_compute_acceleration_from_velocity_trajectory_values()
        real(real64) :: velocity(4), expected(3)

        velocity = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        expected = [1.0_real64, 1.0_real64, 1.0_real64]
        call check_acceleration(5, velocity, expected, &
                                "test_compute_acceleration_from_velocity_trajectory_values: steady rise")

        velocity = [1.0_real64, 3.0_real64, 0.0_real64, -2.0_real64]
        expected = [2.0_real64, -3.0_real64, -2.0_real64]
        call check_acceleration(5, velocity, expected, &
                                "test_compute_acceleration_from_velocity_trajectory_values: turning")

        velocity = [-2.0_real64, 3.0_real64, -3.0_real64, 4.0_real64]
        expected = [5.0_real64, -6.0_real64, 7.0_real64]
        call check_acceleration(5, velocity, expected, &
                                "test_compute_acceleration_from_velocity_trajectory_values: second difference")
    end subroutine test_compute_acceleration_from_velocity_trajectory_values

    !> Three positions, the fewest with an acceleration: the velocity [2, -1] gives [-3].
    subroutine test_compute_acceleration_from_velocity_trajectory_three()
        real(real64) :: velocity(2), expected(1)

        velocity = [2.0_real64, -1.0_real64]
        expected = -3.0_real64
        call check_acceleration(3, velocity, expected, "test_compute_acceleration_from_velocity_trajectory_three")
    end subroutine test_compute_acceleration_from_velocity_trajectory_three

    !> Two positions: one velocity [6] and no acceleration, which is fine and writes nothing.
    subroutine test_compute_acceleration_from_velocity_trajectory_two()
        real(real64) :: velocity(1), expected(0)

        velocity = 6.0_real64
        call check_acceleration(2, velocity, expected, "test_compute_acceleration_from_velocity_trajectory_two")
    end subroutine test_compute_acceleration_from_velocity_trajectory_two

    !> One position: no velocity and no acceleration, which is fine and writes nothing.
    subroutine test_compute_acceleration_from_velocity_trajectory_one()
        real(real64) :: acceleration(1)
        real(real64) :: velocity(0)
        integer(int32) :: ierr

        acceleration = UNWRITTEN
        call compute_acceleration_from_velocity_trajectory(velocity, acceleration, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_acceleration_from_velocity_trajectory_one: ierr")
        call assert_equal_real(acceleration(1), UNWRITTEN, 0.0_real64, &
                               "test_compute_acceleration_from_velocity_trajectory_one: nothing written")
    end subroutine test_compute_acceleration_from_velocity_trajectory_one

    !> `n_timepoints` (argument 3) has the minimum 1, as for the velocity procedures (FES): 0 and
    !| -1 are ERR_INVALID_INPUT, although both size an empty velocity and acceleration; 1 is in
    !| `test_compute_acceleration_from_velocity_trajectory_one`.
    subroutine test_compute_acceleration_from_velocity_trajectory_dimensions()
        real(real64) :: velocity(1), acceleration(1)
        integer(int32) :: ierr

        velocity = 1.0_real64
        call compute_acceleration_from_velocity_trajectory(velocity, acceleration, 0, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_acceleration_from_velocity_trajectory_dimensions: 0", 3)
        call compute_acceleration_from_velocity_trajectory(velocity, acceleration, -1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_acceleration_from_velocity_trajectory_dimensions: -1", 3)
    end subroutine test_compute_acceleration_from_velocity_trajectory_dimensions

    !> NaN, +Inf and -Inf in `velocity` (argument 1) are ERR_NAN_INF; each sits in the last
    !| element, so the check must cover every one.
    subroutine test_compute_acceleration_from_velocity_trajectory_nan_inf()
        real(real64) :: velocity(3), acceleration(2)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            velocity = [1.0_real64, 2.0_real64, 3.0_real64]
            velocity(3) = non_finite_value(i_value)
            call compute_acceleration_from_velocity_trajectory(velocity, acceleration, 4, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_acceleration_from_velocity_trajectory_nan_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_compute_acceleration_from_velocity_trajectory_nan_inf

end module mod_test_compute_acceleration_from_velocity_trajectory
