!> The `compute_acceleration_from_velocity` cases: the first differences over time of every
!| (factor, sample) velocity series of `velocity(n_timepoints - 1, n_factors, n_samples)`,
!| returned as `acceleration(n_timepoints - 2, n_factors, n_samples)`, the three smallest sizes,
!| and the dimension and NaN/Inf checks. Every difference is exact.
module mod_test_compute_acceleration_from_velocity
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_acceleration_from_velocity
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_acceleration_from_velocity

contains

    !> Every case of `compute_acceleration_from_velocity`.
    function get_all_tests_compute_acceleration_from_velocity() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(6))
        all_tests(1) = test_case("test_compute_acceleration_from_velocity_values", &
                                 test_compute_acceleration_from_velocity_values)
        all_tests(2) = test_case("test_compute_acceleration_from_velocity_three_timepoints", &
                                 test_compute_acceleration_from_velocity_three_timepoints)
        all_tests(3) = test_case("test_compute_acceleration_from_velocity_two_timepoints", &
                                 test_compute_acceleration_from_velocity_two_timepoints)
        all_tests(4) = test_case("test_compute_acceleration_from_velocity_one_timepoint", &
                                 test_compute_acceleration_from_velocity_one_timepoint)
        all_tests(5) = test_case("test_compute_acceleration_from_velocity_dimensions", &
                                 test_compute_acceleration_from_velocity_dimensions)
        all_tests(6) = test_case("test_compute_acceleration_from_velocity_rejects_nan_inf", &
                                 test_compute_acceleration_from_velocity_rejects_nan_inf)
    end function get_all_tests_compute_acceleration_from_velocity

    !> Four time points, so three velocities and two accelerations per series, two factors in
    !| two samples:
    !| sample 1: factor 1 [1, 3, 0] -> [2, -3],   factor 2 [0, 0, 0] -> [0, 0];
    !| sample 2: factor 1 [-2, 2, -2] -> [4, -4], factor 2 [5, 1, 6] -> [-4, 5].
    !| Every series differs, so an acceleration in the wrong place shows. The output buffer is
    !| one element longer than the result and starts at UNWRITTEN, so a write past the end shows.
    subroutine test_compute_acceleration_from_velocity_values()
        integer(int32), parameter :: n_factors = 2, n_samples = 2, n_timepoints = 4
        integer(int32), parameter :: n_accelerations = (n_timepoints - 2)*n_factors*n_samples
        real(real64) :: velocity(n_timepoints - 1, n_factors, n_samples), acceleration(n_accelerations + 1)
        real(real64) :: expected(n_timepoints - 2, n_factors, n_samples)
        integer(int32) :: ierr

        velocity(:, 1, 1) = [1.0_real64, 3.0_real64, 0.0_real64]
        velocity(:, 2, 1) = [0.0_real64, 0.0_real64, 0.0_real64]
        velocity(:, 1, 2) = [-2.0_real64, 2.0_real64, -2.0_real64]
        velocity(:, 2, 2) = [5.0_real64, 1.0_real64, 6.0_real64]
        expected(:, 1, 1) = [2.0_real64, -3.0_real64]
        expected(:, 2, 1) = [0.0_real64, 0.0_real64]
        expected(:, 1, 2) = [4.0_real64, -4.0_real64]
        expected(:, 2, 2) = [-4.0_real64, 5.0_real64]
        acceleration = UNWRITTEN

        call compute_acceleration_from_velocity(velocity, acceleration, n_factors, n_samples, n_timepoints, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_acceleration_from_velocity_values: ierr")
        call assert_equal_array_real(acceleration, expected, n_accelerations, 0.0_real64, &
                                     "test_compute_acceleration_from_velocity_values: acceleration", n_timepoints - 2)
        call assert_equal_real(acceleration(n_accelerations + 1), UNWRITTEN, 0.0_real64, &
                               "test_compute_acceleration_from_velocity_values: no write past the end")
    end subroutine test_compute_acceleration_from_velocity_values

    !> Three time points, the fewest with an acceleration: one factor in two samples, velocities
    !| [2, -1] -> [-3] and [0, 4] -> [4].
    subroutine test_compute_acceleration_from_velocity_three_timepoints()
        real(real64) :: velocity(2, 1, 2), acceleration(1, 1, 2), expected(1, 1, 2)
        integer(int32) :: ierr

        velocity(:, 1, 1) = [2.0_real64, -1.0_real64]
        velocity(:, 1, 2) = [0.0_real64, 4.0_real64]
        expected(1, 1, :) = [-3.0_real64, 4.0_real64]
        acceleration = UNWRITTEN

        call compute_acceleration_from_velocity(velocity, acceleration, 1, 2, 3, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_acceleration_from_velocity_three_timepoints: ierr")
        call assert_equal_array_real(acceleration, expected, 2, 0.0_real64, &
                                     "test_compute_acceleration_from_velocity_three_timepoints: acceleration")
    end subroutine test_compute_acceleration_from_velocity_three_timepoints

    !> Two time points: one velocity per series and no acceleration, which is fine; the buffer
    !| behind the empty result stays as it was.
    subroutine test_compute_acceleration_from_velocity_two_timepoints()
        real(real64) :: velocity(1, 2, 1), acceleration(1)
        integer(int32) :: ierr

        velocity(1, :, 1) = [3.0_real64, -1.0_real64]
        acceleration = UNWRITTEN
        call compute_acceleration_from_velocity(velocity, acceleration, 2, 1, 2, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_acceleration_from_velocity_two_timepoints: ierr")
        call assert_equal_real(acceleration(1), UNWRITTEN, 0.0_real64, &
                               "test_compute_acceleration_from_velocity_two_timepoints: nothing written")
    end subroutine test_compute_acceleration_from_velocity_two_timepoints

    !> One time point: no velocity and no acceleration, which is fine and writes nothing.
    subroutine test_compute_acceleration_from_velocity_one_timepoint()
        real(real64) :: velocity(1), acceleration(1)
        integer(int32) :: ierr

        velocity = UNWRITTEN
        acceleration = UNWRITTEN
        call compute_acceleration_from_velocity(velocity, acceleration, 2, 1, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_acceleration_from_velocity_one_timepoint: ierr")
        call assert_equal_real(acceleration(1), UNWRITTEN, 0.0_real64, &
                               "test_compute_acceleration_from_velocity_one_timepoint: nothing written")
    end subroutine test_compute_acceleration_from_velocity_one_timepoint

    !> `n_factors` (argument 3) and `n_samples` (4) are ERR_EMPTY_INPUT at 0 and
    !| ERR_INVALID_INPUT at -1, each with the other valid.
    subroutine test_compute_acceleration_from_velocity_dimensions()
        real(real64) :: velocity(2, 2, 2), acceleration(1, 2, 2)
        integer(int32) :: ierr

        velocity = 1.0_real64
        call compute_acceleration_from_velocity(velocity, acceleration, 0, 2, 3, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_acceleration_from_velocity_dimensions: n_factors = 0", 3)
        call compute_acceleration_from_velocity(velocity, acceleration, -1, 2, 3, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_acceleration_from_velocity_dimensions: n_factors = -1", 3)
        call compute_acceleration_from_velocity(velocity, acceleration, 2, 0, 3, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_acceleration_from_velocity_dimensions: n_samples = 0", 4)
        call compute_acceleration_from_velocity(velocity, acceleration, 2, -1, 3, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_acceleration_from_velocity_dimensions: n_samples = -1", 4)
    end subroutine test_compute_acceleration_from_velocity_dimensions

    !> NaN, +Inf and -Inf in `velocity` (argument 1) are ERR_NAN_INF; each sits in the last
    !| element, so the check must cover the whole array.
    subroutine test_compute_acceleration_from_velocity_rejects_nan_inf()
        real(real64) :: velocity(2, 2, 2), acceleration(1, 2, 2)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            velocity = 1.0_real64
            velocity(2, 2, 2) = non_finite_value(i_value)
            call compute_acceleration_from_velocity(velocity, acceleration, 2, 2, 3, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_acceleration_from_velocity_rejects_nan_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_compute_acceleration_from_velocity_rejects_nan_inf

end module mod_test_compute_acceleration_from_velocity
