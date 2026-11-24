!> The `compute_velocity_trajectories` cases: the first differences over time of every
!| (factor, sample) series of `trajectories(n_factors, n_samples, n_timepoints)`, returned
!| time-first as `velocity(n_timepoints - 1, n_factors, n_samples)`, the two smallest sizes, and
!| the dimension and NaN/Inf checks. Every difference is exact.
module mod_test_compute_velocity_trajectories
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_velocity_trajectories
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_velocity_trajectories

contains

    !> Every case of `compute_velocity_trajectories`.
    function get_all_tests_compute_velocity_trajectories() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(5))
        all_tests(1) = test_case("test_compute_velocity_trajectories_values", test_compute_velocity_trajectories_values)
        all_tests(2) = test_case("test_compute_velocity_trajectories_two_timepoints", &
                                 test_compute_velocity_trajectories_two_timepoints)
        all_tests(3) = test_case("test_compute_velocity_trajectories_single_timepoint", &
                                 test_compute_velocity_trajectories_single_timepoint)
        all_tests(4) = test_case("test_compute_velocity_trajectories_dimensions", &
                                 test_compute_velocity_trajectories_dimensions)
        all_tests(5) = test_case("test_compute_velocity_trajectories_rejects_nan_and_inf", &
                                 test_compute_velocity_trajectories_rejects_nan_and_inf)
    end function get_all_tests_compute_velocity_trajectories

    !> Two factors in three samples over four time points, each series its own shape:
    !| sample 1: factor 1 [0, 1, 3, 6] -> [1, 2, 3],      factor 2 [5, 4, 4, 2] -> [-1, 0, -2];
    !| sample 2: factor 1 [2, 2, 2, 2] -> [0, 0, 0],      factor 2 [-1, 1, -1, 1] -> [2, -2, 2];
    !| sample 3: factor 1 [10, 20, 40, 80] -> [10, 20, 40], factor 2 [3, 0, -3, -6] -> [-3, -3, -3].
    !| Every series differs, so a velocity in the wrong place of the time-first layout shows. The
    !| output buffer is one element longer than the result and starts at UNWRITTEN, so a write
    !| past the end shows too.
    subroutine test_compute_velocity_trajectories_values()
        integer(int32), parameter :: n_factors = 2, n_samples = 3, n_timepoints = 4
        integer(int32), parameter :: n_velocities = (n_timepoints - 1)*n_factors*n_samples
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints), velocity(n_velocities + 1)
        real(real64) :: expected(n_timepoints - 1, n_factors, n_samples)
        integer(int32) :: ierr

        trajectories(1, 1, :) = [0.0_real64, 1.0_real64, 3.0_real64, 6.0_real64]
        trajectories(2, 1, :) = [5.0_real64, 4.0_real64, 4.0_real64, 2.0_real64]
        trajectories(1, 2, :) = [2.0_real64, 2.0_real64, 2.0_real64, 2.0_real64]
        trajectories(2, 2, :) = [-1.0_real64, 1.0_real64, -1.0_real64, 1.0_real64]
        trajectories(1, 3, :) = [10.0_real64, 20.0_real64, 40.0_real64, 80.0_real64]
        trajectories(2, 3, :) = [3.0_real64, 0.0_real64, -3.0_real64, -6.0_real64]
        expected(:, 1, 1) = [1.0_real64, 2.0_real64, 3.0_real64]
        expected(:, 2, 1) = [-1.0_real64, 0.0_real64, -2.0_real64]
        expected(:, 1, 2) = [0.0_real64, 0.0_real64, 0.0_real64]
        expected(:, 2, 2) = [2.0_real64, -2.0_real64, 2.0_real64]
        expected(:, 1, 3) = [10.0_real64, 20.0_real64, 40.0_real64]
        expected(:, 2, 3) = [-3.0_real64, -3.0_real64, -3.0_real64]
        velocity = UNWRITTEN

        call compute_velocity_trajectories(trajectories, velocity, n_factors, n_samples, n_timepoints, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_velocity_trajectories_values: ierr")
        call assert_equal_array_real(velocity, expected, n_velocities, 0.0_real64, &
                                     "test_compute_velocity_trajectories_values: velocity", n_timepoints - 1)
        call assert_equal_real(velocity(n_velocities + 1), UNWRITTEN, 0.0_real64, &
                               "test_compute_velocity_trajectories_values: no write past the end")
    end subroutine test_compute_velocity_trajectories_values

    !> n_timepoints = 2, the smallest with a velocity: factor 1 [1, 4] -> [3], factor 2
    !| [7, 2] -> [-5], one sample.
    subroutine test_compute_velocity_trajectories_two_timepoints()
        real(real64) :: trajectories(2, 1, 2), velocity(1, 2, 1), expected(1, 2, 1)
        integer(int32) :: ierr

        trajectories(1, 1, :) = [1.0_real64, 4.0_real64]
        trajectories(2, 1, :) = [7.0_real64, 2.0_real64]
        expected(1, :, 1) = [3.0_real64, -5.0_real64]
        velocity = UNWRITTEN

        call compute_velocity_trajectories(trajectories, velocity, 2, 1, 2, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_velocity_trajectories_two_timepoints: ierr")
        call assert_equal_array_real(velocity, expected, 2, 0.0_real64, &
                                     "test_compute_velocity_trajectories_two_timepoints: velocity")
    end subroutine test_compute_velocity_trajectories_two_timepoints

    !> n_timepoints = 1, the smallest valid size: no step, so the velocity is empty, and the
    !| buffer behind it stays as it was.
    subroutine test_compute_velocity_trajectories_single_timepoint()
        real(real64) :: trajectories(2, 2, 1), velocity(1)
        integer(int32) :: ierr

        trajectories = 1.0_real64
        velocity = UNWRITTEN
        call compute_velocity_trajectories(trajectories, velocity, 2, 2, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_velocity_trajectories_single_timepoint: ierr")
        call assert_equal_real(velocity(1), UNWRITTEN, 0.0_real64, &
                               "test_compute_velocity_trajectories_single_timepoint: nothing written")
    end subroutine test_compute_velocity_trajectories_single_timepoint

    !> Each dimension on its own, the others valid: `n_factors` (argument 3), `n_samples` (4) and
    !| `n_timepoints` (5) are ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1.
    subroutine test_compute_velocity_trajectories_dimensions()
        integer(int32), parameter :: SIZE_VALUES(2) = [0, -1]
        integer(int32), parameter :: EXPECTED_CODES(2) = [ERR_EMPTY_INPUT, ERR_INVALID_INPUT]
        character(*), parameter :: DIM_NAMES(3) = ["n_factors   ", "n_samples   ", "n_timepoints"]
        real(real64) :: trajectories(2, 2, 3), velocity(2, 2, 2)
        integer(int32) :: dims(3), i_dim, i_value, ierr

        trajectories = 1.0_real64
        do i_dim = 1, 3
            do i_value = 1, 2
                dims = [2, 2, 3]
                dims(i_dim) = SIZE_VALUES(i_value)
                call compute_velocity_trajectories(trajectories, velocity, dims(1), dims(2), dims(3), ierr)
                call assert_err(ierr, EXPECTED_CODES(i_value), "test_compute_velocity_trajectories_dimensions: " &
                                //trim(DIM_NAMES(i_dim))//" = "//merge(" 0", "-1", i_value == 1), i_dim + 2)
            end do
        end do
    end subroutine test_compute_velocity_trajectories_dimensions

    !> NaN, +Inf and -Inf in `trajectories` (argument 1) are ERR_NAN_INF; each sits in the last
    !| element, so the check must cover the whole array.
    subroutine test_compute_velocity_trajectories_rejects_nan_and_inf()
        real(real64) :: trajectories(2, 2, 3), velocity(2, 2, 2)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            trajectories = 1.0_real64
            trajectories(2, 2, 3) = non_finite_value(i_value)
            call compute_velocity_trajectories(trajectories, velocity, 2, 2, 3, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_velocity_trajectories_rejects_nan_and_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_compute_velocity_trajectories_rejects_nan_and_inf

end module mod_test_compute_velocity_trajectories
