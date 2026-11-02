!> Unit test suite for trajectory_normalization routines.
module mod_test_trajectory_normalization
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_normalization
    use tox_errors
    use test_suite, only: test_case, TEST_N_OMP_THREADS
    implicit none

    real(real64), parameter :: TOL = 1.0e-10_real64

contains

    function get_all_tests_trajectory_normalization() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(4))
        all_tests(1) = test_case("test_normalize_variable_timeseries", test_normalize_variable_timeseries)
        all_tests(2) = test_case("test_normalize_single_trajectory", test_normalize_single_trajectory)
        all_tests(3) = test_case("test_normalize_all_trajectories", test_normalize_all_trajectories)
        all_tests(4) = test_case("test_normalize_all_trajectories_many_series", test_normalize_all_trajectories_many_series)
    end function get_all_tests_trajectory_normalization

    !> Test the normalization of variable timeseries.
    subroutine test_normalize_variable_timeseries()
        integer(int32), parameter :: n_points = 5
        real(real64) :: v(n_points), v_norm(n_points), expected(n_points)
        integer(int32) :: ierr, status

        v = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        expected = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]

        call normalize_variable_timeseries(v, v_norm, n_points, status, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_variable_timeseries: ierr")
        call assert_equal_array_real(v_norm, expected, n_points, TOL, "test_normalize_variable_timeseries: values")
    end subroutine test_normalize_variable_timeseries

    !> Test the normalization of a single trajectory.
    subroutine test_normalize_single_trajectory()
        integer(int32), parameter :: n_factors = 2, n_timepoints = 4
        real(real64) :: trajectory(n_timepoints, n_factors), trajectory_norm(n_timepoints, n_factors)
        real(real64) :: expected(n_timepoints)
        integer(int32) :: ierr, status(n_factors), i

        trajectory(:, 1) = [10.0_real64, 11.0_real64, 12.0_real64, 13.0_real64]
        trajectory(:, 2) = [20.0_real64, 21.0_real64, 22.0_real64, 23.0_real64]
        expected = [0.0_real64, 1.0_real64/3.0_real64, 2.0_real64/3.0_real64, 1.0_real64]

        call normalize_single_trajectory(trajectory, trajectory_norm, n_factors, n_timepoints, status, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_single_trajectory: ierr")
        do i = 1, n_factors
            call assert_equal_array_real(trajectory_norm(:, i), expected, n_timepoints, TOL, "test_normalize_single_trajectory: factor in iteration " // i)
            call assert_equal_int(status(i), ERR_OK, "test_normalize_single_trajectory: status in iteration " // i)
        end do
    end subroutine test_normalize_single_trajectory

    !> Test the normalization of all trajectories in a 3D array.
    subroutine test_normalize_all_trajectories()
        integer(int32), parameter :: n_factors = 2, n_samples = 2, n_timepoints = 4
        real(real64) :: trajectories(n_factors, n_samples, n_timepoints)
        real(real64) :: trajectories_norm(n_factors, n_samples, n_timepoints)
        integer(int32) :: ierr, status(n_factors, n_samples), i_factor, i_sample

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        trajectories(2, 1, :) = [10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64]
        trajectories(1, 2, :) = [5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64]
        trajectories(2, 2, :) = [2.0_real64, 4.0_real64, 6.0_real64, 8.0_real64]

        call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, status, &
                                        TEST_N_OMP_THREADS, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_all_trajectories: ierr")
        do i_sample = 1, n_samples
            do i_factor = 1, n_factors
                call assert_equal_real(minval(trajectories_norm(i_factor, i_sample, :)), 0.0_real64, TOL, "test_normalize_all_trajectories: min in iteration i_sample=" // i_sample // ", i_factor=" // i_factor)
                call assert_equal_real(maxval(trajectories_norm(i_factor, i_sample, :)), 1.0_real64, TOL, "test_normalize_all_trajectories: max in iteration i_sample=" // i_sample // ", i_factor=" // i_factor)
                call assert_equal_int(status(i_factor, i_sample), ERR_OK, "test_normalize_all_trajectories: status in iteration i_sample=" // i_sample // ", i_factor=" // i_factor)
            end do
        end do
    end subroutine test_normalize_all_trajectories

    !> Enough series for the parallel loop to run in an `--openmp` build: every series but those of
    !> the first factor is a ramp, whatever its slope normalized to the same values, and the first
    !> factor is constant, which its status says in every sample.
    subroutine test_normalize_all_trajectories_many_series()
        integer(int32), parameter :: n_factors = 80, n_samples = 100, n_timepoints = 4
        real(real64), allocatable :: trajectories(:, :, :), trajectories_norm(:, :, :)
        integer(int32), allocatable :: status(:, :)
        integer(int32) :: ierr, i_factor, i_sample, i_timepoint
        real(real64) :: expected

        allocate (trajectories(n_factors, n_samples, n_timepoints), trajectories_norm(n_factors, n_samples, n_timepoints))
        allocate (status(n_factors, n_samples))
        do i_timepoint = 1, n_timepoints
            do i_sample = 1, n_samples
                trajectories(1, i_sample, i_timepoint) = 7.0_real64
                do i_factor = 2, n_factors
                    trajectories(i_factor, i_sample, i_timepoint) = real((i_factor + i_sample)*i_timepoint, real64)
                end do
            end do
        end do

        call normalize_all_trajectories(trajectories, trajectories_norm, n_factors, n_samples, n_timepoints, status, &
                                        TEST_N_OMP_THREADS, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_all_trajectories_many_series: ierr")
        do i_sample = 1, n_samples
            call assert_equal_int(get_err_code(status(1, i_sample)), ERR_DIVISION_BY_ZERO, "test_normalize_all_trajectories_many_series: status of a constant series")
            do i_factor = 2, n_factors
                call assert_equal_int(get_err_code(status(i_factor, i_sample)), ERR_OK, "test_normalize_all_trajectories_many_series: status of a ramp")
            end do
            do i_timepoint = 1, n_timepoints
                expected = real(i_timepoint - 1, real64)/real(n_timepoints - 1, real64)
                call assert_equal_real(trajectories_norm(1, i_sample, i_timepoint), 0.0_real64, TOL, "test_normalize_all_trajectories_many_series: a constant series")
                do i_factor = 2, n_factors
                    call assert_equal_real(trajectories_norm(i_factor, i_sample, i_timepoint), expected, TOL, "test_normalize_all_trajectories_many_series: a ramp")
                end do
            end do
        end do
    end subroutine test_normalize_all_trajectories_many_series

end module mod_test_trajectory_normalization
