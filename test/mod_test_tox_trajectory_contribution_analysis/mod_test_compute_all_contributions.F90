!> The `compute_all_contributions` cases: the contributions of every selected (factor, dependent)
!| pair in every sample, laid out as `(n_timepoints, n_selected_factors, n_selected_dependents,
!| n_samples)`, in the order the indices are given, repeated indices, each baseline mode, a single
!| time point, the `_expert` tier, and the dimension, index, mode and NaN/Inf checks.
!|
!| The values are small integers and the one mean divides by 4, so every expected value is exact.
module mod_test_compute_all_contributions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_all_contributions, compute_all_contributions_expert
    use tox_trajectory_contribution_analysis_impl, only: MODE_BASELINE_RAW, MODE_BASELINE_MIN, MODE_BASELINE_MEAN
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_all_contributions

    ! The selection of `test_compute_all_contributions_values`: three factors in two samples over
    ! three time points, two of them selected as factors and two as dependents.
    integer(int32), parameter :: N_FACTORS = 3, N_SAMPLES = 2, N_TIMEPOINTS = 3
    integer(int32), parameter :: N_SELECTED_FACTORS = 2, N_SELECTED_DEPENDENTS = 2

contains

    !> Every case of `compute_all_contributions`.
    function get_all_tests_compute_all_contributions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(10))
        all_tests(1) = test_case("test_compute_all_contributions_values", test_compute_all_contributions_values)
        all_tests(2) = test_case("test_compute_all_contributions_expert", test_compute_all_contributions_expert)
        all_tests(3) = test_case("test_compute_all_contributions_mean", test_compute_all_contributions_mean)
        all_tests(4) = test_case("test_compute_all_contributions_min", test_compute_all_contributions_min)
        all_tests(5) = test_case("test_compute_all_contributions_repeated_indices", &
                                 test_compute_all_contributions_repeated_indices)
        all_tests(6) = test_case("test_compute_all_contributions_single_timepoint", &
                                 test_compute_all_contributions_single_timepoint)
        all_tests(7) = test_case("test_compute_all_contributions_dimensions", &
                                 test_compute_all_contributions_dimensions)
        all_tests(8) = test_case("test_compute_all_contributions_index_bounds", &
                                 test_compute_all_contributions_index_bounds)
        all_tests(9) = test_case("test_compute_all_contributions_invalid_mode", &
                                 test_compute_all_contributions_invalid_mode)
        all_tests(10) = test_case("test_compute_all_contributions_rejects_nan_and_inf", &
                                  test_compute_all_contributions_rejects_nan_and_inf)
    end function get_all_tests_compute_all_contributions

    !> Three factors x1, x2, x3 in two samples over three time points:
    !| sample 1: x1 [1, 2, 3], x2 [2, 0, -1], x3 [0, 1, 2];
    !| sample 2: x1 [-1, 0, 2], x2 [1, 3, 1], x3 [2, -2, 1].
    subroutine fill_three_factors(trajectories)
        real(real64), intent(out) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2, 1, :) = [2.0_real64, 0.0_real64, -1.0_real64]
        trajectories(3, 1, :) = [0.0_real64, 1.0_real64, 2.0_real64]
        trajectories(1, 2, :) = [-1.0_real64, 0.0_real64, 2.0_real64]
        trajectories(2, 2, :) = [1.0_real64, 3.0_real64, 1.0_real64]
        trajectories(3, 2, :) = [2.0_real64, -2.0_real64, 1.0_real64]
    end subroutine fill_three_factors

    !> The raw contributions of factors [3, 1] against dependents [2, 3] in the data of
    !| `fill_three_factors`, each point the product of the two values:
    !| sample 1: (x3, x2) [0, 0, -2] total -2; (x1, x2) [2, 0, -3] total -1;
    !|           (x3, x3) [0, 1, 4] total 5;   (x1, x3) [0, 2, 6] total 8;
    !| sample 2: (x3, x2) [2, -6, 1] total -3; (x1, x2) [-1, 0, 2] total 1;
    !|           (x3, x3) [4, 4, 1] total 9;   (x1, x3) [-2, 0, 2] total 0.
    !| x3 is a factor and a dependent at once, and every block differs from the others, so a
    !| block in the wrong place, or an index taken in another order, shows.
    subroutine fill_expected_three_factors(expected_local, expected_total)
        real(real64), intent(out) :: expected_local(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64), intent(out) :: expected_total(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)

        expected_local(:, 1, 1, 1) = [0.0_real64, 0.0_real64, -2.0_real64]
        expected_local(:, 2, 1, 1) = [2.0_real64, 0.0_real64, -3.0_real64]
        expected_local(:, 1, 2, 1) = [0.0_real64, 1.0_real64, 4.0_real64]
        expected_local(:, 2, 2, 1) = [0.0_real64, 2.0_real64, 6.0_real64]
        expected_local(:, 1, 1, 2) = [2.0_real64, -6.0_real64, 1.0_real64]
        expected_local(:, 2, 1, 2) = [-1.0_real64, 0.0_real64, 2.0_real64]
        expected_local(:, 1, 2, 2) = [4.0_real64, 4.0_real64, 1.0_real64]
        expected_local(:, 2, 2, 2) = [-2.0_real64, 0.0_real64, 2.0_real64]
        expected_total(:, 1, 1) = [-2.0_real64, -1.0_real64]
        expected_total(:, 2, 1) = [5.0_real64, 8.0_real64]
        expected_total(:, 1, 2) = [-3.0_real64, 1.0_real64]
        expected_total(:, 2, 2) = [9.0_real64, 0.0_real64]
    end subroutine fill_expected_three_factors

    !> Factors [3, 1] against dependents [2, 3] in two samples, raw: see
    !| `fill_expected_three_factors`. Indices 1 and n_factors are both selected, so both ends of
    !| their range are accepted.
    subroutine test_compute_all_contributions_values()
        real(real64) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)
        real(real64) :: local_contributions(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: total_contributions(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: expected_local(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: expected_total(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        integer(int32) :: factor_indices(N_SELECTED_FACTORS), dependent_indices(N_SELECTED_DEPENDENTS), ierr

        call fill_three_factors(trajectories)
        call fill_expected_three_factors(expected_local, expected_total)
        factor_indices = [3, 1]
        dependent_indices = [2, 3]
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN

        call compute_all_contributions(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, &
                                       N_SELECTED_FACTORS, dependent_indices, N_SELECTED_DEPENDENTS, &
                                       MODE_BASELINE_RAW, local_contributions, total_contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_all_contributions_values: ierr")
        call assert_equal_array_real(local_contributions, expected_local, size(expected_local, kind=int32), &
                                     0.0_real64, "test_compute_all_contributions_values: local_contributions", &
                                     N_TIMEPOINTS)
        call assert_equal_array_real(total_contributions, expected_total, size(expected_total, kind=int32), &
                                     0.0_real64, "test_compute_all_contributions_values: total_contributions", &
                                     N_SELECTED_FACTORS)
    end subroutine test_compute_all_contributions_values

    !> The `_expert` tier on the data of `test_compute_all_contributions_values`, with work
    !| arrays that start out full of a stale value: it must give the same contributions, so it
    !| may not read a work array before writing it.
    subroutine test_compute_all_contributions_expert()
        real(real64) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)
        real(real64) :: local_contributions(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: total_contributions(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: expected_local(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: expected_total(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: tmp_factors(N_TIMEPOINTS, N_SELECTED_FACTORS), tmp_dependent(N_TIMEPOINTS)
        integer(int32) :: factor_indices(N_SELECTED_FACTORS), dependent_indices(N_SELECTED_DEPENDENTS), ierr

        call fill_three_factors(trajectories)
        call fill_expected_three_factors(expected_local, expected_total)
        factor_indices = [3, 1]
        dependent_indices = [2, 3]
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN
        tmp_factors = 1.0e300_real64
        tmp_dependent = -1.0e300_real64

        call compute_all_contributions_expert(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, &
                                              N_SELECTED_FACTORS, dependent_indices, N_SELECTED_DEPENDENTS, &
                                              MODE_BASELINE_RAW, local_contributions, total_contributions, &
                                              tmp_factors, tmp_dependent, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_all_contributions_expert: ierr")
        call assert_equal_array_real(local_contributions, expected_local, size(expected_local, kind=int32), &
                                     0.0_real64, "test_compute_all_contributions_expert: local_contributions", &
                                     N_TIMEPOINTS)
        call assert_equal_array_real(total_contributions, expected_total, size(expected_total, kind=int32), &
                                     0.0_real64, "test_compute_all_contributions_expert: total_contributions", &
                                     N_SELECTED_FACTORS)
    end subroutine test_compute_all_contributions_expert

    !> Calls the plain tier for one factor against one dependent in a single sample and checks
    !| the result exactly.
    subroutine check_single_pair(n_factors, n_timepoints, trajectories, factor_index, dependent_index, &
                                 baseline_mode, expected_local, expected_total, msg)
        integer(int32), intent(in) :: n_factors, n_timepoints, factor_index, dependent_index, baseline_mode
        real(real64), intent(in) :: trajectories(n_factors, 1, n_timepoints)
        real(real64), intent(in) :: expected_local(n_timepoints), expected_total
        character(*), intent(in) :: msg
        real(real64) :: local_contributions(n_timepoints, 1, 1, 1), total_contributions(1, 1, 1)
        integer(int32) :: factor_indices(1), dependent_indices(1), ierr

        factor_indices = factor_index
        dependent_indices = dependent_index
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN
        call compute_all_contributions(trajectories, n_factors, 1, n_timepoints, factor_indices, 1, &
                                       dependent_indices, 1, baseline_mode, local_contributions, &
                                       total_contributions, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(local_contributions, expected_local, n_timepoints, 0.0_real64, &
                                     msg//": local_contributions")
        call assert_equal_real(total_contributions(1, 1, 1), expected_total, 0.0_real64, msg//": total_contributions")
    end subroutine check_single_pair

    !> MODE_BASELINE_MEAN reaches the pair: factor [1, 2, 3, 6] (mean 3, deviations
    !| [-2, -1, 0, 3]) against dependent [8, 2, 4, 10] (mean 6, deviations [2, -4, -2, 4]) is
    !| [-4, 4, 0, 12], total 12.
    subroutine test_compute_all_contributions_mean()
        real(real64) :: trajectories(2, 1, 4), expected_local(4)

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64, 3.0_real64, 6.0_real64]
        trajectories(2, 1, :) = [8.0_real64, 2.0_real64, 4.0_real64, 10.0_real64]
        expected_local = [-4.0_real64, 4.0_real64, 0.0_real64, 12.0_real64]
        call check_single_pair(2, 4, trajectories, 1, 2, MODE_BASELINE_MEAN, expected_local, 12.0_real64, &
                               "test_compute_all_contributions_mean")
    end subroutine test_compute_all_contributions_mean

    !> MODE_BASELINE_MIN reaches the pair (kept from the old suite): factor [2, 4, 6] (minimum 2)
    !| against dependent [1, 3, 5] (minimum 1) is [0*0, 2*2, 4*4] = [0, 4, 16], total 20.
    subroutine test_compute_all_contributions_min()
        real(real64) :: trajectories(2, 1, 3), expected_local(3)

        trajectories(1, 1, :) = [2.0_real64, 4.0_real64, 6.0_real64]
        trajectories(2, 1, :) = [1.0_real64, 3.0_real64, 5.0_real64]
        expected_local = [0.0_real64, 4.0_real64, 16.0_real64]
        call check_single_pair(2, 3, trajectories, 1, 2, MODE_BASELINE_MIN, expected_local, 20.0_real64, &
                               "test_compute_all_contributions_min")
    end subroutine test_compute_all_contributions_min

    !> An index may be selected more than once, and each selection gets its own block: factors
    !| [2, 2] against dependents [1, 1, 1], raw, x1 [1, 2, 3] and x2 [4, 5, 6], fill all six
    !| blocks with x2 * x1 = [4, 10, 18], total 32.
    subroutine test_compute_all_contributions_repeated_indices()
        real(real64) :: trajectories(2, 1, 3), local_contributions(3, 2, 3, 1), total_contributions(2, 3, 1)
        real(real64) :: expected_local(3, 2, 3, 1), expected_total(2, 3, 1)
        integer(int32) :: factor_indices(2), dependent_indices(3), i_factor, i_dependent, ierr

        trajectories(1, 1, :) = [1.0_real64, 2.0_real64, 3.0_real64]
        trajectories(2, 1, :) = [4.0_real64, 5.0_real64, 6.0_real64]
        factor_indices = 2
        dependent_indices = 1
        do i_dependent = 1, 3
            do i_factor = 1, 2
                expected_local(:, i_factor, i_dependent, 1) = [4.0_real64, 10.0_real64, 18.0_real64]
            end do
        end do
        expected_total = 32.0_real64
        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN

        call compute_all_contributions(trajectories, 2, 1, 3, factor_indices, 2, dependent_indices, 3, &
                                       MODE_BASELINE_RAW, local_contributions, total_contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_all_contributions_repeated_indices: ierr")
        call assert_equal_array_real(local_contributions, expected_local, 18, 0.0_real64, &
                                     "test_compute_all_contributions_repeated_indices: local_contributions", 3)
        call assert_equal_array_real(total_contributions, expected_total, 6, 0.0_real64, &
                                     "test_compute_all_contributions_repeated_indices: total_contributions", 2)
    end subroutine test_compute_all_contributions_repeated_indices

    !> n_timepoints = 1 in two samples: x1 is 3 and 5, x2 is -2 and 4. Raw, x1 against x2 is
    !| -6 and 20; the mean of a single point is the point, so with MODE_BASELINE_MEAN both are 0.
    subroutine test_compute_all_contributions_single_timepoint()
        real(real64) :: trajectories(2, 2, 1), local_contributions(1, 1, 1, 2), total_contributions(1, 1, 2)
        real(real64) :: expected(2)
        integer(int32) :: factor_indices(1), dependent_indices(1), ierr

        trajectories(:, 1, 1) = [3.0_real64, -2.0_real64]
        trajectories(:, 2, 1) = [5.0_real64, 4.0_real64]
        factor_indices = 1
        dependent_indices = 2

        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN
        call compute_all_contributions(trajectories, 2, 2, 1, factor_indices, 1, dependent_indices, 1, &
                                       MODE_BASELINE_RAW, local_contributions, total_contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_all_contributions_single_timepoint: raw, ierr")
        expected = [-6.0_real64, 20.0_real64]
        call assert_equal_array_real(local_contributions, expected, 2, 0.0_real64, &
                                     "test_compute_all_contributions_single_timepoint: raw, local_contributions")
        call assert_equal_array_real(total_contributions, expected, 2, 0.0_real64, &
                                     "test_compute_all_contributions_single_timepoint: raw, total_contributions")

        local_contributions = UNWRITTEN
        total_contributions = UNWRITTEN
        call compute_all_contributions(trajectories, 2, 2, 1, factor_indices, 1, dependent_indices, 1, &
                                       MODE_BASELINE_MEAN, local_contributions, total_contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_all_contributions_single_timepoint: mean, ierr")
        expected = 0.0_real64
        call assert_equal_array_real(local_contributions, expected, 2, 0.0_real64, &
                                     "test_compute_all_contributions_single_timepoint: mean, local_contributions")
        call assert_equal_array_real(total_contributions, expected, 2, 0.0_real64, &
                                     "test_compute_all_contributions_single_timepoint: mean, total_contributions")
    end subroutine test_compute_all_contributions_single_timepoint

    !> Each dimension on its own, the others valid: `n_factors` (argument 2), `n_samples` (3),
    !| `n_timepoints` (4), `n_selected_factors` (6) and `n_selected_dependents` (8) are
    !| ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1.
    subroutine test_compute_all_contributions_dimensions()
        real(real64) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)
        real(real64) :: local_contributions(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: total_contributions(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        integer(int32) :: factor_indices(N_SELECTED_FACTORS), dependent_indices(N_SELECTED_DEPENDENTS)
        integer(int32) :: dims(5), positions(5), i_dim, size_value, ierr
        integer(int32), parameter :: SIZE_VALUES(2) = [0, -1]
        integer(int32), parameter :: EXPECTED_CODES(2) = [ERR_EMPTY_INPUT, ERR_INVALID_INPUT]
        character(*), parameter :: DIM_NAMES(5) = ["n_factors            ", "n_samples            ", &
                                                   "n_timepoints         ", "n_selected_factors   ", &
                                                   "n_selected_dependents"]
        integer :: i_value

        call fill_three_factors(trajectories)
        factor_indices = [1, 2]
        dependent_indices = [2, 3]
        positions = [2, 3, 4, 6, 8]

        do i_dim = 1, 5
            do i_value = 1, 2
                size_value = SIZE_VALUES(i_value)
                dims = [N_FACTORS, N_SAMPLES, N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS]
                dims(i_dim) = size_value
                call compute_all_contributions(trajectories, dims(1), dims(2), dims(3), factor_indices, dims(4), &
                                               dependent_indices, dims(5), MODE_BASELINE_RAW, local_contributions, &
                                               total_contributions, ierr)
                call assert_err(ierr, EXPECTED_CODES(i_value), "test_compute_all_contributions_dimensions: " &
                                //trim(DIM_NAMES(i_dim))//" = "//merge(" 0", "-1", size_value == 0), &
                                positions(i_dim))
            end do
        end do
    end subroutine test_compute_all_contributions_dimensions

    !> `factor_indices` (argument 5) and `dependent_indices` (argument 7) lie in [1, n_factors]:
    !| 0 and n_factors + 1, the values just outside, are ERR_INVALID_INPUT. Each sits in the last
    !| element, so every selected index is checked. 1 and n_factors are accepted in
    !| `test_compute_all_contributions_values`.
    subroutine test_compute_all_contributions_index_bounds()
        real(real64) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)
        real(real64) :: local_contributions(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: total_contributions(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        integer(int32) :: factor_indices(N_SELECTED_FACTORS), dependent_indices(N_SELECTED_DEPENDENTS), ierr
        integer(int32) :: outside(2), i_outside
        character(*), parameter :: OUTSIDE_NAMES(2) = ["0            ", "n_factors + 1"]

        call fill_three_factors(trajectories)
        outside = [0, N_FACTORS + 1]

        do i_outside = 1, 2
            factor_indices(1) = 1
            factor_indices(2) = outside(i_outside)
            dependent_indices = [2, 3]
            call compute_all_contributions(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, &
                                           N_SELECTED_FACTORS, dependent_indices, N_SELECTED_DEPENDENTS, &
                                           MODE_BASELINE_RAW, local_contributions, total_contributions, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_all_contributions_index_bounds: factor index " &
                            //trim(OUTSIDE_NAMES(i_outside)), 5)

            factor_indices = [1, 2]
            dependent_indices(1) = 2
            dependent_indices(2) = outside(i_outside)
            call compute_all_contributions(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, &
                                           N_SELECTED_FACTORS, dependent_indices, N_SELECTED_DEPENDENTS, &
                                           MODE_BASELINE_RAW, local_contributions, total_contributions, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_all_contributions_index_bounds: dependent index " &
                            //trim(OUTSIDE_NAMES(i_outside)), 7)
        end do
    end subroutine test_compute_all_contributions_index_bounds

    !> `baseline_mode` (argument 9) is one of the three MODE_BASELINE_* values (1 to 3); 0 and 4,
    !| just outside, are ERR_INVALID_INPUT.
    subroutine test_compute_all_contributions_invalid_mode()
        real(real64) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)
        real(real64) :: local_contributions(N_TIMEPOINTS, N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        real(real64) :: total_contributions(N_SELECTED_FACTORS, N_SELECTED_DEPENDENTS, N_SAMPLES)
        integer(int32) :: factor_indices(N_SELECTED_FACTORS), dependent_indices(N_SELECTED_DEPENDENTS), ierr

        call fill_three_factors(trajectories)
        factor_indices = [1, 2]
        dependent_indices = [2, 3]
        call compute_all_contributions(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, &
                                       N_SELECTED_FACTORS, dependent_indices, N_SELECTED_DEPENDENTS, &
                                       0_int32, local_contributions, total_contributions, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_all_contributions_invalid_mode: 0", 9)
        call compute_all_contributions(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, &
                                       N_SELECTED_FACTORS, dependent_indices, N_SELECTED_DEPENDENTS, &
                                       4_int32, local_contributions, total_contributions, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_all_contributions_invalid_mode: 4", 9)
    end subroutine test_compute_all_contributions_invalid_mode

    !> NaN, +Inf and -Inf in `trajectories` (argument 1) are ERR_NAN_INF, even in a series that
    !| is not selected: the value sits in x2 of the last sample at the last time point, while
    !| only x1 is selected, as factor and as dependent.
    subroutine test_compute_all_contributions_rejects_nan_and_inf()
        real(real64) :: trajectories(N_FACTORS, N_SAMPLES, N_TIMEPOINTS)
        real(real64) :: local_contributions(N_TIMEPOINTS, 1, 1, N_SAMPLES), total_contributions(1, 1, N_SAMPLES)
        integer(int32) :: factor_indices(1), dependent_indices(1), ierr, i_value

        factor_indices = 1
        dependent_indices = 1
        do i_value = 1, N_NON_FINITE
            call fill_three_factors(trajectories)
            trajectories(2, N_SAMPLES, N_TIMEPOINTS) = non_finite_value(i_value)
            call compute_all_contributions(trajectories, N_FACTORS, N_SAMPLES, N_TIMEPOINTS, factor_indices, 1, &
                                           dependent_indices, 1, MODE_BASELINE_RAW, local_contributions, &
                                           total_contributions, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_all_contributions_rejects_nan_and_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_compute_all_contributions_rejects_nan_and_inf

end module mod_test_compute_all_contributions
