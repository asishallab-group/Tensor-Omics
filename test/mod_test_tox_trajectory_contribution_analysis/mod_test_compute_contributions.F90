!> The `compute_contributions` cases: the per-point contributions (factor - its baseline) *
!| (dependent - its baseline) and their total, in each baseline mode, for series that move
!| together and apart, a single point, a constant series, and the mode, dimension and NaN/Inf
!| checks.
!|
!| Every value is a small integer or a multiple of 1/4 and every mean divides by a power of two,
!| so every expected value is exact and needs no tolerance.
module mod_test_compute_contributions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_contributions
    use tox_trajectory_contribution_analysis_impl, only: MODE_BASELINE_RAW, MODE_BASELINE_MIN, MODE_BASELINE_MEAN
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_contributions

contains

    !> Every case of `compute_contributions`.
    function get_all_tests_compute_contributions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(8))
        all_tests(1) = test_case("test_compute_contributions_raw", test_compute_contributions_raw)
        all_tests(2) = test_case("test_compute_contributions_min", test_compute_contributions_min)
        all_tests(3) = test_case("test_compute_contributions_mean", test_compute_contributions_mean)
        all_tests(4) = test_case("test_compute_contributions_single_timepoint", &
                                 test_compute_contributions_single_timepoint)
        all_tests(5) = test_case("test_compute_contributions_constant_series", &
                                 test_compute_contributions_constant_series)
        all_tests(6) = test_case("test_compute_contributions_invalid_mode", test_compute_contributions_invalid_mode)
        all_tests(7) = test_case("test_compute_contributions_dimensions", test_compute_contributions_dimensions)
        all_tests(8) = test_case("test_compute_contributions_rejects_nan_and_inf", &
                                 test_compute_contributions_rejects_nan_and_inf)
    end function get_all_tests_compute_contributions

    !> Calls the procedure and checks `ierr` OK and both outputs exactly. The outputs start at
    !| UNWRITTEN, so an element the procedure never writes shows.
    subroutine check_contributions(n_dims, factor, dependent, baseline_mode, expected_local, expected_total, msg)
        integer(int32), intent(in) :: n_dims, baseline_mode
        real(real64), intent(in) :: factor(n_dims), dependent(n_dims), expected_local(n_dims), expected_total
        character(*), intent(in) :: msg
        real(real64) :: local_contributions(n_dims), total_contribution
        integer(int32) :: ierr

        local_contributions = UNWRITTEN
        total_contribution = UNWRITTEN
        call compute_contributions(factor, dependent, n_dims, baseline_mode, local_contributions, &
                                   total_contribution, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(local_contributions, expected_local, n_dims, 0.0_real64, &
                                     msg//": local_contributions")
        call assert_equal_real(total_contribution, expected_total, 0.0_real64, msg//": total_contribution")
    end subroutine check_contributions

    !> MODE_BASELINE_RAW, baselines 0: factor [1, 2, 3, 4] times dependent [2, 1, 0, -3] is
    !| [2, 2, 0, -12], total -8.
    subroutine test_compute_contributions_raw()
        real(real64) :: factor(4), dependent(4), expected_local(4)

        factor = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent = [2.0_real64, 1.0_real64, 0.0_real64, -3.0_real64]
        expected_local = [2.0_real64, 2.0_real64, 0.0_real64, -12.0_real64]
        call check_contributions(4, factor, dependent, MODE_BASELINE_RAW, expected_local, -8.0_real64, &
                                 "test_compute_contributions_raw")
    end subroutine test_compute_contributions_raw

    !> MODE_BASELINE_MIN: factor [3, 5, 2, 4] has the minimum 2, dependent [1, 2, 0, -1] the
    !| minimum -1, so the deviations are [1, 3, 0, 2] and [2, 3, 1, 0]: [2, 9, 0, 0], total 11.
    !| Each series has its own minimum, at a different point, so a shared baseline would show.
    subroutine test_compute_contributions_min()
        real(real64) :: factor(4), dependent(4), expected_local(4)

        factor = [3.0_real64, 5.0_real64, 2.0_real64, 4.0_real64]
        dependent = [1.0_real64, 2.0_real64, 0.0_real64, -1.0_real64]
        expected_local = [2.0_real64, 9.0_real64, 0.0_real64, 0.0_real64]
        call check_contributions(4, factor, dependent, MODE_BASELINE_MIN, expected_local, 11.0_real64, &
                                 "test_compute_contributions_min")
    end subroutine test_compute_contributions_min

    !> MODE_BASELINE_MEAN, series that move apart: factor [1, 2, 3, 4] and dependent [4, 3, 2, 1]
    !| both have the mean 2.5, so the deviations are [-1.5, -0.5, 0.5, 1.5] and
    !| [1.5, 0.5, -0.5, -1.5]: [-2.25, -0.25, -0.25, -2.25], total -5. Each point is negative,
    !| as the two always lie on opposite sides of their means.
    subroutine test_compute_contributions_mean()
        real(real64) :: factor(4), dependent(4), expected_local(4)

        factor = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        dependent = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        expected_local = [-2.25_real64, -0.25_real64, -0.25_real64, -2.25_real64]
        call check_contributions(4, factor, dependent, MODE_BASELINE_MEAN, expected_local, -5.0_real64, &
                                 "test_compute_contributions_mean")
    end subroutine test_compute_contributions_mean

    !> n_dims = 1, the smallest valid size: factor [3], dependent [-2]. Raw it is 3 * -2 = -6;
    !| the minimum and the mean of one point are the point, so both other modes give 0.
    subroutine test_compute_contributions_single_timepoint()
        real(real64) :: factor(1), dependent(1), expected_local(1)

        factor = 3.0_real64
        dependent = -2.0_real64
        expected_local = -6.0_real64
        call check_contributions(1, factor, dependent, MODE_BASELINE_RAW, expected_local, -6.0_real64, &
                                 "test_compute_contributions_single_timepoint: raw")
        expected_local = 0.0_real64
        call check_contributions(1, factor, dependent, MODE_BASELINE_MIN, expected_local, 0.0_real64, &
                                 "test_compute_contributions_single_timepoint: min")
        call check_contributions(1, factor, dependent, MODE_BASELINE_MEAN, expected_local, 0.0_real64, &
                                 "test_compute_contributions_single_timepoint: mean")
    end subroutine test_compute_contributions_single_timepoint

    !> A series that is constant over time never deviates from its minimum or mean, so every
    !| contribution with it is 0, whether it is the factor or the dependent: [5, 5, 5, 5] against
    !| [1, 2, 3, 4] (kept from the old suite's "all elements equal" case).
    subroutine test_compute_contributions_constant_series()
        real(real64) :: constant(4), varying(4), expected_local(4)

        constant = 5.0_real64
        varying = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        expected_local = 0.0_real64
        call check_contributions(4, constant, varying, MODE_BASELINE_MIN, expected_local, 0.0_real64, &
                                 "test_compute_contributions_constant_series: constant factor, min")
        call check_contributions(4, constant, varying, MODE_BASELINE_MEAN, expected_local, 0.0_real64, &
                                 "test_compute_contributions_constant_series: constant factor, mean")
        call check_contributions(4, varying, constant, MODE_BASELINE_MIN, expected_local, 0.0_real64, &
                                 "test_compute_contributions_constant_series: constant dependent, min")
        call check_contributions(4, varying, constant, MODE_BASELINE_MEAN, expected_local, 0.0_real64, &
                                 "test_compute_contributions_constant_series: constant dependent, mean")
    end subroutine test_compute_contributions_constant_series

    !> `baseline_mode` (argument 4) is one of the three MODE_BASELINE_* values (1 to 3); 0 and 4,
    !| just outside, are ERR_INVALID_INPUT.
    subroutine test_compute_contributions_invalid_mode()
        real(real64) :: factor(2), dependent(2), local_contributions(2), total_contribution
        integer(int32) :: ierr

        factor = [1.0_real64, 2.0_real64]
        dependent = [3.0_real64, 4.0_real64]
        call compute_contributions(factor, dependent, 2, 0_int32, local_contributions, total_contribution, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_contributions_invalid_mode: 0", 4)
        call compute_contributions(factor, dependent, 2, 4_int32, local_contributions, total_contribution, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_contributions_invalid_mode: 4", 4)
    end subroutine test_compute_contributions_invalid_mode

    !> `n_dims` (argument 3) is ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1; 1 is in
    !| `test_compute_contributions_single_timepoint`.
    subroutine test_compute_contributions_dimensions()
        real(real64) :: factor(2), dependent(2), local_contributions(2), total_contribution
        integer(int32) :: ierr

        factor = [1.0_real64, 2.0_real64]
        dependent = [3.0_real64, 4.0_real64]
        call compute_contributions(factor, dependent, 0, MODE_BASELINE_MEAN, local_contributions, &
                                   total_contribution, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_contributions_dimensions: 0", 3)
        call compute_contributions(factor, dependent, -1, MODE_BASELINE_MEAN, local_contributions, &
                                   total_contribution, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_contributions_dimensions: -1", 3)
    end subroutine test_compute_contributions_dimensions

    !> NaN, +Inf and -Inf in `factor` (argument 1) or `dependent` (argument 2) are ERR_NAN_INF.
    !| Each sits in the last point, so the check must cover all n_dims values; the mode is
    !| MODE_BASELINE_RAW, whose baseline does not look at the values, so only the input check
    !| can refuse them.
    subroutine test_compute_contributions_rejects_nan_and_inf()
        real(real64) :: factor(3), dependent(3), local_contributions(3), total_contribution
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            factor = [1.0_real64, 2.0_real64, 3.0_real64]
            factor(3) = non_finite_value(i_value)
            dependent = [3.0_real64, 4.0_real64, 5.0_real64]
            call compute_contributions(factor, dependent, 3, MODE_BASELINE_RAW, local_contributions, &
                                       total_contribution, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_contributions_rejects_nan_and_inf: factor " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)

            factor = [1.0_real64, 2.0_real64, 3.0_real64]
            dependent = [3.0_real64, 4.0_real64, 5.0_real64]
            dependent(3) = non_finite_value(i_value)
            call compute_contributions(factor, dependent, 3, MODE_BASELINE_RAW, local_contributions, &
                                       total_contribution, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_compute_contributions_rejects_nan_and_inf: dependent " &
                            //trim(NON_FINITE_NAMES(i_value)), 2)
        end do
    end subroutine test_compute_contributions_rejects_nan_and_inf

end module mod_test_compute_contributions
