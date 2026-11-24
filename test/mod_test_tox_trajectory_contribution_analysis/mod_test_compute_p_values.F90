!> The `compute_p_values` cases: the fraction of permutations whose contribution is at least the
!| observed one, per time point and for the total, with ties, no permutation as high, negative
!| contributions, a single permutation and time point, and the dimension and NaN/Inf checks.
!|
!| Every case has one or four permutations, so every p-value is a multiple of 1/4 and exact.
module mod_test_compute_p_values
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_trajectory_contribution_analysis, only: compute_p_values
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_trajectory_contribution_analysis_fixtures, only: UNWRITTEN, N_NON_FINITE, NON_FINITE_NAMES, &
                                                                     non_finite_value
    implicit none
    private
    public :: get_all_tests_compute_p_values

contains

    !> Every case of `compute_p_values`.
    function get_all_tests_compute_p_values() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(9))
        all_tests(1) = test_case("test_compute_p_values_values", test_compute_p_values_values)
        all_tests(2) = test_case("test_compute_p_values_ties_count", test_compute_p_values_ties_count)
        all_tests(3) = test_case("test_compute_p_values_none_as_high", test_compute_p_values_none_as_high)
        all_tests(4) = test_case("test_compute_p_values_negative_contributions", &
                                 test_compute_p_values_negative_contributions)
        all_tests(5) = test_case("test_compute_p_values_total_from_totals", test_compute_p_values_total_from_totals)
        all_tests(6) = test_case("test_compute_p_values_single_permutation", &
                                 test_compute_p_values_single_permutation)
        all_tests(7) = test_case("test_compute_p_values_single_timepoint", test_compute_p_values_single_timepoint)
        all_tests(8) = test_case("test_compute_p_values_dimensions", test_compute_p_values_dimensions)
        all_tests(9) = test_case("test_compute_p_values_rejects_nan_and_inf", &
                                 test_compute_p_values_rejects_nan_and_inf)
    end function get_all_tests_compute_p_values

    !> Calls the procedure and checks `ierr` OK and both p-values exactly. The outputs start at
    !| UNWRITTEN, so an element the procedure never writes shows.
    subroutine check_p_values(n_timepoints, n_permutations, local_observed, total_observed, local_perm, &
                              total_perm, expected_local_p, expected_total_p, msg)
        integer(int32), intent(in) :: n_timepoints, n_permutations
        real(real64), intent(in) :: local_observed(n_timepoints), total_observed
        real(real64), intent(in) :: local_perm(n_timepoints, n_permutations), total_perm(n_permutations)
        real(real64), intent(in) :: expected_local_p(n_timepoints), expected_total_p
        character(*), intent(in) :: msg
        real(real64) :: local_p_values(n_timepoints), total_p_value
        integer(int32) :: ierr

        local_p_values = UNWRITTEN
        total_p_value = UNWRITTEN
        call compute_p_values(local_observed, total_observed, local_perm, total_perm, n_timepoints, n_permutations, &
                              local_p_values, total_p_value, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_array_real(local_p_values, expected_local_p, n_timepoints, 0.0_real64, msg//": local_p_values")
        call assert_equal_real(total_p_value, expected_total_p, 0.0_real64, msg//": total_p_value")
    end subroutine check_p_values

    !> Kept from the old suite: observed [2, 0, 2], total 4, against four permutations
    !| [1, 0, 1] (2), [2, 0, 2] (4), [3, 1, 3] (7), [0, 0, 0] (0). At least as high at point 1:
    !| permutations 2 and 3 (2/4); at point 2 all four (4/4); at point 3 two (2/4); the total:
    !| 4 and 7 (2/4).
    subroutine test_compute_p_values_values()
        real(real64) :: local_observed(3), local_perm(3, 4), total_perm(4), expected_local_p(3)

        local_observed = [2.0_real64, 0.0_real64, 2.0_real64]
        local_perm(:, 1) = [1.0_real64, 0.0_real64, 1.0_real64]
        local_perm(:, 2) = [2.0_real64, 0.0_real64, 2.0_real64]
        local_perm(:, 3) = [3.0_real64, 1.0_real64, 3.0_real64]
        local_perm(:, 4) = [0.0_real64, 0.0_real64, 0.0_real64]
        total_perm = [2.0_real64, 4.0_real64, 7.0_real64, 0.0_real64]
        expected_local_p = [0.5_real64, 1.0_real64, 0.5_real64]
        call check_p_values(3, 4, local_observed, 4.0_real64, local_perm, total_perm, expected_local_p, 0.5_real64, &
                            "test_compute_p_values_values")
    end subroutine test_compute_p_values_values

    !> A permutation equal to the observed contribution counts as at least as high: four
    !| permutations that all equal the observed [1, -2], total -1, give p = 1 everywhere.
    !| Counting only the strictly higher ones would give 0.
    subroutine test_compute_p_values_ties_count()
        real(real64) :: local_observed(2), local_perm(2, 4), total_perm(4), expected_local_p(2)
        integer(int32) :: i_perm

        local_observed = [1.0_real64, -2.0_real64]
        do i_perm = 1, 4
            local_perm(:, i_perm) = local_observed
        end do
        total_perm = -1.0_real64
        expected_local_p = 1.0_real64
        call check_p_values(2, 4, local_observed, -1.0_real64, local_perm, total_perm, expected_local_p, 1.0_real64, &
                            "test_compute_p_values_ties_count")
    end subroutine test_compute_p_values_ties_count

    !> No permutation as high as the observed [5, 5], total 10: the p-values are exactly 0, as
    !| the documented fraction gives (no pseudo-count is added).
    subroutine test_compute_p_values_none_as_high()
        real(real64) :: local_observed(2), local_perm(2, 4), total_perm(4), expected_local_p(2)

        local_observed = [5.0_real64, 5.0_real64]
        local_perm(:, 1) = [4.0_real64, 4.0_real64]
        local_perm(:, 2) = [-1.0_real64, 0.0_real64]
        local_perm(:, 3) = [4.5_real64, -5.0_real64]
        local_perm(:, 4) = [0.0_real64, 4.75_real64]
        total_perm = [8.0_real64, -1.0_real64, -0.5_real64, 4.75_real64]
        expected_local_p = 0.0_real64
        call check_p_values(2, 4, local_observed, 10.0_real64, local_perm, total_perm, expected_local_p, 0.0_real64, &
                            "test_compute_p_values_none_as_high")
    end subroutine test_compute_p_values_none_as_high

    !> Negative contributions compare by value, not by magnitude: observed [-3, -1], total -4,
    !| against [-4, -1] (-5), [-2, -5] (-7), [-3, 0] (-3), [-6, 2] (-4). At least -3 at point 1:
    !| -2 and -3 (2/4); at least -1 at point 2: -1, 0 and 2 (3/4); at least -4 in total: -3 and
    !| -4 (2/4).
    subroutine test_compute_p_values_negative_contributions()
        real(real64) :: local_observed(2), local_perm(2, 4), total_perm(4), expected_local_p(2)

        local_observed = [-3.0_real64, -1.0_real64]
        local_perm(:, 1) = [-4.0_real64, -1.0_real64]
        local_perm(:, 2) = [-2.0_real64, -5.0_real64]
        local_perm(:, 3) = [-3.0_real64, 0.0_real64]
        local_perm(:, 4) = [-6.0_real64, 2.0_real64]
        total_perm = [-5.0_real64, -7.0_real64, -3.0_real64, -4.0_real64]
        expected_local_p = [0.5_real64, 0.75_real64]
        call check_p_values(2, 4, local_observed, -4.0_real64, local_perm, total_perm, expected_local_p, 0.5_real64, &
                            "test_compute_p_values_negative_contributions")
    end subroutine test_compute_p_values_negative_contributions

    !> The total p-value compares the totals it is given, not the sums of the per-point
    !| contributions: the observed total is 0 although its points sum to 2, and the permuted
    !| totals [1, -1, 0, 2] (3 of 4 at least 0) differ from their points' sums [0, 0, 0, 0].
    subroutine test_compute_p_values_total_from_totals()
        real(real64) :: local_observed(2), local_perm(2, 4), total_perm(4), expected_local_p(2)
        integer(int32) :: i_perm

        local_observed = [1.0_real64, 1.0_real64]
        do i_perm = 1, 4
            local_perm(:, i_perm) = [0.0_real64, 0.0_real64]
        end do
        total_perm = [1.0_real64, -1.0_real64, 0.0_real64, 2.0_real64]
        expected_local_p = 0.0_real64
        call check_p_values(2, 4, local_observed, 0.0_real64, local_perm, total_perm, expected_local_p, 0.75_real64, &
                            "test_compute_p_values_total_from_totals")
    end subroutine test_compute_p_values_total_from_totals

    !> n_permutations = 1, the smallest valid count: each p-value is 0 or 1. Observed [1, 1],
    !| total 2, against [1, 0], total 3: [1, 0] and 1.
    subroutine test_compute_p_values_single_permutation()
        real(real64) :: local_observed(2), local_perm(2, 1), total_perm(1), expected_local_p(2)

        local_observed = [1.0_real64, 1.0_real64]
        local_perm(:, 1) = [1.0_real64, 0.0_real64]
        total_perm = 3.0_real64
        expected_local_p = [1.0_real64, 0.0_real64]
        call check_p_values(2, 1, local_observed, 2.0_real64, local_perm, total_perm, expected_local_p, 1.0_real64, &
                            "test_compute_p_values_single_permutation")
    end subroutine test_compute_p_values_single_permutation

    !> n_timepoints = 1, the smallest valid size: observed [3], total 3, against [4], [3], [2],
    !| [1] with the totals [1, 2, 3, 4]: 2/4 at the point, 2/4 for the total.
    subroutine test_compute_p_values_single_timepoint()
        real(real64) :: local_observed(1), local_perm(1, 4), total_perm(4), expected_local_p(1)

        local_observed = 3.0_real64
        local_perm(1, :) = [4.0_real64, 3.0_real64, 2.0_real64, 1.0_real64]
        total_perm = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        expected_local_p = 0.5_real64
        call check_p_values(1, 4, local_observed, 3.0_real64, local_perm, total_perm, expected_local_p, 0.5_real64, &
                            "test_compute_p_values_single_timepoint")
    end subroutine test_compute_p_values_single_timepoint

    !> `n_timepoints` (argument 5) and `n_permutations` (6) are ERR_EMPTY_INPUT at 0 and
    !| ERR_INVALID_INPUT at -1, each with the other valid.
    subroutine test_compute_p_values_dimensions()
        real(real64) :: local_observed(2), local_perm(2, 2), total_perm(2), local_p_values(2), total_p_value
        integer(int32) :: ierr

        local_observed = 1.0_real64
        local_perm = 1.0_real64
        total_perm = 2.0_real64
        call compute_p_values(local_observed, 2.0_real64, local_perm, total_perm, 0, 2, local_p_values, &
                              total_p_value, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_p_values_dimensions: n_timepoints = 0", 5)
        call compute_p_values(local_observed, 2.0_real64, local_perm, total_perm, -1, 2, local_p_values, &
                              total_p_value, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_p_values_dimensions: n_timepoints = -1", 5)
        call compute_p_values(local_observed, 2.0_real64, local_perm, total_perm, 2, 0, local_p_values, &
                              total_p_value, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_p_values_dimensions: n_permutations = 0", 6)
        call compute_p_values(local_observed, 2.0_real64, local_perm, total_perm, 2, -1, local_p_values, &
                              total_p_value, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_p_values_dimensions: n_permutations = -1", 6)
    end subroutine test_compute_p_values_dimensions

    !> NaN, +Inf and -Inf are ERR_NAN_INF in each input: `local_contributions_observed`
    !| (argument 1), `total_contribution_observed` (2), `local_contributions_perm_test` (3) and
    !| `total_contributions_perm_test` (4). In the arrays the value sits in the last element.
    subroutine test_compute_p_values_rejects_nan_and_inf()
        real(real64) :: local_observed(2), local_perm(2, 2), total_perm(2), local_p_values(2), total_p_value
        integer(int32) :: ierr, i_value, i_argument
        character(*), parameter :: ARGUMENT_NAMES(4) = ["local_contributions_observed ", &
                                                        "total_contribution_observed  ", &
                                                        "local_contributions_perm_test", &
                                                        "total_contributions_perm_test"]
        real(real64) :: total_observed

        do i_argument = 1, 4
            do i_value = 1, N_NON_FINITE
                local_observed = 1.0_real64
                total_observed = 2.0_real64
                local_perm = 1.0_real64
                total_perm = 2.0_real64
                select case (i_argument)
                case (1)
                    local_observed(2) = non_finite_value(i_value)
                case (2)
                    total_observed = non_finite_value(i_value)
                case (3)
                    local_perm(2, 2) = non_finite_value(i_value)
                case (4)
                    total_perm(2) = non_finite_value(i_value)
                end select
                call compute_p_values(local_observed, total_observed, local_perm, total_perm, 2, 2, local_p_values, &
                                      total_p_value, ierr)
                call assert_err(ierr, ERR_NAN_INF, "test_compute_p_values_rejects_nan_and_inf: " &
                                //trim(ARGUMENT_NAMES(i_argument))//" "//trim(NON_FINITE_NAMES(i_value)), i_argument)
            end do
        end do
    end subroutine test_compute_p_values_rejects_nan_and_inf

end module mod_test_compute_p_values
