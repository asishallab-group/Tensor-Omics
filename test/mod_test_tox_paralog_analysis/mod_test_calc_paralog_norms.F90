!> The `calc_paralog_norms` cases: the euclidean norm of each gene's vector, the producer the plain
!| `detect_subfunctionalization` computes its norms with.
!|
!| The vectors are Pythagorean ([3, 4] is 5 long, [2, 3, 6] 7) times powers of two, so every
!| expected norm is exact and every comparison has tolerance 0, also where the squares themselves
!| would overflow or underflow: the norm is f42's scaled one.
module mod_test_calc_paralog_norms
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use tox_paralog_analysis, only: calc_paralog_norms
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_paralog_analysis_fixtures, only: N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    implicit none
    private
    public :: get_all_tests_calc_paralog_norms

    !> What a norm holds before a call, so that one the procedure never writes shows.
    real(real64), parameter :: UNWRITTEN = -777.0_real64

contains

    !> Every case of `calc_paralog_norms`.
    function get_all_tests_calc_paralog_norms() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(5))
        all_tests(1) = test_case("test_calc_paralog_norms_values", test_calc_paralog_norms_values)
        all_tests(2) = test_case("test_calc_paralog_norms_one_of_each", test_calc_paralog_norms_one_of_each)
        all_tests(3) = test_case("test_calc_paralog_norms_extreme_magnitudes", test_calc_paralog_norms_extreme_magnitudes)
        all_tests(4) = test_case("test_calc_paralog_norms_dimensions", test_calc_paralog_norms_dimensions)
        all_tests(5) = test_case("test_calc_paralog_norms_rejects_nan_and_inf", test_calc_paralog_norms_rejects_nan_and_inf)
    end function get_all_tests_calc_paralog_norms

    !> Four genes in three dimensions, one norm each, in gene order: [3, 4, 0] is 5, [-2, 3, -6] is
    !| 7, [0, 0, 0] is 0, and [0.5, 0, 0] is 0.5.
    subroutine test_calc_paralog_norms_values()
        real(real64) :: genes(3, 4), norms(4), expected(4)
        integer(int32) :: ierr

        genes = reshape([3.0_real64, 4.0_real64, 0.0_real64, -2.0_real64, 3.0_real64, -6.0_real64, &
                         0.0_real64, 0.0_real64, 0.0_real64, 0.5_real64, 0.0_real64, 0.0_real64], [3, 4])
        expected = [5.0_real64, 7.0_real64, 0.0_real64, 0.5_real64]
        norms = UNWRITTEN
        call calc_paralog_norms(genes, 3, 4, norms, ierr)
        call assert_err(ierr, ERR_OK, "test_calc_paralog_norms_values: ierr")
        call assert_equal_array_real(norms, expected, 4, 0.0_real64, "test_calc_paralog_norms_values")
    end subroutine test_calc_paralog_norms_values

    !> The smallest problem, one gene in one dimension: [-2.5] is 2.5 long.
    subroutine test_calc_paralog_norms_one_of_each()
        real(real64) :: genes(1, 1), norms(1)
        integer(int32) :: ierr

        genes = -2.5_real64
        norms = UNWRITTEN
        call calc_paralog_norms(genes, 1, 1, norms, ierr)
        call assert_err(ierr, ERR_OK, "test_calc_paralog_norms_one_of_each: ierr")
        call assert_equal_real(norms(1), 2.5_real64, 0.0_real64, "test_calc_paralog_norms_one_of_each")
    end subroutine test_calc_paralog_norms_one_of_each

    !> Magnitudes whose squares leave the normal range, each norm still exact:
    !| - [3, 4]*2**600 is 5*2**600 (the squares would pass `huge`);
    !| - [3, 4]*2**-600 is 5*2**-600 (the squares would underflow to 0);
    !| - [huge, 0] is `huge`;
    !| - [huge, huge], whose true norm sqrt(2)*huge exceeds `huge`, is +Inf, as f42's `norm`
    !|   documents.
    subroutine test_calc_paralog_norms_extreme_magnitudes()
        real(real64), parameter :: LARGE = 2.0_real64**600, SMALL = 2.0_real64**(-600), HUGE_VALUE = huge(1.0_real64)
        real(real64) :: genes(2, 4), norms(4), expected(4)
        integer(int32) :: ierr

        genes(:, 1) = [3.0_real64, 4.0_real64]*LARGE
        genes(:, 2) = [3.0_real64, 4.0_real64]*SMALL
        genes(:, 3) = [HUGE_VALUE, 0.0_real64]
        genes(:, 4) = [HUGE_VALUE, HUGE_VALUE]
        expected = [5.0_real64*LARGE, 5.0_real64*SMALL, HUGE_VALUE, ieee_value(1.0_real64, ieee_positive_inf)]
        norms = UNWRITTEN
        call calc_paralog_norms(genes, 2, 4, norms, ierr)
        call assert_err(ierr, ERR_OK, "test_calc_paralog_norms_extreme_magnitudes: ierr")
        call assert_equal_array_real(norms, expected, 4, 0.0_real64, "test_calc_paralog_norms_extreme_magnitudes")
    end subroutine test_calc_paralog_norms_extreme_magnitudes

    !> `n_dims` (argument 2) and `n_genes` (argument 3) are ERR_EMPTY_INPUT at 0 and
    !| ERR_INVALID_INPUT at -1. 1 is in `test_calc_paralog_norms_one_of_each`.
    subroutine test_calc_paralog_norms_dimensions()
        real(real64) :: genes(1, 1), norms(1)
        integer(int32) :: ierr

        genes = 1.0_real64
        call calc_paralog_norms(genes, 0, 1, norms, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_calc_paralog_norms_dimensions: n_dims 0", 2)
        call calc_paralog_norms(genes, -1, 1, norms, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_paralog_norms_dimensions: n_dims -1", 2)
        call calc_paralog_norms(genes, 1, 0, norms, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_calc_paralog_norms_dimensions: n_genes 0", 3)
        call calc_paralog_norms(genes, 1, -1, norms, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_paralog_norms_dimensions: n_genes -1", 3)
    end subroutine test_calc_paralog_norms_dimensions

    !> NaN, +Inf and -Inf in `genes` (argument 1), in the last element, are ERR_NAN_INF.
    subroutine test_calc_paralog_norms_rejects_nan_and_inf()
        real(real64) :: genes(2, 2), norms(2)
        integer(int32) :: ierr, i_value

        do i_value = 1, N_NON_FINITE
            genes = 1.0_real64
            genes(2, 2) = non_finite_value(i_value)
            call calc_paralog_norms(genes, 2, 2, norms, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_calc_paralog_norms_rejects_nan_and_inf: " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
        end do
    end subroutine test_calc_paralog_norms_rejects_nan_and_inf

end module mod_test_calc_paralog_norms
