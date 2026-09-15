!> The `relative_axes_expression_from_expression_vector` cases: hand-derived shares of an
!| expression vector, an axis without expression, the zero vector (ERR_DIVISION_BY_ZERO), the
!| n_axes bound and the rejection of NaN and Inf.
!|
!| The procedure hands the vector to `compute_relative_axis_contributions`, whose cases cover the
!| arithmetic in depth; these establish that the entry point delivers it, and check its own input
!| validation. Test names drop the procedure's `relative_axes_` prefix to stay within 63 characters:
!| `test_expression_from_expression_vector_<case>`.
module mod_test_relative_axes_expression_from_expression_vector
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
    use tox_relative_axis_plane_tools
    use test_suite, only: test_case
    use tox_errors
    implicit none
    public

contains

    !> Get array of all available tests.
    function get_all_tests_relative_axes_expression_from_expression_vector() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(4))

        all_tests(1) = test_case("test_expression_from_expression_vector_values", &
                                 test_expression_from_expression_vector_values)
        all_tests(2) = test_case("test_expression_from_expression_vector_zero_vector", &
                                 test_expression_from_expression_vector_zero_vector)
        all_tests(3) = test_case("test_expression_from_expression_vector_dimensions", &
                                 test_expression_from_expression_vector_dimensions)
        all_tests(4) = test_case("test_expression_from_expression_vector_rejects_nan_and_inf", &
                                 test_expression_from_expression_vector_rejects_nan_and_inf)
    end function get_all_tests_relative_axes_expression_from_expression_vector

    !> The RAP-projected expression vector [-1/2, 1/4, 1/4] sums to 0 and its absolute values to 1,
    !| so the shares are the magnitudes themselves: 1/2, 1/4, 1/4. In [6, 0, -2] the second axis
    !| carries no expression and gets no share; the absolute values sum to 8, giving 6/8, 0, 2/8.
    subroutine test_expression_from_expression_vector_values()
        integer(int32) :: ierr
        real(real64), dimension(3) :: expression, contributions, expected

        expression = [-0.5d0, 0.25d0, 0.25d0]
        expected = [0.5d0, 0.25d0, 0.25d0]
        call relative_axes_expression_from_expression_vector(expression, 3, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_expression_from_expression_vector_values: ierr for [-1/2, 1/4, 1/4]")
        call assert_equal_array_real(contributions, expected, 3, 0d0, &
                                     "test_expression_from_expression_vector_values: shares of [-1/2, 1/4, 1/4]")

        expression = [6d0, 0d0, -2d0]
        expected = [0.75d0, 0d0, 0.25d0]
        call relative_axes_expression_from_expression_vector(expression, 3, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_expression_from_expression_vector_values: ierr for [6, 0, -2]")
        call assert_equal_array_real(contributions, expected, 3, 0d0, &
                                     "test_expression_from_expression_vector_values: shares of [6, 0, -2]")
    end subroutine test_expression_from_expression_vector_values

    !> A zero expression vector has nothing to share out: ERR_DIVISION_BY_ZERO.
    subroutine test_expression_from_expression_vector_zero_vector()
        integer(int32) :: ierr
        real(real64), dimension(3) :: expression, contributions

        expression = [0d0, 0d0, 0d0]
        call relative_axes_expression_from_expression_vector(expression, 3, contributions, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_expression_from_expression_vector_zero_vector: [0, 0, 0]")
    end subroutine test_expression_from_expression_vector_zero_vector

    !> n_axes, argument 2: 1 is the smallest valid value, [-7] giving [1]; 0 is ERR_EMPTY_INPUT and
    !| -1 ERR_INVALID_INPUT.
    subroutine test_expression_from_expression_vector_dimensions()
        integer(int32) :: ierr
        real(real64), dimension(1) :: expression, contributions, expected

        expression = [-7d0]
        expected = [1d0]

        call relative_axes_expression_from_expression_vector(expression, 1, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_expression_from_expression_vector_dimensions: n_axes = 1")
        call assert_equal_array_real(contributions, expected, 1, 0d0, &
                                     "test_expression_from_expression_vector_dimensions: share of [-7]")

        call relative_axes_expression_from_expression_vector(expression, 0, contributions, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_expression_from_expression_vector_dimensions: n_axes = 0", 2)
        call relative_axes_expression_from_expression_vector(expression, -1, contributions, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_expression_from_expression_vector_dimensions: n_axes = -1", 2)
    end subroutine test_expression_from_expression_vector_dimensions

    !> NaN, Inf and -Inf in `vec`, argument 1, are rejected up front: TOX takes no NaN or Inf.
    subroutine test_expression_from_expression_vector_rejects_nan_and_inf()
        integer(int32) :: ierr, i_bad
        real(real64), dimension(3) :: expression, contributions
        real(real64) :: bad(3)
        character(len=*), parameter :: bad_names(3) = ["NaN ", "Inf ", "-Inf"]

        bad(1) = ieee_value(1d0, ieee_quiet_nan)
        bad(2) = ieee_value(1d0, ieee_positive_inf)
        bad(3) = ieee_value(1d0, ieee_negative_inf)
        expression(1:2) = [1d0, -1d0]

        do i_bad = 1, size(bad)
            expression(3) = bad(i_bad)
            call relative_axes_expression_from_expression_vector(expression, 3, contributions, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_expression_from_expression_vector_rejects_nan_and_inf: must reject "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_expression_from_expression_vector_rejects_nan_and_inf

end module mod_test_relative_axes_expression_from_expression_vector
