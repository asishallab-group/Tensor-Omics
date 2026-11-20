!> The `relative_axes_changes_from_shift_vector` cases: hand-derived shares of a shift vector,
!| a reversed shift sharing out the same way, the zero shift (ERR_DIVISION_BY_ZERO), the n_axes
!| bound and the rejection of NaN and Inf.
!|
!| The procedure hands the vector to `compute_relative_axis_contributions`, whose cases cover the
!| arithmetic in depth; these establish that the entry point delivers it, and check its own input
!| validation. Test names drop the procedure's `relative_axes_` prefix to stay within 63 characters:
!| `test_changes_from_shift_vector_<case>`.
module mod_test_relative_axes_changes_from_shift_vector
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
    function get_all_tests_relative_axes_changes_from_shift_vector() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(4))

        all_tests(1) = test_case("test_changes_from_shift_vector_values", test_changes_from_shift_vector_values)
        all_tests(2) = test_case("test_changes_from_shift_vector_zero_shift", test_changes_from_shift_vector_zero_shift)
        all_tests(3) = test_case("test_changes_from_shift_vector_dimensions", test_changes_from_shift_vector_dimensions)
        all_tests(4) = test_case("test_changes_from_shift_vector_rejects_nan_and_inf", &
                                 test_changes_from_shift_vector_rejects_nan_and_inf)
    end function get_all_tests_relative_axes_changes_from_shift_vector

    !> The shift [-4, 1, 1, 2] between two RAP-projected vectors sums to 0; its absolute values sum
    !| to 8, so axis 1 carries 4/8 of the change, axes 2 and 3 1/8 each and axis 4 2/8. The
    !| reversed shift [4, -1, -1, -2] moves along the same axes by the same amounts, so it shares
    !| out identically.
    subroutine test_changes_from_shift_vector_values()
        integer(int32) :: ierr
        real(real64), dimension(4) :: shift, contributions, expected

        expected = [0.5d0, 0.125d0, 0.125d0, 0.25d0]

        shift = [-4d0, 1d0, 1d0, 2d0]
        call relative_axes_changes_from_shift_vector(shift, 4, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_changes_from_shift_vector_values: ierr for [-4, 1, 1, 2]")
        call assert_equal_array_real(contributions, expected, 4, 0d0, &
                                     "test_changes_from_shift_vector_values: shares of [-4, 1, 1, 2]")

        shift = [4d0, -1d0, -1d0, -2d0]
        call relative_axes_changes_from_shift_vector(shift, 4, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_changes_from_shift_vector_values: ierr for [4, -1, -1, -2]")
        call assert_equal_array_real(contributions, expected, 4, 0d0, &
                                     "test_changes_from_shift_vector_values: shares of the reversed shift [4, -1, -1, -2]")
    end subroutine test_changes_from_shift_vector_values

    !> A gene whose origin and target coincide has not moved: the zero shift has no change to share
    !| out, which is ERR_DIVISION_BY_ZERO.
    subroutine test_changes_from_shift_vector_zero_shift()
        integer(int32) :: ierr
        real(real64), dimension(3) :: shift, contributions

        shift = [0d0, 0d0, 0d0]
        call relative_axes_changes_from_shift_vector(shift, 3, contributions, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_changes_from_shift_vector_zero_shift: [0, 0, 0]")
    end subroutine test_changes_from_shift_vector_zero_shift

    !> n_axes, argument 2: 1 is the smallest valid value, [2] giving [1]; 0 is ERR_EMPTY_INPUT and
    !| -1 ERR_INVALID_INPUT.
    subroutine test_changes_from_shift_vector_dimensions()
        integer(int32) :: ierr
        real(real64), dimension(1) :: shift, contributions, expected

        shift = [2d0]
        expected = [1d0]

        call relative_axes_changes_from_shift_vector(shift, 1, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_changes_from_shift_vector_dimensions: n_axes = 1")
        call assert_equal_array_real(contributions, expected, 1, 0d0, "test_changes_from_shift_vector_dimensions: share of [2]")

        call relative_axes_changes_from_shift_vector(shift, 0, contributions, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_changes_from_shift_vector_dimensions: n_axes = 0", 2)
        call relative_axes_changes_from_shift_vector(shift, -1, contributions, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_changes_from_shift_vector_dimensions: n_axes = -1", 2)
    end subroutine test_changes_from_shift_vector_dimensions

    !> NaN, Inf and -Inf in `vec`, argument 1, are rejected up front: TOX takes no NaN or Inf.
    subroutine test_changes_from_shift_vector_rejects_nan_and_inf()
        integer(int32) :: ierr, i_bad
        real(real64), dimension(3) :: shift, contributions
        real(real64) :: bad(3)
        character(len=*), parameter :: bad_names(3) = ["NaN ", "Inf ", "-Inf"]

        bad(1) = ieee_value(1d0, ieee_quiet_nan)
        bad(2) = ieee_value(1d0, ieee_positive_inf)
        bad(3) = ieee_value(1d0, ieee_negative_inf)
        shift(1:2) = [1d0, -1d0]

        do i_bad = 1, size(bad)
            shift(3) = bad(i_bad)
            call relative_axes_changes_from_shift_vector(shift, 3, contributions, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_changes_from_shift_vector_rejects_nan_and_inf: must reject "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_changes_from_shift_vector_rejects_nan_and_inf

end module mod_test_relative_axes_changes_from_shift_vector
