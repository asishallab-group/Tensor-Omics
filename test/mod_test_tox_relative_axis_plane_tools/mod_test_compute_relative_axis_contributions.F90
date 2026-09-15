!> The `compute_relative_axis_contributions` cases: hand-derived shares |v_i|/sum|v| for vectors
!| whose components are exact in binary, the sign and scale of a vector changing nothing, a single
!| axis, all-equal components, many axes, the zero vector (ERR_DIVISION_BY_ZERO), the n_axes
!| bound and the rejection of NaN and Inf.
!|
!| The share depends only on the direction of `vec`, not on its length, so the cases use short
!| integer vectors rather than normalized ones: normalizing would only round the input. Most sum
!| their components to zero, as a RAP-projected vector does, and their absolute values to a power
!| of two, so every share is exact.
module mod_test_compute_relative_axis_contributions
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
    function get_all_tests_compute_relative_axis_contributions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(9))

        all_tests(1) = test_case("test_compute_relative_axis_contributions_values", &
                                 test_compute_relative_axis_contributions_values)
        all_tests(2) = test_case("test_compute_relative_axis_contributions_thirds_and_sixths", &
                                 test_compute_relative_axis_contributions_thirds_and_sixths)
        all_tests(3) = test_case("test_compute_relative_axis_contributions_signs", &
                                 test_compute_relative_axis_contributions_signs)
        all_tests(4) = test_case("test_compute_relative_axis_contributions_scale", &
                                 test_compute_relative_axis_contributions_scale)
        all_tests(5) = test_case("test_compute_relative_axis_contributions_one_nonzero_axis", &
                                 test_compute_relative_axis_contributions_one_nonzero_axis)
        all_tests(6) = test_case("test_compute_relative_axis_contributions_all_equal", &
                                 test_compute_relative_axis_contributions_all_equal)
        all_tests(7) = test_case("test_compute_relative_axis_contributions_zero_vector", &
                                 test_compute_relative_axis_contributions_zero_vector)
        all_tests(8) = test_case("test_compute_relative_axis_contributions_dimensions", &
                                 test_compute_relative_axis_contributions_dimensions)
        all_tests(9) = test_case("test_compute_relative_axis_contributions_rejects_nan_and_inf", &
                                 test_compute_relative_axis_contributions_rejects_nan_and_inf)
    end function get_all_tests_compute_relative_axis_contributions

    !> A RAP-projected vector: [1, 3, -4] sums to 0, and its absolute values to 1 + 3 + 4 = 8, so
    !| the shares are 1/8, 3/8 and 4/8, all exact.
    subroutine test_compute_relative_axis_contributions_values()
        integer(int32) :: ierr
        real(real64), dimension(3) :: vec, contributions, expected

        vec = [1d0, 3d0, -4d0]
        expected = [0.125d0, 0.375d0, 0.5d0]

        call compute_relative_axis_contributions(vec, 3, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_values: ierr")
        call assert_equal_array_real(contributions, expected, 3, 0d0, &
                                     "test_compute_relative_axis_contributions_values: shares of [1, 3, -4]")
    end subroutine test_compute_relative_axis_contributions_values

    !> A total that is no power of two: [1, 2, -3] sums to 0, its absolute values to 6, so the
    !| shares are 1/6, 2/6 = 1/3 and 3/6 = 1/2.
    subroutine test_compute_relative_axis_contributions_thirds_and_sixths()
        integer(int32) :: ierr
        real(real64), dimension(3) :: vec, contributions, expected

        vec = [1d0, 2d0, -3d0]
        expected = [1d0/6d0, 1d0/3d0, 0.5d0]

        call compute_relative_axis_contributions(vec, 3, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_thirds_and_sixths: ierr")
        ! A few ulps: an optimizer may divide by 6 as a multiplication with the rounded 1/6, which
        ! is one rounding more than the correctly rounded quotient the expectation holds.
        call assert_equal_array_real(contributions, expected, 3, 1d-15, &
                                     "test_compute_relative_axis_contributions_thirds_and_sixths: shares of [1, 2, -3]")
    end subroutine test_compute_relative_axis_contributions_thirds_and_sixths

    !> Only the magnitude of a component counts: [1, 3, 4], its negation [-1, -3, -4] and the mixed
    !| [1, -3, 4] all have absolute values summing to 8 and give 1/8, 3/8, 4/8.
    subroutine test_compute_relative_axis_contributions_signs()
        integer(int32) :: ierr, i_vec
        real(real64), dimension(3) :: contributions, expected
        real(real64), dimension(3, 3) :: vecs
        character(len=*), parameter :: vec_names(3) = ["[1, 3, 4]   ", "[-1, -3, -4]", "[1, -3, 4]  "]

        vecs(:, 1) = [1d0, 3d0, 4d0]
        vecs(:, 2) = [-1d0, -3d0, -4d0]
        vecs(:, 3) = [1d0, -3d0, 4d0]
        expected = [0.125d0, 0.375d0, 0.5d0]

        do i_vec = 1, 3
            call compute_relative_axis_contributions(vecs(:, i_vec), 3, contributions, ierr)
            call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_signs: ierr for "//trim(vec_names(i_vec)))
            call assert_equal_array_real(contributions, expected, 3, 0d0, &
                                         "test_compute_relative_axis_contributions_signs: shares of "//trim(vec_names(i_vec)))
        end do
    end subroutine test_compute_relative_axis_contributions_signs

    !> Scaling a vector changes no share, down to 2**-30 and up to 2**1000: [s, -3s] gives 1/4
    !| and 3/4 for every power of two s, and 4*2**1000 is still far below huge.
    subroutine test_compute_relative_axis_contributions_scale()
        integer(int32) :: ierr, i_scale
        real(real64), dimension(2) :: vec, contributions, expected
        real(real64) :: scales(3)

        scales = [2d0**(-30), 1d0, 2d0**1000]
        expected = [0.25d0, 0.75d0]

        do i_scale = 1, size(scales)
            vec(1) = scales(i_scale)
            vec(2) = -3d0*scales(i_scale)
            call compute_relative_axis_contributions(vec, 2, contributions, ierr)
            call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_scale: ierr")
            call assert_equal_array_real(contributions, expected, 2, 0d0, &
                                         "test_compute_relative_axis_contributions_scale: shares of [s, -3s]")
        end do
    end subroutine test_compute_relative_axis_contributions_scale

    !> A vector along a single axis puts everything on it: [0, 5, 0] and [0, -5, 0] give [0, 1, 0].
    subroutine test_compute_relative_axis_contributions_one_nonzero_axis()
        integer(int32) :: ierr
        real(real64), dimension(3) :: vec, contributions, expected

        expected = [0d0, 1d0, 0d0]

        vec = [0d0, 5d0, 0d0]
        call compute_relative_axis_contributions(vec, 3, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_one_nonzero_axis: ierr for [0, 5, 0]")
        call assert_equal_array_real(contributions, expected, 3, 0d0, &
                                     "test_compute_relative_axis_contributions_one_nonzero_axis: shares of [0, 5, 0]")

        vec = [0d0, -5d0, 0d0]
        call compute_relative_axis_contributions(vec, 3, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_one_nonzero_axis: ierr for [0, -5, 0]")
        call assert_equal_array_real(contributions, expected, 3, 0d0, &
                                     "test_compute_relative_axis_contributions_one_nonzero_axis: shares of [0, -5, 0]")
    end subroutine test_compute_relative_axis_contributions_one_nonzero_axis

    !> Components of equal magnitude share equally: the normalized RAP vector [1/2, 1/2, -1/2, -1/2]
    !| gives 1/4 each, and so does [2, -2, 2, -2]; 128 alternating +1 and -1 give 1/128 each (a sum
    !| of 128 ones is exact, and so is 1/128).
    subroutine test_compute_relative_axis_contributions_all_equal()
        integer(int32), parameter :: n_many = 128
        integer(int32) :: ierr, i_axis
        real(real64), dimension(4) :: vec, contributions, expected
        real(real64), dimension(n_many) :: many_vec, many_contributions, many_expected

        expected = 0.25d0

        vec = [0.5d0, 0.5d0, -0.5d0, -0.5d0]
        call compute_relative_axis_contributions(vec, 4, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_all_equal: ierr for [1/2, 1/2, -1/2, -1/2]")
        call assert_equal_array_real(contributions, expected, 4, 0d0, &
                                     "test_compute_relative_axis_contributions_all_equal: shares of [1/2, 1/2, -1/2, -1/2]")

        vec = [2d0, -2d0, 2d0, -2d0]
        call compute_relative_axis_contributions(vec, 4, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_all_equal: ierr for [2, -2, 2, -2]")
        call assert_equal_array_real(contributions, expected, 4, 0d0, &
                                     "test_compute_relative_axis_contributions_all_equal: shares of [2, -2, 2, -2]")

        do i_axis = 1, n_many
            many_vec(i_axis) = merge(1d0, -1d0, mod(i_axis, 2) == 1)
        end do
        many_expected = 1d0/128d0
        call compute_relative_axis_contributions(many_vec, n_many, many_contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_all_equal: ierr for 128 axes")
        call assert_equal_array_real(many_contributions, many_expected, n_many, 0d0, &
                                     "test_compute_relative_axis_contributions_all_equal: shares of 128 axes")
    end subroutine test_compute_relative_axis_contributions_all_equal

    !> A zero vector has no direction to share out: ERR_DIVISION_BY_ZERO, for [0, 0, 0], for signed
    !| zeros [-0, 0, -0] (whose magnitudes are 0 as well) and for a single axis [0].
    subroutine test_compute_relative_axis_contributions_zero_vector()
        integer(int32) :: ierr
        real(real64), dimension(3) :: vec, contributions
        real(real64), dimension(1) :: single_vec, single_contributions

        vec = [0d0, 0d0, 0d0]
        call compute_relative_axis_contributions(vec, 3, contributions, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_compute_relative_axis_contributions_zero_vector: [0, 0, 0]")

        vec = [-0d0, 0d0, -0d0]
        call compute_relative_axis_contributions(vec, 3, contributions, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_compute_relative_axis_contributions_zero_vector: [-0, 0, -0]")

        single_vec = [0d0]
        call compute_relative_axis_contributions(single_vec, 1, single_contributions, ierr)
        call assert_err(ierr, ERR_DIVISION_BY_ZERO, "test_compute_relative_axis_contributions_zero_vector: [0]")
    end subroutine test_compute_relative_axis_contributions_zero_vector

    !> n_axes, argument 2: 1 is the smallest valid value, and a single nonzero axis holds the whole
    !| share, [-3] giving [1]; 0 is ERR_EMPTY_INPUT and -1 ERR_INVALID_INPUT.
    subroutine test_compute_relative_axis_contributions_dimensions()
        integer(int32) :: ierr
        real(real64), dimension(1) :: vec, contributions, expected

        vec = [-3d0]
        expected = [1d0]

        call compute_relative_axis_contributions(vec, 1, contributions, ierr)
        call assert_err(ierr, ERR_OK, "test_compute_relative_axis_contributions_dimensions: n_axes = 1")
        call assert_equal_array_real(contributions, expected, 1, 0d0, &
                                     "test_compute_relative_axis_contributions_dimensions: share of [-3]")

        call compute_relative_axis_contributions(vec, 0, contributions, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_compute_relative_axis_contributions_dimensions: n_axes = 0", 2)
        call compute_relative_axis_contributions(vec, -1, contributions, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_compute_relative_axis_contributions_dimensions: n_axes = -1", 2)
    end subroutine test_compute_relative_axis_contributions_dimensions

    !> NaN, Inf and -Inf in `vec`, argument 1, are rejected up front rather than carried into the
    !| shares: TOX takes no NaN or Inf. The bad value sits last, behind valid ones.
    subroutine test_compute_relative_axis_contributions_rejects_nan_and_inf()
        integer(int32) :: ierr, i_bad
        real(real64), dimension(3) :: vec, contributions
        real(real64) :: bad(3)
        character(len=*), parameter :: bad_names(3) = ["NaN ", "Inf ", "-Inf"]

        bad(1) = ieee_value(1d0, ieee_quiet_nan)
        bad(2) = ieee_value(1d0, ieee_positive_inf)
        bad(3) = ieee_value(1d0, ieee_negative_inf)
        vec(1:2) = [1d0, -1d0]

        do i_bad = 1, size(bad)
            vec(3) = bad(i_bad)
            call compute_relative_axis_contributions(vec, 3, contributions, ierr)
            call assert_err(ierr, ERR_NAN_INF, &
                            "test_compute_relative_axis_contributions_rejects_nan_and_inf: must reject "//trim(bad_names(i_bad)), 1)
        end do
    end subroutine test_compute_relative_axis_contributions_rejects_nan_and_inf

end module mod_test_compute_relative_axis_contributions
