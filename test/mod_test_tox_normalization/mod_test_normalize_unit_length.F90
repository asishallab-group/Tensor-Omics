!> The `normalize_unit_length` cases: hand-derived unit vectors, the only norm it rejects (exactly
!| zero), magnitudes whose squares leave the real64 range, subnormal inputs, and the input checks.
module mod_test_normalize_unit_length
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use tox_normalization, only: normalize_unit_length
    use mod_test_tox_normalization_fixtures, only: subnormal_from_bits, flushes_subnormals, SMALLEST_SUBNORMAL_BITS, &
                                                   NINE_TENTHS_OF_TINY_BITS
    use tox_errors
    use test_suite, only: test_case
    implicit none

    real(real64), parameter :: TOL = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests_normalize_unit_length() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(11))
        all_tests(1) = test_case("test_normalize_unit_length_values", test_normalize_unit_length_values)
        all_tests(2) = test_case("test_normalize_unit_length_already_unit", test_normalize_unit_length_already_unit)
        all_tests(3) = test_case("test_normalize_unit_length_zero_vector", test_normalize_unit_length_zero_vector)
        all_tests(4) = test_case("test_normalize_unit_length_tiny_vector", test_normalize_unit_length_tiny_vector)
        all_tests(5) = test_case("test_normalize_unit_length_extreme_magnitudes", test_normalize_unit_length_extreme_magnitudes)
        all_tests(6) = test_case("test_normalize_unit_length_rejects_nan_and_inf", test_normalize_unit_length_rejects_nan_and_inf)
        all_tests(7) = test_case("test_normalize_unit_length_dimensions", test_normalize_unit_length_dimensions)
        all_tests(8) = test_case("test_normalize_unit_length_huge_entries", test_normalize_unit_length_huge_entries)
        all_tests(9) = test_case("test_normalize_unit_length_subnormal_squares", test_normalize_unit_length_subnormal_squares)
        all_tests(10) = test_case("test_normalize_unit_length_smallest_subnormal", &
                                  test_normalize_unit_length_smallest_subnormal)
        all_tests(11) = test_case("test_normalize_unit_length_subnormal_entry", test_normalize_unit_length_subnormal_entry)
    end function get_all_tests_normalize_unit_length

    !> [3, 4, -12] has norm 13, so it becomes [3, 4, -12]/13; a single entry keeps only its sign.
    subroutine test_normalize_unit_length_values()
        integer(int32) :: ierr
        real(real64) :: vector(3), expected(3), single(1)

        vector = [3.0_real64, 4.0_real64, -12.0_real64]
        expected = [3.0_real64/13.0_real64, 4.0_real64/13.0_real64, -12.0_real64/13.0_real64]
        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_values: [3, 4, -12] ierr")
        call assert_equal_array_real(vector, expected, 3, 4*TOL, "test_normalize_unit_length_values: [3, 4, -12]/13")

        single = [-5.0_real64]
        call normalize_unit_length(single, 1, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_values: [-5] ierr")
        call assert_equal_real(single(1), -1.0_real64, 0.0_real64, "test_normalize_unit_length_values: [-5] becomes [-1]")
    end subroutine test_normalize_unit_length_values

    !> A vector of norm 1 comes back unchanged, to rounding.
    subroutine test_normalize_unit_length_already_unit()
        integer(int32) :: ierr
        real(real64) :: vector(3), expected(3)

        ! 0.6164770879765119...^2 + 0.42^2 + 0.666^2 = 1
        vector = [0.61647708797651190566181131895072_real64, -0.42_real64, 0.666_real64]
        expected = vector
        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_already_unit: ierr")
        call assert_equal_array_real(vector, expected, 3, 2*TOL, "test_normalize_unit_length_already_unit: must stay unchanged")
    end subroutine test_normalize_unit_length_already_unit

    !> The zero vector has no direction: ERR_DIVISION_BY_ZERO, and the vector is left as it was.
    subroutine test_normalize_unit_length_zero_vector()
        integer(int32) :: ierr
        real(real64) :: vector(3), expected(3)

        vector = 0.0_real64
        expected = vector
        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_DIVISION_BY_ZERO, "test_normalize_unit_length_zero_vector: ierr")
        call assert_equal_array_real(vector, expected, 3, 0.0_real64, "test_normalize_unit_length_zero_vector: must stay unchanged")
    end subroutine test_normalize_unit_length_zero_vector

    !> Only an exactly zero norm is rejected: a norm of 1e-13 is small, not absent, and the vector
    !| still has a direction.
    subroutine test_normalize_unit_length_tiny_vector()
        integer(int32) :: ierr
        real(real64) :: vector(3), expected(3)

        vector = [1.0e-13_real64, 0.0_real64, 0.0_real64]
        expected = [1.0_real64, 0.0_real64, 0.0_real64]

        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_tiny_vector: a norm of 1e-13 is not zero")
        call assert_equal_array_real(vector, expected, 3, 0.0_real64, &
                                     "test_normalize_unit_length_tiny_vector: must normalize to [1, 0, 0]")
    end subroutine test_normalize_unit_length_tiny_vector

    !> Squaring an entry must neither overflow nor underflow the norm: 1e200 squared is past the
    !| largest real64, and 1e-170 squared is below the smallest, which made a non-zero vector's norm
    !| infinite or exactly zero. Both vectors have a direction and normalize to [1, 0, 0].
    subroutine test_normalize_unit_length_extreme_magnitudes()
        integer(int32) :: ierr
        real(real64) :: vector(3), expected(3)

        expected = [1.0_real64, 0.0_real64, 0.0_real64]

        vector = [1.0e200_real64, 0.0_real64, 0.0_real64]
        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_extreme_magnitudes: 1e200 must not overflow")
        call assert_equal_array_real(vector, expected, 3, 0.0_real64, &
                                     "test_normalize_unit_length_extreme_magnitudes: 1e200 must normalize to [1, 0, 0]")

        vector = [1.0e-170_real64, 0.0_real64, 0.0_real64]
        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_extreme_magnitudes: 1e-170 must not underflow")
        call assert_equal_array_real(vector, expected, 3, 0.0_real64, &
                                     "test_normalize_unit_length_extreme_magnitudes: 1e-170 must normalize to [1, 0, 0]")

        ! a 3-4-5 triangle far past the overflow point
        vector = [3.0e200_real64, 4.0e200_real64, 0.0_real64]
        expected = [0.6_real64, 0.8_real64, 0.0_real64]
        call normalize_unit_length(vector, 3, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_extreme_magnitudes: [3e200, 4e200] ierr")
        call assert_equal_array_real(vector, expected, 3, 4*TOL, &
                                     "test_normalize_unit_length_extreme_magnitudes: [3e200, 4e200] must normalize to [0.6, 0.8, 0]")
    end subroutine test_normalize_unit_length_extreme_magnitudes

    !> NaN and Inf are rejected up front: TOX writes no NaN.
    subroutine test_normalize_unit_length_rejects_nan_and_inf()
        integer(int32) :: ierr, i_bad
        real(real64) :: vector(3), bad(2)

        bad = [ieee_value(1.0_real64, ieee_quiet_nan), ieee_value(1.0_real64, ieee_positive_inf)]
        do i_bad = 1, size(bad)
            vector = [0.6_real64, -0.42_real64, bad(i_bad)]
            call normalize_unit_length(vector, 3, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_NAN_INF, &
                                  "test_normalize_unit_length_rejects_nan_and_inf: must reject "//merge("NaN", "Inf", i_bad == 1))
        end do
    end subroutine test_normalize_unit_length_rejects_nan_and_inf

    !> An empty vector is ERR_EMPTY_INPUT, a negative length ERR_INVALID_INPUT.
    subroutine test_normalize_unit_length_dimensions()
        integer(int32) :: ierr
        real(real64) :: vector(1)

        vector = 1.0_real64
        call normalize_unit_length(vector, 0, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, "test_normalize_unit_length_dimensions: n_dims = 0")
        call normalize_unit_length(vector, -1, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, "test_normalize_unit_length_dimensions: n_dims = -1")
    end subroutine test_normalize_unit_length_dimensions

    !> Entries whose norm exceeds huge still have a direction: [huge, huge] becomes
    !| [1/sqrt(2), 1/sqrt(2)], as the division happens in scaled coordinates.
    subroutine test_normalize_unit_length_huge_entries()
        integer(int32) :: ierr
        real(real64) :: vector(2), expected(2)

        vector = huge(1.0_real64)
        expected = 1.0_real64/sqrt(2.0_real64)
        call normalize_unit_length(vector, 2, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_huge_entries: [huge, huge] ierr")
        call assert_equal_array_real(vector, expected, 2, 2*TOL, &
                                     "test_normalize_unit_length_huge_entries: [huge, huge] becomes [1, 1]/sqrt(2)")
    end subroutine test_normalize_unit_length_huge_entries

    !> [v, v] for normal v whose square is subnormal (below about 1.5e-154 squared falls under
    !| tiny): the result is [1, 1]/sqrt(2) and has length 1 in every build, whether or not it
    !| flushes subnormals, as the scaled squares stay normal.
    subroutine test_normalize_unit_length_subnormal_squares()
        integer(int32) :: ierr, i_value
        real(real64) :: vector(2), expected(2), values(3)
        character(len=96) :: label

        values = [1.572e-162_real64, 2.5e-162_real64, 2.722e-162_real64]
        expected = 0.7071067811865476_real64
        do i_value = 1, size(values)
            write (label, '(a, es10.4, a)') "test_normalize_unit_length_subnormal_squares: [", values(i_value), &
                ", same]"
            vector = values(i_value)
            call normalize_unit_length(vector, 2, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, trim(label)//" ierr")
            call assert_equal_array_real(vector, expected, 2, 2*TOL, trim(label)//" becomes [1, 1]/sqrt(2)")
            call assert_equal_real(sqrt(vector(1)**2 + vector(2)**2), 1.0_real64, 2*TOL, trim(label)//" has length 1")
        end do
    end subroutine test_normalize_unit_length_subnormal_squares

    !> [s, s] for the smallest subnormal s: [1, 1]/sqrt(2), where the build keeps subnormals. Where
    !| it flushes them, s reads as 0 and [s, s] is the zero vector: ERR_DIVISION_BY_ZERO.
    subroutine test_normalize_unit_length_smallest_subnormal()
        integer(int32) :: ierr
        real(real64) :: vector(2), expected(2)

        vector = subnormal_from_bits(SMALLEST_SUBNORMAL_BITS)
        call normalize_unit_length(vector, 2, ierr)
        if (flushes_subnormals()) then
            call assert_equal_int(get_err_code(ierr), ERR_DIVISION_BY_ZERO, &
                                  "test_normalize_unit_length_smallest_subnormal: flushed [s, s] is the zero vector")
        else
            expected = 0.7071067811865476_real64
            call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_smallest_subnormal: [s, s] ierr")
            call assert_equal_array_real(vector, expected, 2, 2*TOL, &
                                         "test_normalize_unit_length_smallest_subnormal: [s, s] becomes [1, 1]/sqrt(2)")
        end if
    end subroutine test_normalize_unit_length_smallest_subnormal

    !> [tiny, 0.9*tiny], whose second entry is subnormal: about [0.74329, 0.66896], of length 1 and
    !| with the entries' ratio, where the build keeps subnormals. Where it flushes them, the second
    !| entry reads as 0 and the vector becomes [1, 0].
    subroutine test_normalize_unit_length_subnormal_entry()
        integer(int32) :: ierr
        real(real64) :: vector(2), expected(2), input_ratio

        vector(1) = tiny(1.0_real64)
        vector(2) = subnormal_from_bits(NINE_TENTHS_OF_TINY_BITS)
        input_ratio = vector(2)/vector(1)
        call normalize_unit_length(vector, 2, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_normalize_unit_length_subnormal_entry: ierr")
        if (flushes_subnormals()) then
            expected = [1.0_real64, 0.0_real64]
            call assert_equal_array_real(vector, expected, 2, 0.0_real64, &
                                         "test_normalize_unit_length_subnormal_entry: flushed, it becomes [1, 0]")
        else
            call assert_equal_real(sqrt(vector(1)**2 + vector(2)**2), 1.0_real64, 2*TOL, &
                                   "test_normalize_unit_length_subnormal_entry: has length 1")
            call assert_equal_real(vector(2)/vector(1), input_ratio, 2*TOL, &
                                   "test_normalize_unit_length_subnormal_entry: keeps the ratio of its entries")
        end if
    end subroutine test_normalize_unit_length_subnormal_entry

end module mod_test_normalize_unit_length
