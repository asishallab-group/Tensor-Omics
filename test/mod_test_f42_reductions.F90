!> The f42 reductions that sum over a set of values: `scaled_length` and `norm` of
!| `f42_vector_impl`, and the power-of-two scaling `scaling_exponent` of `f42_math_impl` they share,
!| which keeps a sum of squares from overflowing or underflowing for any finite input.
module mod_test_f42_reductions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use f42_math_impl, only: scaling_exponent
    use f42_vector_impl, only: norm, scaled_length
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: EPS = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests_f42_reductions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(4))
        all_tests(1) = test_case("test_scaling_exponent_range_edges", test_scaling_exponent_range_edges)
        all_tests(2) = test_case("test_scaled_length_zero_vector", test_scaled_length_zero_vector)
        all_tests(3) = test_case("test_norm_values", test_norm_values)
        all_tests(4) = test_case("test_norm_extreme_magnitudes", test_norm_extreme_magnitudes)
    end function get_all_tests_f42_reductions

    !> Values whose largest magnitude lies in [2**-470, 2**496] are not scaled (exponent 0); one ulp
    !| above the range they are scaled down by 2**600, one ulp below it up by 2**600. Checked on
    !| `scaling_exponent` and on the exponent `scaled_length` reports, which comes from it.
    subroutine test_scaling_exponent_range_edges()
        real(real64), parameter :: LOWER_EDGE = 2.0_real64**(-470), UPPER_EDGE = 2.0_real64**496
        real(real64) :: edges(4), vector(2), scaled_norm
        integer(int32) :: expected(4), exponent, i_edge
        character(len=*), parameter :: names(4) = [character(len=21) :: "2**-470", "one ulp below 2**-470", &
                                                    "2**496", "one ulp above 2**496"]

        edges = [LOWER_EDGE, nearest(LOWER_EDGE, -1.0_real64), UPPER_EDGE, nearest(UPPER_EDGE, 1.0_real64)]
        expected = [0_int32, 600_int32, 0_int32, -600_int32]
        do i_edge = 1, size(edges)
            ! the largest magnitude decides, whatever its sign and whatever the smaller entries
            vector(1) = 0.5_real64*edges(i_edge)
            vector(2) = -edges(i_edge)
            call assert_equal_int(scaling_exponent(vector), expected(i_edge), &
                                  "test_scaling_exponent_range_edges: scaling_exponent, largest "//trim(names(i_edge)))
            call scaled_length(2, vector, exponent, scaled_norm)
            call assert_equal_int(exponent, expected(i_edge), &
                                  "test_scaling_exponent_range_edges: scaled_length, largest "//trim(names(i_edge)))
        end do
    end subroutine test_scaling_exponent_range_edges

    !> The zero vector, and only it, has a scaled length of exactly 0; the smallest normal vector
    !| still has a positive one.
    subroutine test_scaled_length_zero_vector()
        real(real64) :: vector(3), scaled_norm
        integer(int32) :: exponent

        vector = 0.0_real64
        call scaled_length(3, vector, exponent, scaled_norm)
        call assert_equal_real(scaled_norm, 0.0_real64, 0.0_real64, "test_scaled_length_zero_vector: the zero vector")

        vector = [tiny(1.0_real64), 0.0_real64, 0.0_real64]
        call scaled_length(3, vector, exponent, scaled_norm)
        call assert_true(scaled_norm > 0.0_real64, "test_scaled_length_zero_vector: [tiny, 0, 0] has a positive length")
        call assert_equal_real(scale(scaled_norm, -exponent), tiny(1.0_real64), 0.0_real64, &
                               "test_scaled_length_zero_vector: [tiny, 0, 0] has length tiny")
    end subroutine test_scaled_length_zero_vector

    !> Lengths of ordinary vectors: [3, 4] has length 5, [1, 2, 2] length 3, and the zero vector 0.
    subroutine test_norm_values()
        real(real64) :: pair(2), triple(3), zeros(4)

        pair = [3.0_real64, 4.0_real64]
        call assert_equal_real(norm(pair), 5.0_real64, 0.0_real64, "test_norm_values: [3, 4]")
        triple = [1.0_real64, -2.0_real64, 2.0_real64]
        call assert_equal_real(norm(triple), 3.0_real64, 0.0_real64, "test_norm_values: [1, -2, 2]")
        zeros = 0.0_real64
        call assert_equal_real(norm(zeros), 0.0_real64, 0.0_real64, "test_norm_values: the zero vector")
    end subroutine test_norm_values

    !> Squares past the real64 range do not overflow the norm: four entries of 1e300 have norm
    !| 2e300. The norm is Inf only where the true length exceeds huge, as for [huge, huge], whose
    !| length is sqrt(2)*huge.
    subroutine test_norm_extreme_magnitudes()
        real(real64) :: large(4), huge_pair(2), length

        large = 1.0e300_real64
        call assert_equal_real(norm(large), 2.0e300_real64, 2*EPS*2.0e300_real64, &
                               "test_norm_extreme_magnitudes: four entries of 1e300")

        huge_pair = huge(1.0_real64)
        length = norm(huge_pair)
        call assert_true(length > huge(1.0_real64) .and. .not. ieee_is_finite(length), &
                         "test_norm_extreme_magnitudes: [huge, huge] has a length past huge, Inf")
    end subroutine test_norm_extreme_magnitudes

end module mod_test_f42_reductions
