!> The f42 reductions that sum over a set of values: `mean` and `std_dev` of `f42_math_impl`,
!| `scaled_length` and `norm` of `f42_vector_impl`, and the power-of-two scaling
!| `scaling_exponent` they share, which keeps a sum or a sum of squares from overflowing or
!| underflowing for any finite input.
module mod_test_f42_reductions
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_support_underflow_control, ieee_get_underflow_mode, &
                                             ieee_set_underflow_mode
    use f42_math_impl, only: scaling_exponent, mean, std_dev
    use f42_vector_impl, only: norm, scaled_length
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: EPS = epsilon(1.0_real64)

contains

    !> Get array of all available tests.
    function get_all_tests_f42_reductions() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(10))
        all_tests(1) = test_case("test_scaling_exponent_range_edges", test_scaling_exponent_range_edges)
        all_tests(2) = test_case("test_mean_values", test_mean_values)
        all_tests(3) = test_case("test_mean_extreme_magnitudes", test_mean_extreme_magnitudes)
        all_tests(4) = test_case("test_std_dev_centered", test_std_dev_centered)
        all_tests(5) = test_case("test_std_dev_equal_values", test_std_dev_equal_values)
        all_tests(6) = test_case("test_std_dev_extreme_magnitudes", test_std_dev_extreme_magnitudes)
        all_tests(7) = test_case("test_std_dev_close_values_flushed", test_std_dev_close_values_flushed)
        all_tests(8) = test_case("test_scaled_length_zero_vector", test_scaled_length_zero_vector)
        all_tests(9) = test_case("test_norm_values", test_norm_values)
        all_tests(10) = test_case("test_norm_extreme_magnitudes", test_norm_extreme_magnitudes)
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

    !> Means of ordinary values: [1, 2, 3, 4] gives exactly 2.5, and [0.1, 0.2, 0.3, 0.4] gives 0.25
    !| to 2 ulps, as none of those four is exact in binary.
    subroutine test_mean_values()
        real(real64) :: values(4)

        values = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64]
        call assert_equal_real(mean(values), 2.5_real64, 0.0_real64, "test_mean_values: [1, 2, 3, 4]")
        values = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        call assert_equal_real(mean(values), 0.25_real64, 2*spacing(0.25_real64), "test_mean_values: [0.1, 0.2, 0.3, 0.4]")
    end subroutine test_mean_values

    !> The sum of values near huge does not overflow the mean: three copies of 1e308 average to
    !| 1e308, three of huge to huge, each to 1 ulp. Summed plainly, both sums are Inf.
    subroutine test_mean_extreme_magnitudes()
        real(real64) :: values(3), average

        values = 1.0e308_real64
        call assert_equal_real(mean(values), 1.0e308_real64, spacing(1.0e308_real64), &
                               "test_mean_extreme_magnitudes: three copies of 1e308")

        values = huge(1.0_real64)
        average = mean(values)
        call assert_true(ieee_is_finite(average), "test_mean_extreme_magnitudes: three copies of huge, finite")
        call assert_equal_real(average, huge(1.0_real64), spacing(huge(1.0_real64)), &
                               "test_mean_extreme_magnitudes: three copies of huge")
    end subroutine test_mean_extreme_magnitudes

    !> The centered form keeps the spread of values far from zero: [b + 1, b + 2, b + 3] has the
    !| population standard deviation sqrt(2/3) and the sample one 1, for b = 1e8 and 1e10 alike,
    !| where the sum of squares minus the squared mean loses those digits to cancellation.
    !| [1e200, 2e200, 3e200], whose squares overflow, has sqrt(2/3)*1e200 and 1e200. All to 4 eps,
    !| relative.
    subroutine test_std_dev_centered()
        real(real64) :: values(3), offsets(2), population, sample
        integer(int32) :: i_offset
        character(len=64) :: label

        offsets = [1.0e8_real64, 1.0e10_real64]
        do i_offset = 1, size(offsets)
            write (label, '(a, es7.1)') "test_std_dev_centered: b = ", offsets(i_offset)
            values = offsets(i_offset) + [1.0_real64, 2.0_real64, 3.0_real64]
            population = sqrt(2.0_real64/3.0_real64)
            call assert_equal_real(std_dev(values), population, 4*EPS*population, trim(label)//", population")
            call assert_equal_real(std_dev(values, do_bessel_correction=.true._c_bool), 1.0_real64, 4*EPS, &
                                   trim(label)//", sample (Bessel)")
        end do

        values = [1.0e200_real64, 2.0e200_real64, 3.0e200_real64]
        population = sqrt(2.0_real64/3.0_real64)*1.0e200_real64
        sample = 1.0e200_real64
        call assert_equal_real(std_dev(values), population, 4*EPS*population, &
                               "test_std_dev_centered: [1e200, 2e200, 3e200], population")
        call assert_equal_real(std_dev(values, do_bessel_correction=.true._c_bool), sample, 4*EPS*sample, &
                               "test_std_dev_centered: [1e200, 2e200, 3e200], sample (Bessel)")
    end subroutine test_std_dev_centered

    !> Values that are all equal have a spread of exactly 0 in both modes, even where their computed
    !| mean is not exactly the value: ten copies of 0.1, three of 1e308.
    subroutine test_std_dev_equal_values()
        real(real64) :: tenths(10), large(3)

        tenths = 0.1_real64
        call assert_equal_real(std_dev(tenths), 0.0_real64, 0.0_real64, "test_std_dev_equal_values: ten copies of 0.1, population")
        call assert_equal_real(std_dev(tenths, do_bessel_correction=.true._c_bool), 0.0_real64, 0.0_real64, &
                               "test_std_dev_equal_values: ten copies of 0.1, sample (Bessel)")

        large = 1.0e308_real64
        call assert_equal_real(std_dev(large), 0.0_real64, 0.0_real64, "test_std_dev_equal_values: three copies of 1e308, population")
        call assert_equal_real(std_dev(large, do_bessel_correction=.true._c_bool), 0.0_real64, 0.0_real64, &
                               "test_std_dev_equal_values: three copies of 1e308, sample (Bessel)")
    end subroutine test_std_dev_equal_values

    !> [huge, -huge] has the population standard deviation huge; with Bessel's correction it is
    !| sqrt(2)*huge, past the largest real64, so Inf.
    subroutine test_std_dev_extreme_magnitudes()
        real(real64) :: values(2), sample

        values = [huge(1.0_real64), -huge(1.0_real64)]
        call assert_equal_real(std_dev(values), huge(1.0_real64), 0.0_real64, &
                               "test_std_dev_extreme_magnitudes: [huge, -huge], population")
        sample = std_dev(values, do_bessel_correction=.true._c_bool)
        call assert_true(sample > huge(1.0_real64) .and. .not. ieee_is_finite(sample), &
                         "test_std_dev_extreme_magnitudes: [huge, -huge], sample (Bessel) is past huge, Inf")
    end subroutine test_std_dev_extreme_magnitudes

    !> Distinct values close together keep their spread where results below `tiny` are flushed to
    !| zero: [2**-470, 2**-470 + 2**-520] has the population standard deviation 2**-521 exactly and
    !| the sample one 2**-520.5, and [2**-465, 2**-465 + 2**-515] has 2**-516 and 2**-515.5 (the
    !| sample ones to 4 eps, relative). Their squared distances lie far below `tiny`. With gradual
    !| underflow those squares are exact subnormals, so the spread comes out right either way; the
    !| run with underflow set to flush to zero, where the host supports switching, is the one that
    !| shows the distances are scaled before they are squared. The default mode must agree.
    subroutine test_std_dev_close_values_flushed()
        integer(int32), parameter :: N_BASES = 2, DEFAULT_MODE = 1, FLUSHED_MODE = 2
        real(real64), parameter :: bases(N_BASES) = [2.0_real64**(-470), 2.0_real64**(-465)]
        real(real64), parameter :: gaps(N_BASES) = [2.0_real64**(-520), 2.0_real64**(-515)]
        character(len=*), parameter :: base_names(N_BASES) = [character(len=7) :: "2**-470", "2**-465"]
        character(len=*), parameter :: mode_names(2) = [character(len=24) :: "default underflow mode", &
                                                        "underflow flushed to 0"]
        real(real64) :: values(2), population(N_BASES, 2), sample(N_BASES, 2)
        real(real64), volatile :: smallest_normal
        real(real64) :: half_of_smallest_normal
        logical :: can_flush, gradual_mode
        integer(int32) :: i_base, i_mode, n_modes
        character(len=96) :: label

        call measure_spreads(DEFAULT_MODE)

        ! Nothing between switching to flush-to-zero and switching back may return early, so the
        ! results are only stored here and checked after the saved mode is restored.
        n_modes = DEFAULT_MODE
        can_flush = ieee_support_underflow_control(1.0_real64)
        if (can_flush) then
            n_modes = FLUSHED_MODE
            smallest_normal = tiny(1.0_real64)
            call ieee_get_underflow_mode(gradual_mode)
            call ieee_set_underflow_mode(.false.)
            call measure_spreads(FLUSHED_MODE)
            half_of_smallest_normal = smallest_normal/2.0_real64
            call ieee_set_underflow_mode(gradual_mode)
            call assert_equal_real(half_of_smallest_normal, 0.0_real64, 0.0_real64, &
                                   "test_std_dev_close_values_flushed: flush-to-zero is in force, tiny/2 reads as 0")
        end if

        do i_mode = 1, n_modes
            do i_base = 1, N_BASES
                label = "test_std_dev_close_values_flushed: "//trim(base_names(i_base))//", "//trim(mode_names(i_mode))
                call assert_equal_real(population(i_base, i_mode), gaps(i_base)/2.0_real64, 0.0_real64, &
                                       trim(label)//", population")
                call assert_equal_real(sample(i_base, i_mode), gaps(i_base)/sqrt(2.0_real64), &
                                       4*EPS*gaps(i_base)/sqrt(2.0_real64), &
                                       trim(label)//", sample (Bessel)")
            end do
        end do

    contains

        !> Both spreads of each pair of values, in the underflow mode now in force.
        subroutine measure_spreads(i_mode_now)
            integer(int32), intent(in) :: i_mode_now
                !! Column of `population` and `sample` the results go to

            integer(int32) :: i_pair

            do i_pair = 1, N_BASES
                values(1) = bases(i_pair)
                values(2) = bases(i_pair) + gaps(i_pair)
                population(i_pair, i_mode_now) = std_dev(values)
                sample(i_pair, i_mode_now) = std_dev(values, do_bessel_correction=.true._c_bool)
            end do
        end subroutine measure_spreads
    end subroutine test_std_dev_close_values_flushed

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
