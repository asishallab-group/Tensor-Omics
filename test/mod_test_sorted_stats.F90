!> Unit test suite for the sorted-slice statistics of f42_stats_impl: the median absolute
!| deviation of an already-sorted slice (calc_sorted_slice_mad) and the equal-width bin counts of
!| already-sorted values (sorted_bin_counts, by binary search, against histogram_bin_counts, by
!| clamping and binning every value).
module mod_test_sorted_stats
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use f42_stats_impl, only: calc_sorted_slice_mad, sorted_bin_counts, histogram_bin_counts, last_non_nan_position, &
                              calc_median, calc_quantile_impl
    use test_suite, only: test_case
    implicit none
    public

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available sorted-slice statistics tests.
    function get_all_tests_sorted_stats() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(8))
        all_tests(1) = test_case("test_sorted_slice_mad_odd_length", test_sorted_slice_mad_odd_length)
        all_tests(2) = test_case("test_sorted_slice_mad_even_length", test_sorted_slice_mad_even_length)
        all_tests(3) = test_case("test_sorted_slice_mad_all_equal", test_sorted_slice_mad_all_equal)
        all_tests(4) = test_case("test_sorted_slice_mad_single_element", test_sorted_slice_mad_single_element)
        all_tests(5) = test_case("test_sorted_slice_mad_offset_slice", test_sorted_slice_mad_offset_slice)
        all_tests(6) = test_case("test_sorted_bin_counts_matches_histogram_bin_counts", &
                                 test_sorted_bin_counts_matches_histogram_bin_counts)
        all_tests(7) = test_case("test_last_non_nan_position", test_last_non_nan_position)
        all_tests(8) = test_case("test_calc_median", test_calc_median)
    end function get_all_tests_sorted_stats

    !> Issue #217 (C1), odd length, and the merge's left-only branch. Sorted slice
    !| [0, 1, 2, 2.25, 2.5] (given unsorted through `perm`): rank 0.5*4+1 = 3, median 2.
    !| Deviations 2, 1, 0, 0.25, 0.5 sort to 0, 0.25, 0.5, 1, 2 -> MAD 0.5. The merge takes 0 (left),
    !| 0.25 and 0.5 (right, which is then exhausted), then 1 from the left side alone.
    subroutine test_sorted_slice_mad_odd_length()
        real(real64) :: median, mad

        call calc_sorted_slice_mad([2.25_real64, 0.0_real64, 2.5_real64, 1.0_real64, 2.0_real64], 5_int32, &
                                   [2_int32, 4_int32, 5_int32, 1_int32, 3_int32], 1_int32, 5_int32, median, mad)
        call assert_equal_real(median, 2.0_real64, TOL, "test_sorted_slice_mad_odd_length: median")
        call assert_equal_real(mad, 0.5_real64, TOL, "test_sorted_slice_mad_odd_length: mad")
    end subroutine test_sorted_slice_mad_odd_length


    !> Issue #217 (C2), even length. [0, 1, 3, 7]: rank 2.5, median 1 + 0.5*(3-1) = 2; deviations
    !| 2, 1, 1, 5 sort to 1, 1, 2, 5 -> MAD 1 + 0.5*(2-1) = 1.5. [3, 7] exercises the merge's
    !| right-only branch: median 5, the tied middle pair takes the left first, then only the right
    !| is left -> deviations 2, 2, MAD 2.
    subroutine test_sorted_slice_mad_even_length()
        real(real64) :: median, mad

        call calc_sorted_slice_mad([0.0_real64, 1.0_real64, 3.0_real64, 7.0_real64], 4_int32, &
                                   [1_int32, 2_int32, 3_int32, 4_int32], 1_int32, 4_int32, median, mad)
        call assert_equal_real(median, 2.0_real64, TOL, "test_sorted_slice_mad_even_length: median of 4")
        call assert_equal_real(mad, 1.5_real64, TOL, "test_sorted_slice_mad_even_length: mad of 4")

        call calc_sorted_slice_mad([7.0_real64, 3.0_real64], 2_int32, [2_int32, 1_int32], 1_int32, 2_int32, median, mad)
        call assert_equal_real(median, 5.0_real64, TOL, "test_sorted_slice_mad_even_length: median of 2")
        call assert_equal_real(mad, 2.0_real64, TOL, "test_sorted_slice_mad_even_length: mad of 2")
    end subroutine test_sorted_slice_mad_even_length


    !> Issue #217 (C3): four equal values have median 5 and MAD 0.
    subroutine test_sorted_slice_mad_all_equal()
        real(real64) :: median, mad

        call calc_sorted_slice_mad([5.0_real64, 5.0_real64, 5.0_real64, 5.0_real64], 4_int32, &
                                   [1_int32, 2_int32, 3_int32, 4_int32], 1_int32, 4_int32, median, mad)
        call assert_equal_real(median, 5.0_real64, TOL, "test_sorted_slice_mad_all_equal: median")
        call assert_equal_real(mad, 0.0_real64, TOL, "test_sorted_slice_mad_all_equal: mad")
    end subroutine test_sorted_slice_mad_all_equal


    !> Issue #217 (C4): a one-element slice (here position 2 of 3) is its own median, MAD 0.
    subroutine test_sorted_slice_mad_single_element()
        real(real64) :: median, mad

        call calc_sorted_slice_mad([9.0_real64, 3.5_real64, -1.0_real64], 3_int32, [3_int32, 2_int32, 1_int32], &
                                   2_int32, 2_int32, median, mad)
        call assert_equal_real(median, 3.5_real64, TOL, "test_sorted_slice_mad_single_element: median")
        call assert_equal_real(mad, 0.0_real64, TOL, "test_sorted_slice_mad_single_element: mad")
    end subroutine test_sorted_slice_mad_single_element


    !> Issue #217 (C5): only positions 2..4 of the sorted order [-1, 0, 2, 3, 100, 200] are read.
    !| Slice [0, 2, 3]: median 2, deviations 2, 0, 1 -> MAD 1 (the outer values would change both).
    subroutine test_sorted_slice_mad_offset_slice()
        real(real64) :: median, mad

        call calc_sorted_slice_mad([100.0_real64, 3.0_real64, -1.0_real64, 200.0_real64, 0.0_real64, 2.0_real64], &
                                   6_int32, [3_int32, 5_int32, 6_int32, 2_int32, 1_int32, 4_int32], 2_int32, 4_int32, &
                                   median, mad)
        call assert_equal_real(median, 2.0_real64, TOL, "test_sorted_slice_mad_offset_slice: median")
        call assert_equal_real(mad, 1.0_real64, TOL, "test_sorted_slice_mad_offset_slice: mad")
    end subroutine test_sorted_slice_mad_offset_slice


    !> `sorted_bin_counts` (binary search on sorted values, inside values only) plus the clamped
    !| tails reproduces `histogram_bin_counts` exactly, for m = 3, 17, 120 and 256. The values mix
    !| pseudo-random ones in [-2, 4] (some below R_low = -1, some above R_high = 3), every bin edge
    !| `R_low + j*w` itself (including R_low and R_high), and two NaNs. A zero-width range puts every
    !| inside value into bin 1.
    subroutine test_sorted_bin_counts_matches_histogram_bin_counts()
        character(len=*), parameter :: name = "test_sorted_bin_counts_matches_histogram_bin_counts: "
        integer(int32), parameter :: m_values(4) = [3_int32, 17_int32, 120_int32, 256_int32]
        real(real64), parameter :: range_low = -1.0_real64, range_high = 3.0_real64
        real(real64), allocatable :: values(:)
        integer(int32), allocatable :: perm(:), counts(:), hist(:)
        integer(int32) :: i_m, m, j, n, n_valid, n_below, n_inside, min_occ
        real(real64) :: width, nan_val
        character(len=8) :: m_text

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)
        do i_m = 1, size(m_values)
            m = m_values(i_m)
            write (m_text, '(i0)') m
            n = 80_int32 + (m + 1_int32) + 2_int32
            allocate (values(n), perm(n), counts(m), hist(m))
            call lcg_values(int(m, int64)*7_int64, -2.0_real64, 4.0_real64, values(1:80))
            width = (range_high - range_low)/real(m, real64)
            do j = 0, m
                values(81 + j) = range_low + real(j, real64)*width
            end do
            values(n - 1:n) = nan_val
            n_valid = n - 2_int32
            call init_perm(perm)
            call sort_array_heapsort(values, perm)

            call sorted_bin_counts(values, n, perm, n, range_low, range_high, m, counts, n_below, n_inside)
            call assert_equal_int(n_below, count(values < range_low), name//"n_below, m="//trim(m_text))
            call assert_equal_int(n_inside, count(values >= range_low .and. values <= range_high), &
                                  name//"n_inside, m="//trim(m_text))
            call assert_equal_int(sum(counts), n_inside, name//"inside counts sum to n_inside, m="//trim(m_text))
            counts(1) = counts(1) + n_below
            counts(m) = counts(m) + (n_valid - n_below - n_inside)
            call histogram_bin_counts(values, perm, n, n_valid, range_low, range_high, m, hist, min_occ)
            call assert_equal_array_int(counts, hist, m, name//"clamped counts == histogram_bin_counts, m="//trim(m_text))
            deallocate (values, perm, counts, hist)
        end do

        ! Zero-width range: every inside value (== 1.0) in bin 1, nothing elsewhere.
        allocate (values(6), perm(6), counts(5))
        values = [0.5_real64, 1.0_real64, 1.0_real64, 2.0_real64, 1.0_real64, nan_val]
        call init_perm(perm)
        call sort_array_heapsort(values, perm)
        call sorted_bin_counts(values, 6_int32, perm, 6_int32, 1.0_real64, 1.0_real64, 5_int32, counts, n_below, n_inside)
        call assert_equal_int(n_below, 1_int32, name//"zero range: one value below")
        call assert_equal_int(n_inside, 3_int32, name//"zero range: three values equal the range")
        call assert_equal_int(counts(1), 3_int32, name//"zero range: all inside values in bin 1")
        call assert_equal_int(sum(counts(2:5)), 0_int32, name//"zero range: other bins empty")
        ! histogram_bin_counts on the same zero-width range: the bin width falls back to 1, so every
        ! non-NaN value, clamped to 1.0, lands in bin 1.
        allocate (hist(5))
        call histogram_bin_counts(values, perm, 6_int32, 5_int32, 1.0_real64, 1.0_real64, 5_int32, hist, min_occ)
        call assert_equal_int(hist(1), 5_int32, name//"zero range, clamped: all five values in bin 1")
        call assert_equal_int(sum(hist(2:5)), 0_int32, name//"zero range, clamped: other bins empty")
        call assert_equal_int(min_occ, 0_int32, name//"zero range, clamped: min occupancy 0")
    end subroutine test_sorted_bin_counts_matches_histogram_bin_counts


    !> Test helper: `size(values)` pseudo-random values in `[low, high)` from a Park-Miller
    !| (minimal standard, multiplier 48271) linear congruential generator, so the fixture is the
    !| same on every compiler.
    subroutine lcg_values(seed, low, high, values)
        integer(int64), intent(in) :: seed
            !! Seed
        real(real64), intent(in) :: low
            !! Lower end
        real(real64), intent(in) :: high
            !! Upper end
        real(real64), intent(out) :: values(:)
            !! Generated values

        integer(int64), parameter :: lcg_mod = 2147483647_int64, lcg_mul = 48271_int64
        integer(int64) :: state
        integer(int32) :: i

        state = mod(seed, lcg_mod - 1_int64) + 1_int64
        do i = 1, size(values)
            state = mod(lcg_mul*state, lcg_mod)
            values(i) = low + (high - low)*real(state, real64)/real(lcg_mod, real64)
        end do
    end subroutine lcg_values
    !> last_non_nan_position on NaN-last sorted permutations: two NaN at the end of five values
    !| (heapsorted, so NaN really sorts last) give 3; no NaN gives the full length; all NaN give 0;
    !| a single non-NaN value gives 1. A NaN placed before a non-NaN position is not the tail and is
    !| not counted away: the scan stops at the first non-NaN value from the end.
    subroutine test_last_non_nan_position()
        character(len=*), parameter :: name = "test_last_non_nan_position: "
        real(real64) :: nan, values(5)
        integer(int32) :: perm(5)

        nan = ieee_value(1.0_real64, ieee_quiet_nan)
        values = [3.0_real64, nan, 1.0_real64, nan, 2.0_real64]
        call init_perm(perm)
        call sort_array_heapsort(values, perm)
        call assert_equal_int(last_non_nan_position(values, 5_int32, perm), 3_int32, name//"two NaN at the end")

        values = [3.0_real64, 5.0_real64, 1.0_real64, 4.0_real64, 2.0_real64]
        call init_perm(perm)
        call sort_array_heapsort(values, perm)
        call assert_equal_int(last_non_nan_position(values, 5_int32, perm), 5_int32, name//"no NaN")

        values = nan
        call init_perm(perm)
        call assert_equal_int(last_non_nan_position(values, 5_int32, perm), 0_int32, name//"all NaN")

        call assert_equal_int(last_non_nan_position([7.0_real64], 1_int32, [1_int32]), 1_int32, name//"one value")

        ! Not a NaN-last order: position 1 holds a NaN, positions 2 and 3 do not.
        call assert_equal_int(last_non_nan_position([nan, 1.0_real64, 2.0_real64], 3_int32, [1_int32, 2_int32, 3_int32]), &
                              3_int32, name//"stops at the first non-NaN from the end")
    end subroutine test_last_non_nan_position

    !> calc_median: odd count [5, 1, 3] -> 3; even count [4, 1, 3, 2] -> 2.5 (the two middle values
    !| interpolated); a single value is its own median; `tmp_perm` holds the ascending order on
    !| return; and the result equals calc_quantile_impl at 0.5 on the same permutation, bit for bit.
    subroutine test_calc_median()
        character(len=*), parameter :: name = "test_calc_median: "
        integer(int32) :: perm3(3), perm4(4), perm1(1)
        real(real64) :: median, reference, values4(4)

        call calc_median([5.0_real64, 1.0_real64, 3.0_real64], 3_int32, perm3, median)
        call assert_equal_real(median, 3.0_real64, TOL, name//"odd count")
        call assert_equal_int(perm3(1), 2_int32, name//"perm ascending (1)")
        call assert_equal_int(perm3(2), 3_int32, name//"perm ascending (2)")
        call assert_equal_int(perm3(3), 1_int32, name//"perm ascending (3)")

        values4 = [4.0_real64, 1.0_real64, 3.0_real64, 2.0_real64]
        call calc_median(values4, 4_int32, perm4, median)
        call assert_equal_real(median, 2.5_real64, TOL, name//"even count")
        call calc_quantile_impl(values4, 4_int32, perm4, 0.5_real64, reference)
        call assert_true(median == reference, name//"bit-identical to calc_quantile_impl at 0.5")

        call calc_median([-2.0_real64], 1_int32, perm1, median)
        call assert_equal_real(median, -2.0_real64, TOL, name//"single value")
    end subroutine test_calc_median

end module mod_test_sorted_stats
