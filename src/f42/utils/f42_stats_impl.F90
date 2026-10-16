#include <src/macros.h>

!> Descriptive statistics: quantiles, empirical distribution functions, and 2-D LOESS smoothing.
!|
!| One of the modules the `f42_utils` family gathers.
module f42_stats_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_sort_impl, only: binary_search_insertion
    use f42_math_impl, only: clamp
    M_IMPLICIT_NONE

    ! Which monotone predicate sorted_partition_point searches for; the `PREDICATE_` prefix
    ! follows from the argument name `predicate_kind`.
    integer(int32), parameter :: PREDICATE_BELOW = 1_int32
        !! [[f42_stats_impl(module):sorted_partition_point(function)]]: `probed_value < threshold`
    integer(int32), parameter :: PREDICATE_AT_MOST = 2_int32
        !! [[f42_stats_impl(module):sorted_partition_point(function)]]: `probed_value <= threshold`
    integer(int32), parameter :: PREDICATE_BIN_AT_MOST = 3_int32
        !! [[f42_stats_impl(module):sorted_partition_point(function)]]: the equal-width bin index
        !! of `probed_value` is at most `max_bin`

contains

    !> summary: Performs LOESS smoothing on a set of data points
    !| AUTHOR_VIVIAN_BASS
    !| Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
    !| The user must pre-filter data and provide only valid indices in indices_used.
    pure subroutine loess_smooth_2d_impl(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                                         kernel_sigma, kernel_cutoff, y_out)
        integer(int32), intent(in) :: n_total
            !! Total number of reference points.
        integer(int32), intent(in) :: n_target
            !! Number of target points to smooth.
        real(real64), intent(in) :: x_ref(n_total)
            !! Reference x-coordinates.
        real(real64), intent(in) :: y_ref(n_total)
            !! Reference y-coordinates (length n_total).
        integer(int32), intent(in) :: n_used
            !! Number of indices actually used for smoothing.
        integer(int32), intent(in) :: indices_used(n_used)
            !! Indices of reference points used for smoothing (only valid indices).
            !! DM_MIN(1_int32)
            !! DM_MAX(n_total)
        real(real64), intent(in) :: x_query(n_target)
            !! Target x-coordinates to smooth.
        real(real64), intent(in) :: kernel_sigma
            !! Bandwidth parameter for the kernel.
            !! DM_MIN(0.0_real64)
        real(real64), intent(in) :: kernel_cutoff
            !! Cutoff for the kernel, not used if zero
            !! DM_MIN(0.0_real64)
        real(real64), intent(out) :: y_out(n_target)
            !! Output smoothed values (length n_target).

        integer(int32) :: i_target, i_used, used_idx
        real(real64) :: query_x, ref_x, delta, sum_weights, weight
        real(real64) :: min_dist
        integer(int32) :: min_idx
        logical(c_bool) :: exact_match_found, use_kernel

        ! Check if we should use kernel smoothing
        use_kernel = (kernel_sigma > 0.0_real64)

        do concurrent(i_target=1:n_target) &
            local(query_x, sum_weights, min_dist, min_idx, exact_match_found, used_idx, ref_x, delta, i_used, weight) &
            shared(y_out, n_used, indices_used, x_ref, x_query, y_ref, use_kernel, kernel_cutoff, kernel_sigma)

            query_x = x_query(i_target)
            sum_weights = 0.0_real64
            y_out(i_target) = 0.0_real64
            min_dist = huge(1.0_real64)
            min_idx = indices_used(1)
            exact_match_found = .false.

            ! Process all reference points
            do i_used = 1, n_used
                used_idx = indices_used(i_used)
                ref_x = x_ref(used_idx)
                delta = abs(query_x - ref_x)

                ! Check for exact match
                if (delta == 0.0_real64) then
                    y_out(i_target) = y_ref(used_idx)
                    exact_match_found = .true.
                    exit
                end if

                ! Track closest point for potential fallback
                if (delta < min_dist) then
                    min_dist = delta
                    min_idx = used_idx
                end if

                ! Apply kernel smoothing if enabled and within cutoff
                if (use_kernel .and. delta <= kernel_cutoff*kernel_sigma) then
                    weight = exp(-(delta/kernel_sigma)**2)
                    sum_weights = sum_weights + weight
                    y_out(i_target) = y_out(i_target) + weight*y_ref(used_idx)
                end if
            end do

            ! Finalize result if no exact match was found
            if (.not. exact_match_found) then
                if (sum_weights > 0.0_real64) then
                    ! We have weighted average from kernel smoothing
                    y_out(i_target) = y_out(i_target)/sum_weights
                else
                    ! Fallback: use nearest neighbor
                    y_out(i_target) = y_ref(min_idx)
                end if
            end if
        end do
    end subroutine loess_smooth_2d_impl


    !> summary: Compute the Empirical Distribution Function (EDF) from a sorted permutation
    !| AUTHOR_JITU_DABA
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf_impl(values, n_values, values_perm, unique_values, cdf_values, n_unique)
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: values_perm(n_values)
            !! Permutation of `values` in ascending order. The allocating entry point builds
            !! and heapsorts it for you; the expert one takes whatever order you supply.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_values)
        real(real64), intent(out) :: unique_values(n_values)
            !! Sorted unique data values.
            !! DM_RESULT_SIZE_IS(n_unique)
        real(real64), intent(out) :: cdf_values(n_values)
            !! Corresponding cumulative frequencies between 0 and 1.
            !! DM_RESULT_SIZE_IS(n_unique)
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)

        integer(int32) :: i_value
        real(real64) :: current_val, cumulative_count

        unique_values = 0.0_real64
        cdf_values = 0.0_real64

        ! Identify unique values and compute cumulative frequencies
        n_unique = 0
        cumulative_count = 0.0_real64
        current_val = -huge(1.0_real64)  ! Initialize to lowest possible value

        do i_value = 1, n_values
            ! Check if this is a new unique value (exact comparison)
            if (values(values_perm(i_value)) /= current_val) then
                ! New unique value found
                current_val = values(values_perm(i_value))
                n_unique = n_unique + 1

                unique_values(n_unique) = current_val
            end if

            ! Update cumulative count and set CDF value
            cumulative_count = cumulative_count + 1.0_real64
            cdf_values(n_unique) = cumulative_count/real(n_values, real64)
        end do
    end subroutine compute_edf_impl



    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculate the fractional index using linear interpolation method
    pure real(real64) function calc_quantile_rank(level, n) result(rank)
        integer(int32), intent(in) :: n
            !! Sample size
        real(real64), intent(in) :: level
            !! quantile level as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)

        rank = level*real(n - 1, real64) + 1.0_real64
    end function calc_quantile_rank

    !> summary: Calculate the quantile of an array at a given level, given a sorted permutation
    !| AUTHOR_AARON_SCHROEDER
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_quantile_impl(array, n_array, array_perm, level, value, n_considered)
        integer(int32), intent(in) :: n_array
            !! number of elements in `array`
        real(real64), intent(in) :: array(n_array)
            !! input array
        integer(int32), intent(in) :: array_perm(n_array)
            !! Permutation of `array` in ascending order. The allocating entry point builds and
            !! heapsorts it for you; the expert one takes whatever order you supply.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_array)
        real(real64), intent(in) :: level
            !! quantile level as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(out) :: value
            !! the quantile: the value a fraction `level` of the considered entries lies at or below
        integer(int32), intent(in), optional :: n_considered
            !! How many leading entries of `array_perm` the quantile is taken over, for a
            !! quantile of a subset -- the trailing entries are ignored rather than sliced
            !! off, so the permutation stays the shape the sort produced. Zero, the default,
            !! considers all `n_array` of them.
            !! DM_DEFAULT(0_int32)
            !! DM_MIN(0_int32)
            !! DM_MAX(n_array)

        integer(int32) :: n, lower_index
        real(real64) :: index, fraction, lower_value, upper_value

        n = n_array
        if (present(n_considered)) then
            if (n_considered > 0) n = n_considered
        end if

        ! Handle single element case
        if (n_array == 1) then
            value = array(1)
            return
        end if

        index = calc_quantile_rank(level, n)
        lower_index = floor(index)
        fraction = index - real(lower_index, real64)

        ! Handle edge cases for indices
        if (lower_index < 1) then
            value = array(array_perm(1))  ! Smallest value in sorted order
        else if (lower_index >= n) then
            value = array(array_perm(n))  ! Largest value in sorted order
        else
            ! Linear interpolation between adjacent values using permuted indices
            lower_value = array(array_perm(lower_index))
            upper_value = array(array_perm(lower_index + 1))
            value = lower_value + fraction*(upper_value - lower_value)
        end if
    end subroutine calc_quantile_impl

    !> AUTHOR_LASZLO_LANG
    !| Median and raw (unscaled) median absolute deviation (MAD) of the values
    !| `values(perm(first_pos:last_pos))`, which `perm` already orders ascending, in O(k) with no
    !| work buffer and no sort. Both medians follow
    !| [[f42_stats_impl(module):calc_quantile_impl(interface)]]'s convention (rank
    !| [[f42_stats_impl(module):calc_quantile_rank(function)]] at level 0.5, `0.5*(k-1) + 1`,
    !| linearly interpolated). The absolute deviations grow outward from the median on both sides,
    !| so a two-pointer merge outward visits them in ascending order and stops at the needed rank.
    !| A single value has itself as median and a MAD of 0.
    !|
    !| Precondition, not checked: `1 <= first_pos <= last_pos <= n_values`, and the slice contains
    !| no NaN.
    pure subroutine calc_sorted_slice_mad(values, n_values, perm, first_pos, last_pos, median, mad)
        integer(int32), intent(in) :: n_values
            !! Number of elements in `values`
        real(real64), intent(in) :: values(n_values)
            !! Values the slice is taken from
        integer(int32), intent(in) :: perm(n_values)
            !! Permutation ordering `values` ascending at least over `first_pos:last_pos`
        integer(int32), intent(in) :: first_pos
            !! First position of the slice in `perm`
        integer(int32), intent(in) :: last_pos
            !! Last position of the slice in `perm`
        real(real64), intent(out) :: median
            !! Median of the slice
        real(real64), intent(out) :: mad
            !! Raw median absolute deviation of the slice from `median`

        integer(int32) :: n_slice, lower_rank, left_pos, right_pos, i_rank
        real(real64) :: rank, fraction, left_dev, right_dev, dev, lower_dev

        n_slice = last_pos - first_pos + 1
        if (n_slice == 1_int32) then
            median = values(perm(first_pos))
            mad = 0.0_real64
            return
        end if

        ! calc_quantile_impl's rank at level 0.5. For n_slice >= 2 it is below n_slice,
        ! so both interpolation neighbors lie inside the slice.
        rank = calc_quantile_rank(0.5_real64, n_slice)
        lower_rank = floor(rank, kind=int32)
        fraction = rank - real(lower_rank, real64)
        left_pos = first_pos + lower_rank - 1_int32
        right_pos = left_pos + 1_int32
        median = values(perm(left_pos)) + fraction*(values(perm(right_pos)) - values(perm(left_pos)))

        ! Sequential two-pointer merge: each step depends on which side the previous one took.
        ! The deviations at ranks lower_rank and lower_rank + 1 are interpolated like the median.
        lower_dev = 0.0_real64
        dev = 0.0_real64
        do i_rank = 1_int32, lower_rank + 1_int32
            if (left_pos >= first_pos .and. right_pos <= last_pos) then
                left_dev = abs(values(perm(left_pos)) - median)
                right_dev = abs(values(perm(right_pos)) - median)
                if (left_dev <= right_dev) then
                    dev = left_dev
                    left_pos = left_pos - 1_int32
                else
                    dev = right_dev
                    right_pos = right_pos + 1_int32
                end if
            else if (left_pos >= first_pos) then
                dev = abs(values(perm(left_pos)) - median)
                left_pos = left_pos - 1_int32
            else
                dev = abs(values(perm(right_pos)) - median)
                right_pos = right_pos + 1_int32
            end if
            if (i_rank == lower_rank) lower_dev = dev
        end do
        mad = lower_dev + fraction*(dev - lower_dev)
    end subroutine calc_sorted_slice_mad


    !> summary: Calculate the empirical upper-tail probability (effect-size measure) of scaled expression distances (RDI)
    !| AUTHOR_VIVIAN_BASS
    !| This is NOT a null-hypothesis-testing p-value: each distance is compared against the
    !| observed distribution it was drawn from, not an independently generated null distribution.
    !| It instead measures how extreme an observed distance is relative to all observed distances.
    !|
    !| Implements:
    !|   T(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !|
    !| Because distances are non-negative, a one-sided upper-tail probability is used.
    !|
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    pure subroutine compute_scaled_distance_tail_probability_impl(n_genes, rdi, sorted_rdi, sorted_rdi_perm, &
                                                                  tail_probability, c_const)
        integer(int32), intent(in) :: n_genes
            !! Number of genes being processed.
        real(real64), intent(in) :: rdi(n_genes)
            !! empirical distribution D
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), intent(in) :: sorted_rdi(n_genes)
            !! empirical distribution D with non negative values
            !! DM_ALLOW_NAN
            !! DM_ALLOW_INFINITE
        real(real64), intent(out) :: tail_probability(n_genes)
            !! Output array to store the computed upper-tail probability for each gene.
        real(real64), intent(in) :: c_const
            !! Constant used in the computation, typically 1
        integer(int32), intent(in) :: sorted_rdi_perm(n_genes)
            !! Permutation of `sorted_rdi` in ascending order. The allocating entry point builds
            !! and heapsorts it for you; the expert one takes whatever order you supply.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_genes)

        integer(int32) :: i, first_ge, count_ge
        real(real64) :: denom, d

        denom = real(n_genes, real64) + c_const
        if (denom <= 0.0_real64) then
            tail_probability = 1.0_real64
            return
        end if

        do i = 1, n_genes
            d = rdi(i)

            ! Invalid / negative => not an outlier: tail probability = 1
            if (d < 0.0_real64) then
                tail_probability(i) = 1.0_real64
                cycle
            end if

            first_ge = binary_search_insertion(sorted_rdi, sorted_rdi_perm, d)

            if (first_ge <= n_genes) then
                count_ge = n_genes - first_ge + 1_int32
            else
                count_ge = 0_int32
            end if

            tail_probability(i) = (real(count_ge, real64) + c_const)/denom
        end do

    end subroutine compute_scaled_distance_tail_probability_impl

    !> AUTHOR_LASZLO_LANG
    !| Bins the non-NaN values `values(perm(1:n_pool))` into `trial_m` equal-width bins over
    !| `[range_low, range_high]`, clamping a value outside the range into the nearer edge bin
    !| (clamp, then bin index `min(trial_m, int((v - range_low)/bin_width) + 1)`). A range of zero
    !| or negative width falls back to a bin width of 1, so every value lands in one bin instead of
    !| a division by zero. [[f42_stats_impl(module):sorted_bin_counts(subroutine)]] gives the same
    !| counts for the inside values by binary search, without the clamping.
    !|
    !| Known limitation, not guarded: if the range width is not finite (values near +-1e308) or so
    !| small that the bin width underflows to 0, the bin index can leave `[1, trial_m]` and the
    !| count is written out of bounds.
    pure subroutine histogram_bin_counts(values, perm, n_values, n_pool, range_low, range_high, trial_m, &
                                         bin_counts, min_occupancy)
        integer(int32), intent(in) :: n_values
            !! Number of elements in `values`
        integer(int32), intent(in) :: n_pool
            !! Number of non-NaN values; only perm(1:n_pool) is used
        real(real64), intent(in) :: values(n_values)
            !! Values to bin
        integer(int32), intent(in) :: perm(n_values)
            !! Sorting permutation for `values`, ascending, NaN last
        real(real64), intent(in) :: range_low
            !! Lower end of the binned range
        real(real64), intent(in) :: range_high
            !! Upper end of the binned range
        integer(int32), intent(in) :: trial_m
            !! Number of equal-width bins
        integer(int32), intent(out) :: bin_counts(trial_m)
            !! Absolute count of values per bin
        integer(int32), intent(out) :: min_occupancy
            !! minval(bin_counts)

        real(real64) :: bin_width, clamped_value
        integer(int32) :: i_pool, bin_idx

        ! A degenerate range falls back to a fixed bin width instead of dividing by zero.
        if (range_high - range_low <= 0.0_real64) then
            bin_width = 1.0_real64
        else
            bin_width = (range_high - range_low)/real(trial_m, real64)
        end if

        bin_counts = 0_int32
        ! Sequential: two values can increment the same bin.
        do i_pool = 1, n_pool
            clamped_value = clamp(values(perm(i_pool)), min_val=range_low, max_val=range_high)
            bin_idx = min(trial_m, int((clamped_value - range_low)/bin_width) + 1_int32)
            bin_counts(bin_idx) = bin_counts(bin_idx) + 1_int32
        end do

        min_occupancy = minval(bin_counts)
    end subroutine histogram_bin_counts

    !> AUTHOR_LASZLO_LANG
    !| Counts `values(perm(1:n_perm))` (ascending, NaN last) into `trial_m` equal-width bins over
    !| `[range_low, range_high]` in O(trial_m log n_perm): the bin index of
    !| [[f42_stats_impl(module):histogram_bin_counts(subroutine)]] never decreases with the value,
    !| so each bin end is a binary search. Only values inside the range are binned (`n_inside` of
    !| them); values below it are counted in `n_below`, values above it and NaN nowhere. For a range
    !| of positive width the counts equal `histogram_bin_counts`' counts without its clamping of the
    !| outside values into the edge bins.
    !|
    !| Precondition, not checked: `perm(1:n_perm)` orders `values` ascending, NaN last.
    pure subroutine sorted_bin_counts(values, n_values, perm, n_perm, range_low, range_high, trial_m, bin_counts, &
                                      n_below, n_inside)
        integer(int32), intent(in) :: n_values
            !! Number of elements in `values`
        real(real64), intent(in) :: values(n_values)
            !! Values to count, NaN allowed
        integer(int32), intent(in) :: n_perm
            !! Number of positions in `perm`
        integer(int32), intent(in) :: perm(n_perm)
            !! Indices into `values`, ordering the counted values ascending, NaN last
        real(real64), intent(in) :: range_low
            !! Lower end of the binned range
        real(real64), intent(in) :: range_high
            !! Upper end of the binned range
        integer(int32), intent(in) :: trial_m
            !! Number of equal-width bins
        integer(int32), intent(out) :: bin_counts(trial_m)
            !! Number of inside values per bin
        integer(int32), intent(out) :: n_below
            !! Number of values below `range_low`
        integer(int32), intent(out) :: n_inside
            !! Number of values inside `[range_low, range_high]`, the sum of `bin_counts`

        real(real64) :: bin_width
        integer(int32) :: n_up_to_high, n_up_to_bin, n_before_bin, i_bin

        ! Same zero-range guard as histogram_bin_counts.
        if (range_high - range_low <= 0.0_real64) then
            bin_width = 1.0_real64
        else
            bin_width = (range_high - range_low)/real(trial_m, real64)
        end if

        ! A comparison with NaN is false, so both predicates stay monotone over a NaN tail.
        n_below = sorted_partition_point(values, n_values, perm, n_perm, 0_int32, n_perm, PREDICATE_BELOW, &
                                         range_low, bin_width, 0_int32)
        n_up_to_high = sorted_partition_point(values, n_values, perm, n_perm, n_below, n_perm, PREDICATE_AT_MOST, &
                                              range_high, bin_width, 0_int32)
        n_inside = n_up_to_high - n_below

        ! Sequential: each bin's search starts where the previous bin ended. The bin predicate is
        ! only evaluated inside the range, which holds no NaN.
        n_before_bin = n_below
        do i_bin = 1, trial_m - 1_int32
            n_up_to_bin = sorted_partition_point(values, n_values, perm, n_perm, n_before_bin, n_up_to_high, &
                                                 PREDICATE_BIN_AT_MOST, range_low, bin_width, i_bin)
            bin_counts(i_bin) = n_up_to_bin - n_before_bin
            n_before_bin = n_up_to_bin
        end do
        bin_counts(trial_m) = n_up_to_high - n_before_bin
    end subroutine sorted_bin_counts

    !> AUTHOR_LASZLO_LANG
    !| Binary search over positions `first_pos+1 .. last_pos` of a sorted permutation: returns the
    !| last position up to which a monotone predicate holds (`first_pos` when it holds nowhere).
    !|
    !| On `probed_value = values(perm(pos))`, `PREDICATE_BELOW` is `probed_value < threshold`,
    !| `PREDICATE_AT_MOST` is `probed_value <= threshold`, and `PREDICATE_BIN_AT_MOST` is
    !| `int((probed_value - threshold)/bin_width) + 1 <= max_bin` (`threshold` = range's lower end).
    !| Unlike [[f42_sort_impl(module):binary_search_insertion(function)]] it takes a permutation
    !| slice shorter than `values` and the bin predicate. Precondition, not checked: the predicate
    !| holds on a leading run of the positions and fails on the rest.
    pure integer(int32) function sorted_partition_point(values, n_values, perm, n_perm, first_pos, last_pos, &
                                                        predicate_kind, threshold, bin_width, max_bin) result(split_pos)
        integer(int32), intent(in) :: n_values
            !! Number of elements in `values`
        real(real64), intent(in) :: values(n_values)
            !! Values the permutation orders
        integer(int32), intent(in) :: n_perm
            !! Number of positions in `perm`
        integer(int32), intent(in) :: perm(n_perm)
            !! Indices into `values`, ascending order
        integer(int32), intent(in) :: first_pos
            !! The predicate is known to hold up to this position (0 for none)
        integer(int32), intent(in) :: last_pos
            !! Last position searched
        integer(int32), intent(in) :: predicate_kind
            !! PREDICATE_BELOW, PREDICATE_AT_MOST or PREDICATE_BIN_AT_MOST
        real(real64), intent(in) :: threshold
            !! Compared value, or the range's lower end for PREDICATE_BIN_AT_MOST
        real(real64), intent(in) :: bin_width
            !! Bin width, used by PREDICATE_BIN_AT_MOST only
        integer(int32), intent(in) :: max_bin
            !! Largest accepted bin index, used by PREDICATE_BIN_AT_MOST only

        integer(int32) :: holds_pos, fails_pos, mid_pos
        real(real64) :: probed_value
        logical(c_bool) :: holds

        ! Invariant: the predicate holds up to holds_pos and fails at fails_pos (last_pos + 1 is a
        ! virtual failing position).
        holds_pos = first_pos
        fails_pos = last_pos + 1_int32
        do while (fails_pos - holds_pos > 1_int32)
            mid_pos = holds_pos + (fails_pos - holds_pos)/2_int32
            probed_value = values(perm(mid_pos))
            select case (predicate_kind)
            case (PREDICATE_BELOW)
                holds = probed_value < threshold
            case (PREDICATE_AT_MOST)
                holds = probed_value <= threshold
            case default ! PREDICATE_BIN_AT_MOST
                holds = int((probed_value - threshold)/bin_width) + 1_int32 <= max_bin
            end select
            if (holds) then
                holds_pos = mid_pos
            else
                fails_pos = mid_pos
            end if
        end do
        split_pos = holds_pos
    end function sorted_partition_point
end module f42_stats_impl
