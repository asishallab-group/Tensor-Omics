#include <src/macros.h>

!> Descriptive statistics: percentiles, empirical distribution functions, and 2-D LOESS smoothing.
!|
!| One of the modules the `f42_utils` family gathers.
module f42_stats_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use f42_sort_impl, only: binary_search_insertion
    M_IMPLICIT_NONE

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
        integer(int32), intent(in) :: indices_used(n_used)
            !! Indices of reference points used for smoothing (only valid indices).
            !! DM_MIN(1_int32)
            !! DM_MAX(n_total)
        integer(int32), intent(in) :: n_used
            !! Number of indices actually used for smoothing.
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
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
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
    pure real(real64) function calc_percentile_rank(percentile, n) result(rank)
        integer(int32), intent(in) :: n
            !! Sample size
        real(real64), intent(in) :: percentile
            !! desired percentile as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)

        rank = percentile*real(n - 1, real64) + 1.0_real64
    end function calc_percentile_rank

    !> summary: Calculate the percentile of an array given a sorted permutation
    !| AUTHOR_AARON_SCHROEDER
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile_impl(array, n_array, array_perm, percentile, value, n_considered)
        real(real64), intent(in) :: array(n_array)
            !! input array
        integer(int32), intent(in) :: n_array
            !! number of elements in `array`
        integer(int32), intent(in) :: array_perm(n_array)
            !! Permutation of `array` in ascending order. The allocating entry point builds and
            !! heapsorts it for you; the expert one takes whatever order you supply.
            !! DM_MIN(1_int32)
            !! DM_MAX(n_array)
        real(real64), intent(in) :: percentile
            !! desired percentile as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
        real(real64), intent(out) :: value
            !! output percentile value
        integer(int32), intent(in), optional :: n_considered
            !! How many leading entries of `array_perm` the percentile is taken over, for a
            !! percentile of a subset -- the trailing entries are ignored rather than sliced
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

        index = calc_percentile_rank(percentile, n)
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
    end subroutine calc_percentile_impl


    !> summary: Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)
    !| AUTHOR_VIVIAN_BASS
    !| This is NOT a null-hypothesis-testing p-value: each distance is compared against the
    !| observed distribution it was drawn from, not an independently generated null distribution.
    !| It instead measures how extreme an observed distance is relative to all observed distances.
    !|
    !| Implements:
    !|   Q(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !|
    !| Because distances are non-negative, a one-sided upper-tail quantile is used.
    !|
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    pure subroutine compute_scaled_distance_quantile_impl(n_genes, rdi, sorted_rdi, sorted_rdi_perm, &
                                                          quantile, c_const)
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
        real(real64), intent(out) :: quantile(n_genes)
            !! Output array to store the computed quantile for each gene.
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
            quantile = 1.0_real64
            return
        end if

        do i = 1, n_genes
            d = rdi(i)

            ! Invalid / negative => not an outlier: quantile = 1
            if (d < 0.0_real64) then
                quantile(i) = 1.0_real64
                cycle
            end if

            first_ge = binary_search_insertion(sorted_rdi, sorted_rdi_perm, d)

            if (first_ge <= n_genes) then
                count_ge = n_genes - first_ge + 1_int32
            else
                count_ge = 0_int32
            end if

            quantile(i) = (real(count_ge, real64) + c_const)/denom
        end do

    end subroutine compute_scaled_distance_quantile_impl
end module f42_stats_impl
