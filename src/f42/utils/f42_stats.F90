#include <src/macros.h>

!> Descriptive statistics: percentiles, empirical distribution functions, and 2-D LOESS smoothing.
!|
!| One of the modules [[f42_utils(module)]] gathers; `use f42_utils` reaches all of them.
module f42_stats
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_errors, only: validate_all_in_range_int, set_ok, set_err, validate_in_range_real
    use tox_errors, only: is_err, validate_in_range_int, validate_dimension_size, validate_all_in_range_real
    use tox_errors, only: validate_all_in_range_int, ERR_ALLOC_FAIL
    use f42_sort_impl, only: binary_search_insertion, init_perm, sort_array_heapsort, sort_real_heapsort
    M_IMPLICIT_NONE

contains

    !> M_EXPORT_C
    !| summary: Performs LOESS smoothing on a set of data points
    !| AUTHOR_VIVIAN_BASS
    !| Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
    !| The user must pre-filter data and provide only valid indices in indices_used.
    pure subroutine loess_smooth_2d(n_total, n_target, x_ref, y_ref, indices_used, n_used, x_query, &
                                    kernel_sigma, kernel_cutoff, y_out, ierr)
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
        integer(int32), intent(in) :: n_used
            !! Number of indices actually used for smoothing.
        real(real64), intent(in) :: x_query(n_target)
            !! Target x-coordinates to smooth.
        real(real64), intent(in) :: kernel_sigma
            !! Bandwidth parameter for the kernel.
        real(real64), intent(in) :: kernel_cutoff
            !! Cutoff for the kernel, not used if zero
        real(real64), intent(out) :: y_out(n_target)
            !! Output smoothed values (length n_target).
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i_target, i_used, used_idx
        real(real64) :: query_x, ref_x, delta, sum_weights, weight
        real(real64) :: min_dist
        integer(int32) :: min_idx
        logical :: exact_match_found, use_kernel

        ! Initialize error code
        call set_ok(ierr)

        ! Input validation
        call validate_dimension_size(n_total, ierr, arg_pos=1_int32)
        call validate_dimension_size(n_target, ierr, arg_pos=2_int32)
        call validate_dimension_size(n_used, ierr, arg_pos=6_int32)
        call validate_in_range_real(kernel_sigma, ierr, min=0.0_real64, arg_pos=8_int32)
        call validate_in_range_real(kernel_cutoff, ierr, min=0.0_real64, arg_pos=9_int32)
        call validate_all_in_range_int(indices_used, n_used, ierr, min=1_int32, max=n_total, arg_pos=5_int32)

        if (is_err(ierr)) return

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
    end subroutine loess_smooth_2d

    !> M_EXPORT_C
    !| summary: Compute the Empirical Distribution Function (EDF) from pre-sorted permutation
    !| AUTHOR_JITU_DABA
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf(values, n_values, perm, unique_values, cdf_values, n_unique, ierr)
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        integer(int32), intent(in) :: perm(n_values)
            !! Pre-sorted permutation indices (must be sorted by values[perm]).
        real(real64), intent(out) :: unique_values(n_values)
            !! Sorted unique data values.
            !! DM_RESULT_SIZE_IS(n_unique)
        real(real64), intent(out) :: cdf_values(n_values)
            !! Corresponding cumulative frequencies between 0 and 1.
            !! DM_RESULT_SIZE_IS(n_unique)
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Initialize error code and outputs
        call set_ok(ierr)

        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_int(perm, n_values, ierr, min=1_int32, max=n_values, arg_pos=3_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)

        if (is_err(ierr)) return

        call compute_edf_helper(values, n_values, perm, unique_values, cdf_values, n_unique)
    end subroutine compute_edf

    !> AUTHOR_JITU_DABA
    !| (no input validation) Compute the Empirical Distribution Function (EDF) from pre-sorted permutation.
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| Assumes `values` is already sorted by `values[perm]`. Caller controls sorting algorithm.
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf_helper(values, n_values, perm, unique_values, cdf_values, n_unique)
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        integer(int32), intent(in) :: perm(n_values)
            !! Pre-sorted permutation indices (must be sorted by values[perm]).
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
            if (values(perm(i_value)) /= current_val) then
                ! New unique value found
                current_val = values(perm(i_value))
                n_unique = n_unique + 1

                unique_values(n_unique) = current_val
            end if

            ! Update cumulative count and set CDF value
            cumulative_count = cumulative_count + 1.0_real64
            cdf_values(n_unique) = cumulative_count/real(n_values, real64)
        end do
    end subroutine compute_edf_helper

    !> M_EXPORT_C
    !| summary: Sorts the values and computes the Empirical Distribution Function (EDF)
    !| AUTHOR_JITU_DABA
    !| Allocates workspace internally and performs sorting before computing EDF.
    !| Use this for convenience; use compute_edf directly for custom sorting.
    pure subroutine compute_edf_alloc(values, n_values, unique_values, cdf_values, n_unique, ierr)
        real(real64), intent(in) :: values(n_values)
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        real(real64), intent(out) :: unique_values(n_values)
            !! Sorted unique data values.
            !! DM_RESULT_SIZE_IS(n_unique)
        real(real64), intent(out) :: cdf_values(n_values)
            !! Corresponding cumulative frequencies between 0 and 1.
            !! DM_RESULT_SIZE_IS(n_unique)
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(int32), intent(out) :: ierr
            !! Error code

        ! Local workspace arrays with explicit size
        integer(int32), dimension(:), allocatable :: perm

        call set_ok(ierr)

        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(perm(n_values))

        call init_perm(perm)

        call sort_array_heapsort(values, perm)

        ! Compute EDF with sorted permutation
        call compute_edf_helper(values, n_values, perm, unique_values, cdf_values, n_unique)
    end subroutine compute_edf_alloc

    !> AUTHOR_AARON_SCHROEDER
    !| Calculate the percentile of an array given a sorted permutation.
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile(array, permutation, percentile, value, ierr)
        real(real64), intent(in) :: array(:)
            !! input array
        integer(int32), intent(in) :: permutation(:)
            !! permutation vector representing sorted order
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)
        real(real64), intent(out) :: value
            !! output percentile value
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: n

        ! Initialize error
        call set_ok(ierr)

        ! Input validation
        n = size(permutation, kind=int32)
        call validate_dimension_size(n, ierr, arg_pos=2_int32)
        call validate_in_range_int(size(array, kind=int32), ierr, min=1_int32, max=n, arg_pos=1_int32)
        call validate_all_in_range_int(permutation, n, ierr, min=1_int32, max=size(array, kind=int32), arg_pos=2_int32)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=100.0_real64, arg_pos=3_int32)

        if (is_err(ierr)) return

        call calc_percentile_helper(array, permutation, percentile, value)
    end subroutine calc_percentile

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Calculate the fractional index using linear interpolation method
    pure real(real64) function calc_percentile_rank(percentile, n) result(rank)
        integer(int32), intent(in) :: n
            !! Sample size
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)

        rank = (percentile/100.0_real64)*real(n - 1, real64) + 1.0_real64
    end function calc_percentile_rank

    !> AUTHOR_AARON_SCHROEDER
    !| (no input validation) Calculate the percentile of an array given a sorted permutation.
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile_helper(array, permutation, percentile, value)
        real(real64), intent(in) :: array(:)
            !! input array
        integer(int32), intent(in) :: permutation(:)
            !! permutation vector representing sorted order
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)
        real(real64), intent(out) :: value
            !! output percentile value

        integer(int32) :: n, lower_index
        real(real64) :: index, fraction, lower_value, upper_value

        n = size(permutation, kind=int32)

        ! Handle single element case
        if (size(array, kind=int32) == 1) then
            value = array(1)
            return
        end if

        index = calc_percentile_rank(percentile, n)
        lower_index = floor(index)
        fraction = index - real(lower_index, real64)

        ! Handle edge cases for indices
        if (lower_index < 1) then
            value = array(permutation(1))  ! Smallest value in sorted order
        else if (lower_index >= n) then
            value = array(permutation(n))  ! Largest value in sorted order
        else
            ! Linear interpolation between adjacent values using permuted indices
            lower_value = array(permutation(lower_index))
            upper_value = array(permutation(lower_index + 1))
            value = lower_value + fraction*(upper_value - lower_value)
        end if
    end subroutine calc_percentile_helper

    !> AUTHOR_AARON_SCHROEDER
    !| Calculate the percentile of an array, allocating necessary arrays when no sorting permutation is given
    !| @note This subroutine uses quicksort internally which may cause a spike in memory usage for large arrays.
    pure subroutine calc_percentile_alloc(array, percentile, value, ierr)
        real(real64), intent(in) :: array(:)
            !! Input array
        real(real64), intent(in) :: percentile
            !! Desired percentile (0-100)
        real(real64), intent(out) :: value
            !! Output percentile value
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: n
        integer(int32), allocatable :: perm(:)

        n = size(array, kind=int32)
        ! Initialize error
        call set_ok(ierr)

        call validate_dimension_size(n, ierr, arg_pos=1_int32)
        call validate_in_range_real(percentile, ierr, min=0.0_real64, max=100.0_real64, arg_pos=2_int32)

        if (is_err(ierr)) return

        M_ALLOCATE(perm(n))

        call init_perm(perm)

        ! Sort the array indirectly
        call sort_real_heapsort(array, perm)

        ! Calculate percentile using sorted permutation
        call calc_percentile_helper(array, perm, percentile, value)
    end subroutine calc_percentile_alloc

    !> M_EXPORT_C
    !| summary: Calculate the empirical quantile (effect-size measure) of scaled expression distances (RDI)
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
    pure subroutine compute_scaled_distance_quantile(n_genes, rdi, sorted_rdi, perm, quantile, c_const)
        integer(int32), intent(in) :: n_genes
            !! Number of genes being processed.
        real(real64), intent(in) :: rdi(n_genes)
            !! empirical distribution D
        real(real64), intent(in) :: sorted_rdi(n_genes)
            !! empirical distribution D with non negative values
        real(real64), intent(out) :: quantile(n_genes)
            !! Output array to store the computed quantile for each gene.
        real(real64), intent(in) :: c_const
            !! Constant used in the computation, typically 1
        integer(int32), intent(in) :: perm(n_genes)
            !! Permutation array with sorted indices for sorted_rdi

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

            first_ge = binary_search_insertion(sorted_rdi, perm, d)

            if (first_ge <= n_genes) then
                count_ge = n_genes - first_ge + 1_int32
            else
                count_ge = 0_int32
            end if

            quantile(i) = (real(count_ge, real64) + c_const)/denom
        end do

    end subroutine compute_scaled_distance_quantile
end module f42_stats
