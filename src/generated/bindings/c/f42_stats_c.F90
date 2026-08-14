#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_stats(module)]]
!| Descriptive statistics: percentiles, empirical distribution functions, and 2-D LOESS smoothing.
!|
!| One of the modules the `f42_utils` family gathers.
module f42_stats_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL, ERR_ALLOC_FAIL
    M_IMPLICIT_NONE
    private

    public :: loess_smooth_2d_c
    public :: compute_edf_c
    public :: compute_edf_expert_c
    public :: calc_percentile_c
    public :: calc_percentile_expert_c
    public :: compute_scaled_distance_quantile_c
    public :: compute_scaled_distance_quantile_expert_c

contains

    !> summary: C-wrapper for [[f42_stats(module):loess_smooth_2d(subroutine)]]
    !| Smooths `y_ref` at `x_query` using reference points `x_ref`, `y_ref`, and kernel parameters.
    !| The user must pre-filter data and provide only valid indices in indices_used.
    subroutine loess_smooth_2d_c(&
            n_total,&
            n_target,&
            x_ref,&
            y_ref,&
            indices_used,&
            n_used,&
            x_query,&
            kernel_sigma,&
            kernel_cutoff,&
            y_out,&
            ierr&
        ) bind(C, name="loess_smooth_2d_c")
        use f42_stats, only: loess_smooth_2d

        integer(c_int), intent(in), target :: n_total
            !! Total number of reference points.
        integer(c_int), intent(in), target :: n_target
            !! Number of target points to smooth.
        integer(c_int), intent(in), target :: n_used
            !! Number of indices actually used for smoothing.
        real(c_double), dimension(n_total), intent(in), target :: x_ref
            !! Reference x-coordinates.
        real(c_double), dimension(n_total), intent(in), target :: y_ref
            !! Reference y-coordinates (length n_total).
        integer(c_int), dimension(n_used), intent(in), target :: indices_used
            !! Indices of reference points used for smoothing (only valid indices).
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_total`.
        real(c_double), dimension(n_target), intent(in), target :: x_query
            !! Target x-coordinates to smooth.
        real(c_double), intent(in), target :: kernel_sigma
            !! Bandwidth parameter for the kernel.
            !! The minimum valid value is `0.0_real64`.
        real(c_double), intent(in), target :: kernel_cutoff
            !! Cutoff for the kernel, not used if zero
            !! The minimum valid value is `0.0_real64`.
        real(c_double), dimension(n_target), intent(out), target :: y_out
            !! Output smoothed values (length n_target).
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_total)
        M_CHECK_NON_NULL(n_target)
        M_CHECK_NON_NULL(n_used)
        M_CHECK_NON_NULL(kernel_sigma)
        M_CHECK_NON_NULL(kernel_cutoff)
        M_CHECK_ARRAY_NON_NULL(x_ref, n_total)
        M_CHECK_ARRAY_NON_NULL(y_ref, n_total)
        M_CHECK_ARRAY_NON_NULL(indices_used, n_used)
        M_CHECK_ARRAY_NON_NULL(x_query, n_target)
        M_CHECK_ARRAY_NON_NULL(y_out, n_target)

        call loess_smooth_2d(&
            n_total = n_total,&
            n_target = n_target,&
            x_ref = x_ref,&
            y_ref = y_ref,&
            indices_used = indices_used,&
            n_used = n_used,&
            x_query = x_query,&
            kernel_sigma = kernel_sigma,&
            kernel_cutoff = kernel_cutoff,&
            y_out = y_out,&
            ierr = ierr&
        )
    end subroutine loess_smooth_2d_c

    !> summary: C-wrapper for [[f42_stats(module):compute_edf(subroutine)]]
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    subroutine compute_edf_c(&
            values,&
            n_values,&
            unique_values,&
            cdf_values,&
            n_unique,&
            ierr&
        ) bind(C, name="compute_edf_c")
        use f42_stats, only: compute_edf

        integer(c_int), intent(in), target :: n_values
            !! Number of values in the input array.
        real(c_double), dimension(n_values), intent(in), target :: values
            !! Array of observed data values (e.g., contributions or spikes).
        real(c_double), dimension(n_values), intent(out), target :: unique_values
            !! Sorted unique data values.
            !! The first `n_unique` elements will hold the results.
        real(c_double), dimension(n_values), intent(out), target :: cdf_values
            !! Corresponding cumulative frequencies between 0 and 1.
            !! The first `n_unique` elements will hold the results.
        integer(c_int), intent(out), target :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_NON_NULL(n_unique)
        M_CHECK_ARRAY_NON_NULL(values, n_values)
        M_CHECK_ARRAY_NON_NULL(unique_values, n_values)
        M_CHECK_ARRAY_NON_NULL(cdf_values, n_values)

        call compute_edf(&
            values = values,&
            n_values = n_values,&
            unique_values = unique_values,&
            cdf_values = cdf_values,&
            n_unique = n_unique,&
            ierr = ierr&
        )
    end subroutine compute_edf_c

    !> summary: C-wrapper for [[f42_stats(module):compute_edf_expert(subroutine)]]
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    subroutine compute_edf_expert_c(&
            values,&
            n_values,&
            values_perm,&
            unique_values,&
            cdf_values,&
            n_unique,&
            ierr&
        ) bind(C, name="compute_edf_expert_c")
        use f42_stats, only: compute_edf_expert

        integer(c_int), intent(in), target :: n_values
            !! Number of values in the input array.
        real(c_double), dimension(n_values), intent(in), target :: values
            !! Array of observed data values (e.g., contributions or spikes).
        integer(c_int), dimension(n_values), intent(in), target :: values_perm
            !! Permutation of `values` in ascending order. The allocating entry point builds
            !! and heapsorts it for you; the expert one takes whatever order you supply.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_values`.
        real(c_double), dimension(n_values), intent(out), target :: unique_values
            !! Sorted unique data values.
            !! The first `n_unique` elements will hold the results.
        real(c_double), dimension(n_values), intent(out), target :: cdf_values
            !! Corresponding cumulative frequencies between 0 and 1.
            !! The first `n_unique` elements will hold the results.
        integer(c_int), intent(out), target :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_values)
        M_CHECK_NON_NULL(n_unique)
        M_CHECK_ARRAY_NON_NULL(values, n_values)
        M_CHECK_ARRAY_NON_NULL(values_perm, n_values)
        M_CHECK_ARRAY_NON_NULL(unique_values, n_values)
        M_CHECK_ARRAY_NON_NULL(cdf_values, n_values)

        call compute_edf_expert(&
            values = values,&
            n_values = n_values,&
            values_perm = values_perm,&
            unique_values = unique_values,&
            cdf_values = cdf_values,&
            n_unique = n_unique,&
            ierr = ierr&
        )
    end subroutine compute_edf_expert_c

    !> summary: C-wrapper for [[f42_stats(module):calc_percentile(subroutine)]]
    !| Uses linear interpolation between adjacent values.
    subroutine calc_percentile_c(&
            array,&
            n_array,&
            percentile,&
            value,&
            n_considered,&
            ierr&
        ) bind(C, name="calc_percentile_c")
        use f42_stats, only: calc_percentile

        integer(c_int), intent(in), target :: n_array
            !! number of elements in `array`
        real(c_double), dimension(n_array), intent(in), target :: array
            !! input array
        real(c_double), intent(in), target :: percentile
            !! desired percentile as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(out), target :: value
            !! output percentile value
        integer(c_int), intent(in), target :: n_considered
            !! How many leading entries of `array_perm` the percentile is taken over, for a
            !! percentile of a subset -- the trailing entries are ignored rather than sliced
            !! off, so the permutation stays the shape the sort produced. Zero, the default,
            !! considers all `n_array` of them.
            !! The default value is `0_int32`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_array`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_array)
        M_CHECK_NON_NULL(percentile)
        M_CHECK_NON_NULL(value)
        M_CHECK_NON_NULL(n_considered)
        M_CHECK_ARRAY_NON_NULL(array, n_array)

        call calc_percentile(&
            array = array,&
            n_array = n_array,&
            percentile = percentile,&
            value = value,&
            n_considered = n_considered,&
            ierr = ierr&
        )
    end subroutine calc_percentile_c

    !> summary: C-wrapper for [[f42_stats(module):calc_percentile_expert(subroutine)]]
    !| Uses linear interpolation between adjacent values.
    subroutine calc_percentile_expert_c(&
            array,&
            n_array,&
            array_perm,&
            percentile,&
            value,&
            n_considered,&
            ierr&
        ) bind(C, name="calc_percentile_expert_c")
        use f42_stats, only: calc_percentile_expert

        integer(c_int), intent(in), target :: n_array
            !! number of elements in `array`
        real(c_double), dimension(n_array), intent(in), target :: array
            !! input array
        integer(c_int), dimension(n_array), intent(in), target :: array_perm
            !! Permutation of `array` in ascending order. The allocating entry point builds and
            !! heapsorts it for you; the expert one takes whatever order you supply.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_array`.
        real(c_double), intent(in), target :: percentile
            !! desired percentile as a fraction in [0,1] (e.g. 0.95 for the 95th percentile)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `1.0_real64`.
        real(c_double), intent(out), target :: value
            !! output percentile value
        integer(c_int), intent(in), target :: n_considered
            !! How many leading entries of `array_perm` the percentile is taken over, for a
            !! percentile of a subset -- the trailing entries are ignored rather than sliced
            !! off, so the permutation stays the shape the sort produced. Zero, the default,
            !! considers all `n_array` of them.
            !! The default value is `0_int32`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_array`.
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_array)
        M_CHECK_NON_NULL(percentile)
        M_CHECK_NON_NULL(value)
        M_CHECK_NON_NULL(n_considered)
        M_CHECK_ARRAY_NON_NULL(array, n_array)
        M_CHECK_ARRAY_NON_NULL(array_perm, n_array)

        call calc_percentile_expert(&
            array = array,&
            n_array = n_array,&
            array_perm = array_perm,&
            percentile = percentile,&
            value = value,&
            n_considered = n_considered,&
            ierr = ierr&
        )
    end subroutine calc_percentile_expert_c

    !> summary: C-wrapper for [[f42_stats(module):compute_scaled_distance_quantile(subroutine)]]
    !| This is NOT a null-hypothesis-testing p-value: each distance is compared against the
    !| observed distribution it was drawn from, not an independently generated null distribution.
    !| It instead measures how extreme an observed distance is relative to all observed distances.
    !|
    !| Implements:
    !| Q(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !|
    !| Because distances are non-negative, a one-sided upper-tail quantile is used.
    !|
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    subroutine compute_scaled_distance_quantile_c(&
            n_genes,&
            rdi,&
            sorted_rdi,&
            quantile,&
            c_const,&
            ierr&
        ) bind(C, name="compute_scaled_distance_quantile_c")
        use f42_stats, only: compute_scaled_distance_quantile

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes being processed.
        real(c_double), dimension(n_genes), intent(in), target :: rdi
            !! empirical distribution D
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(in), target :: sorted_rdi
            !! empirical distribution D with non negative values
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(out), target :: quantile
            !! Output array to store the computed quantile for each gene.
        real(c_double), intent(in), target :: c_const
            !! Constant used in the computation, typically 1
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(c_const)
        M_CHECK_ARRAY_NON_NULL(rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(quantile, n_genes)

        call compute_scaled_distance_quantile(&
            n_genes = n_genes,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            quantile = quantile,&
            c_const = c_const,&
            ierr = ierr&
        )
    end subroutine compute_scaled_distance_quantile_c

    !> summary: C-wrapper for [[f42_stats(module):compute_scaled_distance_quantile_expert(subroutine)]]
    !| This is NOT a null-hypothesis-testing p-value: each distance is compared against the
    !| observed distribution it was drawn from, not an independently generated null distribution.
    !| It instead measures how extreme an observed distance is relative to all observed distances.
    !|
    !| Implements:
    !| Q(d) = ( #{di in D | di >= d} + c ) / ( |D| + c )
    !|
    !| Because distances are non-negative, a one-sided upper-tail quantile is used.
    !|
    !| Assumptions / preconditions:
    !| - sorted_rdi(1:n_genes) contains the empirical distribution D.
    !| - If invalid RDIs exist (negative), they should already be mapped to 0 in the distribution
    subroutine compute_scaled_distance_quantile_expert_c(&
            n_genes,&
            rdi,&
            sorted_rdi,&
            sorted_rdi_perm,&
            quantile,&
            c_const,&
            ierr&
        ) bind(C, name="compute_scaled_distance_quantile_expert_c")
        use f42_stats, only: compute_scaled_distance_quantile_expert

        integer(c_int), intent(in), target :: n_genes
            !! Number of genes being processed.
        real(c_double), dimension(n_genes), intent(in), target :: rdi
            !! empirical distribution D
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        real(c_double), dimension(n_genes), intent(in), target :: sorted_rdi
            !! empirical distribution D with non negative values
            !! NaN is permitted for this value.
            !! Infinite values are permitted for this value.
        integer(c_int), dimension(n_genes), intent(in), target :: sorted_rdi_perm
            !! Permutation of `sorted_rdi` in ascending order. The allocating entry point builds
            !! and heapsorts it for you; the expert one takes whatever order you supply.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_genes`.
        real(c_double), dimension(n_genes), intent(out), target :: quantile
            !! Output array to store the computed quantile for each gene.
        real(c_double), intent(in), target :: c_const
            !! Constant used in the computation, typically 1
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_genes)
        M_CHECK_NON_NULL(c_const)
        M_CHECK_ARRAY_NON_NULL(rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi, n_genes)
        M_CHECK_ARRAY_NON_NULL(sorted_rdi_perm, n_genes)
        M_CHECK_ARRAY_NON_NULL(quantile, n_genes)

        call compute_scaled_distance_quantile_expert(&
            n_genes = n_genes,&
            rdi = rdi,&
            sorted_rdi = sorted_rdi,&
            sorted_rdi_perm = sorted_rdi_perm,&
            quantile = quantile,&
            c_const = c_const,&
            ierr = ierr&
        )
    end subroutine compute_scaled_distance_quantile_expert_c

end module f42_stats_c
#endif
