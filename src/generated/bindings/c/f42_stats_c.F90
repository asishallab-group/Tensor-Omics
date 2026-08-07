#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_stats(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module f42_stats_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_double, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: compute_edf_c
    public :: compute_edf_expert_c
    public :: calc_percentile_c
    public :: calc_percentile_expert_c

contains

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
            !! desired percentile (0-100)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
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
            !! desired percentile (0-100)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
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

end module f42_stats_c
#endif
