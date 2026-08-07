#include <src/macros.h>

!> summary: Wrappers for [[f42_stats_impl(module)]]
!| Generated from the implementation; do not edit -- regenerate instead.
module f42_stats
    use f42_stats_impl, only: calc_percentile_impl, compute_edf_impl
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use f42_sort_impl, only: init_perm, sort_array_heapsort
    use tox_errors, only: set_ok, is_err, ERR_ALLOC_FAIL, set_err
    use tox_errors, only: validate_all_in_range_int, validate_all_in_range_real, validate_dimension_size, validate_in_range_int
    use tox_errors, only: validate_in_range_real
    M_IMPLICIT_NONE
    private

    public :: compute_edf
    public :: compute_edf_expert
    public :: calc_percentile
    public :: calc_percentile_expert

contains

    !> summary: Validates its inputs, prepares what [[f42_stats_impl(module):compute_edf_impl]] needs, then calls it. The entry point to reach for first; see [[f42_stats(module):compute_edf_expert]] to prepare it yourself.
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf(&
            values,&
            n_values,&
            unique_values,&
            cdf_values,&
            n_unique,&
            ierr&
        )
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        real(real64), dimension(n_values), intent(in) :: values
            !! Array of observed data values (e.g., contributions or spikes).
        real(real64), dimension(n_values), intent(out) :: unique_values
            !! Sorted unique data values.
            !! The first `n_unique` elements will hold the results.
        real(real64), dimension(n_values), intent(out) :: cdf_values
            !! Corresponding cumulative frequencies between 0 and 1.
            !! The first `n_unique` elements will hold the results.
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: values_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(values_perm(n_values))
        call init_perm(values_perm)
        call sort_array_heapsort(values, values_perm)

        call compute_edf_impl(&
            values = values,&
            n_values = n_values,&
            values_perm = values_perm,&
            unique_values = unique_values,&
            cdf_values = cdf_values,&
            n_unique = n_unique&
        )
    end subroutine compute_edf

    !> summary: Validates its inputs, then calls [[f42_stats_impl(module):compute_edf_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_stats(module):compute_edf]] does both.
    !| Returns the sorted unique values and their cumulative frequencies in [0,1].
    !| The number of unique values can be determined by finding the last non-zero cdf_value.
    pure subroutine compute_edf_expert(&
            values,&
            n_values,&
            values_perm,&
            unique_values,&
            cdf_values,&
            n_unique,&
            ierr&
        )
        integer(int32), intent(in) :: n_values
            !! Number of values in the input array.
        real(real64), dimension(n_values), intent(in) :: values
            !! Array of observed data values (e.g., contributions or spikes).
        integer(int32), dimension(n_values), intent(in) :: values_perm
            !! Permutation of `values` in ascending order. The allocating entry point builds
            !! and heapsorts it for you; the expert one takes whatever order you supply.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_values`.
        real(real64), dimension(n_values), intent(out) :: unique_values
            !! Sorted unique data values.
            !! The first `n_unique` elements will hold the results.
        real(real64), dimension(n_values), intent(out) :: cdf_values
            !! Corresponding cumulative frequencies between 0 and 1.
            !! The first `n_unique` elements will hold the results.
        integer(int32), intent(out) :: n_unique
            !! Number of unique values found (actual size of output arrays)
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_values, ierr, arg_pos=2_int32)
        call validate_all_in_range_real(values, n_values, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(values_perm, n_values, ierr, arg_pos=3_int32, min=1_int32, max=n_values)
        if (is_err(ierr)) return
#endif

        call compute_edf_impl(&
            values = values,&
            n_values = n_values,&
            values_perm = values_perm,&
            unique_values = unique_values,&
            cdf_values = cdf_values,&
            n_unique = n_unique&
        )
    end subroutine compute_edf_expert

    !> summary: Validates its inputs, prepares what [[f42_stats_impl(module):calc_percentile_impl]] needs, then calls it. The entry point to reach for first; see [[f42_stats(module):calc_percentile_expert]] to prepare it yourself.
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile(&
            array,&
            n_array,&
            percentile,&
            value,&
            n_considered,&
            ierr&
        )
        integer(int32), intent(in) :: n_array
            !! number of elements in `array`
        real(real64), dimension(n_array), intent(in) :: array
            !! input array
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
        real(real64), intent(out) :: value
            !! output percentile value
        integer(int32), intent(in), optional :: n_considered
            !! How many leading entries of `array_perm` the percentile is taken over, for a
            !! percentile of a subset -- the trailing entries are ignored rather than sliced
            !! off, so the permutation stays the shape the sort produced. Zero, the default,
            !! considers all `n_array` of them.
            !! The default value is `0_int32`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_array`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.
        integer(int32), dimension(:), allocatable :: array_perm

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_array, ierr, arg_pos=2_int32)
        call validate_in_range_real(percentile, ierr, arg_pos=3_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_int(n_considered, ierr, arg_pos=5_int32, min=0_int32, max=n_array)
        call validate_all_in_range_real(array, n_array, ierr, arg_pos=1_int32)
        if (is_err(ierr)) return
#endif

        M_ALLOCATE(array_perm(n_array))
        call init_perm(array_perm)
        call sort_array_heapsort(array, array_perm)

        call calc_percentile_impl(&
            array = array,&
            n_array = n_array,&
            array_perm = array_perm,&
            percentile = percentile,&
            value = value,&
            n_considered = n_considered&
        )
    end subroutine calc_percentile

    !> summary: Validates its inputs, then calls [[f42_stats_impl(module):calc_percentile_impl]] with what you supply. The expert entry point: it allocates nothing and prepares nothing; [[f42_stats(module):calc_percentile]] does both.
    !| Uses linear interpolation between adjacent values.
    pure subroutine calc_percentile_expert(&
            array,&
            n_array,&
            array_perm,&
            percentile,&
            value,&
            n_considered,&
            ierr&
        )
        integer(int32), intent(in) :: n_array
            !! number of elements in `array`
        real(real64), dimension(n_array), intent(in) :: array
            !! input array
        integer(int32), dimension(n_array), intent(in) :: array_perm
            !! Permutation of `array` in ascending order. The allocating entry point builds and
            !! heapsorts it for you; the expert one takes whatever order you supply.
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_array`.
        real(real64), intent(in) :: percentile
            !! desired percentile (0-100)
            !! The minimum valid value is `0.0_real64`.
            !! The maximum valid value is `100.0_real64`.
        real(real64), intent(out) :: value
            !! output percentile value
        integer(int32), intent(in), optional :: n_considered
            !! How many leading entries of `array_perm` the percentile is taken over, for a
            !! percentile of a subset -- the trailing entries are ignored rather than sliced
            !! off, so the permutation stays the shape the sort produced. Zero, the default,
            !! considers all `n_array` of them.
            !! The default value is `0_int32`.
            !! The minimum valid value is `0_int32`.
            !! The maximum valid value is `n_array`.
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_dimension_size(n_array, ierr, arg_pos=2_int32)
        call validate_in_range_real(percentile, ierr, arg_pos=4_int32, min=0.0_real64, max=100.0_real64)
        call validate_in_range_int(n_considered, ierr, arg_pos=6_int32, min=0_int32, max=n_array)
        call validate_all_in_range_real(array, n_array, ierr, arg_pos=1_int32)
        call validate_all_in_range_int(array_perm, n_array, ierr, arg_pos=3_int32, min=1_int32, max=n_array)
        if (is_err(ierr)) return
#endif

        call calc_percentile_impl(&
            array = array,&
            n_array = n_array,&
            array_perm = array_perm,&
            percentile = percentile,&
            value = value,&
            n_considered = n_considered&
        )
    end subroutine calc_percentile_expert

end module f42_stats
