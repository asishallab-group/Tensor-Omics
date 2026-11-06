#include <src/macros.h>

!> summary: Implementation fixture for the extents the bindings compute and the ones they check
!| author: A Developer
!| An output sized by an integer division that Python must truncate as Fortran does; and every
!| form of input axis the shape checks compare (#221 A-shape): a constant one, two the C wrapper
!| compares because the declaration computes them -- one reading another array's contents, one of
!| an optional array --, an optional array sharing an extent, and an input sized by a mask's
!| count. Compiled and run by `test_end_to_end_fortran.py`, and through C, Python and R by
!| `test_end_to_end_impl.py`.
module fx_extents_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    M_IMPLICIT_NONE
    private

    public :: truncated_count_impl, pair_sums_impl, shifted_sum_impl, grouped_sums_impl
    public :: weighted_total_impl, selected_total_impl, optional_shift_impl

contains

    !> summary: The numbers one to `(n - 10)/3 + 5`
    !| author: A Developer
    !| An edge case: an output extent with an integer division of an operand that may be
    !| negative, which truncates toward zero in Fortran and C and floors in Python's `//`.
    pure subroutine truncated_count_impl(n, counted)
        integer(int32), intent(in) :: n
            !! decides the length
        integer(int32), dimension((n - 10)/3 + 5), intent(out) :: counted
            !! one to its own length
        integer(int32) :: i
        counted = [(i, i = 1, size(counted))]
    end subroutine truncated_count_impl

    !> summary: The sum of each row's two values
    !| author: A Developer
    !| An edge case: an input axis the declaration fixes at a constant, which the bindings compare.
    pure subroutine pair_sums_impl(n, pairs, sums)
        integer(int32), intent(in) :: n
            !! number of rows
        real(real64), dimension(n, 2), intent(in) :: pairs
            !! two values per row
        real(real64), dimension(n), intent(out) :: sums
            !! each row's sum
        sums = pairs(:, 1) + pairs(:, 2)
    end subroutine pair_sums_impl

    !> summary: The sum of `y` and of `x`, which has one element more
    !| author: A Developer
    !| An edge case: an input axis the declaration computes (`n + 1`), whose actual extent C is
    !| given and the C wrapper compares.
    pure subroutine shifted_sum_impl(n, y, x, total)
        integer(int32), intent(in) :: n
            !! number of elements of `y`
        real(real64), dimension(n), intent(in) :: y
            !! the first values
        real(real64), dimension(n + 1), intent(in) :: x
            !! the second values, one more
        real(real64), intent(out) :: total
            !! the sum of both
        total = sum(y) + sum(x)
    end subroutine shifted_sum_impl

    !> summary: The sum of each group of values
    !| author: A Developer
    !| An edge case: an input axis computed from another array's contents, like
    !| `calc_tiss_avg`'s `expr(sum(reps_per_tissue), n_genes)`, which the C wrapper may evaluate
    !| only once that array is known not to be null.
    pure subroutine grouped_sums_impl(n_groups, group_sizes, values, sums)
        integer(int32), intent(in) :: n_groups
            !! number of groups
        integer(int32), dimension(n_groups), intent(in) :: group_sizes
            !! how many values each group has, in order
        real(real64), dimension(sum(group_sizes)), intent(in) :: values
            !! the values, group after group
        real(real64), dimension(n_groups), intent(out) :: sums
            !! each group's sum
        integer(int32) :: group, first
        first = 1
        do group = 1, n_groups
            sums(group) = sum(values(first:first + group_sizes(group) - 1))
            first = first + group_sizes(group)
        end do
    end subroutine grouped_sums_impl

    !> summary: The sum of the values, weighted when weights are given
    !| author: A Developer
    !| An edge case: an optional array sharing an extent, compared only when it is given.
    pure subroutine weighted_total_impl(n, values, total, weights)
        integer(int32), intent(in) :: n
            !! number of values
        real(real64), dimension(n), intent(in) :: values
            !! the values
        real(real64), intent(out) :: total
            !! their (weighted) sum
        real(real64), dimension(n), intent(in), optional :: weights
            !! one weight per value
        if (present(weights)) then
            total = sum(values * weights)
        else
            total = sum(values)
        end if
    end subroutine weighted_total_impl

    !> summary: The weighted sum of the selected values
    !| author: A Developer
    !| An edge case: an input sized by a mask's count, which no other input carries.
    pure subroutine selected_total_impl(n_values, values, values_mask, n_selected_values, weights, &
                                        total)
        integer(int32), intent(in) :: n_values
            !! number of values
        real(real64), dimension(n_values), intent(in) :: values
            !! the values
        logical(c_bool), dimension(n_values), intent(in) :: values_mask
            !! which values to sum
        integer(int32), intent(in) :: n_selected_values
            !! how many are selected
        real(real64), dimension(n_selected_values), intent(in) :: weights
            !! one weight per selected value
        real(real64), intent(out) :: total
            !! the weighted sum of the selected values
        total = sum(pack(values, logical(values_mask)) * weights)
    end subroutine selected_total_impl

    !> summary: The sum of `y`, plus that of `x` when given, which has one element more
    !| author: A Developer
    !| An edge case: an optional input axis the declaration computes, compared by the C wrapper
    !| only when the array is present.
    pure subroutine optional_shift_impl(n, y, total, x)
        integer(int32), intent(in) :: n
            !! number of elements of `y`
        real(real64), dimension(n), intent(in) :: y
            !! the first values
        real(real64), intent(out) :: total
            !! the sum of both
        real(real64), dimension(n + 1), intent(in), optional :: x
            !! the second values, one more
        total = sum(y)
        if (present(x)) total = total + sum(x)
    end subroutine optional_shift_impl

end module fx_extents_impl
