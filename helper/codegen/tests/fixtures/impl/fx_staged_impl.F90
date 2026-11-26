#include <src/macros.h>

!> summary: Implementation fixture for bounds naming values set after the top of the wrapper
!| author: A Developer
!| A bound or an equality naming a value a producer computes is checked after the
!| producers, and one naming a value the prologue writes after the prologue; what the caller's
!| values decide is still checked at the top, so the prologue never reads an argument unchecked
!| (rule V). Compiled and run by `test_end_to_end_fortran.py`.
module fx_staged_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: fx_handled_impl, fx_handled_guard, fx_index_impl, fx_index_guard
    public :: fx_two_points_impl, fx_two_points_guard, fx_feeds_impl, fx_refine_guard, fx_refined_impl

contains

    !> summary: Derives a level, and answers a single item itself
    !| author: A Developer
    pure subroutine fx_handled_guard(n_items, x, y, level, total, handled, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        real(real64), dimension(n_items), intent(in) :: x
            !! values
        real(real64), intent(in) :: y
            !! the offset
        real(real64), intent(out) :: level
            !! the level
        real(real64), intent(out) :: total
            !! the answer for a single item
        logical, intent(out) :: handled
            !! a single item
        integer(int32), intent(inout) :: ierr
            !! Error code, left as it is
        level = maxval(x)
        handled = n_items == 1
        total = 0.0_real64
        if (handled) total = x(1) + y
        if (ierr /= 0) return
    end subroutine fx_handled_guard

    !> summary: A bound naming a value the prologue writes, where the prologue may answer the call
    !| author: A Developer
    !| DM_PROLOGUE(fx_handled_guard, fx_staged_impl)
    pure subroutine fx_handled_impl(n_items, x, level, y, total, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        real(real64), dimension(n_items), intent(in) :: x
            !! values
        real(real64), intent(in) :: level
            !! a level
        real(real64), intent(in) :: y
            !! at most the level
            !! DM_MAX(level)
        real(real64), intent(out) :: total
            !! the result
        integer(int32), intent(out) :: ierr
            !! Error code
        ierr = 0
        total = sum(x) + y + 0.0_real64*level
    end subroutine fx_handled_impl

    !> summary: Derives the largest index, and reads the value at k
    !| author: A Developer
    pure subroutine fx_index_guard(n_items, x, k, level, total, handled, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        real(real64), dimension(n_items), intent(in) :: x
            !! values
        integer(int32), intent(in) :: k
            !! an index
        integer(int32), intent(out) :: level
            !! the largest index
        real(real64), intent(out) :: total
            !! zero when the picked value is zero
        logical, intent(out) :: handled
            !! the picked value is zero
        integer(int32), intent(inout) :: ierr
            !! Error code, left as it is
        level = n_items
        handled = x(k) == 0.0_real64
        total = 0.0_real64
        if (ierr /= 0) return
    end subroutine fx_index_guard

    !> summary: An index bounded by a value the prologue writes, which the prologue reads
    !| author: A Developer
    !| DM_PROLOGUE(fx_index_guard, fx_staged_impl)
    pure subroutine fx_index_impl(n_items, x, level, k, total, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        real(real64), dimension(n_items), intent(in) :: x
            !! values
        integer(int32), intent(in) :: level
            !! the largest index
        integer(int32), intent(in) :: k
            !! an index
            !! DM_MIN(1_int32)
            !! DM_MAX(level)
        real(real64), intent(out) :: total
            !! the picked value
        integer(int32), intent(out) :: ierr
            !! Error code
        ierr = 0
        total = x(k)
    end subroutine fx_index_impl

    !> summary: Derives a level
    !| author: A Developer
    pure subroutine fx_two_points_guard(n_items, level, handled, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        integer(int32), intent(out) :: level
            !! the level
        logical, intent(out) :: handled
            !! never
        integer(int32), intent(inout) :: ierr
            !! Error code, left as it is
        level = 100*n_items
        handled = .false.
        if (ierr /= 0) return
    end subroutine fx_two_points_guard

    !> summary: Bounds naming a value a producer computes and one the prologue writes
    !| author: A Developer
    !| DM_PROLOGUE(fx_two_points_guard, fx_staged_impl)
    pure subroutine fx_two_points_impl(n_items, factor, limit, level, y, w, m, total, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        integer(int32), intent(in) :: factor
            !! factor
        integer(int32), intent(in) :: limit
            !! a limit
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
        integer(int32), intent(in) :: level
            !! a level
        integer(int32), intent(in) :: y
            !! between the limit and the level
            !! DM_MIN(limit)
            !! DM_MAX(level)
        integer(int32), intent(in) :: w
            !! at least the limit
            !! DM_MIN(limit)
        integer(int32), intent(in) :: m
            !! the limit itself
            !! DM_EQUALS(limit)
        integer(int32), intent(out) :: total
            !! the result
        integer(int32), intent(out) :: ierr
            !! Error code
        ierr = 0
        total = y + w + m + 0*limit + 0*level + 0*factor + 0*n_items
    end subroutine fx_two_points_impl

    !> summary: A bound naming a value a producer computes, on an input the producer reads
    !| author: A Developer
    pure subroutine fx_feeds_impl(n_items, factor, limit, total)
        integer(int32), intent(in) :: n_items
            !! items
        integer(int32), intent(in) :: factor
            !! factor, at least 2
            !! DM_MIN(2_int32)
            !! DM_MAX(limit)
        integer(int32), intent(in) :: limit
            !! a limit
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
        integer(int32), intent(out) :: total
            !! the result
        total = factor + limit + 0*n_items
    end subroutine fx_feeds_impl

    !> summary: Halves the limit a producer computed
    !| author: A Developer
    pure subroutine fx_refine_guard(limit, handled, ierr)
        integer(int32), intent(inout) :: limit
            !! halved
        logical, intent(out) :: handled
            !! never
        integer(int32), intent(inout) :: ierr
            !! Error code, left as it is
        limit = limit/2_int32
        handled = .false.
        if (ierr /= 0) return
    end subroutine fx_refine_guard

    !> summary: A prologue refining a value a producer computes, and a bound naming it
    !| author: A Developer
    !| The prologue refines `limit` (`intent(inout)`), which it may (item 5); `y`'s bound names
    !| it, so `foo` checks it after the prologue, against the refined value (rule V).
    !| DM_PROLOGUE(fx_refine_guard, fx_staged_impl)
    pure subroutine fx_refined_impl(n_items, factor, limit, y, total, ierr)
        integer(int32), intent(in) :: n_items
            !! items
        integer(int32), intent(in) :: factor
            !! factor
        integer(int32), intent(in) :: limit
            !! a limit, halved by the prologue
            !! DM_OUTPUT_FROM(scaled_size, fx_scaled_size, fx_producers, AUTO)
        integer(int32), intent(in) :: y
            !! at most the limit
            !! DM_MAX(limit)
        integer(int32), intent(out) :: total
            !! the limit, plus `y`
        integer(int32), intent(out) :: ierr
            !! Error code
        ierr = 0
        total = limit + y + 0*n_items + 0*factor
    end subroutine fx_refined_impl
end module fx_staged_impl
