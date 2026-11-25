#include <src/macros.h>

!> summary: Implementation fixture for a value a producer computes that sizes only an input
!| author: A Developer
!| `n_cap` sizes the input `y` and nothing the bindings allocate, so the Python and R expert
!| reads it off `y`, as Fortran's and C's take it from the caller, instead of calling its
!| producer. Its plain tier keeps `n_cap` and computes it, which the R shim of that tier cannot
!| compile yet, so this is built and run through Python alone, by `test_end_to_end_impl.py`.
module fx_carrier_impl
    use, intrinsic :: iso_fortran_env, only: real64, int32
    M_IMPLICIT_NONE
    private

    public :: fx_room_for_impl, fx_capacity_total_impl

contains

    !> summary: Room for a number of items, which may be at most five
    !| author: A Developer
    pure subroutine fx_room_for_impl(n_count, room)
        integer(int32), intent(in) :: n_count
            !! the number of items
            !! DM_MAX(5_int32)
        integer(int32), intent(out) :: room
            !! one more than `n_count`
        room = n_count + 1_int32
    end subroutine fx_room_for_impl

    !> summary: The smallest value, plus the sum of an input as long as a computed capacity
    !| author: A Developer
    pure subroutine fx_capacity_total_impl(n_items, x, x_perm, n_cap, y, total)
        integer(int32), intent(in) :: n_items
            !! number of values
        real(real64), dimension(n_items), intent(in) :: x
            !! the values
        integer(int32), dimension(n_items), intent(in) :: x_perm
            !! the order that sorts `x` ascending
        integer(int32), intent(in) :: n_cap
            !! the capacity
            !! DM_OUTPUT_FROM(room, fx_room_for, fx_carrier, AUTO)
            !!
            !! | Producer input | Supplied by |
            !! |----------------|-------------|
            !! | n_count        | n_items     |
        real(real64), dimension(n_cap), intent(in) :: y
            !! one value per slot of the capacity
        real(real64), intent(out) :: total
            !! the smallest of `x`, plus the sum of `y`
        total = x(x_perm(1)) + sum(y)
    end subroutine fx_capacity_total_impl
end module fx_carrier_impl
