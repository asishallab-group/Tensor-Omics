#include "macros.h"

!> This module implements heaps in Fortran, these are:
!|
!| - min-heaps
!| - max-heaps
!| - top-k-heaps (using min-heap under the hood)
!| - bottom-k-heaps (using max-heap under the hood)
module f42_heaps
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_is_nan
    use tox_errors, only: ERR_EMPTY_INPUT, set_ok, set_err, is_err
    use f42_utils, only: swap_real
    implicit none

contains

    !> Initializes a new array as top-k-heap by filling it with -Inf. This ensures correct heap construction for new non-NaN values
    pure subroutine init_top_k_heap(heap, heap_size)
        integer(int32), intent(in)    :: heap_size
            !! Size of `heap`
        real(real64), intent(out) :: heap(heap_size)
            !! The top-k-heap

        heap = M_NEG_INF
    end subroutine init_top_k_heap

    !> Pushes a `value` to a top-k-heap if it is better than the current best/min (replace best and sift down to keep heap condition)
    pure subroutine top_k_heap_push(heap, heap_size, value)
        integer(int32), intent(in)    :: heap_size
            !! Size of `heap`
        real(real64), intent(inout) :: heap(heap_size)
            !! The top-k-heap array
        real(real64), intent(in)    :: value
            !! Value to push

        if (heap_size <= 0) return

        ! add value only if it is higher than the best value
        if (value > heap(1)) then
            heap(1) = value
            call minheap_sift_down(heap, heap_size)
        end if
    end subroutine top_k_heap_push

    !> Pushes a value to a min-heap if the heap is not full
    pure subroutine minheap_push(heap, max_heap_size, n, value)
        integer(int32), intent(in)    :: max_heap_size
            !! Maximum size of `heap`
        real(real64), intent(inout) :: heap(max_heap_size)
            !! The min-heap array
        integer(int32), intent(inout) :: n
            !! Number of elements already pushed to `heap`, `0 <= n <= max_heap_size`
        real(real64), intent(in)    :: value
            !! Value to push

        integer :: i, parent

        ! If heap is not empty, insert at the end
        if (n < max_heap_size .and. max_heap_size > 0 .and. .not. ieee_is_nan(value)) then
            n = n + 1
            i = n

            ! Bubble-up using hole method
            do while (i > 1)
                parent = i / 2

                ! If the new value belongs above the parent, move parent down
                if (value < heap(parent)) then
                    heap(i) = heap(parent)
                    i = parent
                else
                    exit
                end if
            end do

            ! Final position is guaranteed to be within bounds
            heap(i) = value
        end if
    end subroutine minheap_push

    !> Pops a value from a min-heap if the heap is not empty, else error
    pure subroutine minheap_pop(heap, max_heap_size, n, value, ierr)
        integer(int32), intent(in) :: max_heap_size
            !! Maxmimum size of `heap`
        integer(int32), intent(inout) :: n
            !! Number of elements already pushed to `heap`, `0 <= n <= max_heap_size`
        real(real64), intent(inout) :: heap(max_heap_size)
            !! The min-heap array
        real(real64), intent(out)   :: value
            !! Popped best value from heap. Unset on error
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        if (n <= 0) call set_err(ierr, ERR_EMPTY_INPUT, arg_pos=3_int32)
        if (max_heap_size <= 0) call set_err(ierr, ERR_EMPTY_INPUT, arg_pos=2_int32)
        if (is_err(ierr)) return

        value = heap(1)
        heap(1) = heap(n)
        n = n - 1

        call minheap_sift_down(heap, n)
    end subroutine minheap_pop

    !> The top-down sifting for a heap -> starts from root and swaps the root element with its children until heap condition is fulfilled
    pure subroutine minheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
            !! Current heap size
        real(real64), intent(inout) :: heap(n)
            !! The min-heap array

        integer(int32) :: i, left, right, smallest

        i = 1
        do
            smallest = 2*i
            right = smallest + 1

            if (smallest > n) exit

            if (right <= n) then
                if (heap(right) < heap(smallest)) smallest = right
            end if

            if (heap(smallest) < heap(i)) then
                call swap_real(heap(i), heap(smallest))
                i = smallest
            else
                exit
            end if
        end do
    end subroutine minheap_sift_down

    !> Initializes a new array as bottom-k-heap by filling it with +Inf. This ensures correct heap construction for new non-NaN values
    pure subroutine init_bottom_k_heap(heap, heap_size)
        integer(int32), intent(in)    :: heap_size
            !! Size of `heap`
        real(real64), intent(out) :: heap(heap_size)
            !! The bottom-k-heap array

        heap = M_POS_INF
    end subroutine init_bottom_k_heap

    !> Pushes a `value` to a bottom-k-heap if it is better than the current best/min (replace best and sift down to keep heap condition)
    pure subroutine bottom_k_heap_push(heap, heap_size, value)
        integer(int32), intent(in)    :: heap_size
            !! Size of `heap`
        real(real64), intent(inout) :: heap(heap_size)
            !! The bottom-k-heap array
        real(real64), intent(in)    :: value
            !! Value to push

        if (heap_size <= 0) return

        ! add value only if it is lower than the best value
        if (value < heap(1)) then
            heap(1) = value
            call maxheap_sift_down(heap, heap_size)
        end if
    end subroutine bottom_k_heap_push

    !> Pushes a value to a max-heap if the heap is not full
    pure subroutine maxheap_push(heap, max_heap_size, n, value)
        integer(int32), intent(in)    :: max_heap_size
            !! Maximum size of `heap`
        real(real64), intent(inout) :: heap(max_heap_size)
            !! The max-heap array
        integer(int32), intent(inout) :: n
            !! Number of elements already pushed to `heap`, `0 <= n <= max_heap_size`
        real(real64), intent(in)    :: value
            !! Value to push

        integer :: i, parent

        ! If heap is not empty, insert at the end
        if (n < max_heap_size .and. max_heap_size > 0 .and. .not. ieee_is_nan(value)) then
            n = n + 1
            i = n

            ! Bubble-up using hole method
            do while (i > 1)
                parent = i / 2

                ! If the new value belongs above the parent, move parent down
                if (value > heap(parent)) then
                    heap(i) = heap(parent)
                    i = parent
                else
                    exit
                end if
            end do

            ! Final position is guaranteed to be within bounds
            heap(i) = value
        end if
    end subroutine maxheap_push

    !> Pops a value from a max-heap if the heap is not empty, else error
    pure subroutine maxheap_pop(heap, max_heap_size, n, value, ierr)
        integer(int32), intent(in) :: max_heap_size
            !! Maxmimum size of `heap`
        integer(int32), intent(inout) :: n
            !! Number of elements already pushed to `heap`, `0 <= n <= max_heap_size`
        real(real64), intent(inout) :: heap(max_heap_size)
            !! The max-heap array
        real(real64), intent(out)   :: value
            !! Popped best value from heap. Unset on error
        integer(int32), intent(out) :: ierr
            !! Error code

        call set_ok(ierr)
        if (n <= 0) call set_err(ierr, ERR_EMPTY_INPUT, arg_pos=3_int32)
        if (max_heap_size <= 0) call set_err(ierr, ERR_EMPTY_INPUT, arg_pos=2_int32)
        if (is_err(ierr)) return

        value = heap(1)
        heap(1) = heap(n)
        n = n - 1

        call maxheap_sift_down(heap, n)
    end subroutine maxheap_pop

    !> The top-down sifting for a heap -> starts from root and swaps the root element with its children until heap condition is fulfilled
    pure subroutine maxheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
            !! Current heap size
        real(real64), intent(inout) :: heap(n)
            !! The max-heap array

        integer(int32) :: i, left, right, largest

        i = 1
        do
            largest = 2*i
            right = largest + 1

            if (largest > n) exit

            if (right <= n) then
                if (heap(right) > heap(largest)) largest = right
            end if

            if (heap(largest) > heap(i)) then
                call swap_real(heap(i), heap(largest))
                i = largest
            else
                exit
            end if
        end do
    end subroutine maxheap_sift_down
end module f42_heaps