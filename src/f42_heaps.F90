#include "macros.h"

module f42_heaps
    use safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, ieee_negative_inf
    use tox_errors, only: ERR_EMPTY_INPUT, set_ok, set_err, is_err
    use f42_utils, only: swap_real
    implicit none

contains

    pure subroutine init_top_k_heap(heap, max_heap_size)
        integer(int32), intent(in)    :: max_heap_size
        real(real64), intent(out) :: heap(max_heap_size)

        heap = M_NEG_INF
    end subroutine init_top_k_heap

    pure subroutine top_k_heap_push(heap, max_heap_size, value)
        integer(int32), intent(in)    :: max_heap_size
        real(real64), intent(inout) :: heap(max_heap_size)
        real(real64), intent(in)    :: value

        if (max_heap_size <= 0) return

        ! add value only if it is higher than the best value
        if (value > heap(1)) then
            heap(1) = value
            call minheap_sift_down(heap, max_heap_size)
        end if
    end subroutine top_k_heap_push

    pure subroutine minheap_push(heap, max_heap_size, n, value)
        integer(int32), intent(in)    :: max_heap_size
        real(real64), intent(inout) :: heap(max_heap_size)
        integer(int32), intent(inout) :: n
        real(real64), intent(in)    :: value
        integer :: i, parent

        ! If heap is not empty, insert at the end
        if (n < max_heap_size .and. max_heap_size > 0) then
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

    pure subroutine minheap_pop(heap, max_heap_size, n, value, ierr)
        integer(int32), intent(in) :: max_heap_size
        integer(int32), intent(inout) :: n
        real(real64), intent(inout) :: heap(max_heap_size)
        real(real64), intent(out)   :: value
        integer(int32), intent(out) :: ierr

        call set_ok(ierr)
        if (n <= 0) call set_err(ierr, ERR_EMPTY_INPUT)
        if (max_heap_size <= 0) call set_err(ierr, ERR_EMPTY_INPUT)
        if (is_err(ierr)) return

        value = heap(1)
        heap(1) = heap(n)
        n = n - 1

        call minheap_sift_down(heap, n)
    end subroutine minheap_pop

    pure subroutine minheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
        real(real64), intent(inout) :: heap(n)

        integer(int32) :: i, left, right, smallest

        i = 1
        do
            smallest = 2*i
            right = smallest + 1

            if (smallest > n) exit

            if (right <= n) then
                if (heap(right) < heap(smallest)) smallest = right
            end if

            call swap_real(heap(i), heap(smallest))
            i = smallest
        end do
    end subroutine minheap_sift_down

    pure subroutine init_bottom_k_heap(heap, max_heap_size)
        integer(int32), intent(in)    :: max_heap_size
        real(real64), intent(out) :: heap(max_heap_size)

        heap = M_POS_INF
    end subroutine init_bottom_k_heap

    pure subroutine bottom_k_heap_push(heap, max_heap_size, value)
        integer(int32), intent(in)    :: max_heap_size
        real(real64), intent(inout) :: heap(max_heap_size)
        real(real64), intent(in)    :: value

        if (max_heap_size <= 0) return

        ! add value only if it is lower than the best value
        if (value < heap(1)) then
            heap(1) = value
            call maxheap_sift_down(heap, max_heap_size)
        end if
    end subroutine bottom_k_heap_push

    pure subroutine maxheap_push(heap, max_heap_size, n, value)
        integer(int32), intent(in)    :: max_heap_size
        real(real64), intent(inout) :: heap(max_heap_size)
        integer(int32), intent(inout) :: n
        real(real64), intent(in)    :: value
        integer :: i, parent

        ! If heap is not empty, insert at the end
        if (n < max_heap_size .and. max_heap_size > 0) then
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
        else
            ! If heap is full, add value only if it is worse than the best value
            n = max_heap_size
            if (value < heap(1)) then
                heap(1) = value
                call maxheap_sift_down(heap, n)
            end if
        end if
    end subroutine maxheap_push

    pure subroutine maxheap_pop(heap, max_heap_size, n, value, ierr)
        integer(int32), intent(in) :: max_heap_size
        integer(int32), intent(inout) :: n
        real(real64), intent(inout) :: heap(max_heap_size)
        real(real64), intent(out)   :: value
        integer(int32), intent(out) :: ierr

        call set_ok(ierr)
        if (n <= 0) call set_err(ierr, ERR_EMPTY_INPUT)
        if (max_heap_size <= 0) call set_err(ierr, ERR_EMPTY_INPUT)
        if (is_err(ierr)) return

        value = heap(1)
        heap(1) = heap(n)
        n = n - 1

        call maxheap_sift_down(heap, n)
    end subroutine maxheap_pop

    pure subroutine maxheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
        real(real64), intent(inout) :: heap(n)

        integer(int32) :: i, left, right, smallest

        i = 1
        do
            smallest = 2*i
            right = smallest + 1

            if (smallest > n) exit

            if (right <= n) then
                if (heap(right) > heap(smallest)) smallest = right
            end if

            call swap_real(heap(i), heap(smallest))
            i = smallest
        end do
    end subroutine maxheap_sift_down
end module f42_heaps