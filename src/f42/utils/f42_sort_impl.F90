#include <src/macros.h>

!> Indirect sorting -- every routine reorders a permutation vector rather than the data -- plus the searches over a sorted array.
!|
!| One of the modules [[f42_utils_impl(module)]] gathers; `use f42_utils_impl` reaches all of them.
module f42_sort_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_errors, only: set_ok, is_err, validate_in_range_int, validate_dimension_size
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use f42_math_impl, only: clamp
    M_IMPLICIT_NONE

    !> Generic indirect (permutation-based) sort, dispatches on the array's element type.
    !| See [[f42_sort_impl(module):sort_real(subroutine)]], [[f42_sort_impl(module):sort_integer(subroutine)]],
    !| [[f42_sort_impl(module):sort_character(subroutine)]] for the implementations. Uses quicksort.
    interface sort_array
        module procedure sort_real, sort_integer, sort_character
    end interface sort_array

    !> Generic indirect (permutation-based) heapsort, dispatches on the array's element type.
    !| Prefer this over [[f42_sort_impl(module):sort_array(interface)]] when a worst-case O(n log n)
    !| guarantee is needed (quicksort's manual-stack partitioning can degrade on adversarial input),
    !| or when the temporary quicksort stack workspace should be avoided.
    interface sort_array_heapsort
        module procedure sort_real_heapsort, sort_integer_heapsort, sort_character_heapsort
        module procedure sort_real_heapsort_expl_size
    end interface sort_array_heapsort

    !> NaN-aware `<` comparison operator for `real(real64)`; see [[f42_sort_impl(module):real_less(function)]].
    interface operator(.lessthan.)
        module procedure real_less
    end interface operator(.lessthan.)

    !> NaN-aware `>` comparison operator for `real(real64)`; see [[f42_sort_impl(module):real_greater(function)]].
    interface operator(.greaterthan.)
        module procedure real_greater
    end interface operator(.greaterthan.)

contains

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Initializes a permutation vector with ascending indices
    pure subroutine init_perm(perm)
        integer(int32), dimension(:), intent(out) :: perm
            !! Permutation vector to initialize with indices

        integer(int32) :: i_perm

        do concurrent(i_perm=1:size(perm, kind=int32)) shared(perm)
            perm(i_perm) = i_perm
        end do
    end subroutine init_perm

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Function to find the position to place a value in a sorted array using binary search
    pure integer(int32) function binary_search_insertion(arr, perm, value, lower_idx, upper_idx) result(idx)
        real(real64), dimension(:), contiguous, intent(in) :: arr
            !! Array of values
        integer(int32), dimension(size(arr, kind=int32)), intent(in) :: perm
            !! Sorting permutation of `arr`, NaNs sorted to the end
        real(real64), intent(in) :: value
            !! Value to find
        integer(int32), intent(in), optional :: lower_idx
            !! The lower index to start searching in the array, default `1` -> searching in `perm(lower_idx:upper_idx)`
        integer(int32), intent(in), optional :: upper_idx
            !! The upper index to stop searching in the array, default `size(arr)` -> searching in `perm(lower_idx:upper_idx)`

        integer(int32) :: actual_lower_idx, actual_upper_idx, mid_idx, n

        n = size(arr, kind=int32)

        ! Guard against an empty array: there is no element to compare against, so
        ! report the (only sensible) insertion position of 1.
        if (n <= 0) then
            idx = 1
            return
        end if

        ! NaN is sorted to the end -> If value and last element is NaN, return n, else it is not found
        if (ieee_is_nan(value)) then
            if (ieee_is_nan(arr(perm(n)))) then
                idx = n
            else
                idx = n + 1
            end if
            return
        end if

        M_DEFAULT_VAL(lower_idx, actual_lower_idx, 1_int32)
        actual_lower_idx = clamp(actual_lower_idx, min_val=1_int32, max_val=n)
        M_DEFAULT_VAL(upper_idx, actual_upper_idx, n + 1)
        actual_upper_idx = clamp(actual_upper_idx, min_val=actual_lower_idx - 1, max_val=n) + 1

        do while (actual_lower_idx < actual_upper_idx)
            mid_idx = (actual_lower_idx + actual_upper_idx)/2

            ! Update bounds: Note that NaN is sorted to the end -> if condition fails because of NaN, the upper bound will be reduced
            if (arr(perm(mid_idx)) < value) then
                actual_lower_idx = mid_idx + 1
            else
                actual_upper_idx = mid_idx
            end if
        end do

        idx = actual_lower_idx
    end function binary_search_insertion

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Function to find a value in a sorted array using binary search. Returns -1 if not found
    pure integer(int32) function binary_search(arr, perm, value) result(idx)
        real(real64), dimension(:), contiguous, intent(in) :: arr
            !! Array of values
        integer(int32), dimension(size(arr, kind=int32)), intent(in) :: perm
            !! Sorting permutation of `arr`, NaNs sorted to the end
        real(real64), intent(in) :: value
            !! Value to find

        integer(int32) :: n
        real(real64) :: found

        n = size(arr, kind=int32)

        ! Nothing found in empty arrays
        if (n <= 0) then
            idx = -1
        else
            idx = binary_search_insertion(arr, perm, value)

            ! nothing found if inserted as new element after last one
            if (idx == n + 1) then
                idx = -1
            else
                found = arr(perm(idx))

                ! allow NaN search -> if value doesn't match found one and both aren't NaN, it is not found
                if (found /= value .and. .not. (ieee_is_nan(found) .and. ieee_is_nan(value))) then
                    idx = -1
                end if
            end if
        end if
    end function binary_search

    !> AUTHOR_LASZLO_LANG
    !| The first block `s` with `block_ends(s) >= idx`, by binary search over the ascending block
    !| ends: the contiguous block holding index `idx`, for an array split into consecutive blocks
    !| `block_ends(s-1)+1 .. block_ends(s)` (empty blocks allowed). Unlike
    !| [[f42_sort_impl(module):binary_search_insertion(function)]] it searches integers directly,
    !| without a permutation. Precondition, not checked: `n_blocks >= 1`, `block_ends` ascends and
    !| `1 <= idx <= block_ends(n_blocks)`; outside it the result is `n_blocks` for an `idx` above
    !| every block end and `1` for `n_blocks = 0`.
    pure integer(int32) function block_containing_index(block_ends, n_blocks, idx) result(block_idx)
        integer(int32), intent(in) :: n_blocks
            !! Number of blocks
        integer(int32), intent(in) :: block_ends(n_blocks)
            !! Last index of each block, ascending
        integer(int32), intent(in) :: idx
            !! The index to locate

        integer(int32) :: low_block, high_block, mid_block

        ! Invariant: the answer lies in low_block..high_block.
        low_block = 1_int32
        high_block = n_blocks
        do while (low_block < high_block)
            mid_block = low_block + (high_block - low_block)/2_int32
            if (block_ends(mid_block) >= idx) then
                high_block = mid_block
            else
                low_block = mid_block + 1_int32
            end if
        end do
        block_idx = low_block
    end function block_containing_index

    !> AUTHOR_VIVIAN_BASS
    !| Sort a real array indirectly using quicksort.
    !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
    pure subroutine sort_real(array, perm, tmp_stack_left, tmp_stack_right)
        real(real64), intent(in) :: array(:)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(:)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(:)
            !! Manual stack of right indices for quicksort recursion
        call quicksort_real(array, perm, int(size(array), int32), tmp_stack_left, tmp_stack_right)
    end subroutine sort_real

    !> AUTHOR_VIVIAN_BASS
    !| Sort an integer array indirectly using quicksort.
    !| Similar to `sort_real`, but for integer input.
    pure subroutine sort_integer(array, perm, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: array(:)
            !! Integer input array to sort
        integer(int32), intent(inout) :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(:)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(:)
            !! Manual stack of right indices for quicksort recursion
        call quicksort_int(array, perm, int(size(array), int32), tmp_stack_left, tmp_stack_right)
    end subroutine sort_integer

    !> AUTHOR_VIVIAN_BASS
    !| Sort a character array indirectly using quicksort.
    !| Uses lexicographic ordering and permutation vector sorting.
    pure subroutine sort_character(array, perm, tmp_stack_left, tmp_stack_right)
        character(len=*), intent(in) :: array(:)
            !! Character input array to sort
        integer(int32), intent(inout) :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(:)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(:)
            !! Manual stack of right indices for quicksort recursion
        call quicksort_char(array, perm, int(size(array), int32), tmp_stack_left, tmp_stack_right)
    end subroutine sort_character

    !> AUTHOR_MOHAMED_AKDI
    !| Sort a real array indirectly using heapsort.
    !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
    pure subroutine sort_real_heapsort(array, perm)
        real(real64), intent(in), contiguous :: array(:)
            !! Real input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        call heapsort_real(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_real_heapsort

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Sort a real explicit-size array indirectly using heapsort.
    !| Creates a sorted version of the array by reordering the `perm` vector. The original data in `array` remains unchanged.
    pure subroutine sort_real_heapsort_expl_size(array, perm, n)
        integer(int32), intent(in) :: n
            !! Size of `array`
        real(real64), intent(in) :: array(n)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        call heapsort_real(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_real_heapsort_expl_size

    !> M_EXPORT_C
    !| summary: Sort a real vector, returning the permutation (ascending, NaN last)
    !| AUTHOR_LASZLO_LANG
    pure subroutine sort_real_get_perm(array, n, perm, ierr)
        integer(int32), intent(in) :: n
            !! elements of `array`
        real(real64), intent(in) :: array(n)
            !! values to sort
        integer(int32), intent(out) :: perm(n)
            !! permutation that sorts `array` ascending, NaN last
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i

        call set_ok(ierr)
        perm = [(i, i = 1, n)]
        call sort_real_heapsort_expl_size(array, perm, n)
    end subroutine sort_real_get_perm

    !> AUTHOR_MOHAMED_AKDI
    !| Sort an integer array indirectly using heapsort.
    !| Similar to `sort_real_heapsort`, but for integer input.
    pure subroutine sort_integer_heapsort(array, perm)
        integer(int32), intent(in), contiguous :: array(:)
            !! Integer input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        call heapsort_integer(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_integer_heapsort

    !> AUTHOR_MOHAMED_AKDI
    !| Sort a character array indirectly using heapsort.
    !| Uses lexicographic ordering and permutation vector sorting.
    pure subroutine sort_character_heapsort(array, perm)
        character(len=*), intent(in), contiguous :: array(:)
            !! Character input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        call heapsort_character(array, perm, size(array, kind=int32), size(perm, kind=int32))
    end subroutine sort_character_heapsort

    !> AUTHOR_VIVIAN_BASS
    !| Internal quicksort implementation for real arrays.
    !| Sorts indirectly using the permutation vector `perm`. Manual stack replaces recursion.
    pure subroutine quicksort_real(array, perm, n, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: n
            !! Size of the array
        real(real64), intent(in) :: array(n)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(n)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(n)
            !! Manual stack of right indices for quicksort recursion

        integer(int32) :: left, right, i, j, top, pivot_idx
        real(real64) :: pivot_val

        if (n == 0) return

        top = 1
        tmp_stack_left(top) = 1
        tmp_stack_right(top) = n

        ! Iterative quicksort using explicit stack
        do while (top > 0)
            left = tmp_stack_left(top)
            right = tmp_stack_right(top)
            top = top - 1

            if (left >= right) cycle

            ! Select pivot and initialize pointers
            pivot_idx = (left + right)/2
            pivot_val = array(perm(pivot_idx))
            i = left
            j = right

            ! Partitioning loop
            do
                do while (array(perm(i)) .lessthan.pivot_val)
                    i = i + 1
                end do
                do while (array(perm(j)) .greaterthan.pivot_val)
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_int(perm(i), perm(j))
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            ! Push new ranges onto stack
            if (left < j) then
                top = top + 1
                tmp_stack_left(top) = left
                tmp_stack_right(top) = j
            end if
            if (i < right) then
                top = top + 1
                tmp_stack_left(top) = i
                tmp_stack_right(top) = right
            end if
        end do
    end subroutine quicksort_real

    !> AUTHOR_VIVIAN_BASS
    !| Internal quicksort implementation for integer arrays.
    !| Indirectly sorts `array` using `perm`, same algorithm as `quicksort_real`.
    pure subroutine quicksort_int(array, perm, n, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: n
            !! Size of the array
        integer(int32), intent(in) :: array(n)
            !! Integer input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(n)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(n)
            !! Manual stack of right indices for quicksort recursion

        integer(int32) :: left, right, i, j, top, pivot_idx
        integer(int32) :: pivot_val

        if (n == 0) return

        top = 1
        tmp_stack_left(top) = 1
        tmp_stack_right(top) = n

        do while (top > 0)
            left = tmp_stack_left(top)
            right = tmp_stack_right(top)
            top = top - 1

            if (left >= right) cycle

            pivot_idx = (left + right)/2
            pivot_val = array(perm(pivot_idx))
            i = left
            j = right

            do
                do while (array(perm(i)) < pivot_val)
                    i = i + 1
                end do
                do while (array(perm(j)) > pivot_val)
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_int(perm(i), perm(j))
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            if (left < j) then
                top = top + 1
                tmp_stack_left(top) = left
                tmp_stack_right(top) = j
            end if
            if (i < right) then
                top = top + 1
                tmp_stack_left(top) = i
                tmp_stack_right(top) = right
            end if
        end do
    end subroutine quicksort_int

    !> AUTHOR_VIVIAN_BASS
    !| Internal quicksort implementation for character arrays.
    !| Lexicographic quicksort using string comparison, indirect via `perm`.
    pure subroutine quicksort_char(array, perm, n, tmp_stack_left, tmp_stack_right)
        integer(int32), intent(in) :: n
            !! Size of the array
        character(len=*), intent(in) :: array(n)
            !! Character input array to sort
        integer(int32), intent(inout) :: perm(n)
            !! Permutation vector that will be sorted
        integer(int32), intent(out) :: tmp_stack_left(n)
            !! Manual stack of left indices for quicksort recursion
        integer(int32), intent(out) :: tmp_stack_right(n)
            !! Manual stack of right indices for quicksort recursion
        integer(int32) :: left, right, i, j, top, pivot_idx
            !! Temporary variables
        character(len=len(array)) :: pivot_val

        if (n == 0) return

        top = 1
        tmp_stack_left(top) = 1
        tmp_stack_right(top) = n

        do while (top > 0)
            left = tmp_stack_left(top)
            right = tmp_stack_right(top)
            top = top - 1

            if (left >= right) cycle

            pivot_idx = (left + right)/2
            pivot_val = array(perm(pivot_idx))
            i = left
            j = right

            do
                do while (array(perm(i)) < pivot_val)
                    i = i + 1
                end do
                do while (array(perm(j)) > pivot_val)
                    j = j - 1
                end do
                if (i <= j) then
                    call swap_int(perm(i), perm(j))
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do

            if (left < j) then
                top = top + 1
                tmp_stack_left(top) = left
                tmp_stack_right(top) = j
            end if
            if (i < right) then
                top = top + 1
                tmp_stack_left(top) = i
                tmp_stack_right(top) = right
            end if
        end do
    end subroutine quicksort_char

    !Internal heapsort implementations for real arrays.
    !TODO future: heapify_real/heapify_integer/heapify_character (and heapsort_real/_integer/_character below) are
    !             near-identical ~30-40 line copies differing only in the element type and comparison operator.
    !             Given real already uses the `.greaterthan.` operator interface for NaN-aware comparisons, the
    !             int/character variants could likely share one templated implementation (or at least the
    !             sift-down logic factored out), reducing triplicated maintenance surface.
    !             According to the current 2026 draft (https://j3-fortran.org/doc/year/26/26-007.pdf, section 16.3),
    !             there is a template structure planned to create e.g. templated subroutines
    !             with generic type specified at compile time.
    !             This may also be suitable for other <procedure>_<type> procedures like array serialization.
    !> AUTHOR_MOHAMED_AKDI
    !| Sorts indirectly using the permutation vector `perm`. Uses `heapify_real` to maintain heap property.
    pure subroutine heapsort_real(array, perm, n_arr, n_perm)
        integer(int32), intent(in) :: n_arr
            !! Size of `arr`
        integer(int32), intent(in) :: n_perm
            !! Size of `perm`
        real(real64), intent(in) :: array(n_arr)
            !! Real input array to sort
        integer(int32), intent(inout) :: perm(n_perm)
            !! Permutation vector that will be sorted
        integer(int32) :: n, i
        n = n_perm

        ! Build max heap
        do i = n/2, 1, -1
            call heapify_real(array, perm, n, i)
        end do
        ! Heap sort
        do i = n, 2, -1
            call swap_int(perm(1), perm(i))
            call heapify_real(array, perm, i - 1, 1)
        end do
    end subroutine heapsort_real

    !> AUTHOR_MOHAMED_AKDI
    !| Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_real(array, perm, heap_size, root)
        real(real64), intent(in), contiguous :: array(:)
            !! Real input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: heap_size
            !! Size of the heap
        integer(int32), intent(in) :: root
            !! Root index

        integer(int32) :: current, next_idx, largest_idx

        current = root

        do
            ! Use `largest_idx` directly as the left-child index: left = 2*current
            largest_idx = 2*current
            ! If there is no left child the subtree is a leaf; we're done
            if (largest_idx > heap_size) exit

            next_idx = largest_idx + 1

            ! Only compare the right child (next_idx) when it actually exists
            if (next_idx <= heap_size) then
                if (array(perm(next_idx)) .greaterthan.array(perm(largest_idx))) then
                    largest_idx = next_idx
                end if
            end if

            ! If the larger child is greater than current, swap permutation entries
            if (array(perm(largest_idx)) .greaterthan.array(perm(current))) then
                call swap_int(perm(current), perm(largest_idx))
                current = largest_idx
            else
                exit
            end if
        end do
    end subroutine heapify_real

    !> AUTHOR_MOHAMED_AKDI
    !| Internal heapsort implementation for integer arrays.
    !| Indirectly sorts `array` using `perm`, same algorithm as `heapsort_real`.
    pure subroutine heapsort_integer(array, perm, n_arr, n_perm)
        integer(int32), intent(in) :: n_arr
            !! Size of `arr`
        integer(int32), intent(in) :: n_perm
            !! Size of `perm`
        integer(int32), intent(in) :: array(n_arr)
            !! Integer input array to sort
        integer(int32), intent(inout) :: perm(n_perm)
            !! Permutation vector that will be sorted

        integer(int32) :: n, i

        n = n_perm

        ! Build max-heap
        do i = n/2, 1, -1
            call heapify_integer(array, perm, n, i)
        end do

        ! Heap sort
        do i = n, 2, -1
            call swap_int(perm(1), perm(i))
            call heapify_integer(array, perm, i - 1, 1)
        end do
    end subroutine heapsort_integer

    !> AUTHOR_MOHAMED_AKDI
    !| Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_integer(array, perm, heap_size, root)
        integer(int32), intent(in), contiguous :: array(:)
            !! Integer input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: heap_size
            !! Size of the heap
        integer(int32), intent(in) :: root
            !! Root index

        integer(int32) :: current, next_idx, largest_idx

        current = root

        do
            ! Compute left-child index (use largest_idx as left to avoid extra var)
            largest_idx = 2*current
            ! If there is no left child the subtree is a leaf; nothing to do
            if (largest_idx > heap_size) exit

            ! Compute right-child index (may be out of heap bounds)
            next_idx = largest_idx + 1

            ! Only compare right-child when it actually exists (guarded access)
            if (next_idx <= heap_size) then
                if (array(perm(next_idx)) > array(perm(largest_idx))) then
                    ! Right child is larger than left child
                    largest_idx = next_idx
                end if
            end if

            ! Compare the selected child with the current node; if the child is
            ! greater, swap permutation indices so the larger value moves up the
            ! heap. We swap entries of `perm` (indices), not the array values.
            if (array(perm(largest_idx)) > array(perm(current))) then
                call swap_int(perm(current), perm(largest_idx))
                current = largest_idx
            else
                exit
            end if
        end do
    end subroutine heapify_integer

    !> AUTHOR_MOHAMED_AKDI
    !| Internal heapsort implementation for character arrays.
    !| Lexicographic heapsort using string comparison, indirect via `perm`.
    pure subroutine heapsort_character(array, perm, n_arr, n_perm)
        integer(int32), intent(in) :: n_arr
            !! Size of `arr`
        integer(int32), intent(in) :: n_perm
            !! Size of `perm`
        character(len=*), intent(in)    :: array(n_arr)
            !! Character input array to sort
        integer(int32), intent(inout) :: perm(n_perm)
            !! Permutation vector that will be sorted

        integer(int32) :: n, i

        n = n_perm

        ! Build max-heap
        do i = n/2, 1, -1
            call heapify_character(array, perm, n, i)
        end do

        ! Heap sort
        do i = n, 2, -1
            call swap_int(perm(1), perm(i))
            call heapify_character(array, perm, i - 1, 1)
        end do
    end subroutine heapsort_character

    !> AUTHOR_MOHAMED_AKDI
    !| Iterative heapify (non-recursive, pure)
    !| Restore the max-heap property for the subtree rooted at `root`.
    !| Heap layout (1-based): left child = 2*i, right child = 2*i+1.
    !| We reorder the permutation vector `perm` (indices) rather than moving
    !| array values. Guard accesses by checking child indices before indexing.
    pure subroutine heapify_character(array, perm, heap_size, root)
        character(len=*), intent(in), contiguous :: array(:)
            !! Character input array to sort
        integer(int32), intent(inout), contiguous :: perm(:)
            !! Permutation vector that will be sorted
        integer(int32), intent(in) :: heap_size
            !! Size of the heap
        integer(int32), intent(in) :: root
            !! Root index

        integer(int32) :: current, next_idx, largest_idx

        current = root

        do
            ! Compute left-child index (use largest_idx directly as left = 2*current)
            largest_idx = 2*current
            ! If there is no left child the subtree is a leaf; nothing to do
            if (largest_idx > heap_size) exit

            ! Potential right child index
            next_idx = largest_idx + 1

            ! Only compare right child when it exists to avoid OOB access
            if (next_idx <= heap_size) then
                if (array(perm(next_idx)) > array(perm(largest_idx))) then
                    largest_idx = next_idx
                end if
            end if

            ! If the selected child is larger than current, swap permutation
            ! indices so the larger element moves up. We swap entries in `perm`.
            if (array(perm(largest_idx)) > array(perm(current))) then
                call swap_int(perm(current), perm(largest_idx))
                current = largest_idx
            else
                exit
            end if
        end do
    end subroutine heapify_character

    !> AUTHOR_LASZLO_LANG
    !| Push a value onto a top-k heap: a min-heap holding the `heap_size` largest values pushed so
    !| far, its root the smallest of them. A value above the root replaces it and sifts down
    !| ([[f42_sort_impl(module):minheap_sift_down(subroutine)]]); any other value is dropped. The
    !| caller seeds every slot (with values that form a min-heap, e.g. all equal), so no
    !| initialization routine exists. A `heap_size` of 0 or less does nothing.
    pure subroutine top_k_heap_push(heap, heap_size, value)
        integer(int32), intent(in) :: heap_size
            !! Size of heap
        real(real64), intent(inout) :: heap(heap_size)
            !! The top-k heap array (min-heap of the k largest values pushed so far)
        real(real64), intent(in) :: value
            !! Value to push

        if (heap_size <= 0_int32) return

        if (value > heap(1)) then
            heap(1) = value
            call minheap_sift_down(heap, heap_size)
        end if
    end subroutine top_k_heap_push

    !> AUTHOR_LASZLO_LANG
    !| Push a value onto a bottom-k heap: a max-heap holding the `heap_size` smallest values pushed
    !| so far, its root the largest of them. A value below the root replaces it and sifts down
    !| ([[f42_sort_impl(module):maxheap_sift_down(subroutine)]]); any other value is dropped. The
    !| caller seeds every slot, as for [[f42_sort_impl(module):top_k_heap_push(subroutine)]]. A
    !| `heap_size` of 0 or less does nothing.
    pure subroutine bottom_k_heap_push(heap, heap_size, value)
        integer(int32), intent(in) :: heap_size
            !! Size of heap
        real(real64), intent(inout) :: heap(heap_size)
            !! The bottom-k heap array (max-heap of the k smallest values pushed so far)
        real(real64), intent(in) :: value
            !! Value to push

        if (heap_size <= 0_int32) return

        if (value < heap(1)) then
            heap(1) = value
            call maxheap_sift_down(heap, heap_size)
        end if
    end subroutine bottom_k_heap_push

    !> AUTHOR_LASZLO_LANG
    !| Top-down sift for a min-heap of values (not a permutation, unlike
    !| [[f42_sort_impl(module):heapify_real(subroutine)]]): swaps the root with its smaller child
    !| until the heap condition holds.
    pure subroutine minheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
            !! Current heap size
        real(real64), intent(inout) :: heap(n)
            !! The min-heap array

        integer(int32) :: i, left, right, smallest
        real(real64) :: temp

        i = 1_int32
        do
            left = 2_int32*i
            right = left + 1_int32
            if (left > n) exit

            smallest = left
            if (right <= n) then
                if (heap(right) < heap(smallest)) smallest = right
            end if

            if (heap(smallest) < heap(i)) then
                temp = heap(i)
                heap(i) = heap(smallest)
                heap(smallest) = temp
                i = smallest
            else
                exit
            end if
        end do
    end subroutine minheap_sift_down

    !> AUTHOR_LASZLO_LANG
    !| Top-down sift for a max-heap of values (not a permutation, unlike
    !| [[f42_sort_impl(module):heapify_real(subroutine)]]): swaps the root with its larger child
    !| until the heap condition holds.
    pure subroutine maxheap_sift_down(heap, n)
        integer(int32), intent(in) :: n
            !! Current heap size
        real(real64), intent(inout) :: heap(n)
            !! The max-heap array

        integer(int32) :: i, left, right, largest
        real(real64) :: temp

        i = 1_int32
        do
            left = 2_int32*i
            right = left + 1_int32
            if (left > n) exit

            largest = left
            if (right <= n) then
                if (heap(right) > heap(largest)) largest = right
            end if

            if (heap(largest) > heap(i)) then
                temp = heap(i)
                heap(i) = heap(largest)
                heap(largest) = temp
                i = largest
            else
                exit
            end if
        end do
    end subroutine maxheap_sift_down

    !> AUTHOR_VIVIAN_BASS
    !| Swap two integer values in-place.
    pure subroutine swap_int(a, b)
        integer(int32), intent(inout) :: a
            !! First integer to swap
        integer(int32), intent(inout) :: b
            !! Second integer to swap

        integer(int32) :: temp
        temp = a; a = b; b = temp
    end subroutine swap_int

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Swap two real values in-place.
    pure subroutine swap_real(a, b)
        real(real64), intent(inout) :: a
            !! First real to swap
        real(real64), intent(inout) :: b
            !! Second real to swap
        real(real64) :: temp
        temp = a; a = b; b = temp
    end subroutine swap_real

    !> AUTHOR_MOHAMED_AKDI
    !| Helper: NaN-aware comparisons for real(real64). Treats NaN as greater as every real number, while `NaN` equals `NaN`.
    pure logical function real_less(lhs, rhs)
        real(real64), intent(in) :: lhs
            !! left side of the comparison `lhs < rhs`
        real(real64), intent(in) :: rhs
            !! right side of the comparison `lhs < rhs`

        if (ieee_is_nan(lhs) .and. ieee_is_nan(rhs)) then
            real_less = .false.
        else if (ieee_is_nan(lhs)) then
            real_less = .false.
        else if (ieee_is_nan(rhs)) then
            real_less = .true.
        else
            real_less = (lhs < rhs)
        end if
    end function real_less

    !> AUTHOR_MOHAMED_AKDI
    !| Helper: NaN-aware comparisons for real(real64). Treats NaN as greater as every real number, while `NaN` equals `NaN`.
    pure logical function real_greater(lhs, rhs)
        real(real64), intent(in) :: lhs
            !! left side of the comparison `lhs > rhs`
        real(real64), intent(in) :: rhs
            !! right side of the comparison `lhs > rhs`

        if (ieee_is_nan(lhs) .and. ieee_is_nan(rhs)) then
            real_greater = .false.
        else if (ieee_is_nan(lhs)) then
            real_greater = .true.
        else if (ieee_is_nan(rhs)) then
            real_greater = .false.
        else
            real_greater = (lhs > rhs)
        end if
    end function real_greater

    !> AUTHOR_VIVIAN_BASS
    !| Finds the indices of the true values in a logical mask.
    pure subroutine which(mask, n, idx_out, m_max, m_out, ierr)
        integer(int32), intent(in) :: m_max
            !! Maximum size of `idx_out`.
        integer(int32), intent(in) :: n
            !! Size of the mask.
        logical(c_bool), intent(in) :: mask(n)
            !! Logical array of size n.
        integer(int32), intent(out) :: idx_out(m_max)
            !! Integer array to store the indices of true values.
        integer(int32), intent(out) :: m_out
            !! Actual size of `idx_out` (number of true values found).
        integer(int32), intent(out) :: ierr
            !! Error code

        integer(int32) :: i, count

        ! Initialize error code
        call set_ok(ierr)

        ! Validate inputs
        call validate_dimension_size(n, ierr, arg_pos=2_int32)
        call validate_in_range_int(m_max, ierr, min=1_int32, arg_pos=4_int32)

        if (is_err(ierr)) return

        count = 0
        idx_out = 0  ! Initialize to avoid garbage values
        do i = 1, n
            if (mask(i)) then
                count = count + 1
                if (count <= m_max) then
                    idx_out(count) = i
                end if
            end if
        end do
        m_out = count
    end subroutine which
end module f42_sort_impl
