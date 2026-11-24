#include <src/macros.h>

!> summary: Bit masks: `n_bits` flags packed 32 to an `integer(int32)` word.
!| AUTHOR_FRANZ_ERIC_SILL
!|
!| A bit mask holds the same information as a `logical(c_bool)` array in an eighth of the memory,
!| and whole-mask operations work on 32 flags at once. `misc/bench/bit_masks/` measures both, and
!| records why this module is a set of procedures on plain words rather than a derived type.
!|
!| **Layout.** Flag `i` (counting from 1) is bit `mod(i - 1, 32)` of word `(i - 1)/32 + 1`, the
!| least significant bit first, so `n_bits` flags take `ceiling(n_bits/32)` words. This is
!| the order numpy's `np.unpackbits(words.view(np.uint8), bitorder="little")` and R's
!| `intToBits(words)` decode.
!|
!| **All 32 bits are used.** A word whose bit 31 is set is negative, so a word is tested against
!| zero with `/= 0`, never with `> 0`.
!|
!| **In R the word with only bit 31 set prints as `NA`.** Its bit pattern is the one R uses for
!| `NA_integer_`, and it is still that word: it passes back to any procedure here unchanged, and
!| `intToBits()` decodes it, while `bitwAnd()` and arithmetic take it for a missing value.
!|
!| **The bits past `n_bits` in the last word are always zero.** Every procedure here keeps it that
!| way, and the counts, [[f42_bit_masks_impl(module):bitmask_all(function)]] and
!| [[f42_bit_masks_impl(module):bitmask_equal(function)]] rely on it.
!|
!| **Several masks** are the columns of a `(n_words, n_masks)` array. Every column starts on a
!| new word, so `bitmasks(:, i_mask)` is a mask of its own and passes to any procedure here
!| without a copy.
!|
!| **Parallel loops write whole words.** Two iterations that set flags in the same word race, so
!| a parallel loop over flags may only read a mask; one that writes runs over words or columns.
!|
!| **Hot loops** use the whole-mask procedures, or step from one set flag to the next with
!| [[f42_bit_masks_impl(module):bitmask_next(function)]]. A procedure call per flag costs more
!| than the flag, because nothing is inlined across modules; where a loop must touch single
!| flags, the bit-mask macros in `src/macros.h` test a flag or locate its word inline, and one
!| of them sizes a mask.
module f42_bit_masks_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, int64
    use, intrinsic :: iso_c_binding, only: c_bool
    M_IMPLICIT_NONE
    private

    public :: bitmask_n_words_impl, bitmask_test_impl
    public :: bitmask_from_logical_impl, bitmask_to_logical_impl
    public :: bitmasks_from_logical_2D_impl, bitmasks_to_logical_2D_impl
    public :: bitmask_set, bitmask_clear
    public :: bitmask_fill, bitmask_copy, bitmask_and, bitmask_or, bitmask_and_not, bitmask_xor, bitmask_not
    public :: bitmask_count, bitmask_and_count, bitmask_any, bitmask_all, bitmask_any_and
    public :: bitmask_equal, bitmask_is_subset
    public :: bitmask_first, bitmask_last, bitmask_next, bitmask_indices

contains

    ! ------------------------------------------------------------------ published

    !> summary: The number of words a mask of `n_bits` flags takes.
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Inside Fortran, a bit-mask macro in `src/macros.h` gives the same number.
    pure subroutine bitmask_n_words_impl(n_bits, n_words)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask; no flags take no words
            !! DM_MIN(0_int32)
        integer(int32), intent(out) :: n_words
            !! number of words the mask takes

        n_words = M_BIT_MASK_N_WORDS(n_bits)
    end subroutine bitmask_n_words_impl

    !> summary: Whether flag `i_bit` of a mask is set.
    !| AUTHOR_FRANZ_ERIC_SILL
    !| Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.
    pure subroutine bitmask_test_impl(n_bits, n_words, bitmask, i_bit, is_set)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_bits))
            !! DM_MAX(M_BIT_MASK_N_WORDS(n_bits))
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask
        integer(int32), intent(in) :: i_bit
            !! index of the flag, counting from 1
            !! DM_MIN(1_int32)
            !! DM_MAX(n_bits)
        logical(c_bool), intent(out) :: is_set
            !! `.true.` if the flag is set

        is_set = .false._c_bool
        if (i_bit <= n_bits) is_set = M_BIT_MASK_TEST(bitmask, i_bit)
    end subroutine bitmask_test_impl

    !> summary: Packs a `logical(c_bool)` array into a bit mask.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine bitmask_from_logical_impl(n_bits, n_words, flags, bitmask)
        integer(int32), intent(in) :: n_bits
            !! number of flags
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_bits))
            !! DM_MAX(M_BIT_MASK_N_WORDS(n_bits))
        logical(c_bool), dimension(n_bits), intent(in) :: flags
            !! the flags to pack
        integer(int32), dimension(n_words), intent(out) :: bitmask
            !! the mask, flag `i` set where `flags(i)` is `.true.`

        integer(int32) :: i_word, first_flag

        do i_word = 1, n_words
            first_flag = (i_word - 1_int32)*M_BIT_MASK_WORD_BITS + 1_int32
            bitmask(i_word) = packed_word(min(M_BIT_MASK_WORD_BITS, n_bits - first_flag + 1_int32), &
                                           flags(first_flag:))
        end do
    end subroutine bitmask_from_logical_impl

    !> summary: Unpacks a bit mask into a `logical(c_bool)` array.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine bitmask_to_logical_impl(n_bits, n_words, bitmask, flags)
        integer(int32), intent(in) :: n_bits
            !! number of flags
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_bits))
            !! DM_MAX(M_BIT_MASK_N_WORDS(n_bits))
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask
        logical(c_bool), dimension(n_bits), intent(out) :: flags
            !! `flags(i)` is `.true.` where flag `i` of the mask is set

        integer(int32) :: i_bit

        do i_bit = 1, n_bits
            flags(i_bit) = M_BIT_MASK_TEST(bitmask, i_bit)
        end do
    end subroutine bitmask_to_logical_impl

    !> summary: Packs every column of a `logical(c_bool)` matrix into its own bit mask.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine bitmasks_from_logical_2D_impl(n_bits, n_masks, n_words, flags, bitmasks)
        integer(int32), intent(in) :: n_bits
            !! number of flags per mask
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_masks
            !! number of masks
        integer(int32), intent(in) :: n_words
            !! number of words per mask
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_bits))
            !! DM_MAX(M_BIT_MASK_N_WORDS(n_bits))
        logical(c_bool), dimension(n_bits, n_masks), intent(in) :: flags
            !! the flags, one mask per column
        integer(int32), dimension(n_words, n_masks), intent(out) :: bitmasks
            !! the masks, column `j` packed from `flags(:, j)`

        integer(int32) :: i_mask

        do i_mask = 1, n_masks
            call bitmask_from_logical_impl(n_bits, n_words, flags(:, i_mask), bitmasks(:, i_mask))
        end do
    end subroutine bitmasks_from_logical_2D_impl

    !> summary: Unpacks every column of a matrix of bit masks into a `logical(c_bool)` column.
    !| AUTHOR_FRANZ_ERIC_SILL
    pure subroutine bitmasks_to_logical_2D_impl(n_bits, n_masks, n_words, bitmasks, flags)
        integer(int32), intent(in) :: n_bits
            !! number of flags per mask
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_masks
            !! number of masks
        integer(int32), intent(in) :: n_words
            !! number of words per mask
            !! DM_MIN(M_BIT_MASK_N_WORDS(n_bits))
            !! DM_MAX(M_BIT_MASK_N_WORDS(n_bits))
        integer(int32), dimension(n_words, n_masks), intent(in) :: bitmasks
            !! the masks, one per column
        logical(c_bool), dimension(n_bits, n_masks), intent(out) :: flags
            !! the flags, column `j` unpacked from `bitmasks(:, j)`

        integer(int32) :: i_mask

        do i_mask = 1, n_masks
            call bitmask_to_logical_impl(n_bits, n_words, bitmasks(:, i_mask), flags(:, i_mask))
        end do
    end subroutine bitmasks_to_logical_2D_impl

    ! ------------------------------------------------------------------ single flags

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Sets flag `i_bit`. Not for a parallel loop over flags: see the module's note on parallel loops.
    pure subroutine bitmask_set(n_words, bitmask, i_bit)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the mask
        integer(int32), intent(in) :: i_bit
            !! index of the flag, `1 <= i_bit <= n_bits`

        bitmask(M_BIT_MASK_WORD(i_bit)) = ibset(bitmask(M_BIT_MASK_WORD(i_bit)), M_BIT_MASK_BIT(i_bit))
    end subroutine bitmask_set

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Clears flag `i_bit`. Not for a parallel loop over flags: see the module's note on parallel loops.
    pure subroutine bitmask_clear(n_words, bitmask, i_bit)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the mask
        integer(int32), intent(in) :: i_bit
            !! index of the flag, `1 <= i_bit <= n_bits`

        bitmask(M_BIT_MASK_WORD(i_bit)) = ibclr(bitmask(M_BIT_MASK_WORD(i_bit)), M_BIT_MASK_BIT(i_bit))
    end subroutine bitmask_clear

    ! ------------------------------------------------------------------ whole masks, in place

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Sets every flag to `state`.
    pure subroutine bitmask_fill(n_bits, n_words, bitmask, state)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask
        integer(int32), intent(in) :: n_words
            !! number of words of the mask, `ceiling(n_bits/32)`
        integer(int32), dimension(n_words), intent(out) :: bitmask
            !! the mask
        logical(c_bool), intent(in) :: state
            !! the value of every flag

        if (state) then
            bitmask = not(0_int32)
            call clear_tail(n_bits, n_words, bitmask)
        else
            bitmask = 0_int32
        end if
    end subroutine bitmask_fill

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `bitmask = source_bitmask`
    pure subroutine bitmask_copy(n_words, source_bitmask, bitmask)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: source_bitmask
            !! the mask to copy
        integer(int32), dimension(n_words), intent(out) :: bitmask
            !! the copy

        bitmask = source_bitmask
    end subroutine bitmask_copy

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `bitmask = bitmask .and. other_bitmask`
    pure subroutine bitmask_and(n_words, other_bitmask, bitmask)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: other_bitmask
            !! the second operand
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the first operand, and the result

        bitmask = iand(bitmask, other_bitmask)
    end subroutine bitmask_and

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `bitmask = bitmask .or. other_bitmask`
    pure subroutine bitmask_or(n_words, other_bitmask, bitmask)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: other_bitmask
            !! the second operand
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the first operand, and the result

        bitmask = ior(bitmask, other_bitmask)
    end subroutine bitmask_or

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `bitmask = bitmask .and. .not. other_bitmask`: clears every flag `other_bitmask` sets.
    pure subroutine bitmask_and_not(n_words, other_bitmask, bitmask)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: other_bitmask
            !! the flags to clear
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the first operand, and the result

        bitmask = iand(bitmask, not(other_bitmask))
    end subroutine bitmask_and_not

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `bitmask = bitmask .neqv. other_bitmask`
    pure subroutine bitmask_xor(n_words, other_bitmask, bitmask)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: other_bitmask
            !! the second operand
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the first operand, and the result

        bitmask = ieor(bitmask, other_bitmask)
    end subroutine bitmask_xor

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `bitmask = .not. bitmask`, the bits past `n_bits` kept zero.
    pure subroutine bitmask_not(n_bits, n_words, bitmask)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask
        integer(int32), intent(in) :: n_words
            !! number of words of the mask, `ceiling(n_bits/32)`
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the mask

        bitmask = not(bitmask)
        call clear_tail(n_bits, n_words, bitmask)
    end subroutine bitmask_not

    ! ------------------------------------------------------------------ queries

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The number of set flags.
    pure integer(int32) function bitmask_count(n_words, bitmask) result(n_set)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask

        integer(int32) :: i_word

        n_set = 0_int32
        do i_word = 1, n_words
            n_set = n_set + int(popcnt(bitmask(i_word)), int32)
        end do
    end function bitmask_count

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `count(left .and. right)`, without building the conjunction.
    pure integer(int32) function bitmask_and_count(n_words, left_bitmask, right_bitmask) result(n_set)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: left_bitmask
            !! the first operand
        integer(int32), dimension(n_words), intent(in) :: right_bitmask
            !! the second operand

        integer(int32) :: i_word

        n_set = 0_int32
        do i_word = 1, n_words
            n_set = n_set + int(popcnt(iand(left_bitmask(i_word), right_bitmask(i_word))), int32)
        end do
    end function bitmask_and_count

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether any flag is set.
    pure logical function bitmask_any(n_words, bitmask) result(any_set)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask

        integer(int32) :: i_word

        any_set = .true.
        do i_word = 1, n_words
            if (bitmask(i_word) /= 0_int32) return
        end do
        any_set = .false.
    end function bitmask_any

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether every flag is set. `.true.` for a mask without flags.
    pure logical function bitmask_all(n_bits, n_words, bitmask) result(all_set)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask
        integer(int32), intent(in) :: n_words
            !! number of words of the mask, `ceiling(n_bits/32)`
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask

        integer(int32) :: i_word

        all_set = .false.
        do i_word = 1, n_words - 1_int32
            if (bitmask(i_word) /= not(0_int32)) return
        end do
        if (n_words > 0_int32) then
            if (bitmask(n_words) /= last_word_flags(n_bits)) return
        end if
        all_set = .true.
    end function bitmask_all

    !> AUTHOR_FRANZ_ERIC_SILL
    !| `any(left .and. right)`: whether the masks share a flag. Stops at the first shared word.
    pure logical function bitmask_any_and(n_words, left_bitmask, right_bitmask) result(any_shared)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: left_bitmask
            !! the first operand
        integer(int32), dimension(n_words), intent(in) :: right_bitmask
            !! the second operand

        integer(int32) :: i_word

        any_shared = .true.
        do i_word = 1, n_words
            if (iand(left_bitmask(i_word), right_bitmask(i_word)) /= 0_int32) return
        end do
        any_shared = .false.
    end function bitmask_any_and

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether both masks set the same flags.
    pure logical function bitmask_equal(n_words, left_bitmask, right_bitmask) result(are_equal)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: left_bitmask
            !! the first mask
        integer(int32), dimension(n_words), intent(in) :: right_bitmask
            !! the second mask

        integer(int32) :: i_word

        are_equal = .false.
        do i_word = 1, n_words
            if (left_bitmask(i_word) /= right_bitmask(i_word)) return
        end do
        are_equal = .true.
    end function bitmask_equal

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Whether every flag `subset_bitmask` sets is also set in `superset_bitmask`.
    pure logical function bitmask_is_subset(n_words, subset_bitmask, superset_bitmask) result(is_subset)
        integer(int32), intent(in) :: n_words
            !! number of words of both masks
        integer(int32), dimension(n_words), intent(in) :: subset_bitmask
            !! the mask that may be contained
        integer(int32), dimension(n_words), intent(in) :: superset_bitmask
            !! the mask that may contain it

        integer(int32) :: i_word

        is_subset = .false.
        do i_word = 1, n_words
            if (iand(subset_bitmask(i_word), not(superset_bitmask(i_word))) /= 0_int32) return
        end do
        is_subset = .true.
    end function bitmask_is_subset

    ! ------------------------------------------------------------------ positions

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The index of the first set flag, or 0 if no flag is set.
    pure integer(int32) function bitmask_first(n_words, bitmask) result(i_bit)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask

        i_bit = bitmask_next(n_words, bitmask, 0_int32)
    end function bitmask_first

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The index of the last set flag, or 0 if no flag is set.
    pure integer(int32) function bitmask_last(n_words, bitmask) result(i_bit)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask

        integer(int32) :: i_word

        do i_word = n_words, 1_int32, -1_int32
            if (bitmask(i_word) /= 0_int32) then
                ! not i_word*32 - leadz: i_word*32 overflows for the last word of a mask of 2**31 - 1 flags
                i_bit = (i_word - 1_int32)*M_BIT_MASK_WORD_BITS &
                        + (M_BIT_MASK_WORD_BITS - int(leadz(bitmask(i_word)), int32))
                return
            end if
        end do
        i_bit = 0_int32
    end function bitmask_last

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The index of the first set flag after `after_bit`, or 0 if there is none. Stepping through
    !| the set flags of a mask:
    !|
    !| ```fortran
    !| i_bit = bitmask_next(n_words, bitmask, 0_int32)
    !| do while (i_bit /= 0_int32)
    !|     ! ... flag i_bit is set
    !|     i_bit = bitmask_next(n_words, bitmask, i_bit)
    !| end do
    !| ```
    pure integer(int32) function bitmask_next(n_words, bitmask, after_bit) result(i_bit)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask
        integer(int32), intent(in) :: after_bit
            !! the search starts at flag `after_bit + 1`; 0, or anything below, searches the whole mask

        integer(int32) :: i_word, word, start_after

        i_bit = 0_int32
        start_after = max(after_bit, 0_int32)
        i_word = start_after/M_BIT_MASK_WORD_BITS + 1_int32
        if (i_word > n_words) return
        ! the flags up to start_after are masked out of the first word searched
        word = iand(bitmask(i_word), shiftl(not(0_int32), mod(start_after, M_BIT_MASK_WORD_BITS)))
        do
            if (word /= 0_int32) then
                i_bit = (i_word - 1_int32)*M_BIT_MASK_WORD_BITS + trailing_zeros(word) + 1_int32
                return
            end if
            i_word = i_word + 1_int32
            if (i_word > n_words) return
            word = bitmask(i_word)
        end do
    end function bitmask_next

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The indices of the set flags, in ascending order. With `n_set` from
    !| [[f42_bit_masks_impl(module):bitmask_count(function)]] that is all of them; a smaller
    !| `n_set` gives the first `n_set`, and a larger one fills the rest of `indices` with 0.
    pure subroutine bitmask_indices(n_words, n_set, bitmask, indices)
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), intent(in) :: n_set
            !! number of indices to collect
        integer(int32), dimension(n_words), intent(in) :: bitmask
            !! the mask
        integer(int32), dimension(n_set), intent(out) :: indices
            !! the indices of the first `n_set` set flags, then 0

        integer(int32) :: i_word, i_index, word, position

        i_index = 0_int32
        do i_word = 1, n_words
            word = bitmask(i_word)
            do while (word /= 0_int32)
                if (i_index == n_set) return
                position = trailing_zeros(word)
                i_index = i_index + 1_int32
                indices(i_index) = (i_word - 1_int32)*M_BIT_MASK_WORD_BITS + position + 1_int32
                word = ibclr(word, position)
            end do
        end do
        indices(i_index + 1_int32:) = 0_int32
    end subroutine bitmask_indices

    ! ------------------------------------------------------------------ helpers

    !> AUTHOR_FRANZ_ERIC_SILL
    !| One word packed from up to 32 flags. Branch-free: a branch per flag mispredicts on random masks.
    pure integer(int32) function packed_word(n_flags, flags) result(word)
        integer(int32), intent(in) :: n_flags
            !! number of flags in this word, at most 32
        logical(c_bool), dimension(n_flags), intent(in) :: flags
            !! the flags, the first one becoming bit 0

        integer(int32) :: i_flag

        word = 0_int32
        do i_flag = 1, n_flags
            word = ior(word, shiftl(merge(1_int32, 0_int32, flags(i_flag)), i_flag - 1_int32))
        end do
    end function packed_word

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The last word of a mask of `n_bits` flags with exactly its valid bits set.
    pure integer(int32) function last_word_flags(n_bits) result(word)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask

        if (mod(n_bits, M_BIT_MASK_WORD_BITS) == 0_int32) then
            word = not(0_int32)
        else
            word = shiftr(not(0_int32), M_BIT_MASK_WORD_BITS - mod(n_bits, M_BIT_MASK_WORD_BITS))
        end if
    end function last_word_flags

    !> AUTHOR_FRANZ_ERIC_SILL
    !| Clears the bits past `n_bits` in the last word.
    pure subroutine clear_tail(n_bits, n_words, bitmask)
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
        integer(int32), dimension(n_words), intent(inout) :: bitmask
            !! the mask

        if (n_words > 0_int32) bitmask(n_words) = iand(bitmask(n_words), last_word_flags(n_bits))
    end subroutine clear_tail

    !> AUTHOR_FRANZ_ERIC_SILL
    !| The number of trailing zero bits of a nonzero word. Not `trailz`, which ifx does not inline
    !| (it calls a runtime routine per word). `iand(not(x), x - 1)` sets exactly the bits below the
    !| lowest set bit; the word is widened first because `x - 1` overflows int32 when only bit 31 is set.
    pure integer(int32) function trailing_zeros(word) result(n_zeros)
        integer(int32), intent(in) :: word
            !! the word, not zero

        integer(int64) :: widened

        widened = int(word, int64)
        n_zeros = int(popcnt(iand(not(widened), widened - 1_int64)), int32)
    end function trailing_zeros

end module f42_bit_masks_impl
