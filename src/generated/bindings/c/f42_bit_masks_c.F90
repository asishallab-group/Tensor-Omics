#ifndef NO_C_BINDING
#include <src/macros.h>

!> summary: C-wrappers for [[f42_bit_masks(module)]]
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
!| a `do concurrent` over flags may only read a mask; one that writes runs over words or columns.
!|
!| **Hot loops** use the whole-mask procedures, or step from one set flag to the next with
!| [[f42_bit_masks_impl(module):bitmask_next(function)]]. A procedure call per flag costs more
!| than the flag, because nothing is inlined across modules; where a loop must touch single
!| flags, the bit-mask macros in `src/macros.h` test a flag or locate its word inline, and one
!| of them sizes a mask.
module f42_bit_masks_c
    use f42_safeguard
    use, intrinsic :: iso_c_binding, only: c_associated, c_bool, c_int, c_loc
    use tox_errors, only: set_ok, set_err, ERR_POINTER_NULL
    M_IMPLICIT_NONE
    private

    public :: bitmask_n_words_c
    public :: bitmask_test_c
    public :: bitmask_from_logical_c
    public :: bitmask_to_logical_c
    public :: bitmasks_from_logical_2D_c
    public :: bitmasks_to_logical_2D_c

contains

    !> summary: C-wrapper for [[f42_bit_masks(module):bitmask_n_words(subroutine)]]
    !| Inside Fortran, a bit-mask macro in `src/macros.h` gives the same number.
    subroutine bitmask_n_words_c(&
            n_bits,&
            n_words,&
            ierr&
        ) bind(C, name="bitmask_n_words_c")
        use f42_bit_masks, only: bitmask_n_words

        integer(c_int), intent(in), target :: n_bits
            !! number of flags in the mask; no flags take no words
            !! The minimum valid value is `0_int32`.
        integer(c_int), intent(out), target :: n_words
            !! number of words the mask takes
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bits)
        M_CHECK_NON_NULL(n_words)

        call bitmask_n_words(&
            n_bits = n_bits,&
            n_words = n_words,&
            ierr = ierr&
        )
    end subroutine bitmask_n_words_c

    !> summary: C-wrapper for [[f42_bit_masks(module):bitmask_test(subroutine)]]
    !| Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.
    subroutine bitmask_test_c(&
            n_bits,&
            n_words,&
            bitmask,&
            i_bit,&
            is_set,&
            ierr&
        ) bind(C, name="bitmask_test_c")
        use f42_bit_masks, only: bitmask_test

        integer(c_int), intent(in), target :: n_words
            !! number of words of the mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        integer(c_int), intent(in), target :: n_bits
            !! number of flags in the mask
            !! The minimum valid value is `1_int32`.
        integer(c_int), dimension(n_words), intent(in), target :: bitmask
            !! the mask
        integer(c_int), intent(in), target :: i_bit
            !! index of the flag, counting from 1
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_bits`.
        logical(c_bool), intent(out), target :: is_set
            !! `.true.` if the flag is set
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bits)
        M_CHECK_NON_NULL(n_words)
        M_CHECK_NON_NULL(i_bit)
        M_CHECK_NON_NULL(is_set)
        M_CHECK_ARRAY_NON_NULL(bitmask, n_words)

        call bitmask_test(&
            n_bits = n_bits,&
            n_words = n_words,&
            bitmask = bitmask,&
            i_bit = i_bit,&
            is_set = is_set,&
            ierr = ierr&
        )
    end subroutine bitmask_test_c

    !> summary: C-wrapper for [[f42_bit_masks(module):bitmask_from_logical(subroutine)]]
    subroutine bitmask_from_logical_c(&
            n_bits,&
            n_words,&
            flags,&
            bitmask,&
            ierr&
        ) bind(C, name="bitmask_from_logical_c")
        use f42_bit_masks, only: bitmask_from_logical

        integer(c_int), intent(in), target :: n_bits
            !! number of flags
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_words
            !! number of words of the mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        logical(c_bool), dimension(n_bits), intent(in), target :: flags
            !! the flags to pack
        integer(c_int), dimension(n_words), intent(out), target :: bitmask
            !! the mask, flag `i` set where `flags(i)` is `.true.`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bits)
        M_CHECK_NON_NULL(n_words)
        M_CHECK_ARRAY_NON_NULL(flags, n_bits)
        M_CHECK_ARRAY_NON_NULL(bitmask, n_words)

        call bitmask_from_logical(&
            n_bits = n_bits,&
            n_words = n_words,&
            flags = flags,&
            bitmask = bitmask,&
            ierr = ierr&
        )
    end subroutine bitmask_from_logical_c

    !> summary: C-wrapper for [[f42_bit_masks(module):bitmask_to_logical(subroutine)]]
    subroutine bitmask_to_logical_c(&
            n_bits,&
            n_words,&
            bitmask,&
            flags,&
            ierr&
        ) bind(C, name="bitmask_to_logical_c")
        use f42_bit_masks, only: bitmask_to_logical

        integer(c_int), intent(in), target :: n_bits
            !! number of flags
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_words
            !! number of words of the mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        integer(c_int), dimension(n_words), intent(in), target :: bitmask
            !! the mask
        logical(c_bool), dimension(n_bits), intent(out), target :: flags
            !! `flags(i)` is `.true.` where flag `i` of the mask is set
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bits)
        M_CHECK_NON_NULL(n_words)
        M_CHECK_ARRAY_NON_NULL(bitmask, n_words)
        M_CHECK_ARRAY_NON_NULL(flags, n_bits)

        call bitmask_to_logical(&
            n_bits = n_bits,&
            n_words = n_words,&
            bitmask = bitmask,&
            flags = flags,&
            ierr = ierr&
        )
    end subroutine bitmask_to_logical_c

    !> summary: C-wrapper for [[f42_bit_masks(module):bitmasks_from_logical_2D(subroutine)]]
    subroutine bitmasks_from_logical_2D_c(&
            n_bits,&
            n_masks,&
            n_words,&
            flags,&
            bitmasks,&
            ierr&
        ) bind(C, name="bitmasks_from_logical_2D_c")
        use f42_bit_masks, only: bitmasks_from_logical_2D

        integer(c_int), intent(in), target :: n_bits
            !! number of flags per mask
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_masks
            !! number of masks
        integer(c_int), intent(in), target :: n_words
            !! number of words per mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        logical(c_bool), dimension(n_bits, n_masks), intent(in), target :: flags
            !! the flags, one mask per column
        integer(c_int), dimension(n_words, n_masks), intent(out), target :: bitmasks
            !! the masks, column `j` packed from `flags(:, j)`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bits)
        M_CHECK_NON_NULL(n_masks)
        M_CHECK_NON_NULL(n_words)
        M_CHECK_ARRAY_NON_NULL(flags, n_bits * n_masks)
        M_CHECK_ARRAY_NON_NULL(bitmasks, n_words * n_masks)

        call bitmasks_from_logical_2D(&
            n_bits = n_bits,&
            n_masks = n_masks,&
            n_words = n_words,&
            flags = flags,&
            bitmasks = bitmasks,&
            ierr = ierr&
        )
    end subroutine bitmasks_from_logical_2D_c

    !> summary: C-wrapper for [[f42_bit_masks(module):bitmasks_to_logical_2D(subroutine)]]
    subroutine bitmasks_to_logical_2D_c(&
            n_bits,&
            n_masks,&
            n_words,&
            bitmasks,&
            flags,&
            ierr&
        ) bind(C, name="bitmasks_to_logical_2D_c")
        use f42_bit_masks, only: bitmasks_to_logical_2D

        integer(c_int), intent(in), target :: n_bits
            !! number of flags per mask
            !! The minimum valid value is `1_int32`.
        integer(c_int), intent(in), target :: n_masks
            !! number of masks
        integer(c_int), intent(in), target :: n_words
            !! number of words per mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        integer(c_int), dimension(n_words, n_masks), intent(in), target :: bitmasks
            !! the masks, one per column
        logical(c_bool), dimension(n_bits, n_masks), intent(out), target :: flags
            !! the flags, column `j` unpacked from `bitmasks(:, j)`
        integer(c_int), intent(out), target :: ierr
            !! Error code; zero on success, non-zero on failure.

        M_CHECK_IERR_NON_NULL
        call set_ok(ierr)
        M_CHECK_NON_NULL(n_bits)
        M_CHECK_NON_NULL(n_masks)
        M_CHECK_NON_NULL(n_words)
        M_CHECK_ARRAY_NON_NULL(bitmasks, n_words * n_masks)
        M_CHECK_ARRAY_NON_NULL(flags, n_bits * n_masks)

        call bitmasks_to_logical_2D(&
            n_bits = n_bits,&
            n_masks = n_masks,&
            n_words = n_words,&
            bitmasks = bitmasks,&
            flags = flags,&
            ierr = ierr&
        )
    end subroutine bitmasks_to_logical_2D_c

end module f42_bit_masks_c
#endif
