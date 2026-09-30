#include <src/macros.h>

!> A bit mask holds the same information as a `logical(c_bool)` array in an eighth of the memory,
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
!| **The bits past `n_bits` in the last word are always zero.** Every procedure here keeps it that
!| way, and the counts, [[f42_bit_masks_impl(module):bit_mask_all(function)]] and
!| [[f42_bit_masks_impl(module):bit_mask_equal(function)]] rely on it.
!|
!| **Several masks** are the columns of a `(n_words, n_masks)` array. Every column starts on a
!| new word, so `bit_masks(:, i_mask)` is a mask of its own and passes to any procedure here
!| without a copy.
!|
!| **Parallel loops write whole words.** Two iterations that set flags in the same word race, so
!| a `do concurrent` over flags may only read a mask; one that writes runs over words or columns.
!|
!| **Hot loops** use the whole-mask procedures, or step from one set flag to the next with
!| [[f42_bit_masks_impl(module):bit_mask_next(function)]]. A procedure call per flag costs more
!| than the flag, because nothing is inlined across modules; where a loop must touch single
!| flags, the bit-mask macros in `src/macros.h` test a flag or locate its word inline, and one
!| of them sizes a mask.
!|
!| Generated from [[f42_bit_masks_impl(module)]]; do not edit -- regenerate instead.
module f42_bit_masks
    use f42_safeguard
    use f42_bit_masks_impl, only: bit_mask_from_logical_impl, bit_mask_n_words_impl, bit_mask_test_impl, bit_mask_to_logical_impl
    use f42_bit_masks_impl, only: bit_masks_from_logical_2D_impl, bit_masks_to_logical_2D_impl
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_errors, only: set_ok, is_err, validate_dimension_size, validate_in_range_int
    M_IMPLICIT_NONE
    private

    public :: bit_mask_n_words
    public :: bit_mask_test
    public :: bit_mask_from_logical
    public :: bit_mask_to_logical
    public :: bit_masks_from_logical_2D
    public :: bit_masks_to_logical_2D

contains

    !> summary: Validates its inputs, then calls [[f42_bit_masks_impl(module):bit_mask_n_words_impl]].
    !| Inside Fortran, a bit-mask macro in `src/macros.h` gives the same number.
    pure subroutine bit_mask_n_words(&
            n_bits,&
            n_words,&
            ierr&
        )
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask; no flags take no words
            !! The minimum valid value is `0_int32`.
        integer(int32), intent(out) :: n_words
            !! number of words the mask takes
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bits, ierr, arg_pos=1_int32, min=0_int32)
        if (is_err(ierr)) return
#endif

        call bit_mask_n_words_impl(&
            n_bits = n_bits,&
            n_words = n_words&
        )
    end subroutine bit_mask_n_words

    !> summary: Validates its inputs, then calls [[f42_bit_masks_impl(module):bit_mask_test_impl]].
    !| Inside Fortran, a loop tests a flag inline with a bit-mask macro from `src/macros.h` instead.
    pure subroutine bit_mask_test(&
            n_bits,&
            n_words,&
            bit_mask,&
            i_bit,&
            is_set,&
            ierr&
        )
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        integer(int32), intent(in) :: n_bits
            !! number of flags in the mask
            !! The minimum valid value is `1_int32`.
        integer(int32), dimension(n_words), intent(in) :: bit_mask
            !! the mask
        integer(int32), intent(in) :: i_bit
            !! index of the flag, counting from 1
            !! The minimum valid value is `1_int32`.
            !! The maximum valid value is `n_bits`.
        logical(c_bool), intent(out) :: is_set
            !! `.true.` if the flag is set; `.false.` for `i_bit > n_bits`, which the bits past
            !! `n_bits` would say too, as long as the mask keeps them zero
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bits, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_words, ierr, arg_pos=2_int32, min=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))), max=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))))
        call validate_in_range_int(i_bit, ierr, arg_pos=4_int32, min=1_int32, max=n_bits)
        if (is_err(ierr)) return
#endif

        call bit_mask_test_impl(&
            n_bits = n_bits,&
            n_words = n_words,&
            bit_mask = bit_mask,&
            i_bit = i_bit,&
            is_set = is_set&
        )
    end subroutine bit_mask_test

    !> summary: Validates its inputs, then calls [[f42_bit_masks_impl(module):bit_mask_from_logical_impl]].
    pure subroutine bit_mask_from_logical(&
            n_bits,&
            n_words,&
            flags,&
            bit_mask,&
            ierr&
        )
        integer(int32), intent(in) :: n_bits
            !! number of flags
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        logical(c_bool), dimension(n_bits), intent(in) :: flags
            !! the flags to pack
        integer(int32), dimension(n_words), intent(out) :: bit_mask
            !! the mask, flag `i` set where `flags(i)` is `.true.`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bits, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_words, ierr, arg_pos=2_int32, min=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))), max=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))))
        if (is_err(ierr)) return
#endif

        call bit_mask_from_logical_impl(&
            n_bits = n_bits,&
            n_words = n_words,&
            flags = flags,&
            bit_mask = bit_mask&
        )
    end subroutine bit_mask_from_logical

    !> summary: Validates its inputs, then calls [[f42_bit_masks_impl(module):bit_mask_to_logical_impl]].
    pure subroutine bit_mask_to_logical(&
            n_bits,&
            n_words,&
            bit_mask,&
            flags,&
            ierr&
        )
        integer(int32), intent(in) :: n_bits
            !! number of flags
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_words
            !! number of words of the mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        integer(int32), dimension(n_words), intent(in) :: bit_mask
            !! the mask
        logical(c_bool), dimension(n_bits), intent(out) :: flags
            !! `flags(i)` is `.true.` where flag `i` of the mask is set
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bits, ierr, arg_pos=1_int32, min=1_int32)
        call validate_in_range_int(n_words, ierr, arg_pos=2_int32, min=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))), max=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))))
        if (is_err(ierr)) return
#endif

        call bit_mask_to_logical_impl(&
            n_bits = n_bits,&
            n_words = n_words,&
            bit_mask = bit_mask,&
            flags = flags&
        )
    end subroutine bit_mask_to_logical

    !> summary: Validates its inputs, then calls [[f42_bit_masks_impl(module):bit_masks_from_logical_2D_impl]].
    pure subroutine bit_masks_from_logical_2D(&
            n_bits,&
            n_masks,&
            n_words,&
            flags,&
            bit_masks,&
            ierr&
        )
        integer(int32), intent(in) :: n_bits
            !! number of flags per mask
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_masks
            !! number of masks
        integer(int32), intent(in) :: n_words
            !! number of words per mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        logical(c_bool), dimension(n_bits, n_masks), intent(in) :: flags
            !! the flags, one mask per column
        integer(int32), dimension(n_words, n_masks), intent(out) :: bit_masks
            !! the masks, column `j` packed from `flags(:, j)`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bits, ierr, arg_pos=1_int32, min=1_int32)
        call validate_dimension_size(n_masks, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_words, ierr, arg_pos=3_int32, min=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))), max=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))))
        if (is_err(ierr)) return
#endif

        call bit_masks_from_logical_2D_impl(&
            n_bits = n_bits,&
            n_masks = n_masks,&
            n_words = n_words,&
            flags = flags,&
            bit_masks = bit_masks&
        )
    end subroutine bit_masks_from_logical_2D

    !> summary: Validates its inputs, then calls [[f42_bit_masks_impl(module):bit_masks_to_logical_2D_impl]].
    pure subroutine bit_masks_to_logical_2D(&
            n_bits,&
            n_masks,&
            n_words,&
            bit_masks,&
            flags,&
            ierr&
        )
        integer(int32), intent(in) :: n_bits
            !! number of flags per mask
            !! The minimum valid value is `1_int32`.
        integer(int32), intent(in) :: n_masks
            !! number of masks
        integer(int32), intent(in) :: n_words
            !! number of words per mask
            !! The minimum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
            !! The maximum valid value is `((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32)))`.
        integer(int32), dimension(n_words, n_masks), intent(in) :: bit_masks
            !! the masks, one per column
        logical(c_bool), dimension(n_bits, n_masks), intent(out) :: flags
            !! the flags, column `j` unpacked from `bit_masks(:, j)`
        integer(int32), intent(out) :: ierr
            !! Error code; zero on success, non-zero on failure.

        call set_ok(ierr)
#ifndef NO_INPUT_VALIDATION
        call validate_in_range_int(n_bits, ierr, arg_pos=1_int32, min=1_int32)
        call validate_dimension_size(n_masks, ierr, arg_pos=2_int32)
        call validate_in_range_int(n_words, ierr, arg_pos=3_int32, min=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))), max=((n_bits)/32_int32 + min(1_int32, mod((n_bits), 32_int32))))
        if (is_err(ierr)) return
#endif

        call bit_masks_to_logical_2D_impl(&
            n_bits = n_bits,&
            n_masks = n_masks,&
            n_words = n_words,&
            bit_masks = bit_masks,&
            flags = flags&
        )
    end subroutine bit_masks_to_logical_2D

end module f42_bit_masks
