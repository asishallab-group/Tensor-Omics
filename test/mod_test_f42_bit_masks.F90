#include <src/macros.h>

!> Unit test suite for f42_bit_masks_impl.
!|
!| Every expected word is derived by hand from the documented layout, never from the module's
!| own formulas: flag `i` (counting from 1) is bit `mod(i - 1, 32)` of word `(i - 1)/32 + 1`,
!| least significant bit first; all 32 bits are used, so a word with bit 31 set is negative;
!| the bits past `n_bits` in the last word are zero; several masks are the columns of an
!| `(n_words, n_masks)` array, each column padded to whole words.
module mod_test_f42_bit_masks
    use f42_bit_masks_impl
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use test_suite, only: test_case
    implicit none
    public

    !> Flag 32 alone: bit 31 of a word, which makes the word -2^31. Not `-huge(0_int32) - 1_int32`,
    !| which `-pedantic` flags as outside the symmetric range the standard implies.
    integer(int32), parameter :: ONLY_BIT_31 = ibset(0_int32, 31)
    !> Every bit of a word set.
    integer(int32), parameter :: ALL_BITS = -1_int32

    !> The sample mask most tests share: 40 flags in 2 words, flags 1, 2, 4, 32 | 33, 36, 40 set.
    integer(int32), parameter :: SAMPLE_N_BITS = 40_int32
    integer(int32), parameter :: SAMPLE_N_WORDS = 2_int32
    integer(int32), parameter :: SAMPLE_SET_FLAGS(7) = [1_int32, 2_int32, 4_int32, 32_int32, 33_int32, 36_int32, 40_int32]
    !| word 1: flags 1, 2, 4 -> bits 0, 1, 3 -> 1 + 2 + 8 = 11; flag 32 -> bit 31 -> 11 - 2^31 = -2147483637
    !| word 2: flags 33, 36, 40 -> bits 0, 3, 7 -> 1 + 8 + 128 = 137
    integer(int32), parameter :: SAMPLE_MASK(2) = [-2147483637_int32, 137_int32]

    !> A second operand for the binary operations: flags 2, 3, 32 | 33, 36.
    !| word 1: bits 1, 2 -> 2 + 4 = 6; bit 31 -> 6 - 2^31 = -2147483642
    !| word 2: bits 0, 3 -> 1 + 8 = 9
    integer(int32), parameter :: OTHER_MASK(2) = [-2147483642_int32, 9_int32]

contains

    !> Get array of all available tests.
    function get_all_tests_f42_bit_masks() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(43))
        all_tests(1) = test_case("test_n_words", test_n_words)
        all_tests(2) = test_case("test_n_words_macro", test_n_words_macro)
        all_tests(3) = test_case("test_flag_is_set", test_flag_is_set)
        all_tests(4) = test_case("test_flag_is_set_macro", test_flag_is_set_macro)
        all_tests(5) = test_case("test_from_logical_40_flags", test_from_logical_40_flags)
        all_tests(6) = test_case("test_from_logical_64_flags", test_from_logical_64_flags)
        all_tests(7) = test_case("test_from_logical_clears_tail", test_from_logical_clears_tail)
        all_tests(8) = test_case("test_to_logical_40_flags", test_to_logical_40_flags)
        all_tests(9) = test_case("test_logical_round_trip", test_logical_round_trip)
        all_tests(10) = test_case("test_from_logical_2D", test_from_logical_2D)
        all_tests(11) = test_case("test_to_logical_2D", test_to_logical_2D)
        all_tests(12) = test_case("test_logical_2D_round_trip", test_logical_2D_round_trip)
        all_tests(13) = test_case("test_set_flags", test_set_flags)
        all_tests(14) = test_case("test_set_keeps_set_flag", test_set_keeps_set_flag)
        all_tests(15) = test_case("test_clear_flags", test_clear_flags)
        all_tests(16) = test_case("test_clear_keeps_clear_flag", test_clear_keeps_clear_flag)
        all_tests(17) = test_case("test_fill_true_clears_tail", test_fill_true_clears_tail)
        all_tests(18) = test_case("test_fill_true_whole_words", test_fill_true_whole_words)
        all_tests(19) = test_case("test_fill_false", test_fill_false)
        all_tests(20) = test_case("test_copy", test_copy)
        all_tests(21) = test_case("test_and", test_and)
        all_tests(22) = test_case("test_or", test_or)
        all_tests(23) = test_case("test_and_not", test_and_not)
        all_tests(24) = test_case("test_xor", test_xor)
        all_tests(25) = test_case("test_not_clears_tail", test_not_clears_tail)
        all_tests(26) = test_case("test_not_twice_restores", test_not_twice_restores)
        all_tests(27) = test_case("test_not_whole_words", test_not_whole_words)
        all_tests(28) = test_case("test_count", test_count)
        all_tests(29) = test_case("test_and_count", test_and_count)
        all_tests(30) = test_case("test_any", test_any)
        all_tests(31) = test_case("test_all", test_all)
        all_tests(32) = test_case("test_any_and", test_any_and)
        all_tests(33) = test_case("test_equal", test_equal)
        all_tests(34) = test_case("test_is_subset", test_is_subset)
        all_tests(35) = test_case("test_first_last_empty", test_first_last_empty)
        all_tests(36) = test_case("test_first_last_flag_32", test_first_last_flag_32)
        all_tests(37) = test_case("test_first_last_several_words", test_first_last_several_words)
        all_tests(38) = test_case("test_next", test_next)
        all_tests(39) = test_case("test_next_stepping_loop", test_next_stepping_loop)
        all_tests(40) = test_case("test_indices_all", test_indices_all)
        all_tests(41) = test_case("test_indices_first_n_set", test_indices_first_n_set)
        all_tests(42) = test_case("test_next_negative_after_bit", test_next_negative_after_bit)
        all_tests(43) = test_case("test_indices_larger_n_set", test_indices_larger_n_set)
    end function get_all_tests_f42_bit_masks

    ! ------------------------------------------------------------------ helpers

    !> `bit_mask_test_impl` as a function, so a test reads as one assertion per flag.
    function flag_is_set(n_bits, n_words, bit_mask, i_bit) result(is_set)
        integer(int32), intent(in) :: n_bits, n_words, i_bit
        integer(int32), intent(in) :: bit_mask(n_words)
        logical(c_bool) :: is_set

        call bit_mask_test_impl(n_bits, n_words, bit_mask, i_bit, is_set)
    end function flag_is_set

    ! ------------------------------------------------------------------ published

    !> Words per mask: none for no flags, at and around the word boundaries, and at huge(0_int32)
    !| without overflow.
    subroutine test_n_words()
        integer(int32) :: n_words

        call bit_mask_n_words_impl(0_int32, n_words)
        call assert_equal_int(n_words, 0_int32, "test_n_words: 0 flags take 0 words")
        call bit_mask_n_words_impl(1_int32, n_words)
        call assert_equal_int(n_words, 1_int32, "test_n_words: 1 flag takes 1 word")
        call bit_mask_n_words_impl(31_int32, n_words)
        call assert_equal_int(n_words, 1_int32, "test_n_words: 31 flags take 1 word")
        call bit_mask_n_words_impl(32_int32, n_words)
        call assert_equal_int(n_words, 1_int32, "test_n_words: 32 flags take 1 word")
        call bit_mask_n_words_impl(33_int32, n_words)
        call assert_equal_int(n_words, 2_int32, "test_n_words: 33 flags take 2 words")
        call bit_mask_n_words_impl(64_int32, n_words)
        call assert_equal_int(n_words, 2_int32, "test_n_words: 64 flags take 2 words")
        call bit_mask_n_words_impl(65_int32, n_words)
        call assert_equal_int(n_words, 3_int32, "test_n_words: 65 flags take 3 words")
        ! huge = 2^31 - 1 = 67108863*32 + 31 -> 67108863 full words and one partial = 2^26
        call bit_mask_n_words_impl(huge(0_int32), n_words)
        call assert_equal_int(n_words, 67108864_int32, "test_n_words: huge(0_int32) flags take 2^26 words")
    end subroutine test_n_words

    !> The macro agrees with the procedure.
    subroutine test_n_words_macro()
        integer(int32) :: n_bits, n_words

        n_bits = 0_int32
        n_words = M_BIT_MASK_N_WORDS(n_bits)
        call assert_equal_int(n_words, 0_int32, "test_n_words_macro: 0 flags take 0 words")
        n_bits = 33_int32
        n_words = M_BIT_MASK_N_WORDS(n_bits)
        call assert_equal_int(n_words, 2_int32, "test_n_words_macro: 33 flags take 2 words")
        n_bits = huge(0_int32)
        n_words = M_BIT_MASK_N_WORDS(n_bits)
        call assert_equal_int(n_words, 67108864_int32, "test_n_words_macro: huge(0_int32) flags take 2^26 words")
    end subroutine test_n_words_macro

    !> Single flags in word 1, at bit 31, first of word 2, and the last flag.
    subroutine test_flag_is_set()
        integer(int32), parameter :: n_bits = 70_int32, n_words = 3_int32
        integer(int32) :: bit_mask(n_words)

        ! word 1: flags 1, 5 -> bits 0, 4 -> 1 + 16 = 17; flag 32 -> bit 31 -> 17 - 2^31 = -2147483631
        ! word 2: flag 33 -> bit 0 -> 1
        ! word 3: flag 70 -> (70 - 1) = 2*32 + 5 -> bit 5 -> 32
        bit_mask = [-2147483631_int32, 1_int32, 32_int32]

        call assert_true(flag_is_set(n_bits, n_words, bit_mask, 1_int32), "test_flag_is_set: flag 1")
        call assert_true(flag_is_set(n_bits, n_words, bit_mask, 5_int32), "test_flag_is_set: flag 5")
        call assert_true(flag_is_set(n_bits, n_words, bit_mask, 32_int32), "test_flag_is_set: flag 32 (bit 31)")
        call assert_true(flag_is_set(n_bits, n_words, bit_mask, 33_int32), "test_flag_is_set: flag 33 (word 2)")
        call assert_true(flag_is_set(n_bits, n_words, bit_mask, 70_int32), "test_flag_is_set: flag 70 (last)")
        call assert_false(flag_is_set(n_bits, n_words, bit_mask, 2_int32), "test_flag_is_set: flag 2 is clear")
        call assert_false(flag_is_set(n_bits, n_words, bit_mask, 31_int32), "test_flag_is_set: flag 31 is clear")
        call assert_false(flag_is_set(n_bits, n_words, bit_mask, 34_int32), "test_flag_is_set: flag 34 is clear")
        call assert_false(flag_is_set(n_bits, n_words, bit_mask, 64_int32), "test_flag_is_set: flag 64 is clear")
        call assert_false(flag_is_set(n_bits, n_words, bit_mask, 69_int32), "test_flag_is_set: flag 69 is clear")
    end subroutine test_flag_is_set

    !> M_BIT_MASK_TEST reads the same flags as the procedure.
    subroutine test_flag_is_set_macro()
        integer(int32), parameter :: n_words = 3_int32
        integer(int32) :: bit_mask(n_words), i_bit

        ! the mask of test_flag_is_set: flags 1, 5, 32 | 33 | 70
        bit_mask = [-2147483631_int32, 1_int32, 32_int32]

        i_bit = 1_int32
        call assert_true(M_BIT_MASK_TEST(bit_mask, i_bit), "test_flag_is_set_macro: flag 1")
        i_bit = 32_int32
        call assert_true(M_BIT_MASK_TEST(bit_mask, i_bit), "test_flag_is_set_macro: flag 32 (bit 31)")
        i_bit = 33_int32
        call assert_true(M_BIT_MASK_TEST(bit_mask, i_bit), "test_flag_is_set_macro: flag 33 (word 2)")
        i_bit = 70_int32
        call assert_true(M_BIT_MASK_TEST(bit_mask, i_bit), "test_flag_is_set_macro: flag 70 (last)")
        i_bit = 31_int32
        call assert_false(M_BIT_MASK_TEST(bit_mask, i_bit), "test_flag_is_set_macro: flag 31 is clear")
        i_bit = 69_int32
        call assert_false(M_BIT_MASK_TEST(bit_mask, i_bit), "test_flag_is_set_macro: flag 69 is clear")
    end subroutine test_flag_is_set_macro

    !> Packing 40 flags (not a multiple of 32) gives the sample words exactly.
    subroutine test_from_logical_40_flags()
        logical(c_bool) :: flags(SAMPLE_N_BITS)
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        flags = .false._c_bool
        flags(SAMPLE_SET_FLAGS) = .true._c_bool
        bit_mask = ALL_BITS  ! every word must be overwritten, the tail included

        call bit_mask_from_logical_impl(SAMPLE_N_BITS, SAMPLE_N_WORDS, flags, bit_mask)

        call assert_equal_array_int(bit_mask, SAMPLE_MASK, SAMPLE_N_WORDS, "test_from_logical_40_flags: words")
    end subroutine test_from_logical_40_flags

    !> Packing 64 flags (a multiple of 32): bit 31 is set in both words.
    subroutine test_from_logical_64_flags()
        integer(int32), parameter :: n_bits = 64_int32, n_words = 2_int32
        logical(c_bool) :: flags(n_bits)
        integer(int32) :: bit_mask(n_words), expected(n_words)

        flags = .false._c_bool
        flags([1_int32, 32_int32, 33_int32, 64_int32]) = .true._c_bool
        ! each word: bits 0 and 31 -> 1 - 2^31 = -2147483647
        expected = [-2147483647_int32, -2147483647_int32]

        call bit_mask_from_logical_impl(n_bits, n_words, flags, bit_mask)

        call assert_equal_array_int(bit_mask, expected, n_words, "test_from_logical_64_flags: words")
    end subroutine test_from_logical_64_flags

    !> Packing 40 set flags leaves the 24 bits past flag 40 zero.
    subroutine test_from_logical_clears_tail()
        logical(c_bool) :: flags(SAMPLE_N_BITS)
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS)

        flags = .true._c_bool
        bit_mask = ALL_BITS
        ! word 1: all 32 bits -> -1; word 2: flags 33..40 -> bits 0..7 -> 255
        expected = [ALL_BITS, 255_int32]

        call bit_mask_from_logical_impl(SAMPLE_N_BITS, SAMPLE_N_WORDS, flags, bit_mask)

        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_from_logical_clears_tail: words")
    end subroutine test_from_logical_clears_tail

    !> Unpacking the sample words gives exactly the sample flags.
    subroutine test_to_logical_40_flags()
        logical(c_bool) :: flags(SAMPLE_N_BITS), expected(SAMPLE_N_BITS)

        expected = .false._c_bool
        expected(SAMPLE_SET_FLAGS) = .true._c_bool
        flags = .true._c_bool  ! every flag must be overwritten

        call bit_mask_to_logical_impl(SAMPLE_N_BITS, SAMPLE_N_WORDS, SAMPLE_MASK, flags)

        call assert_equal_array_logical(flags, expected, SAMPLE_N_BITS, "test_to_logical_40_flags: flags")
    end subroutine test_to_logical_40_flags

    !> Pack then unpack 70 flags (3 words, the last partial) gives the flags back.
    subroutine test_logical_round_trip()
        integer(int32), parameter :: n_bits = 70_int32, n_words = 3_int32
        logical(c_bool) :: flags(n_bits), unpacked(n_bits)
        integer(int32) :: bit_mask(n_words), i_bit

        do i_bit = 1_int32, n_bits
            flags(i_bit) = logical(mod(i_bit, 3_int32) == 0_int32 .or. i_bit == 32_int32 &
                                   .or. i_bit == 33_int32 .or. i_bit == n_bits, c_bool)
        end do

        call bit_mask_from_logical_impl(n_bits, n_words, flags, bit_mask)
        call bit_mask_to_logical_impl(n_bits, n_words, bit_mask, unpacked)

        call assert_equal_array_logical(unpacked, flags, n_bits, "test_logical_round_trip: flags")
    end subroutine test_logical_round_trip

    !> Three columns of 33 flags: each column packs into its own 2 words.
    subroutine test_from_logical_2D()
        integer(int32), parameter :: n_bits = 33_int32, n_masks = 3_int32, n_words = 2_int32
        logical(c_bool) :: flags(n_bits, n_masks)
        integer(int32) :: bit_masks(n_words, n_masks), expected(n_words, n_masks)

        flags = .false._c_bool
        flags(1, 1) = .true._c_bool       ! column 1: flags 1, 33
        flags(33, 1) = .true._c_bool
        flags(32, 2) = .true._c_bool      ! column 2: flag 32 only
        flags(:, 3) = .true._c_bool       ! column 3: all 33 flags
        bit_masks = ALL_BITS

        ! column 1: [bit 0 -> 1, bit 0 -> 1]
        ! column 2: [bit 31 -> -2^31, 0]
        ! column 3: [all 32 bits -> -1, flag 33 alone in word 2 -> 1]
        expected(:, 1) = [1_int32, 1_int32]
        expected(:, 2) = [ONLY_BIT_31, 0_int32]
        expected(:, 3) = [ALL_BITS, 1_int32]

        call bit_masks_from_logical_2D_impl(n_bits, n_masks, n_words, flags, bit_masks)

        call assert_equal_array_int(bit_masks, expected, n_words*n_masks, "test_from_logical_2D: words", n_rows=n_words)
    end subroutine test_from_logical_2D

    !> Unpacking the words of test_from_logical_2D gives its flags back, column by column.
    subroutine test_to_logical_2D()
        integer(int32), parameter :: n_bits = 33_int32, n_masks = 3_int32, n_words = 2_int32
        logical(c_bool) :: flags(n_bits, n_masks), expected(n_bits, n_masks)
        integer(int32) :: bit_masks(n_words, n_masks)

        bit_masks(:, 1) = [1_int32, 1_int32]
        bit_masks(:, 2) = [ONLY_BIT_31, 0_int32]
        bit_masks(:, 3) = [ALL_BITS, 1_int32]
        expected = .false._c_bool
        expected(1, 1) = .true._c_bool
        expected(33, 1) = .true._c_bool
        expected(32, 2) = .true._c_bool
        expected(:, 3) = .true._c_bool
        flags = .true._c_bool

        call bit_masks_to_logical_2D_impl(n_bits, n_masks, n_words, bit_masks, flags)

        call assert_equal_array_logical(flags, expected, n_bits*n_masks, "test_to_logical_2D: flags", n_rows=n_bits)
    end subroutine test_to_logical_2D

    !> Pack then unpack a (33, 3) flag matrix gives it back.
    subroutine test_logical_2D_round_trip()
        integer(int32), parameter :: n_bits = 33_int32, n_masks = 3_int32, n_words = 2_int32
        logical(c_bool) :: flags(n_bits, n_masks), unpacked(n_bits, n_masks)
        integer(int32) :: bit_masks(n_words, n_masks), i_bit, i_mask

        do i_mask = 1_int32, n_masks
            do i_bit = 1_int32, n_bits
                flags(i_bit, i_mask) = logical(mod(i_bit + i_mask, 4_int32) == 0_int32, c_bool)
            end do
        end do

        call bit_masks_from_logical_2D_impl(n_bits, n_masks, n_words, flags, bit_masks)
        call bit_masks_to_logical_2D_impl(n_bits, n_masks, n_words, bit_masks, unpacked)

        call assert_equal_array_logical(unpacked, flags, n_bits*n_masks, "test_logical_2D_round_trip: flags", &
                                        n_rows=n_bits)
    end subroutine test_logical_2D_round_trip

    ! ------------------------------------------------------------------ single flags

    !> Setting flag 32, then the first flag of word 2, then flag 1, each leaving its neighbours alone.
    subroutine test_set_flags()
        integer(int32), parameter :: n_words = 2_int32
        integer(int32) :: bit_mask(n_words), expected(n_words)

        bit_mask = 0_int32

        call bit_mask_set(n_words, bit_mask, 32_int32)
        expected = [ONLY_BIT_31, 0_int32]                  ! bit 31 of word 1
        call assert_equal_array_int(bit_mask, expected, n_words, "test_set_flags: after flag 32")

        call bit_mask_set(n_words, bit_mask, 33_int32)
        expected = [ONLY_BIT_31, 1_int32]                  ! plus bit 0 of word 2
        call assert_equal_array_int(bit_mask, expected, n_words, "test_set_flags: after flag 33")

        call bit_mask_set(n_words, bit_mask, 1_int32)
        expected = [-2147483647_int32, 1_int32]            ! plus bit 0 of word 1: 1 - 2^31
        call assert_equal_array_int(bit_mask, expected, n_words, "test_set_flags: after flag 1")
    end subroutine test_set_flags

    !> Setting a flag that is already set changes nothing.
    subroutine test_set_keeps_set_flag()
        integer(int32), parameter :: n_words = 2_int32
        integer(int32) :: bit_mask(n_words)

        bit_mask = SAMPLE_MASK
        call bit_mask_set(n_words, bit_mask, 4_int32)
        call bit_mask_set(n_words, bit_mask, 32_int32)
        call bit_mask_set(n_words, bit_mask, 33_int32)

        call assert_equal_array_int(bit_mask, SAMPLE_MASK, n_words, "test_set_keeps_set_flag: words")
    end subroutine test_set_keeps_set_flag

    !> Clearing flag 32, then the first flag of word 2, then flag 1, each leaving its neighbours alone.
    subroutine test_clear_flags()
        integer(int32), parameter :: n_words = 2_int32
        integer(int32) :: bit_mask(n_words), expected(n_words)

        bit_mask = [ALL_BITS, ALL_BITS]

        call bit_mask_clear(n_words, bit_mask, 32_int32)
        expected = [huge(0_int32), ALL_BITS]               ! bits 0..30 of word 1: 2^31 - 1
        call assert_equal_array_int(bit_mask, expected, n_words, "test_clear_flags: after flag 32")

        call bit_mask_clear(n_words, bit_mask, 33_int32)
        expected = [huge(0_int32), -2_int32]               ! word 2 without bit 0: -1 - 1
        call assert_equal_array_int(bit_mask, expected, n_words, "test_clear_flags: after flag 33")

        call bit_mask_clear(n_words, bit_mask, 1_int32)
        expected = [2147483646_int32, -2_int32]            ! word 1 without bits 0 and 31: 2^31 - 2
        call assert_equal_array_int(bit_mask, expected, n_words, "test_clear_flags: after flag 1")
    end subroutine test_clear_flags

    !> Clearing a flag that is already clear changes nothing.
    subroutine test_clear_keeps_clear_flag()
        integer(int32), parameter :: n_words = 2_int32
        integer(int32) :: bit_mask(n_words)

        bit_mask = SAMPLE_MASK
        call bit_mask_clear(n_words, bit_mask, 3_int32)
        call bit_mask_clear(n_words, bit_mask, 31_int32)
        call bit_mask_clear(n_words, bit_mask, 34_int32)

        call assert_equal_array_int(bit_mask, SAMPLE_MASK, n_words, "test_clear_keeps_clear_flag: words")
    end subroutine test_clear_keeps_clear_flag

    ! ------------------------------------------------------------------ whole masks, in place

    !> Filling 40 flags: word 1 all ones, word 2 only its 8 valid bits.
    subroutine test_fill_true_clears_tail()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS)

        bit_mask = [0_int32, ALL_BITS]
        expected = [ALL_BITS, 255_int32]                   ! flags 33..40 -> bits 0..7 -> 255

        call bit_mask_fill(SAMPLE_N_BITS, SAMPLE_N_WORDS, bit_mask, .true._c_bool)

        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_fill_true_clears_tail: words")
    end subroutine test_fill_true_clears_tail

    !> Filling 64 flags: both words all ones.
    subroutine test_fill_true_whole_words()
        integer(int32), parameter :: n_bits = 64_int32, n_words = 2_int32
        integer(int32) :: bit_mask(n_words), expected(n_words)

        bit_mask = 0_int32
        expected = [ALL_BITS, ALL_BITS]

        call bit_mask_fill(n_bits, n_words, bit_mask, .true._c_bool)

        call assert_equal_array_int(bit_mask, expected, n_words, "test_fill_true_whole_words: words")
    end subroutine test_fill_true_whole_words

    !> Filling with .false. clears every word.
    subroutine test_fill_false()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS)

        bit_mask = [ALL_BITS, 255_int32]
        expected = 0_int32

        call bit_mask_fill(SAMPLE_N_BITS, SAMPLE_N_WORDS, bit_mask, .false._c_bool)

        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_fill_false: words")
    end subroutine test_fill_false

    !> A copy has the source's words, bit 31 included.
    subroutine test_copy()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        bit_mask = 0_int32

        call bit_mask_copy(SAMPLE_N_WORDS, SAMPLE_MASK, bit_mask)

        call assert_equal_array_int(bit_mask, SAMPLE_MASK, SAMPLE_N_WORDS, "test_copy: words")
    end subroutine test_copy

    !> sample .and. other: flags 2, 32 | 33, 36.
    subroutine test_and()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS)

        bit_mask = SAMPLE_MASK
        ! word 1: bits {0,1,3,31} .and. {1,2,31} = {1,31} -> 2 - 2^31 = -2147483646
        ! word 2: bits {0,3,7} .and. {0,3} = {0,3} -> 9
        expected = [-2147483646_int32, 9_int32]

        call bit_mask_and(SAMPLE_N_WORDS, OTHER_MASK, bit_mask)

        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_and: words")
    end subroutine test_and

    !> sample .or. other: flags 1, 2, 3, 4, 32 | 33, 36, 40.
    subroutine test_or()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS)

        bit_mask = SAMPLE_MASK
        ! word 1: bits {0,1,3,31} .or. {1,2,31} = {0,1,2,3,31} -> 15 - 2^31 = -2147483633
        ! word 2: bits {0,3,7} .or. {0,3} = {0,3,7} -> 137
        expected = [-2147483633_int32, 137_int32]

        call bit_mask_or(SAMPLE_N_WORDS, OTHER_MASK, bit_mask)

        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_or: words")
    end subroutine test_or

    !> sample .and. .not. other: flags 1, 4 | 40; and bit 31 survives when other lacks it.
    subroutine test_and_not()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS), other(SAMPLE_N_WORDS)

        bit_mask = SAMPLE_MASK
        ! word 1: bits {0,1,3,31} minus {1,2,31} = {0,3} -> 9
        ! word 2: bits {0,3,7} minus {0,3} = {7} -> 128
        expected = [9_int32, 128_int32]
        call bit_mask_and_not(SAMPLE_N_WORDS, OTHER_MASK, bit_mask)
        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_and_not: sample minus other")

        bit_mask = [-2147483647_int32, 0_int32]             ! flags 1, 32: 1 - 2^31
        other = [1_int32, 0_int32]                          ! flag 1
        expected = [ONLY_BIT_31, 0_int32]                   ! flag 32 is left
        call bit_mask_and_not(SAMPLE_N_WORDS, other, bit_mask)
        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_and_not: bit 31 survives")
    end subroutine test_and_not

    !> sample .neqv. other: flags 1, 3, 4 | 40; and bit 31 set when only one side has it.
    subroutine test_xor()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS), other(SAMPLE_N_WORDS)

        bit_mask = SAMPLE_MASK
        ! word 1: bits {0,1,3,31} xor {1,2,31} = {0,2,3} -> 1 + 4 + 8 = 13
        ! word 2: bits {0,3,7} xor {0,3} = {7} -> 128
        expected = [13_int32, 128_int32]
        call bit_mask_xor(SAMPLE_N_WORDS, OTHER_MASK, bit_mask)
        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_xor: sample xor other")

        bit_mask = [ONLY_BIT_31, 0_int32]                   ! flag 32
        other = [1_int32, 0_int32]                          ! flag 1
        expected = [-2147483647_int32, 0_int32]             ! flags 1, 32: 1 - 2^31
        call bit_mask_xor(SAMPLE_N_WORDS, other, bit_mask)
        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_xor: bit 31 from one side")
    end subroutine test_xor

    !> Negating 40 flags: the 24 tail bits of word 2 stay zero.
    subroutine test_not_clears_tail()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), expected(SAMPLE_N_WORDS)

        bit_mask = SAMPLE_MASK
        ! word 1: every bit but {0,1,3,31} -> (2^31 - 1) - 11 = 2147483636
        ! word 2: bits 0..7 but {0,3,7} = {1,2,4,5,6} -> 2 + 4 + 16 + 32 + 64 = 118
        expected = [2147483636_int32, 118_int32]

        call bit_mask_not(SAMPLE_N_BITS, SAMPLE_N_WORDS, bit_mask)

        call assert_equal_array_int(bit_mask, expected, SAMPLE_N_WORDS, "test_not_clears_tail: words")
    end subroutine test_not_clears_tail

    !> Negating twice gives the mask back.
    subroutine test_not_twice_restores()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        bit_mask = SAMPLE_MASK

        call bit_mask_not(SAMPLE_N_BITS, SAMPLE_N_WORDS, bit_mask)
        call bit_mask_not(SAMPLE_N_BITS, SAMPLE_N_WORDS, bit_mask)

        call assert_equal_array_int(bit_mask, SAMPLE_MASK, SAMPLE_N_WORDS, "test_not_twice_restores: words")
    end subroutine test_not_twice_restores

    !> Negating 64 flags (no tail): every bit flips, bit 31 included.
    subroutine test_not_whole_words()
        integer(int32), parameter :: n_bits = 64_int32, n_words = 2_int32
        integer(int32) :: bit_mask(n_words), expected(n_words)

        bit_mask = [ONLY_BIT_31, 0_int32]
        expected = [huge(0_int32), ALL_BITS]               ! bits 0..30 -> 2^31 - 1; all 32 bits -> -1

        call bit_mask_not(n_bits, n_words, bit_mask)

        call assert_equal_array_int(bit_mask, expected, n_words, "test_not_whole_words: words")
    end subroutine test_not_whole_words

    ! ------------------------------------------------------------------ queries

    !> Set flags counted across words, bit 31 and all-ones words included.
    subroutine test_count()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        call assert_equal_int(bit_mask_count(SAMPLE_N_WORDS, SAMPLE_MASK), 7_int32, &
                              "test_count: the 7 sample flags")

        bit_mask = [ONLY_BIT_31, 0_int32]
        call assert_equal_int(bit_mask_count(SAMPLE_N_WORDS, bit_mask), 1_int32, "test_count: bit 31 alone")

        bit_mask = [ALL_BITS, ALL_BITS]
        call assert_equal_int(bit_mask_count(SAMPLE_N_WORDS, bit_mask), 64_int32, "test_count: all-ones words")

        bit_mask = 0_int32
        call assert_equal_int(bit_mask_count(SAMPLE_N_WORDS, bit_mask), 0_int32, "test_count: no flag")
    end subroutine test_count

    !> count(sample .and. other) = flags 2, 32 | 33, 36 = 4.
    subroutine test_and_count()
        integer(int32) :: left(SAMPLE_N_WORDS), right(SAMPLE_N_WORDS)

        call assert_equal_int(bit_mask_and_count(SAMPLE_N_WORDS, SAMPLE_MASK, OTHER_MASK), 4_int32, &
                              "test_and_count: sample and other share 4 flags")

        left = [5_int32, 0_int32]                           ! flags 1, 3
        right = [2_int32, 1_int32]                          ! flags 2, 33
        call assert_equal_int(bit_mask_and_count(SAMPLE_N_WORDS, left, right), 0_int32, &
                              "test_and_count: disjoint masks")
    end subroutine test_and_count

    !> Any flag set: not for a zero or empty mask; yes for bit 31 alone and for the last word alone.
    subroutine test_any()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS), empty_mask(0)

        bit_mask = 0_int32
        call assert_false(bit_mask_any(SAMPLE_N_WORDS, bit_mask), "test_any: zero mask")
        call assert_false(bit_mask_any(0_int32, empty_mask), "test_any: mask of no words")

        bit_mask = [ONLY_BIT_31, 0_int32]
        call assert_true(bit_mask_any(SAMPLE_N_WORDS, bit_mask), "test_any: only bit 31 (negative word)")

        bit_mask = [0_int32, 1_int32]
        call assert_true(bit_mask_any(SAMPLE_N_WORDS, bit_mask), "test_any: only flag 33")
    end subroutine test_any

    !> Every flag set: full masks of 40 and 64 flags, one missing flag anywhere, no flags.
    subroutine test_all()
        integer(int32) :: bit_mask(2), empty_mask(0)

        bit_mask = [ALL_BITS, 255_int32]                    ! 40 flags, all set
        call assert_true(bit_mask_all(40_int32, 2_int32, bit_mask), "test_all: 40 flags all set")

        bit_mask = [ALL_BITS, 127_int32]                    ! flag 40 (bit 7 of word 2) missing
        call assert_false(bit_mask_all(40_int32, 2_int32, bit_mask), "test_all: 40 flags, flag 40 missing")

        bit_mask = [huge(0_int32), 255_int32]               ! flag 32 (bit 31 of word 1) missing
        call assert_false(bit_mask_all(40_int32, 2_int32, bit_mask), "test_all: 40 flags, flag 32 missing")

        bit_mask = [ALL_BITS, ALL_BITS]
        call assert_true(bit_mask_all(64_int32, 2_int32, bit_mask), "test_all: 64 flags all set")

        bit_mask = [ALL_BITS, huge(0_int32)]                ! flag 64 (bit 31 of word 2) missing
        call assert_false(bit_mask_all(64_int32, 2_int32, bit_mask), "test_all: 64 flags, flag 64 missing")

        call assert_true(bit_mask_all(0_int32, 0_int32, empty_mask), "test_all: a mask without flags")
    end subroutine test_all

    !> Whether two masks share a flag, bit 31 and word 2 included.
    subroutine test_any_and()
        integer(int32) :: left(SAMPLE_N_WORDS), right(SAMPLE_N_WORDS)

        call assert_true(bit_mask_any_and(SAMPLE_N_WORDS, SAMPLE_MASK, OTHER_MASK), &
                         "test_any_and: sample and other share flags")

        left = [ONLY_BIT_31, 0_int32]
        right = [ONLY_BIT_31, 0_int32]
        call assert_true(bit_mask_any_and(SAMPLE_N_WORDS, left, right), "test_any_and: share only flag 32")

        left = [0_int32, 1_int32]
        right = [0_int32, 1_int32]
        call assert_true(bit_mask_any_and(SAMPLE_N_WORDS, left, right), "test_any_and: share only flag 33")

        left = [5_int32, 0_int32]                           ! flags 1, 3
        right = [2_int32, 1_int32]                          ! flags 2, 33
        call assert_false(bit_mask_any_and(SAMPLE_N_WORDS, left, right), "test_any_and: disjoint masks")
    end subroutine test_any_and

    !> Equal masks, and masks differing only in bit 31 or only in word 2.
    subroutine test_equal()
        integer(int32) :: right(SAMPLE_N_WORDS)

        right = SAMPLE_MASK
        call assert_true(bit_mask_equal(SAMPLE_N_WORDS, SAMPLE_MASK, right), "test_equal: same flags")

        right = [11_int32, 137_int32]                       ! the sample without flag 32
        call assert_false(bit_mask_equal(SAMPLE_N_WORDS, SAMPLE_MASK, right), "test_equal: differ in bit 31")

        right = [-2147483637_int32, 9_int32]                ! the sample without flag 40
        call assert_false(bit_mask_equal(SAMPLE_N_WORDS, SAMPLE_MASK, right), "test_equal: differ in word 2")
    end subroutine test_equal

    !> Subsets: a true subset, the reverse, equal masks, the empty mask, one flag outside.
    subroutine test_is_subset()
        integer(int32) :: subset(SAMPLE_N_WORDS), empty_mask(SAMPLE_N_WORDS)

        subset = [-2147483647_int32, 8_int32]               ! flags 1, 32 (1 - 2^31) | 36 (bit 3)
        empty_mask = 0_int32

        call assert_true(bit_mask_is_subset(SAMPLE_N_WORDS, subset, SAMPLE_MASK), &
                         "test_is_subset: flags 1, 32, 36 are in the sample")
        call assert_false(bit_mask_is_subset(SAMPLE_N_WORDS, SAMPLE_MASK, subset), &
                          "test_is_subset: the sample is not in flags 1, 32, 36")
        call assert_true(bit_mask_is_subset(SAMPLE_N_WORDS, SAMPLE_MASK, SAMPLE_MASK), &
                         "test_is_subset: a mask is a subset of itself")
        call assert_true(bit_mask_is_subset(SAMPLE_N_WORDS, empty_mask, SAMPLE_MASK), &
                         "test_is_subset: the empty mask is in the sample")
        call assert_true(bit_mask_is_subset(SAMPLE_N_WORDS, empty_mask, empty_mask), &
                         "test_is_subset: the empty mask is in the empty mask")

        subset = [4_int32, 0_int32]                         ! flag 3, which the sample lacks
        call assert_false(bit_mask_is_subset(SAMPLE_N_WORDS, subset, SAMPLE_MASK), &
                          "test_is_subset: flag 3 is not in the sample")
    end subroutine test_is_subset

    ! ------------------------------------------------------------------ positions

    !> No flag set: first, last and next all give 0.
    subroutine test_first_last_empty()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        bit_mask = 0_int32

        call assert_equal_int(bit_mask_first(SAMPLE_N_WORDS, bit_mask), 0_int32, "test_first_last_empty: first")
        call assert_equal_int(bit_mask_last(SAMPLE_N_WORDS, bit_mask), 0_int32, "test_first_last_empty: last")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, bit_mask, 0_int32), 0_int32, &
                              "test_first_last_empty: next after 0")
    end subroutine test_first_last_empty

    !> Flag 32 alone (a negative word): it is both the first and the last flag.
    subroutine test_first_last_flag_32()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        bit_mask = [ONLY_BIT_31, 0_int32]

        call assert_equal_int(bit_mask_first(SAMPLE_N_WORDS, bit_mask), 32_int32, "test_first_last_flag_32: first")
        call assert_equal_int(bit_mask_last(SAMPLE_N_WORDS, bit_mask), 32_int32, "test_first_last_flag_32: last")
    end subroutine test_first_last_flag_32

    !> First and last flags spread over three words, the first word empty.
    subroutine test_first_last_several_words()
        integer(int32), parameter :: n_words = 3_int32
        integer(int32) :: bit_mask(n_words)

        ! word 2 bit 3 -> flag 32 + 4 = 36; word 3 bit 31 -> flag 64 + 32 = 96
        bit_mask = [0_int32, 8_int32, ONLY_BIT_31]
        call assert_equal_int(bit_mask_first(n_words, bit_mask), 36_int32, "test_first_last_several_words: first 36")
        call assert_equal_int(bit_mask_last(n_words, bit_mask), 96_int32, "test_first_last_several_words: last 96")

        ! word 1 bit 1 -> flag 2; word 3 bit 0 -> flag 65
        bit_mask = [2_int32, 0_int32, 1_int32]
        call assert_equal_int(bit_mask_first(n_words, bit_mask), 2_int32, "test_first_last_several_words: first 2")
        call assert_equal_int(bit_mask_last(n_words, bit_mask), 65_int32, "test_first_last_several_words: last 65")
    end subroutine test_first_last_several_words

    !> The next set flag of the sample (1, 2, 4, 32 | 33, 36, 40) from various starting points.
    subroutine test_next()
        integer(int32) :: bit_mask(SAMPLE_N_WORDS)

        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 0_int32), 1_int32, "test_next: after 0")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 1_int32), 2_int32, "test_next: after 1")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 2_int32), 4_int32, "test_next: after 2")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 4_int32), 32_int32, "test_next: after 4")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 31_int32), 32_int32, "test_next: after 31")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 32_int32), 33_int32, "test_next: after 32")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 36_int32), 40_int32, "test_next: after 36")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 40_int32), 0_int32, &
                              "test_next: after the last set flag")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 32_int32*SAMPLE_N_WORDS), 0_int32, &
                              "test_next: after_bit = 32*n_words, past the end")

        bit_mask = [-2147483647_int32, 0_int32]             ! flags 1, 32
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, bit_mask, 1_int32), 32_int32, &
                              "test_next: skips to bit 31 of the same word")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, bit_mask, 31_int32), 32_int32, &
                              "test_next: after 31 finds 32")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, bit_mask, 32_int32), 0_int32, &
                              "test_next: nothing after 32")
    end subroutine test_next

    !> The documented stepping loop visits exactly the sample's set flags, in order.
    subroutine test_next_stepping_loop()
        integer(int32) :: visited(size(SAMPLE_SET_FLAGS)), n_visited, i_bit

        visited = 0_int32
        n_visited = 0_int32
        i_bit = bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, 0_int32)
        do while (i_bit /= 0_int32)
            n_visited = n_visited + 1_int32
            if (n_visited > size(visited, kind=int32)) exit  ! a next that never reaches 0 must not hang the test
            visited(n_visited) = i_bit
            i_bit = bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, i_bit)
        end do

        call assert_equal_int(n_visited, size(SAMPLE_SET_FLAGS, kind=int32), "test_next_stepping_loop: number visited")
        call assert_equal_array_int(visited, SAMPLE_SET_FLAGS, size(SAMPLE_SET_FLAGS, kind=int32), &
                                    "test_next_stepping_loop: flags visited")
    end subroutine test_next_stepping_loop

    !> With n_set from bit_mask_count, every set flag's index, 32 and 33 included.
    subroutine test_indices_all()
        integer(int32) :: indices(size(SAMPLE_SET_FLAGS)), n_set

        n_set = bit_mask_count(SAMPLE_N_WORDS, SAMPLE_MASK)
        call assert_equal_int(n_set, size(SAMPLE_SET_FLAGS, kind=int32), "test_indices_all: count")

        call bit_mask_indices(SAMPLE_N_WORDS, n_set, SAMPLE_MASK, indices)

        call assert_equal_array_int(indices, SAMPLE_SET_FLAGS, n_set, "test_indices_all: indices")
    end subroutine test_indices_all

    !> A smaller n_set gives only the first n_set indices.
    subroutine test_indices_first_n_set()
        integer(int32), parameter :: n_set = 4_int32
        integer(int32) :: indices(n_set), expected(n_set)

        expected = [1_int32, 2_int32, 4_int32, 32_int32]

        call bit_mask_indices(SAMPLE_N_WORDS, n_set, SAMPLE_MASK, indices)

        call assert_equal_array_int(indices, expected, n_set, "test_indices_first_n_set: indices")
    end subroutine test_indices_first_n_set

    !> A negative after_bit searches the whole mask, like 0.
    subroutine test_next_negative_after_bit()
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, -1_int32), 1_int32, "test_next_negative_after_bit: after -1")
        call assert_equal_int(bit_mask_next(SAMPLE_N_WORDS, SAMPLE_MASK, -huge(0_int32)), 1_int32, &
                              "test_next_negative_after_bit: after -huge")
    end subroutine test_next_negative_after_bit

    !> A larger n_set gives every index, then fills the rest with 0.
    subroutine test_indices_larger_n_set()
        integer(int32), parameter :: n_set = 9_int32
        integer(int32) :: indices(n_set), expected(n_set)

        ! the sample's seven set flags, then two zeros
        expected = [1_int32, 2_int32, 4_int32, 32_int32, 33_int32, 36_int32, 40_int32, 0_int32, 0_int32]
        indices = -1_int32

        call bit_mask_indices(SAMPLE_N_WORDS, n_set, SAMPLE_MASK, indices)

        call assert_equal_array_int(indices, expected, n_set, "test_indices_larger_n_set: indices")
    end subroutine test_indices_larger_n_set

end module mod_test_f42_bit_masks
