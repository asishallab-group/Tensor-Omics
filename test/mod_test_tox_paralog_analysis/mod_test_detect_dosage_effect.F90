!> The `detect_dosage_effect` cases, and its expert entry `detect_dosage_effect_expert`.
!|
!| The detection grows subsets of the kept genes one gene at a time, in index order: it starts
!| from each kept gene but the last, and extends a subset only by genes after its last one. A
!| subset of two or more genes whose summed vector is within `max_angle` of the ancestor is a
!| result when its length is at least (1 + `gain_gamma`) times the ancestor's, and is extended
!| further otherwise; one past `max_angle` is dropped. A result is not extended, and a single gene
!| is never a result. So whether a subset S is found depends on S without its last gene: it must
!| have stayed below the gain, every step of the way.
!|
!| The genes in most cases are multiples of one direction, whose angles to the ancestor are 0
!| exactly, and whose lengths are multiples of a power of two, so every comparison is exact. The
!| expected subsets are bit mask words (bit `i - 1` for gene `i`), compared in any order.
module mod_test_detect_dosage_effect
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_paralog_analysis, only: detect_dosage_effect, detect_dosage_effect_expert, &
                                    filter_paralogs_by_pattern_dosage_effect, calc_paralog_subsets_bitmasks_size
    use tox_errors
    use test_suite, only: test_case
    use f42_math_impl, only: PI
    use f42_vector_impl, only: angle_between
    use mod_test_tox_paralog_analysis_fixtures, only: BIT_31, UNWRITTEN_COUNT, N_NON_FINITE, NON_FINITE_NAMES, &
                                                      non_finite_value, check_subsets, MIXED_N_DIMS, MIXED_N_GENES, &
                                                      MIXED_ANCESTOR, MIXED_GENES, MIXED_ANGLES
    implicit none
    private
    public :: get_all_tests_detect_dosage_effect

    !> The binary-weights family of `test_detect_dosage_effect_binary_weights`: gene i is
    !| 2**(i - 1) times [3, 4], so the length of a subset's sum is 5 times the subset's bit mask
    !| word, and the ancestor [15, 20] (length 25) with gain 1 asks for a length of at least 50:
    !| a word of at least 10.
    integer(int32), parameter :: BINARY_N_GENES = 5, BINARY_N_SLOTS = 15
    real(real64), parameter :: BINARY_ANCESTOR(2) = [15.0_real64, 20.0_real64]
    real(real64), parameter :: BINARY_GENES(2, BINARY_N_GENES) = reshape([ &
                               3.0_real64, 4.0_real64, 6.0_real64, 8.0_real64, 12.0_real64, 16.0_real64, &
                               24.0_real64, 32.0_real64, 48.0_real64, 64.0_real64], [2, BINARY_N_GENES])
    !> Its results with every gene kept and subsets of up to four or more genes, pairs first, then
    !| triples, then quadruples (see `test_detect_dosage_effect_binary_weights`).
    integer(int32), parameter :: BINARY_N_RESULTS = 15
    integer(int32), parameter :: BINARY_RESULTS(1, BINARY_N_RESULTS) = reshape([17, 10, 18, 12, 20, 24, 11, 19, 13, 21, &
                                                                                25, 14, 22, 15, 23], [1, BINARY_N_RESULTS])

contains

    !> Every case of `detect_dosage_effect`.
    function get_all_tests_detect_dosage_effect() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(15))
        all_tests(1) = test_case("test_detect_dosage_effect_binary_weights", test_detect_dosage_effect_binary_weights)
        all_tests(2) = test_case("test_detect_dosage_effect_angle_bound", test_detect_dosage_effect_angle_bound)
        all_tests(3) = test_case("test_detect_dosage_effect_defaults", test_detect_dosage_effect_defaults)
        all_tests(4) = test_case("test_detect_dosage_effect_multiword", test_detect_dosage_effect_multiword)
        all_tests(5) = test_case("test_detect_dosage_effect_size_bound", test_detect_dosage_effect_size_bound)
        all_tests(6) = test_case("test_detect_dosage_effect_no_candidates", test_detect_dosage_effect_no_candidates)
        all_tests(7) = test_case("test_detect_dosage_effect_after_filter", test_detect_dosage_effect_after_filter)
        all_tests(8) = test_case("test_detect_dosage_effect_early_errors", test_detect_dosage_effect_early_errors)
        all_tests(9) = test_case("test_detect_dosage_effect_dimensions", test_detect_dosage_effect_dimensions)
        all_tests(10) = test_case("test_detect_dosage_effect_value_ranges", test_detect_dosage_effect_value_ranges)
        all_tests(11) = test_case("test_detect_dosage_effect_rejects_nan_and_inf", test_detect_dosage_effect_rejects_nan_and_inf)
        all_tests(12) = test_case("test_detect_dosage_effect_expert", test_detect_dosage_effect_expert)
        all_tests(13) = test_case("test_detect_dosage_effect_any_scale", test_detect_dosage_effect_any_scale)
        all_tests(14) = test_case("test_angle_between", test_angle_between)
        all_tests(15) = test_case("test_detect_dosage_effect_starts_in_the_array", test_detect_dosage_effect_starts_in_the_array)
    end function get_all_tests_detect_dosage_effect

    !> The binary-weights family (see BINARY_GENES), all five genes kept, `max_angle` 0 -- every
    !| sum points along the ancestor, at angle 0, which the bound includes. A subset is a result
    !| when its word is at least 10 and its word without its highest bit is below 10:
    !| - pairs: {1, 5} 17, {2, 4} 10 (the gain exactly), {2, 5} 18, {3, 4} 12, {3, 5} 20,
    !|   {4, 5} 24; {1, 2} 3, {1, 3} 5, {1, 4} 9 and {2, 3} 6 grow on;
    !| - triples: {1, 2, 4} 11, {1, 2, 5} 19, {1, 3, 4} 13, {1, 3, 5} 21, {1, 4, 5} 25,
    !|   {2, 3, 4} 14, {2, 3, 5} 22; {1, 2, 3} 7 grows on;
    !| - quadruples: {1, 2, 3, 4} 15, {1, 2, 3, 5} 23. Nothing is left to grow to five.
    !| `max_subset_size` 2 finds the 6 pairs, 3 the 13 up to triples, and 4 and 5 all 15. The work
    !| array holds 15 subsets, the worst case for five genes. Note {2, 4} and {1, 2, 4} are both
    !| results: a result's superset is found when its own prefix stayed below the gain.
    !|
    !| In place of the old Python and R dosage checks (three kept genes, results [6, 3, 5]) and of
    !| the old suite's `test_detect_patterns_dosage_effect` (two genes along the ancestor, gain
    !| 0.3 of 0.2: mask 3), and with them of the old helper cases, whose bookkeeping this drives:
    !| results, candidates and new candidates share the one array.
    subroutine test_detect_dosage_effect_binary_weights()
        integer(int32), parameter :: n_expected(5) = [0, 6, 13, 15, 15]
        integer(int32) :: mask(1), subsets(1, BINARY_N_SLOTS), n_results, ierr, max_subset_size
        character(:), allocatable :: label

        mask = 31
        do max_subset_size = 2, 5
            label = "test_detect_dosage_effect_binary_weights: max_subset_size "//max_subset_size
            n_results = UNWRITTEN_COUNT
            call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, max_subset_size, &
                                      subsets, BINARY_N_SLOTS, max_angle=0.0_real64, gain_gamma=1.0_real64, ierr=ierr)
            call assert_err(ierr, ERR_OK, label//": ierr")
            call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, n_expected(max_subset_size), &
                               BINARY_RESULTS(:, 1:n_expected(max_subset_size)), label)
        end do
    end subroutine test_detect_dosage_effect_binary_weights

    !> `max_angle` decides which sums count, and a sum past it is dropped, not grown. Ancestor
    !| [1, 0], genes [1, 1], [1, -1] and [2, 0], gain 2 (a length of at least 3):
    !| - `max_angle` 0.25: {1, 2} sums to [2, 0], at angle 0 but of length 2, and grows; {1, 3}
    !|   and {2, 3}, [3, 1] and [3, -1], are long enough (sqrt(10)) but at angle atan(1/3), about
    !|   0.322, and are dropped; {1, 2, 3}, [4, 0], is the one result, the word 7;
    !| - `max_angle` pi/4 lets {1, 3} and {2, 3} in: the words 5, 6 and 7.
    !| In place of the old suite's `test_detect_patterns_dosage_effect_near_angle_margin`
    !| (a pair at about 0.03 degrees against a bound of 3 degrees) and of the old
    !| `test_angle_between`, whose angles of 0 and about atan(1/3) this reaches through the
    !| detection.
    subroutine test_detect_dosage_effect_angle_bound()
        integer(int32), parameter :: n_genes = 3, n_slots = 3
        real(real64), parameter :: ancestor(2) = [1.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, n_genes) = reshape([1.0_real64, 1.0_real64, 1.0_real64, -1.0_real64, &
                                                                2.0_real64, 0.0_real64], [2, n_genes])
        integer(int32), parameter :: narrow(1, 1) = reshape([7], [1, 1]), wide(1, 3) = reshape([5, 6, 7], [1, 3])
        integer(int32) :: mask(1), subsets(1, n_slots), n_results, ierr

        mask = 7

        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, 3, subsets, n_slots, &
                                  max_angle=0.25_real64, gain_gamma=2.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_angle_bound: 0.25, ierr")
        call check_subsets(1, n_slots, subsets, n_results, 1, narrow, "test_detect_dosage_effect_angle_bound: 0.25")

        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, 3, subsets, n_slots, &
                                  max_angle=PI/4, gain_gamma=2.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_angle_bound: pi/4, ierr")
        call check_subsets(1, n_slots, subsets, n_results, 3, wide, "test_detect_dosage_effect_angle_bound: pi/4")
    end subroutine test_detect_dosage_effect_angle_bound

    !> Without `max_angle` and `gain_gamma`, the documented defaults hold: pi, so no direction is
    !| ruled out, and 0.1. Ancestor [1, 0], genes [-0.5, 0], [-0.5625, 0] and [-0.6171875, 0]
    !| (-79/128), all pointing away from it, at angle pi:
    !| - {1, 2} is 1.0625 long, short of 1.1, and grows; {1, 3} (1.1171875) and {2, 3} (1.1796875)
    !|   are results, and so is {1, 2, 3} (1.6796875): the words 5, 6 and 7.
    !| A default angle below pi would find nothing; a default gain of 0.0625 or less would find
    !| {1, 2}, one above 0.1171875 -- 0.125, say -- would miss {1, 3}.
    subroutine test_detect_dosage_effect_defaults()
        integer(int32), parameter :: n_genes = 3, n_slots = 3
        real(real64), parameter :: ancestor(2) = [1.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, n_genes) = reshape([-0.5_real64, 0.0_real64, -0.5625_real64, 0.0_real64, &
                                                                -0.6171875_real64, 0.0_real64], [2, n_genes])
        integer(int32), parameter :: expected(1, 3) = reshape([5, 6, 7], [1, 3])
        integer(int32) :: mask(1), subsets(1, n_slots), n_results, ierr

        mask = 7
        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, 3, subsets, n_slots, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_defaults: ierr")
        call check_subsets(1, n_slots, subsets, n_results, 3, expected, "test_detect_dosage_effect_defaults")
    end subroutine test_detect_dosage_effect_defaults

    !> Masks of more than one word, in f42's layout:
    !| - kept from the old suite (`test_detect_patterns_last_gene_on_word_boundary`): 100 genes in
    !|   four words, genes 5, 40 and 64 kept -- gene 64 is the top bit of the second word, where
    !|   before #207 the next gene after {5, 64} was taken to be 38. Every gene is 0.375 times the
    !|   ancestor [1, 0, 0], gain 0.1: no pair reaches 1.1, the triple (1.125) does, and is the
    !|   one result, the words [16, 128 + BIT_31, 0, 0];
    !| - 33 genes in two words, genes 32 and 33 kept, each 0.75 times the ancestor [1, 0]: the
    !|   pair is the one result, the words [BIT_31, 1] -- bit 31 of the first word, bit 0 of the
    !|   second.
    subroutine test_detect_dosage_effect_multiword()
        integer(int32), parameter :: n_wide = 100, n_narrow = 33
        integer(int32), parameter :: expected_four(4, 1) = reshape([16_int32, 128_int32 + BIT_31, 0_int32, 0_int32], [4, 1])
        integer(int32), parameter :: expected_two(2, 1) = reshape([BIT_31, 1_int32], [2, 1])
        real(real64) :: wide_ancestor(3), wide_genes(3, n_wide), narrow_ancestor(2), narrow_genes(2, n_narrow)
        integer(int32) :: four_words(4), four_word_subsets(4, 3), two_words(2), two_word_subsets(2, 1), n_results, ierr

        wide_ancestor = [1.0_real64, 0.0_real64, 0.0_real64]
        wide_genes = 0.0_real64
        wide_genes(1, :) = 0.375_real64
        four_words = [16_int32, 128_int32 + BIT_31, 0_int32, 0_int32]
        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(wide_ancestor, wide_genes, n_wide, 3, four_words, 4, n_results, 3, four_word_subsets, 3, &
                                  max_angle=0.0_real64, gain_gamma=0.1_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_multiword: genes 5, 40 and 64, ierr")
        call check_subsets(4, 3, four_word_subsets, n_results, 1, expected_four, &
                           "test_detect_dosage_effect_multiword: genes 5, 40 and 64")

        narrow_ancestor = [1.0_real64, 0.0_real64]
        narrow_genes = 0.0_real64
        narrow_genes(1, :) = 0.75_real64
        two_words = [BIT_31, 1_int32]
        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(narrow_ancestor, narrow_genes, n_narrow, 2, two_words, 2, n_results, 2, two_word_subsets, 1, &
                                  max_angle=0.0_real64, gain_gamma=0.1_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_multiword: genes 32 and 33, ierr")
        call check_subsets(2, 1, two_word_subsets, n_results, 1, expected_two, "test_detect_dosage_effect_multiword: genes 32 and 33")
    end subroutine test_detect_dosage_effect_multiword

    !> `max_subset_size` bounds the search, after issue #64's T7 (16 genes, size 15), smaller:
    !| ten genes [1, 0], the ancestor [8, 0], gain 0.25, so only all ten together reach the
    !| length 10. Size 9 finds nothing, size 10 the word 1023. The work array holds 511 subsets,
    !| 2**9 - 1, the worst case for ten genes.
    subroutine test_detect_dosage_effect_size_bound()
        integer(int32), parameter :: n_genes = 10, n_slots = 511
        integer(int32), parameter :: everything(1, 1) = reshape([1023], [1, 1])
        real(real64) :: ancestor(2), genes(2, n_genes)
        integer(int32) :: mask(1), n_results, ierr
        integer(int32), allocatable :: subsets(:, :)

        ancestor = [8.0_real64, 0.0_real64]
        genes(1, :) = 1.0_real64
        genes(2, :) = 0.0_real64
        mask = 1023
        allocate (subsets(1, n_slots))

        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, 9, subsets, n_slots, &
                                  max_angle=0.0_real64, gain_gamma=0.25_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_size_bound: size 9, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_dosage_effect_size_bound: size 9 finds nothing")

        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, 10, subsets, n_slots, &
                                  max_angle=0.0_real64, gain_gamma=0.25_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_size_bound: size 10, ierr")
        call check_subsets(1, n_slots, subsets, n_results, 1, everything, "test_detect_dosage_effect_size_bound: size 10")
    end subroutine test_detect_dosage_effect_size_bound

    !> Nothing to find, while every gene is four times the ancestor [1, 0], so any pair of them
    !| would be a result:
    !| - no gene kept (the word 0), and only gene 2 kept (the word 2): no pair to form;
    !| - all three kept, `max_subset_size` 1: a single gene is never a result;
    !| - all three kept, `max_subset_size` 0: the early exit, nothing checked;
    !| - a problem of one gene;
    !| - no gene kept and only gene 2 kept, with a work array of no subsets, which is what the
    !|   sizing gives for fewer than two kept genes (and what the Python and R bindings pass on).
    !| `n_results` is 0 each time. The work array holds three subsets, room for the starting genes,
    !| but in the last case.
    subroutine test_detect_dosage_effect_no_candidates()
        integer(int32), parameter :: n_genes = 3, n_slots = 3
        integer(int32), parameter :: n_cases = 4
        real(real64), parameter :: ancestor(2) = [1.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, n_genes) = reshape([4.0_real64, 0.0_real64, 4.0_real64, 0.0_real64, &
                                                                4.0_real64, 0.0_real64], [2, n_genes])
        integer(int32), parameter :: case_masks(n_cases) = [0, 2, 7, 7], case_sizes(n_cases) = [3, 3, 1, 0]
        character(*), parameter :: case_names(n_cases) = ["no gene kept     ", "one gene kept    ", "max_subset_size 1", &
                                                          "max_subset_size 0"]
        integer(int32) :: mask(1), subsets(1, n_slots), n_results, ierr, i_case

        do i_case = 1, n_cases
            mask = case_masks(i_case)
            n_results = UNWRITTEN_COUNT
            call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, case_sizes(i_case), subsets, n_slots, &
                                      ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_no_candidates: "//trim(case_names(i_case))//", ierr")
            call assert_equal_int(n_results, 0_int32, "test_detect_dosage_effect_no_candidates: "//trim(case_names(i_case)))
        end do

        mask = 1
        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(ancestor, genes(:, 1), 1, 2, mask, 1, n_results, 1, subsets, 1, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_no_candidates: one gene, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_dosage_effect_no_candidates: one gene")

        do i_case = 1, 2
            mask = case_masks(i_case)
            n_results = UNWRITTEN_COUNT
            call detect_dosage_effect(ancestor, genes, n_genes, 2, mask, 1, n_results, 3, subsets, 0, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_no_candidates: "//trim(case_names(i_case)) &
                            //", no subsets, ierr")
            call assert_equal_int(n_results, 0_int32, "test_detect_dosage_effect_no_candidates: "//trim(case_names(i_case)) &
                                  //", no subsets")
        end do
    end subroutine test_detect_dosage_effect_no_candidates

    !> The documented chain, filter -> size -> detect, on the mixed family (see the fixtures), after
    !| the old suite's `test_detect_patterns_mixed_results`: the dosage-effect filter at pi/8 keeps
    !| genes 3 and 4 (the word 12), the sizing gives (2, 1), and with gain 0.2 and `max_angle` 0.1
    !| the pair {3, 4}, 1.25 times the ancestor, is the result 12.
    subroutine test_detect_dosage_effect_after_filter()
        integer(int32) :: gene_to_fam(MIXED_N_GENES), family_bitmasks(1, 1), max_subset_size, work_array_size, n_results, ierr
        integer(int32), parameter :: expected(1, 1) = reshape([12], [1, 1])
        integer(int32), allocatable :: subsets(:, :)

        gene_to_fam = 1
        call filter_paralogs_by_pattern_dosage_effect(MIXED_ANGLES, PI/8, MIXED_N_GENES, 1, gene_to_fam, family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_after_filter: filter, ierr")
        call assert_equal_int(family_bitmasks(1, 1), 12_int32, "test_detect_dosage_effect_after_filter: filter keeps 3 and 4")

        max_subset_size = MIXED_N_GENES
        call calc_paralog_subsets_bitmasks_size(max_subset_size, MIXED_N_GENES, work_array_size, family_bitmasks(:, 1), 1, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_after_filter: sizing, ierr")
        call assert_equal_int(max_subset_size, 2_int32, "test_detect_dosage_effect_after_filter: two genes kept")
        call assert_equal_int(work_array_size, 1_int32, "test_detect_dosage_effect_after_filter: one subset")

        allocate (subsets(1, work_array_size))
        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect(MIXED_ANCESTOR, MIXED_GENES, MIXED_N_GENES, MIXED_N_DIMS, family_bitmasks(:, 1), 1, &
                                  n_results, max_subset_size, subsets, work_array_size, max_angle=0.1_real64, &
                                  gain_gamma=0.2_real64, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_after_filter: detection, ierr")
        call check_subsets(1, work_array_size, subsets, n_results, 1, expected, "test_detect_dosage_effect_after_filter")
    end subroutine test_detect_dosage_effect_after_filter

    !> The two errors the detection itself finds, past the validation:
    !| - a kept bit past `n_genes` (gene 4 of 3, the word 9) is ERR_INVALID_INPUT; it names no
    !|   argument today, though the mask is at fault, so only the code is checked;
    !| - a work array too small for the search (5 subsets for the binary-weights family, which
    !|   needs 15) is ERR_SIZE_MISMATCH, at no argument: the internal bookkeeping finds it, and
    !|   its positions number its own arguments, not the caller's.
    subroutine test_detect_dosage_effect_early_errors()
        integer(int32) :: mask(1), subsets(1, 5), n_results, ierr

        mask = 9
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 5, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_early_errors: a bit past n_genes")

        mask = 31
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 5, subsets, 5, &
                                  max_angle=0.0_real64, gain_gamma=1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_detect_dosage_effect_early_errors: a work array too small", 0)
    end subroutine test_detect_dosage_effect_early_errors

    !> The extents, each refused with its position:
    !| - `n_genes` (argument 3) and `n_dims` (argument 4): ERR_EMPTY_INPUT at 0, ERR_INVALID_INPUT
    !|   at -1;
    !| - `n_mask_chunks` (argument 6) must hold `n_genes` bits: 0 words for 3 genes, and 1 word
    !|   for 33, are ERR_INVALID_INPUT;
    !| - `max_subset_size` (argument 8) is ERR_INVALID_INPUT at -1 (0 is the early exit in
    !|   `test_detect_dosage_effect_no_candidates`);
    !| - `n_paralog_subsets` (argument 10): ERR_INVALID_INPUT at -1; 0 passes where nothing is
    !|   to be stored (`test_detect_dosage_effect_no_candidates`) and is ERR_SIZE_MISMATCH where
    !|   something is (`test_detect_dosage_effect_starts_in_the_array`).
    subroutine test_detect_dosage_effect_dimensions()
        real(real64) :: genes_33(2, 33)
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr

        mask = 3
        genes_33 = 1.0_real64

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 0, 2, mask, 1, n_results, 2, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_dosage_effect_dimensions: n_genes 0", 3)
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, -1, 2, mask, 1, n_results, 2, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_dimensions: n_genes -1", 3)

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 0, mask, 1, n_results, 2, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_dosage_effect_dimensions: n_dims 0", 4)
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, -1, mask, 1, n_results, 2, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_dimensions: n_dims -1", 4)

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 0, n_results, 2, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_dimensions: no word for 3 genes", 6)
        call detect_dosage_effect(BINARY_ANCESTOR, genes_33, 33, 2, mask, 1, n_results, 2, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_dimensions: one word for 33 genes", 6)

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, -1, subsets, 3, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_dimensions: max_subset_size -1", 8)

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 2, subsets, -1, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_dimensions: n_paralog_subsets -1", 10)
    end subroutine test_detect_dosage_effect_dimensions

    !> The documented ranges of the thresholds, kept from the old suite's
    !| `test_detect_patterns_input_validation` and completed:
    !| - `max_angle` (argument 11) runs over [0, pi]: 0 passes in
    !|   `test_detect_dosage_effect_binary_weights`, pi here; -epsilon and the next number above pi
    !|   are ERR_INVALID_INPUT;
    !| - `gain_gamma` (argument 12) must be above 0: `tiny` passes, 0 and -1 are ERR_INVALID_INPUT.
    subroutine test_detect_dosage_effect_value_ranges()
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr

        mask = 7

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, max_angle=PI, ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_value_ranges: max_angle pi")
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                  max_angle=-epsilon(1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_value_ranges: max_angle below 0", 11)
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                  max_angle=nearest(PI, 4.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_value_ranges: max_angle above pi", 11)

        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                  gain_gamma=tiny(1.0_real64), ierr=ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_value_ranges: gain_gamma tiny")
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                  gain_gamma=0.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_value_ranges: gain_gamma 0", 12)
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                  gain_gamma=-1.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_value_ranges: gain_gamma -1", 12)
    end subroutine test_detect_dosage_effect_value_ranges

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `ancestor` (argument 1) and `genes` (argument 2), each
    !| in the last element, and in `max_angle` (argument 11) and `gain_gamma` (argument 12). Kept
    !| from the old suite's `test_detect_patterns_input_validation` (+Inf in the ancestor and in a
    !| gene), now with every non-finite value.
    subroutine test_detect_dosage_effect_rejects_nan_and_inf()
        real(real64) :: ancestor(2), genes(2, 3)
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr, i_value
        character(:), allocatable :: value_name

        mask = 7
        do i_value = 1, N_NON_FINITE
            value_name = trim(NON_FINITE_NAMES(i_value))

            ancestor = BINARY_ANCESTOR
            genes = BINARY_GENES(:, 1:3)
            ancestor(2) = non_finite_value(i_value)
            call detect_dosage_effect(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, ierr=ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_dosage_effect_rejects_nan_and_inf: ancestor "//value_name, 1)

            ancestor = BINARY_ANCESTOR
            genes(2, 3) = non_finite_value(i_value)
            call detect_dosage_effect(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, ierr=ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_dosage_effect_rejects_nan_and_inf: genes "//value_name, 2)

            genes = BINARY_GENES(:, 1:3)
            call detect_dosage_effect(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                      max_angle=non_finite_value(i_value), ierr=ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_dosage_effect_rejects_nan_and_inf: max_angle "//value_name, 11)
            call detect_dosage_effect(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                      gain_gamma=non_finite_value(i_value), ierr=ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_dosage_effect_rejects_nan_and_inf: gain_gamma "//value_name, 12)
        end do
    end subroutine test_detect_dosage_effect_rejects_nan_and_inf

    !> The expert entry takes the two scratch arrays from the caller and allocates nothing: it
    !| finds the same 15 subsets in the binary-weights family, and names its arguments by its own
    !| list, where `max_angle` is 13 and `gain_gamma` 14 (kept from the old suite's
    !| `test_detect_patterns_input_validation`).
    subroutine test_detect_dosage_effect_expert()
        integer(int32) :: mask(1), subsets(1, BINARY_N_SLOTS), tmp_active_bitmask(1), n_results, ierr
        real(real64) :: tmp_paralog_vector(2)

        mask = 31
        n_results = UNWRITTEN_COUNT
        call detect_dosage_effect_expert(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 5, subsets, &
                                         BINARY_N_SLOTS, tmp_active_bitmask, tmp_paralog_vector, 0.0_real64, 1.0_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_expert: ierr")
        call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, BINARY_N_RESULTS, BINARY_RESULTS, "test_detect_dosage_effect_expert")

        call detect_dosage_effect_expert(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 5, subsets, &
                                         BINARY_N_SLOTS, tmp_active_bitmask, tmp_paralog_vector, nearest(PI, 4.0_real64), &
                                         1.0_real64, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_expert: max_angle above pi", 13)
        call detect_dosage_effect_expert(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 5, subsets, &
                                         BINARY_N_SLOTS, tmp_active_bitmask, tmp_paralog_vector, 0.0_real64, 0.0_real64, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_dosage_effect_expert: gain_gamma 0", 14)
    end subroutine test_detect_dosage_effect_expert

    !> BUG (fails today): the binary-weights family scaled by a power of two finds the same 15
    !| subsets. Scaling every vector by 2**k changes no angle and no ratio of lengths, and is exact,
    !| so the result may not change; issue #64's T12 asks for this scale invariance. The angle of a
    !| subset's sum to the ancestor is taken from unscaled sums of squares, and a product of the two
    !| lengths of at most 1e-12 counts as zero (f42's `angle_between`, issue #223):
    !| - by 2**-30 (lengths about 1e-8) that product is about 1e-16 times the word: ERR_DIVISION_BY_ZERO,
    !|   at no argument, for the whole family;
    !| - by 2**520 (lengths about 1e157) the squares pass `huge` and the cosine is Inf/Inf, NaN:
    !|   gfortran at -O0 (--diagnostics) then drops every subset and finds nothing, while gfortran
    !|   at -O3 and ifx happen to find all 15. The outcome is the build's, not the method's.
    subroutine test_detect_dosage_effect_any_scale()
        integer(int32), parameter :: n_scales = 2
        real(real64), parameter :: scales(n_scales) = [2.0_real64**(-30), 2.0_real64**520]
        character(*), parameter :: scale_names(n_scales) = ["2**-30", "2**520"]
        real(real64) :: ancestor(2), genes(2, BINARY_N_GENES)
        integer(int32) :: mask(1), subsets(1, BINARY_N_SLOTS), n_results, ierr, i_scale

        mask = 31
        do i_scale = 1, n_scales
            ancestor = BINARY_ANCESTOR*scales(i_scale)
            genes = BINARY_GENES*scales(i_scale)
            n_results = UNWRITTEN_COUNT
            call detect_dosage_effect(ancestor, genes, BINARY_N_GENES, 2, mask, 1, n_results, 5, subsets, BINARY_N_SLOTS, &
                                      max_angle=0.0_real64, gain_gamma=1.0_real64, ierr=ierr)
            call assert_err(ierr, ERR_OK, "test_detect_dosage_effect_any_scale: by "//scale_names(i_scale)//", ierr")
            call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, BINARY_N_RESULTS, BINARY_RESULTS, &
                               "test_detect_dosage_effect_any_scale: by "//scale_names(i_scale))
        end do
    end subroutine test_detect_dosage_effect_any_scale

    !> f42_vector's `angle_between`, which the dosage effect's angle bound is built on and which
    !| has no suite of its own yet; it moves to the f42 suites with the f42 batch. Kept from the old
    !| suite's `test_angle_between`, with values written by hand:
    !| - [1, 0] to itself and to [2, 0]: the cosine is exactly 1, the angle exactly 0;
    !| - [1, 0] to [-2, 0]: the cosine is exactly -1, the angle pi (the dosage effect's default
    !|   bound relies on acos(-1) being `PI` exactly, see `test_detect_dosage_effect_defaults`);
    !| - [1, 0] to [0, 3]: the cosine is exactly 0, the angle pi/2; and [1, 0] to [1, 1]: the
    !|   cosine is 1/sqrt(2), rounded, the angle pi/4. Both within 2 ulps: acos comes from the
    !|   compiler's library, which need not round correctly, and 1/sqrt(2) is itself rounded;
    !| - scaling both vectors by 2**20 or 2**-20 changes no angle, to the last bit: every product
    !|   and root scales exactly ([3, 4] to [5, 12], an angle of about 0.284).
    !| The old case's 45-degree pair [0.7071, 0.7071] matched pi/4 only to 1e-5, under a tolerance
    !| of epsilon; zero and huge vectors are #223's (see `test_detect_dosage_effect_any_scale`).
    subroutine test_angle_between()
        real(real64), parameter :: TWO_ULPS = 2*epsilon(1.0_real64)
        real(real64), parameter :: scales(2) = [2.0_real64**20, 2.0_real64**(-20)]
        real(real64) :: x_axis(2), other(2), angle, unscaled_angle, scaled_angle
        integer(int32) :: ierr, i_scale

        x_axis = [1.0_real64, 0.0_real64]

        other = x_axis
        call angle_between(x_axis, other, 2, angle, ierr)
        call assert_err(ierr, ERR_OK, "test_angle_between: to itself, ierr")
        call assert_equal_real(angle, 0.0_real64, 0.0_real64, "test_angle_between: to itself")
        other = [2.0_real64, 0.0_real64]
        call angle_between(x_axis, other, 2, angle, ierr)
        call assert_equal_real(angle, 0.0_real64, 0.0_real64, "test_angle_between: to twice itself")

        other = [-2.0_real64, 0.0_real64]
        call angle_between(x_axis, other, 2, angle, ierr)
        call assert_err(ierr, ERR_OK, "test_angle_between: opposite, ierr")
        call assert_equal_real(angle, PI, 0.0_real64, "test_angle_between: opposite")

        other = [0.0_real64, 3.0_real64]
        call angle_between(x_axis, other, 2, angle, ierr)
        call assert_equal_real(angle, PI/2, TWO_ULPS, "test_angle_between: perpendicular")
        other = [1.0_real64, 1.0_real64]
        call angle_between(x_axis, other, 2, angle, ierr)
        call assert_equal_real(angle, PI/4, TWO_ULPS, "test_angle_between: 45 degrees")

        x_axis = [3.0_real64, 4.0_real64]
        other = [5.0_real64, 12.0_real64]
        call angle_between(x_axis, other, 2, unscaled_angle, ierr)
        do i_scale = 1, 2
            x_axis = [3.0_real64, 4.0_real64]*scales(i_scale)
            other = [5.0_real64, 12.0_real64]*scales(i_scale)
            call angle_between(x_axis, other, 2, scaled_angle, ierr)
            call assert_err(ierr, ERR_OK, "test_angle_between: scaled, ierr")
            call assert_equal_real(scaled_angle, unscaled_angle, 0.0_real64, "test_angle_between: scaled by "//scales(i_scale))
        end do
    end subroutine test_angle_between

    !> The detection starts from each kept gene but the last, one subset each, and these must fit
    !| the work array too: three kept genes start two subsets, so a work array of one is
    !| ERR_SIZE_MISMATCH, naming `n_paralog_subsets` (argument 10), before anything is written;
    !| so is a work array of none.
    !| The caller's array has a second column the call is not given, holding 12345; it must stay
    !| so. (The starting subsets used to be written unchecked: past the array's end, into that
    !| column, or a bounds-check abort in a checking build.)
    subroutine test_detect_dosage_effect_starts_in_the_array()
        integer(int32) :: mask(1), subsets(1, 2), n_results, ierr

        mask = 7
        subsets = 12345
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 1, &
                                  max_angle=0.0_real64, gain_gamma=100.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_detect_dosage_effect_starts_in_the_array: refused", 10)
        call assert_equal_int(subsets(1, 2), 12345_int32, "test_detect_dosage_effect_starts_in_the_array: column 2 untouched")

        subsets = 12345
        call detect_dosage_effect(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 0, &
                                  max_angle=0.0_real64, gain_gamma=100.0_real64, ierr=ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_detect_dosage_effect_starts_in_the_array: no subsets", 10)
        call assert_equal_int(subsets(1, 1), 12345_int32, "test_detect_dosage_effect_starts_in_the_array: column 1 untouched")
    end subroutine test_detect_dosage_effect_starts_in_the_array

end module mod_test_detect_dosage_effect
