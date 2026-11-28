!> The `calc_paralog_subsets_bitmasks_size` cases: how many subsets a detection's work array
!| must hold, and the largest subset size it may search, for the genes a filter kept.
!|
!| The expected sizes are counted by hand, not with the routine's own loop. With F genes kept and
!| subsets of at most m genes (m is `max_subset_size`, cut to F), the worst case keeps every
!| subset of the first F - 1 kept genes as a candidate and turns every subset that adds the last
!| kept gene into a result (`test_calc_size_enough_for_worst_case` builds exactly that). At the
!| end of size m the array holds, side by side, the results of every size, the subsets of the
!| first F - 1 genes with 2 to m genes plus the last one, and the m-gene candidates without it:
!| that is one subset of the first F - 1 genes for each of sizes 1 to m, so the size is
!|
!|     sum over k = 1 .. min(m, F - 1) of C(F - 1, k),
!|
!| which is 2**(F - 1) - 1 when m = F. So F = 0 and F = 1 need no array, F = 2 needs 1, F = 5
!| needs 15. With `max_subset_size` 0 there is nothing to size, and the answer is 0.
module mod_test_calc_paralog_subsets_bitmasks_size
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_paralog_analysis, only: calc_paralog_subsets_bitmasks_size, detect_subfunctionalization_expert
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_paralog_analysis_fixtures, only: BIT_31, UNWRITTEN_COUNT
    implicit none
    private
    public :: get_all_tests_calc_paralog_subsets_bitmasks_size

contains

    !> Every case of `calc_paralog_subsets_bitmasks_size`.
    function get_all_tests_calc_paralog_subsets_bitmasks_size() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(7))
        all_tests(1) = test_case("test_calc_size_small_families", test_calc_size_small_families)
        all_tests(2) = test_case("test_calc_size_zero_subset_size", test_calc_size_zero_subset_size)
        all_tests(3) = test_case("test_calc_size_counts_kept_bits", test_calc_size_counts_kept_bits)
        all_tests(4) = test_case("test_calc_size_largest_uncapped", test_calc_size_largest_uncapped)
        all_tests(5) = test_case("test_calc_size_enough_for_worst_case", test_calc_size_enough_for_worst_case)
        all_tests(6) = test_case("test_calc_size_bounds", test_calc_size_bounds)
        all_tests(7) = test_case("test_calc_size_caps_at_largest_fit", test_calc_size_caps_at_largest_fit)
    end function get_all_tests_calc_paralog_subsets_bitmasks_size

    !> Calls the routine and checks `ierr` OK, the subset size it hands back and the work array size.
    !| Both outputs start at UNWRITTEN_COUNT or the request, so an output never written shows.
    subroutine check_size(n_genes, n_mask_chunks, filtered_paralogs_bitmask, requested, expected_subset_size, &
                          expected_work_array_size, msg)
        integer(int32), intent(in) :: n_genes, n_mask_chunks, requested, expected_subset_size, expected_work_array_size
        integer(int32), intent(in) :: filtered_paralogs_bitmask(n_mask_chunks)
        character(*), intent(in) :: msg
        integer(int32) :: max_subset_size, work_array_size, ierr

        max_subset_size = requested
        work_array_size = UNWRITTEN_COUNT
        call calc_paralog_subsets_bitmasks_size(max_subset_size, n_genes, work_array_size, filtered_paralogs_bitmask, &
                                                n_mask_chunks, ierr)
        call assert_err(ierr, ERR_OK, msg//": ierr")
        call assert_equal_int(max_subset_size, expected_subset_size, msg//": max_subset_size")
        call assert_equal_int(work_array_size, expected_work_array_size, msg//": work_array_size")
    end subroutine check_size

    !> Five genes, F of them kept, any subset size (100 is cut to F):
    !| - F = 0 (the word 0): no subset, (0, 0); F = 1 (gene 5 alone, 16): a single gene starts no
    !|   pair, (1, 0); F = 2 (genes 2 and 5, 18): (2, 1); F = 3 (genes 1 to 3, 7): (3, 3); F = 4
    !|   (15): (4, 7); F = 5 (31): (5, 15).
    !| - All five kept, m = 1 to 4: C(4, 1) = 4, + C(4, 2) = 10, + C(4, 3) = 14, + C(4, 4) = 15.
    !| In place of the old Python and R check of three kept genes at size 3 (work array size 3).
    subroutine test_calc_size_small_families()
        integer(int32), parameter :: n_genes = 5, n_cases = 6
        integer(int32) :: masks(1, n_cases), expected_sizes(n_cases), expected_partial(4), i_case, subset_size

        masks = reshape([0, 16, 18, 7, 15, 31], [1, n_cases])
        expected_sizes = [0, 0, 1, 3, 7, 15]
        do i_case = 1, n_cases
            call check_size(n_genes, 1, masks(:, i_case), 100, i_case - 1, expected_sizes(i_case), &
                            "test_calc_size_small_families: F = "//(i_case - 1))
        end do

        expected_partial = [4, 10, 14, 15]
        do subset_size = 1, 4
            call check_size(n_genes, 1, masks(:, n_cases), subset_size, subset_size, expected_partial(subset_size), &
                            "test_calc_size_small_families: F = 5, m = "//subset_size)
        end do
    end subroutine test_calc_size_small_families

    !> `max_subset_size` 0 is in range and the early exit: nothing to size, the work array size is
    !| 0 and the subset size stays 0, however many genes are kept (all five here).
    subroutine test_calc_size_zero_subset_size()
        integer(int32) :: mask(1)

        mask = 31
        call check_size(5, 1, mask, 0, 0, 0, "test_calc_size_zero_subset_size")
    end subroutine test_calc_size_zero_subset_size

    !> Only the number of kept genes counts, wherever they are, in any word:
    !| - genes 2 and 5 of 5 are F = 2, as genes 1 and 2 would be: (2, 1);
    !| - genes 1, 32 and 40 of 40, in two words [1 + BIT_31, 128] -- gene 32 is bit 31 of the
    !|   first word, gene 40 bit 7 of the second -- are F = 3: (3, 3).
    subroutine test_calc_size_counts_kept_bits()
        integer(int32) :: one_word(1), two_words(2)

        one_word = 18
        call check_size(5, 1, one_word, 5, 2, 1, "test_calc_size_counts_kept_bits: genes 2 and 5")

        two_words = [1_int32 + BIT_31, 128_int32]
        call check_size(40, 2, two_words, 40, 3, 3, "test_calc_size_counts_kept_bits: genes 1, 32 and 40")
    end subroutine test_calc_size_counts_kept_bits

    !> All 29 genes kept (the word 2**29 - 1), any size: 2**28 - 1 = 268,435,455 subsets, not
    !| allocated here. Until the sizing counted in int64 this was the largest family it sized
    !| without cutting the subset size; families up to 32 genes fit whole
    !| (`test_calc_size_caps_at_largest_fit`).
    subroutine test_calc_size_largest_uncapped()
        integer(int32) :: mask(1)

        mask = 536870911
        call check_size(29, 1, mask, 29, 29, 268435455, "test_calc_size_largest_uncapped: F = 29")
    end subroutine test_calc_size_largest_uncapped

    !> In place of the old suite's `test_calc_paralog_subsets_bitmasks_size` (which compared the
    !| slots a detection used with the size), and of the old helper cases (`take_active_mask_helper`,
    !| `add_to_results_helper`, `add_new_active_mask_helper`), which this drives through the
    !| published routine: the worst case fills the work array to the last slot, and one slot fewer
    !| is ERR_SIZE_MISMATCH.
    !|
    !| F kept genes, F = 3 to 6, subset sizes m = 2 to F. The ancestor is [3, 4], genes 1 to F - 1
    !| are zero and gene F is the ancestor, so a subset with gene F has residual 0 and is a result
    !| (threshold 0), and one without it has residual 5, which a later gene can always lower: the
    !| smallest norm after it is 0, or 5 when only gene F follows, never above the residual. So no
    !| subset is pruned, as the size assumes. The results are the subsets of the first F - 1 genes
    !| with 1 to m - 1 genes, plus gene F: sum over k = 1 .. m - 1 of C(F - 1, k).
    !|
    !| m = 1 is left out: one slot fewer cannot hold the F - 1 starting genes, and the detection
    !| writes those without checking the array's size, past its end.
    subroutine test_calc_size_enough_for_worst_case()
        integer(int32), parameter :: n_dims = 2, max_genes = 6
        integer(int32), parameter :: binomial(0:max_genes - 1, 0:max_genes - 1) = reshape([ &
                                     1, 1, 1, 1, 1, 1, &
                                     0, 1, 2, 3, 4, 5, &
                                     0, 0, 1, 3, 6, 10, &
                                     0, 0, 0, 1, 4, 10, &
                                     0, 0, 0, 0, 1, 5, &
                                     0, 0, 0, 0, 0, 1], [max_genes, max_genes])
            !! binomial(n, k), C(n, k)
        real(real64) :: ancestor(n_dims), genes(n_dims, max_genes), paralog_norms(max_genes), tmp_work_array(max_genes)
        real(real64) :: tmp_paralog_vector(n_dims)
        integer(int32) :: sorted_paralog_norms_perm(max_genes), mask(1), tmp_active_bitmask(1)
        integer(int32), allocatable :: subsets(:, :)
        integer(int32) :: n_genes, subset_size, max_subset_size, work_array_size, n_results, expected_results, ierr, k
        character(:), allocatable :: label

        ancestor = [3.0_real64, 4.0_real64]
        do n_genes = 3, max_genes
            genes = 0.0_real64
            genes(:, n_genes) = ancestor
            paralog_norms = 0.0_real64
            paralog_norms(n_genes) = 5.0_real64
            ! ascending: the zero genes in index order, then gene F
            do k = 1, n_genes
                sorted_paralog_norms_perm(k) = k
            end do
            mask = 2**n_genes - 1

            do subset_size = 2, n_genes
                label = "test_calc_size_enough_for_worst_case: F = "//n_genes//", m = "//subset_size
                max_subset_size = subset_size
                call calc_paralog_subsets_bitmasks_size(max_subset_size, n_genes, work_array_size, mask, 1, ierr)
                call assert_err(ierr, ERR_OK, label//": sizing")

                expected_results = 0
                do k = 1, subset_size - 1
                    expected_results = expected_results + binomial(n_genes - 1, k)
                end do

                allocate (subsets(1, work_array_size))
                n_results = UNWRITTEN_COUNT
                call detect_subfunctionalization_expert(ancestor, genes, n_genes, n_dims, mask, 1, n_results, max_subset_size, &
                                                        subsets, work_array_size, tmp_active_bitmask, tmp_paralog_vector, &
                                                        0.0_real64, paralog_norms, sorted_paralog_norms_perm, tmp_work_array, &
                                                        ierr)
                call assert_err(ierr, ERR_OK, label//": the size suffices")
                call assert_equal_int(n_results, expected_results, label//": n_results")
                deallocate (subsets)

                allocate (subsets(1, work_array_size - 1))
                call detect_subfunctionalization_expert(ancestor, genes, n_genes, n_dims, mask, 1, n_results, max_subset_size, &
                                                        subsets, work_array_size - 1, tmp_active_bitmask, tmp_paralog_vector, &
                                                        0.0_real64, paralog_norms, sorted_paralog_norms_perm, tmp_work_array, &
                                                        ierr)
                call assert_err(ierr, ERR_SIZE_MISMATCH, label//": one slot fewer", 0)
                deallocate (subsets)
            end do
        end do
    end subroutine test_calc_size_enough_for_worst_case

    !> The extents:
    !| - `max_subset_size` (argument 1) is ERR_INVALID_INPUT at -1;
    !| - `n_genes` (argument 2) runs from 1 to 32 bits per word: 32 genes in one word pass, 33
    !|   are ERR_INVALID_INPUT, and so are 0 and -1;
    !| - `n_mask_chunks` (argument 5) at 0 or -1 leaves no bit for even one gene, so it is
    !|   `n_genes` (argument 2) that is refused first, as ERR_INVALID_INPUT; the extent's own check
    !|   comes after it and cannot be reached with a valid `n_genes`.
    subroutine test_calc_size_bounds()
        integer(int32) :: mask(1), max_subset_size, work_array_size, ierr

        mask = 1

        max_subset_size = -1
        call calc_paralog_subsets_bitmasks_size(max_subset_size, 1, work_array_size, mask, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_size_bounds: max_subset_size -1", 1)

        call check_size(32, 1, mask, 2, 1, 0, "test_calc_size_bounds: 32 genes in one word")

        max_subset_size = 2
        call calc_paralog_subsets_bitmasks_size(max_subset_size, 33, work_array_size, mask, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_size_bounds: 33 genes in one word", 2)
        call calc_paralog_subsets_bitmasks_size(max_subset_size, 0, work_array_size, mask, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_size_bounds: n_genes 0", 2)
        call calc_paralog_subsets_bitmasks_size(max_subset_size, -1, work_array_size, mask, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_size_bounds: n_genes -1", 2)

        call calc_paralog_subsets_bitmasks_size(max_subset_size, 1, work_array_size, mask, 0, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_size_bounds: n_mask_chunks 0", 2)
        call calc_paralog_subsets_bitmasks_size(max_subset_size, 1, work_array_size, mask, -1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_size_bounds: n_mask_chunks -1", 2)
    end subroutine test_calc_size_bounds

    !> Regression: the subset size is cut to "the maximum valid size", the largest whose work
    !| array holds at most `huge` subsets. It used to be cut where an intermediate product of the
    !| routine's count, C(F, s - 1)*(F - s + 1), passed `huge`, far below. The largest valid sizes,
    !| by the count in this module's header (int32's `huge` is 2,147,483,647), and what came back
    !| before the sizing counted in int64:
    !| - F = 30, m = 30: 2**29 - 1 = 536,870,911 fits; it used to give (14, 268,435,455);
    !| - F = 32, m = 32: 2**31 - 1, `huge` itself, fits; it used to give (11, 160,645,503);
    !| - F = 33, m = 33: 15 is the largest size that fits, 1,846,943,452 (16 would be
    !|   2,448,023,842); it used to give (11, 236,618,692);
    !| - F = 100, m = 16 (the old suite's overflow case, which only checked that the size changed):
    !|   6 fits, 1,195,978,575 (7 would be 16,083,010,119); it used to give (5, 75,449,319).
    subroutine test_calc_size_caps_at_largest_fit()
        integer(int32) :: one_word(1), two_words(2), four_words(4)

        one_word = 1073741823
        call check_size(30, 1, one_word, 30, 30, 536870911, "test_calc_size_caps_at_largest_fit: F = 30")

        one_word = -1
        call check_size(32, 1, one_word, 32, 32, huge(0_int32), "test_calc_size_caps_at_largest_fit: F = 32")

        two_words = [-1, 1]
        call check_size(33, 2, two_words, 33, 15, 1846943452, "test_calc_size_caps_at_largest_fit: F = 33")

        four_words = [-1, -1, -1, 15]
        call check_size(100, 4, four_words, 16, 6, 1195978575, "test_calc_size_caps_at_largest_fit: F = 100")
    end subroutine test_calc_size_caps_at_largest_fit

end module mod_test_calc_paralog_subsets_bitmasks_size
