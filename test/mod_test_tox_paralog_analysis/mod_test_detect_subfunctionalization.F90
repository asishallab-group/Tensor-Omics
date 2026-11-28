!> The `detect_subfunctionalization` cases, and its expert entry `detect_subfunctionalization_expert`.
!|
!| The detection grows subsets of the kept genes one gene at a time, in index order, as the dosage
!| effect's does (see `mod_test_detect_dosage_effect`). A subset of two or more genes is a result
!| when its residual, the ancestor minus the genes' sum, is at most `rdi_threshold` long, and is
!| not grown further. Otherwise it grows only if a gene follows its last one and the smallest
!| norm among the genes after its last one (all of them, kept or not, by `paralog_norms`) is at
!| most the residual's length; issue #64's pruning rule.
!|
!| The genes in most cases lie along one axis, with lengths that are small integers or multiples
!| of a power of two, so every residual and every comparison is exact. The expected subsets are
!| bit mask words (bit `i - 1` for gene `i`), compared in any order.
module mod_test_detect_subfunctionalization
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_paralog_analysis, only: detect_subfunctionalization, detect_subfunctionalization_expert, &
                                    filter_paralogs_by_pattern_subfunctionalization, calc_paralog_subsets_bitmasks_size
    use tox_paralog_analysis_impl, only: detect_patterns_impl, MODE_SUBFUNC_PATTERN
    use tox_errors
    use test_suite, only: test_case
    use f42_math_impl, only: PI
    use mod_test_tox_paralog_analysis_fixtures, only: UNWRITTEN_COUNT, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value, &
                                                      check_subsets, MIXED_N_DIMS, MIXED_N_GENES, MIXED_ANCESTOR, &
                                                      MIXED_GENES, MIXED_ANGLES
    implicit none
    private
    public :: get_all_tests_detect_subfunctionalization

    !> The binary-weights family of `test_detect_subfunc_binary_weights`: gene i is [2**(i - 1), 0],
    !| so a subset's sum is [w, 0] for its bit mask word w; the ancestor is [10, 0], and a subset's
    !| residual is |10 - w| long.
    integer(int32), parameter :: BINARY_N_GENES = 5, BINARY_N_SLOTS = 15
    real(real64), parameter :: BINARY_ANCESTOR(2) = [10.0_real64, 0.0_real64]
    real(real64), parameter :: BINARY_GENES(2, BINARY_N_GENES) = reshape([ &
                               1.0_real64, 0.0_real64, 2.0_real64, 0.0_real64, 4.0_real64, 0.0_real64, &
                               8.0_real64, 0.0_real64, 16.0_real64, 0.0_real64], [2, BINARY_N_GENES])
    real(real64), parameter :: BINARY_NORMS(BINARY_N_GENES) = [1.0_real64, 2.0_real64, 4.0_real64, 8.0_real64, 16.0_real64]
    integer(int32), parameter :: BINARY_PERM(BINARY_N_GENES) = [1, 2, 3, 4, 5]
    !> Its results at `rdi_threshold` 1 with subsets of three or more genes, pairs first.
    integer(int32), parameter :: BINARY_N_RESULTS = 3
    integer(int32), parameter :: BINARY_RESULTS(1, BINARY_N_RESULTS) = reshape([9, 10, 11], [1, BINARY_N_RESULTS])

contains

    !> The expert tier with the three scratch arrays it needs made here, so a case reads like a call
    !| of the plain tier with `paralog_norms` and `paralog_norms_perm` added.
    subroutine subfunc_expert(ancestor, genes, n_genes, n_dims, mask, n_mask_chunks, n_results, max_subset_size, subsets, &
                              n_slots, rdi_threshold, norms, perm, ierr)
        integer(int32), intent(in) :: n_genes, n_dims, n_mask_chunks, max_subset_size, n_slots
        real(real64), intent(in) :: ancestor(n_dims), genes(n_dims, n_genes), rdi_threshold, norms(n_genes)
        integer(int32), intent(in) :: mask(n_mask_chunks), perm(n_genes)
        integer(int32), intent(out) :: n_results, subsets(n_mask_chunks, n_slots), ierr
        integer(int32) :: tmp_active_bitmask(n_mask_chunks)
        real(real64) :: tmp_paralog_vector(n_dims), tmp_work_array(n_genes)

        call detect_subfunctionalization_expert(ancestor, genes, n_genes, n_dims, mask, n_mask_chunks, n_results, &
                                                max_subset_size, subsets, n_slots, tmp_active_bitmask, tmp_paralog_vector, &
                                                rdi_threshold, norms, perm, tmp_work_array, ierr)
    end subroutine subfunc_expert

    !> Every case of `detect_subfunctionalization`.
    function get_all_tests_detect_subfunctionalization() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(14))
        all_tests(1) = test_case("test_detect_subfunc_binary_weights", test_detect_subfunc_binary_weights)
        all_tests(2) = test_case("test_detect_subfunc_threshold_inclusive", test_detect_subfunc_threshold_inclusive)
        all_tests(3) = test_case("test_detect_subfunc_pruning_bound", test_detect_subfunc_pruning_bound)
        all_tests(4) = test_case("test_detect_subfunc_suffix_minimum", test_detect_subfunc_suffix_minimum)
        all_tests(5) = test_case("test_detect_subfunc_no_candidates", test_detect_subfunc_no_candidates)
        all_tests(6) = test_case("test_detect_subfunc_after_filter", test_detect_subfunc_after_filter)
        all_tests(7) = test_case("test_detect_subfunc_early_errors", test_detect_subfunc_early_errors)
        all_tests(8) = test_case("test_detect_subfunc_dimensions", test_detect_subfunc_dimensions)
        all_tests(9) = test_case("test_detect_subfunc_value_ranges", test_detect_subfunc_value_ranges)
        all_tests(10) = test_case("test_detect_subfunc_rejects_nan_and_inf", test_detect_subfunc_rejects_nan_and_inf)
        all_tests(11) = test_case("test_detect_subfunc_expert", test_detect_subfunc_expert)
        all_tests(12) = test_case("test_detect_subfunc_order_is_the_callers", test_detect_subfunc_order_is_the_callers)
        all_tests(13) = test_case("test_detect_subfunc_suffix_minima", test_detect_subfunc_suffix_minima)
        all_tests(14) = test_case("test_detect_subfunc_plain_matches_expert", test_detect_subfunc_plain_matches_expert)
    end function get_all_tests_detect_subfunctionalization

    !> The binary-weights family (see BINARY_GENES), all five genes kept, `rdi_threshold` 1; the
    !| plain tier computes the genes' norms and their order. The smallest norm after gene i is gene
    !| i + 1's, 2**i:
    !| - pairs: {1, 4} (9, residual 1) and {2, 4} (10, residual 0) are results; {1, 2} (3,
    !|   residual 7) grows, as gene 3's norm 4 is at most 7; {1, 3} (residual 5), {2, 3} (4) and
    !|   {3, 4} (2) are pruned, the norms after them (8, 8, 16) being longer; pairs with gene 5
    !|   have no gene after them;
    !| - triples: {1, 2, 4} (11, residual 1) is a result; {1, 2, 3} (residual 3) is pruned by
    !|   gene 4's norm 8.
    !| `max_subset_size` 2 finds the words 9 and 10, and 3 to 5 also 11. As in the dosage effect,
    !| {2, 4} and its superset {1, 2, 4} are both results, and another gene order finds others
    !| (issue #244).
    subroutine test_detect_subfunc_binary_weights()
        integer(int32), parameter :: n_expected(5) = [0, 2, 3, 3, 3]
        integer(int32) :: mask(1), subsets(1, BINARY_N_SLOTS), n_results, ierr, max_subset_size
        character(:), allocatable :: label

        mask = 31
        do max_subset_size = 2, 5
            label = "test_detect_subfunc_binary_weights: max_subset_size "//max_subset_size
            n_results = UNWRITTEN_COUNT
            call detect_subfunctionalization(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, &
                                             max_subset_size, subsets, BINARY_N_SLOTS, 1.0_real64, ierr)
            call assert_err(ierr, ERR_OK, label//": ierr")
            call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, n_expected(max_subset_size), &
                               BINARY_RESULTS(:, 1:n_expected(max_subset_size)), label)
        end do
    end subroutine test_detect_subfunc_binary_weights

    !> A residual of exactly `rdi_threshold` is a result:
    !| - ancestor [1.75, 2], genes [1, 0] and [0, 1]: the pair's residual [0.75, 1] is 1.25 long;
    !|   `rdi_threshold` 1.25 finds the word 3, the next number below 1.25 nothing;
    !| - kept from the old suite (`test_detect_patterns_perfect_subfunc_split`, issue #64's T1):
    !|   ancestor [1, 1, 0], genes along the three axes, `rdi_threshold` 0: {1, 2} rebuilds the
    !|   ancestor exactly, the word 3; {1, 3} and {2, 3} are sqrt(2) off.
    !| In place of the old `test_detect_patterns_subfunc_at_angle_margin` (residuals 0.0008 and
    !| 0.002 against 1e-3) and `test_detect_patterns_subfunc_floating_point_epsilon` (a residual
    !| of 1e-12 against 1e-9, found only through a made-up norm of 1 for a gene 1 + 1e-12 long).
    subroutine test_detect_subfunc_threshold_inclusive()
        real(real64), parameter :: pair_ancestor(2) = [1.75_real64, 2.0_real64]
        real(real64), parameter :: pair_genes(2, 2) = reshape([1.0_real64, 0.0_real64, 0.0_real64, 1.0_real64], [2, 2])
        real(real64), parameter :: split_ancestor(3) = [1.0_real64, 1.0_real64, 0.0_real64]
        real(real64), parameter :: split_genes(3, 3) = reshape([1.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64, &
                                                                0.0_real64, 0.0_real64, 0.0_real64, 1.0_real64], [3, 3])
        integer(int32), parameter :: pair(1, 1) = reshape([3], [1, 1])
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr

        mask = 3
        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(pair_ancestor, pair_genes, 2, 2, mask, 1, n_results, 2, subsets, 1, 1.25_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_threshold_inclusive: at 1.25, ierr")
        call check_subsets(1, 1, subsets, n_results, 1, pair, "test_detect_subfunc_threshold_inclusive: at 1.25")

        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(pair_ancestor, pair_genes, 2, 2, mask, 1, n_results, 2, subsets, 1, &
                                         nearest(1.25_real64, -1.0_real64), ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_threshold_inclusive: below 1.25, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_threshold_inclusive: below 1.25 finds nothing")

        mask = 7
        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(split_ancestor, split_genes, 3, 3, mask, 1, n_results, 3, subsets, 3, 0.0_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_threshold_inclusive: perfect split, ierr")
        call check_subsets(1, 3, subsets, n_results, 1, pair, "test_detect_subfunc_threshold_inclusive: perfect split")
    end subroutine test_detect_subfunc_threshold_inclusive

    !> The pruning bound itself: a subset grows while the smallest norm after it is at most its
    !| residual's length. Ancestor [4, 0], genes [1.5, 1] and [1.5, -1], whose pair leaves the
    !| residual [1, 0], and a third gene; `rdi_threshold` 0.25:
    !| Through the expert tier, with the norms and their order given.
    !| - third gene [1, 0], norm 1, equal to the residual: {1, 2} grows, and {1, 2, 3} leaves
    !|   nothing, the word 7;
    !| - third gene [1.25, 0], norm 1.25, longer than the residual: {1, 2} is pruned, so nothing is
    !|   found -- although {1, 2, 3} would leave [-0.25, 0], within the threshold. The rule is a
    !|   heuristic, as issue #64 states it.
    subroutine test_detect_subfunc_pruning_bound()
        real(real64), parameter :: ancestor(2) = [4.0_real64, 0.0_real64]
        integer(int32), parameter :: ascending(3) = [3, 1, 2]
        integer(int32), parameter :: triple(1, 1) = reshape([7], [1, 1])
        real(real64) :: genes(2, 3), norms(3)
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr

        mask = 7
        genes = reshape([1.5_real64, 1.0_real64, 1.5_real64, -1.0_real64, 1.0_real64, 0.0_real64], [2, 3])
        norms = [sqrt(3.25_real64), sqrt(3.25_real64), 1.0_real64]
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, 0.25_real64, norms, ascending, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_pruning_bound: norm equal to the residual, ierr")
        call check_subsets(1, 3, subsets, n_results, 1, triple, "test_detect_subfunc_pruning_bound: norm equal to the residual")

        genes(1, 3) = 1.25_real64
        norms(3) = 1.25_real64
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, 0.25_real64, norms, ascending, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_pruning_bound: norm above the residual, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_pruning_bound: norm above the residual prunes")
    end subroutine test_detect_subfunc_pruning_bound

    !> The bound takes the smallest norm of all the genes after a subset's last gene, not of the
    !| next one only (issue #64's example of the minimum array, [50, 75, 0, 100, 25]). Ancestor
    !| [6, 0], genes [3, 0], [2, 0], [8, 0] and [1, 0], norms [3, 2, 8, 1], `rdi_threshold` 0:
    !| {1, 2} leaves 1; the next gene is 8 long, but gene 4 is 1 long, so {1, 2} grows, and
    !| {1, 2, 4} rebuilds the ancestor, the word 11. It is the one result: no pair or other triple
    !| sums to 6, and {1, 2, 3, 4} sums to 14. Through the expert tier, with the order given.
    subroutine test_detect_subfunc_suffix_minimum()
        integer(int32), parameter :: n_genes = 4, n_slots = 7
        real(real64), parameter :: ancestor(2) = [6.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, n_genes) = reshape([3.0_real64, 0.0_real64, 2.0_real64, 0.0_real64, 8.0_real64, &
                                                                0.0_real64, 1.0_real64, 0.0_real64], [2, n_genes])
        real(real64), parameter :: norms(n_genes) = [3.0_real64, 2.0_real64, 8.0_real64, 1.0_real64]
        integer(int32), parameter :: ascending(n_genes) = [4, 2, 1, 3]
        integer(int32), parameter :: expected(1, 1) = reshape([11], [1, 1])
        integer(int32) :: mask(1), subsets(1, n_slots), n_results, ierr

        mask = 15
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(ancestor, genes, n_genes, 2, mask, 1, n_results, n_genes, subsets, n_slots, 0.0_real64, norms, &
                            ascending, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_suffix_minimum: ierr")
        call check_subsets(1, n_slots, subsets, n_results, 1, expected, "test_detect_subfunc_suffix_minimum")
    end subroutine test_detect_subfunc_suffix_minimum

    !> Nothing to find, while any pair of the three genes [1, 0] rebuilds the ancestor [2, 0]:
    !| - no gene kept (the word 0), and only gene 2 kept (the word 2): no pair to form;
    !| - all three kept, `max_subset_size` 1: a single gene is never a result;
    !| - all three kept, `max_subset_size` 0: the early exit, nothing checked;
    !| - a problem of one gene;
    !| - no gene kept and only gene 2 kept, with a work array of no subsets, which is what the
    !|   sizing gives for fewer than two kept genes (and what the Python and R bindings pass on).
    !| `n_results` is 0 each time. The work array holds three subsets, room for the starting genes,
    !| but in the last case.
    subroutine test_detect_subfunc_no_candidates()
        integer(int32), parameter :: n_genes = 3, n_slots = 3, n_cases = 4
        real(real64), parameter :: ancestor(2) = [2.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, n_genes) = reshape([1.0_real64, 0.0_real64, 1.0_real64, 0.0_real64, &
                                                                1.0_real64, 0.0_real64], [2, n_genes])
        integer(int32), parameter :: case_masks(n_cases) = [0, 2, 7, 7], case_sizes(n_cases) = [3, 3, 1, 0]
        character(*), parameter :: case_names(n_cases) = ["no gene kept     ", "one gene kept    ", "max_subset_size 1", &
                                                          "max_subset_size 0"]
        integer(int32) :: mask(1), subsets(1, n_slots), n_results, ierr, i_case

        do i_case = 1, n_cases
            mask = case_masks(i_case)
            n_results = UNWRITTEN_COUNT
            call detect_subfunctionalization(ancestor, genes, n_genes, 2, mask, 1, n_results, case_sizes(i_case), subsets, &
                                             n_slots, 0.0_real64, ierr)
            call assert_err(ierr, ERR_OK, "test_detect_subfunc_no_candidates: "//trim(case_names(i_case))//", ierr")
            call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_no_candidates: "//trim(case_names(i_case)))
        end do

        mask = 1
        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(ancestor, genes(:, 1), 1, 2, mask, 1, n_results, 1, subsets, 1, 0.0_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_no_candidates: one gene, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_no_candidates: one gene")

        do i_case = 1, 2
            mask = case_masks(i_case)
            n_results = UNWRITTEN_COUNT
            call detect_subfunctionalization(ancestor, genes, n_genes, 2, mask, 1, n_results, 3, subsets, 0, 0.0_real64, ierr)
            call assert_err(ierr, ERR_OK, "test_detect_subfunc_no_candidates: "//trim(case_names(i_case))//", no subsets, ierr")
            call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_no_candidates: "//trim(case_names(i_case)) &
                                  //", no subsets")
        end do
    end subroutine test_detect_subfunc_no_candidates

    !> The documented chain, filter -> size -> detect, on the mixed family (see the fixtures), after
    !| the old suite's `test_detect_patterns_mixed_results`: the subfunctionalization filter at pi/8
    !| keeps genes 1 and 2 (the word 3), the sizing gives (2, 1), and with `rdi_threshold` 1e-3 the
    !| pair {1, 2}, which rebuilds the ancestor, is the result 3. The plain tier computes the norms
    !| and their order itself.
    subroutine test_detect_subfunc_after_filter()
        integer(int32), parameter :: expected(1, 1) = reshape([3], [1, 1])
        integer(int32) :: gene_to_fam(MIXED_N_GENES), family_bitmasks(1, 1), max_subset_size, work_array_size, n_results, ierr
        integer(int32), allocatable :: subsets(:, :)

        gene_to_fam = 1
        call filter_paralogs_by_pattern_subfunctionalization(MIXED_ANGLES, PI/8, MIXED_N_GENES, 1, gene_to_fam, &
                                                             family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_after_filter: filter, ierr")
        call assert_equal_int(family_bitmasks(1, 1), 3_int32, "test_detect_subfunc_after_filter: filter keeps 1 and 2")

        max_subset_size = MIXED_N_GENES
        call calc_paralog_subsets_bitmasks_size(max_subset_size, MIXED_N_GENES, work_array_size, family_bitmasks(:, 1), 1, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_after_filter: sizing, ierr")
        call assert_equal_int(max_subset_size, 2_int32, "test_detect_subfunc_after_filter: two genes kept")
        call assert_equal_int(work_array_size, 1_int32, "test_detect_subfunc_after_filter: one subset")

        allocate (subsets(1, work_array_size))
        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(MIXED_ANCESTOR, MIXED_GENES, MIXED_N_GENES, MIXED_N_DIMS, family_bitmasks(:, 1), 1, &
                                         n_results, max_subset_size, subsets, work_array_size, 1.0e-3_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_after_filter: detection, ierr")
        call check_subsets(1, work_array_size, subsets, n_results, 1, expected, "test_detect_subfunc_after_filter")
    end subroutine test_detect_subfunc_after_filter

    !> The errors the detection itself finds, past the validation:
    !| - a kept bit past `n_genes` (gene 4 of 3, the word 9) is ERR_INVALID_INPUT; it names no
    !|   argument today, though the mask is at fault, so only the code is checked;
    !| - a work array too small for the search is ERR_SIZE_MISMATCH, at no argument: the worst
    !|   case of four genes (zero genes and the ancestor, see `test_calc_size_enough_for_worst_case`)
    !|   needs 7 subsets and gets 3, room for its three starting genes only;
    !| - kept from the old suite (`test_detect_patterns_drops_its_helpers_positions`): an entry of
    !|   `paralog_norms_perm` past `n_genes`, which the generated validation refuses first,
    !|   reaches the implementation when it is called directly; the internal routine that finds it
    !|   numbers its own arguments, so the error comes back at no position.
    subroutine test_detect_subfunc_early_errors()
        real(real64), parameter :: worst_ancestor(2) = [3.0_real64, 4.0_real64]
        real(real64), parameter :: worst_genes(2, 4) = reshape([0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, 0.0_real64, &
                                                                0.0_real64, 3.0_real64, 4.0_real64], [2, 4])
        integer(int32), parameter :: past_n_genes(3) = [1, 2, 4]
        real(real64) :: tmp_paralog_vector(2), tmp_work_array(3)
        integer(int32) :: mask(1), subsets(1, 3), tmp_active_bitmask(1), n_results, ierr

        mask = 9
        call detect_subfunctionalization(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_early_errors: a bit past n_genes")

        mask = 15
        call detect_subfunctionalization(worst_ancestor, worst_genes, 4, 2, mask, 1, n_results, 4, subsets, 3, 0.0_real64, ierr)
        call assert_err(ierr, ERR_SIZE_MISMATCH, "test_detect_subfunc_early_errors: a work array too small", 0)

        mask = 7
        call detect_patterns_impl(ancestor=BINARY_ANCESTOR, genes=BINARY_GENES, n_genes=3, n_dims=2, &
                                  pattern_mode=MODE_SUBFUNC_PATTERN, filtered_paralogs_bitmask=mask, n_mask_chunks=1, &
                                  n_results=n_results, max_subset_size=3, paralog_subsets_bitmasks=subsets, &
                                  n_paralog_subsets=3, tmp_active_bitmask=tmp_active_bitmask, &
                                  tmp_paralog_vector=tmp_paralog_vector, rdi_threshold=1.0_real64, &
                                  paralog_norms=BINARY_NORMS, paralog_norms_perm=past_n_genes, &
                                  tmp_work_array=tmp_work_array, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_early_errors: a permutation entry past n_genes", 0)
    end subroutine test_detect_subfunc_early_errors

    !> The extents, each refused with its position:
    !| - `n_genes` (argument 3) and `n_dims` (argument 4): ERR_EMPTY_INPUT at 0, ERR_INVALID_INPUT
    !|   at -1;
    !| - `n_mask_chunks` (argument 6) must hold `n_genes` bits: 0 words for 3 genes, and 1 word
    !|   for 33, are ERR_INVALID_INPUT;
    !| - `max_subset_size` (argument 8) is ERR_INVALID_INPUT at -1 (0 is the early exit in
    !|   `test_detect_subfunc_no_candidates`);
    !| - `n_paralog_subsets` (argument 10): ERR_INVALID_INPUT at -1; 0 passes where nothing is
    !|   to be stored (`test_detect_subfunc_no_candidates`).
    subroutine test_detect_subfunc_dimensions()
        real(real64) :: genes_33(2, 33)
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr

        mask = 3
        genes_33 = 1.0_real64

        call subfunc(BINARY_GENES, 0, 2, 1, 2, 3)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_subfunc_dimensions: n_genes 0", 3)
        call subfunc(BINARY_GENES, -1, 2, 1, 2, 3)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_dimensions: n_genes -1", 3)

        call subfunc(BINARY_GENES, 3, 0, 1, 2, 3)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_subfunc_dimensions: n_dims 0", 4)
        call subfunc(BINARY_GENES, 3, -1, 1, 2, 3)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_dimensions: n_dims -1", 4)

        call subfunc(BINARY_GENES, 3, 2, 0, 2, 3)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_dimensions: no word for 3 genes", 6)
        call subfunc(genes_33, 33, 2, 1, 2, 3)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_dimensions: one word for 33 genes", 6)

        call subfunc(BINARY_GENES, 3, 2, 1, -1, 3)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_dimensions: max_subset_size -1", 8)

        call subfunc(BINARY_GENES, 3, 2, 1, 2, -1)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_dimensions: n_paralog_subsets -1", 10)
    contains
        !> The detection on the binary family's ancestor with the given extents, `rdi_threshold` 1.
        subroutine subfunc(genes, n_genes, n_dims, n_mask_chunks, max_subset_size, n_paralog_subsets)
            real(real64), intent(in) :: genes(2, *)
            integer(int32), intent(in) :: n_genes, n_dims, n_mask_chunks, max_subset_size, n_paralog_subsets

            call detect_subfunctionalization(BINARY_ANCESTOR, genes, n_genes, n_dims, mask, n_mask_chunks, n_results, &
                                             max_subset_size, subsets, n_paralog_subsets, 1.0_real64, ierr)
        end subroutine subfunc
    end subroutine test_detect_subfunc_dimensions

    !> The documented ranges, kept from the old suite's `test_detect_patterns_input_validation` and
    !| completed, each bad value in the last element:
    !| - `rdi_threshold` (argument 11) must be at least 0: 0 passes in
    !|   `test_detect_subfunc_suffix_minimum`, -epsilon is ERR_INVALID_INPUT;
    !| - the expert tier's `paralog_norms` (its argument 14) must be at least 0: -epsilon is
    !|   ERR_INVALID_INPUT;
    !| - the expert tier's `paralog_norms_perm` (its argument 15) runs from 1 to `n_genes`: 0 and
    !|   `n_genes` + 1 are ERR_INVALID_INPUT.
    subroutine test_detect_subfunc_value_ranges()
        real(real64) :: norms(3)
        integer(int32) :: perm(3), mask(1), subsets(1, 3), n_results, ierr

        mask = 7
        norms = BINARY_NORMS(1:3)
        perm = BINARY_PERM(1:3)

        call detect_subfunctionalization(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                         -epsilon(1.0_real64), ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_value_ranges: rdi_threshold below 0", 11)

        norms(3) = -epsilon(1.0_real64)
        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, norms, perm, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_value_ranges: a norm below 0", 14)
        norms = BINARY_NORMS(1:3)

        perm(3) = 0
        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, norms, perm, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_value_ranges: a permutation entry 0", 15)
        perm(3) = 4
        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, norms, perm, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_value_ranges: a permutation entry n_genes + 1", 15)
    end subroutine test_detect_subfunc_value_ranges

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `ancestor` (argument 1) and `genes` (argument 2), each
    !| in the last element, in `rdi_threshold` (argument 11), and in the expert tier's
    !| `paralog_norms` (its argument 14), in the last element. Kept from the old suite's `test_detect_patterns_input_validation` (+Inf
    !| in the ancestor and a gene, NaN in a norm), now with every non-finite value.
    subroutine test_detect_subfunc_rejects_nan_and_inf()
        real(real64) :: ancestor(2), genes(2, 3), norms(3)
        integer(int32) :: perm(3), mask(1), subsets(1, 3), n_results, ierr, i_value
        character(:), allocatable :: value_name

        mask = 7
        perm = BINARY_PERM(1:3)
        do i_value = 1, N_NON_FINITE
            value_name = trim(NON_FINITE_NAMES(i_value))
            call reset()
            ancestor(2) = non_finite_value(i_value)
            call detect()
            call assert_err(ierr, ERR_NAN_INF, "test_detect_subfunc_rejects_nan_and_inf: ancestor "//value_name, 1)

            call reset()
            genes(2, 3) = non_finite_value(i_value)
            call detect()
            call assert_err(ierr, ERR_NAN_INF, "test_detect_subfunc_rejects_nan_and_inf: genes "//value_name, 2)

            call reset()
            call detect_subfunctionalization(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, &
                                             non_finite_value(i_value), ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_subfunc_rejects_nan_and_inf: rdi_threshold "//value_name, 11)

            call reset()
            norms(3) = non_finite_value(i_value)
            call subfunc_expert(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, norms, perm, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_subfunc_rejects_nan_and_inf: paralog_norms "//value_name, 14)
        end do
    contains
        subroutine reset()
            ancestor = BINARY_ANCESTOR
            genes = BINARY_GENES(:, 1:3)
            norms = BINARY_NORMS(1:3)
        end subroutine reset

        subroutine detect()
            call detect_subfunctionalization(ancestor, genes, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, ierr)
        end subroutine detect
    end subroutine test_detect_subfunc_rejects_nan_and_inf

    !> The expert entry takes the three scratch arrays from the caller and allocates nothing: it
    !| finds the same subsets in the binary-weights family, and names its arguments by its own
    !| list, where `rdi_threshold` is 13, `paralog_norms` 14 and `paralog_norms_perm` 15
    !| (kept from the old suite's `test_detect_patterns_input_validation`).
    subroutine test_detect_subfunc_expert()
        real(real64) :: tmp_paralog_vector(2), tmp_work_array(BINARY_N_GENES), norms(BINARY_N_GENES)
        integer(int32) :: perm(BINARY_N_GENES), mask(1), subsets(1, BINARY_N_SLOTS), tmp_active_bitmask(1), n_results, ierr

        mask = 31
        norms = BINARY_NORMS
        perm = BINARY_PERM
        n_results = UNWRITTEN_COUNT
        call expert(1.0_real64)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_expert: ierr")
        call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, BINARY_N_RESULTS, BINARY_RESULTS, "test_detect_subfunc_expert")

        call expert(-1.0_real64)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_expert: rdi_threshold -1", 13)

        norms(BINARY_N_GENES) = -1.0_real64
        call expert(1.0_real64)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_expert: a norm of -1", 14)
        norms = BINARY_NORMS

        perm(BINARY_N_GENES) = 0
        call expert(1.0_real64)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_expert: a permutation entry 0", 15)
    contains
        subroutine expert(rdi_threshold)
            real(real64), intent(in) :: rdi_threshold

            call detect_subfunctionalization_expert(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 3, &
                                                    subsets, BINARY_N_SLOTS, tmp_active_bitmask, tmp_paralog_vector, &
                                                    rdi_threshold, norms, perm, tmp_work_array, ierr)
        end subroutine expert
    end subroutine test_detect_subfunc_expert

    !> The expert tier checks only the range of `paralog_norms_perm`, so no index is out of bounds;
    !| the order is the caller's job (FES, 2026-11-25: "if not it's wanted"). On the first three
    !| binary-weights genes (norms [1, 2, 4], ancestor [10, 0], `rdi_threshold` 1):
    !| - a permutation that is not ascending, [3, 2, 1], is accepted with `ierr` OK: the minima are
    !|   wrong (4 for every gene), which prunes nothing here, and no subset of 1, 2 and 4 sums to
    !|   within 1 of 10, so none is found;
    !| - an order with a repeated entry, [1, 1, 1], is accepted too, which is pinned with
    !|   `max_subset_size` 1: the validation runs, and the search stops before any pair is formed,
    !|   so the minima the order would leave unwritten are never read. With entries out of range,
    !|   [1, 1, 4], the same call is refused, naming the order (the expert's argument 15).
    subroutine test_detect_subfunc_order_is_the_callers()
        integer(int32), parameter :: repeated(3) = [1, 1, 1], descending(3) = [3, 2, 1], past_n_genes(3) = [1, 1, 4]
        integer(int32) :: mask(1), subsets(1, 3), n_results, ierr

        mask = 7
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 3, subsets, 3, 1.0_real64, BINARY_NORMS, &
                            descending, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_order_is_the_callers: a descending order, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_order_is_the_callers: a descending order")

        n_results = UNWRITTEN_COUNT
        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 1, subsets, 3, 1.0_real64, BINARY_NORMS, &
                            repeated, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_order_is_the_callers: a repeated entry, ierr")
        call assert_equal_int(n_results, 0_int32, "test_detect_subfunc_order_is_the_callers: a repeated entry")

        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, 3, 2, mask, 1, n_results, 1, subsets, 3, 1.0_real64, BINARY_NORMS, &
                            past_n_genes, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_subfunc_order_is_the_callers: an entry past n_genes", 15)
    end subroutine test_detect_subfunc_order_is_the_callers

    !> The plain tier computes the norms (`calc_paralog_norms`) and sorts their order ascending,
    !| and finds what the expert finds given the true norms and an ascending order:
    !| - the suffix-minimum family (see `test_detect_subfunc_suffix_minimum`): norms [3, 2, 8, 1],
    !|   order [4, 2, 1, 3], the word 11. The order matters here: in index order, the minima would
    !|   be the norms themselves, gene 3's 8 would prune {1, 2}, and nothing would be found;
    !| - the binary-weights family: norms [1, 2, 4, 8, 16], order [1, 2, 3, 4, 5], the words 9,
    !|   10 and 11.
    subroutine test_detect_subfunc_plain_matches_expert()
        real(real64), parameter :: ancestor(2) = [6.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, 4) = reshape([3.0_real64, 0.0_real64, 2.0_real64, 0.0_real64, 8.0_real64, &
                                                          0.0_real64, 1.0_real64, 0.0_real64], [2, 4])
        real(real64), parameter :: norms(4) = [3.0_real64, 2.0_real64, 8.0_real64, 1.0_real64]
        integer(int32), parameter :: ascending(4) = [4, 2, 1, 3]
        integer(int32), parameter :: expected(1, 1) = reshape([11], [1, 1])
        integer(int32) :: mask(1), subsets(1, BINARY_N_SLOTS), n_results, ierr

        mask = 15
        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(ancestor, genes, 4, 2, mask, 1, n_results, 4, subsets, 7, 0.0_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_plain_matches_expert: suffix minimum, plain, ierr")
        call check_subsets(1, 7, subsets, n_results, 1, expected, "test_detect_subfunc_plain_matches_expert: suffix minimum, plain")
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(ancestor, genes, 4, 2, mask, 1, n_results, 4, subsets, 7, 0.0_real64, norms, ascending, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_plain_matches_expert: suffix minimum, expert, ierr")
        call check_subsets(1, 7, subsets, n_results, 1, expected, &
                           "test_detect_subfunc_plain_matches_expert: suffix minimum, expert")

        mask = 31
        n_results = UNWRITTEN_COUNT
        call detect_subfunctionalization(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 3, subsets, &
                                         BINARY_N_SLOTS, 1.0_real64, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_plain_matches_expert: binary weights, plain, ierr")
        call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, BINARY_N_RESULTS, BINARY_RESULTS, &
                           "test_detect_subfunc_plain_matches_expert: binary weights, plain")
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(BINARY_ANCESTOR, BINARY_GENES, BINARY_N_GENES, 2, mask, 1, n_results, 3, subsets, BINARY_N_SLOTS, &
                            1.0_real64, BINARY_NORMS, BINARY_PERM, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_plain_matches_expert: binary weights, expert, ierr")
        call check_subsets(1, BINARY_N_SLOTS, subsets, n_results, BINARY_N_RESULTS, BINARY_RESULTS, &
                           "test_detect_subfunc_plain_matches_expert: binary weights, expert")
    end subroutine test_detect_subfunc_plain_matches_expert

    !> The smallest norm after each gene, gene by gene, as the pruning reads it from the norms'
    !| ascending order. Ancestor [5, 0], genes [2, 0], [2, 0], [1, 0] and [4, 0], norms [2, 2, 1, 4],
    !| ascending order [3, 1, 2, 4], `rdi_threshold` 0. Walking the order: gene 3 (norm 1) is the
    !| smallest of genes 1 to 3; genes 1 and 2 come before it and change nothing; gene 4 (norm 4)
    !| is the smallest of gene 4 alone. So the smallest norm after gene 1 is 1, after gene 2 is 1,
    !| after gene 3 is 4. Then:
    !| - pairs: {3, 4} sums to 5, residual 0, a result (the word 12); {1, 2} leaves 1 and grows, as
    !|   gene 3's norm 1 follows it; {1, 3} and {2, 3} leave 2 and are pruned, only gene 4 (norm 4)
    !|   following; pairs with gene 4 leave 1 and have no gene after them;
    !| - triples: {1, 2, 3} sums to 5, a result (the word 7); {1, 2, 4} leaves -3, and stops.
    !| The words 12 and 7. The minima per gene must read [1, 1, 1, 4]; a fill that let a later,
    !| larger entry of the order overwrite earlier minima would read [1, 2, 4, 4], see 4 after
    !| gene 2, prune {1, 2} and find 12 only.
    subroutine test_detect_subfunc_suffix_minima()
        integer(int32), parameter :: n_genes = 4, n_slots = 7
        real(real64), parameter :: ancestor(2) = [5.0_real64, 0.0_real64]
        real(real64), parameter :: genes(2, n_genes) = reshape([2.0_real64, 0.0_real64, 2.0_real64, 0.0_real64, 1.0_real64, &
                                                                0.0_real64, 4.0_real64, 0.0_real64], [2, n_genes])
        real(real64), parameter :: norms(n_genes) = [2.0_real64, 2.0_real64, 1.0_real64, 4.0_real64]
        integer(int32), parameter :: ascending(n_genes) = [3, 1, 2, 4]
        integer(int32), parameter :: expected(1, 2) = reshape([12, 7], [1, 2])
        integer(int32) :: mask(1), subsets(1, n_slots), n_results, ierr

        mask = 15
        n_results = UNWRITTEN_COUNT
        call subfunc_expert(ancestor, genes, n_genes, 2, mask, 1, n_results, n_genes, subsets, n_slots, 0.0_real64, norms, &
                            ascending, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_subfunc_suffix_minima: ierr")
        call check_subsets(1, n_slots, subsets, n_results, 2, expected, "test_detect_subfunc_suffix_minima")
    end subroutine test_detect_subfunc_suffix_minima

end module mod_test_detect_subfunctionalization
