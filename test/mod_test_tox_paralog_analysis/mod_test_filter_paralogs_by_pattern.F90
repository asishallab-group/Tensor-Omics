!> The `filter_paralogs_by_pattern_dosage_effect` and `filter_paralogs_by_pattern_subfunctionalization`
!| cases. Both mark the genes of each family in that family's bit mask, a column of
!| `family_bitmasks`: the dosage-effect filter keeps a gene whose angle is at most the threshold,
!| the subfunctionalization filter one whose angle is at least the threshold. The two share a
!| child because they share everything but that comparison.
!|
!| The angles are multiples of 0.25, so each comparison is exact, and the expected masks are
!| written out as words (bit `i - 1` for gene `i`, see the fixtures).
module mod_test_filter_paralogs_by_pattern
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use tox_paralog_analysis, only: filter_paralogs_by_pattern_dosage_effect, filter_paralogs_by_pattern_subfunctionalization
    use tox_errors
    use test_suite, only: test_case
    use f42_math_impl, only: PI
    use mod_test_tox_paralog_analysis_fixtures, only: BIT_31, UNWRITTEN_WORD, N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    implicit none
    private
    public :: get_all_tests_filter_paralogs_by_pattern

contains

    !> Every case of the two filters.
    function get_all_tests_filter_paralogs_by_pattern() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(8))
        all_tests(1) = test_case("test_filter_paralogs_one_family", test_filter_paralogs_one_family)
        all_tests(2) = test_case("test_filter_paralogs_families", test_filter_paralogs_families)
        all_tests(3) = test_case("test_filter_paralogs_threshold_extremes", test_filter_paralogs_threshold_extremes)
        all_tests(4) = test_case("test_filter_paralogs_words", test_filter_paralogs_words)
        all_tests(5) = test_case("test_filter_paralogs_one_gene", test_filter_paralogs_one_gene)
        all_tests(6) = test_case("test_filter_paralogs_dimensions", test_filter_paralogs_dimensions)
        all_tests(7) = test_case("test_filter_paralogs_value_ranges", test_filter_paralogs_value_ranges)
        all_tests(8) = test_case("test_filter_paralogs_rejects_nan_and_inf", test_filter_paralogs_rejects_nan_and_inf)
    end function get_all_tests_filter_paralogs_by_pattern

    !> In place of the old suite's `test_filter_paralogs_by_pattern` (which counted the kept genes
    !| with the filters' own comparisons) and the old Python and R value checks (angles 0.1 to
    !| 0.9, threshold 0.6, masks 7 and 24): one family, angles [0, 0.25, 0.5, 0.75, 1], threshold
    !| 0.5. The dosage-effect filter keeps genes 1 to 3, the word 7; the subfunctionalization
    !| filter genes 3 to 5, the word 28. Gene 3 sits at the threshold and both keep it. The mask
    !| starts at UNWRITTEN_WORD, so a word the filter never writes shows.
    subroutine test_filter_paralogs_one_family()
        integer(int32), parameter :: n_genes = 5, n_families = 1, n_mask_chunks = 1
        real(real64) :: gene_angles(n_genes)
        integer(int32) :: gene_to_fam(n_genes), family_bitmasks(n_mask_chunks, n_families), ierr

        gene_angles = [0.0_real64, 0.25_real64, 0.5_real64, 0.75_real64, 1.0_real64]
        gene_to_fam = 1

        family_bitmasks = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                      family_bitmasks, n_mask_chunks, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_one_family: dosage effect, ierr")
        call assert_equal_int(family_bitmasks(1, 1), 7_int32, "test_filter_paralogs_one_family: dosage effect keeps genes 1 to 3")

        family_bitmasks = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_subfunctionalization(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                             family_bitmasks, n_mask_chunks, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_one_family: subfunctionalization, ierr")
        call assert_equal_int(family_bitmasks(1, 1), 28_int32, &
                              "test_filter_paralogs_one_family: subfunctionalization keeps genes 3 to 5")
    end subroutine test_filter_paralogs_one_family

    !> Each gene is marked in its own family's column, at its own index. Five genes in three
    !| families, `gene_to_fam` [2, 1, 2, 3, 1], angles [0.25, 1, 0.5, 2, 0.25], threshold 0.5:
    !| - dosage effect keeps genes 1, 3 and 5: family 1 holds gene 5 (the word 16), family 2
    !|   genes 1 and 3 (the word 5), and family 3 nothing (0);
    !| - subfunctionalization keeps genes 2, 3 and 4: family 1 holds gene 2 (the word 2), family
    !|   2 gene 3 (the word 4), and family 3 gene 4 (the word 8).
    subroutine test_filter_paralogs_families()
        integer(int32), parameter :: n_genes = 5, n_families = 3, n_mask_chunks = 1
        real(real64) :: gene_angles(n_genes)
        integer(int32) :: gene_to_fam(n_genes), family_bitmasks(n_mask_chunks, n_families), ierr
        integer(int32) :: expected_dosage(n_mask_chunks, n_families), expected_subfunc(n_mask_chunks, n_families)

        gene_angles = [0.25_real64, 1.0_real64, 0.5_real64, 2.0_real64, 0.25_real64]
        gene_to_fam = [2, 1, 2, 3, 1]
        expected_dosage = reshape([16, 5, 0], [n_mask_chunks, n_families])
        expected_subfunc = reshape([2, 4, 8], [n_mask_chunks, n_families])

        family_bitmasks = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                      family_bitmasks, n_mask_chunks, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_families: dosage effect, ierr")
        call assert_equal_array_int(family_bitmasks, expected_dosage, n_mask_chunks*n_families, &
                                    "test_filter_paralogs_families: dosage effect", n_mask_chunks)

        family_bitmasks = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_subfunctionalization(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                             family_bitmasks, n_mask_chunks, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_families: subfunctionalization, ierr")
        call assert_equal_array_int(family_bitmasks, expected_subfunc, n_mask_chunks*n_families, &
                                    "test_filter_paralogs_families: subfunctionalization", n_mask_chunks)
    end subroutine test_filter_paralogs_families

    !> The angles span their whole range [0, pi], and the threshold has no range but finiteness.
    !| Angles [0, 1, pi]:
    !| - threshold -1, below every angle: dosage effect keeps nothing (0), subfunctionalization
    !|   everything (7);
    !| - threshold 4, above every angle: dosage effect keeps everything, subfunctionalization
    !|   nothing;
    !| - threshold pi: dosage effect keeps everything, subfunctionalization gene 3 (the word 4),
    !|   whose angle is pi itself; threshold 0 the other way round, gene 1 for dosage effect.
    subroutine test_filter_paralogs_threshold_extremes()
        integer(int32), parameter :: n_genes = 3, n_families = 1, n_mask_chunks = 1, n_thresholds = 4
        real(real64) :: gene_angles(n_genes), thresholds(n_thresholds)
        integer(int32) :: gene_to_fam(n_genes), family_bitmasks(n_mask_chunks, n_families), ierr, i_threshold
        integer(int32) :: expected_dosage(n_thresholds), expected_subfunc(n_thresholds)

        gene_angles = [0.0_real64, 1.0_real64, PI]
        gene_to_fam = 1
        thresholds = [-1.0_real64, 4.0_real64, PI, 0.0_real64]
        expected_dosage = [0, 7, 7, 1]
        expected_subfunc = [7, 0, 4, 7]

        do i_threshold = 1, n_thresholds
            family_bitmasks = UNWRITTEN_WORD
            call filter_paralogs_by_pattern_dosage_effect(gene_angles, thresholds(i_threshold), n_genes, n_families, &
                                                          gene_to_fam, family_bitmasks, n_mask_chunks, ierr)
            call assert_err(ierr, ERR_OK, "test_filter_paralogs_threshold_extremes: dosage effect, ierr")
            call assert_equal_int(family_bitmasks(1, 1), expected_dosage(i_threshold), &
                                  "test_filter_paralogs_threshold_extremes: dosage effect at "//thresholds(i_threshold))

            family_bitmasks = UNWRITTEN_WORD
            call filter_paralogs_by_pattern_subfunctionalization(gene_angles, thresholds(i_threshold), n_genes, n_families, &
                                                                 gene_to_fam, family_bitmasks, n_mask_chunks, ierr)
            call assert_err(ierr, ERR_OK, "test_filter_paralogs_threshold_extremes: subfunctionalization, ierr")
            call assert_equal_int(family_bitmasks(1, 1), expected_subfunc(i_threshold), &
                                  "test_filter_paralogs_threshold_extremes: subfunctionalization at "//thresholds(i_threshold))
        end do
    end subroutine test_filter_paralogs_threshold_extremes

    !> Masks of more than one word, in f42's layout. Every gene's angle is 0 or 1, threshold 0.5,
    !| one family:
    !| - 33 genes need two words; genes 1 and 32 kept by dosage effect are the words
    !|   [1 + BIT_31, 0] -- gene 32 is bit 31, which makes the word negative -- and genes 32 and 33
    !|   [BIT_31, 1];
    !| - 64 genes fill two words exactly; gene 64 alone is [0, BIT_31];
    !| - a third word, more than 33 genes need, is in range and comes back 0.
    subroutine test_filter_paralogs_words()
        real(real64) :: angles_33(33), angles_64(64)
        integer(int32) :: to_fam_33(33), to_fam_64(64), two_words(2, 1), three_words(3, 1), ierr
        integer(int32) :: expected_two(2, 1), expected_three(3, 1)

        to_fam_33 = 1
        to_fam_64 = 1

        angles_33 = 1.0_real64
        angles_33([1, 32]) = 0.0_real64
        expected_two = reshape([1_int32 + BIT_31, 0_int32], [2, 1])
        two_words = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(angles_33, 0.5_real64, 33, 1, to_fam_33, two_words, 2, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_words: genes 1 and 32 of 33, ierr")
        call assert_equal_array_int(two_words, expected_two, 2, "test_filter_paralogs_words: genes 1 and 32 of 33")

        angles_33 = 1.0_real64
        angles_33([32, 33]) = 0.0_real64
        expected_two = reshape([BIT_31, 1_int32], [2, 1])
        two_words = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(angles_33, 0.5_real64, 33, 1, to_fam_33, two_words, 2, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_words: genes 32 and 33 of 33, ierr")
        call assert_equal_array_int(two_words, expected_two, 2, "test_filter_paralogs_words: genes 32 and 33 of 33")

        angles_64 = 1.0_real64
        angles_64(64) = 0.0_real64
        expected_two = reshape([0_int32, BIT_31], [2, 1])
        two_words = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(angles_64, 0.5_real64, 64, 1, to_fam_64, two_words, 2, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_words: gene 64 of 64, ierr")
        call assert_equal_array_int(two_words, expected_two, 2, "test_filter_paralogs_words: gene 64 of 64")

        expected_three = reshape([BIT_31, 1_int32, 0_int32], [3, 1])
        three_words = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(angles_33, 0.5_real64, 33, 1, to_fam_33, three_words, 3, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_words: three words for 33 genes, ierr")
        call assert_equal_array_int(three_words, expected_three, 3, "test_filter_paralogs_words: three words for 33 genes")
    end subroutine test_filter_paralogs_words

    !> The smallest problem, one gene in one family: angle 0 is kept by dosage effect (the word 1)
    !| and not by subfunctionalization at threshold 0.5 (the word 0, as when every gene is
    !| filtered out).
    subroutine test_filter_paralogs_one_gene()
        real(real64) :: gene_angles(1)
        integer(int32) :: gene_to_fam(1), family_bitmasks(1, 1), ierr

        gene_angles = 0.0_real64
        gene_to_fam = 1

        family_bitmasks = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_dosage_effect(gene_angles, 0.5_real64, 1, 1, gene_to_fam, family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_one_gene: dosage effect, ierr")
        call assert_equal_int(family_bitmasks(1, 1), 1_int32, "test_filter_paralogs_one_gene: dosage effect keeps it")

        family_bitmasks = UNWRITTEN_WORD
        call filter_paralogs_by_pattern_subfunctionalization(gene_angles, 0.5_real64, 1, 1, gene_to_fam, family_bitmasks, &
                                                             1, ierr)
        call assert_err(ierr, ERR_OK, "test_filter_paralogs_one_gene: subfunctionalization, ierr")
        call assert_equal_int(family_bitmasks(1, 1), 0_int32, "test_filter_paralogs_one_gene: subfunctionalization drops it")
    end subroutine test_filter_paralogs_one_gene

    !> The extents, through both filters:
    !| - `n_genes` (argument 3): ERR_EMPTY_INPUT at 0, ERR_INVALID_INPUT at -1;
    !| - `n_families` (argument 4) runs from 1 to `n_genes`: 0 and `n_genes` + 1 are
    !|   ERR_INVALID_INPUT; `n_genes` itself passes, a family per gene, each column holding its
    !|   own gene;
    !| - `n_mask_chunks` (argument 7) must hold `n_genes` bits: one word is too few for 33 genes,
    !|   and 0 too few for three. Exactly enough passes in `test_filter_paralogs_words`, whose 64
    !|   genes fill two words.
    subroutine test_filter_paralogs_dimensions()
        integer(int32), parameter :: n_genes = 3
        real(real64) :: gene_angles(n_genes), angles_33(33)
        integer(int32) :: gene_to_fam(n_genes), to_fam_33(33), family_bitmasks(1, n_genes), ierr, i_pattern
        integer(int32) :: expected(1, n_genes)

        gene_angles = 0.25_real64
        gene_to_fam = 1
        angles_33 = 0.25_real64
        to_fam_33 = 1

        do i_pattern = 1, 2
            call filter(gene_angles, 0, 1, gene_to_fam, 1)
            call assert_err(ierr, ERR_EMPTY_INPUT, "test_filter_paralogs_dimensions: n_genes 0"//pattern_name(), 3)
            call filter(gene_angles, -1, 1, gene_to_fam, 1)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_dimensions: n_genes -1"//pattern_name(), 3)

            call filter(gene_angles, n_genes, 0, gene_to_fam, 1)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_dimensions: n_families 0"//pattern_name(), 4)
            call filter(gene_angles, n_genes, n_genes + 1, gene_to_fam, 1)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_dimensions: n_families n_genes + 1" &
                            //pattern_name(), 4)

            call filter(angles_33, 33, 1, to_fam_33, 1)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_dimensions: one word for 33 genes"//pattern_name(), 7)
            call filter(gene_angles, n_genes, 1, gene_to_fam, 0)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_dimensions: no word"//pattern_name(), 7)
        end do

        ! a family per gene: gene i alone in column i, by both filters at threshold 0.25
        gene_to_fam = [1, 2, 3]
        expected = reshape([1, 2, 4], [1, n_genes])
        do i_pattern = 1, 2
            family_bitmasks = UNWRITTEN_WORD
            call filter(gene_angles, n_genes, n_genes, gene_to_fam, 1)
            call assert_err(ierr, ERR_OK, "test_filter_paralogs_dimensions: n_families n_genes, ierr"//pattern_name())
            call assert_equal_array_int(family_bitmasks, expected, n_genes, &
                                        "test_filter_paralogs_dimensions: n_families n_genes"//pattern_name())
        end do
    contains
        !> Calls the filter of pattern `i_pattern` at threshold 0.25.
        subroutine filter(angles, n_genes_passed, n_families_passed, to_fam, n_mask_chunks_passed)
            real(real64), contiguous, intent(in) :: angles(:)
            integer(int32), contiguous, intent(in) :: to_fam(:)
            integer(int32), intent(in) :: n_genes_passed, n_families_passed, n_mask_chunks_passed

            if (i_pattern == 1) then
                call filter_paralogs_by_pattern_dosage_effect(angles, 0.25_real64, n_genes_passed, n_families_passed, to_fam, &
                                                              family_bitmasks, n_mask_chunks_passed, ierr)
            else
                call filter_paralogs_by_pattern_subfunctionalization(angles, 0.25_real64, n_genes_passed, n_families_passed, &
                                                                     to_fam, family_bitmasks, n_mask_chunks_passed, ierr)
            end if
        end subroutine filter

        function pattern_name() result(name)
            character(:), allocatable :: name

            if (i_pattern == 1) then
                name = " (dosage effect)"
            else
                name = " (subfunctionalization)"
            end if
        end function pattern_name
    end subroutine test_filter_paralogs_dimensions

    !> One step past each documented bound is ERR_INVALID_INPUT, naming the argument: an angle
    !| (`gene_angles`, argument 1) just above pi or below 0 (-epsilon, a normal number, so no
    !| build reads it as 0); a family id (`gene_to_fam`, argument 5) of 0 -- these filters take no
    !| unassigned genes -- or `n_families` + 1. Each bad value sits in the last element.
    subroutine test_filter_paralogs_value_ranges()
        integer(int32), parameter :: n_genes = 3, n_families = 2
        real(real64) :: gene_angles(n_genes)
        integer(int32) :: gene_to_fam(n_genes), family_bitmasks(1, n_families), ierr

        call reset()
        gene_angles(n_genes) = nearest(PI, 4.0_real64)
        call filter_paralogs_by_pattern_dosage_effect(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                      family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_value_ranges: angle above pi", 1)

        call reset()
        gene_angles(n_genes) = -epsilon(1.0_real64)
        call filter_paralogs_by_pattern_subfunctionalization(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                             family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_value_ranges: angle below 0", 1)

        call reset()
        gene_to_fam(n_genes) = 0
        call filter_paralogs_by_pattern_dosage_effect(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                      family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_value_ranges: family id 0", 5)

        call reset()
        gene_to_fam(n_genes) = n_families + 1
        call filter_paralogs_by_pattern_subfunctionalization(gene_angles, 0.5_real64, n_genes, n_families, gene_to_fam, &
                                                             family_bitmasks, 1, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_filter_paralogs_value_ranges: family id n_families + 1", 5)
    contains
        subroutine reset()
            gene_angles = 0.25_real64
            gene_to_fam = [1, 2, 1]
        end subroutine reset
    end subroutine test_filter_paralogs_value_ranges

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `gene_angles` (argument 1, in the last element) and
    !| in `threshold` (argument 2), through both filters.
    subroutine test_filter_paralogs_rejects_nan_and_inf()
        integer(int32), parameter :: n_genes = 3
        real(real64) :: gene_angles(n_genes)
        integer(int32) :: gene_to_fam(n_genes), family_bitmasks(1, 1), ierr, i_value

        gene_to_fam = 1
        do i_value = 1, N_NON_FINITE
            gene_angles = 0.25_real64
            gene_angles(n_genes) = non_finite_value(i_value)
            call filter_paralogs_by_pattern_dosage_effect(gene_angles, 0.5_real64, n_genes, 1, gene_to_fam, &
                                                          family_bitmasks, 1, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_filter_paralogs_rejects_nan_and_inf: dosage effect, gene_angles " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)
            call filter_paralogs_by_pattern_subfunctionalization(gene_angles, 0.5_real64, n_genes, 1, gene_to_fam, &
                                                                 family_bitmasks, 1, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_filter_paralogs_rejects_nan_and_inf: subfunctionalization, gene_angles " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)

            gene_angles = 0.25_real64
            call filter_paralogs_by_pattern_dosage_effect(gene_angles, non_finite_value(i_value), n_genes, 1, gene_to_fam, &
                                                          family_bitmasks, 1, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_filter_paralogs_rejects_nan_and_inf: dosage effect, threshold " &
                            //trim(NON_FINITE_NAMES(i_value)), 2)
            call filter_paralogs_by_pattern_subfunctionalization(gene_angles, non_finite_value(i_value), n_genes, 1, &
                                                                 gene_to_fam, family_bitmasks, 1, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_filter_paralogs_rejects_nan_and_inf: subfunctionalization, threshold " &
                            //trim(NON_FINITE_NAMES(i_value)), 2)
        end do
    end subroutine test_filter_paralogs_rejects_nan_and_inf

end module mod_test_filter_paralogs_by_pattern
