!> The `detect_neofunctionalization` cases: a gene is flagged on an axis where it differs from
!| its family's ancestor by more than that axis' threshold, `|ancestor - gene| > threshold`, and
!| never when it is in no family (`gene_to_fam` 0).
!|
!| Every value is a multiple of a power of two (0.125, 0.25, ...), so each difference is exact and
!| every flag is decided by an exact comparison; the cases at the threshold itself check that the
!| comparison is strict. The output is (genes, axes), the transpose of `genes`, so the cases use
!| more genes than axes and check the orientation too.
module mod_test_detect_neofunctionalization
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_paralog_analysis, only: detect_neofunctionalization
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_paralog_analysis_fixtures, only: N_NON_FINITE, NON_FINITE_NAMES, non_finite_value
    implicit none
    private
    public :: get_all_tests_detect_neofunctionalization

    logical(c_bool), parameter :: T = .true., F = .false.

contains

    !> Every case of `detect_neofunctionalization`.
    function get_all_tests_detect_neofunctionalization() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(7))
        all_tests(1) = test_case("test_detect_neofunc_threshold_is_strict", test_detect_neofunc_threshold_is_strict)
        all_tests(2) = test_case("test_detect_neofunc_family_ids", test_detect_neofunc_family_ids)
        all_tests(3) = test_case("test_detect_neofunc_range_ends", test_detect_neofunc_range_ends)
        all_tests(4) = test_case("test_detect_neofunc_one_of_each", test_detect_neofunc_one_of_each)
        all_tests(5) = test_case("test_detect_neofunc_dimensions", test_detect_neofunc_dimensions)
        all_tests(6) = test_case("test_detect_neofunc_value_ranges", test_detect_neofunc_value_ranges)
        all_tests(7) = test_case("test_detect_neofunc_rejects_nan_and_inf", test_detect_neofunc_rejects_nan_and_inf)
    end function get_all_tests_detect_neofunctionalization

    !> In place of the old suite's `test_detect_neofunctionalization` (genes equal to their
    !| ancestor, then shifted by a threshold times the family id, whose flags at the threshold
    !| hung on rounding) and the old Python and R value checks: one family, ancestor
    !| [0.5, -0.25, 0], thresholds 0.25 on every axis.
    !| - gene 1, [0.75, -0.25, 0.5]: differences [0.25, 0, 0.5], flags [no, no, yes] -- 0.25 is the
    !|   threshold itself, and the comparison is strict;
    !| - gene 2, [0.125, 0.25, -0.25]: differences [0.375, 0.5, 0.25], flags [yes, yes, no].
    subroutine test_detect_neofunc_threshold_is_strict()
        integer(int32), parameter :: n_axes = 3, n_families = 1, n_genes = 2
        real(real64) :: ancestors(n_axes, n_families), genes(n_axes, n_genes), thresholds(n_axes)
        integer(int32) :: gene_to_fam(n_genes), ierr
        logical(c_bool) :: neofunc(n_genes, n_axes), expected(n_genes, n_axes)

        ancestors(:, 1) = [0.5_real64, -0.25_real64, 0.0_real64]
        genes(:, 1) = [0.75_real64, -0.25_real64, 0.5_real64]
        genes(:, 2) = [0.125_real64, 0.25_real64, -0.25_real64]
        gene_to_fam = 1
        thresholds = 0.25_real64
        expected = reshape([F, T, &
                            F, T, &
                            T, F], [n_genes, n_axes])

        call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_neofunc_threshold_is_strict: ierr")
        call assert_equal_array_logical(neofunc, expected, n_genes*n_axes, "test_detect_neofunc_threshold_is_strict: neofunc", &
                                        n_genes)
    end subroutine test_detect_neofunc_threshold_is_strict

    !> Each gene is compared with its own family's ancestor, and a gene in no family is never
    !| flagged. Ancestors [0.5, 0.5] (family 1) and [-0.5, 0] (family 2), thresholds 0.25:
    !| - gene 1 is family 1's ancestor exactly, but in family 2: differences [1, 0.5], flagged on
    !|   both axes (compared with family 1 it would be on neither);
    !| - gene 2, [0.5, 0], family 1: differences [0, 0.5], flagged on axis 2;
    !| - gene 3, [-1, 1], in no family: never flagged, however far it is from either ancestor;
    !| - gene 4, [-0.5, 0.125], family 2: differences [0, 0.125], not flagged.
    !| The output starts all `.true.`, so a flag the procedure never writes shows.
    subroutine test_detect_neofunc_family_ids()
        integer(int32), parameter :: n_axes = 2, n_families = 2, n_genes = 4
        real(real64) :: ancestors(n_axes, n_families), genes(n_axes, n_genes), thresholds(n_axes)
        integer(int32) :: gene_to_fam(n_genes), ierr
        logical(c_bool) :: neofunc(n_genes, n_axes), expected(n_genes, n_axes)

        ancestors(:, 1) = [0.5_real64, 0.5_real64]
        ancestors(:, 2) = [-0.5_real64, 0.0_real64]
        genes(:, 1) = [0.5_real64, 0.5_real64]
        genes(:, 2) = [0.5_real64, 0.0_real64]
        genes(:, 3) = [-1.0_real64, 1.0_real64]
        genes(:, 4) = [-0.5_real64, 0.125_real64]
        gene_to_fam = [2, 1, 0, 2]
        thresholds = 0.25_real64
        expected = reshape([T, F, F, F, &
                            T, T, F, F], [n_genes, n_axes])

        neofunc = .true.
        call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_neofunc_family_ids: ierr")
        call assert_equal_array_logical(neofunc, expected, n_genes*n_axes, "test_detect_neofunc_family_ids: neofunc", n_genes)
    end subroutine test_detect_neofunc_family_ids

    !> Ancestors, genes and thresholds at both ends of their range [-1, 1]. Ancestor [1, -1],
    !| thresholds [1, -1]:
    !| - gene 1, [-1, -1]: differences [2, 0]; 2 > 1, and 0 > -1, so both axes are flagged -- a
    !|   threshold of -1 flags every gene in a family on that axis, however close;
    !| - gene 2, [0, 1]: differences [1, 2]; 1 is not above the threshold 1, 2 > -1.
    subroutine test_detect_neofunc_range_ends()
        integer(int32), parameter :: n_axes = 2, n_families = 1, n_genes = 2
        real(real64) :: ancestors(n_axes, n_families), genes(n_axes, n_genes), thresholds(n_axes)
        integer(int32) :: gene_to_fam(n_genes), ierr
        logical(c_bool) :: neofunc(n_genes, n_axes), expected(n_genes, n_axes)

        ancestors(:, 1) = [1.0_real64, -1.0_real64]
        genes(:, 1) = [-1.0_real64, -1.0_real64]
        genes(:, 2) = [0.0_real64, 1.0_real64]
        gene_to_fam = 1
        thresholds = [1.0_real64, -1.0_real64]
        expected = reshape([T, F, &
                            T, T], [n_genes, n_axes])

        call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_neofunc_range_ends: ierr")
        call assert_equal_array_logical(neofunc, expected, n_genes*n_axes, "test_detect_neofunc_range_ends: neofunc", n_genes)
    end subroutine test_detect_neofunc_range_ends

    !> The smallest problem, one family, one axis and one gene: [0] against the ancestor [0.5]
    !| with threshold 0.25 is flagged; the same gene in no family is not. And a larger problem
    !| whose genes are all in no family flags nothing.
    subroutine test_detect_neofunc_one_of_each()
        real(real64) :: ancestor(1, 1), gene(1, 1), threshold(1), ancestors(2, 2), genes(2, 3), thresholds(2)
        integer(int32) :: gene_to_fam(1), unassigned(3), ierr
        logical(c_bool) :: neofunc(1, 1), all_neofunc(3, 2)

        ancestor = 0.5_real64
        gene = 0.0_real64
        threshold = 0.25_real64

        gene_to_fam = 1
        neofunc = .false.
        call detect_neofunctionalization(ancestor, 1, gene, 1, gene_to_fam, 1, threshold, neofunc, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_neofunc_one_of_each: one gene, ierr")
        call assert_true(neofunc(1, 1), "test_detect_neofunc_one_of_each: one gene is flagged")

        gene_to_fam = 0
        neofunc = .true.
        call detect_neofunctionalization(ancestor, 1, gene, 1, gene_to_fam, 1, threshold, neofunc, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_neofunc_one_of_each: one unassigned gene, ierr")
        call assert_false(neofunc(1, 1), "test_detect_neofunc_one_of_each: one unassigned gene is not flagged")

        ancestors = 0.5_real64
        genes = -0.5_real64
        thresholds = 0.0_real64
        unassigned = 0
        all_neofunc = .true.
        call detect_neofunctionalization(ancestors, 2, genes, 2, unassigned, 3, thresholds, all_neofunc, ierr)
        call assert_err(ierr, ERR_OK, "test_detect_neofunc_one_of_each: all unassigned, ierr")
        call assert_false(any(all_neofunc), "test_detect_neofunc_one_of_each: all unassigned, nothing flagged")
    end subroutine test_detect_neofunc_one_of_each

    !> `n_families` (argument 2), `n_axes` (argument 4) and `n_genes` (argument 6) are
    !| ERR_EMPTY_INPUT at 0 and ERR_INVALID_INPUT at -1. 1 is in `test_detect_neofunc_one_of_each`.
    subroutine test_detect_neofunc_dimensions()
        real(real64) :: ancestors(1, 1), genes(1, 1), thresholds(1)
        integer(int32) :: gene_to_fam(1), ierr
        logical(c_bool) :: neofunc(1, 1)

        ancestors = 0.5_real64
        genes = 0.5_real64
        thresholds = 0.25_real64
        gene_to_fam = 1

        call detect_neofunctionalization(ancestors, 0, genes, 1, gene_to_fam, 1, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_neofunc_dimensions: n_families 0", 2)
        call detect_neofunctionalization(ancestors, -1, genes, 1, gene_to_fam, 1, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_dimensions: n_families -1", 2)

        call detect_neofunctionalization(ancestors, 1, genes, 0, gene_to_fam, 1, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_neofunc_dimensions: n_axes 0", 4)
        call detect_neofunctionalization(ancestors, 1, genes, -1, gene_to_fam, 1, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_dimensions: n_axes -1", 4)

        call detect_neofunctionalization(ancestors, 1, genes, 1, gene_to_fam, 0, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_EMPTY_INPUT, "test_detect_neofunc_dimensions: n_genes 0", 6)
        call detect_neofunctionalization(ancestors, 1, genes, 1, gene_to_fam, -1, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_dimensions: n_genes -1", 6)
    end subroutine test_detect_neofunc_dimensions

    !> One step past each documented bound is ERR_INVALID_INPUT, naming the argument:
    !| `ancestors` (argument 1), `genes` (argument 3) and `thresholds` (argument 7) just above 1
    !| and just below -1; `gene_to_fam` (argument 5) at -1 and at `n_families` + 1. The bounds
    !| themselves pass in `test_detect_neofunc_range_ends`, and the family id 0 in
    !| `test_detect_neofunc_family_ids`. Each bad value sits in the last element.
    subroutine test_detect_neofunc_value_ranges()
        integer(int32), parameter :: n_axes = 2, n_families = 2, n_genes = 3
        real(real64) :: ancestors(n_axes, n_families), genes(n_axes, n_genes), thresholds(n_axes), beyond(2)
        integer(int32) :: gene_to_fam(n_genes), ierr, i_side
        logical(c_bool) :: neofunc(n_genes, n_axes)

        beyond = [nearest(1.0_real64, 2.0_real64), nearest(-1.0_real64, -2.0_real64)]
        do i_side = 1, 2
            call reset()
            ancestors(n_axes, n_families) = beyond(i_side)
            call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_value_ranges: ancestors "//beyond(i_side), 1)

            call reset()
            genes(n_axes, n_genes) = beyond(i_side)
            call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_value_ranges: genes "//beyond(i_side), 3)

            call reset()
            thresholds(n_axes) = beyond(i_side)
            call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
            call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_value_ranges: thresholds "//beyond(i_side), 7)
        end do

        call reset()
        gene_to_fam(n_genes) = -1
        call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_value_ranges: gene_to_fam -1", 5)

        call reset()
        gene_to_fam(n_genes) = n_families + 1
        call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_detect_neofunc_value_ranges: gene_to_fam n_families + 1", 5)
    contains
        subroutine reset()
            ancestors = 0.5_real64
            genes = 0.25_real64
            thresholds = 0.125_real64
            gene_to_fam = [1, 2, 1]
        end subroutine reset
    end subroutine test_detect_neofunc_value_ranges

    !> NaN, +Inf and -Inf are ERR_NAN_INF in `ancestors` (argument 1), `genes` (argument 3) and
    !| `thresholds` (argument 7), each in the last element, so the check must cover them all.
    !| Kept from the old suite (`test_detect_neofunctionalization_input_validation`, NaN in the
    !| ancestors and +Inf in the genes), now with every non-finite value in every argument.
    subroutine test_detect_neofunc_rejects_nan_and_inf()
        integer(int32), parameter :: n_axes = 3, n_families = 1, n_genes = 2
        real(real64) :: ancestors(n_axes, n_families), genes(n_axes, n_genes), thresholds(n_axes)
        integer(int32) :: gene_to_fam(n_genes), ierr, i_value
        logical(c_bool) :: neofunc(n_genes, n_axes)

        gene_to_fam = 1
        do i_value = 1, N_NON_FINITE
            call reset()
            ancestors(n_axes, n_families) = non_finite_value(i_value)
            call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_neofunc_rejects_nan_and_inf: ancestors " &
                            //trim(NON_FINITE_NAMES(i_value)), 1)

            call reset()
            genes(n_axes, n_genes) = non_finite_value(i_value)
            call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_neofunc_rejects_nan_and_inf: genes " &
                            //trim(NON_FINITE_NAMES(i_value)), 3)

            call reset()
            thresholds(n_axes) = non_finite_value(i_value)
            call detect_neofunctionalization(ancestors, n_families, genes, n_axes, gene_to_fam, n_genes, thresholds, neofunc, ierr)
            call assert_err(ierr, ERR_NAN_INF, "test_detect_neofunc_rejects_nan_and_inf: thresholds " &
                            //trim(NON_FINITE_NAMES(i_value)), 7)
        end do
    contains
        subroutine reset()
            ancestors = 0.5_real64
            genes = 0.25_real64
            thresholds = 0.125_real64
        end subroutine reset
    end subroutine test_detect_neofunc_rejects_nan_and_inf

end module mod_test_detect_neofunctionalization
