!> The `apply_unassigned_mask` cases: the genes `mask` keeps, and only those, remain in
!| `gene_ids`, `expression_vectors` (by column) and `gene_to_fam`, in their order, the arrays
!| reallocated to that many; `n_genes_kept` is that count. All genes kept, none kept (refused), a
!| mask that disagrees with `gene_to_fam`, genes with no sample, and a mask longer than the genes.
!|
!| Expected arrays are written out by hand; the values are integers, so comparisons are exact.
module mod_test_apply_unassigned_mask
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_data_tools, only: apply_unassigned_mask
    use tox_errors
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_apply_unassigned_mask

contains

    !> Every case of `apply_unassigned_mask`.
    function get_all_tests_apply_unassigned_mask() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(6))
        all_tests(1) = test_case("test_apply_unassigned_mask_keeps_the_masked", test_apply_unassigned_mask_keeps_the_masked)
        all_tests(2) = test_case("test_apply_unassigned_mask_keeps_all", test_apply_unassigned_mask_keeps_all)
        all_tests(3) = test_case("test_apply_unassigned_mask_refuses_none_kept", &
                                 test_apply_unassigned_mask_refuses_none_kept)
        all_tests(4) = test_case("test_apply_unassigned_mask_follows_the_mask", test_apply_unassigned_mask_follows_the_mask)
        all_tests(5) = test_case("test_apply_unassigned_mask_no_sample", test_apply_unassigned_mask_no_sample)
        all_tests(6) = test_case("test_apply_unassigned_mask_checks_the_mask_size", &
                                 test_apply_unassigned_mask_checks_the_mask_size)
    end function get_all_tests_apply_unassigned_mask

    !> Four genes g1..g4 in families [1, 0, 2, 0], with two samples: gene j's column is
    !| [10*j + 1, 10*j + 2].
    subroutine four_genes(gene_ids, expression_vectors, gene_to_fam)
        character(len=4), allocatable, intent(out) :: gene_ids(:)
        real(real64), allocatable, intent(out) :: expression_vectors(:, :)
        integer(int32), allocatable, intent(out) :: gene_to_fam(:)

        allocate (gene_ids(4), expression_vectors(2, 4), gene_to_fam(4))
        gene_ids = ["g1", "g2", "g3", "g4"]
        expression_vectors = reshape([11.0_real64, 12.0_real64, 21.0_real64, 22.0_real64, &
                                      31.0_real64, 32.0_real64, 41.0_real64, 42.0_real64], [2, 4])
        gene_to_fam = [1, 0, 2, 0]
    end subroutine four_genes

    !> Kept from the old suite's setup (which filtered the shipped data and checked only ierr):
    !| the mask [T, F, T, F] of the families above keeps g1 and g3, their columns and their
    !| families, in that order.
    subroutine test_apply_unassigned_mask_keeps_the_masked()
        character(len=4), allocatable :: gene_ids(:)
        real(real64), allocatable :: expression_vectors(:, :)
        integer(int32), allocatable :: gene_to_fam(:)
        logical(c_bool) :: mask(4)
        character(len=4) :: expected_ids(2)
        real(real64) :: expected_expression(2, 2)
        integer(int32) :: expected_gene_to_fam(2), n_genes_kept, ierr

        call four_genes(gene_ids, expression_vectors, gene_to_fam)
        mask = [.true., .false., .true., .false.]
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)
        call assert_err(ierr, ERR_OK, "test_apply_unassigned_mask_keeps_the_masked: ierr")
        call assert_equal_int(n_genes_kept, 2, "test_apply_unassigned_mask_keeps_the_masked: n_genes_kept")
        call assert_equal_int(size(gene_ids), 2, "test_apply_unassigned_mask_keeps_the_masked: size of gene_ids")
        call assert_equal_int(size(expression_vectors, 1), 2, "test_apply_unassigned_mask_keeps_the_masked: samples")
        call assert_equal_int(size(expression_vectors, 2), 2, "test_apply_unassigned_mask_keeps_the_masked: genes")
        call assert_equal_int(size(gene_to_fam), 2, "test_apply_unassigned_mask_keeps_the_masked: size of gene_to_fam")
        if (size(gene_ids) /= 2 .or. size(expression_vectors) /= 4 .or. size(gene_to_fam) /= 2) return

        expected_ids = ["g1", "g3"]
        call assert_equal_array_char(gene_ids, expected_ids, 4, 2, "test_apply_unassigned_mask_keeps_the_masked: gene_ids")
        expected_expression = reshape([11.0_real64, 12.0_real64, 31.0_real64, 32.0_real64], [2, 2])
        call assert_equal_array_real(expression_vectors, expected_expression, 4, 0.0_real64, &
                                     "test_apply_unassigned_mask_keeps_the_masked: expression_vectors", 2)
        expected_gene_to_fam = [1, 2]
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 2, &
                                    "test_apply_unassigned_mask_keeps_the_masked: gene_to_fam")
    end subroutine test_apply_unassigned_mask_keeps_the_masked

    !> A mask of all genes keeps the arrays as they were.
    subroutine test_apply_unassigned_mask_keeps_all()
        character(len=4), allocatable :: gene_ids(:), expected_ids(:)
        real(real64), allocatable :: expression_vectors(:, :), expected_expression(:, :)
        integer(int32), allocatable :: gene_to_fam(:), expected_gene_to_fam(:)
        logical(c_bool) :: mask(4)
        integer(int32) :: n_genes_kept, ierr

        call four_genes(gene_ids, expression_vectors, gene_to_fam)
        call four_genes(expected_ids, expected_expression, expected_gene_to_fam)
        mask = .true.
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)
        call assert_err(ierr, ERR_OK, "test_apply_unassigned_mask_keeps_all: ierr")
        call assert_equal_int(n_genes_kept, 4, "test_apply_unassigned_mask_keeps_all: n_genes_kept")
        if (size(gene_ids) /= 4 .or. size(expression_vectors) /= 8 .or. size(gene_to_fam) /= 4) then
            call assert_true(.false., "test_apply_unassigned_mask_keeps_all: the sizes changed")
            return
        end if
        call assert_equal_array_char(gene_ids, expected_ids, 4, 4, "test_apply_unassigned_mask_keeps_all: gene_ids")
        call assert_equal_array_real(expression_vectors, expected_expression, 8, 0.0_real64, &
                                     "test_apply_unassigned_mask_keeps_all: expression_vectors", 2)
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 4, "test_apply_unassigned_mask_keeps_all: gene_to_fam")
    end subroutine test_apply_unassigned_mask_keeps_all

    !> A mask of no gene is refused with ERR_INVALID_INPUT, and the arrays are left as they were;
    !| so is a data set of no gene at all.
    subroutine test_apply_unassigned_mask_refuses_none_kept()
        character(len=4), allocatable :: gene_ids(:), expected_ids(:)
        real(real64), allocatable :: expression_vectors(:, :), expected_expression(:, :)
        integer(int32), allocatable :: gene_to_fam(:), expected_gene_to_fam(:)
        logical(c_bool) :: mask(4), no_mask(0)
        integer(int32) :: n_genes_kept, ierr

        call four_genes(gene_ids, expression_vectors, gene_to_fam)
        call four_genes(expected_ids, expected_expression, expected_gene_to_fam)
        mask = .false.
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_apply_unassigned_mask_refuses_none_kept: ierr")
        call assert_equal_int(n_genes_kept, 0, "test_apply_unassigned_mask_refuses_none_kept: n_genes_kept")
        call assert_equal_array_char(gene_ids, expected_ids, 4, 4, "test_apply_unassigned_mask_refuses_none_kept: gene_ids")
        call assert_equal_array_real(expression_vectors, expected_expression, 8, 0.0_real64, &
                                     "test_apply_unassigned_mask_refuses_none_kept: expression_vectors", 2)
        call assert_equal_array_int(gene_to_fam, expected_gene_to_fam, 4, &
                                    "test_apply_unassigned_mask_refuses_none_kept: gene_to_fam")

        deallocate (gene_ids, expression_vectors, gene_to_fam)
        allocate (gene_ids(0), expression_vectors(2, 0), gene_to_fam(0))
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, no_mask, n_genes_kept, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_apply_unassigned_mask_refuses_none_kept: no gene")
    end subroutine test_apply_unassigned_mask_refuses_none_kept

    !> Only the mask decides, not `gene_to_fam`: [F, T, F, F] keeps g2, whose family is 0.
    subroutine test_apply_unassigned_mask_follows_the_mask()
        character(len=4), allocatable :: gene_ids(:)
        real(real64), allocatable :: expression_vectors(:, :)
        integer(int32), allocatable :: gene_to_fam(:)
        logical(c_bool) :: mask(4)
        real(real64) :: expected_expression(2, 1)
        integer(int32) :: n_genes_kept, ierr

        call four_genes(gene_ids, expression_vectors, gene_to_fam)
        mask = [.false., .true., .false., .false.]
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)
        call assert_err(ierr, ERR_OK, "test_apply_unassigned_mask_follows_the_mask: ierr")
        call assert_equal_int(n_genes_kept, 1, "test_apply_unassigned_mask_follows_the_mask: n_genes_kept")
        if (size(gene_ids) /= 1 .or. size(expression_vectors) /= 2 .or. size(gene_to_fam) /= 1) then
            call assert_true(.false., "test_apply_unassigned_mask_follows_the_mask: sizes")
            return
        end if
        call assert_string_equal(gene_ids(1), "g2", "test_apply_unassigned_mask_follows_the_mask: gene_ids")
        expected_expression = reshape([21.0_real64, 22.0_real64], [2, 1])
        call assert_equal_array_real(expression_vectors, expected_expression, 2, 0.0_real64, &
                                     "test_apply_unassigned_mask_follows_the_mask: expression_vectors")
        call assert_equal_int(gene_to_fam(1), 0, "test_apply_unassigned_mask_follows_the_mask: gene_to_fam")
    end subroutine test_apply_unassigned_mask_follows_the_mask

    !> Genes with no sample: the expression matrix keeps no row and its masked columns.
    subroutine test_apply_unassigned_mask_no_sample()
        character(len=4), allocatable :: gene_ids(:)
        real(real64), allocatable :: expression_vectors(:, :)
        integer(int32), allocatable :: gene_to_fam(:)
        logical(c_bool) :: mask(3)
        integer(int32) :: n_genes_kept, ierr

        allocate (gene_ids(3), expression_vectors(0, 3), gene_to_fam(3))
        gene_ids = ["g1", "g2", "g3"]
        gene_to_fam = [1, 0, 1]
        mask = [.true., .false., .true.]
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)
        call assert_err(ierr, ERR_OK, "test_apply_unassigned_mask_no_sample: ierr")
        call assert_equal_int(n_genes_kept, 2, "test_apply_unassigned_mask_no_sample: n_genes_kept")
        call assert_equal_int(size(expression_vectors, 1), 0, "test_apply_unassigned_mask_no_sample: samples")
        call assert_equal_int(size(expression_vectors, 2), 2, "test_apply_unassigned_mask_no_sample: genes")
    end subroutine test_apply_unassigned_mask_no_sample

    !> BUG (fails today): nothing checks that `mask` has an entry per gene. Here it has one more,
    !| and it keeps that one: n_genes_kept is 3 of 2 genes, and the third kept gene is
    !| uninitialized memory, with ERR_OK. (A mask shorter than the genes reads past its end and
    !| cannot be a case here: it aborts a bounds-checked run.) A mask of the wrong size has to be
    !| refused, ERR_SIZE_MISMATCH as the module's other size checks say.
    subroutine test_apply_unassigned_mask_checks_the_mask_size()
        character(len=4), allocatable :: gene_ids(:)
        real(real64), allocatable :: expression_vectors(:, :)
        integer(int32), allocatable :: gene_to_fam(:)
        logical(c_bool) :: mask(3)
        integer(int32) :: n_genes_kept, ierr

        allocate (gene_ids(2), expression_vectors(1, 2), gene_to_fam(2))
        gene_ids = ["g1", "g2"]
        expression_vectors = 1.0_real64
        gene_to_fam = [1, 1]
        mask = .true.
        call apply_unassigned_mask(gene_ids, expression_vectors, gene_to_fam, mask, n_genes_kept, ierr)
        call assert_true(is_err(ierr), "test_apply_unassigned_mask_checks_the_mask_size: a mask of 3 for 2 genes, kept " &
                         //n_genes_kept//" with ERR_OK")
    end subroutine test_apply_unassigned_mask_checks_the_mask_size
end module mod_test_apply_unassigned_mask
