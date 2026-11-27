!> The `get_unassigned_mask` cases: a gene is kept when its family index is at least 1, so the
!| sentinel 0 (M_GENE_TO_FAM_SENTINEL) and anything below it are dropped; `n_genes_kept` counts
!| the kept genes. Both sides of the bound, all and none kept, one gene and none.
module mod_test_get_unassigned_mask
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use tox_data_tools, only: get_unassigned_mask
    use test_suite, only: test_case
    implicit none
    private
    public :: get_all_tests_get_unassigned_mask

contains

    !> Every case of `get_unassigned_mask`.
    function get_all_tests_get_unassigned_mask() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(3))
        all_tests(1) = test_case("test_get_unassigned_mask_bound", test_get_unassigned_mask_bound)
        all_tests(2) = test_case("test_get_unassigned_mask_all_or_none", test_get_unassigned_mask_all_or_none)
        all_tests(3) = test_case("test_get_unassigned_mask_one_or_no_gene", test_get_unassigned_mask_one_or_no_gene)
    end function get_all_tests_get_unassigned_mask

    !> Kept from the R suite's `test_3_filter_unassigned_genes` (the mask is gene_to_fam > 0):
    !| the indices 1 and huge are kept, 0, -1 and -huge are not, wherever they stand.
    subroutine test_get_unassigned_mask_bound()
        integer(int32) :: gene_to_fam(6), n_genes_kept
        logical(c_bool) :: mask(6), expected(6)

        gene_to_fam = [0, 1, -1, huge(1_int32), 3, -huge(1_int32)]
        call get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
        expected = [.false., .true., .false., .true., .true., .false.]
        call assert_equal_array_logical(mask, expected, 6, "test_get_unassigned_mask_bound: mask")
        call assert_equal_int(n_genes_kept, 3, "test_get_unassigned_mask_bound: n_genes_kept")
    end subroutine test_get_unassigned_mask_bound

    !> Every gene in a family keeps all of them; none in a family keeps none.
    subroutine test_get_unassigned_mask_all_or_none()
        integer(int32) :: gene_to_fam(3), n_genes_kept
        logical(c_bool) :: mask(3), expected(3)

        gene_to_fam = [2, 1, 2]
        call get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
        expected = .true.
        call assert_equal_array_logical(mask, expected, 3, "test_get_unassigned_mask_all_or_none: all, mask")
        call assert_equal_int(n_genes_kept, 3, "test_get_unassigned_mask_all_or_none: all, n_genes_kept")

        gene_to_fam = 0
        call get_unassigned_mask(gene_to_fam, mask, n_genes_kept)
        expected = .false.
        call assert_equal_array_logical(mask, expected, 3, "test_get_unassigned_mask_all_or_none: none, mask")
        call assert_equal_int(n_genes_kept, 0, "test_get_unassigned_mask_all_or_none: none, n_genes_kept")
    end subroutine test_get_unassigned_mask_all_or_none

    !> One gene, in a family and not; and no gene at all, of which none is kept.
    subroutine test_get_unassigned_mask_one_or_no_gene()
        integer(int32) :: one(1), none(0), n_genes_kept
        logical(c_bool) :: one_mask(1), no_mask(0)

        one = 1
        call get_unassigned_mask(one, one_mask, n_genes_kept)
        call assert_true(one_mask(1), "test_get_unassigned_mask_one_or_no_gene: one in a family, mask")
        call assert_equal_int(n_genes_kept, 1, "test_get_unassigned_mask_one_or_no_gene: one in a family, n_genes_kept")

        one = 0
        call get_unassigned_mask(one, one_mask, n_genes_kept)
        call assert_false(one_mask(1), "test_get_unassigned_mask_one_or_no_gene: one in none, mask")
        call assert_equal_int(n_genes_kept, 0, "test_get_unassigned_mask_one_or_no_gene: one in none, n_genes_kept")

        n_genes_kept = -1
        call get_unassigned_mask(none, no_mask, n_genes_kept)
        call assert_equal_int(n_genes_kept, 0, "test_get_unassigned_mask_one_or_no_gene: no gene")
    end subroutine test_get_unassigned_mask_one_or_no_gene
end module mod_test_get_unassigned_mask
