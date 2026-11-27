!> The `get_family_for_gene_index` cases: a gene's family index is its entry of `gene_to_fam`, 0
!| for a gene in no family; a gene index outside 1..size(gene_to_fam) is refused with
!| ERR_INVALID_INPUT and gives family 0.
module mod_test_get_family_for_gene_index
    use asserts
    use, intrinsic :: iso_fortran_env, only: int32
    use tox_data_accessors, only: get_family_for_gene_index
    use tox_errors
    use test_suite, only: test_case
    use mod_test_tox_data_accessors_fixtures, only: N_GENES, FIXTURE_GENE_TO_FAM
    implicit none
    private
    public :: get_all_tests_get_family_for_gene_index

contains

    !> Every case of `get_family_for_gene_index`.
    function get_all_tests_get_family_for_gene_index() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(2))
        all_tests(1) = test_case("test_get_family_for_gene_index_values", test_get_family_for_gene_index_values)
        all_tests(2) = test_case("test_get_family_for_gene_index_bound", test_get_family_for_gene_index_bound)
    end function get_all_tests_get_family_for_gene_index

    !> Kept from the old `test_data_accessors` (gene 2 of the shipped data is in family 48): each
    !| gene of the fixture gives its entry, [2, 0, 1], the 0 of the gene in no family included.
    subroutine test_get_family_for_gene_index_values()
        integer(int32) :: gene_to_fam(N_GENES), expected(N_GENES), i_gene, family_idx, ierr

        gene_to_fam = FIXTURE_GENE_TO_FAM
        expected = [2, 0, 1]
        do i_gene = 1, N_GENES
            family_idx = -1
            call get_family_for_gene_index(i_gene, gene_to_fam, family_idx, ierr)
            call assert_err(ierr, ERR_OK, "test_get_family_for_gene_index_values: ierr of gene "//i_gene)
            call assert_equal_int(family_idx, expected(i_gene), "test_get_family_for_gene_index_values: gene "//i_gene)
        end do
    end subroutine test_get_family_for_gene_index_values

    !> Both sides of both ends: 1 and 3 are genes, 0 and 4 are refused with ERR_INVALID_INPUT
    !| and give 0; with no gene at all, 1 is refused too.
    subroutine test_get_family_for_gene_index_bound()
        integer(int32) :: gene_to_fam(N_GENES), no_gene(0), family_idx, ierr

        gene_to_fam = FIXTURE_GENE_TO_FAM
        call get_family_for_gene_index(1, gene_to_fam, family_idx, ierr)
        call assert_err(ierr, ERR_OK, "test_get_family_for_gene_index_bound: 1")
        call get_family_for_gene_index(N_GENES, gene_to_fam, family_idx, ierr)
        call assert_err(ierr, ERR_OK, "test_get_family_for_gene_index_bound: 3")

        family_idx = -1
        call get_family_for_gene_index(0, gene_to_fam, family_idx, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_family_for_gene_index_bound: 0")
        call assert_equal_int(family_idx, 0, "test_get_family_for_gene_index_bound: 0 gives family 0")
        family_idx = -1
        call get_family_for_gene_index(N_GENES + 1, gene_to_fam, family_idx, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_family_for_gene_index_bound: 4")
        call assert_equal_int(family_idx, 0, "test_get_family_for_gene_index_bound: 4 gives family 0")
        call get_family_for_gene_index(1, no_gene, family_idx, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_get_family_for_gene_index_bound: 1 of no gene")
    end subroutine test_get_family_for_gene_index_bound
end module mod_test_get_family_for_gene_index
