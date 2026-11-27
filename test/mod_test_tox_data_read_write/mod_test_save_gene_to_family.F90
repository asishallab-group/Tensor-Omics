!> The `save_gene_to_family` cases: the gene-to-family mapping, a 1-d int32 array, written in the
!| published layout byte for byte, and a longer file in the way replaced. The checks are the suite
!| fixtures' `check_*` for a 1-d int32 array, shared with every `save_*` of that kind.
module mod_test_save_gene_to_family
    use tox_data_read_write, only: save_gene_to_family
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_save_replaces_integers
    use mod_test_tox_data_read_write_fixtures, only: check_saves_integers
    implicit none
    private
    public :: get_all_tests_save_gene_to_family

contains

    !> Every case of `save_gene_to_family`.
    function get_all_tests_save_gene_to_family() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(2))
        all_tests(1) = test_case("test_save_gene_to_family_layout", test_save_gene_to_family_layout)
        all_tests(2) = test_case("test_save_gene_to_family_replaces_a_file", &
                                 test_save_gene_to_family_replaces_a_file)
    end function get_all_tests_save_gene_to_family

    !> In place of the old `test_read_write_gene_to_fam`, which round-tripped the shipped mapping:
    !| the gene-to-family mapping are written in the published layout, byte for byte.
    subroutine test_save_gene_to_family_layout()
        call check_saves_integers(save_gene_to_family, "test_save_gene_to_family_layout")
    end subroutine test_save_gene_to_family_layout

    !> A longer file in the way is replaced, not overwritten in part.
    subroutine test_save_gene_to_family_replaces_a_file()
        call check_save_replaces_integers(save_gene_to_family, "test_save_gene_to_family_replaces_a_file")
    end subroutine test_save_gene_to_family_replaces_a_file
end module mod_test_save_gene_to_family
