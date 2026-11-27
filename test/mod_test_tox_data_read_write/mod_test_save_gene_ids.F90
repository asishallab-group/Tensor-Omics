!> The `save_gene_ids` cases: gene ids, a 1-d array of strings, written in the published layout byte
!| for byte, and a longer file in the way replaced. The checks are the suite fixtures' `check_*` for
!| a 1-d array of strings, shared with every `save_*` of that kind.
module mod_test_save_gene_ids
    use tox_data_read_write, only: save_gene_ids
    use test_suite, only: test_case
    use mod_test_tox_data_read_write_fixtures, only: check_save_replaces_strings
    use mod_test_tox_data_read_write_fixtures, only: check_saves_strings
    implicit none
    private
    public :: get_all_tests_save_gene_ids

contains

    !> Every case of `save_gene_ids`.
    function get_all_tests_save_gene_ids() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        allocate (all_tests(2))
        all_tests(1) = test_case("test_save_gene_ids_layout", test_save_gene_ids_layout)
        all_tests(2) = test_case("test_save_gene_ids_replaces_a_file", test_save_gene_ids_replaces_a_file)
    end function get_all_tests_save_gene_ids

    !> In place of the old `test_read_write_gene_ids`, which round-tripped the shipped gene ids and
    !| compared only the first 128 characters of each 256: gene ids are written in the published
    !| layout, byte for byte.
    subroutine test_save_gene_ids_layout()
        call check_saves_strings(save_gene_ids, "test_save_gene_ids_layout")
    end subroutine test_save_gene_ids_layout

    !> A longer file in the way is replaced, not overwritten in part.
    subroutine test_save_gene_ids_replaces_a_file()
        call check_save_replaces_strings(save_gene_ids, "test_save_gene_ids_replaces_a_file")
    end subroutine test_save_gene_ids_replaces_a_file
end module mod_test_save_gene_ids
