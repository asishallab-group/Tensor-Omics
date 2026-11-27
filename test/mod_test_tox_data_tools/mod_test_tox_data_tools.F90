!> The `tox_data_tools` suite: one suite for the module, as #194 asks, with a child module per
!| published procedure, `mod_test_<procedure>.F90`, and the fixtures they share. The whole suite
!| lives in `test/mod_test_tox_data_tools/`. This module only gathers the children's cases;
!| `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_data_tools
    use test_suite, only: test_case
    use mod_test_read_expression_vectors_tsv, only: get_all_tests_read_expression_vectors_tsv
    use mod_test_read_gene_ids_from_tsv_file, only: get_all_tests_read_gene_ids_from_tsv_file
    use mod_test_read_orthofinder_file, only: get_all_tests_read_orthofinder_file
    use mod_test_split_string, only: get_all_tests_split_string
    use mod_test_get_unassigned_mask, only: get_all_tests_get_unassigned_mask
    use mod_test_apply_unassigned_mask, only: get_all_tests_apply_unassigned_mask
    implicit none
    private
    public :: get_all_tests_tox_data_tools

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_data_tools() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_read_expression_vectors_tsv(), &
                     get_all_tests_read_gene_ids_from_tsv_file(), &
                     get_all_tests_read_orthofinder_file(), &
                     get_all_tests_split_string(), &
                     get_all_tests_get_unassigned_mask(), &
                     get_all_tests_apply_unassigned_mask()]
    end function get_all_tests_tox_data_tools
end module mod_test_tox_data_tools
