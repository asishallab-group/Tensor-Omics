!> The `tox_data_read_write` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`, and the layout and checks they share. The
!| whole suite lives in `test/mod_test_tox_data_read_write/`. This module only gathers the
!| children's cases; `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_data_read_write
    use test_suite, only: test_case
    use mod_test_save_gene_ids, only: get_all_tests_save_gene_ids
    use mod_test_save_family_ids, only: get_all_tests_save_family_ids
    use mod_test_save_gene_to_family, only: get_all_tests_save_gene_to_family
    use mod_test_save_expression_vectors, only: get_all_tests_save_expression_vectors
    use mod_test_save_family_centroids, only: get_all_tests_save_family_centroids
    use mod_test_save_shift_vectors, only: get_all_tests_save_shift_vectors
    use mod_test_load_gene_ids, only: get_all_tests_load_gene_ids
    use mod_test_load_family_ids, only: get_all_tests_load_family_ids
    use mod_test_load_gene_to_family, only: get_all_tests_load_gene_to_family
    use mod_test_load_expression_vectors, only: get_all_tests_load_expression_vectors
    use mod_test_load_family_centroids, only: get_all_tests_load_family_centroids
    use mod_test_load_shift_vectors, only: get_all_tests_load_shift_vectors
    implicit none
    private
    public :: get_all_tests_tox_data_read_write

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_data_read_write() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_save_gene_ids(), &
                     get_all_tests_save_family_ids(), &
                     get_all_tests_save_gene_to_family(), &
                     get_all_tests_save_expression_vectors(), &
                     get_all_tests_save_family_centroids(), &
                     get_all_tests_save_shift_vectors(), &
                     get_all_tests_load_gene_ids(), &
                     get_all_tests_load_family_ids(), &
                     get_all_tests_load_gene_to_family(), &
                     get_all_tests_load_expression_vectors(), &
                     get_all_tests_load_family_centroids(), &
                     get_all_tests_load_shift_vectors()]
    end function get_all_tests_tox_data_read_write
end module mod_test_tox_data_read_write
