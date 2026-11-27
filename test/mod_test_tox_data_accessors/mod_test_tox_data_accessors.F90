!> The `tox_data_accessors` suite: one suite for the module, as #194 asks, with a child module per
!| published procedure, `mod_test_<procedure>.F90`, and the data set they share. The whole suite
!| lives in `test/mod_test_tox_data_accessors/`. This module only gathers the children's cases;
!| `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_data_accessors
    use test_suite, only: test_case
    use mod_test_get_family_for_gene_index, only: get_all_tests_get_family_for_gene_index
    use mod_test_get_expression_vector, only: get_all_tests_get_expression_vector
    use mod_test_get_family_centroid, only: get_all_tests_get_family_centroid
    use mod_test_get_shift_components, only: get_all_tests_get_shift_components
    use mod_test_get_gene_index, only: get_all_tests_get_gene_index
    use mod_test_get_family_index, only: get_all_tests_get_family_index
    implicit none
    private
    public :: get_all_tests_tox_data_accessors

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_data_accessors() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_get_family_for_gene_index(), &
                     get_all_tests_get_expression_vector(), &
                     get_all_tests_get_family_centroid(), &
                     get_all_tests_get_shift_components(), &
                     get_all_tests_get_gene_index(), &
                     get_all_tests_get_family_index()]
    end function get_all_tests_tox_data_accessors
end module mod_test_tox_data_accessors
