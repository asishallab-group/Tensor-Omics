!> The `tox_gene_centroids` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`. The whole suite lives in
!| `test/mod_test_tox_gene_centroids/`. This module only gathers the children's cases;
!| `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_gene_centroids
    use test_suite, only: test_case
    use mod_test_mean_vector, only: get_all_tests_mean_vector
    use mod_test_group_centroid, only: get_all_tests_group_centroid
    implicit none
    private
    public :: get_all_tests_tox_gene_centroids

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_gene_centroids() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_mean_vector(), &
                     get_all_tests_group_centroid()]
    end function get_all_tests_tox_gene_centroids
end module mod_test_tox_gene_centroids
