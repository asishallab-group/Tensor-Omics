!> The `tox_clustering` suite: one suite for the module, as #194 asks, with a child module per
!| published procedure, `mod_test_<procedure>.F90`, and the fixtures they share. The expert entry
!| of `cluster_factor_trajectories_k_means` is tested in its child. The internal helpers (the
!| k-means assignment and update, the linkage updates and minimum search) are tested through the
!| published procedures. The whole suite lives in `test/mod_test_tox_clustering/`. This module
!| only gathers the children's cases; `run_tests` sees one suite and never knows the children
!| exist.
module mod_test_tox_clustering
    use test_suite, only: test_case
    use mod_test_cluster_factor_trajectories_k_means, only: get_all_tests_cluster_factor_trajectories_k_means
    use mod_test_k_means_clustering, only: get_all_tests_k_means_clustering
    use mod_test_linkage_clustering, only: get_all_tests_linkage_clustering
    implicit none
    private
    public :: get_all_tests_tox_clustering

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_clustering() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_cluster_factor_trajectories_k_means(), &
                     get_all_tests_k_means_clustering(), &
                     get_all_tests_linkage_clustering()]
    end function get_all_tests_tox_clustering
end module mod_test_tox_clustering
