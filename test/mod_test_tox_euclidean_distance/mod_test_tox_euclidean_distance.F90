!> The `tox_euclidean_distance` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`, and the fixtures they share. The whole
!| suite lives in `test/mod_test_tox_euclidean_distance/`. This module only gathers the children's
!| cases; `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_euclidean_distance
    use test_suite, only: test_case
    use mod_test_euclidean_distance, only: get_all_tests_euclidean_distance
    use mod_test_distance_to_centroid, only: get_all_tests_distance_to_centroid
    implicit none
    private
    public :: get_all_tests_tox_euclidean_distance

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_euclidean_distance() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_euclidean_distance(), &
                     get_all_tests_distance_to_centroid()]
    end function get_all_tests_tox_euclidean_distance
end module mod_test_tox_euclidean_distance
