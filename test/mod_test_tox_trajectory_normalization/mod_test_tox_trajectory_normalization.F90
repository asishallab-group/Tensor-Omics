!> The `tox_trajectory_normalization` suite: one suite for the module, as #194 asks, with a
!| child module per published procedure, `mod_test_<procedure>.F90`. The whole suite lives in
!| `test/mod_test_tox_trajectory_normalization/`. This module only gathers the children's
!| cases; `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_trajectory_normalization
    use test_suite, only: test_case
    use mod_test_normalize_variable_timeseries, only: get_all_tests_normalize_variable_timeseries
    use mod_test_normalize_single_trajectory, only: get_all_tests_normalize_single_trajectory
    use mod_test_normalize_all_trajectories, only: get_all_tests_normalize_all_trajectories
    implicit none
    private
    public :: get_all_tests_tox_trajectory_normalization

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_trajectory_normalization() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_normalize_variable_timeseries(), &
                     get_all_tests_normalize_single_trajectory(), &
                     get_all_tests_normalize_all_trajectories()]
    end function get_all_tests_tox_trajectory_normalization
end module mod_test_tox_trajectory_normalization
