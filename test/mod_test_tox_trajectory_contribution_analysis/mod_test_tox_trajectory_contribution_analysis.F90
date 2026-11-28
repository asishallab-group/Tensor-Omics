!> The `tox_trajectory_contribution_analysis` suite: one suite for the module, as #194 asks, with
!| a child module per published procedure, `mod_test_<procedure>.F90`, and the fixtures they
!| share. The whole suite lives in `test/mod_test_tox_trajectory_contribution_analysis/`. This
!| module only gathers the children's cases; `run_tests` sees one suite and never knows the
!| children exist.
module mod_test_tox_trajectory_contribution_analysis
    use test_suite, only: test_case
    use mod_test_perform_permutation_test, only: get_all_tests_perform_permutation_test
    use mod_test_compute_p_values, only: get_all_tests_compute_p_values
    use mod_test_compute_contributions, only: get_all_tests_compute_contributions
    use mod_test_compute_all_contributions, only: get_all_tests_compute_all_contributions
    use mod_test_compute_baselines_factor_dependent, only: get_all_tests_compute_baselines_factor_dependent
    use mod_test_compute_velocity_trajectory, only: get_all_tests_compute_velocity_trajectory
    use mod_test_compute_acceleration_from_velocity_trajectory, only: &
        get_all_tests_compute_acceleration_from_velocity_trajectory
    use mod_test_compute_velocity_trajectories, only: get_all_tests_compute_velocity_trajectories
    use mod_test_compute_acceleration_from_velocity, only: get_all_tests_compute_acceleration_from_velocity
    use mod_test_compute_velocity_acceleration_contributions, only: &
        get_all_tests_compute_velocity_acceleration_contributions
    implicit none
    private
    public :: get_all_tests_tox_trajectory_contribution_analysis

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_trajectory_contribution_analysis() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_perform_permutation_test(), &
                     get_all_tests_compute_p_values(), &
                     get_all_tests_compute_contributions(), &
                     get_all_tests_compute_all_contributions(), &
                     get_all_tests_compute_baselines_factor_dependent(), &
                     get_all_tests_compute_velocity_trajectory(), &
                     get_all_tests_compute_acceleration_from_velocity_trajectory(), &
                     get_all_tests_compute_velocity_trajectories(), &
                     get_all_tests_compute_acceleration_from_velocity(), &
                     get_all_tests_compute_velocity_acceleration_contributions()]
    end function get_all_tests_tox_trajectory_contribution_analysis
end module mod_test_tox_trajectory_contribution_analysis
