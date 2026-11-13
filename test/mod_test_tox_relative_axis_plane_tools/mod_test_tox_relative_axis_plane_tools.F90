!> The `tox_relative_axis_plane_tools` suite: one suite for the module, as #194 asks, with a
!| child module per published procedure, `mod_test_<procedure>.F90`. The whole suite lives in
!| `test/mod_test_tox_relative_axis_plane_tools/`. This module only gathers the children's
!| cases; `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_relative_axis_plane_tools
    use test_suite, only: test_case
    use mod_test_omics_vector_RAP_projection, only: get_all_tests_omics_vector_RAP_projection
    use mod_test_omics_field_RAP_projection, only: get_all_tests_omics_field_RAP_projection
    use mod_test_clock_hand_angle_between_vectors, only: get_all_tests_clock_hand_angle_between_vectors
    use mod_test_clock_hand_angles_for_shift_vectors, only: get_all_tests_clock_hand_angles_for_shift_vectors
    use mod_test_compute_relative_axis_contributions, only: get_all_tests_compute_relative_axis_contributions
    use mod_test_relative_axes_changes_from_shift_vector, only: get_all_tests_relative_axes_changes_from_shift_vector
    use mod_test_relative_axes_expression_from_expression_vector, &
        only: get_all_tests_relative_axes_expression_from_expression_vector
    implicit none
    private
    public :: get_all_tests_tox_relative_axis_plane_tools

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_relative_axis_plane_tools() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_omics_vector_RAP_projection(), &
                     get_all_tests_omics_field_RAP_projection(), &
                     get_all_tests_clock_hand_angle_between_vectors(), &
                     get_all_tests_clock_hand_angles_for_shift_vectors(), &
                     get_all_tests_compute_relative_axis_contributions(), &
                     get_all_tests_relative_axes_changes_from_shift_vector(), &
                     get_all_tests_relative_axes_expression_from_expression_vector()]
    end function get_all_tests_tox_relative_axis_plane_tools
end module mod_test_tox_relative_axis_plane_tools
