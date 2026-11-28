!> The `tox_shift_vectors` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`. The whole suite lives in
!| `test/mod_test_tox_shift_vectors/`. This module only gathers the children's cases;
!| `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_shift_vectors
    use test_suite, only: test_case
    use mod_test_compute_shift_vector_field, only: get_all_tests_compute_shift_vector_field
    implicit none
    private
    public :: get_all_tests_tox_shift_vectors

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_shift_vectors() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = get_all_tests_compute_shift_vector_field()
    end function get_all_tests_tox_shift_vectors
end module mod_test_tox_shift_vectors
