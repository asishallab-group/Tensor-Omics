!> The `tox_tissue_versatility` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`. The whole suite lives in
!| `test/mod_test_tox_tissue_versatility/`. This module only gathers the children's cases;
!| `run_tests` sees one suite and never knows the children exist.
module mod_test_tox_tissue_versatility
    use test_suite, only: test_case
    use mod_test_compute_tissue_versatility, only: get_all_tests_compute_tissue_versatility
    implicit none
    private
    public :: get_all_tests_tox_tissue_versatility

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_tissue_versatility() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = get_all_tests_compute_tissue_versatility()
    end function get_all_tests_tox_tissue_versatility
end module mod_test_tox_tissue_versatility
