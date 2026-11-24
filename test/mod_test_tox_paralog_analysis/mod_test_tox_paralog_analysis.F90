!> The `tox_paralog_analysis` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`, and the fixtures they share. The two
!| filters share one child, as they differ only in the comparison; each detection's expert entry
!| is tested in its child. The whole suite lives in `test/mod_test_tox_paralog_analysis/`. This
!| module only gathers the children's cases; `run_tests` sees one suite and never knows the
!| children exist.
module mod_test_tox_paralog_analysis
    use test_suite, only: test_case
    use mod_test_detect_neofunctionalization, only: get_all_tests_detect_neofunctionalization
    use mod_test_detect_dosage_effect, only: get_all_tests_detect_dosage_effect
    use mod_test_detect_subfunctionalization, only: get_all_tests_detect_subfunctionalization
    use mod_test_filter_paralogs_by_pattern, only: get_all_tests_filter_paralogs_by_pattern
    use mod_test_calc_paralog_subsets_bitmasks_size, only: get_all_tests_calc_paralog_subsets_bitmasks_size
    implicit none
    private
    public :: get_all_tests_tox_paralog_analysis

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_paralog_analysis() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_detect_neofunctionalization(), &
                     get_all_tests_detect_dosage_effect(), &
                     get_all_tests_detect_subfunctionalization(), &
                     get_all_tests_filter_paralogs_by_pattern(), &
                     get_all_tests_calc_paralog_subsets_bitmasks_size()]
    end function get_all_tests_tox_paralog_analysis
end module mod_test_tox_paralog_analysis
