!> The `tox_normalization` suite: one suite for the module, as #194 asks, with a child module
!| per published procedure, `mod_test_<procedure>.F90`, plus the fixtures they
!| share. The whole suite lives in `test/mod_test_tox_normalization/`. This module only gathers
!| the children's cases; `run_tests` sees one suite and never knows the children exist.
!|
!| A suite can have a directory of its own because `asserts` and `test_suite` are a library, the
!| test_framework package: fpm lets a test source use only library modules and modules in its own
!| directory or below, so `run_tests` reaches down to this module, but nothing here could reach
!| framework modules sitting in `test/` itself.
module mod_test_tox_normalization
    use test_suite, only: test_case
    use mod_test_calc_fchange, only: get_all_tests_calc_fchange
    use mod_test_calc_tiss_avg, only: get_all_tests_calc_tiss_avg
    use mod_test_log2_transformation, only: get_all_tests_log2_transformation
    use mod_test_normalization_pipeline, only: get_all_tests_normalization_pipeline
    use mod_test_normalize_by_std_dev, only: get_all_tests_normalize_by_std_dev
    use mod_test_normalize_unit_length, only: get_all_tests_normalize_unit_length
    use mod_test_quantile_normalization, only: get_all_tests_quantile_normalization
    use mod_test_root_mean_sq_normalization, only: get_all_tests_root_mean_sq_normalization
    implicit none
    private
    public :: get_all_tests_tox_normalization

contains

    !> Every case of every child, in the order of the module's procedures.
    function get_all_tests_tox_normalization() result(all_tests)
        type(test_case), allocatable :: all_tests(:)

        all_tests = [get_all_tests_calc_fchange(), &
                     get_all_tests_calc_tiss_avg(), &
                     get_all_tests_log2_transformation(), &
                     get_all_tests_normalization_pipeline(), &
                     get_all_tests_normalize_by_std_dev(), &
                     get_all_tests_normalize_unit_length(), &
                     get_all_tests_quantile_normalization(), &
                     get_all_tests_root_mean_sq_normalization()]
    end function get_all_tests_tox_normalization
end module mod_test_tox_normalization
