#include "macros.h"

!> Unit test suite for tox_data_integration routine.
module mod_test_data_integration
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_negative_inf, ieee_positive_inf
    use tox_data_integration_jsd
    use tox_data_integration_preprocessing
    use tox_data_integration_js_comp_test
    use tox_errors
    use f42_utils, only: above, below, init_random, shuffle_vector
    implicit none

    ! Abstract interface for all test procedures
    abstract interface
        subroutine test_interface()
        end subroutine test_interface
    end interface

    ! Type to hold test name and procedure pointer
    type :: test_case
        character(len=128) :: name
        procedure(test_interface), pointer, nopass :: test_proc => null()
    end type test_case

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests() result(all_tests)
        type(test_case) :: all_tests(26)
        ! jsd
        all_tests(22) = test_case("test_build_residual_histograms", test_build_residual_histograms)
        all_tests(16) = test_case("test_compute_divergence_per_reference_point", test_compute_divergence_per_reference_point)
        all_tests(17) = test_case("test_compute_weighted_global_divergence", test_compute_weighted_global_divergence)

        ! js_comp_test
        all_tests(15) = test_case("test_determine_shared_residual_range", test_determine_shared_residual_range)
        all_tests(18) = test_case("test_estimate_bin_count", test_estimate_bin_count)
        all_tests(19) = test_case("test_compute_fractional_overlap", test_compute_fractional_overlap)
        all_tests(20) = test_case("test_test_neighborhood_overlaps", test_test_neighborhood_overlaps)
        all_tests(21) = test_case("test_create_mean_pmf_helpers", test_create_mean_pmf_helpers)
        all_tests(23) = test_case("test_test_mean_pmf_min_counts_helper", test_test_mean_pmf_min_counts_helper)
        all_tests(24) = test_case("test_bootstrap_histogram_helper", test_bootstrap_histogram_helper)
        all_tests(25) = test_case("test_gjct_permutation_test_helper", test_gjct_permutation_test_helper)
        all_tests(26) = test_case("test_check_plateau_condition_helper", test_check_plateau_condition_helper)

        ! preprocessing
        all_tests(1) = test_case("test_compute_gene_means_basic", test_compute_gene_means_basic)
        all_tests(2) = test_case("test_compute_gene_means_with_nan", test_compute_gene_means_with_nan)
        all_tests(3) = test_case("test_compute_gene_means_all_nan", test_compute_gene_means_all_nan)
        all_tests(4) = test_case("test_compute_gene_means_invalid_input", test_compute_gene_means_invalid_input)
        
        all_tests(5) = test_case("test_compute_residuals_basic", test_compute_residuals_basic)
        all_tests(6) = test_case("test_compute_residuals_with_nan", test_compute_residuals_with_nan)
        all_tests(7) = test_case("test_compute_residuals_all_nan", test_compute_residuals_all_nan)
        all_tests(8) = test_case("test_compute_residuals_invalid_input", test_compute_residuals_invalid_input)
        
        all_tests(9) = test_case("test_pool_means_alloc_basic", test_pool_means_alloc_basic)
        all_tests(10) = test_case("test_pool_means_alloc_with_nan", test_pool_means_alloc_with_nan)
        all_tests(11) = test_case("test_pool_means_alloc_single_study", test_pool_means_alloc_single_study)
        all_tests(12) = test_case("test_pool_means_alloc_invalid_input", test_pool_means_alloc_invalid_input)
        
        all_tests(13) = test_case("test_construct_neighborhoods_basic", test_construct_neighborhoods_basic)
        all_tests(14) = test_case("test_construct_neighborhoods_nan_handling", test_construct_neighborhoods_nan_handling)

        ! ! per-family
        ! all_tests(21) = test_case("test_fjct", test_fjct)
        ! all_tests(22) = test_case("test_fjct_compute_contribution_scores", test_fjct_compute_contribution_scores)
    end function get_all_tests

    !> Run all tox_data_integration tests.
    subroutine run_all_tests_tox_data_integration
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i

        all_tests = get_all_tests()

        do i = 1, size(all_tests)
            call all_tests(i)%test_proc()
            print "(' ',A,' passed.')", trim(all_tests(i)%name)
        end do
        print *, "All tox_data_integration tests passed successfully."
    end subroutine run_all_tests_tox_data_integration

    !> Run specific tox_data_integration tests by name.
    subroutine run_named_tests_tox_data_integration(test_names)
        character(len=*), intent(in) :: test_names(:)
        type(test_case), allocatable :: all_tests(:)
        integer(int32) :: i, j
        logical :: found

        all_tests = get_all_tests()

        do i = 1, size(test_names)
            found = .false.
            do j = 1, size(all_tests)
                if (trim(test_names(i)) == trim(all_tests(j)%name)) then
                    call all_tests(j)%test_proc()
                    print "(' ',A,' passed.')", trim(test_names(i))
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                print *, "Unknown test: ", trim(test_names(i))
            end if
        end do
    end subroutine run_named_tests_tox_data_integration

    subroutine test_check_plateau_condition_helper()
        integer(int32), parameter :: n_studies = 3
        real(real64) :: confidence_interval(2, n_studies)
        real(real64) :: best_ci(2, n_studies)
        integer(int32) :: best_i_point, best_i_neighbor, best_exceeded
        logical :: plateau_found
        real(real64) :: succeeding_ci_overlap
        integer(int32) :: join_method

        succeeding_ci_overlap = 0.5_real64

        ! ============================================================
        ! JOIN_MIN: success case
        ! All overlaps >= 0.5 → min(overlaps) >= 0.5 → plateau_found = .true.
        !
        ! best_ci (all studies): [0.0, 1.0]
        ! confidence_interval:
        !   study 1: [0.2, 0.8]  → overlap = 0.6 / 0.6 = 1.0
        !   study 2: [0.0, 0.5]  → overlap = 0.5 / 0.5 = 1.0
        !   study 3: [0.5, 1.0]  → overlap = 0.5 / 0.5 = 1.0
        ! min(overlaps) = 1.0 ≥ 0.5 → plateau condition met
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            0.2, 0.8, 0.0, &   ! mins
            0.5, 0.1, 0.2  &    ! maxs
        ], shape(confidence_interval))

        best_i_point    = 1
        best_i_neighbor = 1
        best_exceeded   = 0

        join_method = JOIN_MIN

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=2, i_neighbor_count=3, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_true(plateau_found, "test_check_plateau_condition_helper: JOIN_MIN success: all overlaps >= threshold")
        call assert_equal_int(best_i_point, 2, "test_check_plateau_condition_helper: JOIN_MIN success: best_i_point updated")
        call assert_equal_int(best_i_neighbor, 3, "test_check_plateau_condition_helper: JOIN_MIN success: best_i_neighbor updated")
        call assert_equal_int(best_exceeded, 3, "test_check_plateau_condition_helper: JOIN_MIN success: all 3 exceeded")

        ! ============================================================
        ! JOIN_MIN: failure case
        ! Not all overlaps >= 0.5 → min(overlaps) < 0.5 → plateau_found = .false.
        !
        ! best_ci (all studies): [0.0, 1.0]
        ! confidence_interval:
        !   study 1: [0.2, 0.8]  → overlap = 1.0
        !   study 2: [0.0, 0.5]  → overlap = 1.0
        !   study 3: [1.1, 1.2]  → overlap = 0.0
        ! overlaps = [1.0, 1.0, 0.0], min = 0.0 < 0.5 → plateau condition NOT met
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            0.2, 0.8, 0.0, &   ! mins
            0.5, 1.1, 1.2  &    ! maxs
        ], shape(confidence_interval))

        best_i_point    = 10
        best_i_neighbor = 20
        best_exceeded   = 0

        join_method = JOIN_MIN

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=4, i_neighbor_count=5, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_false(plateau_found, "test_check_plateau_condition_helper: JOIN_MIN failure: not all overlaps >= threshold")
        call assert_equal_int(best_exceeded, 2, "test_check_plateau_condition_helper: JOIN_MIN failure: only 2 exceeded")

        ! ============================================================
        ! JOIN_MAX: success case
        ! At least one overlap >= 0.5 → max(overlaps) >= 0.5 → plateau_found = .true.
        !
        ! best_ci: [0.0, 1.0]
        ! confidence_interval:
        !   study 1: [1.1, 1.2] → overlap = 0.0
        !   study 2: [1.3, 1.4] → overlap = 0.0
        !   study 3: [0.2, 0.8] → overlap = 1.0
        ! max(overlaps) = 1.0 ≥ 0.5 → plateau condition met
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            1.1, 1.2, 1.3, &   ! mins
            1.4, 0.2, 0.8   &   ! maxs
        ], shape(confidence_interval))

        best_i_point    = 1
        best_i_neighbor = 1
        best_exceeded   = 0

        join_method = JOIN_MAX

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=6, i_neighbor_count=7, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_true(plateau_found, "test_check_plateau_condition_helper: JOIN_MAX success: at least one overlap >= threshold")
        call assert_equal_int(best_exceeded, 1, "test_check_plateau_condition_helper: JOIN_MAX success: exactly 1 exceeded")

        ! ============================================================
        ! JOIN_MAX: failure case
        ! No overlap >= 0.5 → max(overlaps) < 0.5 → plateau_found = .false.
        !
        ! best_ci: [0.0, 1.0]
        ! confidence_interval:
        !   study 1: [1.1, 1.2] → overlap = 0.0
        !   study 2: [1.3, 1.4] → overlap = 0.0
        !   study 3: [1.5, 1.6] → overlap = 0.0
        ! max(overlaps) = 0.0 < 0.5 → plateau condition NOT met
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            1.1, 1.2, 1.3, &   ! mins
            1.4, 1.4, 1.6  &    ! maxs
        ], shape(confidence_interval))

        best_i_point    = 1
        best_i_neighbor = 1
        best_exceeded   = 0

        join_method = JOIN_MAX

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=8, i_neighbor_count=9, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_false(plateau_found, "test_check_plateau_condition_helper: JOIN_MAX failure: no overlap >= threshold")
        call assert_equal_int(best_exceeded, 0, "test_check_plateau_condition_helper: JOIN_MAX failure: none exceeded")

        ! ============================================================
        ! JOIN_MEDIAN: success case
        ! Majority (>= 2 of 3) overlaps >= 0.5 → median(overlaps) >= 0.5 → plateau_found = .true.
        !
        ! best_ci: [0.0, 1.0]
        ! confidence_interval:
        !   study 1: [0.2, 0.8] → overlap = 1.0
        !   study 2: [0.1, 0.9] → overlap = 1.0
        !   study 3: [1.1, 1.2] → overlap = 0.0
        ! overlaps = [1.0, 1.0, 0.0], median = 1.0 ≥ 0.5 → plateau condition met
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            0.2, 0.8, 0.1, &   ! mins
            0.9, 1.1, 1.2  &    ! maxs
        ], shape(confidence_interval))

        best_i_point    = 1
        best_i_neighbor = 1
        best_exceeded   = 0

        join_method = JOIN_MEDIAN

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=11, i_neighbor_count=12, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_true(plateau_found, "test_check_plateau_condition_helper: JOIN_MEDIAN success: majority overlaps >= threshold")
        call assert_equal_int(best_exceeded, 2, "test_check_plateau_condition_helper: JOIN_MEDIAN success: 2 exceeded")

        ! ============================================================
        ! JOIN_MEDIAN: failure case
        ! Only 1 of 3 overlaps >= 0.5 → median(overlaps) < 0.5 → plateau_found = .false.
        !
        ! best_ci: [0.0, 1.0]
        ! confidence_interval:
        !   study 1: [0.2, 0.8] → overlap = 1.0
        !   study 2: [1.1, 1.2] → overlap = 0.0
        !   study 3: [1.3, 1.4] → overlap = 0.0
        ! overlaps = [1.0, 0.0, 0.0], median = 0.0 < 0.5 → plateau condition NOT met
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            0.2, 0.8, 1.1, &   ! mins
            1.2, 1.3, 1.4 &     ! maxs
        ], shape(confidence_interval))

        best_i_point    = 1
        best_i_neighbor = 1
        best_exceeded   = 0

        join_method = JOIN_MEDIAN

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=13, i_neighbor_count=14, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_false(plateau_found, "test_check_plateau_condition_helper: JOIN_MEDIAN failure: minority overlaps >= threshold")
        call assert_equal_int(best_exceeded, 1, "test_check_plateau_condition_helper: JOIN_MEDIAN failure: only 1 exceeded")

        ! ============================================================
        ! Extra: worse candidate than previous best → early plateau
        ! exceeds_min_CI_overlap < best_params_exceeded_CI_overlap
        ! → plateau_found = .true., best_* NOT updated
        ! JOIN_MIN condition alone would fail
        ! ============================================================
        best_ci = reshape([ &
            0.0, 1.0, 0.0, &   ! mins
            1.0, 0.0, 1.0&      ! maxs
        ], shape(best_ci))

        confidence_interval = reshape([ &
            1.1, 1.2, 0.2, &   ! mins (only 2nd,3rd overlap)
            0.3, 0.8, 0.9&      ! maxs
        ], shape(confidence_interval))

        best_i_point        = 99
        best_i_neighbor     = 88
        best_exceeded       = 3   ! previous best had all 3 exceeding
        join_method         = JOIN_MIN

        call check_plateau_condition_helper( &
            confidence_interval, best_ci, n_studies, &
            best_i_point, best_i_neighbor, best_exceeded, &
            i_point_count=15, i_neighbor_count=16, &
            join_method=join_method, &
            succeeding_ci_overlap=succeeding_ci_overlap, &
            plateau_found=plateau_found)

        call assert_true(plateau_found, "test_check_plateau_condition_helper: Early plateau: new candidate worse than previous best")
        call assert_equal_int(best_i_point, 99, "test_check_plateau_condition_helper: Early plateau: best_i_point unchanged")
        call assert_equal_int(best_i_neighbor, 88, "test_check_plateau_condition_helper: Early plateau: best_i_neighbor unchanged")
        call assert_equal_int(best_exceeded, 3, "test_check_plateau_condition_helper: Early plateau: best_exceeded unchanged")

    end subroutine test_check_plateau_condition_helper

    subroutine test_gjct_permutation_test_helper()
        integer(int32), parameter :: n_bins = 3, n_points = 2, n_studies = 2, n_permutations = 50
        integer(int32) :: i_permutation, i_study, max_possible_seen, n_counts, i_point, i_bin, i_jsd

        integer(int32) :: mean_pmf_counts(n_bins, n_points), tmp_mean_pmf_counts(n_bins, n_points)
        integer(int32) :: mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)

        integer(int32) :: tmp_counts(n_bins, n_points), expected_tmp_counts(n_bins, n_points)
        real(real64) :: tmp_pmfs(n_bins, n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points)
        real(real64) :: tmp_js_divergences(n_points, n_studies)
        real(real64) :: tmp_weights(n_points, n_studies)
        real(real64) :: tmp_global_js_divergence(n_studies)
        real(real64) :: p_values(n_studies), global_jsd_observed(n_studies)

        ! ------------------------------------------------------------
        ! Setup deterministic small test data
        ! ------------------------------------------------------------
        mean_pmf_counts = reshape([ &
            10, 10, &   ! point 1
            10, 10, &   ! point 2
            10, 10&      ! point 3
        ], shape(mean_pmf_counts))

        mean_pmf = real(mean_pmf_counts, real64) / real(n_studies, real64)

        mean_pmf_included_n_reps = [30_int32, 30_int32]
        included_n_reps = reshape([ &
            15, 15, &   ! study 1
            15, 15&      ! study 2
        ], shape(included_n_reps))

        do i_jsd = 0, 1
            global_jsd_observed = real(i_jsd, real64)
            do i_permutation = 1, n_permutations
                call gjct_permutation_test_helper(&
                    mean_pmf_counts, mean_pmf, mean_pmf_included_n_reps, n_points, n_bins,&
                    included_n_reps, n_studies, global_jsd_observed, p_values,&
                    tmp_pmfs, tmp_counts, tmp_mean_pmf_counts, tmp_js_divergences, tmp_weights, tmp_global_js_divergence, i_permutation, random_seed=12345_int32&
                )

                do i_point = 1, n_points
                    do i_bin = 1, n_bins
                        associate (&
                            pmf_count => tmp_counts(i_bin, i_point),&
                            pmf_included_n_reps => included_n_reps(i_point, n_studies)&
                        )
                            call assert_true(pmf_count <= pmf_included_n_reps, "test_gjct_permutation_test_helper: resampled bins should never exceed the total count of the reference point (no replacement)")
                        end associate
                    end do
                end do
            end do

            call assert_equal_int(count(p_values>abs(real(1 - i_jsd, real64)-TOL)), n_studies, "test_gjct_permutation_test_helper: p_values should be either 1 or 0 if determined correctly for jsd_obs=0 and jsd_obs=1")
        end do
    end subroutine test_gjct_permutation_test_helper

    subroutine test_bootstrap_histogram_helper()
        integer(int32), parameter :: n_bins = 3, n_points = 2, n_studies = 2, n_bootstraps = 50
        integer(int32), parameter :: n_bootstrapping_top_k_jsds = n_bootstraps
        integer(int32) :: i_bootstrap, i_study, max_possible_seen, n_counts, i_point, i_bin

        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)

        real(real64) :: confidence_interval(2, n_studies)
        real(real64) :: bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies)

        integer(int32) :: tmp_counts(n_bins, n_points), expected_tmp_counts(n_bins, n_points)
        real(real64) :: tmp_pmfs(n_bins, n_points, n_studies)
        real(real64) :: tmp_mean_pmf(n_bins, n_points)
        real(real64) :: tmp_js_divergences(n_points, n_studies)
        real(real64) :: tmp_weights(n_points, n_studies)
        real(real64) :: tmp_global_js_divergence(n_studies)

        integer(int32) :: i_saw_greater, i_saw_lower
        real(real64) :: p

        ! ------------------------------------------------------------
        ! Setup deterministic small test data
        ! ------------------------------------------------------------
        mean_pmf_counts = reshape([ &
            10, 10, &   ! point 1
            10, 10, &   ! point 2
            10, 10&      ! point 3
        ], shape(mean_pmf_counts))

        mean_pmf_included_n_reps = [30_int32, 30_int32]
        included_n_reps = reshape([ &
            15, 15, &   ! study 1
            15, 15&      ! study 2
        ], shape(included_n_reps))

        confidence_interval = 0.0_real64

        i_saw_greater = 0
        i_saw_lower   = 0
        ! ------------------------------------------------------------
        ! Test each bootstrap iteration independently
        ! ------------------------------------------------------------

        do i_bootstrap = 1, n_bootstraps
            confidence_interval(1, :) = M_POS_INF
            confidence_interval(2, :) = M_NEG_INF
            bootstrapping_top_k_jsds = 0.0_real64

            call bootstrap_histogram_helper( &
                i_bootstrap, included_n_reps, n_bins, n_points, n_studies, &
                mean_pmf_counts, mean_pmf_included_n_reps, confidence_interval, &
                tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                bootstrapping_top_k_jsds, n_bootstrapping_top_k_jsds, &
                tmp_counts, tmp_pmfs, tmp_mean_pmf, random_seed=12345_int32)

            ! --------------------------------------------------------
            ! 1. Confidence interval updated correctly
            ! --------------------------------------------------------
            do i_study = 1, n_studies
                call assert_equal_int(count(ieee_is_finite(bootstrapping_top_k_jsds(:, 1, i_study))), i_bootstrap, "test_bootstrap_histogram_helper: should always push into the top k")
                call assert_equal_int(count(ieee_is_finite(bootstrapping_top_k_jsds(:, 2, i_study))), i_bootstrap, "test_bootstrap_histogram_helper: should always push into the bottom k")
                call assert_equal_real(maxval(bootstrapping_top_k_jsds(:, 1, i_study)), confidence_interval(1, i_study), 0.0_real64, "test_bootstrap_histogram_helper: confidence_interval should have minimum value")
                call assert_equal_real(minval(bootstrapping_top_k_jsds(:, 2, i_study)), confidence_interval(2, i_study), 0.0_real64, "test_bootstrap_histogram_helper: confidence_interval should have maximum value")
            end do

            ! --------------------------------------------------------
            ! 2. tmp_mean_pmf equals mean of tmp_pmfs across studies
            ! --------------------------------------------------------
            call assert_equal_array_real(tmp_mean_pmf, sum(tmp_pmfs, dim=3) / n_studies, size(tmp_mean_pmf, kind=int32), TOL, "test_bootstrap_histogram_helper: tmp_mean_pmf equals mean of tmp_pmfs")

            do i_point = 1, n_points
                do i_bin = 1, n_bins
                    associate (&
                        mean_pmf_count => real(mean_pmf_counts(i_bin, i_point), real64),&
                        pmf_count => real(tmp_counts(i_bin, i_point), real64),&
                        mpmf_included_n_reps => real(mean_pmf_included_n_reps(i_point), real64),&
                        pmf_included_n_reps => real(included_n_reps(i_point, n_studies), real64)&
                    )
                        p = mean_pmf_count / mpmf_included_n_reps
                        i_saw_lower = i_saw_lower + merge(1, 0, pmf_count < p * pmf_included_n_reps)
                        i_saw_greater = i_saw_greater + merge(1, 0, pmf_count > p * pmf_included_n_reps)
                    end associate
                end do
            end do
        end do

        expected_tmp_counts = reshape([ &
            4, 4, &   ! point 1
            7, 5, &   ! point 2
            5, 5&      ! point 3
        ], shape(expected_tmp_counts))
        n_counts = size(expected_tmp_counts, kind=int32)

        call assert_equal_array_int(tmp_counts, expected_tmp_counts, n_counts, "test_bootstrap_histogram_helper: reproducibility")
        max_possible_seen = n_counts * n_bootstraps

        ! All mean pmf counts are equal -> equal probability for greater/lower. Equality is less likely
        call assert_true(i_saw_greater > max_possible_seen / (n_counts / 2), "test_bootstrap_histogram_helper: Multinomial sampling produced many greater counts")
        call assert_true(i_saw_lower > max_possible_seen / (n_counts / 2), "test_bootstrap_histogram_helper: Multinomial sampling produced many lower counts")
    end subroutine test_bootstrap_histogram_helper

    subroutine test_test_mean_pmf_min_counts_helper()
        integer(int32), parameter :: n_bins = 3, n_points = 2
        integer(int32) :: min
        integer(int32) :: counts_ok(n_bins, n_points)
        integer(int32) :: counts_bad(n_bins, n_points)
        logical :: result

        ! Case 1: all bins meet the minimum
        min = 5
        counts_ok = reshape([ &
            5, 6, 7, &   ! point 1
            8, 5, 9&     ! point 2
        ], shape(counts_ok))

        result = test_mean_pmf_min_counts_helper(counts_ok, n_bins, n_points, min)
        call assert_true(result, "test_test_mean_pmf_min_counts_helper: All bins should meet minimum count")

        ! Case 2: at least one bin violates the minimum
        counts_bad = counts_ok
        counts_bad(2,1) = 3   ! violate minimum

        result = test_mean_pmf_min_counts_helper(counts_bad, n_bins, n_points, min)
        call assert_false(result, "test_test_mean_pmf_min_counts_helper: A bin below minimum should cause failure")

    end subroutine test_test_mean_pmf_min_counts_helper

    subroutine test_create_mean_pmf_helpers()
        integer(int32), parameter :: n_bins = 2, n_points = 1, n_studies = 2

        ! Inputs
        real(real64) :: pmfs(n_bins, n_points, n_studies)
        integer(int32) :: counts(n_bins, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)

        ! Outputs
        real(real64) :: mean_pmf(n_bins, n_points)
        real(real64) :: mean_pmf_only(n_bins, n_points)
        integer(int32) :: mean_pmf_included_n_reps(n_points)
        integer(int32) :: mean_pmf_counts(n_bins, n_points)

        ! Expected values
        real(real64) :: expected_mean_pmf(n_bins)
        integer(int32) :: expected_counts(n_bins)
        integer(int32) :: expected_included

        ! ------------------------------------------------------------
        ! Define tiny deterministic test data
        ! ------------------------------------------------------------
        !
        ! Study 1 PMF: [0.2, 0.8]
        ! Study 2 PMF: [0.4, 0.6]
        !
        ! Expected mean PMF = [0.3, 0.7]
        !
        pmfs(:,1,1) = [0.2_real64, 0.8_real64]
        pmfs(:,1,2) = [0.4_real64, 0.6_real64]

        counts(:,1,1) = [2, 8]
        counts(:,1,2) = [4, 6]

        included_n_reps(1,1) = 10
        included_n_reps(1,2) = 12

        expected_mean_pmf = [0.3_real64, 0.7_real64]
        expected_counts   = [6, 14]
        expected_included = 22

        ! ------------------------------------------------------------
        ! Test full helper
        ! ------------------------------------------------------------
        call create_mean_pmf_helper( pmfs, counts, n_bins, n_points, n_studies, &
                                     included_n_reps, mean_pmf, &
                                     mean_pmf_included_n_reps, mean_pmf_counts )

        call assert_equal_array_real(mean_pmf(:,1), expected_mean_pmf, n_bins, TOL, &
             "test_create_mean_pmf_helpers: mean_pmf mismatch")

        call assert_equal_array_int(mean_pmf_counts(:,1), expected_counts, n_bins, &
             "test_create_mean_pmf_helpers: mean_pmf_counts mismatch")

        call assert_equal_int(mean_pmf_included_n_reps(1), expected_included, &
             "test_create_mean_pmf_helpers: mean_pmf_included_n_reps mismatch")

        ! ------------------------------------------------------------
        ! Test PMF-only helper
        ! ------------------------------------------------------------
        call create_mean_pmf_only_helper(pmfs, n_bins, n_points, n_studies, mean_pmf_only)

        call assert_equal_array_real(mean_pmf_only(:,1), expected_mean_pmf, n_bins, TOL, &
             "test_create_mean_pmf_helpers: mean_pmf_only mismatch")

        ! ------------------------------------------------------------
        ! Cross-check: both helpers must produce identical mean_pmf
        ! ------------------------------------------------------------
        call assert_equal_array_real(mean_pmf_only(:,1), mean_pmf(:,1), n_bins, TOL, &
             "test_create_mean_pmf_helpers: mean_pmf_only != mean_pmf")

    end subroutine test_create_mean_pmf_helpers

    subroutine test_test_neighborhood_overlaps()
        integer(int32), parameter :: n_points = 5
        real(real64)   :: min_overlap
        logical        :: ok
        integer(int32) :: nr(2,n_points)

        !===============================================================
        ! 1) Simple valid case: all consecutive neighborhoods overlap
        !===============================================================
        nr(:,1) = [1, 5] ! 2/4
        nr(:,2) = [3, 7] ! 1/4
        nr(:,3) = [6, 9] ! 1/3
        nr(:,4) = [8,12]
        min_overlap = 1.0_real64 / 8.0_real64

        ok = test_neighborhood_overlaps_helper(nr, 4_int32, min_overlap)
        call assert_true(ok, "test_test_neighborhood_overlaps: simple valid case: all overlaps large enough")

        !===============================================================
        ! 2) Exact-boundary overlap (touching at one point)
        !===============================================================
        nr(:,1) = [1, 6]
        nr(:,2) = [5,10]
        min_overlap = 1.0_real64 / 5.0_real64
        ok = test_neighborhood_overlaps_helper(nr, 2_int32, min_overlap)
        call assert_true(ok, "test_test_neighborhood_overlaps: exact boundary overlap")

        ok = test_neighborhood_overlaps_helper(nr, 2_int32, above(min_overlap))
        call assert_false(ok, "test_test_neighborhood_overlaps: exact boundary overlap if min_overlap slightly to high")

        !===============================================================
        ! 3) Zero overlap
        !===============================================================
        nr(:,1) = [1, 4]
        nr(:,2) = [6, 9]
        min_overlap = above(0.0_real64)
        ok = test_neighborhood_overlaps_helper(nr, 2_int32, min_overlap)
        call assert_false(ok, "test_test_neighborhood_overlaps: zero overlap should fail")

        !===============================================================
        ! 7) Duplicate identical ranges → full overlap
        !===============================================================
        nr(:,1) = [2,8]
        nr(:,2) = [2,8]
        nr(:,3) = [2,8]
        nr(:,4) = [2,8]
        min_overlap = 1.0_real64
        ok = test_neighborhood_overlaps_helper(nr, 4_int32, min_overlap)
        call assert_true(ok, "test_test_neighborhood_overlaps: identical ranges satisfy overlap")

        !===============================================================
        ! 8) Single-point case (n_points = 1)
        !    → trivially true (no consecutive pairs)
        !===============================================================
        min_overlap = 1.0_real64
        ok = test_neighborhood_overlaps_helper(nr, 1_int32, min_overlap)
        call assert_true(ok, "test_test_neighborhood_overlaps: single-point case should always pass")

        !===============================================================
        ! 9) Large ranges with partial overlaps
        !===============================================================
        nr(:,1) = [1,100] ! 50/99 = 0.5050505050505051
        nr(:,2) = [50,150] ! 30/100 = 0.3
        nr(:,3) = [120,200] ! 20/80 = 0.25
        nr(:,4) = [180,250] ! 10/70 = 0.14285714285714285
        nr(:,5) = [240,300]

        min_overlap = 10.0_real64 / 70.0_real64
        ok = test_neighborhood_overlaps_helper(nr, 5_int32, min_overlap)
        call assert_true(ok, "test_test_neighborhood_overlaps: large ranges: all overlaps large enough")

        min_overlap = above(min_overlap)
        ok = test_neighborhood_overlaps_helper(nr, 5_int32, min_overlap)
        call assert_false(ok, "test_test_neighborhood_overlaps: large ranges: insufficient overlap for 0.2")

        !===============================================================
        ! 11) min_overlap=0 succeeds
        !===============================================================
        nr(:,1) = [1,10]
        nr(:,2) = [100,200]
        min_overlap = 0.0_real64
        ok = test_neighborhood_overlaps_helper(nr, 2, min_overlap)
        call assert_true(ok, "test_test_neighborhood_overlaps: touching ranges satisfy min_overlap=0")

        min_overlap = above(0.0_real64)
        ok = test_neighborhood_overlaps_helper(nr, 2, min_overlap)
        call assert_false(ok, "test_test_neighborhood_overlaps: touching ranges insufficient for min_overlap=1")

        !===============================================================
        ! 12) Very large min_overlap (impossible)
        !===============================================================
        nr(:,1) = [1,10]
        nr(:,2) = [5,15]
        min_overlap = M_POS_INF
        ok = test_neighborhood_overlaps_helper(nr, 2, min_overlap)
        call assert_false(ok, "test_test_neighborhood_overlaps: min_overlap too large → fail")
    end subroutine test_test_neighborhood_overlaps

    subroutine test_compute_fractional_overlap()
        real(real64) :: a(2), b(2), overlap, overlap_rev

        ! Touching intervals a1 a2=b1 b2
        a = [-1.0_real64, 666.123_real64]
        b = [a(2), huge(1.0_real64)]
        overlap = compute_fractional_overlap_helper(a(1), a(2), b(1), b(2))
        overlap_rev = compute_fractional_overlap_helper(b(1), b(2), a(1), a(2))
        call assert_equal_real(overlap, 0.0_real64, TOL, "test_compute_fractional_overlap: 1. touching intervals a,b")
        call assert_equal_real(overlap_rev, 0.0_real64, TOL, "test_compute_fractional_overlap: 1. touching intervals b,a")

        ! Partial Overlap a1 b1 a2 b2
        a = [-1.0_real64, 666.123_real64]
        b = [0.123_real64, 667.0_real64]
        overlap = compute_fractional_overlap_helper(a(1), a(2), b(1), b(2))
        overlap_rev = compute_fractional_overlap_helper(b(1), b(2), a(1), a(2))
        call assert_equal_real(overlap, 666.0_real64 / 667.123_real64, TOL, "test_compute_fractional_overlap: 2. partial overlap a,b")
        call assert_equal_real(overlap_rev, 666.0_real64 / 666.877_real64, TOL, "test_compute_fractional_overlap: 2. partial overlap b,a")

        ! Inclusion, a1 b1 b2 a2
        a = [-1.0_real64, 667.0_real64]
        b = [0.123_real64, 666.123_real64]
        overlap = compute_fractional_overlap_helper(a(1), a(2), b(1), b(2))
        overlap_rev = compute_fractional_overlap_helper(b(1), b(2), a(1), a(2))
        call assert_equal_real(overlap, 666.0_real64 / 668.0_real64, TOL, "test_compute_fractional_overlap: 3. inclusion a,b")
        call assert_equal_real(overlap_rev, 1.0_real64, TOL, "test_compute_fractional_overlap: 3. inclusion b,a")

        ! Full Overlap, a1=b1 b2=a2
        a = [-1.0_real64, 667.123_real64]
        b = a
        overlap = compute_fractional_overlap_helper(a(1), a(2), b(1), b(2))
        overlap_rev = compute_fractional_overlap_helper(b(1), b(2), a(1), a(2))
        call assert_equal_real(overlap, 1.0_real64, TOL, "test_compute_fractional_overlap: 4. inclusion a,b")
        call assert_equal_real(overlap_rev, 1.0_real64, TOL, "test_compute_fractional_overlap: 4. inclusion b,a")

        ! No Overlap, a1 a2 b1 b2
        a = [-1.0_real64, 667.123_real64]
        b = [2 * a(2), 12 * a(2)]
        overlap = compute_fractional_overlap_helper(a(1), a(2), b(1), b(2))
        overlap_rev = compute_fractional_overlap_helper(b(1), b(2), a(1), a(2))
        call assert_equal_real(overlap, 0.0_real64, TOL, "test_compute_fractional_overlap: 5. inclusion a,b")
        call assert_equal_real(overlap_rev, 0.0_real64, TOL, "test_compute_fractional_overlap: 5. inclusion b,a")

    end subroutine test_compute_fractional_overlap

    subroutine test_estimate_bin_count()
        ! Inputs
        integer(int32) :: n_neighbors, n_residuals, max_n_reps_all_studies
        real(real64), allocatable :: residuals(:)
        real(real64) :: shared_residual_range
        integer(int32) :: n_bins, ierr, n_bins_ref, i

        !===============================
        ! Start of tests
        !===============================
        !----------------------------------------------------
        ! 1) Freedman-Diaconis case: Freedman-Diaconis max
        !----------------------------------------------------
        n_neighbors             = 5_int32
        n_residuals             = 8_int32; allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        shared_residual_range   = 3.0_real64

        ! IQR = (2-0.5)=1.5 -> Freed = 3 / (1.5/(40^(1/3))) = 3/0.43860266073192994 = nint(6.839903786706787) = 7
        ! Sturges = 1 + nint(ln(40)/ln(2)) = 1 + 5 = 6
        residuals = [0.0_real64, 0.0_real64, 1.0_real64, 1.0_real64, &
                   2.0_real64, 2.0_real64, 3.0_real64, 3.0_real64]

        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_estimate_bin_count: Case 1: freedman case, ierr")
        call assert_equal_int(n_bins, 7_int32, "test_estimate_bin_count: Case 1: freedman case, n_bins")


        !----------------------------------------------------
        ! 2) Very small sample size (n_residuals = 1)
        !    -> Freedman bin width is zero -> Sturges fallback = 1 + log(1) / log(2) = 1 + 2 = 3
        !----------------------------------------------------
        n_neighbors             = 3_int32
        n_residuals             = 1_int32; deallocate(residuals); allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        residuals = 0.0_real64
        shared_residual_range   = 0.0_real64

        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_estimate_bin_count: Case 2: one residual case, ierr")
        call assert_equal_int(n_bins, 3_int32, "test_estimate_bin_count: Case 2: one residual case, n_bins")


        !----------------------------------------------------
        ! 3) Invalid parameters: non-positive n_neighbors or max_n_reps_all_studies
        !----------------------------------------------------
        n_neighbors             = 0_int32
        n_residuals             = 5_int32; deallocate(residuals); allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        residuals = [(-2.0_real64), -1.0_real64, 0.0_real64, 1.0_real64, 2.0_real64]
        shared_residual_range   = 4.0_real64

        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, 5_int32), "test_estimate_bin_count: Case 3: expected error for invalid neighbor count")

        n_neighbors             = 3_int32
        max_n_reps_all_studies  = 0_int32

        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 2_int32), "test_estimate_bin_count: Case 3: expected error for invalid max_n_reps_all_studies value (zero)")

        !----------------------------------------------------
        ! 4) Invalid shared_residual_range (<= 0) with non-constant data
        !----------------------------------------------------
        n_neighbors             = 4_int32
        max_n_reps_all_studies  = 20_int32
        n_residuals             = 4_int32; deallocate(residuals); allocate(residuals(n_residuals))
        residuals = [0.0_real64, 1.0_real64, 2.0_real64, 3.0_real64]
        shared_residual_range   = -1.0_real64

        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, 6_int32), "test_estimate_bin_count: Case 4: expected error for negative shared_residual_range")


        !----------------------------------------------------
        ! 5) NaNs in residuals -> NaNs ignored
        !----------------------------------------------------
        ! NaN not ignored:
        !     Sturges = 1+log(5*50)/log(2) = 3
        !     Freed = NaN  (interplation with NaN)
        ! NaN ignored:
        !     Sturges = 1+log(5*50)/log(2) = 3
        !     Freed: IQR = 2.25 - 0.75 = 1.5 -> nint(4 / ( 1.5/((5*50)**(1/3)) )) = 17
        n_neighbors             = 50_int32
        n_residuals             = 5_int32; deallocate(residuals); allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        shared_residual_range   = 4.0_real64
        residuals = [0.0_real64, 1.0_real64, M_NAN, &
                   2.0_real64, 3.0_real64]

        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)
        call assert_equal_int(n_bins, 17_int32, "test_estimate_bin_count: Case 5: NaN case, n_bins")

        !----------------------------------------------------
        ! 6) Monotonicity in effective sample size:
        !    increasing max_n_reps_all_studies * n_neighbors
        !    should not decrease n_bins
        !----------------------------------------------------
        n_residuals             = 20_int32; deallocate(residuals); allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        residuals = [(real(i, real64), i=1,n_residuals)]
        shared_residual_range   = real(n_residuals-1, real64)

        n_neighbors             = 5_int32
        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins_ref, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_estimate_bin_count: Case 6: monotonicity base: ierr == 0")

        n_residuals             = 200_int32; deallocate(residuals); allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        residuals = [(real(i, real64), i=1,n_residuals)]
        shared_residual_range   = real(n_residuals-1, real64)
        n_neighbors             = 5_int32
        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_estimate_bin_count: Case 6: monotonicity larger neighborhood size: ierr == 0")
        call assert_true(n_bins >= n_bins_ref, "test_estimate_bin_count: Case 6: n_bins should not decrease with larger neighborhood size")


        !----------------------------------------------------
        ! 7) Very large effective sample size: check that Sturges grows ~ log2(n_eff)
        !----------------------------------------------------
        n_residuals             = 1000000_int32; deallocate(residuals); allocate(residuals(n_residuals))
        max_n_reps_all_studies  = n_residuals
        residuals = [(real(i, real64), i=1,n_residuals)]
        shared_residual_range   = 0.0_real64 ! -> Freedman is zero

        ! Sturges = nint(1+log(10000000)/log(2)) = 24
        n_neighbors             = 10_int32
        call estimate_bin_count_alloc(residuals, max_n_reps_all_studies, 1_int32, 1_int32, &
                                    n_neighbors, shared_residual_range, n_bins, ierr)
        call assert_equal_int(ierr, 0_int32, "test_estimate_bin_count: Case 7: large n_eff: ierr == 0")

        call assert_equal_int(n_bins, 24_int32, "test_estimate_bin_count: Case 7: large n_eff: n_bins >= Sturges")
    end subroutine test_estimate_bin_count

    subroutine test_determine_shared_residual_range
        integer(int32), parameter :: n_reps_S1 = 4, n_reps_S2 = 3, n_neighbors = 2, n_points = 2, n_studies = 2, max_n_reps_all_studies = max(n_reps_S1, n_reps_S2)
        real(real64), dimension(max_n_reps_all_studies, n_neighbors, n_points, 2), target :: S
        real(real64), dimension(:), pointer :: S_flat
        real(real64) :: R
        integer(int32) :: ierr, n_S
        real(real64) :: q
        n_S = size(S, kind=int32)

        ! ============================================================
        ! Test 1 — Basic correctness with simple values
        ! ============================================================
        !
        S = reshape([ &
            1.0_real64,2.0_real64,3.0_real64,4.0_real64,&
            5.0_real64,6.0_real64,-7.0_real64,8.0_real64,&
            9.0_real64,10.0_real64,11.0_real64,12.0_real64,&
            1.0_real64,1.0_real64,1.0_real64,1.0_real64,&
            2.0_real64,-4.0_real64,6.0_real64,M_NAN,&
            8.0_real64,1.0_real64,3.0_real64,M_NAN,&
            5.0_real64,7.0_real64,9.0_real64,M_NAN,&
            0.0_real64,1.0_real64,2.0_real64,M_NAN  ], shape(S))
        S_flat(1:n_S) => S

        ! Sorted pool:
        !   [0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11, 12]
        !
        ! 95% quantile → 0.95 * (28 - 1) + 1 = 26.65
        ! sorted(26) = 10
        ! sorted(27) = 11
        !
        ! Expected R = 10 + 0.65 * (11-10) = 10.65
        !

        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 1: ierr should be OK")
        call assert_equal_real(R, 10.65_real64, TOL, "test_determine_shared_residual_range: Test 1: R should be 10.65")

        ! ============================================================
        ! Test 2 — Custom quantile (50%)
        ! ============================================================
        !
        ! Median of sorted array above = 0.5 * (sorted(12) + sorted(13)) = 5.5
        !
        q = 50.0_real64
        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr, q)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 2: ierr should be OK")
        call assert_equal_real(R, 4.0_real64, TOL, "test_determine_shared_residual_range: Test 2: R should be 4.0")

        ! ============================================================
        ! Test 3 — Quantile < 0 → error
        ! ============================================================
        q = below(0.0_real64)
        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr, q)
        call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, 5_int32), "test_determine_shared_residual_range: Test 3: ierr should be INVALID_INPUT")

        ! ============================================================
        ! Test 4 — Quantile > 100 → error
        ! ============================================================
        q = above(100.0_real64)
        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr, q)
        call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, 5_int32), "test_determine_shared_residual_range: Test 4: ierr should be INVALID_INPUT")

        ! ============================================================
        ! Test 5 — NaNs must be ignored
        ! ============================================================
        !
        ! Replace some values with NaN; remaining values should determine R.
        !
        S = reshape([ &
            M_NAN, 2.0_real64, 3.0_real64, 4.0_real64, &
            5.0_real64, 6.0_real64, -7.0_real64, 8.0_real64, &
            9.0_real64, 10.0_real64, -11.0_real64, 12.0_real64, &
            1.0_real64,1.0_real64,1.0_real64,1.0_real64, &
            -1.0_real64, 2.0_real64, 3.0_real64, M_NAN,&
            4.0_real64, 5.0_real64, 6.0_real64, M_NAN,&
            -7.0_real64, 8.0_real64, 9.0_real64, M_NAN,&
            10.0_real64, -11.0_real64, M_NAN, M_NAN ], shape(S))

        ! Pool now excludes two NaNs → 26 values
        ! sorted = [1, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12]
        ! 95% quantile → 0.95*(26-1)+1=24.75
        ! sorted(24) = 11
        ! sorted(25) = 11
        ! -> R = 11
        !
        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 5: ierr should be OK")
        call assert_equal_real(R, 11.0_real64, TOL, "test_determine_shared_residual_range: Test 5: R should ignore NaNs")

        ! ============================================================
        ! Test 6 — All zeros
        ! ============================================================
        S = 0.0_real64
        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 6: ierr should be OK")
        call assert_equal_real(R, 0.0_real64, TOL, "test_determine_shared_residual_range: Test 6: R should be zero")

        ! ============================================================
        ! Test 7 — Single residual (n_reps_S1=1, n_reps_S2=1, n_neighbors=1, n_points=1)
        ! ============================================================
        S_flat(1) = 3.0_real64
        S_flat(2) = -4.0_real64
        ! sorted = [3, 4]
        ! rank = 0.95 * (2-1) + 1 = 1.95
        ! R = 3 + (4-3)*0.95 = 3.95
        call determine_shared_residual_range_alloc(S_flat, 2_int32, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 7: ierr should be OK")
        call assert_equal_real(R, 3.95_real64, TOL, "test_determine_shared_residual_range: Test 7: R should be 3.95")

        ! ============================================================
        ! Test 8 — Edge case, only NaN residuals
        ! ============================================================
        S_flat(1) = M_NAN
        call determine_shared_residual_range_alloc(S_flat, 1_int32, 1_int32, 1_int32, R, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 8: ierr should be OK")
        call assert_equal_real(R, 0.0_real64, 0.0_real64, "test_determine_shared_residual_range: Test 8: R should be 0.0")

        ! ============================================================
        ! Test 9 — Edge case, min/max quantile
        ! ============================================================
        S = reshape([ &
            M_NAN, 2.0_real64, 3.0_real64, 4.0_real64, &
            5.0_real64, 6.0_real64, -7.0_real64, 8.0_real64, &
            9.0_real64, 10.0_real64, -13.0_real64, 12.0_real64, &
            1.0_real64,1.0_real64,1.0_real64,1.0_real64, &
            -1.0_real64, 2.0_real64, 3.0_real64, M_NAN,&
            4.0_real64, 5.0_real64, 6.0_real64, M_NAN,&
            -7.0_real64, 8.0_real64, 9.0_real64, M_NAN,&
            10.0_real64, -11.0_real64, M_NAN, M_NAN ], shape(S))
        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr, 100.0_real64)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 9: max quantile ierr should be OK")
        call assert_equal_real(R, 13.0_real64, 0.0_real64, "test_determine_shared_residual_range: Test 9: max quantile should be max abs value")

        call determine_shared_residual_range_alloc(S_flat, n_S, 1_int32, 1_int32, R, ierr, 0.0_real64)
        call assert_equal_int(ierr, ERR_OK, "test_determine_shared_residual_range: Test 9: min quantile ierr should be OK")
        call assert_equal_real(R, 1.0_real64, 0.0_real64, "test_determine_shared_residual_range: Test 9: min quantile should be min abs value")

    end subroutine test_determine_shared_residual_range

    subroutine test_build_residual_histograms
        integer(int32), parameter :: n_reps = 3, n_neighbors = 2, n_points = 3
        integer(int32), parameter :: n_bins = 4
        integer(int32), dimension(n_neighbors, n_points) :: E
        real(real64), dimension(n_reps, n_neighbors * n_points) :: residuals
        real(real64), dimension(n_bins, n_points) :: pmf, expected_pmf
        integer(int32), dimension(n_bins, n_points) :: counts, expected_counts
        integer(int32), dimension(n_points) :: included
        real(real64) :: R
        integer(int32) :: ierr

        ! ============================================================
        ! Test 1 — Simple symmetric case, no NaNs
        ! ============================================================
        !
        ! R = 2, M = 4 → bin width w = 1
        ! Bins: [-2,-1), [-1,0), [0,1), [1,2]
        !
        R = 2.0_real64

        residuals(:, 1) = [-2.0, -0.5, 0.2]
        residuals(:, 2) = [1.7, 0.9, -1.2]
        residuals(:, 3) = 0.0_real64
        residuals(:, 4) = [2.5, -3.0, 1.2] ! (clamping applies -> [2,-2,1.2])
        residuals(:, 5) = [0.4, -0.1, 0.0]

        E(1,1) = 1
        E(2,1) = 2
        E(:,2) = 3
        E(1,3) = 4
        E(2,3) = 5

        call build_residual_histograms(E, n_neighbors, n_points, residuals, n_reps, 5_int32, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 1: ierr should be OK")

        ! point 1
        ! Values fall into bins:
        ! [-2,-1): -2, -1.2 → 2
        ! [-1,0): -0.5 → 1
        ! [0,1): 0.2, 0.9 → 2
        ! [1,2]: 1.7 → 1
        ! 
        ! point 2 — all zeros → all in bin [0,1)
        ! 
        ! point 3 — clamping:
        ! 2.5 → 2
        ! -3 → -2
        ! bins:
        ! [-2,-1): -2 → 1
        ! [-1,0): -0.1 → 1
        ! [0,1): 0.4, 0.0 → 2
        ! [1,2]: 1.2, 2 → 2
        expected_counts = reshape([&
            2, 1, 2, 1,&
            0, 0, 6, 0,&
            1, 1, 2, 2&
        ], [n_bins, n_points])
        expected_pmf = reshape([&
            0.3333333333333333_real64, 0.16666666666666666_real64, 0.3333333333333333_real64, 0.16666666666666666_real64,&
            0.0_real64, 0.0_real64, 1.0_real64, 0.0_real64,&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.3333333333333333_real64, 0.3333333333333333_real64&
        ], [n_bins, n_points])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 1: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 1: pmf don't match")
        call assert_equal_int(included(1), 6, "test_build_residual_histograms: Test 1: included row 1")
        call assert_equal_int(included(2), 6, "test_build_residual_histograms: Test 1: included row 2")
        call assert_equal_int(included(3), 6, "test_build_residual_histograms: Test 1: included row 3")

        ! ============================================================
        ! Test 2 — NaNs must be ignored
        ! ============================================================
        residuals(:, 1) = 0.0_real64
        residuals(:, 2) = [M_NAN, 0.0_real64, M_NAN]
        residuals(:, 3) = [0.0_real64, M_NAN, 0.0_real64]
        residuals(:, 4) = [0.0_real64, 0.0_real64, M_NAN]

        E = 1
        E(1, 1) = 2
        E(1, 2) = 3
        E(2, 3) = 4

        call build_residual_histograms(E, n_neighbors, n_points, residuals, n_reps, n_neighbors * n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 2: ierr should be OK")

        expected_counts = reshape([&
            0, 0, 4, 0,&
            0, 0, 5, 0,&
            0, 0, 5, 0&
        ], [n_bins, n_points])
        expected_pmf = reshape([&
            0.0_real64,  0.0_real64, 1.0_real64, 0.0_real64,&
            0.0_real64,  0.0_real64, 1.0_real64, 0.0_real64,&
            0.0_real64,  0.0_real64, 1.0_real64, 0.0_real64&
        ], [n_bins, n_points])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 2: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 2: pmf don't match")
        call assert_equal_int(included(1), 4, "test_build_residual_histograms: Test 2: included row 1")
        call assert_equal_int(included(2), 5, "test_build_residual_histograms: Test 2: included row 2")
        call assert_equal_int(included(3), 5, "test_build_residual_histograms: Test 2: included row 3")

        ! ============================================================
        ! Test 3 — All NaN → pmf = 0, counts = 0, included = 0
        ! ============================================================
        residuals = M_NAN
        E = 1

        call build_residual_histograms(E, n_neighbors, n_points, residuals, n_reps, n_neighbors * n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 3: ierr should be OK")

        expected_counts = 0.0_real64
        expected_pmf = 0.0_real64
        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 3: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 3: pmf don't match")
        call assert_equal_int(included(1), 0, "test_build_residual_histograms: Test 2: included row 1")
        call assert_equal_int(included(2), 0, "test_build_residual_histograms: Test 2: included row 2")
        call assert_equal_int(included(3), 0, "test_build_residual_histograms: Test 2: included row 3")

        ! ============================================================
        ! Test 4 — Residuals exactly on boundaries
        ! ============================================================
        !
        ! R = 2, bins = 4, width = 1
        ! Values: -2, -1, 0, 1, 2
        !
        residuals(:, 1) = [-2.0, -1.0, 0.0]
        residuals(:, 2) = [1.0, 2.0, 0.0]

        E(1, :) = 1
        E(2, :) = 2

        call build_residual_histograms(E, n_neighbors, n_points, residuals, n_reps, n_neighbors * n_points, R, n_bins, &
                                       counts, pmf, included, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_build_residual_histograms: Test 4: ierr should be OK")

        ! Expected binning:
        ! -2 → bin 1
        ! -1 → bin 2
        !  0 → bin 3
        !  1 → bin 4
        !  2 → bin 4 (right boundary included)
        !
        expected_counts = reshape([&
            1, 1, 2, 2,&
            1, 1, 2, 2,&
            1, 1, 2, 2&
        ], [n_bins, n_points])
        expected_pmf = reshape([&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.3333333333333333_real64, 0.3333333333333333_real64,&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.3333333333333333_real64, 0.3333333333333333_real64,&
            0.16666666666666666_real64, 0.16666666666666666_real64, 0.3333333333333333_real64, 0.3333333333333333_real64&
        ], [n_bins, n_points])

        call assert_equal_array_int(counts, expected_counts, size(counts, kind=int32), "test_build_residual_histograms: Test 3: counts don't match")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, "test_build_residual_histograms: Test 3: pmf don't match")
        call assert_equal_int(included(1), 6, "test_build_residual_histograms: Test 3: included row 1")
        call assert_equal_int(included(2), 6, "test_build_residual_histograms: Test 3: included row 2")
        call assert_equal_int(included(3), 6, "test_build_residual_histograms: Test 3: included row 3")

    end subroutine test_build_residual_histograms

    subroutine test_compute_divergence_per_reference_point
        integer(int32), parameter :: n_points = 3, n_bins = 4
        real(real64), dimension(n_points, n_bins) :: p, q
        real(real64), dimension(n_points) :: jsd, expected_jsd
        integer(int32) :: ierr
        real(real64) :: tol

        ! ============================================================
        ! Test 1 — Identical PMFs → JSD = 0
        ! ============================================================
        p = reshape([0.1, 0.2, 0.3, 0.4, &
                     0.25,0.25,0.25,0.25, &
                     1.0, 0.0, 0.0, 0.0], [n_points,n_bins])
        q = p

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 1: ierr OK")

        expected_jsd = 0.0_real64

        expected_jsd = expected_jsd / log(2.0)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 1: identical PMFs → JSD=0")

        ! ============================================================
        ! Test 2 — Completely disjoint PMFs → JSD = log(2)
        ! ============================================================
        !
        ! p = [1,0,0,0]
        ! q = [0,1,0,0]
        !
        ! mix = [0.5,0.5,0,0]
        !
        ! KL(p||mix) = 1 * log(1/0.5) = log(2)
        ! KL(q||mix) = log(2)
        !
        ! JSD = 0.5*(log2 + log2) = log(2)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,1) = 1.0
        q(1,2) = 1.0

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 2: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(1) = log(2.0_real64)
        expected_jsd = expected_jsd / log(2.0)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 2: disjoint PMFs → JSD=log(2)")
        call assert_equal_real(jsd(2), 0.0_real64, TOL, "test_compute_divergence_per_reference_point: Test 2: rows 2,3 are zero PMFs → JSD=0")
        call assert_equal_real(jsd(3), 0.0_real64, TOL, "test_compute_divergence_per_reference_point: Test 2: rows 2,3 are zero PMFs → JSD=0")

        ! ============================================================
        ! Test 3 — Partially overlapping PMFs (analytic check)
        ! ============================================================
        !
        ! p = [0.5, 0.5, 0, 0]
        ! q = [0.0, 1.0, 0, 0]
        !
        ! mix = [0.25, 0.75, 0, 0]
        !
        ! KL(p||mix) = 0.5*log(0.5/0.25) + 0.5*log(0.5/0.75)
        !            = 0.5*log(2) + 0.5*log(2/3)
        !
        ! KL(q||mix) = 1.0*log(1.0/0.75)
        !
        ! JSD = 0.5*(KL_p + KL_q)
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,:) = [0.5, 0.5, 0.0, 0.0]
        q(1,:) = [0.0, 1.0, 0.0, 0.0]

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 3: ierr OK")

        ! Compute expected value analytically
        expected_jsd = 0.0_real64
        expected_jsd(1) = 0.5_real64 * ( &
            0.5_real64*log(2.0_real64) + 0.5_real64*log(2.0_real64/3.0_real64) &
            + log(1.0_real64/0.75_real64) )

        expected_jsd = expected_jsd / log(2.0)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 3: analytic partial-overlap JSD")

        ! ============================================================
        ! Test 4 — Zero-probability bins handled correctly
        ! ============================================================
        !
        ! p = [1,0,0,0]
        ! q = [0,0,1,0]
        !
        ! mix = [0.5,0,0.5,0]
        !
        ! KL(p||mix) = log(1/0.5) = log(2)
        ! KL(q||mix) = log(1/0.5) = log(2)
        !
        ! JSD = log(2) -> 1
        !
        p = 0.0_real64
        q = 0.0_real64
        p(1,1) = 1.0
        q(1,3) = 1.0

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 4: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(1) = log(2.0_real64)
        expected_jsd = expected_jsd / log(2.0)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 4: zero-probability bins handled correctly")

        ! ============================================================
        ! Test 5 — Multiple neighbors, mixed patterns
        ! ============================================================
        p = reshape([ &
            0.2,0.3,0.5,&
            0.0,1.0,0.0,&
            0.0,0.0,0.25,&
            0.25,0.25,0.25&
        ], [n_points,n_bins])

        q = reshape([ &
            0.2,0.3,0.5,&
            0.0,0.0,1.0,&
            0.0,0.0,0.25,&
            0.25,0.25,0.25&
        ], [n_points,n_bins])

        call compute_divergence_per_reference_point(p, q, n_points, n_bins, jsd, ierr)
        call assert_equal_int(ierr, ERR_OK, "test_compute_divergence_per_reference_point: Test 5: ierr OK")

        expected_jsd = 0.0_real64
        expected_jsd(2) = 0.5 * (1.0_real64 * log(1.0_real64 / 0.5_real64))
        expected_jsd(3) = 0.5 * (1.0_real64 * log(1.0_real64 / 0.5_real64))
        expected_jsd = expected_jsd / log(2.0)
        call assert_equal_array_real(jsd, expected_jsd, size(jsd, kind=int32), TOL, "test_compute_divergence_per_reference_point: Test 5: mixed patterns")
    end subroutine test_compute_divergence_per_reference_point

    subroutine test_compute_weighted_global_divergence
        integer(int32), parameter :: n_points = 4
        real(real64), dimension(n_points) :: jsd
        integer(int32), dimension(n_points) :: n1, n2
        real(real64), dimension(n_points) :: w, expected_weights
        real(real64) :: global_jsd
        integer(int32) :: ierr
        real(real64) :: expected, tol

        ! ============================================================
        ! Test 1 — Simple case: equal sample counts → uniform weights
        ! ============================================================
        !
        ! jsd = [0.1, 0.2, 0.3, 0.4]
        ! n1 = [5,5,5,5]
        ! n2 = [5,5,5,5]
        !
        ! n_j = [10,10,10,10]
        ! T = 40
        ! w = [0.25,0.25,0.25,0.25]
        !
        ! global_jsd = 0.25*(0.1+0.2+0.3+0.4) = 0.25
        !
        jsd = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        n1  = [5_int32,5_int32,5_int32,5_int32]
        n2  = [5_int32,5_int32,5_int32,5_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        expected_weights = 0.25_real64
        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 1: ierr OK")

        ! Strangely this array comparison signals DENORMAL in gfortran if a test fails, seems to be an optimization bug
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 1: uniform weights")
        call assert_equal_real(global_jsd, 0.25_real64, TOL, "test_compute_weighted_global_divergence: Test 1: global JSD")

        ! ============================================================
        ! Test 2 — Unequal sample counts → weighted average
        ! ============================================================
        !
        ! jsd = [1.0, 2.0, 3.0, 4.0]
        ! n1 = [10, 20, 30, 40]
        ! n2 = [ 0, 10, 10, 10]
        !
        ! n_j = [10, 30, 40, 50]
        ! T = 130
        ! w = [10/130, 30/130, 40/130, 50/130]
        !
        ! global_jsd = sum_j w(j)*jsd(j)
        !
        jsd = [0.1_real64, 0.2_real64, 0.3_real64, 0.4_real64]
        n1  = [10_int32,20_int32,30_int32,40_int32]
        n2  = [ 0_int32,10_int32,10_int32,10_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 2: ierr OK")
        
        expected_weights = [1.0_real64, 3.0_real64, 4.0_real64, 5.0_real64] / 13.0_real64
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 2: wweights")

        expected = (0.1_real64*(10.0_real64/130.0_real64) + &
                    0.2_real64*(30.0_real64/130.0_real64) + &
                    0.3_real64*(40.0_real64/130.0_real64) + &
                    0.4_real64*(50.0_real64/130.0_real64))

        call assert_equal_real(global_jsd, expected, TOL, "test_compute_weighted_global_divergence: Test 2: weighted global JSD")

        ! ============================================================
        ! Test 3 — Some neighborhoods have zero samples → weight = 0
        ! ============================================================
        !
        ! jsd = [0.5, 1.0, 2.0, 4.0]
        ! n1 = [0, 10, 0, 5]
        ! n2 = [0,  0, 0, 5]
        !
        ! n_j = [0,10,0,10]
        ! T = 20
        ! w = [0, 10/20, 0, 10/20] = [0,0.5,0,0.5]
        !
        ! global_jsd = 0.5*1.0 + 0.5*0.5 = 0.75
        !
        jsd = [0.5_real64, 1.0_real64, 0.2_real64, 0.5_real64]
        n1  = [0_int32,10_int32,0_int32,5_int32]
        n2  = [0_int32, 0_int32,0_int32,5_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 3: ierr OK")

        expected_weights = [0.0_real64,0.5_real64,0.0_real64,0.5_real64]
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 3: weights with zero-sample neighborhoods")
        call assert_equal_real(global_jsd, 0.75_real64, TOL, "test_compute_weighted_global_divergence: Test 3: weighted global JSD")

        ! ============================================================
        ! Test 4 — All neighborhoods have zero samples → weights=0, global JSD=0
        ! ============================================================
        !
        ! jsd = [1,2,3,4]
        ! n1 = [0,0,0,0]
        ! n2 = [0,0,0,0]
        !
        ! n_j = [0,0,0,0]
        ! T = 0 → weights = 0
        ! global_jsd = 0
        !
        jsd = [0.1_real64,0.2_real64,0.3_real64,0.4_real64]
        n1  = 0
        n2  = 0

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        expected_weights = 0.0_real64
        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 4: ierr OK")
        call assert_no_nan_real(w, size(w, kind=int32), "test_compute_weighted_global_divergence: Test 4: all weights not NaN")
        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 4: all weights zero")
        call assert_equal_real(global_jsd, 0.0_real64, TOL, "test_compute_weighted_global_divergence: Test 4: global JSD zero when no samples")

        ! ============================================================
        ! Test 5 — Mixed jsd values, mixed sample counts
        ! ============================================================
        !
        ! jsd = [0.0, 0.5, 1.0, 2.0]
        ! n1 = [5, 0, 10, 5]
        ! n2 = [5, 5,  0, 5]
        !
        ! n_j = [10,5,10,10]
        ! T = 35
        ! w = [10/35, 5/35, 10/35, 10/35]
        !
        ! global_jsd = sum_j w(j)*jsd(j)
        !
        jsd = [0.0_real64, 0.5_real64, 1.0_real64, 0.2_real64]
        n1  = [5_int32,0_int32,10_int32,5_int32]
        n2  = [5_int32,5_int32, 0_int32,5_int32]

        call compute_weighted_global_divergence(jsd, n_points, n1, n2, &
                                                global_jsd, w, ierr)

        call assert_equal_int(ierr, ERR_OK, "test_compute_weighted_global_divergence: Test 5: ierr OK")

        expected_weights = [10.0_real64, 5.0_real64, 10.0_real64, 10.0_real64] / 35.0_real64

        call assert_equal_array_real(w, expected_weights, size(w, kind=int32), TOL, "test_compute_weighted_global_divergence: Test 5: all weights zero")

        expected = (0.0_real64*(10.0_real64/35.0_real64) + &
                    0.5_real64*( 5.0_real64/35.0_real64) + &
                    1.0_real64*(10.0_real64/35.0_real64) + &
                    0.2_real64*(10.0_real64/35.0_real64))

        call assert_equal_real(global_jsd, expected, TOL, "test_compute_weighted_global_divergence: Test 5: weighted global JSD")

    end subroutine test_compute_weighted_global_divergence

    ! --------------------------------------------------------------------------
    ! Test Cases for compute_gene_means
    ! --------------------------------------------------------------------------

    ! Test case 1: Basic compute_gene_means functionality.
    subroutine test_compute_gene_means_basic()
        integer, parameter :: n_genes = 4, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes)
        real(real64) :: expected_means(n_genes)
        integer(int32) :: ierr
        
        ! Test data
        expr = reshape([1.0, 2.0, 3.0,    &   ! Gene 1: mean = 2.0
                        4.0, 5.0, 6.0,    &   ! Gene 2: mean = 5.0
                        10.0, 20.0, 30.0, &   ! Gene 3: mean = 20.0
                        0.0, 0.0, 0.0],   &   ! Gene 4: mean = 0.0
                       [n_reps, n_genes])
        
        expected_means = [2.0, 5.0, 20.0, 0.0]
        
        call compute_gene_means(expr, n_genes, n_reps, means, n_genes, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_basic: should succeed")
        call assert_equal_array_real(means, expected_means, n_genes, &
                                        TOL, "test_compute_gene_means_basic: means")
    end subroutine test_compute_gene_means_basic

    ! Test case 2: compute_gene_means with NaN values.
    subroutine test_compute_gene_means_with_nan()
        integer, parameter :: n_genes = 3, n_reps = 4
        real(real64) :: expr(n_reps, n_genes), means(n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0_real64, 2.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 3.0_real64]  ! mean = (1+2+3)/3 = 2.0
        expr(:, 2) = [ieee_value(1.0_real64, ieee_quiet_nan), 5.0_real64, 7.0_real64, 9.0_real64]  ! mean = (5+7+9)/3 = 7.0
        expr(:, 3) = [10.0, 20.0, 30.0, 40.0]  ! mean = 25.0
        
        call compute_gene_means(expr, n_genes, n_reps, means, n_genes, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_with_nan: should succeed")
        call assert_equal_real(means(1), 2.0_real64, TOL, "test_compute_gene_means_with_nan: gene 1 mean")
        call assert_equal_real(means(2), 7.0_real64, TOL, "test_compute_gene_means_with_nan: gene 2 mean")
        call assert_equal_real(means(3), 25.0_real64, TOL, "test_compute_gene_means_with_nan: gene 3 mean")
    end subroutine test_compute_gene_means_with_nan

    ! Test case 3: compute_gene_means with all NaN values for a gene.
    subroutine test_compute_gene_means_all_nan()
        integer, parameter :: n_genes = 2, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0, 2.0, 3.0]  ! Normal gene
        expr(:, 2) = ieee_value(1.0_real64, ieee_quiet_nan)  ! All NaN gene
        
        call compute_gene_means(expr, n_genes, n_reps, means, n_genes, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_gene_means_all_nan: should succeed")
        call assert_equal_real(means(1), 2.0_real64, TOL, "test_compute_gene_means_all_nan: gene 1 mean")
        call assert_true(ieee_is_nan(means(2)), "test_compute_gene_means_all_nan: gene 2 should be NaN")
    end subroutine test_compute_gene_means_all_nan

    ! Test case 4: compute_gene_means with invalid input.
    subroutine test_compute_gene_means_invalid_input()
        integer, parameter :: n_genes = 5, n_reps = 3, n_genes_neg = -1
        real(real64) :: expr(3, 1), means(1)
        integer(int32) :: ierr
        
        ! Test with zero genes
        call compute_gene_means(expr, 0_int32, n_reps, means, n_genes, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 2_int32), "test_compute_gene_means_invalid_input: zero genes should fail")

        call compute_gene_means(expr, n_genes, n_reps, means, 0_int32, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 5_int32), "test_compute_gene_means_invalid_input: zero max_n_genes should fail")
        
        ! Test with negative genes
        call compute_gene_means(expr, n_genes_neg, n_reps, means, n_genes_neg, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_INVALID_INPUT, 2_int32), "test_compute_gene_means_invalid_input: negative genes should fail")

        ! Test with zero replicates
        call compute_gene_means(expr, n_genes, 0_int32, means, n_genes, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 3_int32), "test_compute_gene_means_invalid_input: zero replicates should fail")
    end subroutine test_compute_gene_means_invalid_input

    ! --------------------------------------------------------------------------
    ! Test Cases for compute_residuals
    ! --------------------------------------------------------------------------

    ! Test case 5: Basic compute_residuals functionality.
    subroutine test_compute_residuals_basic()
        integer, parameter :: n_genes = 4, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
        real(real64) :: expected_resid(n_reps, n_genes)
        integer(int32) :: ierr
        
        expr = reshape([1.0, 2.0, 3.0,    &   ! Gene 1
                        4.0, 5.0, 6.0,    &   ! Gene 2
                        10.0, 20.0, 30.0, &   ! Gene 3
                        0.0, 0.0, 0.0],   &   ! Gene 4
                       [n_reps, n_genes])
        
        means = [2.0, 5.0, 20.0, 0.0]
        expected_resid = reshape([-1.0, 0.0, 1.0,     &   ! Gene 1 residuals
                                    -1.0, 0.0, 1.0,     &   ! Gene 2 residuals
                                    -10.0, 0.0, 10.0,   &   ! Gene 3 residuals
                                    0.0, 0.0, 0.0],     &   ! Gene 4 residuals
                                 [n_reps, n_genes])
        
        call compute_residuals(expr, n_genes, n_reps, means, n_genes, n_reps, resid, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_basic: should succeed")
        call assert_equal_array_real(resid, &
                                        expected_resid, &
                                        n_reps*n_genes, TOL, &
                                        "test_compute_residuals_basic: residuals")
    end subroutine test_compute_residuals_basic

    ! Test case 6: compute_residuals with NaN values.
    subroutine test_compute_residuals_with_nan()
        integer, parameter :: n_genes = 2, n_reps = 4
        real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0_real64, 2.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 3.0_real64]
        expr(:, 2) = [ieee_value(1.0_real64, ieee_quiet_nan), 5.0_real64, 7.0_real64, 9.0_real64]
        means = [2.0_real64, 7.0_real64]
        
        call compute_residuals(expr, n_genes, n_reps, means, n_genes, n_reps, resid, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_with_nan: should succeed")
        ! Check specific values
        call assert_equal_real(resid(1, 1), -1.0_real64, TOL, "test_compute_residuals_with_nan: resid(1,1)")
        call assert_equal_real(resid(2, 1), 0.0_real64, TOL, "test_compute_residuals_with_nan: resid(2,1)")
        call assert_true(ieee_is_nan(resid(3, 1)), "test_compute_residuals_with_nan: resid(3,1) should be NaN")
        call assert_equal_real(resid(4, 1), 1.0_real64, TOL, "test_compute_residuals_with_nan: resid(4,1)")
        
        call assert_true(ieee_is_nan(resid(1, 2)), "test_compute_residuals_with_nan: resid(1,2) should be NaN")
        call assert_equal_real(resid(2, 2), -2.0_real64, TOL, "test_compute_residuals_with_nan: resid(2,2)")
        call assert_equal_real(resid(3, 2), 0.0_real64, TOL, "test_compute_residuals_with_nan: resid(3,2)")
        call assert_equal_real(resid(4, 2), 2.0_real64, TOL, "test_compute_residuals_with_nan: resid(4,2)")
    end subroutine test_compute_residuals_with_nan

    ! Test case 7: compute_residuals with all NaN values.
    subroutine test_compute_residuals_all_nan()
        integer, parameter :: n_genes = 2, n_reps = 3
        real(real64) :: expr(n_reps, n_genes), means(n_genes), resid(n_reps, n_genes)
        integer(int32) :: ierr
        
        expr(:, 1) = [1.0, 2.0, 3.0]
        expr(:, 2) = ieee_value(1.0_real64, ieee_quiet_nan)  ! All NaN
        means = [2.0, 0.0]  ! Second mean is irrelevant
        
        call compute_residuals(expr, n_genes, n_reps, means, n_genes, n_reps, resid, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_compute_residuals_all_nan: should succeed")
        ! All residuals for gene 2 should be NaN
        call assert_true(all(ieee_is_nan(resid(:, 2))), "test_compute_residuals_all_nan: all residuals for NaN gene should be NaN")
    end subroutine test_compute_residuals_all_nan

    ! Test case 8: compute_residuals with invalid input.
    subroutine test_compute_residuals_invalid_input()
        integer, parameter :: n_genes = 0, n_reps = 3
        real(real64) :: expr(3, 1), means(1), resid(3, 1)
        integer(int32) :: ierr
        
        ! Test with zero genes
        call compute_residuals(expr, n_genes, n_reps, means, n_genes, n_reps, resid, ierr)
        call assert_not_equal_int(ierr, ERR_OK, "test_compute_residuals_invalid_input: zero genes should fail")

    end subroutine test_compute_residuals_invalid_input

    ! --------------------------------------------------------------------------
    ! Test Cases for pool_means_alloc
    ! --------------------------------------------------------------------------

    ! Test case 9: Basic pool_means_alloc functionality.
    subroutine test_pool_means_alloc_basic()
        integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 5, n_points = 3, max_n_genes_all_studies = max(n_genes_S1, n_genes_S2), n_studies = 2
        real(real64), target :: means(max_n_genes_all_studies, n_studies), x_star(n_points)
        integer(int32) :: N_pool, ierr
        real(real64), dimension(:), pointer :: means_flat

        means_flat(1:size(means, kind=int32)) => means
        
        means(:, 1) = [1.0, 3.0, 5.0, 7.0, 9.0]
        means(:, 2) = [2.0, 4.0, 6.0, 8.0, 10.0]
        
        call pool_means_alloc(means_flat, n_studies, max_n_genes_all_studies, n_points, N_pool, x_star, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_basic: should succeed")
        call assert_equal_int(N_pool, 10, "test_pool_means_alloc_basic: N_pool should be 10")
        
        ! Check that x_star contains quantiles from pooled data
        ! Pooled data: [1,2,3,4,5,6,7,8,9,10]
        ! For n_points=3, quantiles at positions: 10/4=2.5, 20/4=5.0, 30/4=7.5
        ! Floored: 2, 5, 7 -> values: 2, 5, 7 -> interpolation to 3.25, 5.5 and 7.75
        call assert_equal_real(x_star(1), 3.25_real64, TOL, "test_pool_means_alloc_basic: first quantile")
        call assert_equal_real(x_star(2), 5.5_real64, TOL, "test_pool_means_alloc_basic: second quantile")
        call assert_equal_real(x_star(3), 7.75_real64, TOL, "test_pool_means_alloc_basic: third quantile")
    end subroutine test_pool_means_alloc_basic

    ! Test case 10: pool_means_alloc with NaN values.
    subroutine test_pool_means_alloc_with_nan()
        integer, parameter :: n_genes_S1 = 4, n_genes_S2 = 4, n_points = 2, max_n_genes_all_studies = max(n_genes_S1, n_genes_S2), n_studies = 2
        real(real64), target :: means(max_n_genes_all_studies, n_studies), x_star(n_points)
        integer(int32) :: N_pool, ierr
        real(real64), dimension(:), pointer :: means_flat

        means_flat(1:size(means, kind=int32)) => means
        
        means(:, 1) = [1.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 3.0_real64, 5.0_real64]
        means(:, 2) = [2.0_real64, 4.0_real64, ieee_value(1.0_real64, ieee_quiet_nan), 6.0_real64]
        
        call pool_means_alloc(means_flat, n_studies, max_n_genes_all_studies, n_points, N_pool, x_star, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_with_nan: should succeed")
        call assert_equal_int(N_pool, 6, "test_pool_means_alloc_with_nan: N_pool should exclude NaN values")
        
        ! Pooled data (excluding NaN): [1,2,3,4,5,6]
        ! Values: 2.666, 4.3333 -> interpolation
        call assert_equal_real(x_star(1), 2.0_real64 + 2.0_real64/3.0_real64, TOL, "test_pool_means_alloc_with_nan: first quantile")
        call assert_equal_real(x_star(2), 4.0_real64 + 1.0_real64/3.0_real64, TOL, "test_pool_means_alloc_with_nan: second quantile")
    end subroutine test_pool_means_alloc_with_nan

    ! Test case 11: pool_means_alloc with single study.
    subroutine test_pool_means_alloc_single_study()
        integer, parameter :: n_genes_S1 = 5, n_points = 3, max_n_genes_all_studies = n_genes_S1, n_studies = 1
        real(real64), target :: means(max_n_genes_all_studies, n_studies), x_star(n_points)
        integer(int32) :: N_pool, ierr
        real(real64), dimension(:), pointer :: means_flat

        means_flat(1:size(means, kind=int32)) => means

        means(:, 1) = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        call pool_means_alloc(means_flat, n_studies, max_n_genes_all_studies, n_points, N_pool, x_star, ierr)
        
        call assert_equal_int(ierr, ERR_OK, "test_pool_means_alloc_single_study: should succeed")
    end subroutine test_pool_means_alloc_single_study

    ! Test case 12: pool_means_alloc with invalid input.
    subroutine test_pool_means_alloc_invalid_input()
        integer, parameter :: n_genes_S1 = 5, n_genes_S2 = 5, n_points = 3, max_n_genes_all_studies = max(n_genes_S1, n_genes_S2), n_studies = 2
        real(real64), target :: means(max_n_genes_all_studies, n_studies), x_star(n_points)
        integer(int32) :: N_pool, ierr
        real(real64), dimension(:), pointer :: means_flat

        means_flat(1:size(means, kind=int32)) => means
        
        means(:, 1) = [1.0, 2.0, 3.0, 4.0, 5.0]
        means(:, 2) = [1.0, 2.0, 3.0, 4.0, 5.0]
        
        ! Test with zero genes in S1
        call pool_means_alloc(means_flat, 0_int32, max_n_genes_all_studies, n_points, N_pool, x_star, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 2_int32), "test_pool_means_alloc_invalid_input: zero studies should fail")
        
        ! Test with zero points
        call pool_means_alloc(means_flat, n_studies, max_n_genes_all_studies, 0_int32, N_pool, x_star, ierr)
        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 4_int32), "test_pool_means_alloc_invalid_input: zero points should fail")
    end subroutine test_pool_means_alloc_invalid_input

    subroutine test_construct_neighborhoods_basic()
        integer(int32), parameter :: n_points    = 3
        integer(int32), parameter :: n_genes_S   = 5
        integer(int32), parameter :: n_neighbors = 2

        integer(int32) :: ierr
        real(real64) :: x_star(n_points)
        real(real64) :: mean_S(n_genes_S)
        integer(int32) :: neighborhood_residuals(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)

        ! -----------------------------
        ! Inputs
        ! -----------------------------
        x_star = [ 0.0_real64, 2.0_real64, 10.0_real64 ]
        mean_S = [ 1.0, 2.5, 2.5, 10.5, 20.0 ]

        ! -----------------------------
        ! Call routine
        ! -----------------------------
        call construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)

        call assert_equal_int(ierr, ERR_OK, "ierr must be ERR_OK")


        ! third gene sorts before second
        call assert_equal_array_int( neighborhood_residuals(:,1), [1,3], n_neighbors, &
            "test_construct_neighborhoods_basic: Incorrect neighborhood residuals for point 1" )
        call assert_equal_array_int( neighborhood_range(:,1), [1,3], n_neighbors, &
            "test_construct_neighborhoods_basic: Incorrect neighborhood range for point 1" )

        ! third gene sorts before second
        call assert_equal_array_int( neighborhood_residuals(:,2), [3,2], n_neighbors, &
            "test_construct_neighborhoods_basic: Incorrect neighborhood residuals for point 2" )
        call assert_equal_array_int( neighborhood_range(:,2), [2,3], n_neighbors, &
            "test_construct_neighborhoods_basic: Incorrect neighborhood range for point 2" )

        ! third gene sorts before second
        call assert_equal_array_int( neighborhood_residuals(:,3), [4,2], n_neighbors, &
            "test_construct_neighborhoods_basic: Incorrect neighborhood residuals for point 3" )
        call assert_equal_array_int( neighborhood_range(:,3), [2,4], n_neighbors, &
            "test_construct_neighborhoods_basic: Incorrect neighborhood range for point 3" )

        call construct_neighborhoods_alloc(0_int32, x_star, n_genes_S, mean_S, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)

        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 1_int32), "test_construct_neighborhoods_basic: for n_points=0 ierr must be ERR_EMPTY_INPUT")

        call construct_neighborhoods_alloc(n_points, x_star, 0_int32, mean_S, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)

        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 3_int32), "test_construct_neighborhoods_basic: for n_genes_S=0 ierr must be ERR_EMPTY_INPUT")

        call construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, neighborhood_residuals, neighborhood_range, 0_int32, ierr)

        call assert_equal_int(ierr, create_err_code(ERR_EMPTY_INPUT, 7_int32), "test_construct_neighborhoods_basic: for n_neighbors=0 ierr must be ERR_EMPTY_INPUT")
    end subroutine test_construct_neighborhoods_basic

    subroutine test_construct_neighborhoods_nan_handling()
        integer(int32), parameter :: n_points    = 2
        integer(int32), parameter :: n_genes_S   = 4
        integer(int32), parameter :: n_neighbors = 2

        integer(int32) :: ierr
        real(real64) :: x_star(n_points)
        real(real64) :: mean_S(n_genes_S)
        integer(int32) :: neighborhood_residuals(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)

        x_star = [ 5.0_real64, M_NAN ]
        mean_S = [ 4.0_real64, M_NAN, 6.0_real64, M_NAN ]

        call construct_neighborhoods_alloc(n_points, x_star, n_genes_S, mean_S, neighborhood_residuals, neighborhood_range, n_neighbors, ierr)

        call assert_equal_int(ierr, ERR_OK, "ierr must be ERR_OK")

        ! Only genes 1 and 3 are valid (non-NaN)
        call assert_equal_array_int( neighborhood_residuals(:,1), [1,3], n_neighbors, &
            "test_construct_neighborhoods_nan_handling: NaN mean handling incorrect" )
        call assert_equal_array_int( neighborhood_range(:,1), [1,2], n_neighbors, &
            "test_construct_neighborhoods_nan_handling: NaN mean handling Incorrect neighborhood range for point 1")

        call assert_equal_array_int( neighborhood_residuals(:,2), [1,3], n_neighbors, &
            "test_construct_neighborhoods_nan_handling: NaN x_star handling incorrect" )
        call assert_equal_array_int( neighborhood_range(:,2), [1,2], n_neighbors, &
            "test_construct_neighborhoods_nan_handling: NaN x star handling Incorrect neighborhood range for point 2")
    end subroutine test_construct_neighborhoods_nan_handling

end module mod_test_data_integration
