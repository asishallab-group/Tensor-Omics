!> Unit test suite for the gJCT/fJCT parameter-search building blocks added for the JSD-Comp-Test
!| port (Issue #126): the ranged kNN neighborhood construction, the N-study shared residual
!| range, the counts-to-pmf step factored out of the histogram builder, the GAMMA-decay
!| candidate-grid generator and its bin-count estimator, the two candidate admissibility gates,
!| the plateau check that decides when the parameter search has converged, the consensus-pmf
!| builder, and the bootstrap confidence interval (with its recommended top-k/bottom-k heap
!| size). Later stages of the port keep adding to this suite.
module mod_test_data_integration_js_comp_test
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_data_integration
    use tox_data_integration_js_comp_test, only: estimate_bin_count, generate_js_comp_test_candidates, &
                                                  generate_js_comp_test_candidates_expert, check_neighborhood_overlaps, &
                                                  check_mean_pmf_min_counts, check_plateau_condition, &
                                                  check_effect_size_plateau_condition, create_mean_pmf, &
                                                  create_mean_pmf_only, bootstrap_histogram, run_js_comp_test, &
                                                  run_js_comp_test_parameter_search
    use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, &
                                                       MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, &
                                                       MODE_PLATEAU_BOTH, calc_js_comp_test_n_top_k_jsds
    use tox_errors
    use f42_math_impl, only: above, below
    use test_suite, only: test_case

    implicit none

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests_data_integration_js_comp_test() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(54))

        all_tests(1) = test_case("test_construct_neighborhoods_ranged_basic", test_construct_neighborhoods_ranged_basic)
        all_tests(2) = test_case("test_construct_neighborhoods_ranged_tie_extends_range", &
                                 test_construct_neighborhoods_ranged_tie_extends_range)
        all_tests(3) = test_case("test_construct_neighborhoods_ranged_all_nan_means", &
                                 test_construct_neighborhoods_ranged_all_nan_means)
        all_tests(4) = test_case("test_construct_neighborhoods_ranged_fewer_genes_than_neighbors", &
                                 test_construct_neighborhoods_ranged_fewer_genes_than_neighbors)
        all_tests(5) = test_case("test_construct_neighborhoods_ranged_plain_and_validation", &
                                 test_construct_neighborhoods_ranged_plain_and_validation)

        all_tests(6) = test_case("test_all_studies_range_matches_two_study_routine", &
                                 test_all_studies_range_matches_two_study_routine)
        all_tests(7) = test_case("test_determine_all_studies_shared_residual_range_three_studies", &
                                 test_determine_all_studies_shared_residual_range_three_studies)
        all_tests(8) = test_case("test_all_studies_range_invalid_quantile", &
                                 test_all_studies_range_invalid_quantile)

        all_tests(9) = test_case("test_calc_pmf_matches_build_residual_histograms", &
                                 test_calc_pmf_matches_build_residual_histograms)
        all_tests(10) = test_case("test_calc_pmf_zero_included_reps", test_calc_pmf_zero_included_reps)
        all_tests(11) = test_case("test_calc_pmf_rejects_negative_counts", test_calc_pmf_rejects_negative_counts)

        all_tests(12) = test_case("test_estimate_bin_count_basic_hand_computed", test_estimate_bin_count_basic_hand_computed)
        all_tests(13) = test_case("test_estimate_bin_count_all_nan_gives_one_bin", &
                                  test_estimate_bin_count_all_nan_gives_one_bin)
        all_tests(14) = test_case("test_estimate_bin_count_clamped_to_max_n_bins", &
                                  test_estimate_bin_count_clamped_to_max_n_bins)

        all_tests(15) = test_case("test_generate_js_comp_test_candidates_collapses_at_8742", &
                                  test_generate_js_comp_test_candidates_collapses_at_8742)
        all_tests(16) = test_case("test_generate_js_comp_test_candidates_has_two_distinct_at_8743", &
                                  test_generate_js_comp_test_candidates_has_two_distinct_at_8743)
        all_tests(17) = test_case("test_generate_js_comp_test_candidates_validation", &
                                  test_generate_js_comp_test_candidates_validation)

        all_tests(18) = test_case("test_neighborhood_overlaps_all_pass", test_neighborhood_overlaps_all_pass)
        all_tests(19) = test_case("test_neighborhood_overlaps_exactly_at_threshold_passes", &
                                  test_neighborhood_overlaps_exactly_at_threshold_passes)
        all_tests(20) = test_case("test_neighborhood_overlaps_below_threshold_fails", &
                                  test_neighborhood_overlaps_below_threshold_fails)

        all_tests(21) = test_case("test_mean_pmf_min_counts_all_pass", test_mean_pmf_min_counts_all_pass)
        all_tests(22) = test_case("test_mean_pmf_min_counts_exactly_at_threshold_passes", &
                                  test_mean_pmf_min_counts_exactly_at_threshold_passes)
        all_tests(23) = test_case("test_mean_pmf_min_counts_below_threshold_fails", &
                                  test_mean_pmf_min_counts_below_threshold_fails)

        all_tests(24) = test_case("test_check_plateau_condition_join_min_requires_all", &
                                  test_check_plateau_condition_join_min_requires_all)
        all_tests(25) = test_case("test_check_plateau_condition_join_max_requires_any", &
                                  test_check_plateau_condition_join_max_requires_any)
        all_tests(26) = test_case("test_check_plateau_condition_join_median_requires_majority", &
                                  test_check_plateau_condition_join_median_requires_majority)
        all_tests(27) = test_case("test_check_plateau_condition_worse_than_previous_short_circuits", &
                                  test_check_plateau_condition_worse_than_previous_short_circuits)

        all_tests(28) = test_case("test_create_mean_pmf_two_studies_hand_computed", &
                                  test_create_mean_pmf_two_studies_hand_computed)
        all_tests(29) = test_case("test_create_mean_pmf_three_studies_includes_self", &
                                  test_create_mean_pmf_three_studies_includes_self)
        all_tests(30) = test_case("test_create_mean_pmf_only_matches_create_mean_pmf", &
                                  test_create_mean_pmf_only_matches_create_mean_pmf)
        all_tests(31) = test_case("test_create_mean_pmf_validation", test_create_mean_pmf_validation)

        all_tests(32) = test_case("test_calc_js_comp_test_n_top_k_jsds_default_hand_computed", &
                                  test_calc_js_comp_test_n_top_k_jsds_default_hand_computed)
        all_tests(33) = test_case("test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level", &
                                  test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level)
        all_tests(34) = test_case("test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one", &
                                  test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one)

        all_tests(35) = test_case("test_bootstrap_histogram_seeded_reproducibility", &
                                  test_bootstrap_histogram_seeded_reproducibility)
        all_tests(36) = test_case("test_bootstrap_histogram_degenerate_single_bin_zero_ci", &
                                  test_bootstrap_histogram_degenerate_single_bin_zero_ci)

        all_tests(37) = test_case("test_gjct_permutation_test_conservation_of_counts", &
                                  test_gjct_permutation_test_conservation_of_counts)
        all_tests(38) = test_case("test_gjct_permutation_test_seeded_reproducibility", &
                                  test_gjct_permutation_test_seeded_reproducibility)
        all_tests(39) = test_case("test_permutation_pvalue_can_be_exactly_zero_known_limitation", &
                                  test_permutation_pvalue_can_be_exactly_zero_known_limitation)

        all_tests(40) = test_case("test_run_js_comp_test_two_studies_hand_traceable", &
                                  test_run_js_comp_test_two_studies_hand_traceable)
        all_tests(41) = test_case("test_run_js_comp_test_three_studies_outlier_has_small_p_value", &
                                  test_run_js_comp_test_three_studies_outlier_has_small_p_value)
        all_tests(42) = test_case("test_param_search_no_plateau_falls_back_to_finest", &
                                  test_param_search_no_plateau_falls_back_to_finest)
        all_tests(43) = test_case("test_param_search_single_candidate_bypasses_plateau", &
                                  test_param_search_single_candidate_bypasses_plateau)
        all_tests(44) = test_case("test_param_search_finds_plateau_mid_grid", &
                                  test_param_search_finds_plateau_mid_grid)

        all_tests(45) = test_case("test_effect_size_plateau_first_candidate_no_delta", &
                                  test_effect_size_plateau_first_candidate_no_delta)
        all_tests(46) = test_case("test_effect_size_plateau_single_transition_insufficient", &
                                  test_effect_size_plateau_single_transition_insufficient)
        all_tests(47) = test_case("test_effect_size_plateau_two_consecutive_transitions", &
                                  test_effect_size_plateau_two_consecutive_transitions)
        all_tests(48) = test_case("test_effect_size_plateau_resets_on_non_qualifying", &
                                  test_effect_size_plateau_resets_on_non_qualifying)
        all_tests(49) = test_case("test_effect_size_plateau_median_max_hand_computed", &
                                  test_effect_size_plateau_median_max_hand_computed)
        all_tests(50) = test_case("test_param_search_effect_size_mode_plateau", &
                                  test_param_search_effect_size_mode_plateau)
        all_tests(51) = test_case("test_param_search_both_mode_uses_earlier_trigger", &
                                  test_param_search_both_mode_uses_earlier_trigger)
        all_tests(52) = test_case("test_param_search_no_plateau_uses_smallest_uncertainty", &
                                  test_param_search_no_plateau_uses_smallest_uncertainty)
        all_tests(53) = test_case("test_param_search_no_plateau_effect_size_falls_back", &
                                  test_param_search_no_plateau_effect_size_falls_back)
        all_tests(54) = test_case("test_param_search_no_plateau_both_falls_back", &
                                  test_param_search_no_plateau_both_falls_back)
    end function get_all_tests_data_integration_js_comp_test

    !> Basic two-reference-point case, computed by hand from a sorted `mean_S`; cross-checked
    !| against `construct_neighborhoods_expert`'s equivalent distance-sort fixture (both find the
    !| same nearest genes, since with no ties the two algorithms agree).
    subroutine test_construct_neighborhoods_ranged_basic()
        integer(int32), parameter :: n_points = 2, n_genes_S = 5, n_neighbors = 2
        real(real64) :: x_star(n_points), mean_S(n_genes_S)
        integer(int32) :: mean_S_perm(n_genes_S)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)
        integer(int32) :: ierr

        x_star = [2.0_real64, 10.0_real64]
        mean_S = [1.0_real64, 2.5_real64, 9.0_real64, 10.5_real64, 20.0_real64]
        mean_S_perm = [1, 2, 3, 4, 5] ! already ascending, no ties

        call construct_neighborhoods_ranged_expert(n_points, x_star, n_genes_S, mean_S, mean_S_perm, n_neighbors, &
                                                    neighborhood_indices, neighborhood_range, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_construct_neighborhoods_ranged_basic: ierr should be OK")

        ! x_star(1)=2.0: distances to [1,2.5,9,10.5,20] are [1,0.5,7,8.5,18] -> nearest are gene 2, then gene 1
        call assert_equal_array_int(neighborhood_indices(:, 1), [2, 1], n_neighbors, &
                                    "test_construct_neighborhoods_ranged_basic: indices for point 1")
        call assert_equal_array_int(neighborhood_range(:, 1), [1, 2], 2_int32, &
                                    "test_construct_neighborhoods_ranged_basic: range for point 1")

        ! x_star(2)=10.0: distances to [1,2.5,9,10.5,20] are [9,7.5,1,0.5,10] -> nearest are gene 4, then gene 3
        call assert_equal_array_int(neighborhood_indices(:, 2), [4, 3], n_neighbors, &
                                    "test_construct_neighborhoods_ranged_basic: indices for point 2")
        call assert_equal_array_int(neighborhood_range(:, 2), [3, 4], 2_int32, &
                                    "test_construct_neighborhoods_ranged_basic: range for point 2")
    end subroutine test_construct_neighborhoods_ranged_basic

    !> Three genes tie on the same mean at the boundary of a single-neighbor neighborhood: only
    !| one of them is picked as the neighbor, but `neighborhood_range` must still be extended to
    !| span all three, because a later admissibility gate reasons about the range, not just the
    !| picked neighbor.
    subroutine test_construct_neighborhoods_ranged_tie_extends_range()
        integer(int32), parameter :: n_points = 1, n_genes_S = 5, n_neighbors = 1
        real(real64) :: x_star(n_points), mean_S(n_genes_S)
        integer(int32) :: mean_S_perm(n_genes_S)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)
        integer(int32) :: ierr

        x_star = [3.0_real64]
        mean_S = [1.0_real64, 3.0_real64, 3.0_real64, 3.0_real64, 5.0_real64]
        mean_S_perm = [1, 2, 3, 4, 5] ! already ascending

        call construct_neighborhoods_ranged_expert(n_points, x_star, n_genes_S, mean_S, mean_S_perm, n_neighbors, &
                                                    neighborhood_indices, neighborhood_range, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_construct_neighborhoods_ranged_tie_extends_range: ierr should be OK")

        ! Binary search lands on the first "3.0" (gene 2); it is closer than gene 1 (distance 0 vs 2), so
        ! gene 2 is the sole neighbor -- but genes 3 and 4 share its mean and must extend the range.
        call assert_equal_array_int(neighborhood_indices(:, 1), [2], n_neighbors, &
                                    "test_construct_neighborhoods_ranged_tie_extends_range: picked neighbor")
        call assert_equal_array_int(neighborhood_range(:, 1), [2, 4], 2_int32, &
                                    "test_construct_neighborhoods_ranged_tie_extends_range: range must include all ties")
    end subroutine test_construct_neighborhoods_ranged_tie_extends_range

    !> When every gene's mean is NaN, the reference point's own value is irrelevant: the routine
    !| falls back to the first `min(n_genes_S, n_neighbors)` positions, and NaN never compares
    !| equal to itself, so the tie-extension step must not extend the range past that fallback.
    subroutine test_construct_neighborhoods_ranged_all_nan_means()
        integer(int32), parameter :: n_points = 1, n_genes_S = 3, n_neighbors = 2
        real(real64) :: x_star(n_points), mean_S(n_genes_S)
        integer(int32) :: mean_S_perm(n_genes_S)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)
        integer(int32) :: ierr
        real(real64) :: nan_val

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)
        x_star = [5.0_real64]
        mean_S = [nan_val, nan_val, nan_val]
        mean_S_perm = [1, 2, 3]

        call construct_neighborhoods_ranged_expert(n_points, x_star, n_genes_S, mean_S, mean_S_perm, n_neighbors, &
                                                    neighborhood_indices, neighborhood_range, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_construct_neighborhoods_ranged_all_nan_means: ierr should be OK")
        call assert_equal_array_int(neighborhood_indices(:, 1), [1, 2], n_neighbors, &
                                    "test_construct_neighborhoods_ranged_all_nan_means: fallback indices")
        call assert_equal_array_int(neighborhood_range(:, 1), [1, 2], 2_int32, &
                                    "test_construct_neighborhoods_ranged_all_nan_means: fallback range")
    end subroutine test_construct_neighborhoods_ranged_all_nan_means

    !> `n_neighbors` exceeding the number of genes is a documented, deliberately-preserved
    !| limitation: the leftover slots are filled with out-of-range sentinel indices
    !| (`n_genes_S+1, n_genes_S+2, ...`) rather than left undefined.
    subroutine test_construct_neighborhoods_ranged_fewer_genes_than_neighbors()
        integer(int32), parameter :: n_points = 1, n_genes_S = 2, n_neighbors = 3
        real(real64) :: x_star(n_points), mean_S(n_genes_S)
        integer(int32) :: mean_S_perm(n_genes_S)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)
        integer(int32) :: ierr

        x_star = [2.0_real64]
        mean_S = [1.0_real64, 2.0_real64]
        mean_S_perm = [1, 2]

        call construct_neighborhoods_ranged_expert(n_points, x_star, n_genes_S, mean_S, mean_S_perm, n_neighbors, &
                                                    neighborhood_indices, neighborhood_range, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_construct_neighborhoods_ranged_fewer_genes_than_neighbors: ierr should be OK")
        ! gene 2 (distance 0), then gene 1 (distance 1), then both pointers are exhausted -> sentinel 3
        call assert_equal_array_int(neighborhood_indices(:, 1), [2, 1, 3], n_neighbors, &
                                    "test_construct_neighborhoods_ranged_fewer_genes_than_neighbors: indices")
        call assert_equal_array_int(neighborhood_range(:, 1), [1, 2], 2_int32, &
                                    "test_construct_neighborhoods_ranged_fewer_genes_than_neighbors: range")
    end subroutine test_construct_neighborhoods_ranged_fewer_genes_than_neighbors

    !> The plain wrapper builds and sorts `mean_S_perm` itself, and the generated validation
    !| rejects a non-positive `n_neighbors` and empty extents.
    subroutine test_construct_neighborhoods_ranged_plain_and_validation()
        integer(int32), parameter :: n_points = 2, n_genes_S = 5, n_neighbors = 2
        real(real64) :: x_star(n_points), mean_S(n_genes_S)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points)
        integer(int32) :: neighborhood_range(2, n_points)
        integer(int32) :: ierr

        x_star = [2.0_real64, 10.0_real64]
        mean_S = [1.0_real64, 2.5_real64, 9.0_real64, 10.5_real64, 20.0_real64]

        call construct_neighborhoods_ranged(n_points, x_star, n_genes_S, mean_S, n_neighbors, &
                                            neighborhood_indices, neighborhood_range, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_construct_neighborhoods_ranged_plain_and_validation: ierr should be OK")
        call assert_equal_array_int(neighborhood_indices(:, 1), [2, 1], n_neighbors, &
                                    "test_construct_neighborhoods_ranged_plain_and_validation: indices for point 1")

        call construct_neighborhoods_ranged(0_int32, x_star, n_genes_S, mean_S, n_neighbors, &
                                            neighborhood_indices, neighborhood_range, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_EMPTY_INPUT, &
                              "test_construct_neighborhoods_ranged_plain_and_validation: n_points=0")

        call construct_neighborhoods_ranged(n_points, x_star, n_genes_S, mean_S, 0_int32, &
                                            neighborhood_indices, neighborhood_range, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_construct_neighborhoods_ranged_plain_and_validation: n_neighbors=0 rejected", &
                        arg_pos=5_int32)
    end subroutine test_construct_neighborhoods_ranged_plain_and_validation

    !> Feeding the same two studies both through `determine_study_shared_residual_range` (with
    !| S2 as-is) and through the N-study routine (S2 padded with NaN up to S1's replicate count,
    !| n_studies=2) must give the exact same range -- this is a regression-safety cross-check
    !| against the already-tested two-study routine, not just a hand-computed number.
    subroutine test_all_studies_range_matches_two_study_routine()
        integer(int32), parameter :: n_reps_S1 = 4, n_reps_S2 = 3, n_neighbors = 2, n_points = 2
        integer(int32), parameter :: n_studies = 2, max_n_reps = 4
        real(real64), dimension(n_reps_S1, n_neighbors, n_points) :: S1
        real(real64), dimension(n_reps_S2, n_neighbors, n_points) :: S2
        real(real64), dimension(max_n_reps, n_neighbors, n_points, n_studies) :: all_studies
        real(real64) :: R_two_study, R_all_studies
        integer(int32) :: ierr

        S1 = reshape([1, 2, 3, 4, 5, 6, -7, 8, 9, 10, 11, 12, 1, 1, 1, 1], shape(S1))
        S2 = reshape([2, -4, 6, 8, 1, 3, 5, 7, 9, 0, 1, 2], shape(S2))

        call determine_study_shared_residual_range(S1, S2, n_reps_S1, n_reps_S2, n_neighbors, n_points, &
                                                    R_two_study, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_all_studies_range_matches_two_study_routine: two-study ierr")

        all_studies = ieee_value(1.0_real64, ieee_quiet_nan)
        all_studies(:, :, :, 1) = S1
        all_studies(1:n_reps_S2, :, :, 2) = S2
        ! all_studies(4, :, :, 2) stays NaN -- padding S2 up to max_n_reps

        call determine_all_studies_shared_residual_range(all_studies, max_n_reps, n_neighbors, n_points, n_studies, &
                                                          R_all_studies, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_all_studies_range_matches_two_study_routine: N-study ierr")
        call assert_equal_real(R_all_studies, R_two_study, TOL, &
                               "test_all_studies_range_matches_two_study_routine: R mismatch")
        call assert_equal_real(R_all_studies, 10.65_real64, TOL, &
                               "test_all_studies_range_matches_two_study_routine: R should be 10.65")
    end subroutine test_all_studies_range_matches_two_study_routine

    !> Three single-replicate studies, hand-computed: pooled absolute residuals sorted are
    !| [3, 4, 5]; the default 95% quantile has rank `0.95*(3-1)+1 = 2.9`, so
    !| `R = sorted(2) + 0.9*(sorted(3)-sorted(2)) = 4 + 0.9*1 = 4.9`.
    subroutine test_determine_all_studies_shared_residual_range_three_studies()
        integer(int32), parameter :: n_neighbors = 1, n_points = 1, n_studies = 3, max_n_reps = 1
        real(real64), dimension(max_n_reps, n_neighbors, n_points, n_studies) :: all_studies
        real(real64) :: R
        integer(int32) :: ierr

        all_studies(1, 1, 1, 1) = 3.0_real64
        all_studies(1, 1, 1, 2) = -4.0_real64
        all_studies(1, 1, 1, 3) = 5.0_real64

        call determine_all_studies_shared_residual_range(all_studies, max_n_reps, n_neighbors, n_points, n_studies, &
                                                          R, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_all_studies_shared_residual_range_three_studies: ierr should be OK")
        call assert_equal_real(R, 4.9_real64, TOL, &
                               "test_determine_all_studies_shared_residual_range_three_studies: R should be 4.9")
    end subroutine test_determine_all_studies_shared_residual_range_three_studies

    !> Bounds on `residual_range_quantile` still apply through the N-study wrapper.
    subroutine test_all_studies_range_invalid_quantile()
        integer(int32), parameter :: n_neighbors = 1, n_points = 1, n_studies = 2, max_n_reps = 1
        real(real64), dimension(max_n_reps, n_neighbors, n_points, n_studies) :: all_studies
        real(real64) :: R, q
        integer(int32) :: ierr

        all_studies(1, 1, 1, 1) = 1.0_real64
        all_studies(1, 1, 1, 2) = 2.0_real64

        q = below(0.0_real64)
        call determine_all_studies_shared_residual_range(all_studies, max_n_reps, n_neighbors, n_points, n_studies, &
                                                          R, ierr=ierr, residual_range_quantile=q)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, &
                              "test_all_studies_range_invalid_quantile: below 0 rejected")

        q = above(1.0_real64)
        call determine_all_studies_shared_residual_range(all_studies, max_n_reps, n_neighbors, n_points, n_studies, &
                                                          R, ierr=ierr, residual_range_quantile=q)
        call assert_equal_int(get_err_code(ierr), ERR_INVALID_INPUT, &
                              "test_all_studies_range_invalid_quantile: above 1 rejected")
    end subroutine test_all_studies_range_invalid_quantile

    !> `calc_pmf` is the counts-to-pmf half of `build_residual_histograms`, factored out; feeding
    !| it the exact `counts`/`included_n_reps` that routine's own basic fixture produces must
    !| reproduce that fixture's expected pmf exactly.
    subroutine test_calc_pmf_matches_build_residual_histograms()
        integer(int32), parameter :: n_points = 3, n_bins = 4
        integer(int32), dimension(n_points, n_bins) :: counts
        integer(int32), dimension(n_points) :: included_n_reps
        real(real64), dimension(n_points, n_bins) :: pmf, expected_pmf
        integer(int32) :: ierr

        counts = reshape([2, 0, 1, 1, 0, 1, 2, 6, 2, 1, 0, 2], [n_points, n_bins])
        included_n_reps = [6, 6, 6]

        expected_pmf = reshape([ &
                               0.3333333333333333_real64, 0.0_real64, 0.16666666666666666_real64, &
                               0.16666666666666666_real64, 0.0_real64, 0.16666666666666666_real64, &
                               0.3333333333333333_real64, 1.0_real64, 0.3333333333333333_real64, &
                               0.16666666666666666_real64, 0.0_real64, 0.3333333333333333_real64 &
                               ], [n_points, n_bins])

        call calc_pmf(counts, included_n_reps, n_points, n_bins, pmf, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_pmf_matches_build_residual_histograms: ierr should be OK")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, &
                                     "test_calc_pmf_matches_build_residual_histograms: pmf mismatch")
    end subroutine test_calc_pmf_matches_build_residual_histograms

    !> A reference point with zero included replicates (every residual there was NaN) must report
    !| an all-zero pmf row rather than dividing by zero.
    subroutine test_calc_pmf_zero_included_reps()
        integer(int32), parameter :: n_points = 2, n_bins = 2
        integer(int32), dimension(n_points, n_bins) :: counts
        integer(int32), dimension(n_points) :: included_n_reps
        real(real64), dimension(n_points, n_bins) :: pmf, expected_pmf
        integer(int32) :: ierr

        counts = reshape([3, 0, 5, 0], [n_points, n_bins])
        included_n_reps = [8, 0]

        expected_pmf = reshape([0.375_real64, 0.0_real64, 0.625_real64, 0.0_real64], [n_points, n_bins])

        call calc_pmf(counts, included_n_reps, n_points, n_bins, pmf, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_calc_pmf_zero_included_reps: ierr should be OK")
        call assert_equal_array_real(pmf, expected_pmf, size(pmf, kind=int32), TOL, &
                                     "test_calc_pmf_zero_included_reps: pmf mismatch")
    end subroutine test_calc_pmf_zero_included_reps

    !> `counts` is documented as non-negative; the generated wrapper must reject a negative entry.
    subroutine test_calc_pmf_rejects_negative_counts()
        integer(int32), parameter :: n_points = 1, n_bins = 1
        integer(int32), dimension(n_points, n_bins) :: counts
        integer(int32), dimension(n_points) :: included_n_reps
        real(real64), dimension(n_points, n_bins) :: pmf
        integer(int32) :: ierr

        counts = reshape([-1], [n_points, n_bins])
        included_n_reps = [1]

        call calc_pmf(counts, included_n_reps, n_points, n_bins, pmf, ierr)

        call assert_err(ierr, ERR_INVALID_INPUT, "test_calc_pmf_rejects_negative_counts: negative count rejected", &
                        arg_pos=1_int32)
    end subroutine test_calc_pmf_rejects_negative_counts

    !> Hand-computed: `residuals = [-4,-3,-2,-1,0,1,2,3,4,5]` (already ascending, no ties needed),
    !| `max_n_reps_all_studies=1, n_neighbors=1` so `n_reps_neighborhood=1`. Sturges gives
    !| `1 + nint(log(1)/LOG_2) = 1`. The 25th/75th percentiles (rank `0.25*9+1=3.25` and
    !| `0.75*9+1=7.75`) interpolate to `-1.75` and `2.75`, so the Freedman-Diaconis half-width is
    !| `(2.75 - (-1.75)) / 1^(1/3) = 4.5`; with `shared_residual_range=9.0`,
    !| `nint(9.0/4.5) = 2`, which beats Sturges, so `n_bins` should be `2`.
    subroutine test_estimate_bin_count_basic_hand_computed()
        integer(int32), parameter :: n_residuals = 10
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, ierr

        residuals = [-4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
                    1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call estimate_bin_count(residuals, n_residuals, 1_int32, 1_int32, 9.0_real64, n_bins, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_estimate_bin_count_basic_hand_computed: ierr should be OK")
        call assert_equal_int(n_bins, 2_int32, "test_estimate_bin_count_basic_hand_computed: n_bins should be 2")
    end subroutine test_estimate_bin_count_basic_hand_computed

    !> When every residual is NaN, the pool is empty and the routine must fall back to a single
    !| bin rather than computing percentiles of nothing.
    subroutine test_estimate_bin_count_all_nan_gives_one_bin()
        integer(int32), parameter :: n_residuals = 4
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, ierr
        real(real64) :: nan_val

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)
        residuals = nan_val

        call estimate_bin_count(residuals, n_residuals, 1_int32, 1_int32, 9.0_real64, n_bins, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_estimate_bin_count_all_nan_gives_one_bin: ierr should be OK")
        call assert_equal_int(n_bins, 1_int32, "test_estimate_bin_count_all_nan_gives_one_bin: n_bins should be 1")
    end subroutine test_estimate_bin_count_all_nan_gives_one_bin

    !> Reusing the basic fixture's residuals but with `shared_residual_range=5000.0`, the
    !| Freedman-Diaconis half-width is unchanged (`4.5`), so the raw estimate is
    !| `nint(5000.0/4.5) = 1111`, far above MAX_N_BINS; the output must clamp to 256.
    subroutine test_estimate_bin_count_clamped_to_max_n_bins()
        integer(int32), parameter :: n_residuals = 10
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, ierr

        residuals = [-4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
                    1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call estimate_bin_count(residuals, n_residuals, 1_int32, 1_int32, 5000.0_real64, n_bins, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_estimate_bin_count_clamped_to_max_n_bins: ierr should be OK")
        call assert_equal_int(n_bins, 256_int32, "test_estimate_bin_count_clamped_to_max_n_bins: n_bins should clamp to 256")
    end subroutine test_estimate_bin_count_clamped_to_max_n_bins

    !> The documented small-N candidate-grid collapse, bracketed at its exact threshold: with
    !| `max_n_genes_all_studies=8742`, `n_points_high = clamp(ceil(4*sqrt(8742)), 300, 1500) = 374`
    !| and `n_points_low = max(300, ceil(0.2*374)) = 300`; after one grid iteration
    !| `n_points_high` becomes `374*0.8 = 299.2 < 300`, so the loop exits before a second distinct
    !| `n_points` value is ever produced -- every candidate pair the grid returns must share the
    !| same `n_points`. This is real, derived behavior the grid depends on, not a bug.
    subroutine test_generate_js_comp_test_candidates_collapses_at_8742()
        integer(int32), parameter :: n_residuals = 5
        real(real64) :: residuals(n_residuals)
        integer(int32) :: candidates(2, 16), n_bins_candidates(16), n_candidates, ierr

        residuals = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call generate_js_comp_test_candidates(8742_int32, residuals, n_residuals, 1_int32, 1.0_real64, &
                                              candidates, n_bins_candidates, n_candidates, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_generate_js_comp_test_candidates_collapses_at_8742: ierr should be OK")
        call assert_true(n_candidates >= 1_int32, &
                         "test_generate_js_comp_test_candidates_collapses_at_8742: at least one candidate")
        call assert_equal_int(count(candidates(1, 1:n_candidates) == candidates(1, 1)), n_candidates, &
                              "test_generate_js_comp_test_candidates_collapses_at_8742: "// &
                              "every candidate must share the same n_points")
        call assert_equal_int(candidates(1, 1), 374_int32, &
                              "test_generate_js_comp_test_candidates_collapses_at_8742: n_points should be 374")
    end subroutine test_generate_js_comp_test_candidates_collapses_at_8742

    !> The other side of the bracket: with `max_n_genes_all_studies=8743`,
    !| `n_points_high = ceil(4*sqrt(8743)) = 375` and `n_points_low = max(300, ceil(0.2*375)) = 300`;
    !| after one iteration `n_points_high` becomes `375*0.8 = 300.0`, which is NOT less than
    !| `n_points_low=300`, so the loop proceeds to a second, distinct `n_points=300` candidate. At
    !| least two distinct `n_points` values must appear in the grid.
    subroutine test_generate_js_comp_test_candidates_has_two_distinct_at_8743()
        integer(int32), parameter :: n_residuals = 5
        real(real64) :: residuals(n_residuals)
        integer(int32) :: candidates(2, 16), n_bins_candidates(16), n_candidates, ierr
        integer(int32) :: n_distinct_n_points, i_candidate

        residuals = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call generate_js_comp_test_candidates(8743_int32, residuals, n_residuals, 1_int32, 1.0_real64, &
                                              candidates, n_bins_candidates, n_candidates, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_generate_js_comp_test_candidates_has_two_distinct_at_8743: ierr should be OK")

        n_distinct_n_points = 0_int32
        do i_candidate = 1, n_candidates
            if (i_candidate == 1) then
                n_distinct_n_points = n_distinct_n_points + 1_int32
            else if (candidates(1, i_candidate) /= candidates(1, i_candidate - 1)) then
                n_distinct_n_points = n_distinct_n_points + 1_int32
            end if
        end do
        call assert_true(n_distinct_n_points >= 2_int32, &
                         "test_generate_js_comp_test_candidates_has_two_distinct_at_8743: "// &
                         "at least two distinct n_points values expected")
        call assert_equal_int(candidates(1, 1), 375_int32, &
                              "test_generate_js_comp_test_candidates_has_two_distinct_at_8743: "// &
                              "first n_points should be 375")
        call assert_equal_int(candidates(1, n_candidates), 300_int32, &
                              "test_generate_js_comp_test_candidates_has_two_distinct_at_8743: "// &
                              "last n_points should collapse to 300")
    end subroutine test_generate_js_comp_test_candidates_has_two_distinct_at_8743

    !> `max_n_genes_all_studies` must be positive, and the expert tier's `residuals_perm` is
    !| bounds-checked exactly like every other permutation argument.
    subroutine test_generate_js_comp_test_candidates_validation()
        integer(int32), parameter :: n_residuals = 5
        real(real64) :: residuals(n_residuals)
        integer(int32) :: residuals_perm(n_residuals)
        integer(int32) :: candidates(2, 16), n_bins_candidates(16), n_candidates, ierr

        residuals = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]
        residuals_perm = [1, 2, 3, 4, 5]

        call generate_js_comp_test_candidates(0_int32, residuals, n_residuals, 1_int32, 1.0_real64, &
                                              candidates, n_bins_candidates, n_candidates, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_generate_js_comp_test_candidates_validation: max_n_genes_all_studies=0 rejected", &
                        arg_pos=1_int32)

        call generate_js_comp_test_candidates_expert(100_int32, residuals, residuals_perm, n_residuals, 1_int32, &
                                                     1.0_real64, candidates, n_bins_candidates, n_candidates, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_generate_js_comp_test_candidates_validation: expert tier accepts a valid permutation")
    end subroutine test_generate_js_comp_test_candidates_validation

    !> Three reference points whose neighborhood ranges fully overlap (identical spans): the
    !| fractional overlap between any consecutive pair is exactly 1.0, so the gate must pass even
    !| at the strictest possible threshold.
    subroutine test_neighborhood_overlaps_all_pass()
        integer(int32), parameter :: n_points = 3
        integer(int32) :: neighborhood_range(2, n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        neighborhood_range(:, 1) = [1, 10]
        neighborhood_range(:, 2) = [1, 10]
        neighborhood_range(:, 3) = [1, 10]

        call check_neighborhood_overlaps(neighborhood_range, n_points, 1.0_real64, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_neighborhood_overlaps_all_pass: ierr should be OK")
        call assert_true(all_pass, "test_neighborhood_overlaps_all_pass: identical ranges should always pass")
    end subroutine test_neighborhood_overlaps_all_pass

    !> Two points with ranges `[1,10]` and `[6,10]`: `left_max=min(10,10)=10`,
    !| `right_min=max(1,6)=6`, so the overlap is exactly `(10-6)/(10-1) = 4/9`. Setting
    !| `min_neighbor_overlap` to that exact fraction must still pass, since the gate is `overlap
    !| >= min_neighbor_overlap`, not a strict `>`.
    subroutine test_neighborhood_overlaps_exactly_at_threshold_passes()
        integer(int32), parameter :: n_points = 2
        integer(int32) :: neighborhood_range(2, n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        neighborhood_range(:, 1) = [1, 10]
        neighborhood_range(:, 2) = [6, 10]

        call check_neighborhood_overlaps(neighborhood_range, n_points, 4.0_real64/9.0_real64, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_neighborhood_overlaps_exactly_at_threshold_passes: ierr should be OK")
        call assert_true(all_pass, "test_neighborhood_overlaps_exactly_at_threshold_passes: "// &
                         "overlap exactly at threshold must pass")
    end subroutine test_neighborhood_overlaps_exactly_at_threshold_passes

    !> The same `4/9` overlap as above, but with `min_neighbor_overlap=0.5 > 4/9`: the gate must
    !| fail.
    subroutine test_neighborhood_overlaps_below_threshold_fails()
        integer(int32), parameter :: n_points = 2
        integer(int32) :: neighborhood_range(2, n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        neighborhood_range(:, 1) = [1, 10]
        neighborhood_range(:, 2) = [6, 10]

        call check_neighborhood_overlaps(neighborhood_range, n_points, 0.5_real64, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_neighborhood_overlaps_below_threshold_fails: ierr should be OK")
        call assert_false(all_pass, "test_neighborhood_overlaps_below_threshold_fails: "// &
                          "overlap below threshold must fail")
    end subroutine test_neighborhood_overlaps_below_threshold_fails

    !> Every bin comfortably exceeds the minimum count.
    subroutine test_mean_pmf_min_counts_all_pass()
        integer(int32), parameter :: n_bins = 2, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        mean_pmf_counts = reshape([10, 10, 10, 10], [n_bins, n_points])

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_points, 5_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_mean_pmf_min_counts_all_pass: ierr should be OK")
        call assert_true(all_pass, "test_mean_pmf_min_counts_all_pass: all bins above minimum should pass")
    end subroutine test_mean_pmf_min_counts_all_pass

    !> Every bin equals the minimum count exactly: the gate is `count >= min`, not a strict `>`,
    !| so this must still pass.
    subroutine test_mean_pmf_min_counts_exactly_at_threshold_passes()
        integer(int32), parameter :: n_bins = 2, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        mean_pmf_counts = reshape([3, 3, 3, 3], [n_bins, n_points])

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_points, 3_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_mean_pmf_min_counts_exactly_at_threshold_passes: ierr should be OK")
        call assert_true(all_pass, "test_mean_pmf_min_counts_exactly_at_threshold_passes: "// &
                         "count exactly at minimum must pass")
    end subroutine test_mean_pmf_min_counts_exactly_at_threshold_passes

    !> One bin (of four) is one below the minimum: the gate must fail.
    subroutine test_mean_pmf_min_counts_below_threshold_fails()
        integer(int32), parameter :: n_bins = 2, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        mean_pmf_counts = reshape([3, 3, 3, 2], [n_bins, n_points])

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_points, 3_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_mean_pmf_min_counts_below_threshold_fails: ierr should be OK")
        call assert_false(all_pass, "test_mean_pmf_min_counts_below_threshold_fails: "// &
                          "one bin below minimum must fail the whole gate")
    end subroutine test_mean_pmf_min_counts_below_threshold_fails

    !> METHOD_JOIN_MIN succeeds only once every study's confidence-interval overlap exceeds the
    !| threshold. With 3 studies, 2-of-3 passing must fail, and 3-of-3 passing must succeed.
    !| `best_exceeded_ci_overlap_count` starts at 0 so the short-circuit never fires here.
    subroutine test_check_plateau_condition_join_min_requires_all()
        integer(int32), parameter :: n_studies = 3
        real(real64) :: confidence_interval(2, n_studies), best_ci(2, n_studies)
        integer(int32) :: best_candidate_index, best_exceeded_count, ierr
        logical(c_bool) :: plateau_found

        ! Study 1,2: candidate fully inside the wide best interval -> overlap 1.0 (passes).
        ! Study 3: candidate far outside the best interval -> overlap 0.0 (fails).
        confidence_interval(:, 1) = [0.4_real64, 0.6_real64]
        confidence_interval(:, 2) = [0.4_real64, 0.6_real64]
        confidence_interval(:, 3) = [-2.0_real64, -1.0_real64]
        best_ci(:, 1) = [0.0_real64, 1.0_real64]
        best_ci(:, 2) = [0.0_real64, 1.0_real64]
        best_ci(:, 3) = [0.0_real64, 1.0_real64]
        best_candidate_index = 0_int32
        best_exceeded_count = 0_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 1_int32, METHOD_JOIN_MIN, 0.9_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_join_min_requires_all: 2-of-3 ierr should be OK")
        call assert_false(plateau_found, &
                          "test_check_plateau_condition_join_min_requires_all: 2-of-3 studies must not plateau")

        ! Now all three studies pass. best_ci must be reset: the first call's else-branch already
        ! overwrote it with the first call's own confidence_interval.
        confidence_interval(:, 3) = [0.4_real64, 0.6_real64]
        best_ci(:, 1) = [0.0_real64, 1.0_real64]
        best_ci(:, 2) = [0.0_real64, 1.0_real64]
        best_ci(:, 3) = [0.0_real64, 1.0_real64]
        best_candidate_index = 0_int32
        best_exceeded_count = 0_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 1_int32, METHOD_JOIN_MIN, 0.9_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_join_min_requires_all: 3-of-3 ierr should be OK")
        call assert_true(plateau_found, &
                         "test_check_plateau_condition_join_min_requires_all: 3-of-3 studies must plateau")
    end subroutine test_check_plateau_condition_join_min_requires_all

    !> METHOD_JOIN_MAX succeeds once any one study's overlap exceeds the threshold. With 0-of-3
    !| passing it must fail, and with exactly 1-of-3 passing it must succeed.
    subroutine test_check_plateau_condition_join_max_requires_any()
        integer(int32), parameter :: n_studies = 3
        real(real64) :: confidence_interval(2, n_studies), best_ci(2, n_studies)
        integer(int32) :: best_candidate_index, best_exceeded_count, ierr
        logical(c_bool) :: plateau_found

        confidence_interval(:, 1) = [-2.0_real64, -1.0_real64]
        confidence_interval(:, 2) = [-2.0_real64, -1.0_real64]
        confidence_interval(:, 3) = [-2.0_real64, -1.0_real64]
        best_ci(:, 1) = [0.0_real64, 1.0_real64]
        best_ci(:, 2) = [0.0_real64, 1.0_real64]
        best_ci(:, 3) = [0.0_real64, 1.0_real64]
        best_candidate_index = 0_int32
        best_exceeded_count = 0_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 1_int32, METHOD_JOIN_MAX, 0.9_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_join_max_requires_any: 0-of-3 ierr should be OK")
        call assert_false(plateau_found, &
                          "test_check_plateau_condition_join_max_requires_any: 0-of-3 studies must not plateau")

        confidence_interval(:, 1) = [0.4_real64, 0.6_real64]
        best_ci(:, 1) = [0.0_real64, 1.0_real64]
        best_ci(:, 2) = [0.0_real64, 1.0_real64]
        best_ci(:, 3) = [0.0_real64, 1.0_real64]
        best_candidate_index = 0_int32
        best_exceeded_count = 0_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 1_int32, METHOD_JOIN_MAX, 0.9_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_join_max_requires_any: 1-of-3 ierr should be OK")
        call assert_true(plateau_found, &
                         "test_check_plateau_condition_join_max_requires_any: 1-of-3 studies must plateau")
    end subroutine test_check_plateau_condition_join_max_requires_any

    !> METHOD_JOIN_MEDIAN succeeds once a majority pass: with 3 studies the threshold is
    !| `(3-1)/2 + 1 = 2`. 1-of-3 must fail, 2-of-3 must succeed.
    subroutine test_check_plateau_condition_join_median_requires_majority()
        integer(int32), parameter :: n_studies = 3
        real(real64) :: confidence_interval(2, n_studies), best_ci(2, n_studies)
        integer(int32) :: best_candidate_index, best_exceeded_count, ierr
        logical(c_bool) :: plateau_found

        confidence_interval(:, 1) = [0.4_real64, 0.6_real64]
        confidence_interval(:, 2) = [-2.0_real64, -1.0_real64]
        confidence_interval(:, 3) = [-2.0_real64, -1.0_real64]
        best_ci(:, 1) = [0.0_real64, 1.0_real64]
        best_ci(:, 2) = [0.0_real64, 1.0_real64]
        best_ci(:, 3) = [0.0_real64, 1.0_real64]
        best_candidate_index = 0_int32
        best_exceeded_count = 0_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 1_int32, METHOD_JOIN_MEDIAN, 0.9_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_join_median_requires_majority: 1-of-3 ierr should be OK")
        call assert_false(plateau_found, &
                          "test_check_plateau_condition_join_median_requires_majority: 1-of-3 must not plateau")

        confidence_interval(:, 2) = [0.4_real64, 0.6_real64]
        best_ci(:, 1) = [0.0_real64, 1.0_real64]
        best_ci(:, 2) = [0.0_real64, 1.0_real64]
        best_ci(:, 3) = [0.0_real64, 1.0_real64]
        best_candidate_index = 0_int32
        best_exceeded_count = 0_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 1_int32, METHOD_JOIN_MEDIAN, 0.9_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_join_median_requires_majority: 2-of-3 ierr should be OK")
        call assert_true(plateau_found, &
                         "test_check_plateau_condition_join_median_requires_majority: 2-of-3 must plateau")
    end subroutine test_check_plateau_condition_join_median_requires_majority

    !> When the new candidate exceeds fewer studies' overlap threshold than the previous best did
    !| (`best_exceeded_ci_overlap_count=2` initially, the new candidate passes 0 studies),
    !| `plateau_found` must be `.true.` immediately -- regardless of `join_method` -- and none of
    !| the `best_*` state may be overwritten, so the caller keeps the previous best candidate.
    subroutine test_check_plateau_condition_worse_than_previous_short_circuits()
        integer(int32), parameter :: n_studies = 2
        real(real64) :: confidence_interval(2, n_studies), best_ci(2, n_studies), expected_best_ci(2, n_studies)
        integer(int32) :: best_candidate_index, best_exceeded_count, ierr
        logical(c_bool) :: plateau_found

        ! New candidate: both studies clearly fail the overlap check against best_ci below.
        confidence_interval(:, 1) = [-5.0_real64, -4.0_real64]
        confidence_interval(:, 2) = [-5.0_real64, -4.0_real64]
        best_ci(:, 1) = [0.1_real64, 0.2_real64]
        best_ci(:, 2) = [0.3_real64, 0.4_real64]
        expected_best_ci = best_ci
        best_candidate_index = 5_int32
        best_exceeded_count = 2_int32

        call check_plateau_condition(confidence_interval, best_ci, n_studies, best_candidate_index, &
                                     best_exceeded_count, 9_int32, METHOD_JOIN_MIN, 0.5_real64, plateau_found, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_check_plateau_condition_worse_than_previous_short_circuits: ierr should be OK")
        call assert_true(plateau_found, &
                         "test_check_plateau_condition_worse_than_previous_short_circuits: "// &
                         "a worse candidate must plateau immediately")
        call assert_equal_int(best_candidate_index, 5_int32, &
                              "test_check_plateau_condition_worse_than_previous_short_circuits: "// &
                              "best_candidate_index must not be overwritten")
        call assert_equal_int(best_exceeded_count, 2_int32, &
                              "test_check_plateau_condition_worse_than_previous_short_circuits: "// &
                              "best_exceeded_ci_overlap_count must not be overwritten")
        call assert_equal_array_real(best_ci, expected_best_ci, size(best_ci, kind=int32), TOL, &
                                     "test_check_plateau_condition_worse_than_previous_short_circuits: "// &
                                     "best_candidate_pair_confidence_interval must not be overwritten")
    end subroutine test_check_plateau_condition_worse_than_previous_short_circuits

    !> `has_previous = .false.` (the very first admissible candidate): no transition exists, so
    !| `delta`/`delta_median`/`delta_max` must all be the `-1.0` sentinel, `n_consecutive_ok` must
    !| stay/reset at 0, and `plateau_found` must be `.false.` -- even though `n_consecutive_ok` is
    !| deliberately seeded non-zero on entry, to prove `has_previous = .false.` resets it rather
    !| than merely leaving it untouched.
    subroutine test_effect_size_plateau_first_candidate_no_delta()
        integer(int32), parameter :: n_studies = 2
        real(real64) :: global_js_divergence(n_studies), prev_global_js_divergence(n_studies), delta(n_studies)
        real(real64) :: delta_median, delta_max
        integer(int32) :: n_consecutive_ok, ierr
        logical(c_bool) :: plateau_found

        global_js_divergence = [0.1_real64, 0.2_real64]
        prev_global_js_divergence = [0.0_real64, 0.0_real64]
        n_consecutive_ok = 5_int32

        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.false., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_effect_size_plateau_first_candidate_no_delta: ierr should be OK")
        call assert_equal_array_real(delta, [-1.0_real64, -1.0_real64], n_studies, TOL, &
                                     "test_effect_size_plateau_first_candidate_no_delta: "// &
                                     "delta is the -1.0 sentinel")
        call assert_equal_real(delta_median, -1.0_real64, TOL, &
                              "test_effect_size_plateau_first_candidate_no_delta: "// &
                              "delta_median is the -1.0 sentinel")
        call assert_equal_real(delta_max, -1.0_real64, TOL, &
                              "test_effect_size_plateau_first_candidate_no_delta: "// &
                              "delta_max is the -1.0 sentinel")
        call assert_equal_int(n_consecutive_ok, 0_int32, &
                              "test_effect_size_plateau_first_candidate_no_delta: "// &
                              "n_consecutive_ok is reset to 0, not left at its seeded 5")
        call assert_false(plateau_found, &
                          "test_effect_size_plateau_first_candidate_no_delta: "// &
                          "no plateau on the very first admissible candidate")
    end subroutine test_effect_size_plateau_first_candidate_no_delta

    !> One qualifying transition (both studies' relative JSD change under threshold) is not enough
    !| on its own: `delta_min_consecutive_transitions=2` requires two IN A ROW, so `plateau_found`
    !| must be `.false.` and `n_consecutive_ok` must be exactly 1 after this single call.
    subroutine test_effect_size_plateau_single_transition_insufficient()
        integer(int32), parameter :: n_studies = 2
        real(real64) :: global_js_divergence(n_studies), prev_global_js_divergence(n_studies), delta(n_studies)
        real(real64) :: delta_median, delta_max
        integer(int32) :: n_consecutive_ok, ierr
        logical(c_bool) :: plateau_found

        ! Both studies: relative change = |0.101-0.100|/0.100 = 0.01, well under both thresholds.
        prev_global_js_divergence = [0.100_real64, 0.100_real64]
        global_js_divergence = [0.101_real64, 0.101_real64]
        n_consecutive_ok = 0_int32

        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_effect_size_plateau_single_transition_insufficient: ierr should be OK")
        call assert_equal_int(n_consecutive_ok, 1_int32, &
                              "test_effect_size_plateau_single_transition_insufficient: "// &
                              "one qualifying transition recorded")
        call assert_false(plateau_found, &
                          "test_effect_size_plateau_single_transition_insufficient: "// &
                          "one transition alone must not plateau")
    end subroutine test_effect_size_plateau_single_transition_insufficient

    !> Two consecutive qualifying transitions -- exactly `delta_min_consecutive_transitions` -- must
    !| declare a plateau on the second call.
    subroutine test_effect_size_plateau_two_consecutive_transitions()
        integer(int32), parameter :: n_studies = 2
        real(real64) :: global_js_divergence(n_studies), prev_global_js_divergence(n_studies), delta(n_studies)
        real(real64) :: delta_median, delta_max
        integer(int32) :: n_consecutive_ok, ierr
        logical(c_bool) :: plateau_found

        prev_global_js_divergence = [0.100_real64, 0.100_real64]
        global_js_divergence = [0.101_real64, 0.101_real64]
        n_consecutive_ok = 0_int32

        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)
        call assert_false(plateau_found, &
                          "test_effect_size_plateau_two_consecutive_transitions: "// &
                          "first of two must not yet plateau")

        prev_global_js_divergence = global_js_divergence
        global_js_divergence = [0.102_real64, 0.102_real64]

        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_effect_size_plateau_two_consecutive_transitions: "// &
                              "ierr should be OK")
        call assert_equal_int(n_consecutive_ok, 2_int32, &
                              "test_effect_size_plateau_two_consecutive_transitions: "// &
                              "two consecutive qualifying transitions recorded")
        call assert_true(plateau_found, &
                         "test_effect_size_plateau_two_consecutive_transitions: "// &
                         "two consecutive qualifying transitions must plateau")
    end subroutine test_effect_size_plateau_two_consecutive_transitions

    !> A qualifying transition followed by a NON-qualifying one must reset the counter to 0, so a
    !| third qualifying transition right after is only the first of a new streak, not the second.
    subroutine test_effect_size_plateau_resets_on_non_qualifying()
        integer(int32), parameter :: n_studies = 2
        real(real64) :: global_js_divergence(n_studies), prev_global_js_divergence(n_studies), delta(n_studies)
        real(real64) :: delta_median, delta_max
        integer(int32) :: n_consecutive_ok, ierr
        logical(c_bool) :: plateau_found

        ! Transition 1: qualifies (relative change 0.01).
        prev_global_js_divergence = [0.100_real64, 0.100_real64]
        global_js_divergence = [0.101_real64, 0.101_real64]
        n_consecutive_ok = 0_int32
        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)
        call assert_equal_int(n_consecutive_ok, 1_int32, &
                              "test_effect_size_plateau_resets_on_non_qualifying: "// &
                              "transition 1 qualifies")

        ! Transition 2: a large jump -- relative change 1.0, well past both thresholds -> resets.
        prev_global_js_divergence = global_js_divergence
        global_js_divergence = [0.202_real64, 0.202_real64]
        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)
        call assert_equal_int(n_consecutive_ok, 0_int32, &
                              "test_effect_size_plateau_resets_on_non_qualifying: "// &
                              "transition 2 does not qualify -> counter reset")
        call assert_false(plateau_found, &
                          "test_effect_size_plateau_resets_on_non_qualifying: "// &
                          "must not plateau right after a reset")

        ! Transition 3: qualifies again, but this is only the FIRST of a new streak.
        prev_global_js_divergence = global_js_divergence
        global_js_divergence = [0.203_real64, 0.203_real64]
        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_effect_size_plateau_resets_on_non_qualifying: "// &
                              "ierr should be OK")
        call assert_equal_int(n_consecutive_ok, 1_int32, &
                              "test_effect_size_plateau_resets_on_non_qualifying: "// &
                              "transition 3 is only the first of a NEW streak")
        call assert_false(plateau_found, &
                          "test_effect_size_plateau_resets_on_non_qualifying: "// &
                          "one transition of a new streak must not plateau")
    end subroutine test_effect_size_plateau_resets_on_non_qualifying

    !> Three studies with hand-computable, distinct relative changes:
    !| `prev=[0.10,0.20,0.40]`, `current=[0.11,0.24,0.60]` ->
    !| `delta = [0.01/0.10, 0.04/0.20, 0.20/0.40] = [0.1, 0.2, 0.5]`.
    !| Median of `[0.1,0.2,0.5]` is `0.2`, max is `0.5` -- both hand-computed, independent of
    !| `calc_percentile_impl`'s own internals.
    subroutine test_effect_size_plateau_median_max_hand_computed()
        integer(int32), parameter :: n_studies = 3
        real(real64) :: global_js_divergence(n_studies), prev_global_js_divergence(n_studies), delta(n_studies)
        real(real64) :: delta_median, delta_max
        integer(int32) :: n_consecutive_ok, ierr
        logical(c_bool) :: plateau_found

        prev_global_js_divergence = [0.10_real64, 0.20_real64, 0.40_real64]
        global_js_divergence = [0.11_real64, 0.24_real64, 0.60_real64]
        n_consecutive_ok = 0_int32

        call check_effect_size_plateau_condition(global_js_divergence, prev_global_js_divergence, n_studies, &
                                                 logical(.true., kind=c_bool), 0.05_real64, 0.10_real64, 1.0e-10_real64, &
                                                 2_int32, n_consecutive_ok, delta, delta_median, delta_max, plateau_found, &
                                                 ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_effect_size_plateau_median_max_hand_computed: ierr should be OK")
        call assert_equal_array_real(delta, [0.1_real64, 0.2_real64, 0.5_real64], n_studies, TOL, &
                                     "test_effect_size_plateau_median_max_hand_computed: "// &
                                     "per-study delta hand-computed")
        call assert_equal_real(delta_median, 0.2_real64, TOL, &
                              "test_effect_size_plateau_median_max_hand_computed: "// &
                              "delta_median hand-computed")
        call assert_equal_real(delta_max, 0.5_real64, TOL, &
                              "test_effect_size_plateau_median_max_hand_computed: "// &
                              "delta_max hand-computed")
        ! Both well above the default thresholds (0.05/0.10) -> does not qualify as a transition.
        call assert_equal_int(n_consecutive_ok, 0_int32, &
                              "test_effect_size_plateau_median_max_hand_computed: "// &
                              "large delta does not qualify as a plateau transition")
    end subroutine test_effect_size_plateau_median_max_hand_computed

    !> Two studies, one reference point, two bins, hand-computed:
    !| S1 pmf=[0.3,0.7] counts=[3,7] included=10; S2 pmf=[0.5,0.5] counts=[5,5] included=10.
    !| `mean_pmf = ([0.3,0.7]+[0.5,0.5])/2 = [0.4,0.6]`, `mean_pmf_counts = [3+5,7+5] = [8,12]`,
    !| `mean_pmf_included_n_reps = 10+10 = 20`.
    subroutine test_create_mean_pmf_two_studies_hand_computed()
        integer(int32), parameter :: n_bins = 2, n_points = 1, n_studies = 2
        real(real64) :: pmfs(n_bins, n_points, n_studies), mean_pmf(n_bins, n_points), expected_mean_pmf(n_bins, n_points)
        integer(int32) :: counts(n_bins, n_points, n_studies), mean_pmf_counts(n_bins, n_points)
        integer(int32) :: expected_mean_pmf_counts(n_bins, n_points)
        integer(int32) :: included_n_reps(n_points, n_studies), mean_pmf_included_n_reps(n_points)
        integer(int32) :: ierr

        pmfs(:, 1, 1) = [0.3_real64, 0.7_real64]
        pmfs(:, 1, 2) = [0.5_real64, 0.5_real64]
        counts(:, 1, 1) = [3, 7]
        counts(:, 1, 2) = [5, 5]
        included_n_reps(1, :) = [10, 10]

        expected_mean_pmf(:, 1) = [0.4_real64, 0.6_real64]
        expected_mean_pmf_counts(:, 1) = [8, 12]

        call create_mean_pmf(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                             mean_pmf_included_n_reps, mean_pmf_counts, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_create_mean_pmf_two_studies_hand_computed: ierr should be OK")
        call assert_equal_array_real(mean_pmf, expected_mean_pmf, size(mean_pmf, kind=int32), TOL, &
                                     "test_create_mean_pmf_two_studies_hand_computed: mean_pmf mismatch")
        call assert_equal_array_int(mean_pmf_counts, expected_mean_pmf_counts, size(mean_pmf_counts, kind=int32), &
                                    "test_create_mean_pmf_two_studies_hand_computed: mean_pmf_counts mismatch")
        call assert_equal_int(mean_pmf_included_n_reps(1), 20_int32, &
                              "test_create_mean_pmf_two_studies_hand_computed: mean_pmf_included_n_reps mismatch")
    end subroutine test_create_mean_pmf_two_studies_hand_computed

    !> Three studies, one reference point, two bins, hand-computed, and this test DOCUMENTS THE
    !| KNOWN LIMITATION rather than a desired property -- do not "fix" it by making this test
    !| expect a leave-one-out average without also deliberately updating the implementation and
    !| this comment together.
    !|
    !| S1 pmf=[0.2,0.8], S2 pmf=[0.4,0.6], S3 pmf=[0.9,0.1] (counts/included analogous, scaled by
    !| 10). `mean_pmf = ([0.2,0.8]+[0.4,0.6]+[0.9,0.1])/3 = [0.5,0.5]` -- computed from ALL THREE
    !| studies. A true leave-one-out background for study 1 would average only S2 and S3:
    !| `([0.4,0.6]+[0.9,0.1])/2 = [0.65,0.35]`, which is NOT what this routine has any way to
    !| produce (it has no i_study/exclude argument, by design -- see the plan). This test asserts
    !| the ALL-STUDIES-INCLUDING-SELF result, not the leave-one-out one.
    subroutine test_create_mean_pmf_three_studies_includes_self()
        integer(int32), parameter :: n_bins = 2, n_points = 1, n_studies = 3
        real(real64) :: pmfs(n_bins, n_points, n_studies), mean_pmf(n_bins, n_points), expected_mean_pmf(n_bins, n_points)
        real(real64) :: leave_one_out_mean_pmf_study1(n_bins)
        integer(int32) :: counts(n_bins, n_points, n_studies), mean_pmf_counts(n_bins, n_points)
        integer(int32) :: expected_mean_pmf_counts(n_bins, n_points)
        integer(int32) :: included_n_reps(n_points, n_studies), mean_pmf_included_n_reps(n_points)
        integer(int32) :: ierr

        pmfs(:, 1, 1) = [0.2_real64, 0.8_real64]
        pmfs(:, 1, 2) = [0.4_real64, 0.6_real64]
        pmfs(:, 1, 3) = [0.9_real64, 0.1_real64]
        counts(:, 1, 1) = [2, 8]
        counts(:, 1, 2) = [4, 6]
        counts(:, 1, 3) = [9, 1]
        included_n_reps(1, :) = [10, 10, 10]

        expected_mean_pmf(:, 1) = [0.5_real64, 0.5_real64]
        expected_mean_pmf_counts(:, 1) = [15, 15]

        call create_mean_pmf(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                             mean_pmf_included_n_reps, mean_pmf_counts, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_create_mean_pmf_three_studies_includes_self: ierr should be OK")
        call assert_equal_array_real(mean_pmf, expected_mean_pmf, size(mean_pmf, kind=int32), TOL, &
                                     "test_create_mean_pmf_three_studies_includes_self: "// &
                                     "mean_pmf must be the mean over ALL n_studies including self "// &
                                     "(known limitation, not leave-one-out)")
        call assert_equal_array_int(mean_pmf_counts, expected_mean_pmf_counts, size(mean_pmf_counts, kind=int32), &
                                    "test_create_mean_pmf_three_studies_includes_self: "// &
                                    "mean_pmf_counts mismatch")
        call assert_equal_int(mean_pmf_included_n_reps(1), 30_int32, &
                              "test_create_mean_pmf_three_studies_includes_self: "// &
                              "mean_pmf_included_n_reps mismatch")

        ! What a true leave-one-out background for study 1 would be, for contrast -- NOT what
        ! mean_pmf above equals, and not asserted against it.
        leave_one_out_mean_pmf_study1 = ([0.4_real64, 0.6_real64] + [0.9_real64, 0.1_real64])/2.0_real64
        call assert_true(abs(mean_pmf(1, 1) - leave_one_out_mean_pmf_study1(1)) > TOL, &
                         "test_create_mean_pmf_three_studies_includes_self: "// &
                         "mean_pmf must differ from the leave-one-out background -- if this ever "// &
                         "fails, the known limitation was fixed without updating this test")
    end subroutine test_create_mean_pmf_three_studies_includes_self

    !> `create_mean_pmf_only` must produce exactly the same `mean_pmf` as `create_mean_pmf` for the
    !| same `pmfs` input, since it is the same averaging computation with the counts/included_n_reps
    !| bookkeeping dropped.
    subroutine test_create_mean_pmf_only_matches_create_mean_pmf()
        integer(int32), parameter :: n_bins = 2, n_points = 1, n_studies = 2
        real(real64) :: pmfs(n_bins, n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points), mean_pmf_only(n_bins, n_points)
        integer(int32) :: counts(n_bins, n_points, n_studies), mean_pmf_counts(n_bins, n_points)
        integer(int32) :: included_n_reps(n_points, n_studies), mean_pmf_included_n_reps(n_points)
        integer(int32) :: ierr

        pmfs(:, 1, 1) = [0.3_real64, 0.7_real64]
        pmfs(:, 1, 2) = [0.5_real64, 0.5_real64]
        counts(:, 1, 1) = [3, 7]
        counts(:, 1, 2) = [5, 5]
        included_n_reps(1, :) = [10, 10]

        call create_mean_pmf(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                             mean_pmf_included_n_reps, mean_pmf_counts, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_create_mean_pmf_only_matches_create_mean_pmf: create_mean_pmf ierr should be OK")

        call create_mean_pmf_only(pmfs, n_bins, n_points, n_studies, mean_pmf_only, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_create_mean_pmf_only_matches_create_mean_pmf: create_mean_pmf_only ierr should be OK")

        call assert_equal_array_real(mean_pmf_only, mean_pmf, size(mean_pmf, kind=int32), TOL, &
                                     "test_create_mean_pmf_only_matches_create_mean_pmf: mean_pmf mismatch")
    end subroutine test_create_mean_pmf_only_matches_create_mean_pmf

    !> `pmfs` is documented as a probability in `[0,1]`, and `counts` as non-negative; the
    !| generated wrapper must reject a violation of either.
    subroutine test_create_mean_pmf_validation()
        integer(int32), parameter :: n_bins = 1, n_points = 1, n_studies = 1
        real(real64) :: pmfs(n_bins, n_points, n_studies), mean_pmf(n_bins, n_points)
        integer(int32) :: counts(n_bins, n_points, n_studies), mean_pmf_counts(n_bins, n_points)
        integer(int32) :: included_n_reps(n_points, n_studies), mean_pmf_included_n_reps(n_points)
        integer(int32) :: ierr

        pmfs(:, 1, 1) = [1.5_real64]
        counts(:, 1, 1) = [1]
        included_n_reps(1, :) = [1]

        call create_mean_pmf(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                             mean_pmf_included_n_reps, mean_pmf_counts, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_create_mean_pmf_validation: pmf above 1.0 rejected", &
                        arg_pos=1_int32)

        pmfs(:, 1, 1) = [0.5_real64]
        counts(:, 1, 1) = [-1]

        call create_mean_pmf(pmfs, counts, n_bins, n_points, n_studies, included_n_reps, mean_pmf, &
                             mean_pmf_included_n_reps, mean_pmf_counts, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, "test_create_mean_pmf_validation: negative count rejected", &
                        arg_pos=2_int32)
    end subroutine test_create_mean_pmf_validation

    !> Hand-computed: `n_bootstraps=1000` with the default significance level (2.5%, not passed)
    !| gives `n_top_k = max(1, floor(0.025 * 1000)) = 25`.
    subroutine test_calc_js_comp_test_n_top_k_jsds_default_hand_computed()
        integer(int32) :: n_top_k

        call calc_js_comp_test_n_top_k_jsds(1000_int32, n_top_k=n_top_k)

        call assert_equal_int(n_top_k, 25_int32, &
                              "test_calc_js_comp_test_n_top_k_jsds_default_hand_computed: "// &
                              "default 2.5% of 1000 bootstraps should give 25")
    end subroutine test_calc_js_comp_test_n_top_k_jsds_default_hand_computed

    !> Hand-computed with an explicit significance level: `n_bootstraps=200`, `significance=5.0`
    !| gives `n_top_k = max(1, floor(0.05 * 200)) = 10`.
    subroutine test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level()
        integer(int32) :: n_top_k

        call calc_js_comp_test_n_top_k_jsds(200_int32, two_sided_bootstrapping_significance_level=5.0_real64, &
                                            n_top_k=n_top_k)

        call assert_equal_int(n_top_k, 10_int32, &
                              "test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level: "// &
                              "5% of 200 bootstraps should give 10")
    end subroutine test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level

    !> `n_bootstraps=10` with the default 2.5% significance level gives `floor(0.025*10) = floor(0.25) = 0`,
    !| which must clamp up to the documented minimum of 1, never 0 (an empty heap).
    subroutine test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one()
        integer(int32) :: n_top_k

        call calc_js_comp_test_n_top_k_jsds(10_int32, n_top_k=n_top_k)

        call assert_equal_int(n_top_k, 1_int32, &
                              "test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one: "// &
                              "floor(0.025*10)=0 must clamp up to 1")
    end subroutine test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one

    !> Calling `bootstrap_histogram` twice with the same `random_seed` and otherwise identical
    !| inputs must produce bit-for-bit identical `confidence_interval` output: the GSL stream is
    !| deterministic once seeded, so nothing here should be able to introduce drift between calls.
    subroutine test_bootstrap_histogram_seeded_reproducibility()
        integer(int32), parameter :: n_bins = 3, n_points = 2, n_studies = 2, n_bootstraps = 15
        integer(int32) :: mean_pmf_counts(n_bins, n_points), mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: ci_first(2, n_studies), ci_second(2, n_studies)
        integer(int32) :: ierr

        mean_pmf_counts(:, 1) = [5, 3, 2]
        mean_pmf_counts(:, 2) = [2, 2, 6]
        mean_pmf_included_n_reps = [10, 10]
        included_n_reps(:, 1) = [4, 4]
        included_n_reps(:, 2) = [6, 6]

        ci_first(:, 1) = [0.2_real64, 0.2_real64]
        ci_first(:, 2) = [0.3_real64, 0.3_real64]

        call bootstrap_histogram(n_bootstraps, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf_included_n_reps, &
                                 included_n_reps, ci_first, ierr=ierr, random_seed=123_int32)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_bootstrap_histogram_seeded_reproducibility: first call ierr should be OK")

        ci_second(:, 1) = [0.2_real64, 0.2_real64]
        ci_second(:, 2) = [0.3_real64, 0.3_real64]

        call bootstrap_histogram(n_bootstraps, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf_included_n_reps, &
                                 included_n_reps, ci_second, ierr=ierr, random_seed=123_int32)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_bootstrap_histogram_seeded_reproducibility: second call ierr should be OK")

        call assert_equal_array_real(ci_first, ci_second, size(ci_first, kind=int32), TOL, &
                                     "test_bootstrap_histogram_seeded_reproducibility: "// &
                                     "same random_seed twice must give identical confidence_interval")
    end subroutine test_bootstrap_histogram_seeded_reproducibility

    !> Degenerate, exactly hand-computable case: a single histogram bin means every resample's pmf
    !| is `[1.0]`, for both a study and the (resampled) consensus mean alike, so the Jensen-Shannon
    !| divergence between them is exactly 0.0 on every single bootstrap draw, regardless of which
    !| counts random_multinomial happens to draw. Seeding `confidence_interval` at `[0.0, 0.0]`
    !| (the true observed value, since the un-bootstrapped JSD is 0.0 too for the same reason)
    !| means neither heap push (`value > heap(1)` / `value < heap(1)`) is ever strictly true, so the
    !| confidence interval must collapse to exactly `[0.0, 0.0]` for every study.
    subroutine test_bootstrap_histogram_degenerate_single_bin_zero_ci()
        integer(int32), parameter :: n_bins = 1, n_points = 1, n_studies = 2, n_bootstraps = 5
        integer(int32) :: mean_pmf_counts(n_bins, n_points), mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: confidence_interval(2, n_studies), expected(2, n_studies)
        integer(int32) :: ierr

        mean_pmf_counts(1, 1) = 10_int32
        mean_pmf_included_n_reps = [10_int32]
        included_n_reps(1, :) = [5_int32, 5_int32]

        confidence_interval = 0.0_real64
        expected = 0.0_real64

        call bootstrap_histogram(n_bootstraps, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf_included_n_reps, &
                                 included_n_reps, confidence_interval, ierr=ierr, random_seed=7_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_bootstrap_histogram_degenerate_single_bin_zero_ci: ierr should be OK")
        call assert_equal_array_real(confidence_interval, expected, size(confidence_interval, kind=int32), TOL, &
                                     "test_bootstrap_histogram_degenerate_single_bin_zero_ci: "// &
                                     "a single-bin histogram always has JSD=0, so the CI must collapse to [0,0]")
    end subroutine test_bootstrap_histogram_degenerate_single_bin_zero_ci

    !> Basic conservation-of-counts: a single reference point, single study, single permutation --
    !| `tmp_counts` (the last permutation's resampled histogram for the last study processed) must
    !| sum to exactly `included_n_reps` for that study/point, since `random_multiv_hypergeom` draws
    !| without replacement and every drawn element must land in exactly one bin.
    subroutine test_gjct_permutation_test_conservation_of_counts()
        integer(int32), parameter :: n_bins = 3, n_points = 1, n_studies = 1, n_permutations = 1
        integer(int32) :: mean_pmf_counts(n_bins, n_points), mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points), global_jsd_observed(n_studies), p_values(n_studies)
        integer(int32) :: tmp_mean_pmf_counts(n_bins, n_points), tmp_counts(n_bins, n_points)
        real(real64) :: tmp_pmfs(n_bins, n_points, n_studies)
        real(real64) :: tmp_js_divergences(n_points, n_studies), tmp_weights(n_points, n_studies)
        real(real64) :: tmp_global_js_divergence(n_studies)
        real(real64) :: tmp_pmf_point_major(n_points, n_bins)
        integer(int32) :: tmp_counts_point_major(n_points, n_bins)
        integer(int32) :: ierr

        mean_pmf_counts(:, 1) = [5, 3, 2]
        mean_pmf_included_n_reps = [10]
        included_n_reps(1, :) = [4]
        mean_pmf(:, 1) = real(mean_pmf_counts(:, 1), real64)/10.0_real64
        global_jsd_observed = 0.0_real64

        call gjct_permutation_test_expert(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                          mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values, &
                                          tmp_mean_pmf_counts, tmp_counts, tmp_pmfs, tmp_js_divergences, tmp_weights, &
                                          tmp_global_js_divergence, tmp_pmf_point_major, tmp_counts_point_major, &
                                          ierr=ierr, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_gjct_permutation_test_conservation_of_counts: ierr should be OK")
        call assert_equal_int(sum(tmp_counts(:, 1)), included_n_reps(1, 1), &
                              "test_gjct_permutation_test_conservation_of_counts: "// &
                              "resampled counts must sum to the study's per-point draw size")
    end subroutine test_gjct_permutation_test_conservation_of_counts

    !> Calling `gjct_permutation_test` twice with the same `random_seed` and otherwise identical
    !| inputs must produce bit-for-bit identical `p_values`: the GSL stream is deterministic once
    !| seeded, so nothing here should be able to introduce drift between calls.
    subroutine test_gjct_permutation_test_seeded_reproducibility()
        integer(int32), parameter :: n_bins = 3, n_points = 2, n_studies = 2, n_permutations = 20
        integer(int32) :: mean_pmf_counts(n_bins, n_points), mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points), global_jsd_observed(n_studies)
        real(real64) :: p_values_first(n_studies), p_values_second(n_studies)
        integer(int32) :: ierr

        mean_pmf_counts(:, 1) = [5, 3, 2]
        mean_pmf_counts(:, 2) = [2, 2, 6]
        mean_pmf_included_n_reps = [10, 10]
        included_n_reps(:, 1) = [4, 4]
        included_n_reps(:, 2) = [6, 6]
        mean_pmf(:, 1) = real(mean_pmf_counts(:, 1), real64)/10.0_real64
        mean_pmf(:, 2) = real(mean_pmf_counts(:, 2), real64)/10.0_real64
        global_jsd_observed = 0.3_real64

        call gjct_permutation_test(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                   mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values_first, &
                                   ierr=ierr, random_seed=123_int32)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_gjct_permutation_test_seeded_reproducibility: first call ierr should be OK")

        call gjct_permutation_test(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                   mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values_second, &
                                   ierr=ierr, random_seed=123_int32)
        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_gjct_permutation_test_seeded_reproducibility: second call ierr should be OK")

        call assert_equal_array_real(p_values_first, p_values_second, size(p_values_first, kind=int32), TOL, &
                                     "test_gjct_permutation_test_seeded_reproducibility: "// &
                                     "same random_seed twice must give identical p_values")
    end subroutine test_gjct_permutation_test_seeded_reproducibility

    !> Documents the deferred limitation from the plan: `gjct_permutation_test_impl` ports
    !| origin/125-stabilize-jscomp's p-value formula exactly, WITHOUT the `(1+count)/(n+1)` Laplace
    !| correction, so `p` can come out exactly `0.0` whenever no permutation's resampled JSD
    !| reaches the observed value. Contrived, deterministically (not merely probabilistically) so:
    !| a single-bin histogram (`n_bins=1`) means every resample's pmf is `[1.0]`, for both the
    !| study and the (unperturbed) consensus mean alike, so the Jensen-Shannon divergence between
    !| them is exactly `0.0` on *every* permutation, regardless of which counts
    !| `random_multiv_hypergeom` happens to draw. Observing a strictly positive
    !| `global_jsd_observed` then guarantees no permutation's JSD can reach it, so
    !| `p_values = anint(0)/n_permutations = 0.0` exactly -- not the `1/(n_permutations+1)` the
    !| Laplace-corrected formula would give.
    subroutine test_permutation_pvalue_can_be_exactly_zero_known_limitation()
        integer(int32), parameter :: n_bins = 1, n_points = 1, n_studies = 1, n_permutations = 10
        integer(int32) :: mean_pmf_counts(n_bins, n_points), mean_pmf_included_n_reps(n_points)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points), global_jsd_observed(n_studies), p_values(n_studies)
        integer(int32) :: ierr

        mean_pmf_counts(1, 1) = 10_int32
        mean_pmf_included_n_reps = [10_int32]
        included_n_reps(1, :) = [5_int32]
        mean_pmf(1, 1) = 1.0_real64
        global_jsd_observed = 0.5_real64 ! strictly greater than the always-0.0 resampled JSD

        call gjct_permutation_test(n_permutations, n_bins, n_points, n_studies, mean_pmf_counts, mean_pmf, &
                                   mean_pmf_included_n_reps, included_n_reps, global_jsd_observed, p_values, &
                                   ierr=ierr, random_seed=7_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_permutation_pvalue_can_be_exactly_zero_known_limitation: ierr should be OK")
        call assert_equal_real(p_values(1), 0.0_real64, TOL, &
                               "test_permutation_pvalue_can_be_exactly_zero_known_limitation: "// &
                               "known limitation: p must be exactly 0.0, not (1+0)/(n_permutations+1)")
    end subroutine test_permutation_pvalue_can_be_exactly_zero_known_limitation

    !> End-to-end, fully closed-form 2-study case run through `run_js_comp_test` with
    !| `n_permutations=0`: GSL's `create_rng` still runs (so `ierr` is genuinely exercised), but the
    !| permutation loop itself never executes, so `gjct_permutation_test_impl` leaves everything
    !| untouched and `run_js_comp_test_impl`'s final re-derivation step reproduces exactly the same
    !| values it already had -- making the whole pipeline deterministic and hand-traceable.
    !|
    !| One reference point, one neighbor per study (`x_star=1.0` sits exactly on gene 1's mean, so
    !| both studies trivially pick gene 1 -- see test_construct_neighborhoods_ranged_basic's own
    !| binary-search trace for why an exact match gives insertion index 1). Study 1's two replicate
    !| residuals are `[-1, 1]`, study 2's are `[-3, 3]`, binned into 4 bins spanning `[-4, 4]` (bin
    !| width 2.0, via `build_residual_histograms_impl`'s own
    !| `bin_idx = min(n_bins, int((clamped+R)/bin_width)+1)`): study 1 -> pmf `[0, 0.5, 0.5, 0]`,
    !| study 2 -> pmf `[0.5, 0, 0, 0.5]`, mean pmf -> uniform `[0.25, 0.25, 0.25, 0.25]`.
    !|
    !| Study 2's pmf is study 1's own pmf under the bin permutation `1<->4, 2<->3`, which also fixes
    !| the uniform mean pmf -- so by that symmetry both studies' JSD against the mean must be
    !| identical. Working the Jensen-Shannon sum out by hand (`compute_divergence_per_reference_point_impl`'s
    !| own formula, `0.5*sum(s1*log(s1/S_mean) + s2*log(s2/S_mean))`, then rescaled by `/LOG_2`)
    !| collapses to the closed form `1.5 - 0.75*log2(3)`: bins 1 and 4 each contribute
    !| `0.25*ln(2)`, bins 2 and 3 each contribute `0.5*ln(4/3) + 0.25*ln(2/3)`, so the raw (pre-rescale)
    !| sum is `0.5*ln2 + ln(4/3) + 0.5*ln(2/3) = 3*ln2 - 1.5*ln3`; halving (the JSD's own `0.5*`) and
    !| dividing by `ln2` gives `1.5 - 0.75*(ln3/ln2) = 1.5 - 0.75*log2(3)`. With `n_points=1` the
    !| single reference point's weight is always exactly 1.0 (its own plus the consensus's included
    !| reps divide out exactly), so `global_js_divergence` equals that same closed form too.
    subroutine test_run_js_comp_test_two_studies_hand_traceable()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 2, max_n_reps_all_studies = 2
        integer(int32), parameter :: n_points = 1, n_neighbors = 1, n_bins = 4
        real(real64), parameter :: LOG2_3 = 1.5849625007211562_real64 ! log2(3) = ln(3)/ln(2)
        real(real64), parameter :: EXPECTED_JSD = 1.5_real64 - 0.75_real64*LOG2_3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        integer(int32) :: gene_means_perms(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        real(real64) :: x_star(n_points)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points, n_studies)
        integer(int32) :: neighborhood_range(2, n_points, n_studies)
        real(real64) :: pmfs(n_bins, n_points, n_studies)
        integer(int32) :: counts(n_bins, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points)
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: mean_pmf_included_n_reps(n_points)
        real(real64) :: js_divergences(n_points, n_studies), weights(n_points, n_studies)
        real(real64) :: global_js_divergence(n_studies), p_values(n_studies)
        integer(int32) :: ierr
        real(real64) :: nan_val

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)

        gene_means(:, 1) = [1.0_real64, 5.0_real64]
        gene_means(:, 2) = [1.0_real64, 5.0_real64]
        gene_means_perms(:, 1) = [1, 2]
        gene_means_perms(:, 2) = [1, 2]

        residuals(:, 1, 1) = [-1.0_real64, 1.0_real64]
        residuals(:, 2, 1) = nan_val
        residuals(:, 1, 2) = [-3.0_real64, 3.0_real64]
        residuals(:, 2, 2) = nan_val

        x_star = [1.0_real64]

        call run_js_comp_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, n_bins, &
                              4.0_real64, gene_means, gene_means_perms, residuals, x_star, neighborhood_indices, &
                              neighborhood_range, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, p_values, &
                              ierr=ierr, n_permutations=0_int32, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_run_js_comp_test_two_studies_hand_traceable: ierr should be OK")

        call assert_equal_array_int(neighborhood_indices(:, 1, 1), [1], n_neighbors, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 1 neighbor is gene 1")
        call assert_equal_array_int(neighborhood_indices(:, 1, 2), [1], n_neighbors, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 2 neighbor is gene 1")

        call assert_equal_array_int(counts(:, 1, 1), [0, 1, 1, 0], n_bins, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 1 counts")
        call assert_equal_array_int(counts(:, 1, 2), [1, 0, 0, 1], n_bins, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 2 counts")
        call assert_equal_int(included_n_reps(1, 1), 2_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: study 1 included_n_reps")
        call assert_equal_int(included_n_reps(1, 2), 2_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: study 2 included_n_reps")

        call assert_equal_array_real(mean_pmf(:, 1), [0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64], n_bins, TOL, &
                                     "test_run_js_comp_test_two_studies_hand_traceable: mean_pmf should be uniform")
        call assert_equal_array_int(mean_pmf_counts(:, 1), [1, 1, 1, 1], n_bins, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: mean_pmf_counts")
        call assert_equal_int(mean_pmf_included_n_reps(1), 4_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: mean_pmf_included_n_reps")

        call assert_equal_real(global_js_divergence(1), EXPECTED_JSD, 1d-9, &
                               "test_run_js_comp_test_two_studies_hand_traceable: study 1 global JSD, closed form")
        call assert_equal_real(global_js_divergence(2), EXPECTED_JSD, 1d-9, &
                               "test_run_js_comp_test_two_studies_hand_traceable: study 2 global JSD, closed form")
        call assert_equal_real(weights(1, 1), 1.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: single reference point weight is 1.0")
        call assert_equal_real(weights(1, 2), 1.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: single reference point weight is 1.0")

        call assert_equal_real(p_values(1), 0.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: n_permutations=0 -> p_values stay 0.0")
        call assert_equal_real(p_values(2), 0.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: n_permutations=0 -> p_values stay 0.0")
    end subroutine test_run_js_comp_test_two_studies_hand_traceable

    !> Three studies, same single-point/single-neighbor topology as the hand-traceable case above.
    !| Studies 1 and 2 are identical to each other (residuals `[-1,-1,1,1]`, pmf `[0, 0.5, 0.5, 0]`);
    !| study 3 is a deliberately constructed outlier whose every replicate lands in the SAME bin
    !| (`[-3.9,-3.9,-3.9,-3.9]`, pmf `[1, 0, 0, 0]`). The consensus mean pmf, averaging all three, is
    !| `[1/3, 1/3, 1/3, 0]` with `mean_pmf_counts=[4, 4, 4, 0]`. Study 3's own pmf puts everything in
    !| a bin that holds only 4 of the pooled pool's 12 replicates, so drawing (without replacement,
    !| `random_multiv_hypergeom`) another 4-for-4 landing entirely in that one bin purely by chance
    !| is exceedingly rare -- its empirical p-value must come out small, while studies 1/2's own
    !| draws, being close to what the consensus itself is built from, should not be nearly as
    !| extreme. A fixed random_seed makes the outcome fully deterministic, so this is not a flaky
    !| probabilistic assertion.
    subroutine test_run_js_comp_test_three_studies_outlier_has_small_p_value()
        integer(int32), parameter :: n_studies = 3, max_n_genes_all_studies = 2, max_n_reps_all_studies = 4
        integer(int32), parameter :: n_points = 1, n_neighbors = 1, n_bins = 4, n_permutations = 500
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        integer(int32) :: gene_means_perms(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        real(real64) :: x_star(n_points)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points, n_studies)
        integer(int32) :: neighborhood_range(2, n_points, n_studies)
        real(real64) :: pmfs(n_bins, n_points, n_studies)
        integer(int32) :: counts(n_bins, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(n_bins, n_points)
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: mean_pmf_included_n_reps(n_points)
        real(real64) :: js_divergences(n_points, n_studies), weights(n_points, n_studies)
        real(real64) :: global_js_divergence(n_studies), p_values(n_studies)
        integer(int32) :: ierr, i_study
        real(real64) :: nan_val

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)

        do i_study = 1, n_studies
            gene_means(:, i_study) = [1.0_real64, 5.0_real64]
            gene_means_perms(:, i_study) = [1, 2]
            residuals(:, 2, i_study) = nan_val
        end do
        residuals(:, 1, 1) = [-1.0_real64, -1.0_real64, 1.0_real64, 1.0_real64]
        residuals(:, 1, 2) = [-1.0_real64, -1.0_real64, 1.0_real64, 1.0_real64]
        residuals(:, 1, 3) = [-3.9_real64, -3.9_real64, -3.9_real64, -3.9_real64]

        x_star = [1.0_real64]

        call run_js_comp_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, n_bins, &
                              4.0_real64, gene_means, gene_means_perms, residuals, x_star, neighborhood_indices, &
                              neighborhood_range, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, p_values, &
                              ierr=ierr, n_permutations=n_permutations, random_seed=42_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_run_js_comp_test_three_studies_outlier_has_small_p_value: ierr should be OK")

        call assert_equal_array_int(mean_pmf_counts(:, 1), [4, 4, 4, 0], n_bins, &
                                    "test_run_js_comp_test_three_studies_outlier_has_small_p_value: mean_pmf_counts")

        call assert_true(p_values(3) <= 0.05_real64, &
                         "test_run_js_comp_test_three_studies_outlier_has_small_p_value: "// &
                         "the outlier study's p-value should be small")
        call assert_true(p_values(3) < p_values(1), &
                         "test_run_js_comp_test_three_studies_outlier_has_small_p_value: "// &
                         "the outlier's p-value should be smaller than study 1's")
        call assert_true(p_values(3) < p_values(2), &
                         "test_run_js_comp_test_three_studies_outlier_has_small_p_value: "// &
                         "the outlier's p-value should be smaller than study 2's")
    end subroutine test_run_js_comp_test_three_studies_outlier_has_small_p_value

    !> Forces `run_js_comp_test_parameter_search`'s second admissibility gate
    !| (`min_count_per_mean_bin`) impossibly high, so no candidate in the grid can ever pass it and
    !| `check_plateau_condition` is never even called: `plateau_found` stays `.false.` for the whole
    !| search. With `max_n_genes_all_studies=2000` the GAMMA-decay grid produces two candidates that
    !| both share `n_points=300` (`ceil(4*sqrt(2000))=179`, clamped up to `MIN_POINTS=300`, and a
    !| single GAMMA step already drops below `n_points_low=300`) with `n_neighbors=26` and `13`
    !| (`floor(2000/(0.25*300))=26`, `max(1,floor(2000/(0.5*300)))=13`) -- both values derived purely
    !| from the GAMMA-decay constants and `max_n_genes_all_studies`, independent of the synthetic
    !| gene/residual data below. Since no candidate ever plateaus, the search must fall back to the
    !| FIRST (finest resolution) candidate, `(n_points, n_neighbors) = (300, 26)`, and reset
    !| `best_candidate_pair_confidence_interval` to `-1.0` throughout.
    subroutine test_param_search_no_plateau_falls_back_to_finest()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 2000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = real(i_gene, real64)
                residuals(:, i_gene, i_study) = [-0.5_real64, 0.0_real64, 0.5_real64]
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 1.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=1000000_int32, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_no_plateau_falls_back_to_finest: ierr should be OK")
        call assert_false(plateau_established, &
                          "test_param_search_no_plateau_falls_back_to_finest: plateau_established is false")
        call assert_equal_int(n_points, 300_int32, &
                              "test_param_search_no_plateau_falls_back_to_finest: "// &
                              "falls back to the finest-resolution n_points")
        call assert_equal_int(n_neighbors, 26_int32, &
                              "test_param_search_no_plateau_falls_back_to_finest: "// &
                              "falls back to the finest-resolution n_neighbors")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_falls_back_to_finest: "// &
                                     "study 1 CI reset to -1.0")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_falls_back_to_finest: "// &
                                     "study 2 CI reset to -1.0")
    end subroutine test_param_search_no_plateau_falls_back_to_finest

    !> With `max_n_genes_all_studies=100`, the GAMMA-decay grid collapses to exactly ONE candidate
    !| (`floor(100/(0.25*300))=1` and `max(1,floor(100/(0.5*300)))=max(1,0)=1` are the same value, so
    !| the second `KX_FACTORS` entry does not add a distinct candidate, and a single GAMMA step
    !| already drops `n_points_high` below `n_points_low=300`): `(n_points, n_neighbors) = (300, 1)`. Both
    !| admissibility gates are relaxed to their most permissive settings
    !| (`min_neighbor_overlap=0.0`, `min_count_per_mean_bin=0`) so the sole candidate exercises the
    !| full pipeline -- both gates pass, `bootstrap_histogram` and `check_plateau_condition` really
    !| do run -- yet the search must still return that one candidate regardless of whatever
    !| `plateau_found` comes out as, because `n_candidates < 2` bypasses the plateau machinery
    !| entirely (the same fallback branch `run_js_comp_test_parameter_search_no_plateau_falls_back_to_finest`
    !| reaches by forcing every candidate to fail a gate instead of there only being one to begin
    !| with).
    subroutine test_param_search_single_candidate_bypasses_plateau()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 100, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = real(i_gene, real64)
                residuals(:, i_gene, i_study) = [-0.5_real64, 0.0_real64, 0.5_real64]
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 1.0_real64, 5_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_single_candidate_bypasses_plateau: ierr should be OK")
        call assert_true(plateau_established, &
                         "test_param_search_single_candidate_bypasses_plateau: plateau_established is true")
        call assert_equal_int(n_points, 300_int32, &
                              "test_param_search_single_candidate_bypasses_plateau: "// &
                              "the sole candidate's n_points")
        call assert_equal_int(n_neighbors, 1_int32, &
                              "test_param_search_single_candidate_bypasses_plateau: "// &
                              "the sole candidate's n_neighbors")
        call assert_true(best_candidate_pair_confidence_interval(1, 1) >= 0.0_real64 .and. &
                         best_candidate_pair_confidence_interval(2, 1) <= 1.0_real64, &
                         "test_param_search_single_candidate_bypasses_plateau: "// &
                         "relaxed gates let bootstrap actually run, giving a real (not -1.0) CI")
    end subroutine test_param_search_single_candidate_bypasses_plateau

    !> Exercises the search's actual "success" path: stopping at a genuine mid-grid plateau, rather
    !| than either the no-plateau fallback (`test_param_search_no_plateau_falls_back_to_finest`) or
    !| the single-candidate bypass (`test_param_search_single_candidate_bypasses_plateau`) above.
    !|
    !| With `max_n_genes_all_studies=10000`, the GAMMA-decay grid produces exactly FOUR candidates:
    !| `n_points_high = clamp(ceil(4*sqrt(10000)), 300, 1500) = 400`,
    !| `n_points_low = max(300, ceil(0.2*400)) = 300`.
    !| - i=1: `n_points=400`, neighbors `floor(10000/(0.25*400))=100` and
    !|   `max(1,floor(10000/(0.5*400)))=50` -> candidates (400,100), (400,50). `n_points_high *= 0.8
    !|   -> 320`.
    !| - i=2: `320 >= 300`, so `n_points=320`, neighbors `floor(10000/(0.25*320))=125` and
    !|   `max(1,floor(10000/(0.5*320)))=62` -> candidates (320,125), (320,62). `n_points_high *= 0.8
    !|   -> 256`.
    !| - i=3: `256 < 300` -> loop exits.
    !| So the grid is `[(400,100), (400,50), (320,125), (320,62)]`.
    !|
    !| Every gene's residual and gene mean is the SAME constant across every study (0.0 and 5.0
    !| respectively). This deterministically produces a plateau at the SECOND candidate, regardless
    !| of `join_method` or `random_seed`:
    !| - Constant gene means make every reference point's neighborhood range collapse to the tie-
    !|   extended full gene range for every candidate, so consecutive neighborhoods always overlap
    !|   by exactly 1.0 -- the first admissibility gate always passes.
    !| - Constant residuals mean every neighborhood's histogram places all mass in the SAME single
    !|   bin (`build_residual_histograms_impl` bins purely on the residual value, independent of
    !|   `n_neighbors`), so every study's pmf and the consensus (mean) pmf are identical delta
    !|   distributions -> the observed JSD is exactly 0.0 for every candidate, and every bootstrap
    !|   resample draws from that same single-bin pool, so no resampled value can ever exceed the
    !|   seeded 0.0 in either direction -> the bootstrapped CI is exactly `[0.0, 0.0]`, for every
    !|   candidate, deterministically (not dependent on random_seed).
    !| - `check_plateau_condition_impl`'s `compute_fractional_overlap` scores a candidate whose CI
    !|   sits inside the running best as "no worse"; comparing `[0.0, 0.0]` against a first "best" of
    !|   `[-1.0, -1.0]` (candidate 1) gives overlap 0.0 (no overlap at all -- fails every join
    !|   method), but comparing `[0.0, 0.0]` against a running best that is ITSELF `[0.0, 0.0]`
    !|   (candidate 2, once candidate 1 became the new "best") gives overlap 1.0 for every study,
    !|   which exceeds the default `succeeding_ci_overlap=0.9` for every study at once -- so
    !|   METHOD_JOIN_MIN/MAX/MEDIAN would all detect the plateau there. The search must therefore
    !|   stop at candidate 2 = `(n_points, n_neighbors) = (400, 50)`, never reaching candidates 3/4,
    !|   and must NOT fall back to candidate 1 or reset the CI to the `-1.0` sentinel.
    !|
    !| Gate thresholds are relaxed to the same permissive settings
    !| `test_param_search_single_candidate_bypasses_plateau` uses to let the gates genuinely pass
    !| (`min_count_per_mean_bin=0`, `min_neighbor_overlap=0.0`) -- not the impossible values
    !| `test_param_search_no_plateau_falls_back_to_finest` uses to deliberately force every
    !| candidate to fail a gate.
    subroutine test_param_search_finds_plateau_mid_grid()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 10000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = 5.0_real64
                residuals(:, i_gene, i_study) = 0.0_real64
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 1.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_finds_plateau_mid_grid: ierr should be OK")
        call assert_true(plateau_established, &
                         "test_param_search_finds_plateau_mid_grid: plateau_established is true")

        ! Both candidates the search actually reached (1 and 2) passed both admissibility gates,
        ! so n_admissible_evaluated == 2 -- candidates 3/4 have no trace entry at all, not a
        ! zero-filled one, since the search stopped before reaching them.
        call assert_equal_int(n_admissible_evaluated, 2_int32, &
                              "test_param_search_finds_plateau_mid_grid: two admissible candidates evaluated")
        call assert_equal_int(trace_n_points(1), 400_int32, &
                              "test_param_search_finds_plateau_mid_grid: trace_n_points(1)")
        call assert_equal_int(trace_n_neighbors(1), 100_int32, &
                              "test_param_search_finds_plateau_mid_grid: trace_n_neighbors(1)")
        call assert_equal_int(trace_n_points(2), 400_int32, &
                              "test_param_search_finds_plateau_mid_grid: trace_n_points(2)")
        call assert_equal_int(trace_n_neighbors(2), 50_int32, &
                              "test_param_search_finds_plateau_mid_grid: trace_n_neighbors(2)")
        call assert_equal_array_real(trace_global_js_divergence(:, 1), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_global_js_divergence(:,1) is 0.0")
        call assert_equal_array_real(trace_global_js_divergence(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_global_js_divergence(:,2) is 0.0")
        call assert_equal_array_real(trace_ci_lower(:, 1), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_ci_lower(:,1) is 0.0")
        call assert_equal_array_real(trace_ci_upper(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_ci_upper(:,2) is 0.0")
        call assert_equal_array_real(trace_ci_width(:, 1), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_ci_width(:,1) is 0.0")
        call assert_equal_array_real(trace_ci_width(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_ci_width(:,2) is 0.0")
        call assert_equal_array_real(trace_ci_width_relative(:, 1), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_ci_width_relative(:,1) is 0.0")
        call assert_equal_array_real(trace_ci_width_relative(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_ci_width_relative(:,2) is 0.0")

        ! Slot 1 is the first admissible candidate -- no predecessor to diff against, so the
        ! effect-size trace is the -1.0 sentinel throughout, never a real (mis-)computed value.
        call assert_equal_array_real(trace_delta(:, 1), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_delta(:,1) is the -1.0 sentinel")
        call assert_equal_real(trace_delta_median(1), -1.0_real64, TOL, &
                              "test_param_search_finds_plateau_mid_grid: trace_delta_median(1) is the -1.0 sentinel")
        call assert_equal_real(trace_delta_max(1), -1.0_real64, TOL, &
                              "test_param_search_finds_plateau_mid_grid: trace_delta_max(1) is the -1.0 sentinel")

        ! Slot 2 is a real transition (candidate 1 -> candidate 2), both JSD == 0.0, so delta == 0.0
        ! for both studies -- a genuine computed value, not the sentinel.
        call assert_equal_array_real(trace_delta(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: trace_delta(:,2) is 0.0")
        call assert_equal_real(trace_delta_median(2), 0.0_real64, TOL, &
                              "test_param_search_finds_plateau_mid_grid: trace_delta_median(2) is 0.0")
        call assert_equal_real(trace_delta_max(2), 0.0_real64, TOL, &
                              "test_param_search_finds_plateau_mid_grid: trace_delta_max(2) is 0.0")

        ! The grid's SECOND candidate, (400, 50) -- not the first (finest, 400/100) and not the
        ! no-plateau fallback's candidate either -- proving the search stopped early at a genuine
        ! plateau rather than running to completion or falling back.
        call assert_equal_int(n_points, 400_int32, &
                              "test_param_search_finds_plateau_mid_grid: "// &
                              "should stop at the second candidate's n_points")
        call assert_equal_int(n_neighbors, 50_int32, &
                              "test_param_search_finds_plateau_mid_grid: "// &
                              "should stop at the second candidate's n_neighbors")

        ! A real, bootstrapped [0.0, 0.0] CI -- not the -1.0 sentinel the no-plateau fallback resets
        ! to -- for both studies.
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: "// &
                                     "study 1 CI should be exactly [0.0, 0.0]")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_finds_plateau_mid_grid: "// &
                                     "study 2 CI should be exactly [0.0, 0.0]")
    end subroutine test_param_search_finds_plateau_mid_grid

    !> Same degenerate data as `test_param_search_finds_plateau_mid_grid` (constant gene means and
    !| residuals -> JSD exactly 0.0 and CI exactly [0.0, 0.0] for every candidate in the 4-candidate
    !| grid `[(400,100), (400,50), (320,125), (320,62)]`), but with `plateau_mode=MODE_PLATEAU_EFFECT_SIZE`.
    !| Under the default `plateau_mode` (CI overlap), that test's search stops at the SECOND
    !| candidate. Under effect size alone, `delta_min_consecutive_transitions=2` (the default) needs
    !| TWO consecutive qualifying transitions before it plateaus: the first (candidate 1 -> 2) is
    !| only the first, so the search must continue past candidate 2 -- proving CI overlap is
    !| genuinely ignored under this mode, not just usually also satisfied -- and stop only once the
    !| second consecutive qualifying transition (candidate 2 -> 3) completes, at the THIRD candidate,
    !| `(320, 125)`.
    subroutine test_param_search_effect_size_mode_plateau()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 10000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = 5.0_real64
                residuals(:, i_gene, i_study) = 0.0_real64
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 1.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               plateau_mode=MODE_PLATEAU_EFFECT_SIZE, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_effect_size_mode_plateau: ierr should be OK")
        call assert_true(plateau_established, &
                         "test_param_search_effect_size_mode_plateau: plateau_established is true")
        call assert_equal_int(n_admissible_evaluated, 3_int32, &
                              "test_param_search_effect_size_mode_plateau: "// &
                              "needs two consecutive qualifying transitions -> three candidates evaluated")
        call assert_equal_int(n_points, 320_int32, &
                              "test_param_search_effect_size_mode_plateau: "// &
                              "stops at the THIRD candidate's n_points, later than CI-overlap mode's second")
        call assert_equal_int(n_neighbors, 125_int32, &
                              "test_param_search_effect_size_mode_plateau: "// &
                              "stops at the third candidate's n_neighbors")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_effect_size_mode_plateau: "// &
                                     "study 1 CI is the triggering candidate's own [0.0, 0.0], from the override")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), [0.0_real64, 0.0_real64], 2_int32, TOL, &
                                     "test_param_search_effect_size_mode_plateau: "// &
                                     "study 2 CI is the triggering candidate's own [0.0, 0.0], from the override")
    end subroutine test_param_search_effect_size_mode_plateau

    !> Same degenerate data and grid as the two tests above. Under `plateau_mode=MODE_PLATEAU_BOTH`,
    !| the search must stop at whichever criterion plateaus FIRST -- here that is CI overlap at the
    !| SECOND candidate (exactly `test_param_search_finds_plateau_mid_grid`'s result), not effect
    !| size's third-candidate result from `test_param_search_effect_size_mode_plateau` above. This is
    !| the one case that actually exercises the `select case` `MODE_PLATEAU_BOTH` branch: with only
    !| one criterion ever selected in the other two tests, BOTH is the only mode where the two
    !| criteria's plateau points genuinely differ and the earlier one must win.
    subroutine test_param_search_both_mode_uses_earlier_trigger()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 10000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = 5.0_real64
                residuals(:, i_gene, i_study) = 0.0_real64
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 1.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               plateau_mode=MODE_PLATEAU_BOTH, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_both_mode_uses_earlier_trigger: ierr should be OK")
        call assert_true(plateau_established, &
                         "test_param_search_both_mode_uses_earlier_trigger: plateau_established is true")
        call assert_equal_int(n_admissible_evaluated, 2_int32, &
                              "test_param_search_both_mode_uses_earlier_trigger: "// &
                              "CI overlap's earlier plateau wins -> only two candidates evaluated")
        call assert_equal_int(n_points, 400_int32, &
                              "test_param_search_both_mode_uses_earlier_trigger: "// &
                              "stops at the second candidate's n_points, same as CI-overlap-only mode")
        call assert_equal_int(n_neighbors, 50_int32, &
                              "test_param_search_both_mode_uses_earlier_trigger: "// &
                              "stops at the second candidate's n_neighbors")
    end subroutine test_param_search_both_mode_uses_earlier_trigger

    !> Exercises the new Issue #178 no-plateau fallback (`plateau_established = .false.` +
    !| smallest-bootstrap-uncertainty candidate selection), which only the tests above never
    !| reach: they all either plateau or bypass the plateau machinery entirely via a
    !| single/zero-admissible-candidate grid. This one needs real, non-degenerate data -- the
    !| all-constant fixture the other tests use gives an exactly-[0,0] CI for every candidate,
    !| which trivially plateaus (perfect overlap) rather than ever exercising this fallback.
    !|
    !| `max_n_genes_all_studies=20000` gives a 6-candidate grid (`[(566,141),(566,70),(453,176),
    !| (453,88),(362,220),(362,110)]`). Both studies' residuals are independent draws from the
    !| same simple linear-congruential pseudo-random sequence (`x_{n+1} = (1103515245*n + 12345)
    !| mod 2^31`, scaled to `[-2, 2]`) -- no systematic per-study difference, so the observed JSD
    !| stays small and the bootstrapped CIs genuinely vary in width from one candidate to the next
    !| without any one of them ever nesting inside the running best closely enough to satisfy
    !| `succeeding_ci_overlap`'s default 90% threshold. This was found empirically (per the
    !| project's own plan for this change, which flagged the exact fixture as needing
    !| construction/iteration, not something derivable on paper) -- confirmed via direct
    !| experimentation that every candidate here is genuinely evaluated, no plateau is ever found
    !| under any of the three `join_method`s, and candidate 5 = `(362, 220)` has a distinctly
    !| smaller confidence-interval width than every other candidate, in both studies (median CI
    !| width per candidate: `0.00167, 0.00250, 0.00121, 0.00218, 0.00110, 0.00200` -- candidate 5's
    !| `0.00110` is the unambiguous minimum). The routine must therefore return candidate 5's own
    !| `(n_points, n_neighbors)` and its real, bootstrapped confidence interval (not `-1.0`), with
    !| `plateau_established = .false.` The exact confidence-interval values below come directly
    !| from running this fixture through the actual implementation (not hand-derived), since they
    !| depend on the real bootstrap resampling -- reproducible bit-for-bit given the fixed
    !| `random_seed` and deterministic input data.
    subroutine test_param_search_no_plateau_uses_smallest_uncertainty()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study, i_rep, k

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = real(i_gene, real64)
                do i_rep = 1, max_n_reps_all_studies
                    k = i_rep + max_n_reps_all_studies*(i_gene - 1) + &
                        (i_study - 1)*max_n_reps_all_studies*max_n_genes_all_studies
                    residuals(i_rep, i_gene, i_study) = 4.0_real64* &
                        (mod(1103515245.0_real64*real(k, real64) + 12345.0_real64, 2147483648.0_real64) &
                         /2147483648.0_real64 - 0.5_real64)
                end do
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 3.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_no_plateau_uses_smallest_uncertainty: ierr should be OK")
        call assert_equal_int(n_admissible_evaluated, 6_int32, &
                              "test_param_search_no_plateau_uses_smallest_uncertainty: "// &
                              "all six candidates were admissible and evaluated")
        call assert_false(plateau_established, &
                          "test_param_search_no_plateau_uses_smallest_uncertainty: "// &
                          "no candidate ever plateaus under this fixture")
        call assert_equal_int(n_points, 362_int32, &
                              "test_param_search_no_plateau_uses_smallest_uncertainty: "// &
                              "smallest-uncertainty candidate's n_points")
        call assert_equal_int(n_neighbors, 220_int32, &
                              "test_param_search_no_plateau_uses_smallest_uncertainty: "// &
                              "smallest-uncertainty candidate's n_neighbors")
        call assert_true(all(best_candidate_pair_confidence_interval /= -1.0_real64), &
                         "test_param_search_no_plateau_uses_smallest_uncertainty: "// &
                         "a real confidence interval is returned, not the -1.0 sentinel")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), &
                                     [5.212653804e-05_real64, 1.149283224e-03_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_uses_smallest_uncertainty: study 1 CI")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), &
                                     [5.212150674e-05_real64, 1.150810589e-03_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_uses_smallest_uncertainty: study 2 CI")
    end subroutine test_param_search_no_plateau_uses_smallest_uncertainty

    !> Same LCG-residual fixture as `test_param_search_no_plateau_uses_smallest_uncertainty`
    !| (guaranteed no CI-overlap plateau across all 6 candidates), but with
    !| `plateau_mode=MODE_PLATEAU_EFFECT_SIZE` -- verifying the smallest-bootstrap-uncertainty
    !| fallback is genuinely gated to `MODE_PLATEAU_CI_OVERLAP` (per the Issue #178 fallback plan's
    !| explicit scope decision), not applied under this mode. If effect size also never plateaus on
    !| this fixture, the routine must fall back to the OLD behavior: candidate 1 (finest,
    !| `(566, 141)`), CI reset to `-1.0`, `plateau_established = .false.` -- not candidate 5's
    !| smallest-uncertainty result.
    subroutine test_param_search_no_plateau_effect_size_falls_back()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study, i_rep, k

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = real(i_gene, real64)
                do i_rep = 1, max_n_reps_all_studies
                    k = i_rep + max_n_reps_all_studies*(i_gene - 1) + &
                        (i_study - 1)*max_n_reps_all_studies*max_n_genes_all_studies
                    residuals(i_rep, i_gene, i_study) = 4.0_real64* &
                        (mod(1103515245.0_real64*real(k, real64) + 12345.0_real64, 2147483648.0_real64) &
                         /2147483648.0_real64 - 0.5_real64)
                end do
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 3.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               plateau_mode=MODE_PLATEAU_EFFECT_SIZE, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_no_plateau_effect_size_falls_back: ierr should be OK")
        call assert_false(plateau_established, &
                          "test_param_search_no_plateau_effect_size_falls_back: "// &
                          "no candidate ever plateaus under effect-size mode on this fixture")
        call assert_equal_int(n_points, 566_int32, &
                              "test_param_search_no_plateau_effect_size_falls_back: "// &
                              "falls back to the finest-resolution n_points, NOT candidate 5's smallest-uncertainty n_points")
        call assert_equal_int(n_neighbors, 141_int32, &
                              "test_param_search_no_plateau_effect_size_falls_back: "// &
                              "falls back to the finest-resolution n_neighbors")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_effect_size_falls_back: "// &
                                     "study 1 CI reset to -1.0, NOT a real bootstrapped CI")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_effect_size_falls_back: "// &
                                     "study 2 CI reset to -1.0")
    end subroutine test_param_search_no_plateau_effect_size_falls_back

    !> Same fixture again, `plateau_mode=MODE_PLATEAU_BOTH`. Since CI overlap never plateaus here
    !| (confirmed by `test_param_search_no_plateau_uses_smallest_uncertainty`) and effect size never
    !| plateaus here either (confirmed by
    !| `test_param_search_no_plateau_effect_size_falls_back`), BOTH mode's
    !| `ci_plateau_found .or. effect_size_plateau_found` is false for every candidate too -- so this
    !| must land on the exact same old-fallback branch (candidate 1, CI reset to -1.0), not the new
    !| smallest-uncertainty selection, confirming the gate checks `plateau_mode ==
    !| MODE_PLATEAU_CI_OVERLAP` specifically rather than merely `/= MODE_PLATEAU_EFFECT_SIZE`.
    subroutine test_param_search_no_plateau_both_falls_back()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, n_bins, ierr, n_admissible_evaluated
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study, i_rep, k

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = real(i_gene, real64)
                do i_rep = 1, max_n_reps_all_studies
                    k = i_rep + max_n_reps_all_studies*(i_gene - 1) + &
                        (i_study - 1)*max_n_reps_all_studies*max_n_genes_all_studies
                    residuals(i_rep, i_gene, i_study) = 4.0_real64* &
                        (mod(1103515245.0_real64*real(k, real64) + 12345.0_real64, 2147483648.0_real64) &
                         /2147483648.0_real64 - 0.5_real64)
                end do
            end do
        end do

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 3.0_real64, 10_int32, METHOD_JOIN_MIN, n_points, n_neighbors, &
                                               n_bins, best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, ierr=ierr, &
                                               min_count_per_mean_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               plateau_mode=MODE_PLATEAU_BOTH, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_no_plateau_both_falls_back: ierr should be OK")
        call assert_false(plateau_established, &
                          "test_param_search_no_plateau_both_falls_back: "// &
                          "neither criterion ever plateaus under this fixture")
        call assert_equal_int(n_points, 566_int32, &
                              "test_param_search_no_plateau_both_falls_back: "// &
                              "falls back to the finest-resolution n_points, NOT candidate 5's smallest-uncertainty n_points")
        call assert_equal_int(n_neighbors, 141_int32, &
                              "test_param_search_no_plateau_both_falls_back: "// &
                              "falls back to the finest-resolution n_neighbors")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_both_falls_back: "// &
                                     "study 1 CI reset to -1.0, NOT a real bootstrapped CI")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_both_falls_back: "// &
                                     "study 2 CI reset to -1.0")
    end subroutine test_param_search_no_plateau_both_falls_back

end module mod_test_data_integration_js_comp_test
