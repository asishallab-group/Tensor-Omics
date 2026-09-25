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
    use tox_data_integration_js_comp_test, only: estimate_bin_count, determine_bin_count_occupancy, &
                                                  determine_bin_count_occupancy_exhaustive, &
                                                  generate_js_comp_test_candidates, check_neighborhood_overlaps, &
                                                  check_mean_pmf_min_counts, check_plateau_condition, &
                                                  check_effect_size_plateau_condition, create_mean_pmf, &
                                                  create_mean_pmf_only, bootstrap_histogram, run_js_comp_test, &
                                                  run_js_comp_test_parameter_search
    use tox_data_integration_js_comp_test_impl, only: METHOD_JOIN_MIN, METHOD_JOIN_MAX, METHOD_JOIN_MEDIAN, &
                                                       MODE_PLATEAU_CI_OVERLAP, MODE_PLATEAU_EFFECT_SIZE, &
                                                       MODE_PLATEAU_BOTH, calc_js_comp_test_n_top_k_jsds, &
                                                       calc_js_comp_test_candidate_bounds, &
                                                       gather_pooled_neighborhood_residuals
    use tox_errors
    use test_suite, only: test_case

    implicit none

    real(real64), parameter :: TOL = 1d-12

contains

    !> Get array of all available tests.
    function get_all_tests_data_integration_js_comp_test() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(78))

        all_tests(1) = test_case("test_construct_neighborhoods_ranged_basic", test_construct_neighborhoods_ranged_basic)
        all_tests(2) = test_case("test_construct_neighborhoods_ranged_tie_extends_range", &
                                 test_construct_neighborhoods_ranged_tie_extends_range)
        all_tests(3) = test_case("test_construct_neighborhoods_ranged_all_nan_means", &
                                 test_construct_neighborhoods_ranged_all_nan_means)
        all_tests(4) = test_case("test_construct_neighborhoods_ranged_fewer_genes_than_neighbors", &
                                 test_construct_neighborhoods_ranged_fewer_genes_than_neighbors)
        all_tests(5) = test_case("test_construct_neighborhoods_ranged_plain_and_validation", &
                                 test_construct_neighborhoods_ranged_plain_and_validation)

        all_tests(6) = test_case("test_calc_pmf_matches_build_residual_histograms", &
                                 test_calc_pmf_matches_build_residual_histograms)
        all_tests(7) = test_case("test_calc_pmf_zero_included_reps", test_calc_pmf_zero_included_reps)
        all_tests(8) = test_case("test_calc_pmf_rejects_negative_counts", test_calc_pmf_rejects_negative_counts)

        all_tests(9) = test_case("test_estimate_bin_count_basic_hand_computed", test_estimate_bin_count_basic_hand_computed)
        all_tests(10) = test_case("test_estimate_bin_count_all_nan_gives_one_bin", &
                                  test_estimate_bin_count_all_nan_gives_one_bin)
        all_tests(11) = test_case("test_estimate_bin_count_clamped_to_max_n_bins", &
                                  test_estimate_bin_count_clamped_to_max_n_bins)

        all_tests(12) = test_case("test_generate_js_comp_test_candidates_collapses_at_8742", &
                                  test_generate_js_comp_test_candidates_collapses_at_8742)
        all_tests(13) = test_case("test_generate_js_comp_test_candidates_has_two_distinct_at_8743", &
                                  test_generate_js_comp_test_candidates_has_two_distinct_at_8743)
        all_tests(14) = test_case("test_generate_js_comp_test_candidates_validation", &
                                  test_generate_js_comp_test_candidates_validation)

        all_tests(15) = test_case("test_neighborhood_overlaps_all_pass", test_neighborhood_overlaps_all_pass)
        all_tests(16) = test_case("test_neighborhood_overlaps_exactly_at_threshold_passes", &
                                  test_neighborhood_overlaps_exactly_at_threshold_passes)
        all_tests(17) = test_case("test_neighborhood_overlaps_below_threshold_fails", &
                                  test_neighborhood_overlaps_below_threshold_fails)

        all_tests(18) = test_case("test_mean_pmf_min_counts_all_pass", test_mean_pmf_min_counts_all_pass)
        all_tests(19) = test_case("test_mean_pmf_min_counts_exactly_at_threshold_passes", &
                                  test_mean_pmf_min_counts_exactly_at_threshold_passes)
        all_tests(20) = test_case("test_mean_pmf_min_counts_below_threshold_fails", &
                                  test_mean_pmf_min_counts_below_threshold_fails)

        all_tests(21) = test_case("test_check_plateau_condition_join_min_requires_all", &
                                  test_check_plateau_condition_join_min_requires_all)
        all_tests(22) = test_case("test_check_plateau_condition_join_max_requires_any", &
                                  test_check_plateau_condition_join_max_requires_any)
        all_tests(23) = test_case("test_check_plateau_condition_join_median_requires_majority", &
                                  test_check_plateau_condition_join_median_requires_majority)
        all_tests(24) = test_case("test_check_plateau_condition_worse_than_previous_short_circuits", &
                                  test_check_plateau_condition_worse_than_previous_short_circuits)

        all_tests(25) = test_case("test_create_mean_pmf_two_studies_hand_computed", &
                                  test_create_mean_pmf_two_studies_hand_computed)
        all_tests(26) = test_case("test_create_mean_pmf_three_studies_includes_self", &
                                  test_create_mean_pmf_three_studies_includes_self)
        all_tests(27) = test_case("test_create_mean_pmf_only_matches_create_mean_pmf", &
                                  test_create_mean_pmf_only_matches_create_mean_pmf)
        all_tests(28) = test_case("test_create_mean_pmf_validation", test_create_mean_pmf_validation)

        all_tests(29) = test_case("test_calc_js_comp_test_n_top_k_jsds_default_hand_computed", &
                                  test_calc_js_comp_test_n_top_k_jsds_default_hand_computed)
        all_tests(30) = test_case("test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level", &
                                  test_calc_js_comp_test_n_top_k_jsds_explicit_significance_level)
        all_tests(31) = test_case("test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one", &
                                  test_calc_js_comp_test_n_top_k_jsds_clamped_to_at_least_one)

        all_tests(32) = test_case("test_bootstrap_histogram_seeded_reproducibility", &
                                  test_bootstrap_histogram_seeded_reproducibility)
        all_tests(33) = test_case("test_bootstrap_histogram_degenerate_single_bin_zero_ci", &
                                  test_bootstrap_histogram_degenerate_single_bin_zero_ci)

        all_tests(34) = test_case("test_gjct_permutation_test_conservation_of_counts", &
                                  test_gjct_permutation_test_conservation_of_counts)
        all_tests(35) = test_case("test_gjct_permutation_test_seeded_reproducibility", &
                                  test_gjct_permutation_test_seeded_reproducibility)
        all_tests(36) = test_case("test_permutation_pvalue_laplace_corrected_never_exactly_zero", &
                                  test_permutation_pvalue_laplace_corrected_never_exactly_zero)

        all_tests(37) = test_case("test_run_js_comp_test_two_studies_hand_traceable", &
                                  test_run_js_comp_test_two_studies_hand_traceable)
        all_tests(38) = test_case("test_run_js_comp_test_three_studies_outlier_has_small_p_value", &
                                  test_run_js_comp_test_three_studies_outlier_has_small_p_value)
        all_tests(39) = test_case("test_param_search_no_plateau_falls_back_to_finest", &
                                  test_param_search_no_plateau_falls_back_to_finest)
        all_tests(40) = test_case("test_param_search_single_candidate_bypasses_plateau", &
                                  test_param_search_single_candidate_bypasses_plateau)
        all_tests(41) = test_case("test_param_search_finds_plateau_mid_grid", &
                                  test_param_search_finds_plateau_mid_grid)

        all_tests(42) = test_case("test_effect_size_plateau_first_candidate_no_delta", &
                                  test_effect_size_plateau_first_candidate_no_delta)
        all_tests(43) = test_case("test_effect_size_plateau_single_transition_insufficient", &
                                  test_effect_size_plateau_single_transition_insufficient)
        all_tests(44) = test_case("test_effect_size_plateau_two_consecutive_transitions", &
                                  test_effect_size_plateau_two_consecutive_transitions)
        all_tests(45) = test_case("test_effect_size_plateau_resets_on_non_qualifying", &
                                  test_effect_size_plateau_resets_on_non_qualifying)
        all_tests(46) = test_case("test_effect_size_plateau_median_max_hand_computed", &
                                  test_effect_size_plateau_median_max_hand_computed)
        all_tests(47) = test_case("test_param_search_effect_size_mode_plateau", &
                                  test_param_search_effect_size_mode_plateau)
        all_tests(48) = test_case("test_param_search_both_mode_uses_earlier_trigger", &
                                  test_param_search_both_mode_uses_earlier_trigger)
        all_tests(49) = test_case("test_param_search_no_plateau_uses_smallest_uncertainty", &
                                  test_param_search_no_plateau_uses_smallest_uncertainty)
        all_tests(50) = test_case("test_param_search_no_plateau_effect_size_smallest_uncertainty", &
                                  test_param_search_no_plateau_effect_size_smallest_uncertainty)
        all_tests(51) = test_case("test_param_search_no_plateau_both_uses_smallest_uncertainty", &
                                  test_param_search_no_plateau_both_uses_smallest_uncertainty)

        all_tests(52) = test_case("test_estimate_bin_count_sturges_wins_when_greater_than_fd", &
                                  test_estimate_bin_count_sturges_wins_when_greater_than_fd)
        all_tests(53) = test_case("test_estimate_bin_count_near_zero_iqr_guard_falls_back_to_sturges", &
                                  test_estimate_bin_count_near_zero_iqr_falls_back_to_sturges)

        all_tests(54) = test_case("test_determine_bin_count_occupancy_finds_valid_below_m_max", &
                                  test_occupancy_finds_valid_below_m_max)
        all_tests(55) = test_case("test_determine_bin_count_occupancy_reaches_m_max_validly", &
                                  test_occupancy_reaches_m_max_validly)
        all_tests(56) = test_case("test_determine_bin_count_occupancy_m_min_itself_invalid_failure", &
                                  test_occupancy_m_min_itself_invalid_failure)
        all_tests(57) = test_case("test_determine_bin_count_occupancy_refinement_picks_above_m_valid", &
                                  test_occupancy_refinement_picks_above_m_valid)
        all_tests(58) = test_case("test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid", &
                                  test_occupancy_refinement_finds_nothing_above_m_valid)
        all_tests(59) = test_case("test_determine_bin_count_occupancy_all_residuals_nan", &
                                  test_occupancy_all_residuals_nan)
        all_tests(60) = test_case("test_determine_bin_count_occupancy_geometric_step_guarantees_progress", &
                                  test_occupancy_geometric_step_guarantees_progress)
        all_tests(61) = test_case("test_determine_bin_count_occupancy_diagnostics_hand_computed", &
                                  test_occupancy_diagnostics_hand_computed)
        all_tests(62) = test_case("test_determine_bin_count_occupancy_defaults_match_issue_suggestions", &
                                  test_occupancy_defaults_match_issue_suggestions)

        all_tests(63) = test_case("test_mean_pmf_min_counts_per_point_bins_differ_all_pass", &
                                  test_mean_pmf_min_counts_per_point_bins_differ_all_pass)
        all_tests(64) = test_case("test_mean_pmf_min_counts_per_point_bins_differ_one_point_fails", &
                                  test_mean_pmf_min_counts_per_point_bins_differ_one_point_fails)

        all_tests(65) = test_case("test_param_search_occupancy_failure_rejects_candidate", &
                                  test_param_search_occupancy_failure_rejects_candidate)
        all_tests(66) = test_case("test_param_search_different_neighborhoods_different_m_j", &
                                  test_param_search_different_neighborhoods_different_m_j)
        all_tests(67) = test_case("test_param_search_final_n_bins_matches_selected_trace_column", &
                                  test_param_search_final_n_bins_matches_selected_trace_column)
        all_tests(68) = test_case("test_run_js_comp_test_occupancy_failed_point_still_contributes", &
                                  test_run_js_comp_test_occupancy_failed_point_still_contributes)
        all_tests(69) = test_case("test_occupancy_min_residuals_per_bin_zero_accepted", &
                                  test_occupancy_min_residuals_per_bin_zero_accepted)
        all_tests(70) = test_case("test_run_js_comp_test_accepts_min_residuals_per_bin_zero", &
                                  test_run_js_comp_test_accepts_min_residuals_per_bin_zero)
        all_tests(71) = test_case("test_occupancy_range_asymmetric_skewed_residuals", &
                                  test_occupancy_range_asymmetric_skewed_residuals)
        all_tests(72) = test_case("test_occupancy_range_hand_computed_percentile", &
                                  test_occupancy_range_hand_computed_percentile)
        all_tests(73) = test_case("test_occupancy_exhaustive_matches_production_on_all_fixtures", &
                                  test_occupancy_exhaustive_matches_production_on_all_fixtures)
        all_tests(74) = test_case("test_occupancy_exhaustive_diverges_on_adversarial_fixture", &
                                  test_occupancy_exhaustive_diverges_on_adversarial_fixture)
        all_tests(75) = test_case("test_occupancy_exhaustive_hand_computed_own_correctness", &
                                  test_occupancy_exhaustive_hand_computed_own_correctness)
        all_tests(76) = test_case("test_gather_pooled_residuals_hand_computed", &
                                  test_gather_pooled_residuals_hand_computed)
        all_tests(77) = test_case("test_gather_pooled_residuals_rejects_empty_dimension", &
                                  test_gather_pooled_residuals_rejects_empty_dimension)
        all_tests(78) = test_case("test_gather_pooled_residuals_rejects_gene_index_oob", &
                                  test_gather_pooled_residuals_rejects_gene_index_oob)
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
    !| `1 + nint(log(1)/LOG_2) = 1`, clamped to `sturges_bins=1`. The 25th/75th percentiles (rank
    !| `0.25*9+1=3.25` and `0.75*9+1=7.75`) interpolate to `-1.75` and `2.75`, so the
    !| Freedman-Diaconis half-width is `(2.75 - (-1.75)) / 1^(1/3) = 4.5`; with
    !| `shared_residual_range=9.0`, `nint(9.0/4.5) = 2`, clamped to `fd_bins=2`. `fd_bins` beats
    !| `sturges_bins`, so `n_bins` should be `2`.
    subroutine test_estimate_bin_count_basic_hand_computed()
        integer(int32), parameter :: n_residuals = 10
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, sturges_bins, fd_bins, ierr

        residuals = [-4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
                    1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call estimate_bin_count(residuals, n_residuals, 1_int32, 1_int32, 9.0_real64, n_bins, sturges_bins, fd_bins, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_estimate_bin_count_basic_hand_computed: ierr should be OK")
        call assert_equal_int(sturges_bins, 1_int32, "test_estimate_bin_count_basic_hand_computed: sturges_bins should be 1")
        call assert_equal_int(fd_bins, 2_int32, "test_estimate_bin_count_basic_hand_computed: fd_bins should be 2")
        call assert_equal_int(n_bins, 2_int32, "test_estimate_bin_count_basic_hand_computed: n_bins should be 2")
    end subroutine test_estimate_bin_count_basic_hand_computed

    !> When every residual is NaN, the pool is empty and the routine must fall back to a single
    !| bin -- for all three outputs, since no Sturges/Freedman-Diaconis computation happens at
    !| all on this early-return path -- rather than computing percentiles of nothing.
    subroutine test_estimate_bin_count_all_nan_gives_one_bin()
        integer(int32), parameter :: n_residuals = 4
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, sturges_bins, fd_bins, ierr
        real(real64) :: nan_val

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)
        residuals = nan_val

        call estimate_bin_count(residuals, n_residuals, 1_int32, 1_int32, 9.0_real64, n_bins, sturges_bins, fd_bins, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_estimate_bin_count_all_nan_gives_one_bin: ierr should be OK")
        call assert_equal_int(sturges_bins, 1_int32, &
                              "test_estimate_bin_count_all_nan_gives_one_bin: sturges_bins should be 1")
        call assert_equal_int(fd_bins, 1_int32, "test_estimate_bin_count_all_nan_gives_one_bin: fd_bins should be 1")
        call assert_equal_int(n_bins, 1_int32, "test_estimate_bin_count_all_nan_gives_one_bin: n_bins should be 1")
    end subroutine test_estimate_bin_count_all_nan_gives_one_bin

    !> Reusing the basic fixture's residuals but with `shared_residual_range=5000.0`, `sturges_bins`
    !| is unchanged at `1`. The Freedman-Diaconis half-width is unchanged (`4.5`), so the raw
    !| estimate is `nint(5000.0/4.5) = 1111`, far above MAX_N_BINS; `fd_bins` (and therefore
    !| `n_bins`, since it beats `sturges_bins`) must clamp to 256.
    subroutine test_estimate_bin_count_clamped_to_max_n_bins()
        integer(int32), parameter :: n_residuals = 10
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, sturges_bins, fd_bins, ierr

        residuals = [-4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
                    1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call estimate_bin_count(residuals, n_residuals, 1_int32, 1_int32, 5000.0_real64, n_bins, sturges_bins, fd_bins, &
                                ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_estimate_bin_count_clamped_to_max_n_bins: ierr should be OK")
        call assert_equal_int(sturges_bins, 1_int32, &
                              "test_estimate_bin_count_clamped_to_max_n_bins: sturges_bins should be 1")
        call assert_equal_int(fd_bins, 256_int32, "test_estimate_bin_count_clamped_to_max_n_bins: fd_bins should clamp to 256")
        call assert_equal_int(n_bins, 256_int32, "test_estimate_bin_count_clamped_to_max_n_bins: n_bins should clamp to 256")
    end subroutine test_estimate_bin_count_clamped_to_max_n_bins

    !> Covers the branch where Sturges' estimate exceeds Freedman-Diaconis, so the combined
    !| `n_bins` equals `sturges_bins`, not `fd_bins` -- no existing test covered this before.
    !| Reuses the basic fixture's residuals (so Q25=-1.75, Q75=2.75, IQR=4.5 as in
    !| `test_estimate_bin_count_basic_hand_computed`), but with `max_n_reps_all_studies=8,
    !| n_neighbors=1` so `n_reps_neighborhood=8`. Sturges gives
    !| `1 + nint(log(8)/LOG_2) = 1 + nint(3.0) = 4`, clamped to `sturges_bins=4`. The
    !| Freedman-Diaconis half-width is `4.5 / 8^(1/3) = 4.5/2.0 = 2.25`; with
    !| `shared_residual_range=2.25`, `nint(2.25/2.25) = 1`, clamped to `fd_bins=1`. `sturges_bins`
    !| (4) beats `fd_bins` (1), so `n_bins` should be `4`.
    subroutine test_estimate_bin_count_sturges_wins_when_greater_than_fd()
        integer(int32), parameter :: n_residuals = 10
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, sturges_bins, fd_bins, ierr

        residuals = [-4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
                    1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64]

        call estimate_bin_count(residuals, n_residuals, 8_int32, 1_int32, 2.25_real64, n_bins, sturges_bins, fd_bins, &
                                ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_estimate_bin_count_sturges_wins_when_greater_than_fd: ierr should be OK")
        call assert_equal_int(sturges_bins, 4_int32, &
                              "test_estimate_bin_count_sturges_wins_when_greater_than_fd: sturges_bins should be 4")
        call assert_equal_int(fd_bins, 1_int32, &
                              "test_estimate_bin_count_sturges_wins_when_greater_than_fd: fd_bins should be 1")
        call assert_equal_int(n_bins, 4_int32, &
                              "test_estimate_bin_count_sturges_wins_when_greater_than_fd: n_bins should be sturges_bins=4")
    end subroutine test_estimate_bin_count_sturges_wins_when_greater_than_fd

    !> Covers the near-zero-IQR guard: four identical residuals give `quartile_75 - quartile_25
    !| = 0` exactly, so `half_bin_width = 0.0` and `is_close(half_bin_width, 0.0_real64)` fires,
    !| skipping the Freedman-Diaconis division entirely (avoiding a divide-by-near-zero blowup).
    !| `max_n_reps_all_studies=8, n_neighbors=1` gives `n_reps_neighborhood=8`, so Sturges gives
    !| `1 + nint(log(8)/LOG_2) = 4`, clamped to `sturges_bins=4`. Since the guard fires, `fd_bins`
    !| falls back to the (clamped) Sturges estimate instead of a division blow-up, so `fd_bins`
    !| should also be `4`, and `n_bins = max(4, 4) = 4`. `shared_residual_range` is irrelevant
    !| here since the Freedman-Diaconis term is never reached.
    subroutine test_estimate_bin_count_near_zero_iqr_falls_back_to_sturges()
        integer(int32), parameter :: n_residuals = 4
        real(real64) :: residuals(n_residuals)
        integer(int32) :: n_bins, sturges_bins, fd_bins, ierr

        residuals = [3.0_real64, 3.0_real64, 3.0_real64, 3.0_real64]

        call estimate_bin_count(residuals, n_residuals, 8_int32, 1_int32, 9.0_real64, n_bins, sturges_bins, fd_bins, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_estimate_bin_count_near_zero_iqr_guard_falls_back_to_sturges: ierr should be OK")
        call assert_equal_int(sturges_bins, 4_int32, &
                              "test_estimate_bin_count_near_zero_iqr_guard_falls_back_to_sturges: sturges_bins should be 4")
        call assert_equal_int(fd_bins, 4_int32, &
                              "test_estimate_bin_count_near_zero_iqr_guard_falls_back_to_sturges: "// &
                              "fd_bins should fall back to sturges_bins=4")
        call assert_equal_int(n_bins, 4_int32, &
                              "test_estimate_bin_count_near_zero_iqr_guard_falls_back_to_sturges: n_bins should be 4")
    end subroutine test_estimate_bin_count_near_zero_iqr_falls_back_to_sturges

    !> The documented small-N candidate-grid collapse, bracketed at its exact threshold: with
    !| `max_n_genes_all_studies=8742`, `n_points_high = clamp(ceil(4*sqrt(8742)), 300, 1500) = 374`
    !| and `n_points_low = max(300, ceil(0.2*374)) = 300`; after one grid iteration
    !| `n_points_high` becomes `374*0.8 = 299.2 < 300`, so the loop exits before a second distinct
    !| `n_points` value is ever produced -- every candidate pair the grid returns must share the
    !| same `n_points`. This is real, derived behavior the grid depends on, not a bug.
    subroutine test_generate_js_comp_test_candidates_collapses_at_8742()
        integer(int32) :: candidates(2, 16), n_candidates, ierr

        call generate_js_comp_test_candidates(8742_int32, candidates, n_candidates, ierr)

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
        integer(int32) :: candidates(2, 16), n_candidates, ierr
        integer(int32) :: n_distinct_n_points, i_candidate

        call generate_js_comp_test_candidates(8743_int32, candidates, n_candidates, ierr)

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

    !> `max_n_genes_all_studies` must be positive. Issue #187, Step 2.7: this routine's `_expert`
    !| tier no longer exists (no `tmp_`/work-array/permutation left in its signature once the
    !| bin-estimate side effect was removed), so this test no longer exercises it.
    subroutine test_generate_js_comp_test_candidates_validation()
        integer(int32) :: candidates(2, 16), n_candidates, ierr

        call generate_js_comp_test_candidates(0_int32, candidates, n_candidates, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_generate_js_comp_test_candidates_validation: max_n_genes_all_studies=0 rejected", &
                        arg_pos=1_int32)
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

    !> Every bin comfortably exceeds the minimum count. `n_bins_per_point` is uniformly `n_bins`
    !| here, so the per-point guard is a no-op and this exercises the same behavior the old
    !| uniform-`n_bins` signature did.
    subroutine test_mean_pmf_min_counts_all_pass()
        integer(int32), parameter :: n_bins = 2, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: n_bins_per_point(n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        mean_pmf_counts = reshape([10, 10, 10, 10], [n_bins, n_points])
        n_bins_per_point = n_bins

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_bins_per_point, n_points, 5_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_mean_pmf_min_counts_all_pass: ierr should be OK")
        call assert_true(all_pass, "test_mean_pmf_min_counts_all_pass: all bins above minimum should pass")
    end subroutine test_mean_pmf_min_counts_all_pass

    !> Every bin equals the minimum count exactly: the gate is `count >= min`, not a strict `>`,
    !| so this must still pass. `n_bins_per_point` is uniformly `n_bins`.
    subroutine test_mean_pmf_min_counts_exactly_at_threshold_passes()
        integer(int32), parameter :: n_bins = 2, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: n_bins_per_point(n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        mean_pmf_counts = reshape([3, 3, 3, 3], [n_bins, n_points])
        n_bins_per_point = n_bins

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_bins_per_point, n_points, 3_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_mean_pmf_min_counts_exactly_at_threshold_passes: ierr should be OK")
        call assert_true(all_pass, "test_mean_pmf_min_counts_exactly_at_threshold_passes: "// &
                         "count exactly at minimum must pass")
    end subroutine test_mean_pmf_min_counts_exactly_at_threshold_passes

    !> One bin (of four) is one below the minimum: the gate must fail. `n_bins_per_point` is
    !| uniformly `n_bins`.
    subroutine test_mean_pmf_min_counts_below_threshold_fails()
        integer(int32), parameter :: n_bins = 2, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: n_bins_per_point(n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        mean_pmf_counts = reshape([3, 3, 3, 2], [n_bins, n_points])
        n_bins_per_point = n_bins

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_bins_per_point, n_points, 3_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_mean_pmf_min_counts_below_threshold_fails: ierr should be OK")
        call assert_false(all_pass, "test_mean_pmf_min_counts_below_threshold_fails: "// &
                          "one bin below minimum must fail the whole gate")
    end subroutine test_mean_pmf_min_counts_below_threshold_fails

    !> Issue #187: two points with DIFFERENT `n_bins_per_point`, where each point's own
    !| valid-range columns satisfy `min_count`, but the padded columns beyond each point's own
    !| count are deliberately set BELOW `min_count`. If the `i_bin <= n_bins_per_point(i_point)`
    !| guard were removed or broken, this exact fixture would fail on the padded columns -- so a
    !| gate that still passes here proves the guard is load-bearing, not incidental.
    !| Point 1 has 3 valid bins (all count=10) plus 1 padded column (count=0, deliberately below
    !| min_count=5); point 2 has only 1 valid bin (count=10) plus 3 padded columns (count=0).
    subroutine test_mean_pmf_min_counts_per_point_bins_differ_all_pass()
        integer(int32), parameter :: n_bins = 4, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: n_bins_per_point(n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        n_bins_per_point = [3, 1]
        ! Point 1: valid bins 1:3 = 10 (>= min_count), padded bin 4 = 0 (< min_count, must be skipped)
        mean_pmf_counts(:, 1) = [10, 10, 10, 0]
        ! Point 2: valid bin 1 = 10 (>= min_count), padded bins 2:4 = 0 (< min_count, must be skipped)
        mean_pmf_counts(:, 2) = [10, 0, 0, 0]

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_bins_per_point, n_points, 5_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_mean_pmf_min_counts_per_point_bins_differ_all_pass: ierr should be OK")
        call assert_true(all_pass, "test_mean_pmf_min_counts_per_point_bins_differ_all_pass: "// &
                         "each point's own valid-range bins pass; padded columns below min_count "// &
                         "must be ignored by the i_bin <= n_bins_per_point(i_point) guard")
    end subroutine test_mean_pmf_min_counts_per_point_bins_differ_all_pass

    !> Issue #187: same per-point-bins-differ setup as above, but point 2's own valid-range bin
    !| genuinely fails to reach `min_count` this time. Confirms all-or-nothing is preserved across
    !| points, just correctly scoped per point now -- a genuine failure inside a point's own valid
    !| range still fails the whole gate, even though point 1 is fine.
    subroutine test_mean_pmf_min_counts_per_point_bins_differ_one_point_fails()
        integer(int32), parameter :: n_bins = 4, n_points = 2
        integer(int32) :: mean_pmf_counts(n_bins, n_points)
        integer(int32) :: n_bins_per_point(n_points)
        logical(c_bool) :: all_pass
        integer(int32) :: ierr

        n_bins_per_point = [3, 1]
        ! Point 1: valid bins 1:3 = 10 (>= min_count), padded bin 4 = 0 (must be skipped)
        mean_pmf_counts(:, 1) = [10, 10, 10, 0]
        ! Point 2: valid bin 1 = 2 (< min_count = 5, a genuine failure), padded bins 2:4 = 0
        mean_pmf_counts(:, 2) = [2, 0, 0, 0]

        call check_mean_pmf_min_counts(mean_pmf_counts, n_bins, n_bins_per_point, n_points, 5_int32, all_pass, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_mean_pmf_min_counts_per_point_bins_differ_one_point_fails: ierr should be OK")
        call assert_false(all_pass, "test_mean_pmf_min_counts_per_point_bins_differ_one_point_fails: "// &
                          "point 2's own valid-range bin genuinely below min_count must fail the "// &
                          "whole gate, even though point 1's own valid range is fine")
    end subroutine test_mean_pmf_min_counts_per_point_bins_differ_one_point_fails

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

    !> Verifies the `(1+count)/(n+1)` Laplace add-one correction in `gjct_permutation_test_impl`:
    !| even in the worst-case, deterministic scenario where `count` is guaranteed to be `0` across
    !| every permutation, `p` must come out as `1/(n_permutations+1)`, never exactly `0.0`.
    !| Contrived, deterministically (not merely probabilistically) so: a single-bin histogram
    !| (`n_bins=1`) means every resample's pmf is `[1.0]`, for both the study and the (unperturbed)
    !| consensus mean alike, so the Jensen-Shannon divergence between them is exactly `0.0` on
    !| *every* permutation, regardless of which counts `random_multiv_hypergeom` happens to draw.
    !| Observing a strictly positive `global_jsd_observed` then guarantees no permutation's JSD can
    !| reach it, so `count=0` for all `n_permutations=10` permutations, and the corrected formula
    !| gives `p_values = anint(0+1)/(10+1) = 1/11`.
    subroutine test_permutation_pvalue_laplace_corrected_never_exactly_zero()
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
                              "test_permutation_pvalue_laplace_corrected_never_exactly_zero: ierr should be OK")
        call assert_equal_real(p_values(1), 1.0_real64/11.0_real64, TOL, &
                               "test_permutation_pvalue_laplace_corrected_never_exactly_zero: "// &
                               "Laplace correction: p must be (0+1)/(n_permutations+1) = 1/11, not 0.0")
    end subroutine test_permutation_pvalue_laplace_corrected_never_exactly_zero

    !> End-to-end 2-study case run through `run_js_comp_test` with `n_permutations=0`: GSL's
    !| `create_rng` still runs (so `ierr` is genuinely exercised), but the permutation loop itself
    !| never executes, so `gjct_permutation_test_impl` leaves everything untouched and
    !| `run_js_comp_test_impl`'s final re-derivation step reproduces exactly the same values it
    !| already had -- making the whole pipeline deterministic and hand-traceable.
    !|
    !| One reference point, one neighbor per study (`x_star=1.0` sits exactly on gene 1's mean, so
    !| both studies trivially pick gene 1 -- see test_construct_neighborhoods_ranged_basic's own
    !| binary-search trace for why an exact match gives insertion index 1). Study 1's two replicate
    !| residuals are `[-1, 1]`, study 2's are `[-3, 3]`, pooled into `[-1, 1, -3, 3]`, sorted
    !| `[-3, -1, 1, 3]` (`N_j=4`).
    !|
    !| Step 3 (per-neighborhood residual range): `shared_residual_range_low`/`_high` are now this
    !| point's own 5th/95th percentile of its own pooled residuals, not a caller-supplied scalar.
    !| `rank(0.05,4)=0.05*3+1=1.15` -> interpolate value(1)=-3, value(2)=-1 at fraction 0.15 ->
    !| R_low = -3 + 0.15*2 = -2.7; `rank(0.95,4)=0.95*3+1=3.85` -> interpolate value(3)=1,
    !| value(4)=3 at fraction 0.85 -> R_high = 1 + 0.85*2 = 2.7 (span 5.4, narrower than the old
    !| hand-picked `R=4.0`).
    !|
    !| Issue #187's occupancy search picks the bin count: `min_residuals_per_bin=1` is passed
    !| explicitly (its default, 10, could never be satisfied by only 4 pooled residuals). Starting
    !| from the default `m_min=3`, bin width is `5.4/3=1.8`, boundaries `[-2.7,-0.9),[-0.9,0.9),
    !| [0.9,2.7]`: `-3`(clamped to `-2.7`) and `-1` both land in bin 1, `1` and `3`(clamped to
    !| `2.7`) both land in bin 3, leaving bin 2 completely EMPTY -- `min_occ=0 < 1`, inadmissible
    !| already at `m_min`. This is the FAILURE case from Issue #187's own pseudocode: even the
    !| minimum resolution is unsupported, so `occupancy_failed(1) = .true.` and
    !| `n_bins_per_point(1)` stays at `m_min = 3` (the routine's own documented behavior: a point
    !| that fails occupancy still gets a real histogram built at `m_min`, and still contributes to
    !| `global_js_divergence` -- see this routine's own doc block).
    !|
    !| At `M=3`, `R_low=-2.7`, `R_high=2.7` (bin width 1.8): study 1's `[-1, 1]` -> bin 1 (`-1`),
    !| bin 3 (`1`) -> counts `[1, 0, 1]`. Study 2's `[-3, 3]`, both clamped to `[-2.7, 2.7]` ->
    !| bin 1 (`-2.7`), bin 3 (`2.7`) -> counts `[1, 0, 1]` -- IDENTICAL to study 1's, an accidental
    !| consequence of this fixture's own symmetric residuals under the new narrower, still-symmetric
    !| derived range (a coincidence of this particular data, not a general property of the
    !| asymmetric-range design). Both studies' pmf is therefore `[0.5, 0, 0.5]`, the mean pmf is the
    !| same `[0.5, 0, 0.5]`, and since each study's own pmf exactly equals the consensus,
    !| `global_js_divergence` is exactly `0.0` for both -- no closed-form derivation needed, unlike
    !| the pre-Step-3 antisymmetric-bin-permutation case this fixture used to exercise.
    subroutine test_run_js_comp_test_two_studies_hand_traceable()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 2, max_n_reps_all_studies = 2
        integer(int32), parameter :: n_points = 1, n_neighbors = 1
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        integer(int32) :: gene_means_perms(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        real(real64) :: x_star(n_points)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points, n_studies)
        integer(int32) :: neighborhood_range(2, n_points, n_studies)
        integer(int32) :: n_bins_per_point(n_points), max_n_bins_per_point
        real(real64) :: shared_residual_range_low(n_points), shared_residual_range_high(n_points)
        logical(c_bool) :: occupancy_failed(n_points)
        integer(int32) :: n_pooled_residuals(n_points), min_bin_occupancy(n_points), max_bin_occupancy(n_points)
        real(real64) :: mean_bin_occupancy(n_points)
        integer(int32) :: sturges_bins(n_points), fd_bins(n_points)
        real(real64) :: pmfs(256, n_points, n_studies)
        integer(int32) :: counts(256, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(256, n_points)
        integer(int32) :: mean_pmf_counts(256, n_points)
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

        call run_js_comp_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                              gene_means, gene_means_perms, residuals, x_star, neighborhood_indices, &
                              neighborhood_range, n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                              max_n_bins_per_point, occupancy_failed, &
                              n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, &
                              sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, p_values, &
                              ierr=ierr, n_permutations=0_int32, random_seed=1_int32, min_residuals_per_bin=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, "test_run_js_comp_test_two_studies_hand_traceable: ierr should be OK")

        call assert_equal_array_int(neighborhood_indices(:, 1, 1), [1], n_neighbors, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 1 neighbor is gene 1")
        call assert_equal_array_int(neighborhood_indices(:, 1, 2), [1], n_neighbors, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 2 neighbor is gene 1")

        call assert_equal_real(shared_residual_range_low(1), -2.7_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: "// &
                               "shared_residual_range_low == 5th percentile of pooled [-3,-1,1,3]")
        call assert_equal_real(shared_residual_range_high(1), 2.7_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: "// &
                               "shared_residual_range_high == 95th percentile of pooled [-3,-1,1,3]")
        call assert_equal_int(n_bins_per_point(1), 3_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: "// &
                              "occupancy search FAILS at m_min=3 (bin 2 empty), so n_bins_per_point stays m_min")
        call assert_equal_int(max_n_bins_per_point, 3_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: max_n_bins_per_point is 3")
        call assert_true(occupancy_failed(1), &
                         "test_run_js_comp_test_two_studies_hand_traceable: "// &
                         "occupancy genuinely fails at the default m_min under the new narrower derived range")

        call assert_equal_array_int(counts(1:3, 1, 1), [1, 0, 1], 3_int32, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 1 counts")
        call assert_equal_array_int(counts(1:3, 1, 2), [1, 0, 1], 3_int32, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: study 2 counts")
        call assert_equal_int(included_n_reps(1, 1), 2_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: study 1 included_n_reps")
        call assert_equal_int(included_n_reps(1, 2), 2_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: study 2 included_n_reps")

        call assert_equal_array_real(mean_pmf(1:3, 1), [0.5_real64, 0.0_real64, 0.5_real64], 3_int32, TOL, &
                                     "test_run_js_comp_test_two_studies_hand_traceable: "// &
                                     "mean_pmf equals both studies' own identical pmf")
        call assert_equal_array_int(mean_pmf_counts(1:3, 1), [2, 0, 2], 3_int32, &
                                    "test_run_js_comp_test_two_studies_hand_traceable: mean_pmf_counts")
        call assert_equal_int(mean_pmf_included_n_reps(1), 4_int32, &
                              "test_run_js_comp_test_two_studies_hand_traceable: mean_pmf_included_n_reps")

        call assert_equal_real(global_js_divergence(1), 0.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: "// &
                               "study 1's pmf equals the consensus exactly -> JSD == 0")
        call assert_equal_real(global_js_divergence(2), 0.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: "// &
                               "study 2's pmf equals the consensus exactly -> JSD == 0")
        call assert_equal_real(weights(1, 1), 1.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: single reference point weight is 1.0")
        call assert_equal_real(weights(1, 2), 1.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: single reference point weight is 1.0")

        call assert_equal_real(p_values(1), 0.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: n_permutations=0 -> p_values stay 0.0")
        call assert_equal_real(p_values(2), 0.0_real64, TOL, &
                               "test_run_js_comp_test_two_studies_hand_traceable: n_permutations=0 -> p_values stay 0.0")
    end subroutine test_run_js_comp_test_two_studies_hand_traceable

    !> Compliance-review fix (Issue #187 cleanup): `min_residuals_per_bin`'s `DM_MIN` used to
    !| disagree between this routine (`1`) and `run_js_comp_test_parameter_search` (`0`) even
    !| though both forward the same argument into the same
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]
    !| occupancy search -- a caller passing `0` was rejected here while accepted (and relied upon)
    !| by the sibling entry point's own tests. Fixed by loosening this routine's bound to `0` too.
    !|
    !| Identical fixture to `test_run_js_comp_test_two_studies_hand_traceable` above (same
    !| `gene_means`/`residuals`/`x_star`, pooled residuals `[-1, 1, -3, 3]`, own
    !| `shared_residual_range_low`/`_high` = -2.7/2.7 -- see that test's docstring for the
    !| derivation), except `min_residuals_per_bin=0_int32` (proving the loosened bound is genuinely
    !| accepted, not just documentation) and an explicit `m_max=4_int32` to bound the occupancy
    !| search's early-return branch deterministically -- with `min_residuals_per_bin=0`, occupancy
    !| is trivially satisfied at every candidate `M` (`min_occ >= 0` always holds), so without a
    !| bounded `m_max` the search would run all the way to the default `m_max=120` instead of
    !| stopping at `M=4`. Unlike the hand-traceable test above (which FAILS at the default `m_min=3`
    !| and never reaches `M=4` at all), forcing `M=4` here reproduces the SAME per-bin split that
    !| test's own pre-Step-3 fixture had (`[0,1,1,0]`/`[1,0,0,1]`): at `M=4` the still-symmetric
    !| `[-2.7, 2.7]` range splits this data into the same relative quarters the old `[-4, 4]` range
    !| did, so every downstream numeric result (`counts`, `mean_pmf`, `global_js_divergence`) is
    !| unchanged from before Step 3 -- only `occupancy_failed`'s reason for being `.false.` differs
    !| (trivially satisfied here vs. genuinely satisfied in a from-scratch search).
    subroutine test_run_js_comp_test_accepts_min_residuals_per_bin_zero()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 2, max_n_reps_all_studies = 2
        integer(int32), parameter :: n_points = 1, n_neighbors = 1
        real(real64), parameter :: LOG2_3 = 1.5849625007211562_real64 ! log2(3) = ln(3)/ln(2)
        real(real64), parameter :: EXPECTED_JSD = 1.5_real64 - 0.75_real64*LOG2_3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        integer(int32) :: gene_means_perms(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        real(real64) :: x_star(n_points)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points, n_studies)
        integer(int32) :: neighborhood_range(2, n_points, n_studies)
        integer(int32) :: n_bins_per_point(n_points), max_n_bins_per_point
        real(real64) :: shared_residual_range_low(n_points), shared_residual_range_high(n_points)
        logical(c_bool) :: occupancy_failed(n_points)
        integer(int32) :: n_pooled_residuals(n_points), min_bin_occupancy(n_points), max_bin_occupancy(n_points)
        real(real64) :: mean_bin_occupancy(n_points)
        integer(int32) :: sturges_bins(n_points), fd_bins(n_points)
        real(real64) :: pmfs(256, n_points, n_studies)
        integer(int32) :: counts(256, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(256, n_points)
        integer(int32) :: mean_pmf_counts(256, n_points)
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

        call run_js_comp_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                              gene_means, gene_means_perms, residuals, x_star, neighborhood_indices, &
                              neighborhood_range, n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                              max_n_bins_per_point, occupancy_failed, &
                              n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, &
                              sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, p_values, &
                              ierr=ierr, n_permutations=0_int32, random_seed=1_int32, min_residuals_per_bin=0_int32, &
                              m_max=4_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: "// &
                              "min_residuals_per_bin=0 must be accepted, not rejected by DM_MIN")
        call assert_equal_real(shared_residual_range_low(1), -2.7_real64, TOL, &
                               "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: shared_residual_range_low")
        call assert_equal_real(shared_residual_range_high(1), 2.7_real64, TOL, &
                               "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: shared_residual_range_high")

        call assert_equal_int(n_bins_per_point(1), 4_int32, &
                              "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: "// &
                              "occupancy search reaches m_max=4 (trivially satisfied throughout)")
        call assert_equal_int(max_n_bins_per_point, 4_int32, &
                              "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: max_n_bins_per_point is 4")
        call assert_false(occupancy_failed(1), &
                          "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: "// &
                          "min_residuals_per_bin=0 trivially satisfies the occupancy criterion")

        call assert_equal_array_int(counts(1:4, 1, 1), [0, 1, 1, 0], 4_int32, &
                                    "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: study 1 counts")
        call assert_equal_array_int(counts(1:4, 1, 2), [1, 0, 0, 1], 4_int32, &
                                    "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: study 2 counts")

        call assert_equal_array_real(mean_pmf(1:4, 1), [0.25_real64, 0.25_real64, 0.25_real64, 0.25_real64], 4_int32, TOL, &
                                     "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: mean_pmf should be uniform")

        call assert_equal_real(global_js_divergence(1), EXPECTED_JSD, 1d-9, &
                               "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: study 1 global JSD, closed form")
        call assert_equal_real(global_js_divergence(2), EXPECTED_JSD, 1d-9, &
                               "test_run_js_comp_test_accepts_min_residuals_per_bin_zero: study 2 global JSD, closed form")
    end subroutine test_run_js_comp_test_accepts_min_residuals_per_bin_zero

    !> Three studies, same single-point/single-neighbor topology as the hand-traceable case above.
    !| Studies 1 and 2 are identical to each other (residuals `[-1,-1,1,1]`, pmf `[0, 0.5, 0.5, 0]`);
    !| study 3 is a deliberately constructed outlier whose every replicate lands in the SAME bin
    !| (`[-3.9,-3.9,-3.9,-3.9]`, pmf `[1, 0, 0, 0]`). `m_min=4` is passed explicitly so Issue #187's
    !| occupancy search reproduces the same forced `n_bins=4` this test hand-traced before Issue
    !| #187, but Step 3's own per-neighborhood range changes WHICH bin ends up empty: pooled
    !| residuals sorted are `[-3.9,-3.9,-3.9,-3.9,-1,-1,-1,-1,1,1,1,1]` (`N=12`), giving
    !| `rank(0.05,12)=1.55` -> `R_low=-3.9` (interpolating between two `-3.9` values) and
    !| `rank(0.95,12)=11.45` -> `R_high=1.0` (interpolating between two `1` values) -- an
    !| asymmetric range, unlike the old caller-supplied symmetric `[-4,4]`. At `M=4` over
    !| `[-3.9, 1.0]` (bin width `1.225`), the SECOND bin (`[-2.675,-1.45)`) is the one that ends up
    !| empty, not the fourth: the `-3.9`s land in bin 1, both studies' `-1`s land in bin 3, and both
    !| studies' `1`s land in bin 4 (clamped/capped to the last bin). Occupancy still fails even at
    !| `m=4` (an empty bin can never reach any positive `min_residuals_per_bin`) regardless of that
    !| threshold's value: `occupancy_failed(1) = .true.` and `n_bins_per_point(1) = m_min = 4` by
    !| construction, exactly the asymmetry `run_js_comp_test_impl`'s own doc block describes -- this
    !| occupancy-failed point still gets a real histogram and still contributes to
    !| `global_js_divergence`/the permutation test below, entirely unaffected by the occupancy
    !| search's own verdict. The consensus mean pmf, averaging all three, is `[1/3, 0, 1/3, 1/3]`
    !| with `mean_pmf_counts=[4, 0, 4, 4]`. Study 3's own pmf puts everything in a bin that holds
    !| only 4 of the pooled pool's 12 replicates, so drawing (without
    !| replacement, `random_multiv_hypergeom`) another 4-for-4 landing entirely in that one bin
    !| purely by chance is exceedingly rare -- its empirical p-value must come out small, while
    !| studies 1/2's own draws, being close to what the consensus itself is built from, should not
    !| be nearly as extreme. A fixed random_seed makes the outcome fully deterministic, so this is
    !| not a flaky probabilistic assertion.
    subroutine test_run_js_comp_test_three_studies_outlier_has_small_p_value()
        integer(int32), parameter :: n_studies = 3, max_n_genes_all_studies = 2, max_n_reps_all_studies = 4
        integer(int32), parameter :: n_points = 1, n_neighbors = 1, n_permutations = 500
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        integer(int32) :: gene_means_perms(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        real(real64) :: x_star(n_points)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points, n_studies)
        integer(int32) :: neighborhood_range(2, n_points, n_studies)
        integer(int32) :: n_bins_per_point(n_points), max_n_bins_per_point
        real(real64) :: shared_residual_range_low(n_points), shared_residual_range_high(n_points)
        logical(c_bool) :: occupancy_failed(n_points)
        integer(int32) :: n_pooled_residuals(n_points), min_bin_occupancy(n_points), max_bin_occupancy(n_points)
        real(real64) :: mean_bin_occupancy(n_points)
        integer(int32) :: sturges_bins(n_points), fd_bins(n_points)
        real(real64) :: pmfs(256, n_points, n_studies)
        integer(int32) :: counts(256, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(256, n_points)
        integer(int32) :: mean_pmf_counts(256, n_points)
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

        call run_js_comp_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                              gene_means, gene_means_perms, residuals, x_star, neighborhood_indices, &
                              neighborhood_range, n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                              max_n_bins_per_point, occupancy_failed, &
                              n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, &
                              sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, p_values, &
                              ierr=ierr, n_permutations=n_permutations, random_seed=42_int32, m_min=4_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_run_js_comp_test_three_studies_outlier_has_small_p_value: ierr should be OK")

        call assert_equal_real(shared_residual_range_low(1), -3.9_real64, TOL, &
                               "test_run_js_comp_test_three_studies_outlier_has_small_p_value: shared_residual_range_low")
        call assert_equal_real(shared_residual_range_high(1), 1.0_real64, TOL, &
                               "test_run_js_comp_test_three_studies_outlier_has_small_p_value: shared_residual_range_high")
        call assert_equal_int(n_bins_per_point(1), 4_int32, &
                              "test_run_js_comp_test_three_studies_outlier_has_small_p_value: n_bins_per_point == m_min")
        call assert_true(occupancy_failed(1), &
                         "test_run_js_comp_test_three_studies_outlier_has_small_p_value: "// &
                         "the pooled pool's 2nd bin is always empty, so occupancy fails even at m_min")

        call assert_equal_array_int(mean_pmf_counts(1:4, 1), [4, 0, 4, 4], 4_int32, &
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

    !> Issue #187's own documented asymmetry between `run_js_comp_test` and
    !| `run_js_comp_test_parameter_search`: THIS routine has no admissibility gate, so a reference
    !| point whose occupancy search fails even at `m_min` still gets a real histogram built and
    !| still contributes to `global_js_divergence` -- it is never rejected or skipped, unlike what
    !| `run_js_comp_test_parameter_search`'s own `check_mean_pmf_min_counts` gate would do to an
    !| entire candidate containing such a point.
    !|
    !| Two reference points, one neighbor per study, 2 studies: point 1 (`x_star=1.0`, gene 1) is
    !| the SAME fixture as the hand-traceable test above (residuals `[-1,1]`/`[-3,3]`,
    !| `min_residuals_per_bin=1` override) -- own `shared_residual_range_low`/`_high` = -2.7/2.7,
    !| own occupancy search FAILS at the default `m_min=3` (bin 2 empty), giving
    !| `n_bins_per_point(1)=3`, `occupancy_failed(1)=.true.`, and (as that test's own docstring
    !| derives) both studies land on the IDENTICAL pmf `[0.5, 0, 0.5]`, so point 1's own JSD is
    !| exactly `0.0` for both studies -- see that test's docstring for the full derivation, not
    !| repeated here. Point 2 (`x_star=5.0`, gene 2) is deliberately data-starved: study 1's two
    !| replicates are both `-3.9`; study 2's are `-3.9`/`3.9` (one of each). Pooled and sorted:
    !| `[-3.9,-3.9,-3.9,3.9]` (`N=4`). `rank(0.05,4)=1.15` -> `R_low=-3.9` (interpolating between
    !| two `-3.9` values); `rank(0.95,4)=3.85` -> `R_high=-3.9+0.85*(3.9-(-3.9))=2.73`. At the
    !| default `m_min=3` over `[-3.9, 2.73]` (bin width `2.21`), the middle bin
    !| (`[-1.69, 0.52)`) is never populated by any of these three-out-of-four-residuals-at--3.9
    !| values, so occupancy fails already at `m_min=3` and `n_bins_per_point(2) = m_min = 3` by
    !| construction (`FAILURE` per Issue #187's own policy) -- coincidentally the same `m=3` this
    !| point's occupancy search would have failed at under the old caller-supplied symmetric range
    !| too, though the range bounds themselves differ.
    !|
    !| At `n_bins_per_point(2)=3`, study 1's own pmf at point 2 is `[1, 0, 0]` (both replicates in
    !| bin 1, `-3.9` sitting exactly at `R_low`), study 2's is `[0.5, 0, 0.5]` (one replicate each
    !| in bins 1 and 3, `3.9` clamped to `R_high=2.73`); their consensus mean pmf is
    !| `[0.75, 0, 0.25]` -- numerically IDENTICAL to what the pre-Step-3 fixture had, since both the
    !| old symmetric and new asymmetric ranges happen to place these particular residuals in the
    !| same relative bins. Working `compute_divergence_per_reference_point_impl`'s own formula out
    !| by hand for each study against that consensus (bin 2 skipped throughout, since `S_mean=0`
    !| there for both) therefore reproduces the exact same pre-Step-3 constants:
    !|
    !| - Study 1: bin 1 contributes `1*ln(1/0.875) + 0.75*ln(0.75/0.875)` (`S_mean=0.875`), bin 3
    !|   contributes `0.25*ln(0.25/0.125)` (`S_mean=0.125`); halved and rescaled by `/LOG_2` gives
    !|   `POINT2_JSD_STUDY1 ~= 0.13792538097002990`.
    !| - Study 2: bin 1 contributes `0.5*ln(0.5/0.625) + 0.75*ln(0.75/0.625)` (`S_mean=0.625`), bin 3
    !|   contributes `0.5*ln(0.5/0.375) + 0.25*ln(0.25/0.375)` (`S_mean=0.375`); halved and rescaled
    !|   gives `POINT2_JSD_STUDY2 ~= 0.048794940695398498`.
    !|
    !| Both values are real, positive, and clearly distinct from each other, from zero, and from
    !| point 1's own (now exactly `0.0`) JSD -- proof this occupancy-failed point is not silently
    !| zeroed, skipped, or coincidentally aliased onto point 1's own value, but genuinely computed
    !| from its own (degenerate) data.
    !|
    !| Both points have `included_n_reps=2` for every study and no NaN anywhere, so
    !| `compute_weighted_global_divergence_impl`'s weights come out equal (`0.5` each: `(2+4)/12`
    !| for either point, `total_sample_count=12`). `global_js_divergence` is therefore exactly
    !| `0.5*0.0 + 0.5*point_2_jsd_study_i` per study -- clearly different from `0.0` alone, which is
    !| what would happen if point 2's occupancy-failed contribution were instead silently rejected
    !| or zero-weighted.
    subroutine test_run_js_comp_test_occupancy_failed_point_still_contributes()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 2, max_n_reps_all_studies = 2
        integer(int32), parameter :: n_points = 2, n_neighbors = 1
        real(real64), parameter :: POINT1_JSD = 0.0_real64 ! see hand-traceable test's own docstring
        real(real64), parameter :: POINT2_JSD_STUDY1 = 0.13792538097002990_real64 ! hand-derived, see doc block above
        real(real64), parameter :: POINT2_JSD_STUDY2 = 0.048794940695398498_real64 ! hand-derived, see doc block above
        real(real64), parameter :: EXPECTED_GLOBAL_JSD_STUDY1 = 0.5_real64*POINT1_JSD + 0.5_real64*POINT2_JSD_STUDY1
        real(real64), parameter :: EXPECTED_GLOBAL_JSD_STUDY2 = 0.5_real64*POINT1_JSD + 0.5_real64*POINT2_JSD_STUDY2
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        integer(int32) :: gene_means_perms(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        real(real64) :: x_star(n_points)
        integer(int32) :: neighborhood_indices(n_neighbors, n_points, n_studies)
        integer(int32) :: neighborhood_range(2, n_points, n_studies)
        integer(int32) :: n_bins_per_point(n_points), max_n_bins_per_point
        real(real64) :: shared_residual_range_low(n_points), shared_residual_range_high(n_points)
        logical(c_bool) :: occupancy_failed(n_points)
        integer(int32) :: n_pooled_residuals(n_points), min_bin_occupancy(n_points), max_bin_occupancy(n_points)
        real(real64) :: mean_bin_occupancy(n_points)
        integer(int32) :: sturges_bins(n_points), fd_bins(n_points)
        real(real64) :: pmfs(256, n_points, n_studies)
        integer(int32) :: counts(256, n_points, n_studies)
        integer(int32) :: included_n_reps(n_points, n_studies)
        real(real64) :: mean_pmf(256, n_points)
        integer(int32) :: mean_pmf_counts(256, n_points)
        integer(int32) :: mean_pmf_included_n_reps(n_points)
        real(real64) :: js_divergences(n_points, n_studies), weights(n_points, n_studies)
        real(real64) :: global_js_divergence(n_studies), p_values(n_studies)
        integer(int32) :: ierr

        gene_means(:, 1) = [1.0_real64, 5.0_real64]
        gene_means(:, 2) = [1.0_real64, 5.0_real64]
        gene_means_perms(:, 1) = [1, 2]
        gene_means_perms(:, 2) = [1, 2]

        ! Gene 1 (point 1's neighbor): well-supported, same fixture as the hand-traceable test.
        residuals(:, 1, 1) = [-1.0_real64, 1.0_real64]
        residuals(:, 1, 2) = [-3.0_real64, 3.0_real64]
        ! Gene 2 (point 2's neighbor): deliberately data-starved -- 3 of the 4 pooled residuals sit
        ! at -3.9, so the pooled pool's middle bin is permanently empty regardless of bin count.
        residuals(:, 2, 1) = [-3.9_real64, -3.9_real64]
        residuals(:, 2, 2) = [-3.9_real64, 3.9_real64]

        x_star = [1.0_real64, 5.0_real64]

        call run_js_comp_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, n_points, n_neighbors, &
                              gene_means, gene_means_perms, residuals, x_star, neighborhood_indices, &
                              neighborhood_range, n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                              max_n_bins_per_point, occupancy_failed, &
                              n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, max_bin_occupancy, &
                              sturges_bins, fd_bins, pmfs, counts, included_n_reps, mean_pmf, mean_pmf_counts, &
                              mean_pmf_included_n_reps, js_divergences, weights, global_js_divergence, p_values, &
                              ierr=ierr, n_permutations=0_int32, random_seed=1_int32, min_residuals_per_bin=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_run_js_comp_test_occupancy_failed_point_still_contributes: ierr should be OK")

        call assert_true(occupancy_failed(1), &
                         "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                         "point 1 now fails occupancy too under the new narrower derived range "// &
                         "(see the hand-traceable test's own docstring)")
        call assert_equal_int(n_bins_per_point(1), 3_int32, &
                              "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                              "point 1's own bin count")

        call assert_true(occupancy_failed(2), &
                         "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                         "point 2's occupancy search fails even at m_min")
        call assert_equal_int(n_bins_per_point(2), 3_int32, &
                              "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                              "point 2's n_bins_per_point == m_min (default 3), per Issue #187's own FAILURE policy")

        call assert_equal_int(max_n_bins_per_point, 3_int32, &
                              "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                              "max_n_bins_per_point: both points now select 3 bins")

        ! The occupancy-failed point still gets a REAL, non-skipped, non-zeroed contribution: its
        ! own weight is positive (both replicates counted, nothing excluded) and its JSD is a real,
        ! hand-derived, distinct-per-study value, not some sentinel or a value coincidentally
        ! aliased onto point 1's own JSD.
        call assert_true(weights(2, 1) > 0.0_real64, &
                         "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                         "point 2 has a positive, non-skipped weight in study 1")
        call assert_true(weights(2, 2) > 0.0_real64, &
                         "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                         "point 2 has a positive, non-skipped weight in study 2")
        call assert_equal_real(js_divergences(2, 1), POINT2_JSD_STUDY1, 1d-9, &
                               "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                               "point 2's own JSD, study 1, is the real, hand-computed value, not a sentinel")
        call assert_equal_real(js_divergences(2, 2), POINT2_JSD_STUDY2, 1d-9, &
                               "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                               "point 2's own JSD, study 2, is the real, hand-computed value, not a sentinel")

        ! The global JSD genuinely reflects point 2's contribution: it clearly differs from point
        ! 1's own JSD alone, which is what would happen if the occupancy-failed point had instead
        ! been silently rejected or zero-weighted.
        call assert_true(abs(global_js_divergence(1) - POINT1_JSD) > 0.01_real64, &
                         "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                         "global JSD, study 1, is measurably changed by point 2's real contribution")
        call assert_true(abs(global_js_divergence(2) - POINT1_JSD) > 0.01_real64, &
                         "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                         "global JSD, study 2, is measurably changed by point 2's real contribution")
        call assert_equal_real(global_js_divergence(1), EXPECTED_GLOBAL_JSD_STUDY1, 1d-9, &
                               "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                               "global JSD, study 1, closed form")
        call assert_equal_real(global_js_divergence(2), EXPECTED_GLOBAL_JSD_STUDY2, 1d-9, &
                               "test_run_js_comp_test_occupancy_failed_point_still_contributes: "// &
                               "global JSD, study 2, closed form")
    end subroutine test_run_js_comp_test_occupancy_failed_point_still_contributes

    !> Forces `run_js_comp_test_parameter_search`'s second admissibility gate
    !| (`min_residuals_per_bin`) impossibly high, so no candidate in the grid can ever pass it and
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
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=1000000_int32, random_seed=1_int32)

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
        ! Issue #187, Step 2.7: generate_js_comp_test_candidates_impl no longer computes any bin
        ! estimate, so this fallback (candidate 1 never snapshotted, since the plateau machinery
        ! never got past gate 2) now broadcasts actual_m_min (default 3) instead of the retired
        ! n_bins_candidates(1) global-pool estimate.
        call assert_equal_int(n_bins_per_point(1), 3_int32, &
                              "test_param_search_no_plateau_falls_back_to_finest: "// &
                              "n_bins_per_point falls back to m_min (default 3)")
    end subroutine test_param_search_no_plateau_falls_back_to_finest

    !> With `max_n_genes_all_studies=100`, the GAMMA-decay grid collapses to exactly ONE candidate
    !| (`floor(100/(0.25*300))=1` and `max(1,floor(100/(0.5*300)))=max(1,0)=1` are the same value, so
    !| the second `KX_FACTORS` entry does not add a distinct candidate, and a single GAMMA step
    !| already drops `n_points_high` below `n_points_low=300`): `(n_points, n_neighbors) = (300, 1)`. Both
    !| admissibility gates are relaxed to their most permissive settings
    !| (`min_neighbor_overlap=0.0`, `min_residuals_per_bin=0`) so the sole candidate exercises the
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
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 5_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
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
    !| (`min_residuals_per_bin=0`, `min_neighbor_overlap=0.0`) -- not the impossible values
    !| `test_param_search_no_plateau_falls_back_to_finest` uses to deliberately force every
    !| candidate to fail a gate.
    subroutine test_param_search_finds_plateau_mid_grid()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 10000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
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
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
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
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
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
    !| smaller confidence-interval width than every other candidate, in both studies. The routine
    !| must therefore return candidate 5's own `(n_points, n_neighbors)` and its real, bootstrapped
    !| confidence interval (not `-1.0`), with `plateau_established = .false.` The exact
    !| confidence-interval values below come directly from running this fixture through the actual
    !| implementation (not hand-derived), since they depend on the real bootstrap resampling --
    !| reproducible bit-for-bit given the fixed `random_seed` and deterministic input data.
    !| Re-verified after Issue #187's per-point occupancy search replaced the old global-pool
    !| Sturges/Freedman-Diaconis bin-count broadcast this routine used to use: with
    !| `min_residuals_per_bin=0` here, every neighborhood's occupancy search is trivially satisfied
    !| at every candidate bin count up to `m_max` (default 120), so it always selects the ceiling
    !| `m_max=120` bins per point -- a real, different bin count from before, which is why the
    !| bootstrapped confidence-interval values below (script-verified via a temporary debug print,
    !| not hand-derived) differ from this test's pre-#187 values, even though the candidate
    !| selection itself (still candidate 5) does not. Re-verified again after Step 3 (per-neighborhood
    !| residual range): `m_max=120` bins is still selected everywhere, unchanged, but each point's
    !| histogram range is now its own percentile-derived `[R_low, R_high]` instead of the old
    !| dataset-wide symmetric range, which shifts the actual bin boundaries -- and therefore the
    !| bootstrapped CI values below -- even though candidate selection (still candidate 5) and every
    !| other structural assertion in this test are unaffected.
    subroutine test_param_search_no_plateau_uses_smallest_uncertainty()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
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
                                     [9.3534521490265630e-04_real64, 1.8753358853084691e-02_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_uses_smallest_uncertainty: study 1 CI")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), &
                                     [9.3509303187972945e-04_real64, 1.8663009912309732e-02_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_uses_smallest_uncertainty: study 2 CI")

    end subroutine test_param_search_no_plateau_uses_smallest_uncertainty

    !> Same LCG-residual fixture as `test_param_search_no_plateau_uses_smallest_uncertainty`
    !| (guaranteed no CI-overlap plateau across all 6 candidates), but with
    !| `plateau_mode=MODE_PLATEAU_EFFECT_SIZE` -- confirming Issue #178's own fallback ("retain the
    !| previously defined fallback of selecting the admissible parameter setting with the smallest
    !| bootstrap uncertainty") now applies to effect-size mode too, not just CI-overlap. Since
    !| effect size also never plateaus on this fixture, the routine must land on the exact same
    !| smallest-uncertainty candidate (candidate 5, `(362, 220)`) and the exact same bootstrapped CI
    !| values as `test_param_search_no_plateau_uses_smallest_uncertainty` -- the ranking computation
    !| doesn't depend on which plateau criterion selected the no-plateau branch, only the gate that
    !| used to restrict its use to CI-overlap, which is now gone.
    subroutine test_param_search_no_plateau_effect_size_smallest_uncertainty()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               plateau_mode=MODE_PLATEAU_EFFECT_SIZE, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_no_plateau_effect_size_smallest_uncertainty: ierr should be OK")
        call assert_false(plateau_established, &
                          "test_param_search_no_plateau_effect_size_smallest_uncertainty: "// &
                          "no candidate ever plateaus under effect-size mode on this fixture")
        call assert_equal_int(n_points, 362_int32, &
                              "test_param_search_no_plateau_effect_size_smallest_uncertainty: "// &
                              "smallest-uncertainty candidate's n_points, same as CI-overlap mode's own fallback")
        call assert_equal_int(n_neighbors, 220_int32, &
                              "test_param_search_no_plateau_effect_size_smallest_uncertainty: "// &
                              "smallest-uncertainty candidate's n_neighbors")
        call assert_true(all(best_candidate_pair_confidence_interval /= -1.0_real64), &
                         "test_param_search_no_plateau_effect_size_smallest_uncertainty: "// &
                         "a real confidence interval is returned, not the -1.0 sentinel")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), &
                                     [9.3534521490265630e-04_real64, 1.8753358853084691e-02_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_effect_size_smallest_uncertainty: study 1 CI")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), &
                                     [9.3509303187972945e-04_real64, 1.8663009912309732e-02_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_effect_size_smallest_uncertainty: study 2 CI")
    end subroutine test_param_search_no_plateau_effect_size_smallest_uncertainty

    !> Same fixture again, `plateau_mode=MODE_PLATEAU_BOTH`. Since CI overlap never plateaus here
    !| (confirmed by `test_param_search_no_plateau_uses_smallest_uncertainty`) and effect size never
    !| plateaus here either (confirmed by
    !| `test_param_search_no_plateau_effect_size_smallest_uncertainty`), BOTH mode's
    !| `ci_plateau_found .or. effect_size_plateau_found` is false for every candidate too -- so this
    !| must land on Issue #178's own smallest-uncertainty fallback exactly like the other two modes
    !| now do, not the old finest-resolution/-1.0 sentinel: the fallback no longer depends on
    !| `plateau_mode` at all, only on whether at least one candidate was ever admissible.
    subroutine test_param_search_no_plateau_both_uses_smallest_uncertainty()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=0_int32, min_neighbor_overlap=0.0_real64, &
                                               plateau_mode=MODE_PLATEAU_BOTH, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_no_plateau_both_uses_smallest_uncertainty: ierr should be OK")
        call assert_false(plateau_established, &
                          "test_param_search_no_plateau_both_uses_smallest_uncertainty: "// &
                          "neither criterion ever plateaus under this fixture")
        call assert_equal_int(n_points, 362_int32, &
                              "test_param_search_no_plateau_both_uses_smallest_uncertainty: "// &
                              "smallest-uncertainty candidate's n_points, same as the other 2 modes' own fallback")
        call assert_equal_int(n_neighbors, 220_int32, &
                              "test_param_search_no_plateau_both_uses_smallest_uncertainty: "// &
                              "smallest-uncertainty candidate's n_neighbors")
        call assert_true(all(best_candidate_pair_confidence_interval /= -1.0_real64), &
                         "test_param_search_no_plateau_both_uses_smallest_uncertainty: "// &
                         "a real confidence interval is returned, not the -1.0 sentinel")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), &
                                     [9.3534521490265630e-04_real64, 1.8753358853084691e-02_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_both_uses_smallest_uncertainty: study 1 CI")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), &
                                     [9.3509303187972945e-04_real64, 1.8663009912309732e-02_real64], 2_int32, TOL, &
                                     "test_param_search_no_plateau_both_uses_smallest_uncertainty: study 2 CI")
    end subroutine test_param_search_no_plateau_both_uses_smallest_uncertainty

    !> Issue #187's occupancy FAILURE policy, end to end through the parameter search: a genuinely
    !| degenerate neighborhood (not merely an artificially huge `min_residuals_per_bin`, contrast
    !| `test_param_search_no_plateau_falls_back_to_finest` above, which forces failure that way)
    !| must make gate 2 (`check_mean_pmf_min_counts_impl`) reject the candidate, so it is never
    !| counted as admissible.
    !|
    !| Same grid as `test_param_search_no_plateau_falls_back_to_finest`
    !| (`max_n_genes_all_studies=2000`, finest candidate `(300, 26)`), and `min_neighbor_overlap=0.0`
    !| so gate 1 always passes trivially (tie-extended full-range neighborhoods from constant gene
    !| means, exactly as that test's own comment explains). Unlike that test, `min_residuals_per_bin`
    !| is left at its DEFAULT (10) here -- the rejection instead comes from making every residual in
    !| the whole dataset exactly `0.0_real64`. A pooled neighborhood's residuals are then ALL identical,
    !| so at any bin count `M >= 2` every bin except the one containing `0.0` has occupancy exactly
    !| `0 < 10`, and `M=1` is never tried (the default `m_min=3`) -- so
    !| `determine_bin_count_occupancy_impl` cannot find any admissible `M` for ANY point, in ANY
    !| candidate: `occupancy_failed=.true.` throughout, `selected_n_bins=m_min=3`, and the resulting
    !| 3-bin mean pmf (also all mass in one bin) then fails gate 2's own `min_residuals_per_bin` check
    !| the same way, for every candidate in the grid -- not just the first one tested. This is a
    !| stronger, more literal reading of "genuine occupancy FAILURE" than an oversized threshold: the
    !| data itself cannot support ANY histogram resolution above 1 bin.
    !|
    !| `n_admissible_evaluated` must therefore be exactly `0` (no `trace_*` columns at all, per Issue
    !| #178's own "a candidate that failed either gate has no trace_* entry at all" convention), and
    !| the routine falls back to the finest-resolution candidate `(300, 26)` with the confidence
    !| interval reset to `-1.0` and `plateau_established=.false.` -- case 3 of the final
    !| candidate-selection block, exactly like `test_param_search_no_plateau_falls_back_to_finest`,
    !| but reached here via genuine data sparsity rather than an artificial threshold. Values below
    !| are script-verified against the actual implementation (not hand-derived beyond the reasoning
    !| above), via a temporary debug print removed before landing this test.
    subroutine test_param_search_occupancy_failure_rejects_candidate()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 2000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
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
                residuals(:, i_gene, i_study) = 0.0_real64
            end do
        end do

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 5_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, &
                                               plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_neighbor_overlap=0.0_real64, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_occupancy_failure_rejects_candidate: ierr should be OK")
        call assert_equal_int(n_admissible_evaluated, 0_int32, &
                              "test_param_search_occupancy_failure_rejects_candidate: "// &
                              "every candidate's genuine occupancy failure makes gate 2 reject it, "// &
                              "so none is ever counted admissible")
        call assert_false(plateau_established, &
                          "test_param_search_occupancy_failure_rejects_candidate: "// &
                          "no candidate was ever admissible, so no plateau could be established")
        call assert_equal_int(n_points, 300_int32, &
                              "test_param_search_occupancy_failure_rejects_candidate: "// &
                              "falls back to the finest-resolution n_points")
        call assert_equal_int(n_neighbors, 26_int32, &
                              "test_param_search_occupancy_failure_rejects_candidate: "// &
                              "falls back to the finest-resolution n_neighbors")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 1), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_occupancy_failure_rejects_candidate: "// &
                                     "study 1 CI reset to -1.0, since gate 2 rejected every candidate")
        call assert_equal_array_real(best_candidate_pair_confidence_interval(:, 2), [-1.0_real64, -1.0_real64], 2_int32, TOL, &
                                     "test_param_search_occupancy_failure_rejects_candidate: "// &
                                     "study 2 CI reset to -1.0")
    end subroutine test_param_search_occupancy_failure_rejects_candidate

    !> Issue #187's central claim -- "different neighborhoods may use different numbers of bins" --
    !| verified through `trace_selected_n_bins`, not just through `determine_bin_count_occupancy`
    !| directly (already covered by the `test_determine_bin_count_occupancy_*` suite below): a
    !| SINGLE candidate whose 300 reference points split into two clearly different residual
    !| dispersion regimes must end up with a genuinely non-constant `trace_selected_n_bins` column.
    !|
    !| `max_n_genes_all_studies=100` collapses the grid to exactly one candidate, `(300, 1)` --
    !| exactly `test_param_search_single_candidate_bypasses_plateau`'s own grid collapse, so with
    !| `n_neighbors=1` every reference point's pooled neighborhood is just ONE gene's residuals,
    !| repeated across `max_n_reps_all_studies=3` replicates and `n_studies=2` studies: exactly 6
    !| pooled residuals per point, regardless of which gene. Genes 1-50 ("clustered") all get the
    !| SAME residual, `0.0_real64`, on every replicate of every study. Genes 51-100 ("spread") get 6
    !| DISTINCT values laid out across `(rep, study)` -- `{-2.5, -1.5, -0.5, 0.5, 1.5, 2.5}`.
    !|
    !| `m_min=1` (overriding the default 3) and `min_residuals_per_bin=1` (overriding the default 10)
    !| are both necessary for this fixture: a clustered point's 6 identical values can never fill more
    !| than 1 bin without leaving another one empty, so at the default `m_min=3` it would be a genuine
    !| occupancy FAILURE (`test_param_search_occupancy_failure_rejects_candidate` above already covers
    !| that path) rather than a valid, admissible, small `M_j` -- `m_min=1` lets a clustered point
    !| legitimately select `M_j=1` instead. Hand-traced through
    !| [[tox_data_integration_js_comp_test_impl(module):determine_bin_count_occupancy_impl(interface)]]'s
    !| own geometric search (`gamma_occupancy` default 1.25); Step 3 makes both points' own
    !| `shared_residual_range_low`/`_high` their own 5th/95th-percentile-derived range, not the old
    !| dataset-wide `[-3, 3]`:
    !| - Clustered (all six residuals `0.0`): every percentile of a constant array is that same
    !|   constant, so `R_low = R_high = 0.0` -- the same zero-range degenerate fallback
    !|   (`bin_width=1.0`, single occupied bin) `histogram_bin_counts` already used before Step 3,
    !|   just reached by a different route (a derived, not caller-supplied, zero span). `M=1`: all
    !|   six land in the one bin, `min_occ=6>=1` -> valid. `M=2` next (`ceil(1.25*1)=2`): the
    !|   fallback `bin_width=1.0` does not depend on `M`, so both bins still resolve to the same
    !|   single occupied index -> `min_occ=0<1` -> invalid immediately. Stage 2 refines the (2,1)
    !|   interval, which is empty (nothing strictly between 1 and 2) -> `best_m=1`, unchanged from
    !|   before Step 3.
    !| - Spread (`{-2.5,-1.5,-0.5,0.5,1.5,2.5}`): `rank(0.05,6)=0.05*5+1=1.25` -> interpolate
    !|   value(1)=-2.5, value(2)=-1.5 at fraction 0.25 -> `R_low=-2.25`; `rank(0.95,6)=5.75` ->
    !|   interpolate value(5)=1.5, value(6)=2.5 at fraction 0.75 -> `R_high=2.25` (span `4.5`,
    !|   narrower than the old `[-3,3]`, so the two extreme values now clamp). Retracing the
    !|   geometric ladder against this narrower span: `M=1` (bin_width 4.5, one bin) -> `min_occ=6`;
    !|   `M=2` (bin_width 2.25) -> `[3,3]`, `min_occ=3`; `M=3` (bin_width 1.5) -> `[2,2,2]`,
    !|   `min_occ=2`; `M=4` (`ceil(1.25*3)=3.75->4`, bin_width 1.125) -> `[2,1,1,2]`, `min_occ=1` --
    !|   all four valid so far. `M=5` (`ceil(1.25*4)=5`, bin_width 0.9, boundaries
    !|   `[-2.25,-1.35),[-1.35,-0.45),[-0.45,0.45),[0.45,1.35),[1.35,2.25]`) -> the clamped/original
    !|   values `{-2.25,-1.5,-0.5,0.5,1.5,2.25}` fall `[2,1,0,1,2]` -- the MIDDLE bin
    !|   `[-0.45,0.45)` is empty (no residual value lands near zero), `min_occ=0<1` -> invalid.
    !|   Stage 2 refines the interval strictly between `m_valid=4` and `m_invalid=5`, which is
    !|   empty (no integer strictly between 4 and 5) -> `best_m=4`, a materially different result
    !|   from the pre-Step-3 `best_m=6` (the old wider `[-3,3]` span never created that empty
    !|   middle bin at any `M` up to 6).
    !| So `trace_selected_n_bins(:, 1)` must contain BOTH `1` (every point whose nearest gene is
    !| clustered) and `4` (every point whose nearest gene is spread) -- never a single constant value
    !| across all 300 rows. `min_neighbor_overlap=0.0` keeps gate 1 trivial (as in the single-candidate
    !| test above), and since every point's own selected `M_j` is, by construction, the LARGEST `M`
    !| its own pooled residuals actually support, gate 2 also passes for every point, so the single
    !| candidate is genuinely admissible (`n_admissible_evaluated=1`) -- script-verified against the
    !| actual implementation, matching this hand trace exactly.
    subroutine test_param_search_different_neighborhoods_different_m_j()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 100, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate, i_gene, i_study
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)

        do i_study = 1, n_studies
            do i_gene = 1, max_n_genes_all_studies
                gene_means(i_gene, i_study) = real(i_gene, real64)
            end do
        end do
        do i_gene = 1, 50
            residuals(:, i_gene, :) = 0.0_real64
        end do
        do i_gene = 51, max_n_genes_all_studies
            residuals(1, i_gene, 1) = -2.5_real64
            residuals(2, i_gene, 1) = -1.5_real64
            residuals(3, i_gene, 1) = -0.5_real64
            residuals(1, i_gene, 2) = 0.5_real64
            residuals(2, i_gene, 2) = 1.5_real64
            residuals(3, i_gene, 2) = 2.5_real64
        end do

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 5_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, &
                                               plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=1_int32, min_neighbor_overlap=0.0_real64, &
                                               m_min=1_int32, random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_different_neighborhoods_different_m_j: ierr should be OK")
        call assert_equal_int(n_admissible_evaluated, 1_int32, &
                              "test_param_search_different_neighborhoods_different_m_j: "// &
                              "the sole (collapsed-grid) candidate is admissible")
        call assert_equal_int(trace_n_points(1), 300_int32, &
                              "test_param_search_different_neighborhoods_different_m_j: trace_n_points(1)")
        call assert_true(logical(minval(trace_selected_n_bins(1:trace_n_points(1), 1)) /= &
                                 maxval(trace_selected_n_bins(1:trace_n_points(1), 1)), kind=c_bool), &
                         "test_param_search_different_neighborhoods_different_m_j: "// &
                         "trace_selected_n_bins(1:trace_n_points(1), 1) is NOT constant across its rows")
        call assert_equal_int(minval(trace_selected_n_bins(1:trace_n_points(1), 1)), 1_int32, &
                              "test_param_search_different_neighborhoods_different_m_j: "// &
                              "the clustered points' own M_j (hand-traced above)")
        call assert_equal_int(maxval(trace_selected_n_bins(1:trace_n_points(1), 1)), 4_int32, &
                              "test_param_search_different_neighborhoods_different_m_j: "// &
                              "the spread points' own M_j (hand-traced above)")
    end subroutine test_param_search_different_neighborhoods_different_m_j

    !> Issue #187's final-selection wiring: the routine's returned `n_bins_per_point` must equal the
    !| `trace_selected_n_bins` COLUMN of whichever candidate was ACTUALLY selected -- not the column
    !| of the first admissible candidate, and not a stale snapshot left over from an earlier
    !| candidate. `test_param_search_finds_plateau_mid_grid` and the no-plateau fallback tests above
    !| already confirm the SELECTED `(n_points, n_neighbors)` is correct; this test additionally
    !| confirms the per-point BIN COUNTS returned alongside it are correct too, in a fixture where the
    !| search genuinely evaluates more than one admissible candidate before stopping (so the snapshot
    !| this checks against is not trivially the only one ever written).
    !|
    !| Same LCG-residual fixture as `test_param_search_no_plateau_uses_smallest_uncertainty`
    !| (`max_n_genes_all_studies=20000`, deterministic pseudo-random residuals, `random_seed=1`), but
    !| with `min_residuals_per_bin=30` instead of that test's `0` -- `0` makes every point's occupancy
    !| search trivially run all the way to `m_max` (every `min_occ>=0` unconditionally), which would
    !| make `trace_selected_n_bins` constant everywhere and this test unable to distinguish "returned
    !| the right column" from "returned any column" (they would all read the same value). `30` is
    !| large enough, relative to this fixture's per-candidate pooled residual counts, to make the
    !| occupancy search actually bind: bin counts vary WITHIN a candidate's own column (real per-point
    !| dispersion differences in the pseudo-random residuals) -- script-verified via a temporary debug
    !| print (swept several thresholds; higher ones, e.g. `70` and above, turned out to make EVERY
    !| point in a column pick the identical bin count under Step 3's own local per-point range --
    !| every point's own pooled residual set is the same size for a fixed candidate regardless of
    !| which genes it pools, and this fixture's residuals are close enough to identically distributed
    !| that the occupancy search converges on one answer everywhere at a high threshold; `30` sits low
    !| enough to still show real point-to-point spread): `n_admissible_evaluated=2` (this fixture now
    !| plateaus one candidate earlier than it did before Step 3 -- a genuine consequence of the new
    !| per-point range changing the observed JSD/CI, not a test bug), `plateau_established=.true.`,
    !| selecting the SECOND admissible candidate `(566, 70)` -- not the first, so the running
    !| best-candidate bin-count snapshot (`tmp_best_n_bins_per_point`/`tmp_best_shared_residual_range_low`/
    !| `_high`, and `trace_selected_n_bins` column bookkeeping) is genuinely overwritten as the loop
    !| progresses, not written once and left alone. That column's own bin counts range over `{10, 12}`
    !| (not constant), while column 1's own range over `{23, 24}` -- both non-constant, but with
    !| clearly different values -- so a wiring bug that returned column 1's snapshot instead would
    !| still be caught by the final elementwise comparison below (exact per-point equality against
    !| column 2, not merely "non-constant"), not silently masked by both columns looking similar.
    subroutine test_param_search_final_n_bins_matches_selected_trace_column()
        integer(int32), parameter :: n_studies = 2, max_n_genes_all_studies = 20000, max_n_reps_all_studies = 3
        real(real64) :: gene_means(max_n_genes_all_studies, n_studies)
        real(real64) :: residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies)
        integer(int32) :: n_points, n_neighbors, ierr, n_admissible_evaluated
        integer(int32) :: max_n_points_candidate, max_n_neighbors_candidate
        integer(int32), allocatable :: n_bins_per_point(:), trace_selected_n_bins(:, :), trace_n_pooled_residuals(:, :)
        real(real64), allocatable :: shared_residual_range_low(:), shared_residual_range_high(:)
        real(real64), allocatable :: trace_shared_residual_range_low(:, :), trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: trace_min_bin_occupancy(:, :), trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: trace_sturges_bins(:, :), trace_fd_bins(:, :)
        real(real64), allocatable :: trace_mean_bin_occupancy(:, :)
        logical(c_bool), allocatable :: trace_occupancy_failed(:, :)
        real(real64) :: best_candidate_pair_confidence_interval(2, n_studies)
        logical(c_bool) :: plateau_established
        integer(int32) :: trace_n_points(16), trace_n_neighbors(16)
        real(real64) :: trace_global_js_divergence(n_studies, 16), trace_ci_lower(n_studies, 16), trace_ci_upper(n_studies, 16)
        real(real64) :: trace_ci_width(n_studies, 16), trace_ci_width_relative(n_studies, 16)
        real(real64) :: trace_delta(n_studies, 16), trace_delta_median(16), trace_delta_max(16)
        integer(int32) :: i_gene, i_study, i_rep, k, t, selected_column

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

        call calc_js_comp_test_candidate_bounds(max_n_genes_all_studies, max_n_points_candidate, max_n_neighbors_candidate)
        allocate (n_bins_per_point(max_n_points_candidate))
        allocate (trace_selected_n_bins(max_n_points_candidate, 16), trace_n_pooled_residuals(max_n_points_candidate, 16))
        allocate (trace_min_bin_occupancy(max_n_points_candidate, 16), trace_max_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_sturges_bins(max_n_points_candidate, 16), trace_fd_bins(max_n_points_candidate, 16))
        allocate (trace_mean_bin_occupancy(max_n_points_candidate, 16))
        allocate (trace_occupancy_failed(max_n_points_candidate, 16))
        allocate (shared_residual_range_low(max_n_points_candidate), shared_residual_range_high(max_n_points_candidate))
        allocate (trace_shared_residual_range_low(max_n_points_candidate, 16), &
                 trace_shared_residual_range_high(max_n_points_candidate, 16))

        call run_js_comp_test_parameter_search(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, gene_means, &
                                               residuals, 10_int32, METHOD_JOIN_MIN, &
                                               max_n_points_candidate, max_n_neighbors_candidate, n_points, n_neighbors, &
                                               n_bins_per_point, shared_residual_range_low, shared_residual_range_high, &
                                               best_candidate_pair_confidence_interval, &
                                               plateau_established, &
                                               n_admissible_evaluated, &
                                               trace_n_points, trace_n_neighbors, trace_global_js_divergence, trace_ci_lower, &
                                               trace_ci_upper, trace_ci_width, trace_ci_width_relative, trace_delta, &
                                               trace_delta_median, trace_delta_max, &
                                               trace_selected_n_bins, trace_occupancy_failed, trace_n_pooled_residuals, &
                                               trace_min_bin_occupancy, trace_mean_bin_occupancy, trace_max_bin_occupancy, &
                                               trace_sturges_bins, trace_fd_bins, &
                                               trace_shared_residual_range_low, trace_shared_residual_range_high, &
                                               ierr=ierr, &
                                               min_residuals_per_bin=30_int32, min_neighbor_overlap=0.0_real64, &
                                               random_seed=1_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_param_search_final_n_bins_matches_selected_trace_column: ierr should be OK")
        call assert_equal_int(n_admissible_evaluated, 2_int32, &
                              "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                              "two candidates evaluated before the plateau -- not the first")
        call assert_true(plateau_established, &
                         "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                         "a genuine CI-overlap plateau is found under this stricter min_residuals_per_bin")
        call assert_equal_int(n_points, 566_int32, &
                              "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                              "selects the SECOND admissible candidate's n_points, not the first")
        call assert_equal_int(n_neighbors, 70_int32, &
                              "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                              "selects the second admissible candidate's n_neighbors")

        ! Find the trace column matching the actually-selected (n_points, n_neighbors) pair, rather
        ! than assuming it is column n_admissible_evaluated -- self-verifying regardless of which
        ! column the routine actually stopped at.
        selected_column = -1_int32
        do t = 1, n_admissible_evaluated
            if (trace_n_points(t) == n_points .and. trace_n_neighbors(t) == n_neighbors) then
                selected_column = t
                exit
            end if
        end do
        call assert_equal_int(selected_column, 2_int32, &
                              "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                              "the selected candidate is trace column 2 (script-verified)")
        call assert_true(logical(minval(trace_selected_n_bins(1:n_points, selected_column)) /= &
                                 maxval(trace_selected_n_bins(1:n_points, selected_column)), kind=c_bool), &
                         "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                         "the selected column's own bin counts are not constant either (10 and 12)")
        call assert_true(all(n_bins_per_point(1:n_points) == trace_selected_n_bins(1:n_points, selected_column)), &
                         "test_param_search_final_n_bins_matches_selected_trace_column: "// &
                         "returned n_bins_per_point equals the ACTUALLY-selected candidate's trace column, "// &
                         "not a stale earlier snapshot")
    end subroutine test_param_search_final_n_bins_matches_selected_trace_column

    !> Residuals are every integer in [-60, 59] (120 values). Unlike the pre-Step-3 dataset-wide
    !| symmetric range, `shared_residual_range_low`/`_high` are now this neighborhood's own
    !| 5th/95th-percentile bounds: with `calc_percentile_rank(q, n) = q*(n-1)+1` and the sorted
    !| values `value(k) = k - 61` for rank `k = 1..120`,
    !|   rank(0.05, 120) = 0.05*119 + 1 = 6.95 -> interpolate value(6)=-55, value(7)=-54 at
    !|     fraction 0.95 -> R_low = -55 + 0.95*1 = -54.05
    !|   rank(0.95, 120) = 0.95*119 + 1 = 114.05 -> interpolate value(114)=53, value(115)=54 at
    !|     fraction 0.05 -> R_high = 53 + 0.05*1 = 53.05
    !| so the range [-54.05, 53.05] (span 107.1) already trims off the most extreme ~5% on each
    !| side -- narrower than the full [-60, 59] data span -- which is why the two edge bins below
    !| end up more populated than an interior one (clamping pushes the trimmed tails into them).
    !| Tracing the default geometric ladder (m_min=3, gamma_occupancy=1.25: 3 -> 4 -> 5 -> 7 -> 9 ->
    !| 12 -> ...) against min_residuals_per_bin=10:
    !|   M          3    4    5    7    9   12
    !|   min_occ   36   27   21   15   12    9
    !| M=9 is the last admissible rung (min_occ=12); M=12 fails (min_occ=9). Stage 2 then refines
    !| 10 and 11: M=10 has min_occ=10 (exactly at the threshold, passes), M=11 has min_occ=9
    !| (fails) -- so best_m becomes 10, the larger of the two survivors -- this fixture is reused,
    !| with the SAME hand-computed range, by test_occupancy_diagnostics_hand_computed below.
    subroutine test_occupancy_finds_valid_below_m_max()
        integer(int32), parameter :: n_residuals = 120
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        do i = 1, n_residuals
            residuals(i) = real(i - 61, real64) ! -60, -59, ..., 59
        end do

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_finds_valid_below_m_max: ierr should be OK")
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_finds_valid_below_m_max: occupancy should not fail")
        call assert_equal_real(shared_residual_range_low, -54.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_finds_valid_below_m_max: "// &
                               "shared_residual_range_low == 5th percentile of [-60,59]")
        call assert_equal_real(shared_residual_range_high, 53.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_finds_valid_below_m_max: "// &
                               "shared_residual_range_high == 95th percentile of [-60,59]")
        call assert_equal_int(selected_n_bins, 10_int32, &
                              "test_determine_bin_count_occupancy_finds_valid_below_m_max: "// &
                              "should select M=10, strictly between m_min=3 and m_max=120")
        call assert_equal_int(n_pooled_residuals, 120_int32, &
                              "test_determine_bin_count_occupancy_finds_valid_below_m_max: n_pooled_residuals should be 120")
        call assert_equal_int(min_bin_occupancy, 10_int32, &
                              "test_determine_bin_count_occupancy_finds_valid_below_m_max: min_bin_occupancy at M=10")
        call assert_equal_int(max_bin_occupancy, 17_int32, &
                              "test_determine_bin_count_occupancy_finds_valid_below_m_max: max_bin_occupancy at M=10")
        call assert_equal_real(mean_bin_occupancy, 12.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_finds_valid_below_m_max: mean_bin_occupancy == 120/10")
    end subroutine test_occupancy_finds_valid_below_m_max

    !> Residuals are every integer in [-1200, 1199] (2400 values). `shared_residual_range_low`/
    !| `_high` are the 5th/95th percentiles of this same data: `rank(0.05,2400)=0.05*2399+1=120.95`
    !| -> interpolate value(120)=-1081, value(121)=-1080 at fraction 0.95 -> R_low = -1081 +
    !| 0.95*1 = -1080.05; `rank(0.95,2400)=0.95*2399+1=2280.05` -> interpolate value(2280)=1079,
    !| value(2281)=1080 at fraction 0.05 -> R_high = 1079 + 0.05*1 = 1079.05. This range
    !| ([-1080.05, 1079.05], span 2159.1) is narrower than the full data span, so the two edge bins
    !| absorb the trimmed ~5%-per-side tails and are far more populated than an interior one --
    !| still comfortably above the default min_residuals_per_bin=10 at every M up to and including
    !| the default m_max=120, so occupancy holds all the way through the default geometric ladder's
    !| last rung. `selected_n_bins == m_max` is reachable ONLY via the early-return branch: stage 2
    !| only ever tests strictly less than m_invalid, and m_invalid can be at most m_max, so a value
    !| of exactly m_max could never come out of stage 2's refinement loop even if it ran -- this
    !| assertion is therefore only satisfiable by the early-return `if (trial_m == actual_m_max)`
    !| branch actually firing.
    subroutine test_occupancy_reaches_m_max_validly()
        integer(int32), parameter :: n_residuals = 2400
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        do i = 1, n_residuals
            residuals(i) = real(i - 1201, real64) ! -1200, -1199, ..., 1199
        end do

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_reaches_m_max_validly: ierr should be OK")
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_reaches_m_max_validly: occupancy should not fail")
        call assert_equal_real(shared_residual_range_low, -1080.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_reaches_m_max_validly: "// &
                               "shared_residual_range_low == 5th percentile of [-1200,1199]")
        call assert_equal_real(shared_residual_range_high, 1079.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_reaches_m_max_validly: "// &
                               "shared_residual_range_high == 95th percentile of [-1200,1199]")
        call assert_equal_int(selected_n_bins, 120_int32, &
                              "test_determine_bin_count_occupancy_reaches_m_max_validly: "// &
                              "should reach m_max=120, only possible via the early-return branch")
        call assert_equal_int(n_pooled_residuals, 2400_int32, &
                              "test_determine_bin_count_occupancy_reaches_m_max_validly: n_pooled_residuals should be 2400")
        call assert_equal_int(min_bin_occupancy, 18_int32, &
                              "test_determine_bin_count_occupancy_reaches_m_max_validly: min_bin_occupancy at M=120")
        call assert_equal_int(max_bin_occupancy, 138_int32, &
                              "test_determine_bin_count_occupancy_reaches_m_max_validly: "// &
                              "max_bin_occupancy at M=120 (the two trimmed-tail edge bins)")
        call assert_equal_real(mean_bin_occupancy, 20.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_reaches_m_max_validly: mean_bin_occupancy == 2400/120")
    end subroutine test_occupancy_reaches_m_max_validly

    !> Only 5 pooled residuals total, so with the default min_residuals_per_bin=10 no bin at any
    !| M can ever reach 10 (there are not even 10 residuals to spread across bins) -- occupancy
    !| fails already at the default m_min=3 itself, the FAILURE case from Issue #187's own
    !| pseudocode ("even the minimum resolution is unsupported"). `shared_residual_range_low`/
    !| `_high` are this fixture's own 5th/95th percentiles: sorted `[-9,-5,0,5,9]`,
    !| `rank(0.05,5)=0.05*4+1=1.2` -> interpolate value(1)=-9, value(2)=-5 at fraction 0.2 -> R_low
    !| = -9 + 0.2*4 = -8.2; `rank(0.95,5)=0.95*4+1=4.8` -> interpolate value(4)=5, value(5)=9 at
    !| fraction 0.8 -> R_high = 5 + 0.8*4 = 8.2. This still FAILS regardless of range, since the
    !| occupancy criterion is bounded by the residual *count*, not by where the bin boundaries
    !| happen to fall.
    subroutine test_occupancy_m_min_itself_invalid_failure()
        integer(int32), parameter :: n_residuals = 5
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        residuals = [-9.0_real64, -5.0_real64, 0.0_real64, 5.0_real64, 9.0_real64]

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: ierr should be OK")
        call assert_true(occupancy_failed, &
                         "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                         "5 residuals can never fill any bin to 10")
        call assert_equal_real(shared_residual_range_low, -8.2_real64, TOL, &
                               "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                               "shared_residual_range_low == 5th percentile of [-9,-5,0,5,9]")
        call assert_equal_real(shared_residual_range_high, 8.2_real64, TOL, &
                               "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                               "shared_residual_range_high == 95th percentile of [-9,-5,0,5,9]")
        call assert_equal_int(selected_n_bins, 3_int32, &
                              "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                              "selected_n_bins should still be set to m_min=3 on FAILURE")
        call assert_equal_int(n_pooled_residuals, 5_int32, &
                              "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: n_pooled_residuals should be 5")
        call assert_equal_int(min_bin_occupancy, 0_int32, &
                              "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                              "min_bin_occupancy should be 0 on FAILURE")
        call assert_equal_int(max_bin_occupancy, 0_int32, &
                              "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                              "max_bin_occupancy should be 0 on FAILURE")
        call assert_equal_real(mean_bin_occupancy, 0.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_m_min_itself_invalid_failure: "// &
                               "mean_bin_occupancy should be 0 on FAILURE")
    end subroutine test_occupancy_m_min_itself_invalid_failure

    !> Compliance-review fix (Issue #187 cleanup): `min_residuals_per_bin`'s `DM_MIN` used to
    !| disagree between this routine (`1`) and `run_js_comp_test_parameter_search` (`0`), even
    !| though both forward the same argument into the same occupancy search -- a caller passing
    !| `0` here was rejected while the same value was accepted (and deliberately relied upon by
    !| several of `run_js_comp_test_parameter_search`'s own tests, e.g.
    !| `test_param_search_single_candidate_bypasses_plateau`) by the sibling entry point. Fixed by
    !| loosening this routine's bound to `0` too, matching the already-correct sibling.
    !|
    !| Reuses `test_occupancy_m_min_itself_invalid_failure`'s exact 5-residual fixture
    !| (`[-9,-5,0,5,9]`, its own R_low/R_high = -8.2/8.2 -- see that test's docstring for the
    !| derivation), which FAILS at the default `min_residuals_per_bin=10` (no bin can ever reach 10
    !| with only 5 residuals) -- but with `min_residuals_per_bin=0_int32` passed explicitly, every
    !| candidate bin count is trivially admissible (`min_occ >= 0` always holds), so the search
    !| never finds an inadmissible candidate. `m_max=3_int32` (equal to the default `m_min`) bounds
    !| the search to its very first candidate via the early-return branch: at `M=3` over
    !| `[-8.2, 8.2]` (bin width 16.4/3 ~= 5.467), the 5 residuals land
    !| `[-9,-5]->bin1 (clamped), [0]->bin2, [5,9]->bin3 (clamped)`, i.e. counts `[2,1,2]` --
    !| numerically the same bin counts as the pre-Step-3 fixture had, since this particular data is
    !| symmetric enough that trimming its own tails to +-8.2 still splits the same way as the old
    !| +-10 range did.
    subroutine test_occupancy_min_residuals_per_bin_zero_accepted()
        integer(int32), parameter :: n_residuals = 5
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        residuals = [-9.0_real64, -5.0_real64, 0.0_real64, 5.0_real64, 9.0_real64]

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr, &
                                           min_residuals_per_bin=0_int32, m_max=3_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_occupancy_min_residuals_per_bin_zero_accepted: "// &
                              "min_residuals_per_bin=0 must be accepted, not rejected by DM_MIN")
        call assert_false(occupancy_failed, &
                          "test_occupancy_min_residuals_per_bin_zero_accepted: "// &
                          "min_residuals_per_bin=0 trivially satisfies the occupancy criterion everywhere")
        call assert_equal_real(shared_residual_range_low, -8.2_real64, TOL, &
                               "test_occupancy_min_residuals_per_bin_zero_accepted: "// &
                               "shared_residual_range_low == 5th percentile of [-9,-5,0,5,9]")
        call assert_equal_real(shared_residual_range_high, 8.2_real64, TOL, &
                               "test_occupancy_min_residuals_per_bin_zero_accepted: "// &
                               "shared_residual_range_high == 95th percentile of [-9,-5,0,5,9]")
        call assert_equal_int(selected_n_bins, 3_int32, &
                              "test_occupancy_min_residuals_per_bin_zero_accepted: "// &
                              "early return at m_max=3 (bounded to keep the search's first candidate deterministic)")
        call assert_equal_int(n_pooled_residuals, 5_int32, &
                              "test_occupancy_min_residuals_per_bin_zero_accepted: n_pooled_residuals should be 5")
        call assert_equal_int(min_bin_occupancy, 1_int32, &
                              "test_occupancy_min_residuals_per_bin_zero_accepted: min bin count at M=3 is 1 (bin 2)")
        call assert_equal_int(max_bin_occupancy, 2_int32, &
                              "test_occupancy_min_residuals_per_bin_zero_accepted: max bin count at M=3 is 2 (bins 1 and 3)")
        call assert_equal_real(mean_bin_occupancy, 5.0_real64/3.0_real64, TOL, &
                               "test_occupancy_min_residuals_per_bin_zero_accepted: mean_bin_occupancy == 5/3")
    end subroutine test_occupancy_min_residuals_per_bin_zero_accepted

    !> Residuals are every integer in [-130, 129] EXCEPT the 16 integers [51, 67) (244 residuals
    !| total: 260 - 16 -- the previous version of this fixture declared `n_residuals = 243`, one
    !| short of the true count and an out-of-bounds write into `residuals(244)`; fixed here as part
    !| of this test's own Step-3 rewrite, not a Step-3 behavior change).
    !|
    !| `shared_residual_range_low`/`_high` are this fixture's own 5th/95th percentiles, not the old
    !| dataset-wide `+-130`: sorted, `rank(0.05,244)=0.05*243+1=13.15` -> interpolate the 13th/14th
    !| smallest values at fraction 0.15 -> R_low = -117.85; `rank(0.95,244)=0.95*243+1=231.85` ->
    !| interpolate the 231st/232nd smallest values at fraction 0.85 -> R_high = 116.85 (span
    !| 234.7). Tracing the default geometric ladder (m_min=3, gamma_occupancy=1.25:
    !| 3 -> 4 -> 5 -> 7 -> 9 -> 12 -> 15 -> 19 -> 24 -> 30) against min_residuals_per_bin=1 (an
    !| explicit override, so "admissible" just means "no bin is completely empty"):
    !|   M         12   15   19   24   30
    !|   min_occ   11    3    1    1    0
    !| M=24 is the last admissible rung (min_occ=1); M=30 fails (min_occ=0, the [51,67) gap finally
    !| empties a bin at this resolution). Stage 2 then exhaustively tests 25..29 -- all FAIL
    !| (min_occ=0 at every one) -- so best_m stays at m_valid=24, strictly greater than the smaller
    !| admissible rungs below it. At M=24 the bin counts range from 1 (min) to 22 (max), summing to
    !| 244 as required.
    subroutine test_occupancy_refinement_picks_above_m_valid()
        integer(int32), parameter :: n_residuals = 244
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i, idx
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        idx = 0
        do i = -130, 129
            if (i < 51 .or. i >= 67) then
                idx = idx + 1
                residuals(idx) = real(i, real64)
            end if
        end do

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: ierr should be OK")
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: occupancy should not fail")
        call assert_equal_real(shared_residual_range_low, -117.85_real64, TOL, &
                               "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: "// &
                               "shared_residual_range_low == 5th percentile of the gap fixture")
        call assert_equal_real(shared_residual_range_high, 116.85_real64, TOL, &
                               "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: "// &
                               "shared_residual_range_high == 95th percentile of the gap fixture")
        call assert_equal_int(n_pooled_residuals, 244_int32, &
                              "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: n_pooled_residuals")
        call assert_equal_int(selected_n_bins, 24_int32, &
                              "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: "// &
                              "best_m should be 24 (last admissible geometric rung; M=25..29 all fail)")
        call assert_equal_int(min_bin_occupancy, 1_int32, &
                              "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: min_bin_occupancy at M=24")
        call assert_equal_int(max_bin_occupancy, 22_int32, &
                              "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: max_bin_occupancy at M=24")
        call assert_equal_real(mean_bin_occupancy, 244.0_real64/24.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_refinement_picks_above_m_valid: mean_bin_occupancy")
    end subroutine test_occupancy_refinement_picks_above_m_valid

    !> Residuals are every integer in [-48, 47] (96 values). `shared_residual_range_low`/`_high`
    !| are this fixture's own 5th/95th percentiles: `rank(0.05,96)=0.05*95+1=5.75` -> interpolate
    !| value(5)=-44, value(6)=-43 at fraction 0.75 -> R_low = -44 + 0.75 = -43.25;
    !| `rank(0.95,96)=0.95*95+1=91.25` -> interpolate value(91)=42, value(92)=43 at fraction 0.25
    !| -> R_high = 42 + 0.25 = 42.25 (span 85.5). With the default min_residuals_per_bin=10, stage
    !| 2 now actually finds an improvement over m_valid on this range (M=8 passes at min_occ=10,
    !| beating m_valid=7) -- the opposite of what this test is named for, so
    !| min_residuals_per_bin=9_int32 is passed explicitly instead, close to the default but chosen
    !| so stage 2 genuinely finds nothing better here. Tracing the default geometric ladder against
    !| min_residuals_per_bin=9:
    !|   M          3   4    5    7    9   12
    !|   min_occ   28  21   17   12    9    7
    !| M=9 is the last admissible rung (m_valid=9, min_occ=9, exactly at the threshold); M=12 fails
    !| (min_occ=7, m_invalid=12). Stage 2 then exhaustively tests 10 and 11 (min_occ 8 and 7
    !| respectively, both < 9), so neither beats m_valid -- best_m stays 9, demonstrating stage 2
    !| running and finding nothing better.
    subroutine test_occupancy_refinement_finds_nothing_above_m_valid()
        integer(int32), parameter :: n_residuals = 96
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        do i = 1, n_residuals
            residuals(i) = real(i - 49, real64) ! -48, -47, ..., 47
        end do

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=9_int32, &
                                           ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: ierr should be OK")
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                          "occupancy should not fail")
        call assert_equal_real(shared_residual_range_low, -43.25_real64, TOL, &
                               "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                               "shared_residual_range_low == 5th percentile of [-48,47]")
        call assert_equal_real(shared_residual_range_high, 42.25_real64, TOL, &
                               "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                               "shared_residual_range_high == 95th percentile of [-48,47]")
        call assert_equal_int(selected_n_bins, 9_int32, &
                              "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                              "best_m should stay at m_valid=9 -- neither M=10 nor M=11 beats it")
        call assert_equal_int(n_pooled_residuals, 96_int32, &
                              "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                              "n_pooled_residuals should be 96")
        call assert_equal_int(min_bin_occupancy, 9_int32, &
                              "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                              "min_bin_occupancy at M=9")
        call assert_equal_int(max_bin_occupancy, 15_int32, &
                              "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                              "max_bin_occupancy at M=9")
        call assert_equal_real(mean_bin_occupancy, 96.0_real64/9.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_refinement_finds_nothing_above_m_valid: "// &
                               "mean_bin_occupancy == 96/9")
    end subroutine test_occupancy_refinement_finds_nothing_above_m_valid

    !> Every pooled residual is NaN: n_pooled_residuals must come out 0, matching
    !| estimate_bin_count_impl's own all-NaN branch (sturges_bins=fd_bins=1), and occupancy must
    !| be reported as FAILURE with all three occupancy diagnostics AND the two range bounds
    !| zeroed, per the early-return branch that never gets as far as computing a percentile or
    !| testing any candidate M.
    subroutine test_occupancy_all_residuals_nan()
        integer(int32), parameter :: n_residuals = 4
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr
        real(real64) :: mean_bin_occupancy, nan_val, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        nan_val = ieee_value(1.0_real64, ieee_quiet_nan)
        residuals = nan_val

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: ierr should be OK")
        call assert_true(occupancy_failed, &
                         "test_determine_bin_count_occupancy_all_residuals_nan: all-NaN pool must FAIL")
        call assert_equal_int(selected_n_bins, 3_int32, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: selected_n_bins should be m_min=3")
        call assert_equal_real(shared_residual_range_low, 0.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_all_residuals_nan: "// &
                               "shared_residual_range_low should be 0 (all-NaN degenerate fallback)")
        call assert_equal_real(shared_residual_range_high, 0.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_all_residuals_nan: "// &
                               "shared_residual_range_high should be 0 (all-NaN degenerate fallback)")
        call assert_equal_int(n_pooled_residuals, 0_int32, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: n_pooled_residuals should be 0")
        call assert_equal_int(min_bin_occupancy, 0_int32, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: min_bin_occupancy should be 0")
        call assert_equal_int(max_bin_occupancy, 0_int32, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: max_bin_occupancy should be 0")
        call assert_equal_real(mean_bin_occupancy, 0.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_all_residuals_nan: mean_bin_occupancy should be 0")
        call assert_equal_int(sturges_bins, 1_int32, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: sturges_bins should fall back to 1")
        call assert_equal_int(fd_bins, 1_int32, &
                              "test_determine_bin_count_occupancy_all_residuals_nan: fd_bins should fall back to 1")
    end subroutine test_occupancy_all_residuals_nan

    !> gamma_occupancy=1.01 (close to 1) forces roughly 117 tiny geometric steps from m_min=3 up to
    !| m_max=120 instead of the default gamma=1.25's ~16 -- reusing
    !| test_occupancy_reaches_m_max_validly's own dense fixture (every integer in [-1200, 1199],
    !| own R_low/R_high = -1080.05/1079.05 -- see that test's docstring for the derivation --
    !| comfortably occupied at every M up to 120) so the ONLY thing that changes is how many times
    !| `next_m = max(trial_m + 1, next_m)`'s
    !| progress guarantee actually fires. Mathematically, for any gamma_occupancy > 1.0 (enforced
    !| by DM_MIN(above(1.0_real64))) and integer trial_m >= 1, `ceiling(gamma_occupancy*trial_m)`
    !| is already guaranteed > trial_m in exact arithmetic, so this specific gamma value cannot
    !| itself demonstrate the guard changing the answer -- but it exercises the exact same
    !| `next_m = max(trial_m + 1, next_m)` line roughly 117 consecutive times (for most trial_m in
    !| this run, ceiling(1.01*trial_m) equals exactly trial_m+1, so the guard IS the value that
    !| gets used at almost every step), which is precisely the kind of sustained exercise that
    !| would surface a regression (e.g. a future rewrite that replaced `ceiling` with `int`
    !| truncation would silently stall at trial_m=99, since int(1.01*99)=int(99.99)=99): a
    !| non-terminating or wrong-answer result here is the visible symptom such a regression would
    !| produce, even though this specific test cannot itself force the mathematically-impossible
    !| "rounds back down" case for a validated gamma_occupancy > 1.
    subroutine test_occupancy_geometric_step_guarantees_progress()
        integer(int32), parameter :: n_residuals = 2400
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        do i = 1, n_residuals
            residuals(i) = real(i - 1201, real64) ! -1200, -1199, ..., 1199
        end do

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, gamma_occupancy=1.01_real64, ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_geometric_step_guarantees_progress: ierr should be OK")
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_geometric_step_guarantees_progress: "// &
                          "occupancy should not fail")
        call assert_equal_int(selected_n_bins, 120_int32, &
                              "test_determine_bin_count_occupancy_geometric_step_guarantees_progress: "// &
                              "the search must still terminate and reach m_max=120, not stall below it")
        call assert_equal_int(min_bin_occupancy, 18_int32, &
                              "test_determine_bin_count_occupancy_geometric_step_guarantees_progress: min_bin_occupancy")
        call assert_equal_int(max_bin_occupancy, 138_int32, &
                              "test_determine_bin_count_occupancy_geometric_step_guarantees_progress: max_bin_occupancy")
    end subroutine test_occupancy_geometric_step_guarantees_progress

    !> Reuses test_occupancy_finds_valid_below_m_max's exact residual fixture (every integer in
    !| [-60, 59], own R_low/R_high = -54.05/53.05 -- see that test's docstring for the derivation),
    !| but with max_n_reps_all_studies=8, n_neighbors=1 so estimate_bin_count_impl's own
    !| diagnostics are non-trivial (n_reps_neighborhood=8). `estimate_bin_count_impl`'s own quartile
    !| arithmetic is unchanged by Step 3 -- it still runs directly on the raw pooled residuals, not
    !| on R_low/R_high -- so its own inputs are exactly as before: sorted ascending, the 120
    !| residuals are `value(i) = i - 61` for rank `i = 1..120`. With
    !| `calc_percentile_rank(q, n) = q*(n-1)+1`:
    !|   rank(0.25, 120) = 0.25*119 + 1 = 30.75 -> interpolate value(30)=-31, value(31)=-30 at
    !|     fraction 0.75 -> quartile_25 = -31 + 0.75*1 = -30.25
    !|   rank(0.75, 120) = 0.75*119 + 1 = 90.25 -> interpolate value(90)=29, value(91)=30 at
    !|     fraction 0.25 -> quartile_75 = 29 + 0.25*1 = 29.25
    !| IQR = 29.25 - (-30.25) = 59.5. half_bin_width = 59.5 / 8**(1/3) = 59.5/2.0 = 29.75. What DOES
    !| change under Step 3 is what `estimate_bin_count_impl` is fed for its own one-sided
    !| `shared_residual_range` argument: half the new asymmetric span, `(53.05 - (-54.05))/2 =
    !| 53.55` (Step 3's own `half_span` local in `determine_bin_count_occupancy_impl`), not the old
    !| dataset-wide `60.0`. fd_raw = nint(53.55/29.75) = nint(1.8...) = 2 -> fd_bins = 2 (same
    !| result as before Step 3, since 53.55 and 60.0 both round to fd_raw=2 against this
    !| half_bin_width -- a coincidence of this particular fixture, not a general guarantee).
    !| sturges_raw = 1 + nint(log(8)/LOG_2) = 1 + nint(3.0) = 4 -> sturges_bins = 4 (untouched by
    !| Step 3 entirely -- Sturges never used shared_residual_range).
    !| These four diagnostics (sturges_bins, fd_bins, n_pooled_residuals, and the occupancy
    !| triple at whatever M the search itself selects) are independently derived here: the
    !| search's own selected_n_bins is asserted too (M=10, matching the sibling test above), but
    !| every diagnostic assertion below stands on its own arithmetic, not on that search result.
    subroutine test_occupancy_diagnostics_hand_computed()
        integer(int32), parameter :: n_residuals = 120
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        do i = 1, n_residuals
            residuals(i) = real(i - 61, real64) ! -60, -59, ..., 59
        end do

        call determine_bin_count_occupancy(residuals, n_residuals, 8_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: ierr should be OK")
        call assert_equal_int(n_pooled_residuals, 120_int32, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: n_pooled_residuals")
        call assert_equal_real(shared_residual_range_low, -54.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_diagnostics_hand_computed: shared_residual_range_low")
        call assert_equal_real(shared_residual_range_high, 53.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_diagnostics_hand_computed: shared_residual_range_high")
        call assert_equal_int(sturges_bins, 4_int32, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: sturges_bins")
        call assert_equal_int(fd_bins, 2_int32, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: fd_bins")
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_diagnostics_hand_computed: occupancy should not fail")
        call assert_equal_int(selected_n_bins, 10_int32, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: search should still pick M=10")
        call assert_equal_int(min_bin_occupancy, 10_int32, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: min_bin_occupancy at M=10")
        call assert_equal_int(max_bin_occupancy, 17_int32, &
                              "test_determine_bin_count_occupancy_diagnostics_hand_computed: max_bin_occupancy at M=10")
        call assert_equal_real(mean_bin_occupancy, 12.0_real64, TOL, &
                               "test_determine_bin_count_occupancy_diagnostics_hand_computed: mean_bin_occupancy == 120/10")
    end subroutine test_occupancy_diagnostics_hand_computed

    !> Four sub-cases, each isolating exactly one optional argument's resolved default by
    !| comparing a call with it absent against a call overriding only that one argument (every
    !| other optional either absent in both calls, or pinned to the SAME override in both, so the
    !| observed difference can only be attributed to the one argument actually being varied). Two
    !| of the four (B, C) needed genuinely new fixtures/parameters for Step 3: with a per-neighborhood
    !| percentile-derived range instead of a caller-supplied dataset-wide one, their old pre-Step-3
    !| fixtures no longer isolate the same default -- see each sub-case's own note below.
    !|
    !| (A) min_residuals_per_bin default=10: 20 residuals (every integer in [-10,9]), own
    !|     R_low/R_high = -9.05/8.05 (`rank(0.05,20)=0.05*19+1=1.95` -> interpolate value(1)=-10,
    !|     value(2)=-9 at fraction 0.95 -> R_low = -10+0.95 = -9.05; `rank(0.95,20)=0.95*19+1=19.05`
    !|     -> interpolate value(19)=8, value(20)=9 at fraction 0.05 -> R_high = 8+0.05 = 8.05), give
    !|     min_occ=6 at the default m_min=3 (bins [7,6,7]) -- 6 < 10 (default) FAILS, but 6 >= 5
    !|     (explicit override) PASSES.
    !| (B) m_min default=3: pre-Step-3 this fixture clustered 5 residuals in the extreme upper end
    !|     of a caller-supplied [-10,10] range to leave the lower bins of small M empty. Under
    !|     Step 3 the range is no longer caller-supplied -- it is derived FROM the clustered data
    !|     itself, so it never spans further than the data does, and the old fixture no longer
    !|     produces an empty bin at any M. Replaced with a fixture that keeps this property under a
    !|     percentile-derived range: 4 residuals tightly clustered near 0 (0.0, 0.1, 0.2, 0.3) plus
    !|     one extreme outlier (100.0) -- the outlier still stretches R_high far beyond the cluster
    !|     even though R_high is itself only the 95th percentile, not the outlier's own value.
    !|     Sorted, `rank(0.05,5)=1.2` -> R_low = 0.0 + 0.2*(0.1-0.0) = 0.02; `rank(0.95,5)=4.8` ->
    !|     R_high = 0.3 + 0.8*(100.0-0.3) = 80.06. At the default m_min=3 (bin_width =
    !|     80.04/3 ~= 26.68): the cluster's 4 points all land in bin 1, the outlier (clamped to
    !|     80.06) lands in bin 3, and bin 2 is completely empty -- FAILURE, with
    !|     min_residuals_per_bin=1 explicit (so "admissible" just means "no empty bin"). Starting
    !|     from an explicit m_min=1 instead: M=1 trivially passes (all 5 in one bin), M=2 also
    !|     passes (bin_width=40.02: the 4 clustered points and the clamped outlier split 4/1, no
    !|     empty bin), but M=3 fails the same way as above -- so the search stops advancing at
    !|     M=2, landing on selected_n_bins=2 (not 1: unlike the old fixture, M=2 is also genuinely
    !|     admissible here).
    !| (C) gamma_occupancy default=1.25: reuses test_occupancy_refinement_picks_above_m_valid's own
    !|     244-residual gap fixture (min_residuals_per_bin=1 explicit in both calls, to isolate
    !|     gamma alone; own R_low/R_high = -117.85/116.85, see that test's docstring). The default
    !|     ladder (3,4,5,7,9,12,15,19,24) gives selected_n_bins=24 (verified there). Pre-Step-3 this
    !|     sub-case contrasted gamma_occupancy=3.0's wider ladder (3,9,27) against the default --
    !|     but under the new percentile-derived range, that 3.0 override lands on the SAME answer
    !|     (24) as the default, no longer demonstrating a difference. An explicit
    !|     gamma_occupancy=1.5_real64 instead walks 3,5,8,12,18 (m_valid=12, m_invalid=18) and its
    !|     refinement window {13,...,17} finds selected_n_bins=17 -- a materially different answer,
    !|     confirming the default is genuinely 1.25, not 1.5.
    !| (D) m_max default=120: reuses test_occupancy_reaches_m_max_validly's own 2400-residual dense
    !|     fixture (own R_low/R_high = -1080.05/1079.05, see that test's docstring), comfortably
    !|     occupied (>= 10 per bin) at every M up to 120. Absent, the search reaches the default
    !|     m_max=120 (early return); with an explicit m_max=50 (still comfortably occupied), the
    !|     ladder instead clamps to and returns 50.
    subroutine test_occupancy_defaults_match_issue_suggestions()
        real(real64) :: residuals_a(20), residuals_b(5), residuals_c(244), residuals_d(2400)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr, i, idx
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        ! (A) min_residuals_per_bin default=10
        do i = 1, 20
            residuals_a(i) = real(i - 11, real64) ! -10, -9, ..., 9
        end do
        call determine_bin_count_occupancy(residuals_a, 20_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)
        call assert_true(occupancy_failed, &
                         "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                         "(A) default min_residuals_per_bin=10 should FAIL this fixture (min_occ=6 < 10)")
        call assert_equal_real(shared_residual_range_low, -9.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                               "(A) shared_residual_range_low == 5th percentile of [-10,9]")
        call assert_equal_real(shared_residual_range_high, 8.05_real64, TOL, &
                               "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                               "(A) shared_residual_range_high == 95th percentile of [-10,9]")
        call determine_bin_count_occupancy(residuals_a, 20_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=5_int32, &
                                           ierr=ierr)
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                          "(A) explicit min_residuals_per_bin=5 should PASS the same fixture (min_occ=6 >= 5)")

        ! (B) m_min default=3 -- cluster + one extreme outlier, see docstring above for why this
        ! fixture replaces the pre-Step-3 one
        residuals_b = [0.0_real64, 0.1_real64, 0.2_real64, 0.3_real64, 100.0_real64]
        call determine_bin_count_occupancy(residuals_b, 5_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           ierr=ierr)
        call assert_true(occupancy_failed, &
                         "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                         "(B) default m_min=3 should FAIL (M=3 has an empty middle bin)")
        call assert_equal_real(shared_residual_range_low, 0.02_real64, TOL, &
                               "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                               "(B) shared_residual_range_low == 5th percentile of the cluster+outlier fixture")
        call assert_equal_real(shared_residual_range_high, 80.06_real64, TOL, &
                               "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                               "(B) shared_residual_range_high == 95th percentile of the cluster+outlier fixture")
        call determine_bin_count_occupancy(residuals_b, 5_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           m_min=1_int32, ierr=ierr)
        call assert_false(occupancy_failed, &
                          "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                          "(B) explicit m_min=1 should PASS, growing as far as M=2")
        call assert_equal_int(selected_n_bins, 2_int32, &
                              "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                              "(B) selected_n_bins should be 2 with m_min=1 (M=1 and M=2 both admissible, M=3 fails)")

        ! (C) gamma_occupancy default=1.25
        idx = 0
        do i = -130, 129
            if (i < 51 .or. i >= 67) then
                idx = idx + 1
                residuals_c(idx) = real(i, real64)
            end if
        end do
        call determine_bin_count_occupancy(residuals_c, 244_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           ierr=ierr)
        call assert_equal_int(selected_n_bins, 24_int32, &
                              "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                              "(C) default gamma_occupancy=1.25 should select M=24")
        call determine_bin_count_occupancy(residuals_c, 244_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           gamma_occupancy=1.5_real64, ierr=ierr)
        call assert_equal_int(selected_n_bins, 17_int32, &
                              "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                              "(C) explicit gamma_occupancy=1.5 should instead select M=17")

        ! (D) m_max default=120
        do i = 1, 2400
            residuals_d(i) = real(i - 1201, real64) ! -1200, -1199, ..., 1199
        end do
        call determine_bin_count_occupancy(residuals_d, 2400_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, ierr=ierr)
        call assert_equal_int(selected_n_bins, 120_int32, &
                              "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                              "(D) default m_max=120 should be reached")
        call determine_bin_count_occupancy(residuals_d, 2400_int32, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, m_max=50_int32, ierr=ierr)
        call assert_equal_int(selected_n_bins, 50_int32, &
                              "test_determine_bin_count_occupancy_defaults_match_issue_suggestions: "// &
                              "(D) explicit m_max=50 should instead be reached")
    end subroutine test_occupancy_defaults_match_issue_suggestions

    !> Step 3 (per-neighborhood residual range): confirms shared_residual_range_low/_high are
    !| genuinely asymmetric -- NOT mirror images of each other -- for a residual distribution whose
    !| core is symmetric around zero but whose tail is not. 12 residuals: -5,-4,-3,-2,-1,0,1,2,3,4,5
    !| (symmetric) plus one extreme positive outlier, 50.0. Sorted, with
    !| `calc_percentile_rank(q,n) = q*(n-1)+1`:
    !|   rank(0.05, 12) = 0.05*11 + 1 = 1.55 -> interpolate value(1)=-5, value(2)=-4 at fraction
    !|     0.55 -> R_low = -5 + 0.55*1 = -4.45
    !|   rank(0.95, 12) = 0.95*11 + 1 = 11.45 -> interpolate value(11)=5, value(12)=50 at fraction
    !|     0.45 -> R_high = 5 + 0.45*45 = 25.25
    !| |R_low| = 4.45 and R_high = 25.25 differ by more than 20 -- clearly not mirror images of one
    !| another -- even though the 11 non-outlier residuals are themselves perfectly symmetric,
    !| demonstrating that a single extreme value on one side is enough to pull that side's bound far
    !| out while leaving the other essentially where a symmetric range would have put it. This is
    !| exactly the design Step 3 replaces the old single symmetric shared_residual_range with.
    subroutine test_occupancy_range_asymmetric_skewed_residuals()
        integer(int32), parameter :: n_residuals = 12
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        residuals = [-5.0_real64, -4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
                    1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 50.0_real64]

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_occupancy_range_asymmetric_skewed_residuals: ierr should be OK")
        call assert_equal_real(shared_residual_range_low, -4.45_real64, TOL, &
                               "test_occupancy_range_asymmetric_skewed_residuals: "// &
                               "shared_residual_range_low == 5th percentile of the skewed fixture")
        call assert_equal_real(shared_residual_range_high, 25.25_real64, TOL, &
                               "test_occupancy_range_asymmetric_skewed_residuals: "// &
                               "shared_residual_range_high == 95th percentile of the skewed fixture")
        call assert_true(abs(shared_residual_range_low) - shared_residual_range_high < -1.0_real64, &
                         "test_occupancy_range_asymmetric_skewed_residuals: "// &
                         "|shared_residual_range_low| and shared_residual_range_high must NOT be mirror images "// &
                         "(the old symmetric [-R,R] design would have made them equal)")
    end subroutine test_occupancy_range_asymmetric_skewed_residuals

    !> Step 3 (per-neighborhood residual range): a direct hand-computed check of
    !| shared_residual_range_low/_high against calc_percentile_impl's own nearest-rank/interpolation
    !| formula, independent of whatever bin count the occupancy search separately selects (that
    !| search's own output is not asserted here at all -- see the other test_occupancy_range_* and
    !| test_occupancy_* tests for that). 8 residuals, every integer 1..8 (chosen so both quantile
    !| ranks land on a clean fractional index):
    !|   rank(0.05, 8) = 0.05*7 + 1 = 1.35 -> interpolate value(1)=1, value(2)=2 at fraction 0.35
    !|     -> R_low = 1 + 0.35*1 = 1.35
    !|   rank(0.95, 8) = 0.95*7 + 1 = 7.65 -> interpolate value(7)=7, value(8)=8 at fraction 0.65
    !|     -> R_high = 7 + 0.65*1 = 7.65
    subroutine test_occupancy_range_hand_computed_percentile()
        integer(int32), parameter :: n_residuals = 8
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, n_pooled_residuals, min_bin_occupancy, max_bin_occupancy
        integer(int32) :: sturges_bins, fd_bins, ierr
        real(real64) :: mean_bin_occupancy, shared_residual_range_low, shared_residual_range_high
        logical(c_bool) :: occupancy_failed

        residuals = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64]

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, selected_n_bins, &
                                           occupancy_failed, shared_residual_range_low, shared_residual_range_high, &
                                           n_pooled_residuals, min_bin_occupancy, mean_bin_occupancy, &
                                           max_bin_occupancy, sturges_bins, fd_bins, min_residuals_per_bin=1_int32, &
                                           ierr=ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_occupancy_range_hand_computed_percentile: ierr should be OK")
        call assert_equal_real(shared_residual_range_low, 1.35_real64, TOL, &
                               "test_occupancy_range_hand_computed_percentile: "// &
                               "shared_residual_range_low == calc_percentile_impl(0.05) nearest-rank interpolation")
        call assert_equal_real(shared_residual_range_high, 7.65_real64, TOL, &
                               "test_occupancy_range_hand_computed_percentile: "// &
                               "shared_residual_range_high == calc_percentile_impl(0.95) nearest-rank interpolation")
    end subroutine test_occupancy_range_hand_computed_percentile

    !> Step 3's last missing piece: `determine_bin_count_occupancy_exhaustive` is a brute-force
    !| reference implementation of the same per-point bin-count search, testing every candidate `M`
    !| in `[m_min, m_max]` independently instead of stopping at the production routine's first
    !| occupancy failure. On "ordinary" data the two should agree: reuses 11 of the 12 existing
    !| `test_occupancy_*` fixtures' own residual data and non-default optional arguments verbatim
    !| (`test_occupancy_diagnostics_hand_computed` is skipped -- its residuals and range are
    !| numerically identical to `test_occupancy_finds_valid_below_m_max`'s own fixture, since
    !| `max_n_reps_all_studies`/`n_neighbors` never enter the occupancy search itself, only
    !| `estimate_bin_count_impl`'s diagnostic Sturges/FD calculation, so it would just repeat that
    !| same comparison). For each fixture: call `determine_bin_count_occupancy` first (as the
    !| original test already does) to get its own `R_low`/`R_high`/`n_pooled_residuals`, then feed
    !| those captured values directly into `determine_bin_count_occupancy_exhaustive` with the same
    !| `m_min`/`m_max`/`min_residuals_per_bin` (the exhaustive routine has no `gamma_occupancy`, so
    !| that one optional never carries over), and assert `selected_n_bins`/`occupancy_failed` agree.
    !| Independently re-verified via a Python port of both algorithms before writing the assertions
    !| below (not merely assumed from the production routine's own tests), since the production
    !| routine's own Stage 2 refinement window only ever tests integers strictly between `m_valid`
    !| and the first-found `m_invalid` -- it never looks past `m_invalid` at all, so agreement here
    !| is a genuine property of these fixtures, not a logical certainty.
    subroutine test_occupancy_exhaustive_matches_production_on_all_fixtures()
        real(real64) :: r_a(120), r_b(2400), r_c(5), r_d(244), r_e(96), r_f(4), r_g(20), r_h(12), r_i(8)
        integer(int32) :: i, idx
        integer(int32) :: prod_n_bins, exh_n_bins, n_pooled, ierr1, ierr2
        integer(int32) :: min_occ, max_occ, sturges_bins, fd_bins
        real(real64) :: r_low, r_high, mean_occ
        logical(c_bool) :: prod_failed, exh_failed

        do i = 1, 120
            r_a(i) = real(i - 61, real64) ! -60, -59, ..., 59
        end do
        do i = 1, 2400
            r_b(i) = real(i - 1201, real64) ! -1200, -1199, ..., 1199
        end do
        r_c = [-9.0_real64, -5.0_real64, 0.0_real64, 5.0_real64, 9.0_real64]
        idx = 0
        do i = -130, 129
            if (i < 51 .or. i >= 67) then
                idx = idx + 1
                r_d(idx) = real(i, real64)
            end if
        end do
        do i = 1, 96
            r_e(i) = real(i - 49, real64) ! -48, -47, ..., 47
        end do
        r_f = ieee_value(1.0_real64, ieee_quiet_nan)
        do i = 1, 20
            r_g(i) = real(i - 11, real64) ! -10, -9, ..., 9
        end do
        r_h = [-5.0_real64, -4.0_real64, -3.0_real64, -2.0_real64, -1.0_real64, 0.0_real64, &
              1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 50.0_real64]
        r_i = [1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64, 7.0_real64, 8.0_real64]

        ! 1. test_occupancy_finds_valid_below_m_max's own fixture, default args.
        call determine_bin_count_occupancy(r_a, 120_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1)
        call determine_bin_count_occupancy_exhaustive(r_a, 120_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 1 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 1 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 1 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 1 selected_n_bins agrees")

        ! 2. test_occupancy_reaches_m_max_validly's own fixture, default args.
        call determine_bin_count_occupancy(r_b, 2400_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1)
        call determine_bin_count_occupancy_exhaustive(r_b, 2400_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 2 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 2 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 2 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 2 selected_n_bins agrees")

        ! 3. test_occupancy_m_min_itself_invalid_failure's own fixture, default args.
        call determine_bin_count_occupancy(r_c, 5_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1)
        call determine_bin_count_occupancy_exhaustive(r_c, 5_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 3 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 3 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 3 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 3 selected_n_bins agrees")

        ! 4. test_occupancy_min_residuals_per_bin_zero_accepted's own fixture (same as 3),
        !    min_residuals_per_bin=0, m_max=3.
        call determine_bin_count_occupancy(r_c, 5_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           min_residuals_per_bin=0_int32, m_max=3_int32)
        call determine_bin_count_occupancy_exhaustive(r_c, 5_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       min_residuals_per_bin=0_int32, m_max=3_int32)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 4 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 4 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 4 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 4 selected_n_bins agrees")

        ! 5. test_occupancy_refinement_picks_above_m_valid's own fixture, min_residuals_per_bin=1.
        call determine_bin_count_occupancy(r_d, 244_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           min_residuals_per_bin=1_int32)
        call determine_bin_count_occupancy_exhaustive(r_d, 244_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       min_residuals_per_bin=1_int32)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 5 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 5 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 5 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 5 selected_n_bins agrees")

        ! 6. test_occupancy_refinement_finds_nothing_above_m_valid's own fixture, min_residuals_per_bin=9.
        call determine_bin_count_occupancy(r_e, 96_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           min_residuals_per_bin=9_int32)
        call determine_bin_count_occupancy_exhaustive(r_e, 96_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       min_residuals_per_bin=9_int32)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 6 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 6 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 6 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 6 selected_n_bins agrees")

        ! 7. test_occupancy_all_residuals_nan's own fixture, default args -- both routines must
        !    agree on the degenerate n_pooled_residuals=0 FAILURE path.
        call determine_bin_count_occupancy(r_f, 4_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1)
        call determine_bin_count_occupancy_exhaustive(r_f, 4_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 7 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 7 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 7 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 7 selected_n_bins agrees")

        ! 8. test_occupancy_geometric_step_guarantees_progress reuses fixture 2's own residuals with
        !    gamma_occupancy=1.01 -- gamma has no counterpart on the exhaustive routine (it tests
        !    every M regardless of any growth factor), so only m_min/m_max/min_residuals_per_bin
        !    (all default here) carry over.
        call determine_bin_count_occupancy(r_b, 2400_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           gamma_occupancy=1.01_real64)
        call determine_bin_count_occupancy_exhaustive(r_b, 2400_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 8 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 8 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 8 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 8 selected_n_bins agrees")

        ! (test_occupancy_diagnostics_hand_computed skipped -- see doc block above for why.)

        ! 9. test_occupancy_defaults_match_issue_suggestions's own sub-case (A) fixture, with the
        !    explicit min_residuals_per_bin=5 override (its second, passing call).
        call determine_bin_count_occupancy(r_g, 20_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           min_residuals_per_bin=5_int32)
        call determine_bin_count_occupancy_exhaustive(r_g, 20_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       min_residuals_per_bin=5_int32)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 9 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 9 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 9 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 9 selected_n_bins agrees")

        ! 10. test_occupancy_range_asymmetric_skewed_residuals's own fixture, min_residuals_per_bin=1.
        call determine_bin_count_occupancy(r_h, 12_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           min_residuals_per_bin=1_int32)
        call determine_bin_count_occupancy_exhaustive(r_h, 12_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       min_residuals_per_bin=1_int32)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 10 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 10 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 10 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 10 selected_n_bins agrees")

        ! 11. test_occupancy_range_hand_computed_percentile's own fixture, min_residuals_per_bin=1.
        call determine_bin_count_occupancy(r_i, 8_int32, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, ierr=ierr1, &
                                           min_residuals_per_bin=1_int32)
        call determine_bin_count_occupancy_exhaustive(r_i, 8_int32, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       min_residuals_per_bin=1_int32)
        call assert_equal_int(get_err_code(ierr1), ERR_OK, "test_occupancy_exhaustive_matches_production: 11 prod ierr")
        call assert_equal_int(get_err_code(ierr2), ERR_OK, "test_occupancy_exhaustive_matches_production: 11 exh ierr")
        call assert_true(prod_failed .eqv. exh_failed, "test_occupancy_exhaustive_matches_production: 11 occupancy_failed agrees")
        call assert_equal_int(exh_n_bins, prod_n_bins, "test_occupancy_exhaustive_matches_production: 11 selected_n_bins agrees")
    end subroutine test_occupancy_exhaustive_matches_production_on_all_fixtures

    !> The actual point of this routine's existence, demonstrated rather than just claimed: a
    !| hand-verified adversarial fixture where the production routine's geometric-search-then-
    !| refinement algorithm provably misses a larger admissible `M` its own ladder never reaches.
    !| Residuals `{0,5,10,16,22,28,29,40,46,52,56,60}` (12 values), range pinned to the data's exact
    !| min/max via `lower_residual_range_quantile=0.0`/`upper_residual_range_quantile=1.0` (so
    !| `R_low=0`, `R_high=60` exactly, not an interpolated percentile), `m_min=3`, `m_max=6`,
    !| `min_residuals_per_bin=2`.
    !|
    !| Direct enumeration (bin_width = 60/M): `M=3` -> counts `[4,3,5]`, `min_occ=3>=2` admissible.
    !| `M=4` (the production ladder's next rung, `ceil(1.25*3)=3.75->4`) -> counts `[3,4,1,4]`,
    !| `min_occ=1<2` -- FIRST FAILURE, so the production routine stops here: `m_valid=3`,
    !| `m_invalid=4`, the refinement interval strictly between them is empty, `selected_n_bins=3`.
    !| It never tries `M=5` or `M=6` at all. But `M=5` -> counts `[3,2,2,2,3]`, `min_occ=2>=2` --
    !| independently admissible, and LARGER than what production returned. (`M=6` -> counts
    !| `[2,2,3,0,2,3]`, `min_occ=0`, correctly inadmissible, so `5` is genuinely the exhaustive
    !| maximum, not just "a" larger admissible value found by luck.) The exhaustive search, which
    !| never stops early, correctly finds `5`.
    subroutine test_occupancy_exhaustive_diverges_on_adversarial_fixture()
        integer(int32), parameter :: n_residuals = 12
        real(real64) :: residuals(n_residuals)
        integer(int32) :: prod_n_bins, exh_n_bins, n_pooled, ierr1, ierr2
        integer(int32) :: min_occ, max_occ, sturges_bins, fd_bins
        real(real64) :: r_low, r_high, mean_occ
        logical(c_bool) :: prod_failed, exh_failed

        residuals = [0.0_real64, 5.0_real64, 10.0_real64, 16.0_real64, 22.0_real64, 28.0_real64, &
                    29.0_real64, 40.0_real64, 46.0_real64, 52.0_real64, 56.0_real64, 60.0_real64]

        call determine_bin_count_occupancy(residuals, n_residuals, 1_int32, 1_int32, prod_n_bins, prod_failed, &
                                           r_low, r_high, n_pooled, min_occ, mean_occ, max_occ, sturges_bins, fd_bins, &
                                           ierr=ierr1, m_min=3_int32, m_max=6_int32, min_residuals_per_bin=2_int32, &
                                           lower_residual_range_quantile=0.0_real64, upper_residual_range_quantile=1.0_real64)

        call assert_equal_int(get_err_code(ierr1), ERR_OK, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: prod ierr")
        call assert_equal_real(r_low, 0.0_real64, TOL, &
                               "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                               "R_low pinned to the data's exact min")
        call assert_equal_real(r_high, 60.0_real64, TOL, &
                               "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                               "R_high pinned to the data's exact max")
        call assert_equal_int(n_pooled, 12_int32, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: n_pooled_residuals")
        call assert_false(prod_failed, &
                          "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                          "production finds SOME admissible M (just not the largest)")
        call assert_equal_int(prod_n_bins, 3_int32, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                              "production stops at M=3, blind to M=5's own admissibility")

        call determine_bin_count_occupancy_exhaustive(residuals, n_residuals, n_pooled, r_low, r_high, exh_n_bins, &
                                                       exh_failed, min_occ, mean_occ, max_occ, ierr=ierr2, &
                                                       m_min=3_int32, m_max=6_int32, min_residuals_per_bin=2_int32)

        call assert_equal_int(get_err_code(ierr2), ERR_OK, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: exh ierr")
        call assert_false(exh_failed, &
                          "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                          "exhaustive search finds M=5 admissible")
        call assert_equal_int(exh_n_bins, 5_int32, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                              "exhaustive correctly finds the true maximum, M=5")
        call assert_equal_int(min_occ, 2_int32, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                              "min_bin_occupancy at M=5")
        call assert_equal_int(max_occ, 3_int32, &
                              "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                              "max_bin_occupancy at M=5")

        call assert_not_equal_int(exh_n_bins, prod_n_bins, &
                                  "test_occupancy_exhaustive_diverges_on_adversarial_fixture: "// &
                                  "this IS the intentional divergence this routine exists to catch, not a bug")
    end subroutine test_occupancy_exhaustive_diverges_on_adversarial_fixture

    !> A fresh, non-reused fixture exercising `determine_bin_count_occupancy_exhaustive` directly,
    !| with no dependence on `determine_bin_count_occupancy` at all: `R_low`/`R_high`/
    !| `n_pooled_residuals` are literal, hand-picked values, not derived from any prior call.
    !| 10 residuals `0.5, 1.5, ..., 9.5`, `R_low=0`, `R_high=10` (bin_width = 10/M), `m_min=2`,
    !| `m_max=5`, `min_residuals_per_bin=3`.
    !|
    !| Direct enumeration: `M=2` -> counts `[5,5]`, `min_occ=5` -- admissible, but this is the
    !| SMALLEST candidate tested, not the answer. `M=3` -> counts `[3,4,3]`, `min_occ=3>=3` --
    !| admissible. `M=4` -> counts `[2,3,2,3]`, `min_occ=2<3` -- inadmissible. `M=5` -> counts
    !| `[2,2,2,2,2]`, `min_occ=2<3` -- inadmissible, and this is the LARGEST candidate tested, not
    !| the answer either. The true answer is `M=3`: the largest bin count that is actually
    !| admissible, distinct from both "smallest passing" and "largest tested regardless of
    !| admissibility" -- specifically catches a candidate/scratch-array mixup bug class (e.g.
    !| accidentally reducing over the wrong loop variable, or reusing a stale `candidate_bin_counts`
    !| from a different `trial_m`), independent of whether `do concurrent` or a plain `do` was
    !| chosen for the search itself.
    subroutine test_occupancy_exhaustive_hand_computed_own_correctness()
        integer(int32), parameter :: n_residuals = 10
        real(real64) :: residuals(n_residuals)
        integer(int32) :: selected_n_bins, min_bin_occupancy, max_bin_occupancy, ierr
        real(real64) :: mean_bin_occupancy
        logical(c_bool) :: occupancy_failed

        residuals = [0.5_real64, 1.5_real64, 2.5_real64, 3.5_real64, 4.5_real64, 5.5_real64, &
                    6.5_real64, 7.5_real64, 8.5_real64, 9.5_real64]

        call determine_bin_count_occupancy_exhaustive(residuals, n_residuals, n_residuals, 0.0_real64, 10.0_real64, &
                                                       selected_n_bins, occupancy_failed, min_bin_occupancy, &
                                                       mean_bin_occupancy, max_bin_occupancy, ierr=ierr, &
                                                       m_min=2_int32, m_max=5_int32, min_residuals_per_bin=3_int32)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_occupancy_exhaustive_hand_computed_own_correctness: ierr should be OK")
        call assert_false(occupancy_failed, &
                          "test_occupancy_exhaustive_hand_computed_own_correctness: M=2 and M=3 are both admissible")
        call assert_equal_int(selected_n_bins, 3_int32, &
                              "test_occupancy_exhaustive_hand_computed_own_correctness: "// &
                              "true answer M=3 -- neither smallest-passing (2) nor largest-tested (5)")
        call assert_equal_int(min_bin_occupancy, 3_int32, &
                              "test_occupancy_exhaustive_hand_computed_own_correctness: min_bin_occupancy at M=3")
        call assert_equal_int(max_bin_occupancy, 4_int32, &
                              "test_occupancy_exhaustive_hand_computed_own_correctness: max_bin_occupancy at M=3")
        call assert_equal_real(mean_bin_occupancy, 10.0_real64/3.0_real64, TOL, &
                               "test_occupancy_exhaustive_hand_computed_own_correctness: mean_bin_occupancy == 10/3")
    end subroutine test_occupancy_exhaustive_hand_computed_own_correctness

    !> `gather_pooled_neighborhood_residuals` is the exact per-point pooling step
    !| `run_js_comp_test`/`run_js_comp_test_parameter_search` perform internally, now exported so a
    !| caller can reconstruct that same input on real data. 2 reps, 3 genes, 2 neighbors, 2 studies,
    !| every residual value distinct so a layout bug (wrong stride/order) cannot hide.
    !|
    !| `residuals(:, :, 1)` (study 1): gene 1 -> [1,2], gene 2 -> [3,4], gene 3 -> [5,6].
    !| `residuals(:, :, 2)` (study 2): gene 1 -> [10,20], gene 2 -> [30,40], gene 3 -> [50,60].
    !| `neighborhood_indices_point`: study 1 neighbors = genes [1, 3]; study 2 neighbors = genes
    !| [2, 1].
    !|
    !| Hand layout, per `n_predecessors = ((i_study-1)*n_neighbors + (i_neighbor-1))*max_n_reps`:
    !| study 1/neighbor 1 (gene 1) -> pooled(1:2) = [1, 2]; study 1/neighbor 2 (gene 3) ->
    !| pooled(3:4) = [5, 6]; study 2/neighbor 1 (gene 2) -> pooled(5:6) = [30, 40]; study
    !| 2/neighbor 2 (gene 1) -> pooled(7:8) = [10, 20].
    subroutine test_gather_pooled_residuals_hand_computed()
        integer(int32), parameter :: max_n_reps = 2, max_n_genes = 3, n_neighbors = 2, n_studies = 2
        real(real64) :: residuals(max_n_reps, max_n_genes, n_studies)
        integer(int32) :: neighborhood_indices_point(n_neighbors, n_studies)
        real(real64) :: pooled_residuals(max_n_reps*n_neighbors*n_studies)
        integer(int32) :: ierr

        residuals(:, :, 1) = reshape([1.0_real64, 2.0_real64, 3.0_real64, 4.0_real64, 5.0_real64, 6.0_real64], [2, 3])
        residuals(:, :, 2) = reshape([10.0_real64, 20.0_real64, 30.0_real64, 40.0_real64, 50.0_real64, 60.0_real64], &
                                     [2, 3])
        neighborhood_indices_point = reshape([1_int32, 3_int32, 2_int32, 1_int32], [2, 2])

        call gather_pooled_neighborhood_residuals(residuals, max_n_reps, max_n_genes, n_neighbors, n_studies, &
                                                   neighborhood_indices_point, pooled_residuals, ierr)

        call assert_equal_int(get_err_code(ierr), ERR_OK, &
                              "test_gather_pooled_residuals_hand_computed: ierr should be OK")
        call assert_equal_array_real(pooled_residuals, &
                                     [1.0_real64, 2.0_real64, 5.0_real64, 6.0_real64, &
                                      30.0_real64, 40.0_real64, 10.0_real64, 20.0_real64], &
                                     max_n_reps*n_neighbors*n_studies, TOL, &
                                     "test_gather_pooled_residuals_hand_computed: pooled layout")
    end subroutine test_gather_pooled_residuals_hand_computed

    !> `n_studies=0` must be rejected by the shared `validate_dimension_size` check at `arg_pos=5`,
    !| before the gather body ever runs (the array contents below are irrelevant to this test).
    subroutine test_gather_pooled_residuals_rejects_empty_dimension()
        integer(int32), parameter :: max_n_reps = 2, max_n_genes = 3, n_neighbors = 2, n_studies = 2
        real(real64) :: residuals(max_n_reps, max_n_genes, n_studies)
        integer(int32) :: neighborhood_indices_point(n_neighbors, n_studies)
        real(real64) :: pooled_residuals(max_n_reps*n_neighbors*n_studies)
        integer(int32) :: ierr

        residuals = 0.0_real64
        neighborhood_indices_point = 1_int32

        call gather_pooled_neighborhood_residuals(residuals, max_n_reps, max_n_genes, n_neighbors, 0_int32, &
                                                   neighborhood_indices_point, pooled_residuals, ierr)

        call assert_err(ierr, ERR_EMPTY_INPUT, &
                        "test_gather_pooled_residuals_rejects_empty_dimension: n_studies=0 rejected", &
                        arg_pos=5_int32)
    end subroutine test_gather_pooled_residuals_rejects_empty_dimension

    !> A gene index of `max_n_genes_all_studies + 1` in `neighborhood_indices_point` is out of
    !| bounds for `residuals`'s own gene extent and must be rejected with `ERR_INVALID_INPUT` at
    !| `arg_pos=6`, not silently read out-of-bounds.
    subroutine test_gather_pooled_residuals_rejects_gene_index_oob()
        integer(int32), parameter :: max_n_reps = 2, max_n_genes = 3, n_neighbors = 2, n_studies = 2
        real(real64) :: residuals(max_n_reps, max_n_genes, n_studies)
        integer(int32) :: neighborhood_indices_point(n_neighbors, n_studies)
        real(real64) :: pooled_residuals(max_n_reps*n_neighbors*n_studies)
        integer(int32) :: ierr

        residuals = 0.0_real64
        neighborhood_indices_point = reshape([1_int32, max_n_genes + 1_int32, 2_int32, 1_int32], [2, 2])

        call gather_pooled_neighborhood_residuals(residuals, max_n_reps, max_n_genes, n_neighbors, n_studies, &
                                                   neighborhood_indices_point, pooled_residuals, ierr)

        call assert_err(ierr, ERR_INVALID_INPUT, &
                        "test_gather_pooled_residuals_rejects_gene_index_oob: "// &
                        "gene index max_n_genes+1 rejected", arg_pos=6_int32)
    end subroutine test_gather_pooled_residuals_rejects_gene_index_oob

end module mod_test_data_integration_js_comp_test
