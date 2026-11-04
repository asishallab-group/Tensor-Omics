!> Unit test suite for the adaptive study-label Permutation-Test of the JS-Comp-Test
!| (tox_data_integration_js_comp_test_permutation_adaptive_impl): one assignment of samples to studies, the whole test
!| with its observed value, null distribution and p-values, and the rerun of the adaptive parameter search.
module mod_test_data_integration_js_comp_test_permutation_adaptive
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use tox_errors
    use test_suite, only: test_case
    use tox_data_integration_preprocessing, only: construct_adaptive_neighborhoods
    use tox_data_integration_js_comp_test_adaptive_search, only: calc_adaptive_js_comp_test_bounds
    use tox_data_integration_js_comp_test_impl, only: MAX_N_BINS, MAX_CANDIDATE_PAIRS
    use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_SYMMETRIC
    use mod_test_data_integration_js_comp_test, only: assert_equal_array_real_nan_exact
    use tox_data_integration_js_comp_test_permutation_adaptive, only: run_js_comp_test_adaptive_on_sample_assignment, &
                                                                      run_js_comp_test_adaptive_permutation
    use tox_data_integration_js_comp_test_permutation_bounds, only: calc_adaptive_js_comp_test_permutation_bounds
    use tox_data_integration_js_comp_test_permutation_adaptive_impl, only: &
                                                                           run_js_comp_test_adaptive_on_sample_assignment_impl, &
                                                                           run_js_comp_test_adaptive_permutation_impl

    use mod_test_data_integration_js_comp_test_permutation_common, only: make_expression, reference_means_residuals, &
                                                                         reference_adaptive, &
                                                                         reference_adaptive_search, &
                                                                         redraw_assignments, identity, &
                                                                         equal_and_not_nan, K_START, N_NULL_P_VALUES, &
                                                                         K_STEP, K_MAX, MIN_NULL_ABOVE_MINIMUM, &
                                                                         option_set_t

    implicit none

contains

    !> Get array of all available tests.
    function get_all_tests_js_comp_test_permutation_adaptive() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(10))
        all_tests(1) = test_case("test_adaptive_identity_matches_plain_pipeline", &
                                 test_adaptive_identity_matches_plain_pipeline)
        all_tests(2) = test_case("test_adaptive_swap_and_within_study_permutation", &
                                 test_adaptive_swap_and_within_study_permutation)
        all_tests(3) = test_case("test_adaptive_on_assignment_rejects_invalid_input", &
                                 test_adaptive_on_assignment_rejects_invalid_input)
        all_tests(4) = test_case("test_adaptive_label_test_reproducible_rows_match_assignment_run", &
                                 test_adaptive_label_test_reproducible_rows_match_assignment_run)
        all_tests(5) = test_case("test_adaptive_label_test_null_and_batch", test_adaptive_label_test_null_and_batch)
        all_tests(6) = test_case("test_adaptive_rerun_parameter_search", test_adaptive_rerun_parameter_search)
        all_tests(7) = test_case("test_adaptive_label_test_zero_permutations_direct_call", &
                                 test_adaptive_label_test_zero_permutations_direct_call)
        all_tests(8) = test_case("test_adaptive_label_test_failed_permutations_bookkeeping", &
                                 test_adaptive_label_test_failed_permutations_bookkeeping)
        all_tests(9) = test_case("test_adaptive_label_test_rejects_invalid_input", &
                                 test_adaptive_label_test_rejects_invalid_input)
        all_tests(10) = test_case("test_adaptive_forwards_construction_options", &
                                  test_adaptive_forwards_construction_options)
    end function get_all_tests_js_comp_test_permutation_adaptive

    !> Test helper: assignment run adaptive through the plain wrapper, fixed knobs (`divergence_mode`
    !| forwarded, absent: the default).
    subroutine adaptive_on_assignment(expression, assignment, global_jsd, ierr, divergence_mode)
        real(real64), intent(in) :: expression(:, :, :)
            !! The data
        integer(int32), intent(in) :: assignment(:)
            !! The assignment
        real(real64), intent(out) :: global_jsd(:)
            !! Global JSD per study
        integer(int32), intent(out) :: ierr
            !! Error code
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to the routine
        integer(int32) :: used_ks, used_kst, used_km, n_points

        call run_js_comp_test_adaptive_on_sample_assignment(size(expression, 3, kind=int32), &
                                                            size(expression, 2, kind=int32), &
                                                            size(expression, 1, kind=int32), expression, &
                                                            size(assignment, kind=int32), assignment, K_START, K_STEP, &
                                                            K_MAX, global_jsd, used_ks, used_kst, used_km, n_points, &
                                                            divergence_mode=divergence_mode, ierr=ierr)
    end subroutine adaptive_on_assignment

    !> Test helper: calls run_js_comp_test_adaptive_on_sample_assignment_impl directly, every work array allocated here.
    subroutine direct_adaptive_on_assignment(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                             expression, n_samples, sample_assignment, k_start, k_step, k_max, &
                                             max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                             global_js_divergence, used_k_start, used_k_step, used_k_max, &
                                             n_points, ierr, rerun_parameter_search, n_bootstraps, random_seed, &
                                             divergence_mode)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: n_samples
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: k_start
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: k_step
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: k_max
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: used_k_start
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: used_k_step
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: used_k_max
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: n_points
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to run_js_comp_test_adaptive_on_sample_assignment_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        logical(c_bool), allocatable :: tmp_sample_mark(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_pooled_neighborhood_range(:, :)
        integer(int32), allocatable :: tmp_n_neighbors_per_point(:, :)
        integer(int32), allocatable :: tmp_stop_reason(:)
        real(real64), allocatable :: tmp_neighborhood_dispersion(:)
        real(real64), allocatable :: tmp_neighborhood_mad(:)
        integer(int32), allocatable :: tmp_point_neighborhood_indices(:, :)
        real(real64), allocatable :: tmp_neighbor_residuals(:, :)
        integer(int32), allocatable :: tmp_point_pool(:)
        integer(int32), allocatable :: tmp_point_draw(:)
        integer(int32), allocatable :: tmp_n_bins_per_point(:)
        real(real64), allocatable :: tmp_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_shared_residual_range_high(:)
        logical(c_bool), allocatable :: tmp_occupancy_failed(:)
        integer(int32), allocatable :: tmp_n_pooled_residuals(:)
        integer(int32), allocatable :: tmp_min_bin_occupancy(:)
        real(real64), allocatable :: tmp_mean_bin_occupancy(:)
        integer(int32), allocatable :: tmp_max_bin_occupancy(:)
        integer(int32), allocatable :: tmp_sturges_bins(:)
        integer(int32), allocatable :: tmp_fd_bins(:)
        real(real64), allocatable :: tmp_pmfs(:, :, :)
        integer(int32), allocatable :: tmp_counts(:, :, :)
        integer(int32), allocatable :: tmp_included_n_reps(:, :)
        real(real64), allocatable :: tmp_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_included_n_reps(:)
        real(real64), allocatable :: tmp_js_divergences(:, :)
        real(real64), allocatable :: tmp_weights(:, :)
        real(real64), allocatable :: tmp_global_js_divergence(:)
        real(real64), allocatable :: tmp_pooled_residuals(:)
        integer(int32), allocatable :: tmp_pooled_residuals_perm(:)
        integer(int32), allocatable :: tmp_residual_to_study_perm(:)
        integer(int32), allocatable :: tmp_bin_counts_search(:)
        integer(int32), allocatable :: tmp_study_block_ends(:)
        integer(int32), allocatable :: tmp_study_block_sizes(:)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_k_start(:)
        integer(int32), allocatable :: tmp_search_trace_k_step(:)
        integer(int32), allocatable :: tmp_search_trace_k_max(:)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        real(real64), allocatable :: tmp_search_trace_global_js_divergence(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_lower(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_upper(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width_relative(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_overlap(:, :)
        real(real64), allocatable :: tmp_search_trace_delta(:, :)
        real(real64), allocatable :: tmp_search_trace_delta_median(:)
        real(real64), allocatable :: tmp_search_trace_delta_max(:)
        integer(int32), allocatable :: tmp_search_trace_selected_n_bins(:, :)
        logical(c_bool), allocatable :: tmp_search_trace_occupancy_failed(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_pooled_residuals(:, :)
        integer(int32), allocatable :: tmp_search_trace_min_bin_occupancy(:, :)
        real(real64), allocatable :: tmp_search_trace_mean_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_sturges_bins(:, :)
        integer(int32), allocatable :: tmp_search_trace_fd_bins(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_low(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: tmp_search_candidate_k_start(:)
        integer(int32), allocatable :: tmp_search_candidate_n_points(:)
        integer(int32), allocatable :: tmp_search_candidate_status(:)
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
        integer(int32), allocatable :: tmp_delta_perm(:)
        real(real64), allocatable :: tmp_best_uncertainty_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_best_n_bins_per_point(:)
        integer(int32), allocatable :: tmp_best_uncertainty_n_bins_per_point(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_high(:)
        allocate (tmp_sample_slots(max_n_reps_all_studies*n_studies))
        allocate (tmp_sample_mark(n_samples))
        allocate (tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means(max_n_genes_all_studies, n_studies))
        allocate (tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        allocate (tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        allocate (tmp_stop_reason(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        allocate (tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        allocate (tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        allocate (tmp_point_pool(MAX_N_BINS))
        allocate (tmp_point_draw(MAX_N_BINS))
        allocate (tmp_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_occupancy_failed(max_n_points_candidate))
        allocate (tmp_n_pooled_residuals(max_n_points_candidate))
        allocate (tmp_min_bin_occupancy(max_n_points_candidate))
        allocate (tmp_mean_bin_occupancy(max_n_points_candidate))
        allocate (tmp_max_bin_occupancy(max_n_points_candidate))
        allocate (tmp_sturges_bins(max_n_points_candidate))
        allocate (tmp_fd_bins(max_n_points_candidate))
        allocate (tmp_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_counts(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_included_n_reps(max_n_points_candidate, n_studies))
        allocate (tmp_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        allocate (tmp_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_weights(max_n_points_candidate, n_studies))
        allocate (tmp_global_js_divergence(n_studies))
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_step(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_median(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_status(MAX_CANDIDATE_PAIRS))
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_delta_perm(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_adaptive_on_sample_assignment_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=sample_assignment, k_start=k_start, &
            k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, global_js_divergence=global_js_divergence, &
            used_k_start=used_k_start, used_k_step=used_k_step, used_k_max=used_k_max, n_points=n_points, &
            tmp_sample_slots=tmp_sample_slots, tmp_sample_mark=tmp_sample_mark, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
            tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
            tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
            tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
            tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_k_start=tmp_search_trace_k_start, tmp_search_trace_k_step=tmp_search_trace_k_step, &
            tmp_search_trace_k_max=tmp_search_trace_k_max, tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
            tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
            tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_search_candidate_k_start=tmp_search_candidate_k_start, &
            tmp_search_candidate_n_points=tmp_search_candidate_n_points, &
            tmp_search_candidate_status=tmp_search_candidate_status, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            divergence_mode=divergence_mode, ierr=ierr)
    end subroutine direct_adaptive_on_assignment

    !> Test helper: calls run_js_comp_test_adaptive_permutation_impl directly, every work array allocated here.
    subroutine direct_adaptive_label_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                          k_start, k_step, k_max, max_n_points_candidate, &
                                          n_bootstrapping_top_k_jsds, n_permutations, &
                                          global_js_divergence_observed, p_values, jsd_null, &
                                          n_failed_permutations, ierr, rerun_parameter_search, n_bootstraps, &
                                          random_seed, &
                                          divergence_mode)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: k_start
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: k_step
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: k_max
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in) :: n_permutations
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(out) :: n_failed_permutations
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to run_js_comp_test_adaptive_permutation_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        integer(int32), allocatable :: tmp_sample_assignment(:)
        real(real64), allocatable :: tmp_permuted_global_js_divergence(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_pooled_neighborhood_range(:, :)
        integer(int32), allocatable :: tmp_n_neighbors_per_point(:, :)
        integer(int32), allocatable :: tmp_stop_reason(:)
        real(real64), allocatable :: tmp_neighborhood_dispersion(:)
        real(real64), allocatable :: tmp_neighborhood_mad(:)
        integer(int32), allocatable :: tmp_point_neighborhood_indices(:, :)
        real(real64), allocatable :: tmp_neighbor_residuals(:, :)
        integer(int32), allocatable :: tmp_point_pool(:)
        integer(int32), allocatable :: tmp_point_draw(:)
        integer(int32), allocatable :: tmp_n_bins_per_point(:)
        real(real64), allocatable :: tmp_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_shared_residual_range_high(:)
        logical(c_bool), allocatable :: tmp_occupancy_failed(:)
        integer(int32), allocatable :: tmp_n_pooled_residuals(:)
        integer(int32), allocatable :: tmp_min_bin_occupancy(:)
        real(real64), allocatable :: tmp_mean_bin_occupancy(:)
        integer(int32), allocatable :: tmp_max_bin_occupancy(:)
        integer(int32), allocatable :: tmp_sturges_bins(:)
        integer(int32), allocatable :: tmp_fd_bins(:)
        real(real64), allocatable :: tmp_pmfs(:, :, :)
        integer(int32), allocatable :: tmp_counts(:, :, :)
        integer(int32), allocatable :: tmp_included_n_reps(:, :)
        real(real64), allocatable :: tmp_mean_pmf(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_counts(:, :)
        integer(int32), allocatable :: tmp_mean_pmf_included_n_reps(:)
        real(real64), allocatable :: tmp_js_divergences(:, :)
        real(real64), allocatable :: tmp_weights(:, :)
        real(real64), allocatable :: tmp_global_js_divergence(:)
        real(real64), allocatable :: tmp_pooled_residuals(:)
        integer(int32), allocatable :: tmp_pooled_residuals_perm(:)
        integer(int32), allocatable :: tmp_residual_to_study_perm(:)
        integer(int32), allocatable :: tmp_bin_counts_search(:)
        integer(int32), allocatable :: tmp_study_block_ends(:)
        integer(int32), allocatable :: tmp_study_block_sizes(:)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_k_start(:)
        integer(int32), allocatable :: tmp_search_trace_k_step(:)
        integer(int32), allocatable :: tmp_search_trace_k_max(:)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        real(real64), allocatable :: tmp_search_trace_global_js_divergence(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_lower(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_upper(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_width_relative(:, :)
        real(real64), allocatable :: tmp_search_trace_ci_overlap(:, :)
        real(real64), allocatable :: tmp_search_trace_delta(:, :)
        real(real64), allocatable :: tmp_search_trace_delta_median(:)
        real(real64), allocatable :: tmp_search_trace_delta_max(:)
        integer(int32), allocatable :: tmp_search_trace_selected_n_bins(:, :)
        logical(c_bool), allocatable :: tmp_search_trace_occupancy_failed(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_pooled_residuals(:, :)
        integer(int32), allocatable :: tmp_search_trace_min_bin_occupancy(:, :)
        real(real64), allocatable :: tmp_search_trace_mean_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_max_bin_occupancy(:, :)
        integer(int32), allocatable :: tmp_search_trace_sturges_bins(:, :)
        integer(int32), allocatable :: tmp_search_trace_fd_bins(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_low(:, :)
        real(real64), allocatable :: tmp_search_trace_shared_residual_range_high(:, :)
        integer(int32), allocatable :: tmp_search_candidate_k_start(:)
        integer(int32), allocatable :: tmp_search_candidate_n_points(:)
        integer(int32), allocatable :: tmp_search_candidate_status(:)
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
        integer(int32), allocatable :: tmp_delta_perm(:)
        real(real64), allocatable :: tmp_best_uncertainty_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_best_n_bins_per_point(:)
        integer(int32), allocatable :: tmp_best_uncertainty_n_bins_per_point(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_best_uncertainty_shared_residual_range_high(:)
        allocate (tmp_sample_slots(max_n_reps_all_studies*n_studies))
        allocate (tmp_sample_assignment(max_n_reps_all_studies*n_studies))
        allocate (tmp_permuted_global_js_divergence(n_studies))
        allocate (tmp_permuted_expression(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means(max_n_genes_all_studies, n_studies))
        allocate (tmp_residuals(max_n_reps_all_studies, max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_neighborhood_range(2, max_n_genes_all_studies*n_studies))
        allocate (tmp_n_neighbors_per_point(n_studies, max_n_genes_all_studies*n_studies))
        allocate (tmp_stop_reason(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_dispersion(max_n_genes_all_studies*n_studies))
        allocate (tmp_neighborhood_mad(max_n_genes_all_studies*n_studies))
        allocate (tmp_point_neighborhood_indices(max_n_genes_all_studies, n_studies))
        allocate (tmp_neighbor_residuals(max_n_reps_all_studies, max_n_genes_all_studies))
        allocate (tmp_point_pool(MAX_N_BINS))
        allocate (tmp_point_draw(MAX_N_BINS))
        allocate (tmp_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_occupancy_failed(max_n_points_candidate))
        allocate (tmp_n_pooled_residuals(max_n_points_candidate))
        allocate (tmp_min_bin_occupancy(max_n_points_candidate))
        allocate (tmp_mean_bin_occupancy(max_n_points_candidate))
        allocate (tmp_max_bin_occupancy(max_n_points_candidate))
        allocate (tmp_sturges_bins(max_n_points_candidate))
        allocate (tmp_fd_bins(max_n_points_candidate))
        allocate (tmp_pmfs(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_counts(max_n_points_candidate, MAX_N_BINS, n_studies))
        allocate (tmp_included_n_reps(max_n_points_candidate, n_studies))
        allocate (tmp_mean_pmf(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_counts(max_n_points_candidate, MAX_N_BINS))
        allocate (tmp_mean_pmf_included_n_reps(max_n_points_candidate))
        allocate (tmp_js_divergences(max_n_points_candidate, n_studies))
        allocate (tmp_weights(max_n_points_candidate, n_studies))
        allocate (tmp_global_js_divergence(n_studies))
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_genes_all_studies*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_step(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_k_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_global_js_divergence(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_lower(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_upper(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_width_relative(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_ci_overlap(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta(n_studies, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_median(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_delta_max(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_selected_n_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_occupancy_failed(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_pooled_residuals(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_min_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_mean_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_max_bin_occupancy(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_sturges_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_fd_bins(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_low(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_shared_residual_range_high(max_n_points_candidate, MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_k_start(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_candidate_status(MAX_CANDIDATE_PAIRS))
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_delta_perm(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_adaptive_permutation_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, k_start=k_start, k_step=k_step, k_max=k_max, &
            max_n_points_candidate=max_n_points_candidate, n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
            n_permutations=n_permutations, global_js_divergence_observed=global_js_divergence_observed, &
            p_values=p_values, jsd_null=jsd_null, n_failed_permutations=n_failed_permutations, &
            tmp_sample_slots=tmp_sample_slots, tmp_sample_assignment=tmp_sample_assignment, &
            tmp_permuted_global_js_divergence=tmp_permuted_global_js_divergence, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, &
            tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
            tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
            tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, &
            tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
            tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
            tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
            tmp_global_js_divergence=tmp_global_js_divergence, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_k_start=tmp_search_trace_k_start, tmp_search_trace_k_step=tmp_search_trace_k_step, &
            tmp_search_trace_k_max=tmp_search_trace_k_max, tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, &
            tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, &
            tmp_search_trace_delta=tmp_search_trace_delta, &
            tmp_search_trace_delta_median=tmp_search_trace_delta_median, &
            tmp_search_trace_delta_max=tmp_search_trace_delta_max, &
            tmp_search_trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
            tmp_search_trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
            tmp_search_trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
            tmp_search_trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
            tmp_search_trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
            tmp_search_trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
            tmp_search_trace_sturges_bins=tmp_search_trace_sturges_bins, &
            tmp_search_trace_fd_bins=tmp_search_trace_fd_bins, &
            tmp_search_trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
            tmp_search_trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
            tmp_search_candidate_k_start=tmp_search_candidate_k_start, &
            tmp_search_candidate_n_points=tmp_search_candidate_n_points, &
            tmp_search_candidate_status=tmp_search_candidate_status, &
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            divergence_mode=divergence_mode, ierr=ierr)
    end subroutine direct_adaptive_label_test

    !> The adaptive identity assignment reproduces the plain construction + run exactly.
    subroutine test_adaptive_identity_matches_plain_pipeline()
        character(*), parameter :: name = "test_adaptive_identity_matches_plain_pipeline: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: reference(3), result(3)
        integer(int32) :: ierr, n_points_reference, used_ks, used_kst, used_km, n_points, n_samples

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 12_int32, &
                             expression)
        expression(:, 7, 2) = ieee_value(1.0_real64, ieee_quiet_nan)
        n_samples = 15_int32

        call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points_reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")
        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, n_samples, &
                                                            identity(n_samples), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call assert_equal_array_real_nan_exact(result, reference, size(reference, kind=int32), name//"identity equals the plain pipeline exactly")
        call assert_no_nan_real(result, 3_int32, name//"no NaN")
        call assert_equal_int(n_points, n_points_reference, name//"n_points")
        call assert_equal_int(used_ks, K_START, name//"used_k_start")
        call assert_equal_int(used_kst, K_STEP, name//"used_k_step")
        call assert_equal_int(used_km, K_MAX, name//"used_k_max")
    end subroutine test_adaptive_identity_matches_plain_pipeline

    !> The adaptive counterpart of the two hand-checkable cases: a swap swaps the values, a
    !| within-study reorder keeps them. The swap uses continuous data: the construction sorts all
    !| studies' means pooled, and dyadic means would tie across studies, where the order among tied
    !| entries follows their (swapped) pooled positions. The swap passes the symmetric divergence
    !| explicitly: the default pairwise divergence gives two studies one shared value.
    subroutine test_adaptive_swap_and_within_study_permutation()
        character(*), parameter :: name = "test_adaptive_swap_and_within_study_permutation: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), swapped(2), permuted(2)
        integer(int32) :: ierr

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 23_int32, expression)
        call adaptive_on_assignment(expression, identity(8_int32), observed, ierr, divergence_mode=MODE_DIVERGENCE_SYMMETRIC)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"observed ierr")
        call adaptive_on_assignment(expression, [5_int32, 6_int32, 7_int32, 8_int32, 1_int32, 2_int32, 3_int32, 4_int32], &
                                    swapped, ierr, divergence_mode=MODE_DIVERGENCE_SYMMETRIC)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"swapped ierr")
        call assert_true(observed(1) /= observed(2), name//"the two studies differ")
        call assert_true(swapped(1) == observed(2), name//"J1' == J2")
        call assert_true(swapped(2) == observed(1), name//"J2' == J1")

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 24_int32, expression, dyadic=.true._c_bool)
        call adaptive_on_assignment(expression, identity(8_int32), observed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"dyadic observed ierr")
        call adaptive_on_assignment(expression, [2_int32, 4_int32, 1_int32, 3_int32, 6_int32, 5_int32, 8_int32, 7_int32], &
                                    permuted, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"permuted ierr")
        call assert_equal_array_real_nan_exact(permuted, observed, size(observed, kind=int32), name//"within-study reorder keeps the values")
    end subroutine test_adaptive_swap_and_within_study_permutation

    !> Invalid inputs of the adaptive assignment run, and a construction with no neighborhood
    !| (`k_start` beyond the data), more than the capacity, or a status other than ok: with
    !| `min_study_neighbors` above every study's share of a neighborhood, the construction builds
    !| neighborhoods within the capacity but reports
    !| `STAT_ADAPTIVE_CONSTRUCTION_EMPTY_STUDY_NEIGHBORHOOD` (checked directly here), which level
    !| (a) rejects with invalid input; so does whole test on the observed data.
    subroutine test_adaptive_on_assignment_rejects_invalid_input()
        character(*), parameter :: name = "test_adaptive_on_assignment_rejects_invalid_input: "
        integer(int32), parameter :: too_many_study_neighbors = 1000_int32
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64), allocatable :: x_star(:), dispersion(:), mads(:)
        integer(int32), allocatable :: ranges(:, :), n_per_point(:, :), stops(:)
        real(real64) :: result(3), observed(3), p(3), null(2, 3)
        integer(int32) :: ierr, used_ks, used_kst, used_km, n_points, cap, top_k, search_cap, max_nn, status, n_failed

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 32_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 3) = ieee_value(1.0_real64, ieee_quiet_nan)
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, K_START, max_n_points_candidate=cap, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)

        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                           [1_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                            9_int32, 10_int32, 12_int32, 12_int32], K_START, K_STEP, K_MAX, cap, top_k, &
                                           result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: repeated value", arg_pos=6_int32)
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 13_int32, identity(13_int32), K_START, &
                                           K_STEP, K_MAX, cap, top_k, result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: n_samples", arg_pos=5_int32)
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, empty_study, 8_int32, identity(8_int32), K_START, &
                                           K_STEP, K_MAX, cap, top_k, result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: empty study", arg_pos=4_int32)
        call calc_adaptive_js_comp_test_bounds(300_int32, 3_int32, search_cap, ierr)
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), K_START, &
                                           K_STEP, K_MAX, search_cap - 1_int32, top_k, result, used_ks, used_kst, used_km, &
                                           n_points, ierr, rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: capacity below the search's", arg_pos=10_int32)
        ! Capacity 1: the construction builds more neighborhoods.
        call direct_adaptive_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), K_START, &
                                           K_STEP, K_MAX, 1_int32, top_k, result, used_ks, used_kst, used_km, n_points, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: more neighborhoods than the capacity", arg_pos=0_int32)
        call assert_true(n_points > 1_int32, name//"n_points reports the construction's count")
        call assert_true(all(ieee_is_nan(result)), name//"result NaN on error")

        ! k_start beyond the 900 pooled means: the construction builds nothing (plain wrapper).
        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                                            identity(12_int32), 1000_int32, 10_int32, 1000_int32, result, &
                                                            used_ks, used_kst, used_km, n_points, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"no neighborhood built", arg_pos=0_int32)
        call assert_equal_int(n_points, 0_int32, name//"n_points 0")
        ! Plain wrapper: a repeated value.
        call adaptive_on_assignment(expression, [1_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                                 9_int32, 9_int32, 11_int32, 12_int32], result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"wrapper: repeated value", arg_pos=0_int32)

        ! A construction status other than ok, with the neighborhood count within the capacity.
        call reference_means_residuals(expression, gene_means, residuals)
        allocate (x_star(900), ranges(2, 900), n_per_point(3, 900), stops(900), dispersion(900), mads(900))
        call construct_adaptive_neighborhoods(3_int32, 300_int32, 4_int32, gene_means, residuals, K_START, K_STEP, K_MAX, &
                                              n_points, x_star, ranges, n_per_point, stops, dispersion, mads, max_nn, &
                                              status, min_study_neighbors=too_many_study_neighbors, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"construction ierr")
        call assert_equal_int(status, STAT_ADAPTIVE_CONSTRUCTION_EMPTY_STUDY_NEIGHBORHOOD, name//"construction status")
        call assert_in_range_int(n_points, 1_int32, cap, name//"neighborhood count within the capacity")
        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                                            identity(12_int32), K_START, K_STEP, K_MAX, result, &
                                                            used_ks, used_kst, used_km, n_points, &
                                                            min_study_neighbors=too_many_study_neighbors, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"construction status not ok rejected", arg_pos=0_int32)
        call assert_true(all(ieee_is_nan(result)), name//"construction status: result NaN")
        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                                            identity(12_int32), K_START, K_STEP, K_MAX, result, &
                                                            used_ks, used_kst, used_km, n_points, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"the same call with the default min_study_neighbors")
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 4_int32, expression, K_START, K_STEP, K_MAX, &
                                                   2_int32, observed, p, null, n_failed, &
                                                   min_study_neighbors=too_many_study_neighbors, ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"whole test: observed construction status not ok", arg_pos=0_int32)
        call assert_true(all(ieee_is_nan(observed)) .and. all(ieee_is_nan(p)), name//"whole test: NaN outputs")
    end subroutine test_adaptive_on_assignment_rejects_invalid_input

    !> The adaptive counterpart of the reproducibility and row checks (same seed, identical null rows;
    !| each row is assignment run on the replayed draw).
    subroutine test_adaptive_label_test_reproducible_rows_match_assignment_run()
        character(*), parameter :: name = "test_adaptive_label_test_reproducible_rows_match_assignment_run: "
        integer(int32), parameter :: n = 5_int32, seed = 8_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), observed2(3), p2(3), null2(n, 3), row(3)
        integer(int32) :: ierr, n_failed, n_failed2, assignments(15, n), i

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 42_int32, &
                             expression)
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 6_int32, expression, K_START, K_STEP, &
                                                   K_MAX, n, observed, p, null, n_failed, random_seed=seed, &
                                                   ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 6_int32, expression, K_START, K_STEP, &
                                                   K_MAX, n, observed2, p2, null2, n_failed2, random_seed=seed, &
                                                   ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"second run ierr")
        call assert_equal_array_real_nan_exact(p, p2, 3_int32, name//"reproducible p-values")
        call assert_equal_array_real_nan_exact(observed, observed2, 3_int32, name//"reproducible observed values")
        call assert_equal_array_real_nan_exact(reshape(null, [3*n]), reshape(null2, [3*n]), 3_int32*n, &
                                               name//"reproducible null rows")
        call assert_equal_int(n_failed, 0_int32, name//"no failed permutation")
        ! The p-values recounted from the null rows (no failed permutation, so n valid rows).
        do i = 1, 3
            call assert_equal_real(p(i), real(count(null(:, i) >= observed(i)) + 1, real64)/real(n + 1, real64), &
                                   1.0e-15_real64, name//"p recounted from the null rows")
        end do

        call adaptive_on_assignment(expression, identity(15_int32), row, ierr)
        call assert_equal_array_real_nan_exact(observed, row, size(row, kind=int32), name//"observed is assignment run on the identity")
        call redraw_assignments(seed, 15_int32, n, assignments)
        do i = 1, n
            call adaptive_on_assignment(expression, assignments(:, i), row, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"assignment run ierr")
            call assert_equal_array_real_nan_exact(null(i, :), row, size(row, kind=int32), &
                                                   name//"row equals assignment run on the drawn assignment")
        end do
    end subroutine test_adaptive_label_test_reproducible_rows_match_assignment_run

    !> The adaptive counterpart of the null/batch check.
    subroutine test_adaptive_label_test_null_and_batch()
        character(*), parameter :: name = "test_adaptive_label_test_null_and_batch: "
        integer(int32), parameter :: n = 19_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), p_min
        integer(int32) :: ierr, n_failed, i_round, n_above

        p_min = 1.0_real64/real(n + 1, real64)
        n_above = 0_int32
        do i_round = 1, 4
            call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 0.5_real64], &
                                 200_int32 + i_round, expression)
            call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, &
                                                       K_STEP, K_MAX, n, observed, p, null, n_failed, &
                                                       random_seed=i_round, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null ierr")
            ! Exact comparison: p and p_min are both (count + 1)/(n + 1) quotients computed alike.
            n_above = n_above + count(p > p_min, kind=int32)
        end do
        call assert_true(n_above >= MIN_NULL_ABOVE_MINIMUM, name//"null p-values mostly above the minimum (got "// &
                         n_above//" of "//N_NULL_P_VALUES//")")

        call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 2.0_real64], 209_int32, &
                             expression)
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, &
                                                   K_MAX, n, observed, p, null, n_failed, random_seed=5_int32, &
                                                   ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"batch ierr")
        call assert_true(p(3) == p_min, name//"batch study at the minimum p-value")
    end subroutine test_adaptive_label_test_null_and_batch

    !> Adaptive rerun option: on, the result equals the adaptive search on the identity data followed
    !| by the construction and run with the knobs it selected; off, the fixed-knob run.
    subroutine test_adaptive_rerun_parameter_search()
        character(*), parameter :: name = "test_adaptive_rerun_parameter_search: "
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: reference(3), result(3), observed(3), p(3), null(1, 3), null4(4, 3), row(3)
        integer(int32) :: ierr, used_ks, used_kst, used_km, n_points, ks, kst, km, n_points_reference, n_failed, i
        integer(int32) :: assignments(13, 4)

        call make_expression(300_int32, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], 52_int32, &
                             expression)
        call reference_means_residuals(expression, gene_means, residuals)
        call reference_adaptive_search(gene_means, residuals, 20_int32, ks, kst, km)
        call reference_adaptive(expression, ks, kst, km, reference, n_points_reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")

        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                            identity(13_int32), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, &
                                                            rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                                            ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"rerun ierr")
        call assert_equal_int(used_ks, ks, name//"used_k_start is the search's")
        call assert_equal_int(used_kst, kst, name//"used_k_step is the search's")
        call assert_equal_int(used_km, km, name//"used_k_max is the search's")
        call assert_equal_int(n_points, n_points_reference, name//"n_points")
        call assert_equal_array_real_nan_exact(result, reference, size(reference, kind=int32), name//"rerun equals search + construction + run")

        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                            identity(13_int32), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, &
                                                            rerun_parameter_search=.false._c_bool, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off ierr")
        call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points_reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off reference ierr")
        call assert_true(equal_and_not_nan(result, reference), name//"off equals the fixed-knob run")
        call assert_equal_int(used_ks, K_START, name//"off: used_k_start")

        call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                            identity(13_int32), K_START, K_STEP, K_MAX, result, used_ks, &
                                                            used_kst, used_km, n_points, &
                                                            rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                                            ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"second rerun ierr")
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, &
                                                   K_MAX, 1_int32, observed, p, null, n_failed, &
                                                   rerun_parameter_search=.true._c_bool, &
                                                   n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"label test with rerun ierr")
        call assert_equal_array_real_nan_exact(observed, result, size(result, kind=int32), name//"label test with rerun: observed value")

        ! The null rows of the shuffled data sets, each of which reruns the search: row i is assignment run with the rerun on
        ! the i-th drawn assignment (seed 5, the same bootstrap seed), or all NaN when that run is invalid input.
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, K_MAX, &
                                                   4_int32, observed, p, null4, n_failed, random_seed=5_int32, &
                                                   rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null rows with rerun: ierr")
        call redraw_assignments(5_int32, 13_int32, 4_int32, assignments)
        do i = 1, 4
            call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, &
                                                                assignments(:, i), K_START, K_STEP, K_MAX, row, used_ks, &
                                                                used_kst, used_km, n_points, random_seed=5_int32, &
                                                                rerun_parameter_search=.true._c_bool, &
                                                                n_bootstraps=20_int32, ierr=ierr)
            if (get_err_code(ierr) == ERR_OK) then
                call assert_equal_array_real_nan_exact(null4(i, :), row, size(row, kind=int32), &
                                                       name//"null row equals assignment run with the rerun")
            else
                call assert_err(ierr, ERR_INVALID_INPUT, name//"a failing assignment run run is invalid input")
                call assert_true(all(ieee_is_nan(null4(i, :))), name//"a failing assignment run run: the null row is NaN")
            end if
        end do
        call assert_true(any(null4(:, 1) /= observed(1) .or. ieee_is_nan(null4(:, 1))), &
                         name//"a shuffled data set gives other values than the observed one")
    end subroutine test_adaptive_rerun_parameter_search

    !> A direct implementation call with `n_permutations = 0` (adaptive): the observed value is
    !| computed, the p-values are NaN, nothing failed.
    subroutine test_adaptive_label_test_zero_permutations_direct_call()
        character(*), parameter :: name = "test_adaptive_label_test_zero_permutations_direct_call: "
        real(real64), allocatable :: expression(:, :, :), null(:, :)
        real(real64) :: observed(3), p(3), reference(3)
        integer(int32) :: ierr, n_failed, max_p, top_k, n_points

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 61_int32, &
                             expression)
        allocate (null(0, 3))
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, K_START, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call direct_adaptive_label_test(3_int32, 300_int32, 4_int32, expression, K_START, K_STEP, K_MAX, max_p, top_k, &
                                        0_int32, observed, p, null, n_failed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive ierr")
        call assert_true(all(ieee_is_nan(p)), name//"adaptive p-values NaN")
        call reference_adaptive(expression, K_START, K_STEP, K_MAX, reference, n_points, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"adaptive reference ierr")
        call assert_equal_array_real_nan_exact(observed, reference, size(reference, kind=int32), name//"adaptive observed")
    end subroutine test_adaptive_label_test_zero_permutations_direct_call

    !> Failed permutations, a natural case: every study has 2 of its 5 samples shifted by +20, so
    !| on the observed data every study's gene means are shifted alike. A shuffle that gives the
    !| studies unequal numbers of shifted samples separates their means, some adaptive neighborhood
    !| then holds no gene of some study, and the run rejects it with invalid input. Such a
    !| permutation's `jsd_null` row is NaN, it is counted in `n_failed_permutations`, and the
    !| p-values count only the valid rows; when every permutation fails, the p-values are NaN.
    subroutine test_adaptive_label_test_failed_permutations_bookkeeping()
        character(*), parameter :: name = "test_adaptive_label_test_failed_permutations_bookkeeping: "
        integer(int32), parameter :: n = 19_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), expected, row(3)
        integer(int32) :: ierr, n_failed, i_study, n_valid, n_at_least, i, assignments(15, n), i_seed
        logical(c_bool) :: is_failed(n), found
        real(real64) :: null_one(1, 3)

        call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 0.5_real64], 71_int32, &
                             expression)
        expression(1:2, :, :) = expression(1:2, :, :) + 20.0_real64
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, K_STEP, &
                                                   K_MAX, n, observed, p, null, n_failed, random_seed=3_int32, &
                                                   ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call assert_no_nan_real(observed, 3_int32, name//"observed data run")
        do i = 1, n
            is_failed(i) = all(ieee_is_nan(null(i, :)))
            call assert_true(is_failed(i) .or. .not. any(ieee_is_nan(null(i, :))), name//"a row is all NaN or none")
        end do
        call assert_equal_int(n_failed, count(is_failed, kind=int32), name//"n_failed counts the NaN rows")
        call assert_true(n_failed > 0_int32 .and. n_failed < n, name//"some but not all permutations failed (got "// &
                         n_failed//")")
        ! A failed row is exactly a assignment run call that fails with invalid input.
        call redraw_assignments(3_int32, 15_int32, n, assignments)
        do i = 1, n
            call adaptive_on_assignment(expression, assignments(:, i), row, ierr)
            if (is_failed(i)) then
                call assert_err(ierr, ERR_INVALID_INPUT, name//"failed row: assignment run rejects the assignment")
            else
                call assert_equal_array_real_nan_exact(null(i, :), row, size(row, kind=int32), &
                                                       name//"valid row equals assignment run")
            end if
        end do
        n_valid = count(.not. is_failed, kind=int32)
        do i_study = 1, 3
            n_at_least = count((.not. is_failed) .and. null(:, i_study) >= observed(i_study), kind=int32)
            expected = real(n_at_least + 1, real64)/real(n_valid + 1, real64)
            call assert_true(p(i_study) == expected, name//"p over the valid permutations only")
        end do
        ! Every permutation failed (n_valid = 0): one permutation, and the first seed whose single
        ! draw fails on this data. The p-values are then NaN, the observed values are not.
        found = .false.
        do i_seed = 1, 40
            call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 5_int32, expression, K_START, &
                                                       K_STEP, K_MAX, 1_int32, observed, p, null_one, &
                                                       n_failed, random_seed=i_seed, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"single permutation ierr")
            if (n_failed == 1_int32) then
                found = .true.
                exit
            end if
        end do
        call assert_true(found, name//"some seed makes the single permutation fail")
        call assert_true(all(ieee_is_nan(p)), name//"all permutations failed: p NaN")
        call assert_true(all(ieee_is_nan(null_one)), name//"all permutations failed: the null row NaN")
        call assert_no_nan_real(observed, 3_int32, name//"all permutations failed: observed kept")
    end subroutine test_adaptive_label_test_failed_permutations_bookkeeping

    !> Whole test of the adaptive Permutation-Test rejects what does not depend on the assignment with its own
    !| positions (a study without samples, a rerun capacity below the search's), and returns an error on the observed
    !| data at once with NaN outputs.
    subroutine test_adaptive_label_test_rejects_invalid_input()
        character(*), parameter :: name = "test_adaptive_label_test_rejects_invalid_input: "
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :)
        real(real64) :: observed(3), p(3), null(2, 3)
        integer(int32) :: ierr, n_failed, max_p, top_k, search_cap

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 81_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        call calc_adaptive_js_comp_test_permutation_bounds(300_int32, 3_int32, K_START, max_n_points_candidate=max_p, &
                                                           n_bootstrapping_top_k_jsds=top_k, ierr=ierr)
        call direct_adaptive_label_test(3_int32, 300_int32, 4_int32, empty_study, K_START, K_STEP, K_MAX, max_p, top_k, &
                                        2_int32, observed, p, null, n_failed, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"adaptive: empty study", arg_pos=4_int32)
        call calc_adaptive_js_comp_test_bounds(300_int32, 3_int32, search_cap, ierr)
        call direct_adaptive_label_test(3_int32, 300_int32, 4_int32, expression, K_START, K_STEP, K_MAX, &
                                        search_cap - 1_int32, top_k, 2_int32, observed, p, null, n_failed, ierr, &
                                        rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"adaptive: capacity", arg_pos=8_int32)

        ! The observed data fail (k_start beyond the data): returned at once, every output NaN.
        call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 4_int32, expression, 1000_int32, &
                                                   10_int32, 1000_int32, 2_int32, observed, p, null, n_failed, &
                                                   ierr=ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"observed data fail", arg_pos=0_int32)
        call assert_true(all(ieee_is_nan(observed)) .and. all(ieee_is_nan(p)) .and. all(ieee_is_nan(null)), &
                         name//"observed failure: NaN outputs")
        call assert_equal_int(n_failed, 0_int32, name//"observed failure: nothing counted")
    end subroutine test_adaptive_label_test_rejects_invalid_input

    !> The options of the adaptive neighborhood construction reach the run through both routines. For each of `tau`,
    !| `mad_distance_factor` and `max_pooled_residuals`, assignment run on the identity and whole test equal construction +
    !| run_js_comp_test_adaptive with that option, and that value differs from the run without it, so the test fails if
    !| the option is dropped or goes to another argument. (`min_study_neighbors` changes no value, it only rejects
    !| neighborhoods; test_adaptive_on_assignment_rejects_invalid_input covers it.)
    subroutine test_adaptive_forwards_construction_options()
        character(*), parameter :: name = "test_adaptive_forwards_construction_options: "
        integer(int32), parameter :: N_SAMPLES = 15_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: expected(3), baseline(3), result(3), observed(3), p(3), null(1, 3)
        type(option_set_t) :: options, no_options
        character(len=24) :: labels(3)
        integer(int32) :: ierr, i_case, used_ks, used_kst, used_km, n_points, n_failed

        labels = [character(len=24) :: "tau", "mad_distance_factor", "max_pooled_residuals"]
        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 41_int32, &
                             expression)
        call reference_adaptive(expression, K_START, K_STEP, K_MAX, baseline, n_points, ierr, options=no_options)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"baseline ierr")
        do i_case = 1, 3
            call reset_options(options)
            select case (i_case)
            case (1)
                options%tau = 0.001_real64
            case (2)
                options%mad_distance_factor = 3.0_real64
            case (3)
                options%max_pooled_residuals = 60_int32
            end select
            call reference_adaptive(expression, K_START, K_STEP, K_MAX, expected, n_points, ierr, options=options)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(labels(i_case))//": reference ierr")
            call assert_true(any(expected /= baseline), name//trim(labels(i_case))//": the option changes the value")
            call run_js_comp_test_adaptive_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, N_SAMPLES, &
                                                                identity(N_SAMPLES), K_START, K_STEP, K_MAX, result, &
                                                                used_ks, used_kst, used_km, n_points, &
                                                                tau=options%tau, mad_distance_factor=options%mad_distance_factor, &
                                                                max_pooled_residuals=options%max_pooled_residuals, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(labels(i_case))//": assignment run ierr")
            call assert_true(equal_and_not_nan(result, expected), name//trim(labels(i_case))//": assignment run")
            call run_js_comp_test_adaptive_permutation(3_int32, 300_int32, 6_int32, expression, K_START, K_STEP, K_MAX, &
                                                       1_int32, observed, p, null, n_failed, tau=options%tau, &
                                                       mad_distance_factor=options%mad_distance_factor, &
                                                       max_pooled_residuals=options%max_pooled_residuals, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//trim(labels(i_case))//": whole test ierr")
            call assert_true(equal_and_not_nan(observed, expected), name//trim(labels(i_case))//": whole test")
        end do

    end subroutine test_adaptive_forwards_construction_options

    !> Test helper: leave every component of `options` unallocated (all defaults).
    subroutine reset_options(options)
        type(option_set_t), intent(inout) :: options
            !! The option set to clear

        if (allocated(options%tau)) deallocate (options%tau)
        if (allocated(options%mad_distance_factor)) deallocate (options%mad_distance_factor)
        if (allocated(options%max_pooled_residuals)) deallocate (options%max_pooled_residuals)
        if (allocated(options%min_study_neighbors)) deallocate (options%min_study_neighbors)
    end subroutine reset_options

end module mod_test_data_integration_js_comp_test_permutation_adaptive
