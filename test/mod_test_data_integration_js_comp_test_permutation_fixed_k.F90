!> Unit test suite for the k-NN study-label Permutation-Test of the JS-Comp-Test
!| (tox_data_integration_js_comp_test_permutation_fixed_k_impl): one assignment of samples to studies, the whole test
!| with its observed value, null distribution and p-values, and the rerun of the fixed-k parameter search.
module mod_test_data_integration_js_comp_test_permutation_fixed_k
    use asserts
    use, intrinsic :: iso_fortran_env, only: real64, int32
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    use tox_errors
    use test_suite, only: test_case
    use f42_random_gsl, only: random_shuffle
    use f42_sort_impl, only: init_perm
    use tox_data_integration_js_comp_test_fixed_k_search, only: generate_js_comp_test_candidates
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: MODE_NEIGHBORHOOD_SIZE_LARGE, &
                                                                     MODE_NEIGHBORHOOD_SIZE_SMALL
    use tox_data_integration_js_comp_test_impl, only: MAX_N_BINS, MAX_CANDIDATE_PAIRS
    use tox_data_integration_jsd_impl, only: MODE_DIVERGENCE_SYMMETRIC
    use mod_test_data_integration_js_comp_test, only: assert_equal_array_real_nan_exact
    use tox_data_integration_js_comp_test_permutation_fixed_k, only: run_js_comp_test_on_sample_assignment, &
                                                                     run_js_comp_test_permutation
    use tox_data_integration_js_comp_test_permutation_bounds, only: calc_js_comp_test_permutation_bounds
    use tox_data_integration_js_comp_test_permutation_fixed_k_impl, only: run_js_comp_test_on_sample_assignment_impl, &
                                                                          run_js_comp_test_permutation_impl

    use mod_test_data_integration_js_comp_test_permutation_common, only: make_expression, reference_means_residuals, &
                                                                         reference_knn, reference_fixed_k_search, &
                                                                         redraw_assignments, identity, &
                                                                         equal_and_not_nan, make_knn_grid_fixture, &
                                                                         N_POINTS_KNN, N_NULL_P_VALUES, &
                                                                         KNN_GRID_GENES, N_NEIGHBORS_KNN, &
                                                                         MIN_NULL_ABOVE_MINIMUM, KNN_GRID_SEED

    implicit none

contains

    !> Get array of all available tests.
    function get_all_tests_js_comp_test_permutation_fixed_k() result(all_tests)
        type(test_case), allocatable :: all_tests(:)
        allocate (all_tests(12))
        all_tests(1) = test_case("test_knn_identity_matches_plain_pipeline", test_knn_identity_matches_plain_pipeline)
        all_tests(2) = test_case("test_knn_swapping_two_studies_swaps_their_values", &
                                 test_knn_swapping_two_studies_swaps_their_values)
        all_tests(3) = test_case("test_knn_within_study_permutation_keeps_values", &
                                 test_knn_within_study_permutation_keeps_values)
        all_tests(4) = test_case("test_knn_on_assignment_rejects_invalid_input", &
                                 test_knn_on_assignment_rejects_invalid_input)
        all_tests(5) = test_case("test_knn_label_test_reproducible_rows_match_assignment_run", &
                                 test_knn_label_test_reproducible_rows_match_assignment_run)
        all_tests(6) = test_case("test_knn_label_test_null_and_batch", test_knn_label_test_null_and_batch)
        all_tests(7) = test_case("test_knn_rerun_parameter_search", test_knn_rerun_parameter_search)
        all_tests(8) = test_case("test_knn_rerun_forwards_neighborhood_size_mode", &
                                 test_knn_rerun_forwards_neighborhood_size_mode)
        all_tests(9) = test_case("test_knn_rerun_forwards_grid_options", test_knn_rerun_forwards_grid_options)
        all_tests(10) = test_case("test_knn_label_test_zero_permutations_direct_call", &
                                 test_knn_label_test_zero_permutations_direct_call)
        all_tests(11) = test_case("test_knn_label_test_rejects_invalid_input", &
                                 test_knn_label_test_rejects_invalid_input)
        all_tests(12) = test_case("test_knn_rerun_two_genes_low_grid_floor", test_knn_rerun_two_genes_low_grid_floor)
    end function get_all_tests_js_comp_test_permutation_fixed_k

    !> Test helper: assignment run k-NN through the plain wrapper, fixed parameters (`divergence_mode`
    !| forwarded, absent: the default).
    subroutine knn_on_assignment(expression, assignment, global_jsd, ierr, divergence_mode)
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
        integer(int32) :: used_np, used_nn

        call run_js_comp_test_on_sample_assignment(size(expression, 3, kind=int32), size(expression, 2, kind=int32), &
                                                   size(expression, 1, kind=int32), expression, &
                                                   size(assignment, kind=int32), assignment, N_POINTS_KNN, &
                                                   N_NEIGHBORS_KNN, global_jsd, used_np, used_nn, &
                                                   divergence_mode=divergence_mode, ierr=ierr)
    end subroutine knn_on_assignment

    !> Test helper: calls run_js_comp_test_on_sample_assignment_impl directly, every work array allocated here.
    subroutine direct_knn_on_assignment(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                        n_samples, sample_assignment, n_points, n_neighbors, &
                                        max_n_points_candidate, max_n_neighbors_candidate, &
                                        n_bootstrapping_top_k_jsds, global_js_divergence, used_n_points, &
                                        used_n_neighbors, ierr, rerun_parameter_search, n_bootstraps, random_seed, &
                                        divergence_mode, grid_points_factor, grid_min_points)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_samples
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_points
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_neighbors
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(out) :: used_n_points
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(out) :: used_n_neighbors
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        real(real64), intent(in), optional :: grid_points_factor
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), intent(in), optional :: grid_min_points
            !! Forwarded to run_js_comp_test_on_sample_assignment_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        logical(c_bool), allocatable :: tmp_sample_mark(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perms(:, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_neighborhood_indices(:, :, :)
        integer(int32), allocatable :: tmp_neighborhood_range(:, :, :)
        real(real64), allocatable :: tmp_neighborhood_residuals_gathered(:, :, :)
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
        integer(int32), allocatable :: tmp_point_n_neighbors(:)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        integer(int32), allocatable :: tmp_search_trace_n_neighbors(:)
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
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
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
        allocate (tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_points_candidate))
        allocate (tmp_neighborhood_indices(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_range(2, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
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
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_point_n_neighbors(n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_neighbors(MAX_CANDIDATE_PAIRS))
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
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_on_sample_assignment_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=sample_assignment, n_points=n_points, &
            n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
            max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, global_js_divergence=global_js_divergence, &
            used_n_points=used_n_points, used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, &
            tmp_sample_mark=tmp_sample_mark, tmp_permuted_expression=tmp_permuted_expression, &
            tmp_gene_means=tmp_gene_means, tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
            tmp_point_pool=tmp_point_pool, tmp_point_draw=tmp_point_draw, &
            tmp_n_bins_per_point=tmp_n_bins_per_point, &
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
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
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
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            divergence_mode=divergence_mode, grid_points_factor=grid_points_factor, grid_min_points=grid_min_points, &
            ierr=ierr)
    end subroutine direct_knn_on_assignment

    !> Test helper: calls run_js_comp_test_permutation_impl directly, every work array allocated here.
    subroutine direct_knn_label_test(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                     n_points, n_neighbors, max_n_points_candidate, max_n_neighbors_candidate, &
                                     n_bootstrapping_top_k_jsds, n_permutations, global_js_divergence_observed, &
                                     p_values, jsd_null, n_failed_permutations, ierr, rerun_parameter_search, &
                                     n_bootstraps, random_seed, &
                                     divergence_mode, grid_points_factor, grid_min_points)
        integer(int32), intent(in) :: n_studies
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Forwarded to run_js_comp_test_permutation_impl
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: n_points
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: n_neighbors
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: max_n_points_candidate
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in) :: n_permutations
            !! Forwarded to run_js_comp_test_permutation_impl
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Forwarded to run_js_comp_test_permutation_impl
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Forwarded to run_js_comp_test_permutation_impl
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(out) :: n_failed_permutations
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(out) :: ierr
            !! Error code of the implementation
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in), optional :: n_bootstraps
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in), optional :: random_seed
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in), optional :: divergence_mode
            !! Forwarded to run_js_comp_test_permutation_impl
        real(real64), intent(in), optional :: grid_points_factor
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), intent(in), optional :: grid_min_points
            !! Forwarded to run_js_comp_test_permutation_impl
        integer(int32), allocatable :: tmp_sample_slots(:)
        integer(int32), allocatable :: tmp_sample_assignment(:)
        real(real64), allocatable :: tmp_permuted_global_js_divergence(:)
        real(real64), allocatable :: tmp_permuted_expression(:, :, :)
        real(real64), allocatable :: tmp_gene_means(:, :)
        real(real64), allocatable :: tmp_residuals(:, :, :)
        integer(int32), allocatable :: tmp_gene_means_perms(:, :)
        integer(int32), allocatable :: tmp_gene_means_perm_all(:)
        real(real64), allocatable :: tmp_x_star(:)
        integer(int32), allocatable :: tmp_neighborhood_indices(:, :, :)
        integer(int32), allocatable :: tmp_neighborhood_range(:, :, :)
        real(real64), allocatable :: tmp_neighborhood_residuals_gathered(:, :, :)
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
        integer(int32), allocatable :: tmp_point_n_neighbors(:)
        integer(int32), allocatable :: tmp_search_n_bins_per_point(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_low(:)
        real(real64), allocatable :: tmp_search_shared_residual_range_high(:)
        real(real64), allocatable :: tmp_search_confidence_interval(:, :)
        integer(int32), allocatable :: tmp_search_trace_n_points(:)
        integer(int32), allocatable :: tmp_search_trace_n_neighbors(:)
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
        real(real64), allocatable :: tmp_candidate_confidence_interval(:, :)
        real(real64), allocatable :: tmp_bootstrapping_top_k_jsds(:, :, :)
        real(real64), allocatable :: tmp_prev_global_js_divergence(:)
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
        allocate (tmp_gene_means_perms(max_n_genes_all_studies, n_studies))
        allocate (tmp_gene_means_perm_all(max_n_genes_all_studies*n_studies))
        allocate (tmp_x_star(max_n_points_candidate))
        allocate (tmp_neighborhood_indices(max_n_neighbors_candidate, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_range(2, max_n_points_candidate, n_studies))
        allocate (tmp_neighborhood_residuals_gathered(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate))
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
        allocate (tmp_pooled_residuals(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_pooled_residuals_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_residual_to_study_perm(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies))
        allocate (tmp_bin_counts_search(MAX_N_BINS))
        allocate (tmp_study_block_ends(n_studies))
        allocate (tmp_study_block_sizes(n_studies))
        allocate (tmp_point_n_neighbors(n_studies))
        allocate (tmp_search_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_search_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_search_confidence_interval(2, n_studies))
        allocate (tmp_search_trace_n_points(MAX_CANDIDATE_PAIRS))
        allocate (tmp_search_trace_n_neighbors(MAX_CANDIDATE_PAIRS))
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
        allocate (tmp_candidate_confidence_interval(2, n_studies))
        allocate (tmp_bootstrapping_top_k_jsds(n_bootstrapping_top_k_jsds, 2, n_studies))
        allocate (tmp_prev_global_js_divergence(n_studies))
        allocate (tmp_best_uncertainty_confidence_interval(2, n_studies))
        allocate (tmp_best_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_uncertainty_n_bins_per_point(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_shared_residual_range_high(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_low(max_n_points_candidate))
        allocate (tmp_best_uncertainty_shared_residual_range_high(max_n_points_candidate))
        call run_js_comp_test_permutation_impl(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_points=n_points, n_neighbors=n_neighbors, &
            max_n_points_candidate=max_n_points_candidate, max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, n_permutations=n_permutations, &
            global_js_divergence_observed=global_js_divergence_observed, p_values=p_values, jsd_null=jsd_null, &
            n_failed_permutations=n_failed_permutations, tmp_sample_slots=tmp_sample_slots, &
            tmp_sample_assignment=tmp_sample_assignment, &
            tmp_permuted_global_js_divergence=tmp_permuted_global_js_divergence, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
            tmp_point_pool=tmp_point_pool, tmp_point_draw=tmp_point_draw, &
            tmp_n_bins_per_point=tmp_n_bins_per_point, &
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
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
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
            tmp_candidate_confidence_interval=tmp_candidate_confidence_interval, &
            tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
            tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, &
            tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
            tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
            tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
            tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
            tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
            tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
            tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
            rerun_parameter_search=rerun_parameter_search, n_bootstraps=n_bootstraps, random_seed=random_seed, &
            divergence_mode=divergence_mode, grid_points_factor=grid_points_factor, grid_min_points=grid_min_points, &
            ierr=ierr)
    end subroutine direct_knn_label_test

    !> The identity assignment reproduces, bit for bit, the plain k-NN pipeline on the observed
    !| data, with unequal sample counts (NaN padding) and gene 7 missing in study 2.
    subroutine test_knn_identity_matches_plain_pipeline()
        character(*), parameter :: name = "test_knn_identity_matches_plain_pipeline: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: reference(3), result(3)
        integer(int32) :: ierr, used_np, used_nn, n_samples

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 11_int32, &
                             expression)
        expression(:, 7, 2) = ieee_value(1.0_real64, ieee_quiet_nan)
        n_samples = 15_int32

        call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")
        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 6_int32, expression, n_samples, &
                                                   identity(n_samples), N_POINTS_KNN, N_NEIGHBORS_KNN, result, used_np, &
                                                   used_nn, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call assert_equal_array_real_nan_exact(result, reference, size(reference, kind=int32), name//"identity equals the plain pipeline exactly")
        call assert_no_nan_real(result, 3_int32, name//"no NaN")
        call assert_equal_int(used_np, N_POINTS_KNN, name//"used_n_points")
        call assert_equal_int(used_nn, N_NEIGHBORS_KNN, name//"used_n_neighbors")
    end subroutine test_knn_identity_matches_plain_pipeline

    !> Two studies of 4 samples each: the assignment that swaps the studies' samples swaps their
    !| values exactly (J1' = J2, J2' = J1). With the symmetric divergence, passed explicitly: the
    !| default pairwise divergence gives two studies one shared value, so there is nothing to swap.
    subroutine test_knn_swapping_two_studies_swaps_their_values()
        character(*), parameter :: name = "test_knn_swapping_two_studies_swaps_their_values: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), swapped(2)
        integer(int32) :: ierr

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 21_int32, expression, dyadic=.true._c_bool)
        call knn_on_assignment(expression, identity(8_int32), observed, ierr, divergence_mode=MODE_DIVERGENCE_SYMMETRIC)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"observed ierr")
        call knn_on_assignment(expression, [5_int32, 6_int32, 7_int32, 8_int32, 1_int32, 2_int32, 3_int32, 4_int32], &
                               swapped, ierr, divergence_mode=MODE_DIVERGENCE_SYMMETRIC)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"swapped ierr")
        call assert_true(observed(1) /= observed(2), name//"the two studies differ")
        call assert_true(swapped(1) == observed(2), name//"J1' == J2")
        call assert_true(swapped(2) == observed(1), name//"J2' == J1")
    end subroutine test_knn_swapping_two_studies_swaps_their_values

    !> An assignment that only reorders the samples inside each study leaves every value unchanged
    !| (dyadic data, so the gene means are exact in any order).
    subroutine test_knn_within_study_permutation_keeps_values()
        character(*), parameter :: name = "test_knn_within_study_permutation_keeps_values: "
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(2), permuted(2)
        integer(int32) :: ierr

        call make_expression(200_int32, [4_int32, 4_int32], [0.5_real64, 0.8_real64], 22_int32, expression, dyadic=.true._c_bool)
        call knn_on_assignment(expression, identity(8_int32), observed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"observed ierr")
        call knn_on_assignment(expression, [3_int32, 1_int32, 4_int32, 2_int32, 8_int32, 7_int32, 6_int32, 5_int32], &
                               permuted, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"permuted ierr")
        call assert_equal_array_real_nan_exact(permuted, observed, size(observed, kind=int32), name//"values unchanged")
    end subroutine test_knn_within_study_permutation_keeps_values

    !> Invalid inputs of the k-NN assignment run: through the plain wrapper (runtime errors carry no
    !| position there) and through the implementation (each names its argument).
    subroutine test_knn_on_assignment_rejects_invalid_input()
        character(*), parameter :: name = "test_knn_on_assignment_rejects_invalid_input: "
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :)
        real(real64) :: result(3)
        integer(int32) :: ierr, used_np, used_nn, max_p, max_k, top_k

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 31_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 2) = ieee_value(1.0_real64, ieee_quiet_nan)

        ! Plain wrapper.
        call knn_on_assignment(expression, [1_int32, 1_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                            9_int32, 10_int32, 11_int32, 12_int32], result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"repeated value rejected", arg_pos=0_int32)
        call assert_true(all(ieee_is_nan(result)), name//"result NaN on error")
        call knn_on_assignment(expression, identity(11_int32), result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_samples below the sample count rejected", arg_pos=0_int32)
        call knn_on_assignment(expression, identity(13_int32), result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"n_samples above the sample count rejected", arg_pos=0_int32)
        call knn_on_assignment(empty_study, identity(8_int32), result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"study without samples rejected", arg_pos=0_int32)
        call knn_on_assignment(expression, [0_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, &
                                            9_int32, 10_int32, 11_int32, 12_int32], result, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"value 0 rejected by the wrapper", arg_pos=6_int32)

        ! Implementation: positions, and the rerun capacity checks.
        call calc_js_comp_test_permutation_bounds(300_int32, N_POINTS_KNN, N_NEIGHBORS_KNN, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, &
                                      [2_int32, 2_int32, 3_int32, 4_int32, 5_int32, 6_int32, 7_int32, 8_int32, 9_int32, &
                                       10_int32, 11_int32, 12_int32], N_POINTS_KNN, N_NEIGHBORS_KNN, max_p, max_k, top_k, &
                                      result, used_np, used_nn, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: repeated value", arg_pos=6_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 11_int32, identity(11_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, max_p, max_k, top_k, result, used_np, used_nn, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: n_samples", arg_pos=5_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, empty_study, 8_int32, identity(8_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, max_p, max_k, top_k, result, used_np, used_nn, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: empty study", arg_pos=4_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, max_p - 1_int32, max_k, top_k, result, used_np, used_nn, ierr, &
                                      rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: point capacity below the search's", arg_pos=9_int32)
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), 4_int32, &
                                      4_int32, max_p, 3_int32, top_k, result, used_np, used_nn, ierr, &
                                      rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"impl: neighbor capacity below the search's", arg_pos=10_int32)
        ! The same small capacities without the rerun are fine.
        call direct_knn_on_assignment(3_int32, 300_int32, 4_int32, expression, 12_int32, identity(12_int32), N_POINTS_KNN, &
                                      N_NEIGHBORS_KNN, N_POINTS_KNN, N_NEIGHBORS_KNN, top_k, result, used_np, used_nn, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"impl: exact capacities without rerun")
    end subroutine test_knn_on_assignment_rejects_invalid_input

    !> k-NN: a fixed seed reproduces the test (identical null rows), another seed gives other null
    !| rows; the p-values lie in `[1/(n+1), 1]`; the observed value is assignment run on the identity; every
    !| row of `jsd_null` is assignment run on the assignment the test drew (replayed here through
    !| `init_perm` + `random_shuffle` on a generator with the same seed, each draw checked to be a
    !| permutation).
    subroutine test_knn_label_test_reproducible_rows_match_assignment_run()
        character(*), parameter :: name = "test_knn_label_test_reproducible_rows_match_assignment_run: "
        integer(int32), parameter :: n = 6_int32, seed = 7_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), observed2(3), p2(3), null2(n, 3), row(3)
        integer(int32) :: ierr, n_failed, n_failed2, assignments(15, n), i

        call make_expression(300_int32, [6_int32, 5_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 41_int32, &
                             expression)
        call run_js_comp_test_permutation(3_int32, 300_int32, 6_int32, expression, N_POINTS_KNN, &
                                          N_NEIGHBORS_KNN, n, observed, p, null, n_failed, random_seed=seed, &
                                          ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"ierr")
        call run_js_comp_test_permutation(3_int32, 300_int32, 6_int32, expression, N_POINTS_KNN, &
                                          N_NEIGHBORS_KNN, n, observed2, p2, null2, n_failed2, random_seed=seed, &
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
        ! Another seed draws other assignments, so other null rows (the observed value stays).
        call run_js_comp_test_permutation(3_int32, 300_int32, 6_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, n, &
                                          observed2, p2, null2, n_failed2, random_seed=seed + 1_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"other seed ierr")
        call assert_equal_array_real_nan_exact(observed, observed2, 3_int32, name//"other seed: same observed values")
        call assert_true(any(null /= null2), name//"other seed: other null rows")

        call knn_on_assignment(expression, identity(15_int32), row, ierr)
        call assert_equal_array_real_nan_exact(observed, row, size(row, kind=int32), name//"observed is assignment run on the identity")
        call redraw_assignments(seed, 15_int32, n, assignments)
        do i = 1, n
            call knn_on_assignment(expression, assignments(:, i), row, ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"assignment run ierr")
            call assert_equal_array_real_nan_exact(null(i, :), row, size(row, kind=int32), &
                                                   name//"row equals assignment run on the drawn assignment")
        end do
        call assert_true(any(assignments(:, 1) /= identity(15_int32)), name//"the draw is not the identity")
    end subroutine test_knn_label_test_reproducible_rows_match_assignment_run

    !> k-NN: on null data (every study drawn from the same distribution) the p-values are mostly
    !| well above the minimum `1/(n+1)`; a study with four times the noise gets exactly `1/(n+1)`.
    subroutine test_knn_label_test_null_and_batch()
        character(*), parameter :: name = "test_knn_label_test_null_and_batch: "
        integer(int32), parameter :: n = 19_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: observed(3), p(3), null(n, 3), p_min
        integer(int32) :: ierr, n_failed, i_round, n_above

        p_min = 1.0_real64/real(n + 1, real64)
        n_above = 0_int32
        do i_round = 1, 4
            call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 0.5_real64], &
                                 100_int32 + i_round, expression)
            call run_js_comp_test_permutation(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, &
                                              N_NEIGHBORS_KNN, n, observed, p, null, n_failed, &
                                              random_seed=i_round, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null ierr")
            ! Exact comparison: p and p_min are both (count + 1)/(n + 1) quotients computed alike.
            n_above = n_above + count(p > p_min, kind=int32)
        end do
        call assert_true(n_above >= MIN_NULL_ABOVE_MINIMUM, name//"null p-values mostly above the minimum (got "// &
                         n_above//" of "//N_NULL_P_VALUES//")")

        call make_expression(300_int32, [5_int32, 5_int32, 5_int32], [0.5_real64, 0.5_real64, 2.0_real64], 109_int32, &
                             expression)
        call run_js_comp_test_permutation(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, &
                                          N_NEIGHBORS_KNN, n, observed, p, null, n_failed, random_seed=5_int32, &
                                          ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"batch ierr")
        call assert_true(p(3) == p_min, name//"batch study at the minimum p-value")
    end subroutine test_knn_label_test_null_and_batch

    !> k-NN rerun option: on, the result equals the search on the identity data followed by the run
    !| with the parameters it selected; off, the fixed-parameter run (test 1). The whole test with
    !| the rerun option reports the same observed value.
    subroutine test_knn_rerun_parameter_search()
        character(*), parameter :: name = "test_knn_rerun_parameter_search: "
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: reference(3), result(3), fixed(3), observed(3), p(3), null(1, 3), null4(4, 3), row(3)
        integer(int32) :: ierr, used_np, used_nn, search_np, search_nn, n_failed, i, assignments(13, 4)

        call make_expression(300_int32, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], 51_int32, &
                             expression)
        call reference_means_residuals(expression, gene_means, residuals)
        call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np, search_nn)
        call reference_knn(expression, search_np, search_nn, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")

        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, identity(13_int32), &
                                                   N_POINTS_KNN, N_NEIGHBORS_KNN, result, used_np, used_nn, &
                                                   rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"rerun ierr")
        call assert_equal_int(used_np, search_np, name//"used_n_points is the search's")
        call assert_equal_int(used_nn, search_nn, name//"used_n_neighbors is the search's")
        call assert_equal_array_real_nan_exact(result, reference, size(reference, kind=int32), name//"rerun equals search + run")

        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, identity(13_int32), &
                                                   N_POINTS_KNN, N_NEIGHBORS_KNN, fixed, used_np, used_nn, &
                                                   rerun_parameter_search=.false._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off ierr")
        call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"off reference ierr")
        call assert_true(equal_and_not_nan(fixed, reference), name//"off equals the fixed-parameter run")
        call assert_equal_int(used_np, N_POINTS_KNN, name//"off: used_n_points")

        call run_js_comp_test_permutation(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, &
                                          N_NEIGHBORS_KNN, 1_int32, observed, p, null, n_failed, &
                                          rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"label test with rerun ierr")
        call assert_equal_array_real_nan_exact(observed, result, size(result, kind=int32), name//"label test with rerun: observed value")

        ! The null rows of the shuffled data sets, each of which reruns the search: row i equals assignment run with the rerun
        ! on the i-th drawn assignment (seed 5, the same bootstrap seed).
        call run_js_comp_test_permutation(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                          4_int32, observed, p, null4, n_failed, random_seed=5_int32, &
                                          rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null rows with rerun: ierr")
        call redraw_assignments(5_int32, 13_int32, 4_int32, assignments)
        do i = 1, 4
            call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, assignments(:, i), &
                                                       N_POINTS_KNN, N_NEIGHBORS_KNN, row, used_np, used_nn, &
                                                       random_seed=5_int32, rerun_parameter_search=.true._c_bool, &
                                                       n_bootstraps=20_int32, ierr=ierr)
            call assert_equal_int(get_err_code(ierr), ERR_OK, name//"null row assignment run ierr")
            call assert_equal_array_real_nan_exact(null4(i, :), row, size(row, kind=int32), &
                                                   name//"null row equals assignment run with the rerun")
        end do
        call assert_true(any(null4(1, :) /= observed), name//"a shuffled data set gives other values than the observed one")
    end subroutine test_knn_rerun_parameter_search

    !> With the rerun option, the k-NN routines forward `neighborhood_size_mode` to the parameter
    !| search. Fixture with a multi-candidate default grid (KNN_GRID_GENES genes, data seed
    !| KNN_GRID_SEED, see there): the large and the small grid have at least three candidates each
    !| and share no `n_neighbors`, so the selected value tells which grid the search walked. For the
    !| default, explicit large and small mode, assignment run with rerun reports exactly the
    !| `(n_points, n_neighbors)` the parameter search selects directly with that mode, and its
    !| value equals the run at that pair; the default equals explicit large, and small selects a
    !| different `n_neighbors` than large (so the test fails if the mode is not forwarded). The
    !| label test with rerun and the small mode reports the same observed value.
    subroutine test_knn_rerun_forwards_neighborhood_size_mode()
        character(*), parameter :: name = "test_knn_rerun_forwards_neighborhood_size_mode: "
        integer(int32), parameter :: MIN_GRID_CANDIDATES = 3_int32
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: default_result(3), large_result(3), small_result(3), reference(3), observed(3), p(3), &
                        null(1, 3)
        integer(int32) :: ierr, used_np, used_nn, n_large, n_small, n_failed, search_np, search_nn, large_nn, n_samples
        integer(int32) :: large_grid(2, MAX_CANDIDATE_PAIRS), small_grid(2, MAX_CANDIDATE_PAIRS)

        call make_knn_grid_fixture(expression)
        n_samples = count(.not. ieee_is_nan(expression(:, 1, :)), kind=int32)
        call reference_means_residuals(expression, gene_means, residuals)
        call generate_js_comp_test_candidates(KNN_GRID_GENES, large_grid, n_large, &
                                              neighborhood_size_mode=MODE_NEIGHBORHOOD_SIZE_LARGE, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"large grid ierr")
        call generate_js_comp_test_candidates(KNN_GRID_GENES, small_grid, n_small, &
                                              neighborhood_size_mode=MODE_NEIGHBORHOOD_SIZE_SMALL, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"small grid ierr")
        call assert_true(n_large >= MIN_GRID_CANDIDATES .and. n_small >= MIN_GRID_CANDIDATES, &
                         name//"both grids have at least three candidates")
        call assert_true(.not. any(spread(large_grid(2, 1:n_large), 1, n_small) == &
                                   spread(small_grid(2, 1:n_small), 2, n_large)), &
                         name//"the two grids share no n_neighbors")

        call run_js_comp_test_on_sample_assignment(3_int32, size(expression, 2, kind=int32), &
                                                   size(expression, 1, kind=int32), expression, n_samples, &
                                                   identity(n_samples), N_POINTS_KNN, N_NEIGHBORS_KNN, default_result, &
                                                   used_np, used_nn, rerun_parameter_search=.true._c_bool, &
                                                   n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"default ierr")
        call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np, search_nn)
        call assert_equal_array_int([used_np, used_nn], [search_np, search_nn], 2_int32, &
                                    name//"default: the search's own selection")
        call assert_true(any(large_grid(2, 1:n_large) == used_nn), name//"default: n_neighbors from the large grid")
        large_nn = used_nn

        call run_js_comp_test_on_sample_assignment(3_int32, size(expression, 2, kind=int32), &
                                                   size(expression, 1, kind=int32), expression, n_samples, &
                                                   identity(n_samples), N_POINTS_KNN, N_NEIGHBORS_KNN, large_result, &
                                                   used_np, used_nn, rerun_parameter_search=.true._c_bool, &
                                                   n_bootstraps=20_int32, &
                                                   neighborhood_size_mode=MODE_NEIGHBORHOOD_SIZE_LARGE, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"large ierr")
        call assert_true(equal_and_not_nan(large_result, default_result), name//"explicit large equals the default")

        call run_js_comp_test_on_sample_assignment(3_int32, size(expression, 2, kind=int32), &
                                                   size(expression, 1, kind=int32), expression, n_samples, &
                                                   identity(n_samples), N_POINTS_KNN, N_NEIGHBORS_KNN, small_result, &
                                                   used_np, used_nn, rerun_parameter_search=.true._c_bool, &
                                                   n_bootstraps=20_int32, &
                                                   neighborhood_size_mode=MODE_NEIGHBORHOOD_SIZE_SMALL, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"small ierr")
        call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np, search_nn, &
                                      neighborhood_size_mode=MODE_NEIGHBORHOOD_SIZE_SMALL)
        call assert_equal_array_int([used_np, used_nn], [search_np, search_nn], 2_int32, &
                                    name//"small: the search's own selection with the small mode")
        call assert_true(any(small_grid(2, 1:n_small) == used_nn), name//"small: n_neighbors from the small grid")
        call assert_true(used_nn /= large_nn, name//"small selects another n_neighbors than large")
        call reference_knn(expression, used_np, used_nn, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"small reference ierr")
        call assert_true(equal_and_not_nan(small_result, reference), name//"small equals the run at the selected pair")

        call run_js_comp_test_permutation(3_int32, size(expression, 2, kind=int32), size(expression, 1, kind=int32), &
                                          expression, N_POINTS_KNN, N_NEIGHBORS_KNN, 1_int32, observed, p, null, n_failed, &
                                          rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                          neighborhood_size_mode=MODE_NEIGHBORHOOD_SIZE_SMALL, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"label test ierr")
        call assert_true(equal_and_not_nan(observed, small_result), name//"label test forwards the mode")
    end subroutine test_knn_rerun_forwards_neighborhood_size_mode

    !> With the rerun option, the two k-NN routines forward `grid_points_factor` and
    !| `grid_min_points` to the parameter search (and size its arrays for that grid). G = 300: the
    !| default grid is the single candidate (200, 6) (`ceiling(4*sqrt(300)) = 70` clamps up to the
    !| floor 200), so the default rerun uses 200 points; with `grid_points_factor = 17.32` and
    !| `grid_min_points = 240` the grid is (300, 4), (240, 5) (`ceiling(17.32*sqrt(300)) = 300`,
    !| `0.8*300 = 240` the floor itself), so the rerun uses 300 or 240 points -- another value than
    !| the default's, and exactly the pair the parameter search selects directly with those
    !| options. Assignment run then equals the run at that pair, and the label test with rerun reports
    !| the same observed value.
    subroutine test_knn_rerun_forwards_grid_options()
        character(*), parameter :: name = "test_knn_rerun_forwards_grid_options: "
        real(real64), parameter :: LARGER_POINTS_FACTOR = 17.32_real64
        integer(int32), parameter :: LARGER_MIN_POINTS = 240_int32, DEFAULT_GRID_POINTS = 200_int32
        real(real64), allocatable :: expression(:, :, :), gene_means(:, :), residuals(:, :, :)
        real(real64) :: default_result(3), result(3), reference(3), observed(3), p(3), null(1, 3)
        integer(int32) :: ierr, used_np, used_nn, search_np, search_nn, n_failed

        call make_expression(300_int32, [5_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.6_real64], 51_int32, &
                             expression)
        call reference_means_residuals(expression, gene_means, residuals)

        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, identity(13_int32), &
                                                   N_POINTS_KNN, N_NEIGHBORS_KNN, default_result, used_np, used_nn, &
                                                   rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"default ierr")
        call assert_equal_int(used_np, DEFAULT_GRID_POINTS, name//"default grid: 200 points")

        call run_js_comp_test_on_sample_assignment(3_int32, 300_int32, 5_int32, expression, 13_int32, identity(13_int32), &
                                                   N_POINTS_KNN, N_NEIGHBORS_KNN, result, used_np, used_nn, &
                                                   rerun_parameter_search=.true._c_bool, n_bootstraps=20_int32, &
                                                   grid_points_factor=LARGER_POINTS_FACTOR, &
                                                   grid_min_points=LARGER_MIN_POINTS, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"grid options ierr")
        call reference_fixed_k_search(gene_means, residuals, 20_int32, search_np, search_nn, &
                                      grid_points_factor=LARGER_POINTS_FACTOR, grid_min_points=LARGER_MIN_POINTS)
        call assert_equal_array_int([used_np, used_nn], [search_np, search_nn], 2_int32, &
                                    name//"grid options: the search's own selection with those options")
        call assert_true(used_np /= DEFAULT_GRID_POINTS, name//"grid options: another n_points than the default grid's")
        call reference_knn(expression, used_np, used_nn, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"reference ierr")
        call assert_true(equal_and_not_nan(result, reference), name//"grid options: equals the run at the selected pair")

        call run_js_comp_test_permutation(3_int32, 300_int32, 5_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, &
                                          1_int32, observed, p, null, n_failed, rerun_parameter_search=.true._c_bool, &
                                          n_bootstraps=20_int32, grid_points_factor=LARGER_POINTS_FACTOR, &
                                          grid_min_points=LARGER_MIN_POINTS, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"label test ierr")
        call assert_true(equal_and_not_nan(observed, result), name//"label test forwards the grid options")
    end subroutine test_knn_rerun_forwards_grid_options

    !> A direct implementation call with `n_permutations = 0` (k-NN): the observed value is
    !| computed, the p-values are NaN, nothing failed.
    subroutine test_knn_label_test_zero_permutations_direct_call()
        character(*), parameter :: name = "test_knn_label_test_zero_permutations_direct_call: "
        real(real64), allocatable :: expression(:, :, :), null(:, :)
        real(real64) :: observed(3), p(3), reference(3)
        integer(int32) :: ierr, n_failed, max_p, max_k, top_k

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 61_int32, &
                             expression)
        allocate (null(0, 3))
        call calc_js_comp_test_permutation_bounds(300_int32, N_POINTS_KNN, N_NEIGHBORS_KNN, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, max_p, max_k, &
                                   top_k, 0_int32, observed, p, null, n_failed, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"knn ierr")
        call assert_true(all(ieee_is_nan(p)), name//"knn p-values NaN")
        call assert_equal_int(n_failed, 0_int32, name//"knn n_failed")
        call reference_knn(expression, N_POINTS_KNN, N_NEIGHBORS_KNN, reference, ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"knn reference ierr")
        call assert_equal_array_real_nan_exact(observed, reference, size(reference, kind=int32), name//"knn observed")
    end subroutine test_knn_label_test_zero_permutations_direct_call

    !> Whole test of the k-NN Permutation-Test rejects what does not depend on the assignment with its own
    !| positions (a study without samples, a rerun capacity below the search's).
    subroutine test_knn_label_test_rejects_invalid_input()
        character(*), parameter :: name = "test_knn_label_test_rejects_invalid_input: "
        real(real64), allocatable :: expression(:, :, :), empty_study(:, :, :)
        real(real64) :: observed(3), p(3), null(2, 3)
        integer(int32) :: ierr, n_failed, max_p, max_k, top_k

        call make_expression(300_int32, [4_int32, 4_int32, 4_int32], [0.5_real64, 0.5_real64, 0.5_real64], 81_int32, &
                             expression)
        empty_study = expression
        empty_study(:, :, 1) = ieee_value(1.0_real64, ieee_quiet_nan)
        call calc_js_comp_test_permutation_bounds(300_int32, N_POINTS_KNN, N_NEIGHBORS_KNN, max_n_points_candidate=max_p, &
                                                  max_n_neighbors_candidate=max_k, n_bootstrapping_top_k_jsds=top_k, &
                                                  ierr=ierr)
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, empty_study, N_POINTS_KNN, N_NEIGHBORS_KNN, max_p, max_k, &
                                   top_k, 2_int32, observed, p, null, n_failed, ierr)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"knn: empty study", arg_pos=4_int32)
        call assert_true(all(ieee_is_nan(p)) .and. all(ieee_is_nan(null)), name//"knn: NaN outputs")
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, expression, N_POINTS_KNN, N_NEIGHBORS_KNN, max_p - 1_int32, &
                                   max_k, top_k, 2_int32, observed, p, null, n_failed, ierr, &
                                   rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"knn: point capacity", arg_pos=7_int32)
        call direct_knn_label_test(3_int32, 300_int32, 4_int32, expression, 4_int32, 4_int32, max_p, 3_int32, top_k, &
                                   2_int32, observed, p, null, n_failed, ierr, rerun_parameter_search=.true._c_bool)
        call assert_err(ierr, ERR_INVALID_INPUT, name//"knn: neighbor capacity", arg_pos=8_int32)
    end subroutine test_knn_label_test_rejects_invalid_input

    !> The capped fixed-k grid (test_fixed_k_grid_caps_n_neighbors_at_gene_count in the main suite)
    !| through the k-NN Permutation-Test with the rerun option: G = 2 genes, two studies of 4
    !| samples, `grid_min_points = 2`, whose grid is (6, 1), (5, 1), (4, 2), (3, 2), (2, 2) (the
    !| last capped from 3). Assignment run on the identity and whole test with 3 permutations rerun the
    !| search on every data set and return without error, the run using at most 2 neighbors.
    !| Before the cap the search built neighborhoods of 3 genes out of 2 (gene index 3, past the
    !| data); a `--diagnostics` build (`-fcheck=all`) is what catches such a read.
    subroutine test_knn_rerun_two_genes_low_grid_floor()
        character(*), parameter :: name = "test_knn_rerun_two_genes_low_grid_floor: "
        integer(int32), parameter :: n = 3_int32, LOW_MIN_POINTS = 2_int32
        real(real64), allocatable :: expression(:, :, :)
        real(real64) :: result(2), observed(2), p(2), null(n, 2)
        integer(int32) :: ierr, used_np, used_nn, n_failed

        call make_expression(2_int32, [4_int32, 4_int32], [0.5_real64, 0.5_real64], 61_int32, expression)
        call run_js_comp_test_on_sample_assignment(2_int32, 2_int32, 4_int32, expression, 8_int32, identity(8_int32), &
                                                   2_int32, 2_int32, result, used_np, used_nn, &
                                                   rerun_parameter_search=.true._c_bool, n_bootstraps=5_int32, &
                                                   grid_min_points=LOW_MIN_POINTS, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"assignment run ierr")
        call assert_true(used_nn >= 1_int32 .and. used_nn <= 2_int32, name//"assignment run: used_n_neighbors <= G")
        call assert_true(used_np >= 2_int32 .and. used_np <= 6_int32, name//"assignment run: used_n_points from the grid")

        call run_js_comp_test_permutation(2_int32, 2_int32, 4_int32, expression, 2_int32, 2_int32, n, observed, p, null, &
                                          n_failed, rerun_parameter_search=.true._c_bool, n_bootstraps=5_int32, &
                                          grid_min_points=LOW_MIN_POINTS, random_seed=3_int32, ierr=ierr)
        call assert_equal_int(get_err_code(ierr), ERR_OK, name//"whole test ierr")
        call assert_true(equal_and_not_nan(observed, result), name//"whole test: observed is assignment run on the identity")
    end subroutine test_knn_rerun_two_genes_low_grid_floor

end module mod_test_data_integration_js_comp_test_permutation_fixed_k
