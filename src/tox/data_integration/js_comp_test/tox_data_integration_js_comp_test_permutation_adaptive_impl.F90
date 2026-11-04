#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test, adaptive
!|
!| The adaptive counterpart of [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module)]], with the same input layout, assignment
!| convention and p-values: one assignment of samples to slots, one adaptive JS-Comp-Test
!| ([[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]), and the whole test
!| ([[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_permutation_impl(interface)]]). The shared helpers are in
!| [[tox_data_integration_js_comp_test_permutation_common_impl(module)]].
module tox_data_integration_js_comp_test_permutation_adaptive_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, ERR_OK, ERR_INVALID_INPUT
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_shuffle
    use f42_sort_impl, only: init_perm
    use tox_data_integration_preprocessing_impl, only: construct_adaptive_neighborhoods_impl
    use tox_data_integration_js_comp_test_adaptive_search_impl, only: run_js_comp_test_adaptive_parameter_search_impl
    use tox_data_integration_js_comp_test_run_impl, only: run_js_comp_test_adaptive_impl
    use tox_data_integration_js_comp_test_permutation_common_impl, only: check_permutation_input, &
                                                                         assign_samples_and_center, &
                                                                         label_permutation_p_values
    M_IMPLICIT_NONE
    private
    public :: run_js_comp_test_adaptive_on_sample_assignment_impl, run_js_comp_test_adaptive_permutation_impl

contains

    !> summary: Run the adaptive JS-Comp-Test on one assignment of samples to studies
    !| AUTHOR_LASZLO_LANG
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]],
    !| with the same input layout and assignment convention (the identity is the observed data).
    !|
    !| Steps: per study, the gene means and residuals of the relabelled data; the adaptive neighborhoods for the
    !| growth knobs `k_start`, `k_step`, `k_max` and the construction options
    !| ([[tox_data_integration_preprocessing_impl(module):construct_adaptive_neighborhoods_impl(interface)]]); then
    !| [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_adaptive_impl(interface)]] on its ranges, with `divergence_mode` as given. `n_points` is the number
    !| of neighborhoods the construction built.
    !|
    !| With `rerun_parameter_search = .true.` the relabelled data first run
    !| [[tox_data_integration_js_comp_test_adaptive_search_impl(module):run_js_comp_test_adaptive_parameter_search_impl(interface)]]
    !| (with its own capacity, so it rejects the same candidates as a direct call), and the construction and run use the knobs
    !| it returns (also its fallback when no candidate was admissible). `k_start`, `k_step` and `k_max` are then not used, except
    !| `k_start` to size the work arrays, and `used_k_start`, `used_k_step`, `used_k_max` report what was used. Every call runs the
    !| whole JS-Comp-Test, with the rerun option the whole parameter search first.
    !|
    !| Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
    !| number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. Also when the construction does not succeed: its status is not ok (too few means, a
    !| neighborhood with too few residuals or with fewer than `min_study_neighbors` entries of some study), it builds no neighborhood, or more
    !| neighborhoods than the work arrays hold. In the Fortran expert tier, with the rerun option, also when the capacity is below
    !| the search's own.
    !| Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
    !| Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
    !| values for it.
    subroutine run_js_comp_test_adaptive_on_sample_assignment_impl(n_studies, max_n_genes_all_studies, &
                                                                   max_n_reps_all_studies, expression, n_samples, &
                                                                   sample_assignment, k_start, k_step, k_max, &
                                                                   max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                                                   global_js_divergence, used_k_start, used_k_step, &
                                                                   used_k_max, n_points, tmp_sample_slots, tmp_sample_mark, &
                                                                   tmp_permuted_expression, tmp_gene_means, tmp_residuals, &
                                                                   tmp_gene_means_perm_all, tmp_x_star, &
                                                                   tmp_pooled_neighborhood_range, tmp_n_neighbors_per_point, &
                                                                   tmp_stop_reason, tmp_neighborhood_dispersion, &
                                                                   tmp_neighborhood_mad, tmp_point_neighborhood_indices, &
                                                                   tmp_neighbor_residuals, tmp_point_pool, tmp_point_draw, &
                                                                   tmp_n_bins_per_point, tmp_shared_residual_range_low, &
                                                                   tmp_shared_residual_range_high, tmp_occupancy_failed, &
                                                                   tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                                   tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                                                   tmp_sturges_bins, tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                                   tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                                   tmp_mean_pmf_included_n_reps, tmp_js_divergences, &
                                                                   tmp_weights, tmp_global_js_divergence, &
                                                                   tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                                   tmp_residual_to_study_perm, tmp_bin_counts_search, &
                                                                   tmp_study_block_ends, tmp_study_block_sizes, &
                                                                   tmp_search_n_bins_per_point, &
                                                                   tmp_search_shared_residual_range_low, &
                                                                   tmp_search_shared_residual_range_high, &
                                                                   tmp_search_confidence_interval, tmp_search_trace_k_start, &
                                                                   tmp_search_trace_k_step, tmp_search_trace_k_max, &
                                                                   tmp_search_trace_n_points, &
                                                                   tmp_search_trace_global_js_divergence, &
                                                                   tmp_search_trace_ci_lower, tmp_search_trace_ci_upper, &
                                                                   tmp_search_trace_ci_width, &
                                                                   tmp_search_trace_ci_width_relative, &
                                                                   tmp_search_trace_ci_overlap, tmp_search_trace_delta, &
                                                                   tmp_search_trace_delta_median, tmp_search_trace_delta_max, &
                                                                   tmp_search_trace_selected_n_bins, &
                                                                   tmp_search_trace_occupancy_failed, &
                                                                   tmp_search_trace_n_pooled_residuals, &
                                                                   tmp_search_trace_min_bin_occupancy, &
                                                                   tmp_search_trace_mean_bin_occupancy, &
                                                                   tmp_search_trace_max_bin_occupancy, &
                                                                   tmp_search_trace_sturges_bins, tmp_search_trace_fd_bins, &
                                                                   tmp_search_trace_shared_residual_range_low, &
                                                                   tmp_search_trace_shared_residual_range_high, &
                                                                   tmp_search_candidate_k_start, &
                                                                   tmp_search_candidate_n_points, &
                                                                   tmp_search_candidate_status, &
                                                                   tmp_candidate_confidence_interval, &
                                                                   tmp_bootstrapping_top_k_jsds, &
                                                                   tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                                   tmp_best_uncertainty_confidence_interval, &
                                                                   tmp_best_n_bins_per_point, &
                                                                   tmp_best_uncertainty_n_bins_per_point, &
                                                                   tmp_best_shared_residual_range_low, &
                                                                   tmp_best_shared_residual_range_high, &
                                                                   tmp_best_uncertainty_shared_residual_range_low, &
                                                                   tmp_best_uncertainty_shared_residual_range_high, &
                                                                   rerun_parameter_search, n_bootstraps, join_method, &
                                                                   min_residuals_per_bin, min_neighbor_overlap, &
                                                                   succeeding_ci_overlap, plateau_mode, &
                                                                   delta_median_threshold, delta_max_threshold, &
                                                                   delta_epsilon, delta_min_consecutive_transitions, m_min, &
                                                                   m_max, gamma_occupancy, lower_residual_range_quantile, &
                                                                   upper_residual_range_quantile, &
                                                                   two_sided_bootstrapping_significance_level, random_seed, &
                                                                   tau, mad_distance_factor, max_pooled_residuals, &
                                                                   min_study_neighbors, bin_count_mode, &
                                                                   delta_abs_tolerance, divergence_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: n_samples
            !! Number of samples, the non-padding slots of `expression`
            !! DM_MIN(1_int32)
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
            !! A permutation of `1..n_samples`: destination sample `k` receives source sample
            !! `sample_assignment(k)`; the identity is the observed data
            !! DM_MIN(1_int32)
            !! DM_MAX(n_samples)
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! DM_MIN(k_start)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays; a construction with more neighborhoods fails
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(1_int32)
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
            !! two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
            !! error
        integer(int32), intent(out) :: used_k_start
            !! `k_start` the construction used (`k_start`, or the search's choice)
        integer(int32), intent(out) :: used_k_step
            !! `k_step` the construction used (`k_step`, or the search's choice)
        integer(int32), intent(out) :: used_k_max
            !! `k_max` the construction used (`k_max`, or the search's choice)
        integer(int32), intent(out) :: n_points
            !! Number of neighborhoods the construction built (0 when it built none)
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        logical(c_bool), dimension(n_samples), intent(out) :: tmp_sample_mark
            !! Working array: the values of `sample_assignment` already seen
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
            !! Working array: the relabelled data set
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search, the construction and the run: sorting permutation of the
            !! pooled gene means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_x_star
            !! Working array of the search and the construction: the reference points
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_neighborhood_range
            !! Working array of the search and the construction: each neighborhood's pooled range,
            !! read by the run
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: tmp_n_neighbors_per_point
            !! Working array of the search, the construction and the run: per-study entry count per
            !! neighborhood
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_stop_reason
            !! Working array of the search and the construction: stop reason per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_dispersion
            !! Working array of the search and the construction: dispersion per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_mad
            !! Working array of the search and the construction: MAD per neighborhood
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array of the search and the run: one reference point's gene indices per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array of the search and the run: one (point, study) pair's gathered residuals
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_start
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_step
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_max
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_k_start
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_n_points
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_status
            !! Working array: the search's candidate log
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array of the search: sorting permutation for the per-study medians
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the adaptive parameter search on the relabelled data first and use its
            !! knobs; `.false.`: use `k_start`, `k_step` and `k_max`
            !! DM_DEFAULT(M_RERUN_PARAMETER_SEARCH_DEFAULT)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (at least half of the studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! DM_DEFAULT(M_METHOD_JOIN_MIN)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_EFFECT_SIZE)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(M_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! M_LOWER_RESIDUAL_RANGE_QUANTILE_ORDER_DOC
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Rerun search: seed of the bootstrap's random number generator
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_STONE_EQUAL)
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): M_DELTA_ABS_TOLERANCE_DOC
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_DELTA_ABS_TOLERANCE_DEFAULT)
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence is computed, in the run (observed and every
            !! shuffled data set) and in the rerun search
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Each study's JSD against the consensus of all studies | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]] |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! DM_DEFAULT(M_MODE_DIVERGENCE_PAIRWISE_K2)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples_found, actual_n_bootstraps, actual_join_method, search_capacity
        logical(c_bool) :: actual_rerun

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, M_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(join_method, actual_join_method, M_METHOD_JOIN_MIN)
        global_js_divergence = M_NAN
        used_k_start = k_start
        used_k_step = k_step
        used_k_max = k_max
        n_points = 0_int32
        search_capacity = 0_int32  ! only set by the checks with the rerun option

        ! Runtime checks (codegen guide Sec 5.14), each naming the argument at fault; with the rerun
        ! option they also set the search's own capacity.
        call check_permutation_input(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                     actual_rerun, max_n_points_candidate, 10_int32, tmp_sample_slots, n_samples_found, &
                                     n_samples=n_samples, sample_assignment=sample_assignment, &
                                     sample_mark=tmp_sample_mark, search_capacity=search_capacity, ierr=ierr)
        if (is_err(ierr)) return

        call run_js_comp_test_adaptive_on_sample_assignment_core(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=sample_assignment, k_start=k_start, &
            k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, search_capacity=search_capacity, &
            global_js_divergence=global_js_divergence, used_k_start=used_k_start, used_k_step=used_k_step, &
            used_k_max=used_k_max, n_points=n_points, tmp_sample_slots=tmp_sample_slots, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_occupancy_failed=tmp_occupancy_failed, &
            tmp_n_pooled_residuals=tmp_n_pooled_residuals, tmp_min_bin_occupancy=tmp_min_bin_occupancy, &
            tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, tmp_max_bin_occupancy=tmp_max_bin_occupancy, &
            tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
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
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, tmp_search_trace_delta=tmp_search_trace_delta, &
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
            actual_rerun=actual_rerun, actual_n_bootstraps=actual_n_bootstraps, actual_join_method=actual_join_method, &
            min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
            succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
            delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
            delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
            m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, &
            two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
            random_seed=random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
            max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
            bin_count_mode=bin_count_mode, delta_abs_tolerance=delta_abs_tolerance, divergence_mode=divergence_mode, &
            ierr=ierr)
    end subroutine run_js_comp_test_adaptive_on_sample_assignment_impl

    !> summary: The compute part of `run_js_comp_test_adaptive_on_sample_assignment_impl`, without its input checks
    !| AUTHOR_LASZLO_LANG
    !| Private. Assumes what the checks of the published routine ensure: `tmp_sample_slots` holds the sample slots found by
    !| `check_permutation_input`, `sample_assignment` is a permutation of `1..n_samples`, and the capacities fit
    !| the rerun search. The published routine calls it once, the whole-test routine for the observed data and every
    !| permutation, so the checks run once per test. The arguments are those of the published routine, except the resolved
    !| `actual_rerun`, `actual_n_bootstraps` and `actual_join_method` (defaults applied), and `search_capacity` (from the checks).
    subroutine run_js_comp_test_adaptive_on_sample_assignment_core(n_studies, max_n_genes_all_studies, &
                                                                   max_n_reps_all_studies, expression, n_samples, &
                                                                   sample_assignment, k_start, k_step, k_max, &
                                                                   max_n_points_candidate, n_bootstrapping_top_k_jsds, &
                                                                   search_capacity, global_js_divergence, &
                                                                   used_k_start, used_k_step, used_k_max, n_points, &
                                                                   tmp_sample_slots, tmp_permuted_expression, &
                                                                   tmp_gene_means, tmp_residuals, &
                                                                   tmp_gene_means_perm_all, tmp_x_star, &
                                                                   tmp_pooled_neighborhood_range, &
                                                                   tmp_n_neighbors_per_point, tmp_stop_reason, &
                                                                   tmp_neighborhood_dispersion, tmp_neighborhood_mad, &
                                                                   tmp_point_neighborhood_indices, &
                                                                   tmp_neighbor_residuals, tmp_point_pool, &
                                                                   tmp_point_draw, tmp_n_bins_per_point, &
                                                                   tmp_shared_residual_range_low, &
                                                                   tmp_shared_residual_range_high, &
                                                                   tmp_occupancy_failed, tmp_n_pooled_residuals, &
                                                                   tmp_min_bin_occupancy, tmp_mean_bin_occupancy, &
                                                                   tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                                   tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                                   tmp_included_n_reps, tmp_mean_pmf, &
                                                                   tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                                                   tmp_js_divergences, tmp_weights, &
                                                                   tmp_global_js_divergence, tmp_pooled_residuals, &
                                                                   tmp_pooled_residuals_perm, &
                                                                   tmp_residual_to_study_perm, tmp_bin_counts_search, &
                                                                   tmp_study_block_ends, tmp_study_block_sizes, &
                                                                   tmp_search_n_bins_per_point, &
                                                                   tmp_search_shared_residual_range_low, &
                                                                   tmp_search_shared_residual_range_high, &
                                                                   tmp_search_confidence_interval, &
                                                                   tmp_search_trace_k_start, tmp_search_trace_k_step, &
                                                                   tmp_search_trace_k_max, tmp_search_trace_n_points, &
                                                                   tmp_search_trace_global_js_divergence, &
                                                                   tmp_search_trace_ci_lower, &
                                                                   tmp_search_trace_ci_upper, &
                                                                   tmp_search_trace_ci_width, &
                                                                   tmp_search_trace_ci_width_relative, &
                                                                   tmp_search_trace_ci_overlap, &
                                                                   tmp_search_trace_delta, &
                                                                   tmp_search_trace_delta_median, &
                                                                   tmp_search_trace_delta_max, &
                                                                   tmp_search_trace_selected_n_bins, &
                                                                   tmp_search_trace_occupancy_failed, &
                                                                   tmp_search_trace_n_pooled_residuals, &
                                                                   tmp_search_trace_min_bin_occupancy, &
                                                                   tmp_search_trace_mean_bin_occupancy, &
                                                                   tmp_search_trace_max_bin_occupancy, &
                                                                   tmp_search_trace_sturges_bins, &
                                                                   tmp_search_trace_fd_bins, &
                                                                   tmp_search_trace_shared_residual_range_low, &
                                                                   tmp_search_trace_shared_residual_range_high, &
                                                                   tmp_search_candidate_k_start, &
                                                                   tmp_search_candidate_n_points, &
                                                                   tmp_search_candidate_status, &
                                                                   tmp_candidate_confidence_interval, &
                                                                   tmp_bootstrapping_top_k_jsds, &
                                                                   tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                                   tmp_best_uncertainty_confidence_interval, &
                                                                   tmp_best_n_bins_per_point, &
                                                                   tmp_best_uncertainty_n_bins_per_point, &
                                                                   tmp_best_shared_residual_range_low, &
                                                                   tmp_best_shared_residual_range_high, &
                                                                   tmp_best_uncertainty_shared_residual_range_low, &
                                                                   tmp_best_uncertainty_shared_residual_range_high, &
                                                                   actual_rerun, actual_n_bootstraps, &
                                                                   actual_join_method, min_residuals_per_bin, &
                                                                   min_neighbor_overlap, succeeding_ci_overlap, &
                                                                   plateau_mode, delta_median_threshold, &
                                                                   delta_max_threshold, delta_epsilon, &
                                                                   delta_min_consecutive_transitions, m_min, m_max, &
                                                                   gamma_occupancy, lower_residual_range_quantile, &
                                                                   upper_residual_range_quantile, &
                                                                   two_sided_bootstrapping_significance_level, &
                                                                   random_seed, tau, mad_distance_factor, &
                                                                   max_pooled_residuals, min_study_neighbors, &
                                                                   bin_count_mode, delta_abs_tolerance, &
                                                                   divergence_mode, ierr)
        integer(int32), intent(in) :: n_studies
        integer(int32), intent(in) :: max_n_genes_all_studies
        integer(int32), intent(in) :: max_n_reps_all_studies
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
        integer(int32), intent(in) :: n_samples
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
        integer(int32), intent(in) :: k_start
        integer(int32), intent(in) :: k_step
        integer(int32), intent(in) :: k_max
        integer(int32), intent(in) :: max_n_points_candidate
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
        integer(int32), intent(in) :: search_capacity
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
        integer(int32), intent(out) :: used_k_start
        integer(int32), intent(out) :: used_k_step
        integer(int32), intent(out) :: used_k_max
        integer(int32), intent(out) :: n_points
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(in) :: tmp_sample_slots
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_x_star
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_neighborhood_range
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: tmp_n_neighbors_per_point
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_stop_reason
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_dispersion
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_mad
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_start
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_step
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_max
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_k_start
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_n_points
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_status
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
        logical(c_bool), intent(in) :: actual_rerun
        integer(int32), intent(in) :: actual_n_bootstraps
        integer(int32), intent(in) :: actual_join_method
        integer(int32), intent(in), optional :: min_residuals_per_bin
        real(real64), intent(in), optional :: min_neighbor_overlap
        real(real64), intent(in), optional :: succeeding_ci_overlap
        integer(int32), intent(in), optional :: plateau_mode
        real(real64), intent(in), optional :: delta_median_threshold
        real(real64), intent(in), optional :: delta_max_threshold
        real(real64), intent(in), optional :: delta_epsilon
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
        integer(int32), intent(in), optional :: m_min
        integer(int32), intent(in), optional :: m_max
        real(real64), intent(in), optional :: gamma_occupancy
        real(real64), intent(in), optional :: lower_residual_range_quantile
        real(real64), intent(in), optional :: upper_residual_range_quantile
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
        integer(int32), intent(in), optional :: random_seed
        real(real64), intent(in), optional :: tau
        real(real64), intent(in), optional :: mad_distance_factor
        integer(int32), intent(in), optional :: max_pooled_residuals
        integer(int32), intent(in), optional :: min_study_neighbors
        integer(int32), intent(in), optional :: bin_count_mode
        real(real64), intent(in), optional :: delta_abs_tolerance
        integer(int32), intent(in), optional :: divergence_mode
        integer(int32), intent(out) :: ierr

        integer(int32) :: callee_ierr, max_n_bins_per_point, search_n_points, search_n_admissible_evaluated
        integer(int32) :: search_n_candidates_tried, max_n_neighbors, construction_status
        logical(c_bool) :: search_plateau_established

        call set_ok(ierr)
        global_js_divergence = M_NAN
        used_k_start = k_start
        used_k_step = k_step
        used_k_max = k_max
        n_points = 0_int32

        call assign_samples_and_center(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, n_samples, &
                                       sample_assignment, tmp_sample_slots, tmp_permuted_expression, tmp_gene_means, &
                                       tmp_residuals)

        ! Every capacity-sized buffer is passed whole to explicit-shape dummies declared at the
        ! sizes actually used (sequence association, no copy).
        if (actual_rerun) then
            call run_js_comp_test_adaptive_parameter_search_impl( &
                n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
                max_n_reps_all_studies=max_n_reps_all_studies, gene_means=tmp_gene_means, residuals=tmp_residuals, &
                n_bootstraps=actual_n_bootstraps, join_method=actual_join_method, max_n_points_candidate=search_capacity, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, k_start=used_k_start, k_step=used_k_step, &
                k_max=used_k_max, n_points=search_n_points, n_bins_per_point=tmp_search_n_bins_per_point, &
                shared_residual_range_low=tmp_search_shared_residual_range_low, &
                shared_residual_range_high=tmp_search_shared_residual_range_high, &
                best_candidate_confidence_interval=tmp_search_confidence_interval, &
                plateau_established=search_plateau_established, n_admissible_evaluated=search_n_admissible_evaluated, &
                trace_k_start=tmp_search_trace_k_start, trace_k_step=tmp_search_trace_k_step, &
                trace_k_max=tmp_search_trace_k_max, trace_n_points=tmp_search_trace_n_points, &
                trace_global_js_divergence=tmp_search_trace_global_js_divergence, &
                trace_ci_lower=tmp_search_trace_ci_lower, trace_ci_upper=tmp_search_trace_ci_upper, &
                trace_ci_width=tmp_search_trace_ci_width, trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
                trace_ci_overlap=tmp_search_trace_ci_overlap, trace_delta=tmp_search_trace_delta, &
                trace_delta_median=tmp_search_trace_delta_median, trace_delta_max=tmp_search_trace_delta_max, &
                trace_selected_n_bins=tmp_search_trace_selected_n_bins, &
                trace_occupancy_failed=tmp_search_trace_occupancy_failed, &
                trace_n_pooled_residuals=tmp_search_trace_n_pooled_residuals, &
                trace_min_bin_occupancy=tmp_search_trace_min_bin_occupancy, &
                trace_mean_bin_occupancy=tmp_search_trace_mean_bin_occupancy, &
                trace_max_bin_occupancy=tmp_search_trace_max_bin_occupancy, &
                trace_sturges_bins=tmp_search_trace_sturges_bins, trace_fd_bins=tmp_search_trace_fd_bins, &
                trace_shared_residual_range_low=tmp_search_trace_shared_residual_range_low, &
                trace_shared_residual_range_high=tmp_search_trace_shared_residual_range_high, &
                n_candidates_tried=search_n_candidates_tried, candidate_k_start=tmp_search_candidate_k_start, &
                candidate_n_points=tmp_search_candidate_n_points, candidate_status=tmp_search_candidate_status, &
                tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
                tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
                tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
                tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
                tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
                tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_n_bins_per_point=tmp_n_bins_per_point, &
                tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
                tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
                tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
                tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
                tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
                tmp_global_js_divergence=tmp_global_js_divergence, tmp_point_pool=tmp_point_pool, &
                tmp_point_draw=tmp_point_draw, tmp_confidence_interval=tmp_candidate_confidence_interval, &
                tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
                tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_delta_perm, &
                tmp_best_uncertainty_confidence_interval=tmp_best_uncertainty_confidence_interval, &
                tmp_pooled_residuals=tmp_pooled_residuals, tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
                tmp_bin_counts_search=tmp_bin_counts_search, tmp_residual_to_study_perm=tmp_residual_to_study_perm, &
                tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
                tmp_occupancy_failed=tmp_occupancy_failed, tmp_n_pooled_residuals=tmp_n_pooled_residuals, &
                tmp_min_bin_occupancy=tmp_min_bin_occupancy, tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, &
                tmp_max_bin_occupancy=tmp_max_bin_occupancy, tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, &
                tmp_best_n_bins_per_point=tmp_best_n_bins_per_point, &
                tmp_best_uncertainty_n_bins_per_point=tmp_best_uncertainty_n_bins_per_point, &
                tmp_best_shared_residual_range_low=tmp_best_shared_residual_range_low, &
                tmp_best_shared_residual_range_high=tmp_best_shared_residual_range_high, &
                tmp_best_uncertainty_shared_residual_range_low=tmp_best_uncertainty_shared_residual_range_low, &
                tmp_best_uncertainty_shared_residual_range_high=tmp_best_uncertainty_shared_residual_range_high, &
                min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
                succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
                delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
                delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
                m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
                lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
                max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
                bin_count_mode=bin_count_mode, &
                delta_abs_tolerance=delta_abs_tolerance, divergence_mode=divergence_mode, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                return
            end if
        end if

        call construct_adaptive_neighborhoods_impl( &
            n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
            max_n_reps_all_studies=max_n_reps_all_studies, gene_means=tmp_gene_means, residuals=tmp_residuals, &
            k_start=used_k_start, k_step=used_k_step, k_max=used_k_max, n_points=n_points, x_star=tmp_x_star, &
            pooled_neighborhood_range=tmp_pooled_neighborhood_range, n_neighbors_per_point=tmp_n_neighbors_per_point, &
            stop_reason=tmp_stop_reason, neighborhood_dispersion=tmp_neighborhood_dispersion, &
            neighborhood_mad=tmp_neighborhood_mad, max_n_neighbors=max_n_neighbors, &
            construction_status=construction_status, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tau=tau, &
            mad_distance_factor=mad_distance_factor, max_pooled_residuals=max_pooled_residuals, &
            min_study_neighbors=min_study_neighbors)
        ! Not one argument's fault: the relabelled data decide how many neighborhoods emerge and
        ! whether the construction succeeds (any status other than ERR_OK is a failed construction,
        ! as the adaptive search counts it).
        if (construction_status /= ERR_OK .or. n_points < 1_int32 .or. n_points > max_n_points_candidate) then
            call set_err_once(ierr, ERR_INVALID_INPUT)
            return
        end if

        call run_js_comp_test_adaptive_impl( &
            n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
            max_n_reps_all_studies=max_n_reps_all_studies, n_points=n_points, gene_means=tmp_gene_means, &
            residuals=tmp_residuals, pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            n_neighbors_per_point=tmp_n_neighbors_per_point, n_bins_per_point=tmp_n_bins_per_point, &
            shared_residual_range_low=tmp_shared_residual_range_low, &
            shared_residual_range_high=tmp_shared_residual_range_high, max_n_bins_per_point=max_n_bins_per_point, &
            occupancy_failed=tmp_occupancy_failed, n_pooled_residuals=tmp_n_pooled_residuals, &
            min_bin_occupancy=tmp_min_bin_occupancy, mean_bin_occupancy=tmp_mean_bin_occupancy, &
            max_bin_occupancy=tmp_max_bin_occupancy, sturges_bins=tmp_sturges_bins, fd_bins=tmp_fd_bins, pmfs=tmp_pmfs, &
            counts=tmp_counts, included_n_reps=tmp_included_n_reps, mean_pmf=tmp_mean_pmf, &
            mean_pmf_counts=tmp_mean_pmf_counts, mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            js_divergences=tmp_js_divergences, weights=tmp_weights, global_js_divergence=global_js_divergence, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, &
            tmp_bin_counts_search=tmp_bin_counts_search, tmp_residual_to_study_perm=tmp_residual_to_study_perm, &
            tmp_study_block_ends=tmp_study_block_ends, tmp_study_block_sizes=tmp_study_block_sizes, &
            min_residuals_per_bin=min_residuals_per_bin, m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, bin_count_mode=bin_count_mode, &
            divergence_mode=divergence_mode, ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence = M_NAN
            return
        end if
    end subroutine run_js_comp_test_adaptive_on_sample_assignment_core

    !> summary: Permutation-Test of the adaptive JS-Comp-Test that permutes the samples' study labels
    !| AUTHOR_LASZLO_LANG
    !| The adaptive counterpart of
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]]:
    !| the same test, each data set run by
    !| [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| with the same data, knobs and options.
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
    !|    `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
    !|    shuffled data. Row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.
    !|
    !| With `rerun_parameter_search = .true.` every data set, the observed one included, first runs the adaptive parameter search. `random_seed` then
    !| also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
    !| draws. The routine is sequential: a caller who wants parallel permutations calls
    !| [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module):run_js_comp_test_adaptive_on_sample_assignment_impl(interface)]]
    !| itself, once per assignment.
    !|
    !| Failures: a permutation that fails with invalid input (a study left without residuals at some reference point, or a construction that does not
    !| succeed: a status other than ok, no neighborhood, or more neighborhoods than the work arrays hold) gets a NaN row in `jsd_null`, is
    !| counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
    !| `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
    !| returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
    !| is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own fails the same way with the rerun option.
    !| Known limitation: a gene missing in a whole study travels with the samples.
    subroutine run_js_comp_test_adaptive_permutation_impl(n_studies, max_n_genes_all_studies, &
                                                          max_n_reps_all_studies, expression, k_start, &
                                                          k_step, k_max, max_n_points_candidate, &
                                                          n_bootstrapping_top_k_jsds, n_permutations, &
                                                          global_js_divergence_observed, p_values, &
                                                          jsd_null, n_failed_permutations, &
                                                          tmp_sample_slots, tmp_sample_assignment, &
                                                          tmp_permuted_global_js_divergence, &
                                                          tmp_permuted_expression, &
                                                          tmp_gene_means, tmp_residuals, &
                                                          tmp_gene_means_perm_all, tmp_x_star, &
                                                          tmp_pooled_neighborhood_range, &
                                                          tmp_n_neighbors_per_point, tmp_stop_reason, &
                                                          tmp_neighborhood_dispersion, &
                                                          tmp_neighborhood_mad, &
                                                          tmp_point_neighborhood_indices, &
                                                          tmp_neighbor_residuals, tmp_point_pool, &
                                                          tmp_point_draw, tmp_n_bins_per_point, &
                                                          tmp_shared_residual_range_low, &
                                                          tmp_shared_residual_range_high, &
                                                          tmp_occupancy_failed, tmp_n_pooled_residuals, &
                                                          tmp_min_bin_occupancy, tmp_mean_bin_occupancy, &
                                                          tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                          tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                          tmp_included_n_reps, tmp_mean_pmf, &
                                                          tmp_mean_pmf_counts, &
                                                          tmp_mean_pmf_included_n_reps, &
                                                          tmp_js_divergences, tmp_weights, &
                                                          tmp_global_js_divergence, tmp_pooled_residuals, &
                                                          tmp_pooled_residuals_perm, &
                                                          tmp_residual_to_study_perm, &
                                                          tmp_bin_counts_search, tmp_study_block_ends, &
                                                          tmp_study_block_sizes, tmp_search_n_bins_per_point, &
                                                          tmp_search_shared_residual_range_low, &
                                                          tmp_search_shared_residual_range_high, &
                                                          tmp_search_confidence_interval, &
                                                          tmp_search_trace_k_start, &
                                                          tmp_search_trace_k_step, tmp_search_trace_k_max, &
                                                          tmp_search_trace_n_points, &
                                                          tmp_search_trace_global_js_divergence, &
                                                          tmp_search_trace_ci_lower, &
                                                          tmp_search_trace_ci_upper, &
                                                          tmp_search_trace_ci_width, &
                                                          tmp_search_trace_ci_width_relative, &
                                                          tmp_search_trace_ci_overlap, &
                                                          tmp_search_trace_delta, &
                                                          tmp_search_trace_delta_median, &
                                                          tmp_search_trace_delta_max, &
                                                          tmp_search_trace_selected_n_bins, &
                                                          tmp_search_trace_occupancy_failed, &
                                                          tmp_search_trace_n_pooled_residuals, &
                                                          tmp_search_trace_min_bin_occupancy, &
                                                          tmp_search_trace_mean_bin_occupancy, &
                                                          tmp_search_trace_max_bin_occupancy, &
                                                          tmp_search_trace_sturges_bins, &
                                                          tmp_search_trace_fd_bins, &
                                                          tmp_search_trace_shared_residual_range_low, &
                                                          tmp_search_trace_shared_residual_range_high, &
                                                          tmp_search_candidate_k_start, &
                                                          tmp_search_candidate_n_points, &
                                                          tmp_search_candidate_status, &
                                                          tmp_candidate_confidence_interval, &
                                                          tmp_bootstrapping_top_k_jsds, &
                                                          tmp_prev_global_js_divergence, tmp_delta_perm, &
                                                          tmp_best_uncertainty_confidence_interval, &
                                                          tmp_best_n_bins_per_point, &
                                                          tmp_best_uncertainty_n_bins_per_point, &
                                                          tmp_best_shared_residual_range_low, &
                                                          tmp_best_shared_residual_range_high, &
                                                          tmp_best_uncertainty_shared_residual_range_low, &
                                                          tmp_best_uncertainty_shared_residual_range_high, &
                                                          rerun_parameter_search, n_bootstraps, &
                                                          join_method, min_residuals_per_bin, &
                                                          min_neighbor_overlap, succeeding_ci_overlap, &
                                                          plateau_mode, delta_median_threshold, &
                                                          delta_max_threshold, delta_epsilon, &
                                                          delta_min_consecutive_transitions, m_min, m_max, &
                                                          gamma_occupancy, lower_residual_range_quantile, &
                                                          upper_residual_range_quantile, &
                                                          two_sided_bootstrapping_significance_level, &
                                                          random_seed, tau, mad_distance_factor, &
                                                          max_pooled_residuals, min_study_neighbors, &
                                                          bin_count_mode, &
                                                          delta_abs_tolerance, divergence_mode, ierr)
        integer(int32), intent(in) :: n_studies
            !! Number of studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_genes_all_studies
            !! Maximum number of genes across all studies
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: max_n_reps_all_studies
            !! Maximum number of replicates (sample slots) across all studies
            !! DM_MIN(1_int32)
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
            !! Expression of every sample slot of every gene in every study (see above); NaN for a
            !! missing value, a missing gene or a padding slot
            !! DM_ALLOW_NAN
        integer(int32), intent(in) :: k_start
            !! Pooled entries every neighborhood takes unconditionally (also sizes the work arrays)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_step
            !! Pooled entries staged per adaptive growth round
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: k_max
            !! Largest number of pooled entries in one neighborhood
            !! DM_MIN(k_start)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays; a construction with more neighborhoods fails
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_adaptive_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_permutations
            !! Number of random assignments (rows of `jsd_null`)
            !! DM_MIN(1_int32)
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence_observed
            !! Weighted global JSD of each study on the observed data (the identity assignment)
        real(real64), dimension(n_studies), intent(out) :: p_values
            !! Empirical p-value per study, `(count + 1)/(n_valid + 1)` over the valid permutations
            !! whose value reached the observed one; NaN when none is valid or on error
        real(real64), dimension(n_permutations, n_studies), intent(out) :: jsd_null
            !! Row `i`: permutation `i`'s global JSD per study; NaN for a failed permutation
        integer(int32), intent(out) :: n_failed_permutations
            !! Number of permutations whose JS-Comp-Test failed with invalid input
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_slots
            !! Working array: the flat slot index of every sample
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(out) :: tmp_sample_assignment
            !! Working array: the current assignment, in its first `n_samples` elements
        real(real64), dimension(n_studies), intent(out) :: tmp_permuted_global_js_divergence
            !! Working array: one permutation's global JSD per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
            !! Working array: the relabelled data set
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
            !! Working array: the relabelled data set's gene means per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
            !! Working array: the relabelled data set's residuals
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search, the construction and the run: sorting permutation of the
            !! pooled gene means
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_x_star
            !! Working array of the search and the construction: the reference points
        integer(int32), dimension(2, max_n_genes_all_studies*n_studies), intent(out) :: tmp_pooled_neighborhood_range
            !! Working array of the search and the construction: each neighborhood's pooled range,
            !! read by the run
        integer(int32), dimension(n_studies, max_n_genes_all_studies*n_studies), intent(out) :: tmp_n_neighbors_per_point
            !! Working array of the search, the construction and the run: per-study entry count per
            !! neighborhood
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_stop_reason
            !! Working array of the search and the construction: stop reason per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_dispersion
            !! Working array of the search and the construction: dispersion per neighborhood
        real(real64), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_neighborhood_mad
            !! Working array of the search and the construction: MAD per neighborhood
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_point_neighborhood_indices
            !! Working array of the search and the run: one reference point's gene indices per study
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies), intent(out) :: tmp_neighbor_residuals
            !! Working array of the search and the run: one (point, study) pair's gathered residuals
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_pool
            !! Working array of the search and the run: one point's pool
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_point_draw
            !! Working array of the search and the run: one study's draw at one point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_bins_per_point
            !! Working array of the search and the run: bin count per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_low
            !! Working array of the search and the run: lower residual-range bound per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_shared_residual_range_high
            !! Working array of the search and the run: upper residual-range bound per reference point
        logical(c_bool), dimension(max_n_points_candidate), intent(out) :: tmp_occupancy_failed
            !! Working array of the search and the run: bin-count search failure flag per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_n_pooled_residuals
            !! Working array of the search and the run: pooled residual count per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_min_bin_occupancy
            !! Working array of the search and the run: minimum bin occupancy per reference point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_mean_bin_occupancy
            !! Working array of the search and the run: mean bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_max_bin_occupancy
            !! Working array of the search and the run: maximum bin occupancy per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_sturges_bins
            !! Working array of the search and the run: Sturges' rule diagnostic per reference point
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_fd_bins
            !! Working array of the search and the run: Freedman-Diaconis diagnostic per reference point
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_pmfs
            !! Working array of the search and the run: every study's pmfs
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS, n_studies), intent(out) :: tmp_counts
            !! Working array of the search and the run: every study's histogram counts
        integer(int32), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_included_n_reps
            !! Working array of the search and the run: every study's included residual counts
        real(real64), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf
            !! Working array of the search and the run: the consensus pmf
        integer(int32), dimension(max_n_points_candidate, M_MAX_N_BINS), intent(out) :: tmp_mean_pmf_counts
            !! Working array of the search and the run: the consensus histogram counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_mean_pmf_included_n_reps
            !! Working array of the search and the run: the consensus included residual counts
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_js_divergences
            !! Working array of the search and the run: per-point JSD values
        real(real64), dimension(max_n_points_candidate, n_studies), intent(out) :: tmp_weights
            !! Working array of the search and the run: per-point weights
        real(real64), dimension(n_studies), intent(out) :: tmp_global_js_divergence
            !! Working array of the search and the run: global JSD scratch
        real(real64), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_genes_all_studies*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_start
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_step
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_k_max
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_global_js_divergence
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_lower
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_upper
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_width_relative
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_ci_overlap
            !! Working array: the search's trace
        real(real64), dimension(n_studies, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_median
            !! Working array: the search's trace
        real(real64), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_delta_max
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_selected_n_bins
            !! Working array: the search's trace
        logical(c_bool), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_occupancy_failed
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_n_pooled_residuals
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_min_bin_occupancy
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_mean_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_max_bin_occupancy
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_sturges_bins
            !! Working array: the search's trace
        integer(int32), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_fd_bins
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_low
            !! Working array: the search's trace
        real(real64), dimension(max_n_points_candidate, M_MAX_CANDIDATE_PAIRS), intent(out) :: &
            tmp_search_trace_shared_residual_range_high
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_k_start
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_n_points
            !! Working array: the search's candidate log
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_candidate_status
            !! Working array: the search's candidate log
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
        integer(int32), dimension(n_studies), intent(out) :: tmp_delta_perm
            !! Working array of the search: sorting permutation for the per-study medians
        real(real64), dimension(2, n_studies), intent(out) :: tmp_best_uncertainty_confidence_interval
            !! Working array of the search: the smallest-uncertainty candidate's confidence interval
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_n_bins_per_point
            !! Working array of the search: the best candidate's bin counts
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_n_bins_per_point
            !! Working array of the search: the smallest-uncertainty candidate's bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_low
            !! Working array of the search: the best candidate's lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_shared_residual_range_high
            !! Working array of the search: the best candidate's upper residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_low
            !! Working array of the search: the smallest-uncertainty candidate's lower bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_best_uncertainty_shared_residual_range_high
            !! Working array of the search: the smallest-uncertainty candidate's upper bounds
        logical(c_bool), intent(in), optional :: rerun_parameter_search
            !! `.true.`: run the adaptive parameter search on the relabelled data first and use its
            !! knobs; `.false.`: use `k_start`, `k_step` and `k_max`
            !! DM_DEFAULT(M_RERUN_PARAMETER_SEARCH_DEFAULT)
        integer(int32), intent(in), optional :: n_bootstraps
            !! Number of bootstraps per candidate of the rerun parameter search
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        integer(int32), intent(in), optional :: join_method
            !! How the rerun search joins the studies' confidence-interval overlaps
            !!
            !! | Method | Value |
            !! |--------|-------|
            !! | Minimum overlap (all studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MIN(variable)]] |
            !! | Maximum overlap (any one study passes) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MAX(variable)]] |
            !! | Median overlap (at least half of the studies must pass) | [[tox_data_integration_js_comp_test_search_impl(module):METHOD_JOIN_MEDIAN(variable)]] |
            !! DM_DEFAULT(M_METHOD_JOIN_MIN)
        integer(int32), intent(in), optional :: min_residuals_per_bin
            !! Minimum number of pooled residuals every bin must reach in the bin-count search (and,
            !! in the rerun search, in every consensus bin)
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_OCCUPANCY_MIN_RESIDUALS_PER_BIN_DEFAULT)
        real(real64), intent(in), optional :: min_neighbor_overlap
            !! Rerun search: minimum fractional overlap of two consecutive neighborhoods (gate 1)
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_MIN_NEIGHBOR_OVERLAP_DEFAULT)
        real(real64), intent(in), optional :: succeeding_ci_overlap
            !! Rerun search: minimum fractional confidence-interval overlap with the best candidate so far
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_SUCCEEDING_CI_OVERLAP_DEFAULT)
        integer(int32), intent(in), optional :: plateau_mode
            !! Rerun search: which plateau criterion stops the search; by default, as in the
            !! searches, the change of the observed JSD values from one admissible candidate to the next
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | CI overlap only | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_CI_OVERLAP(variable)]] |
            !! | Relative-effect-size stability only: the change of the observed JSD values (default) | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_EFFECT_SIZE(variable)]] |
            !! | Either criterion | [[tox_data_integration_js_comp_test_search_impl(module):MODE_PLATEAU_BOTH(variable)]] |
            !! DM_DEFAULT(M_MODE_PLATEAU_EFFECT_SIZE)
        real(real64), intent(in), optional :: delta_median_threshold
            !! Rerun search: upper bound on the median relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MEDIAN_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_max_threshold
            !! Rerun search: upper bound on the largest relative JSD change for a transition to qualify
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_MAX_THRESHOLD_DEFAULT)
        real(real64), intent(in), optional :: delta_epsilon
            !! Rerun search: floor of the denominators of the relative JSD change and CI width
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_DELTA_EPSILON_DEFAULT)
        integer(int32), intent(in), optional :: delta_min_consecutive_transitions
            !! Rerun search: number of consecutive qualifying transitions for an effect-size plateau
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_DELTA_MIN_CONSECUTIVE_TRANSITIONS_DEFAULT)
        integer(int32), intent(in), optional :: m_min
            !! Smallest bin count the bin-count search tests (M_min)
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MIN_DEFAULT)
        integer(int32), intent(in), optional :: m_max
            !! Largest bin count the bin-count search tests (M_max); raised to `m_min` if smaller
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_N_BINS)
            !! DM_DEFAULT(M_OCCUPANCY_M_MAX_DEFAULT)
        real(real64), intent(in), optional :: gamma_occupancy
            !! Geometric growth factor of the occupancy search; unused in the Stone modes
            !! DM_MIN(above(1.0_real64))
            !! DM_DEFAULT(M_OCCUPANCY_GAMMA_DEFAULT)
        real(real64), intent(in), optional :: lower_residual_range_quantile
            !! Quantile in [0,1] of each reference point's lower residual-range bound
            !! M_LOWER_RESIDUAL_RANGE_QUANTILE_ORDER_DOC
            !! DM_MIN(0.0_real64)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_LOWER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        ! The lower bound below (lower_residual_range_quantile) is inclusive on purpose: above() must not be
        ! used, because the wrapper would evaluate above() with an absent optional when a Fortran
        ! caller omits lower_residual_range_quantile.
        real(real64), intent(in), optional :: upper_residual_range_quantile
            !! Quantile in [0,1] of each reference point's upper residual-range bound
            !! DM_MIN(lower_residual_range_quantile)
            !! DM_MAX(1.0_real64)
            !! DM_DEFAULT(M_OCCUPANCY_UPPER_RESIDUAL_RANGE_QUANTILE_DEFAULT)
        real(real64), intent(in), optional :: two_sided_bootstrapping_significance_level
            !! Rerun search: two-sided significance level of the bootstrap interval, in percent per end
            !! DM_MIN(0.0_real64)
            !! DM_MAX(100.0_real64)
            !! DM_DEFAULT(M_TWO_SIDED_BOOTSTRAPPING_SIGNIFICANCE_LEVEL_DEFAULT)
        integer(int32), intent(in), optional :: random_seed
            !! Seed of the random number generator of the label draws, and of the bootstrap of
            !! every rerun parameter search (a separate generator)
            !! DM_DEFAULT(M_RANDOM_SEED_DEFAULT)
        real(real64), intent(in), optional :: tau
            !! Largest relative dispersion increase a growth round may cause and still be committed
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_TAU_DEFAULT)
        real(real64), intent(in), optional :: mad_distance_factor
            !! Multiple of a neighborhood's MAD that the next reference point's target lies beyond it
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_ADAPTIVE_MAD_DISTANCE_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: max_pooled_residuals
            !! Cap on a neighborhood's non-NaN residuals in the growth rounds; 0 for no cap
            !! DM_MIN(0_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MAX_POOLED_RESIDUALS_DEFAULT)
        integer(int32), intent(in), optional :: min_study_neighbors
            !! Fewest entries of each study every neighborhood must have (construction status; in
            !! the rerun search, candidate admissibility)
            !! DM_MIN(1_int32)
            !! DM_DEFAULT(M_ADAPTIVE_MIN_STUDY_NEIGHBORS_DEFAULT)
        integer(int32), intent(in), optional :: bin_count_mode
            !! How each reference point's histogram bin count is chosen; the Stone modes test every
            !! bin count from `m_min` to `m_max`
            !! ([[tox_data_integration_js_comp_test_binning_impl(module):determine_bin_count_stone_impl(interface)]]),
            !! with `min_residuals_per_bin` as the floor
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Occupancy search: the largest bin count whose pooled histogram has every bin at or above `min_residuals_per_bin` | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_OCCUPANCY(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's own score (default) | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_EQUAL(variable)]] |
            !! | Stone's cross-validation: minimise the sum of every study's score multiplied by its residual count inside the range | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_SIZE_WEIGHTED(variable)]] |
            !! | Stone's cross-validation: minimise the score of all studies' residuals pooled | [[tox_data_integration_js_comp_test_binning_impl(module):MODE_BIN_COUNT_STONE_POOLED(variable)]] |
            !! DM_DEFAULT(M_MODE_BIN_COUNT_STONE_EQUAL)
        real(real64), intent(in), optional :: delta_abs_tolerance
            !! Rerun search (used only with `rerun_parameter_search`): M_DELTA_ABS_TOLERANCE_DOC
            !! DM_MIN(0.0_real64)
            !! DM_DEFAULT(M_DELTA_ABS_TOLERANCE_DEFAULT)
        integer(int32), intent(in), optional :: divergence_mode
            !! How each study's per-point divergence is computed, in the run (observed and every
            !! shuffled data set) and in the rerun search
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Each study's JSD against the consensus of all studies | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_SYMMETRIC(variable)]] |
            !! | Two studies: their pairwise JSD, the same for both; any other number of studies: as the symmetric mode (default) | [[tox_data_integration_jsd_impl(module):MODE_DIVERGENCE_PAIRWISE_K2(variable)]] |
            !! DM_DEFAULT(M_MODE_DIVERGENCE_PAIRWISE_K2)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples, i_permutation, actual_random_seed, callee_ierr, actual_n_bootstraps, actual_join_method
        integer(int32) :: used_k_start, used_k_step, used_k_max, n_points, search_capacity
        logical(c_bool) :: actual_rerun
        type(rng_t) :: rng

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, M_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(random_seed, actual_random_seed, M_RANDOM_SEED_DEFAULT)
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(join_method, actual_join_method, M_METHOD_JOIN_MIN)
        global_js_divergence_observed = M_NAN
        p_values = M_NAN
        jsd_null = M_NAN
        n_failed_permutations = 0_int32
        search_capacity = 0_int32  ! only set by the checks with the rerun option

        ! The checks that do not depend on the assignment run once here, with this routine's argument
        ! positions. Every data set below runs the core, which does not repeat them.
        call check_permutation_input(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                     actual_rerun, max_n_points_candidate, 8_int32, tmp_sample_slots, n_samples, &
                                     search_capacity=search_capacity, ierr=ierr)
        if (is_err(ierr)) return

        ! The observed data: the identity assignment. An error here ends the test at once.
        call init_perm(tmp_sample_assignment(1:n_samples))
        call run_js_comp_test_adaptive_on_sample_assignment_core(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, k_start=k_start, &
            k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, search_capacity=search_capacity, &
            global_js_divergence=global_js_divergence_observed, used_k_start=used_k_start, used_k_step=used_k_step, &
            used_k_max=used_k_max, n_points=n_points, tmp_sample_slots=tmp_sample_slots, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_pooled_neighborhood_range=tmp_pooled_neighborhood_range, &
            tmp_n_neighbors_per_point=tmp_n_neighbors_per_point, tmp_stop_reason=tmp_stop_reason, &
            tmp_neighborhood_dispersion=tmp_neighborhood_dispersion, tmp_neighborhood_mad=tmp_neighborhood_mad, &
            tmp_point_neighborhood_indices=tmp_point_neighborhood_indices, &
            tmp_neighbor_residuals=tmp_neighbor_residuals, tmp_point_pool=tmp_point_pool, &
            tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
            tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
            tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_occupancy_failed=tmp_occupancy_failed, &
            tmp_n_pooled_residuals=tmp_n_pooled_residuals, tmp_min_bin_occupancy=tmp_min_bin_occupancy, &
            tmp_mean_bin_occupancy=tmp_mean_bin_occupancy, tmp_max_bin_occupancy=tmp_max_bin_occupancy, &
            tmp_sturges_bins=tmp_sturges_bins, tmp_fd_bins=tmp_fd_bins, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
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
            tmp_search_trace_ci_lower=tmp_search_trace_ci_lower, tmp_search_trace_ci_upper=tmp_search_trace_ci_upper, &
            tmp_search_trace_ci_width=tmp_search_trace_ci_width, &
            tmp_search_trace_ci_width_relative=tmp_search_trace_ci_width_relative, &
            tmp_search_trace_ci_overlap=tmp_search_trace_ci_overlap, tmp_search_trace_delta=tmp_search_trace_delta, &
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
            actual_rerun=actual_rerun, actual_n_bootstraps=actual_n_bootstraps, actual_join_method=actual_join_method, &
            min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
            succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
            delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
            delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
            m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, &
            two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
            random_seed=actual_random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
            max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
            bin_count_mode=bin_count_mode, delta_abs_tolerance=delta_abs_tolerance, divergence_mode=divergence_mode, &
            ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence_observed = M_NAN
            return
        end if

        rng = create_rng(callee_ierr, actual_random_seed)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            call destroy_rng(rng)
            return
        end if
        ! Sequential: one generator draws every assignment in turn. It exists, and n_samples >= 1 after
        ! the check above (every study has a sample).
        ! Each assignment is a shuffle of the identity, so it depends only on the generator state and not on the
        ! previous assignment. The reset costs n_samples steps, negligible next to one test run.
        do i_permutation = 1, n_permutations
            call init_perm(tmp_sample_assignment(1:n_samples))
            call random_shuffle(rng, n_samples, tmp_sample_assignment(1:n_samples))
            call run_js_comp_test_adaptive_on_sample_assignment_core(n_studies=n_studies, &
                max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
                expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, k_start=k_start, &
                k_step=k_step, k_max=k_max, max_n_points_candidate=max_n_points_candidate, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, search_capacity=search_capacity, &
                global_js_divergence=tmp_permuted_global_js_divergence, used_k_start=used_k_start, &
                used_k_step=used_k_step, used_k_max=used_k_max, n_points=n_points, tmp_sample_slots=tmp_sample_slots, &
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
                actual_rerun=actual_rerun, actual_n_bootstraps=actual_n_bootstraps, &
                actual_join_method=actual_join_method, min_residuals_per_bin=min_residuals_per_bin, &
                min_neighbor_overlap=min_neighbor_overlap, succeeding_ci_overlap=succeeding_ci_overlap, &
                plateau_mode=plateau_mode, delta_median_threshold=delta_median_threshold, &
                delta_max_threshold=delta_max_threshold, delta_epsilon=delta_epsilon, &
                delta_min_consecutive_transitions=delta_min_consecutive_transitions, m_min=m_min, m_max=m_max, &
                gamma_occupancy=gamma_occupancy, lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=actual_random_seed, tau=tau, mad_distance_factor=mad_distance_factor, &
                max_pooled_residuals=max_pooled_residuals, min_study_neighbors=min_study_neighbors, &
                bin_count_mode=bin_count_mode, delta_abs_tolerance=delta_abs_tolerance, &
                divergence_mode=divergence_mode, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                ! Invalid input of a relabelled data set (e.g. a study without residuals at some
                ! reference point) fails this permutation only; anything else ends the test.
                if (get_err_code(callee_ierr) == ERR_INVALID_INPUT) then
                    n_failed_permutations = n_failed_permutations + 1_int32
                    cycle
                end if
                call set_err_once(ierr, get_err_code(callee_ierr))
                exit
            end if
            jsd_null(i_permutation, :) = tmp_permuted_global_js_divergence
        end do
        call destroy_rng(rng)
        if (is_err(ierr)) return

        call label_permutation_p_values(n_permutations, n_studies, global_js_divergence_observed, jsd_null, p_values)
    end subroutine run_js_comp_test_adaptive_permutation_impl
end module tox_data_integration_js_comp_test_permutation_adaptive_impl
