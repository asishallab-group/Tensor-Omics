#include <src/macros.h>
#include <src/tox/data_integration/macros.h>

!> # Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Study-Label Permutation-Test, k-NN
!|
!| Asks whether the assignment of samples to studies (batches) influences the JS-Comp-Test values
!| at all. Null hypothesis: it does not. The test shuffles the study labels of the samples
!| (without replacement), reruns the whole JS-Comp-Test on every shuffled data set, and compares
!| each study's observed global JSD with the values the shuffled data sets give. Study `i`'s p-value is
!| `(count + 1)/(n_valid + 1)`, where `count` is how many of the `n_valid` valid shuffled data sets give
!| study `i` a global JSD at least as large as the observed one. The +1 counts the observed data itself.
!|
!| **Input layout.** All studies arrive in one array `expression(max_n_reps_all_studies,
!| max_n_genes_all_studies, n_studies)`: gene `g` is the same gene in every study, a gene missing in
!| a study is NaN there, and a replicate slot `(r, s)` whose values are all NaN is padding (a study
!| with fewer samples than `max_n_reps_all_studies`). Every other slot is a sample. Samples are
!| numbered in slot order, replicate fastest, then study.
!|
!| **Entry points.** One assignment of samples to slots, one k-NN JS-Comp-Test (the fixed
!| `(n_points, n_neighbors)` of `run_js_comp_test`):
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]].
!| The whole test, `n_permutations` random assignments plus the observed one:
!| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_permutation_impl(interface)]].
!| The adaptive counterparts are in [[tox_data_integration_js_comp_test_permutation_adaptive_impl(module)]], the shared helpers in [[tox_data_integration_js_comp_test_permutation_common_impl(module)]].
!| `calc_js_comp_test_permutation_bounds` and `calc_adaptive_js_comp_test_permutation_bounds`
!| ([[tox_data_integration_js_comp_test_permutation_bounds_impl(module)]]) size the work arrays for a caller that allocates its own.
!|
!| The adaptive routines are the recommended way, as for the JS-Comp-Test itself; the k-NN ones
!| stay available.
module tox_data_integration_js_comp_test_permutation_fixed_k_impl
    use f42_safeguard
    use, intrinsic :: iso_fortran_env, only: int32, real64
    use, intrinsic :: iso_c_binding, only: c_bool
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use tox_errors, only: set_ok, set_err_once, is_err, get_err_code, ERR_INVALID_INPUT
    use f42_random_gsl, only: rng_t, create_rng, destroy_rng, random_shuffle
    use f42_sort_impl, only: init_perm
    use tox_data_integration_preprocessing_impl, only: pool_means_impl, sort_gene_means_per_study, &
                                                       sort_pooled_gene_means
    use tox_data_integration_js_comp_test_fixed_k_search_impl, only: run_js_comp_test_parameter_search_impl
    use tox_data_integration_js_comp_test_run_impl, only: run_js_comp_test_impl
    use tox_data_integration_js_comp_test_permutation_common_impl, only: check_permutation_input, &
                                                                         assign_samples_and_center, &
                                                                         label_permutation_p_values
    M_IMPLICIT_NONE
    private
    public :: run_js_comp_test_on_sample_assignment_impl, run_js_comp_test_permutation_impl

contains

    !> summary: Run the k-NN JS-Comp-Test on one assignment of samples to studies
    !| AUTHOR_LASZLO_LANG
    !| One building block of the study-label Permutation-Test: relabels the samples by `sample_assignment` and runs
    !| the whole JS-Comp-Test on the relabelled data set. `expression` holds every study (input layout as in the module
    !| description). `sample_assignment` is a permutation of `1..n_samples`: destination sample `k` receives all genes of source
    !| sample `sample_assignment(k)`, so every study keeps its slots, sample count and padding. The identity is the observed
    !| data, and its result is exactly that of a caller doing the steps below on the observed data.
    !|
    !| Steps: per study, the gene means and residuals of the relabelled data
    !| ([[tox_data_integration_preprocessing_impl(module):compute_gene_means_impl(interface)]], [[tox_data_integration_preprocessing_impl(module):compute_residuals_impl(interface)]]); the reference points from
    !| the pooled means ([[tox_data_integration_preprocessing_impl(module):pool_means_impl(interface)]]); then [[tox_data_integration_js_comp_test_run_impl(module):run_js_comp_test_impl(interface)]] with the
    !| bin-count options and `divergence_mode` as given.
    !|
    !| With `rerun_parameter_search = .true.` the relabelled data first run
    !| [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):run_js_comp_test_parameter_search_impl(interface)]]
    !| with the search options below, and the run uses the `(n_points, n_neighbors)` it returns (also its fallback when no
    !| candidate was admissible). `n_points` and `n_neighbors` then only size the work arrays, and `used_n_points` and
    !| `used_n_neighbors` report what the run used. Every call runs the whole JS-Comp-Test, with the rerun option the
    !| whole parameter search first.
    !|
    !| Errors (`ierr`): invalid input (`ERR_INVALID_INPUT`) when some study has no sample at all, when `n_samples` is not the
    !| number of samples in `expression`, or when `sample_assignment` is not a permutation of `1..n_samples`. In the Fortran expert tier, with the rerun option, also when a capacity is below the
    !| search's own bound (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`).
    !| Errors of the search and the run are returned with their code, and `global_js_divergence` is then NaN.
    !| Known limitation: a gene missing in a whole study travels with the samples, so after relabelling that study may have
    !| values for it.
    subroutine run_js_comp_test_on_sample_assignment_impl(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                          expression, n_samples, sample_assignment, n_points, n_neighbors, &
                                                          max_n_points_candidate, max_n_neighbors_candidate, &
                                                          n_bootstrapping_top_k_jsds, global_js_divergence, used_n_points, &
                                                          used_n_neighbors, tmp_sample_slots, tmp_sample_mark, &
                                                          tmp_permuted_expression, tmp_gene_means, tmp_residuals, &
                                                          tmp_gene_means_perms, tmp_gene_means_perm_all, tmp_x_star, &
                                                          tmp_neighborhood_indices, tmp_neighborhood_range, &
                                                          tmp_neighborhood_residuals_gathered, tmp_point_pool, tmp_point_draw, &
                                                          tmp_n_bins_per_point, tmp_shared_residual_range_low, &
                                                          tmp_shared_residual_range_high, tmp_occupancy_failed, &
                                                          tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                          tmp_mean_bin_occupancy, tmp_max_bin_occupancy, tmp_sturges_bins, &
                                                          tmp_fd_bins, tmp_pmfs, tmp_counts, tmp_included_n_reps, tmp_mean_pmf, &
                                                          tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps, &
                                                          tmp_js_divergences, tmp_weights, tmp_global_js_divergence, &
                                                          tmp_pooled_residuals, tmp_pooled_residuals_perm, &
                                                          tmp_residual_to_study_perm, tmp_bin_counts_search, &
                                                          tmp_study_block_ends, tmp_study_block_sizes, tmp_point_n_neighbors, &
                                                          tmp_search_n_bins_per_point, tmp_search_shared_residual_range_low, &
                                                          tmp_search_shared_residual_range_high, &
                                                          tmp_search_confidence_interval, tmp_search_trace_n_points, &
                                                          tmp_search_trace_n_neighbors, tmp_search_trace_global_js_divergence, &
                                                          tmp_search_trace_ci_lower, tmp_search_trace_ci_upper, &
                                                          tmp_search_trace_ci_width, tmp_search_trace_ci_width_relative, &
                                                          tmp_search_trace_ci_overlap, tmp_search_trace_delta, &
                                                          tmp_search_trace_delta_median, tmp_search_trace_delta_max, &
                                                          tmp_search_trace_selected_n_bins, tmp_search_trace_occupancy_failed, &
                                                          tmp_search_trace_n_pooled_residuals, &
                                                          tmp_search_trace_min_bin_occupancy, &
                                                          tmp_search_trace_mean_bin_occupancy, &
                                                          tmp_search_trace_max_bin_occupancy, tmp_search_trace_sturges_bins, &
                                                          tmp_search_trace_fd_bins, tmp_search_trace_shared_residual_range_low, &
                                                          tmp_search_trace_shared_residual_range_high, &
                                                          tmp_candidate_confidence_interval, tmp_bootstrapping_top_k_jsds, &
                                                          tmp_prev_global_js_divergence, &
                                                          tmp_best_uncertainty_confidence_interval, tmp_best_n_bins_per_point, &
                                                          tmp_best_uncertainty_n_bins_per_point, &
                                                          tmp_best_shared_residual_range_low, &
                                                          tmp_best_shared_residual_range_high, &
                                                          tmp_best_uncertainty_shared_residual_range_low, &
                                                          tmp_best_uncertainty_shared_residual_range_high, &
                                                          rerun_parameter_search, n_bootstraps, join_method, &
                                                          min_residuals_per_bin, min_neighbor_overlap, succeeding_ci_overlap, &
                                                          plateau_mode, delta_median_threshold, delta_max_threshold, &
                                                          delta_epsilon, delta_min_consecutive_transitions, m_min, m_max, &
                                                          gamma_occupancy, lower_residual_range_quantile, &
                                                          upper_residual_range_quantile, &
                                                          two_sided_bootstrapping_significance_level, random_seed, &
                                                          bin_count_mode, &
                                                          neighborhood_size_mode, &
                                                          delta_abs_tolerance, divergence_mode, grid_points_factor, grid_min_points, ierr)
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
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(n_points)
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(n_neighbors)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(1_int32)
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
            !! Weighted global JSD of each study, the one `divergence_mode` selects (by default for exactly
            !! two studies their pairwise JSD, otherwise against the consensus), for this assignment; NaN on
            !! error
        integer(int32), intent(out) :: used_n_points
            !! Number of reference points the run used (`n_points`, or the search's choice)
        integer(int32), intent(out) :: used_n_neighbors
            !! Number of neighbors per neighborhood the run used (`n_neighbors`, or the search's choice)
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
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array of the search and the run: each study's sorting permutation of its gene means
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search and the run: sorting permutation of the pooled gene means
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array of the search and the run: the reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
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
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array of the search and the run: per-study neighbor counts, then sort scratch
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_neighbors
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
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
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
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood size the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
            !! default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Large neighborhoods, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods, `KX_FACTORS(2)` = 0.5 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! DM_DEFAULT(M_MODE_NEIGHBORHOOD_SIZE_LARGE)
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
        real(real64), intent(in), optional :: grid_points_factor
            !! Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
            !! grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_GRID_POINTS_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: grid_min_points
            !! Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
            !! floor, its last and coarsest `n_points` candidate
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_POINTS)
            !! DM_DEFAULT(M_GRID_MIN_POINTS_DEFAULT)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples_found, actual_n_bootstraps, actual_join_method
        logical(c_bool) :: actual_rerun

        call set_ok(ierr)
        M_DEFAULT_VAL(rerun_parameter_search, actual_rerun, M_RERUN_PARAMETER_SEARCH_DEFAULT)
        M_DEFAULT_VAL(n_bootstraps, actual_n_bootstraps, M_PERMUTATION_N_BOOTSTRAPS_DEFAULT)
        M_DEFAULT_VAL(join_method, actual_join_method, M_METHOD_JOIN_MIN)
        global_js_divergence = M_NAN
        used_n_points = n_points
        used_n_neighbors = n_neighbors

        ! Runtime checks (codegen guide Sec 5.14), each naming the argument at fault.
        call check_permutation_input(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                     actual_rerun, max_n_points_candidate, 9_int32, tmp_sample_slots, n_samples_found, &
                                     n_samples=n_samples, sample_assignment=sample_assignment, &
                                     sample_mark=tmp_sample_mark, max_n_neighbors_candidate=max_n_neighbors_candidate, &
                                     neighbors_arg_pos=10_int32, grid_points_factor=grid_points_factor, &
                                     grid_min_points=grid_min_points, ierr=ierr)
        if (is_err(ierr)) return

        call run_js_comp_test_on_sample_assignment_core(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=sample_assignment, n_points=n_points, &
            n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
            max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, global_js_divergence=global_js_divergence, &
            used_n_points=used_n_points, used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, tmp_point_pool=tmp_point_pool, &
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
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
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
            actual_rerun=actual_rerun, actual_n_bootstraps=actual_n_bootstraps, actual_join_method=actual_join_method, &
            min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
            succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
            delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
            delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
            m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, &
            two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
            random_seed=random_seed, bin_count_mode=bin_count_mode, neighborhood_size_mode=neighborhood_size_mode, &
            delta_abs_tolerance=delta_abs_tolerance, divergence_mode=divergence_mode, &
            grid_points_factor=grid_points_factor, grid_min_points=grid_min_points, ierr=ierr)
    end subroutine run_js_comp_test_on_sample_assignment_impl

    !> summary: The compute part of `run_js_comp_test_on_sample_assignment_impl`, without its input checks
    !| AUTHOR_LASZLO_LANG
    !| Private. Assumes what the checks of the published routine ensure: `tmp_sample_slots` holds the sample slots found by
    !| `check_permutation_input`, `sample_assignment` is a permutation of `1..n_samples`, and the capacities fit
    !| the rerun search. The published routine calls it once, the whole-test routine for the observed data and every
    !| permutation, so the checks run once per test. The arguments are those of the published routine, except the resolved
    !| `actual_rerun`, `actual_n_bootstraps` and `actual_join_method` (defaults applied).
    subroutine run_js_comp_test_on_sample_assignment_core(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, &
                                                          expression, n_samples, sample_assignment, n_points, &
                                                          n_neighbors, max_n_points_candidate, &
                                                          max_n_neighbors_candidate, n_bootstrapping_top_k_jsds, &
                                                          global_js_divergence, used_n_points, used_n_neighbors, &
                                                          tmp_sample_slots, tmp_permuted_expression, tmp_gene_means, &
                                                          tmp_residuals, tmp_gene_means_perms, &
                                                          tmp_gene_means_perm_all, tmp_x_star, &
                                                          tmp_neighborhood_indices, tmp_neighborhood_range, &
                                                          tmp_neighborhood_residuals_gathered, tmp_point_pool, &
                                                          tmp_point_draw, tmp_n_bins_per_point, &
                                                          tmp_shared_residual_range_low, &
                                                          tmp_shared_residual_range_high, tmp_occupancy_failed, &
                                                          tmp_n_pooled_residuals, tmp_min_bin_occupancy, &
                                                          tmp_mean_bin_occupancy, tmp_max_bin_occupancy, &
                                                          tmp_sturges_bins, tmp_fd_bins, tmp_pmfs, tmp_counts, &
                                                          tmp_included_n_reps, tmp_mean_pmf, tmp_mean_pmf_counts, &
                                                          tmp_mean_pmf_included_n_reps, tmp_js_divergences, &
                                                          tmp_weights, tmp_global_js_divergence, tmp_pooled_residuals, &
                                                          tmp_pooled_residuals_perm, tmp_residual_to_study_perm, &
                                                          tmp_bin_counts_search, tmp_study_block_ends, &
                                                          tmp_study_block_sizes, tmp_point_n_neighbors, &
                                                          tmp_search_n_bins_per_point, &
                                                          tmp_search_shared_residual_range_low, &
                                                          tmp_search_shared_residual_range_high, &
                                                          tmp_search_confidence_interval, tmp_search_trace_n_points, &
                                                          tmp_search_trace_n_neighbors, &
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
                                                          tmp_candidate_confidence_interval, &
                                                          tmp_bootstrapping_top_k_jsds, tmp_prev_global_js_divergence, &
                                                          tmp_best_uncertainty_confidence_interval, &
                                                          tmp_best_n_bins_per_point, &
                                                          tmp_best_uncertainty_n_bins_per_point, &
                                                          tmp_best_shared_residual_range_low, &
                                                          tmp_best_shared_residual_range_high, &
                                                          tmp_best_uncertainty_shared_residual_range_low, &
                                                          tmp_best_uncertainty_shared_residual_range_high, &
                                                          actual_rerun, actual_n_bootstraps, actual_join_method, &
                                                          min_residuals_per_bin, min_neighbor_overlap, &
                                                          succeeding_ci_overlap, plateau_mode, delta_median_threshold, &
                                                          delta_max_threshold, delta_epsilon, &
                                                          delta_min_consecutive_transitions, m_min, m_max, &
                                                          gamma_occupancy, lower_residual_range_quantile, &
                                                          upper_residual_range_quantile, &
                                                          two_sided_bootstrapping_significance_level, random_seed, &
                                                          bin_count_mode, neighborhood_size_mode, delta_abs_tolerance, &
                                                          divergence_mode, grid_points_factor, grid_min_points, ierr)
        integer(int32), intent(in) :: n_studies
        integer(int32), intent(in) :: max_n_genes_all_studies
        integer(int32), intent(in) :: max_n_reps_all_studies
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(in) :: expression
        integer(int32), intent(in) :: n_samples
        integer(int32), dimension(n_samples), intent(in) :: sample_assignment
        integer(int32), intent(in) :: n_points
        integer(int32), intent(in) :: n_neighbors
        integer(int32), intent(in) :: max_n_points_candidate
        integer(int32), intent(in) :: max_n_neighbors_candidate
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
        real(real64), dimension(n_studies), intent(out) :: global_js_divergence
        integer(int32), intent(out) :: used_n_points
        integer(int32), intent(out) :: used_n_neighbors
        integer(int32), dimension(max_n_reps_all_studies*n_studies), intent(in) :: tmp_sample_slots
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: &
            tmp_permuted_expression
        real(real64), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means
        real(real64), dimension(max_n_reps_all_studies, max_n_genes_all_studies, n_studies), intent(out) :: tmp_residuals
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
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
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_neighbors
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
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
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
        integer(int32), intent(in), optional :: bin_count_mode
        integer(int32), intent(in), optional :: neighborhood_size_mode
        real(real64), intent(in), optional :: delta_abs_tolerance
        integer(int32), intent(in), optional :: divergence_mode
        real(real64), intent(in), optional :: grid_points_factor
        integer(int32), intent(in), optional :: grid_min_points
        integer(int32), intent(out) :: ierr

        integer(int32) :: n_pool, callee_ierr, max_n_bins_per_point, search_n_admissible_evaluated
        logical(c_bool) :: search_plateau_established

        call set_ok(ierr)
        global_js_divergence = M_NAN
        used_n_points = n_points
        used_n_neighbors = n_neighbors

        call assign_samples_and_center(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, n_samples, &
                                       sample_assignment, tmp_sample_slots, tmp_permuted_expression, tmp_gene_means, &
                                       tmp_residuals)

        ! Every capacity-sized buffer is passed whole to explicit-shape dummies declared at the
        ! sizes actually used (sequence association, no copy), as the fixed-k search does.
        if (actual_rerun) then
            call run_js_comp_test_parameter_search_impl( &
                n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
                max_n_reps_all_studies=max_n_reps_all_studies, gene_means=tmp_gene_means, residuals=tmp_residuals, &
                n_bootstraps=actual_n_bootstraps, join_method=actual_join_method, &
                max_n_points_candidate=max_n_points_candidate, max_n_neighbors_candidate=max_n_neighbors_candidate, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, n_points=used_n_points, &
                n_neighbors=used_n_neighbors, n_bins_per_point=tmp_search_n_bins_per_point, &
                shared_residual_range_low=tmp_search_shared_residual_range_low, &
                shared_residual_range_high=tmp_search_shared_residual_range_high, &
                best_candidate_pair_confidence_interval=tmp_search_confidence_interval, &
                plateau_established=search_plateau_established, n_admissible_evaluated=search_n_admissible_evaluated, &
                trace_n_points=tmp_search_trace_n_points, trace_n_neighbors=tmp_search_trace_n_neighbors, &
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
                tmp_gene_means_perms=tmp_gene_means_perms, tmp_gene_means_perm_all=tmp_gene_means_perm_all, &
                tmp_x_star=tmp_x_star, tmp_neighborhood_indices_all_studies=tmp_neighborhood_indices, &
                tmp_neighborhood_range=tmp_neighborhood_range, &
                tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, tmp_point_pool=tmp_point_pool, &
                tmp_point_draw=tmp_point_draw, tmp_n_bins_per_point=tmp_n_bins_per_point, &
                tmp_shared_residual_range_low=tmp_shared_residual_range_low, &
                tmp_shared_residual_range_high=tmp_shared_residual_range_high, tmp_pmfs=tmp_pmfs, tmp_counts=tmp_counts, &
                tmp_included_n_reps=tmp_included_n_reps, tmp_mean_pmf=tmp_mean_pmf, &
                tmp_mean_pmf_counts=tmp_mean_pmf_counts, tmp_mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
                tmp_js_divergences=tmp_js_divergences, tmp_weights=tmp_weights, &
                tmp_global_js_divergence=tmp_global_js_divergence, &
                tmp_confidence_interval=tmp_candidate_confidence_interval, &
                tmp_bootstrapping_top_k_jsds=tmp_bootstrapping_top_k_jsds, &
                tmp_prev_global_js_divergence=tmp_prev_global_js_divergence, tmp_delta_perm=tmp_point_n_neighbors, &
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
                random_seed=random_seed, bin_count_mode=bin_count_mode, &
                neighborhood_size_mode=neighborhood_size_mode, &
                delta_abs_tolerance=delta_abs_tolerance, divergence_mode=divergence_mode, &
                grid_points_factor=grid_points_factor, grid_min_points=grid_min_points, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                call set_err_once(ierr, get_err_code(callee_ierr))
                return
            end if
        end if

        ! The reference points as the fixed-k search computes them, each study's order as the plain
        ! run_js_comp_test's prologue builds it.
        call sort_pooled_gene_means(tmp_gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perm_all)
        call pool_means_impl(tmp_gene_means, tmp_gene_means_perm_all, max_n_genes_all_studies*n_studies, used_n_points, &
                             n_pool, tmp_x_star)
        call sort_gene_means_per_study(tmp_gene_means, max_n_genes_all_studies, n_studies, tmp_gene_means_perms)

        call run_js_comp_test_impl( &
            n_studies=n_studies, max_n_genes_all_studies=max_n_genes_all_studies, &
            max_n_reps_all_studies=max_n_reps_all_studies, n_points=used_n_points, n_neighbors=used_n_neighbors, &
            gene_means=tmp_gene_means, gene_means_perms=tmp_gene_means_perms, residuals=tmp_residuals, x_star=tmp_x_star, &
            neighborhood_indices=tmp_neighborhood_indices, neighborhood_range=tmp_neighborhood_range, &
            n_bins_per_point=tmp_n_bins_per_point, shared_residual_range_low=tmp_shared_residual_range_low, &
            shared_residual_range_high=tmp_shared_residual_range_high, max_n_bins_per_point=max_n_bins_per_point, &
            occupancy_failed=tmp_occupancy_failed, n_pooled_residuals=tmp_n_pooled_residuals, &
            min_bin_occupancy=tmp_min_bin_occupancy, mean_bin_occupancy=tmp_mean_bin_occupancy, &
            max_bin_occupancy=tmp_max_bin_occupancy, sturges_bins=tmp_sturges_bins, fd_bins=tmp_fd_bins, pmfs=tmp_pmfs, &
            counts=tmp_counts, included_n_reps=tmp_included_n_reps, mean_pmf=tmp_mean_pmf, &
            mean_pmf_counts=tmp_mean_pmf_counts, mean_pmf_included_n_reps=tmp_mean_pmf_included_n_reps, &
            js_divergences=tmp_js_divergences, weights=tmp_weights, global_js_divergence=global_js_divergence, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, &
            tmp_pooled_residuals=tmp_pooled_residuals, &
            tmp_pooled_residuals_perm=tmp_pooled_residuals_perm, tmp_bin_counts_search=tmp_bin_counts_search, &
            tmp_residual_to_study_perm=tmp_residual_to_study_perm, tmp_study_block_ends=tmp_study_block_ends, &
            tmp_study_block_sizes=tmp_study_block_sizes, tmp_point_n_neighbors=tmp_point_n_neighbors, &
            min_residuals_per_bin=min_residuals_per_bin, m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, bin_count_mode=bin_count_mode, &
            divergence_mode=divergence_mode, ierr=callee_ierr)
        if (is_err(callee_ierr)) then
            call set_err_once(ierr, get_err_code(callee_ierr))
            global_js_divergence = M_NAN
            return
        end if
    end subroutine run_js_comp_test_on_sample_assignment_core

    !> summary: Permutation-Test of the k-NN JS-Comp-Test that permutes the samples' study labels
    !| AUTHOR_LASZLO_LANG
    !| Tests whether the assignment of samples to studies influences each study's JS-Comp-Test value (null
    !| hypothesis: it does not), with the data, parameters and options of
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]].
    !|
    !| 1. The observed data (the identity assignment) give `global_js_divergence_observed`.
    !| 2. `n_permutations` times, the study labels are shuffled without replacement (a uniformly random permutation of
    !|    `1..n_samples`, GSL's `gsl_ran_shuffle`, one generator seeded with `random_seed`), and the whole JS-Comp-Test runs on the
    !|    shuffled data. Row `i` of `jsd_null` holds its values.
    !| 3. Study `i`'s p-value is `(#{valid permutations with J_perm >= J_obs} + 1) / (n_valid + 1)`, so it is never 0.
    !|
    !| With `rerun_parameter_search = .true.` every data set, the observed one included, first runs the fixed-k parameter search. `random_seed` then
    !| also seeds the bootstrap of every search: each search starts from the same seed, so every data set gets the same bootstrap
    !| draws. The routine is sequential: a caller who wants parallel permutations calls
    !| [[tox_data_integration_js_comp_test_permutation_fixed_k_impl(module):run_js_comp_test_on_sample_assignment_impl(interface)]]
    !| itself, once per assignment.
    !|
    !| Failures: a permutation that fails with invalid input (a safeguard: the relabelling keeps every study's sample count, so the k-NN test does not normally trigger it) gets a NaN row in `jsd_null`, is
    !| counted in `n_failed_permutations` and is not in `n_valid`. No valid permutation gives NaN p-values (also for
    !| `n_permutations = 0` in a direct call; the published entry points reject 0). Any other error ends the test and is
    !| returned, with NaN p-values and the rows run so far kept. An error on the observed data, or a study without any sample,
    !| is returned at once with everything NaN. In the Fortran expert tier a capacity below the search's own bound
    !| (`calc_js_comp_test_candidate_bounds_impl` for the same `grid_points_factor` and `grid_min_points`) fails the same way with the rerun option.
    !| Known limitation: a gene missing in a whole study travels with the samples.
    subroutine run_js_comp_test_permutation_impl(n_studies, max_n_genes_all_studies, &
                                                 max_n_reps_all_studies, expression, n_points, &
                                                 n_neighbors, max_n_points_candidate, &
                                                 max_n_neighbors_candidate, n_bootstrapping_top_k_jsds, &
                                                 n_permutations, global_js_divergence_observed, p_values, &
                                                 jsd_null, n_failed_permutations, tmp_sample_slots, &
                                                 tmp_sample_assignment, tmp_permuted_global_js_divergence, &
                                                 tmp_permuted_expression, tmp_gene_means, &
                                                 tmp_residuals, tmp_gene_means_perms, &
                                                 tmp_gene_means_perm_all, tmp_x_star, &
                                                 tmp_neighborhood_indices, tmp_neighborhood_range, &
                                                 tmp_neighborhood_residuals_gathered, tmp_point_pool, &
                                                 tmp_point_draw, tmp_n_bins_per_point, &
                                                 tmp_shared_residual_range_low, &
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
                                                 tmp_point_n_neighbors, tmp_search_n_bins_per_point, &
                                                 tmp_search_shared_residual_range_low, &
                                                 tmp_search_shared_residual_range_high, &
                                                 tmp_search_confidence_interval, &
                                                 tmp_search_trace_n_points, tmp_search_trace_n_neighbors, &
                                                 tmp_search_trace_global_js_divergence, &
                                                 tmp_search_trace_ci_lower, tmp_search_trace_ci_upper, &
                                                 tmp_search_trace_ci_width, &
                                                 tmp_search_trace_ci_width_relative, &
                                                 tmp_search_trace_ci_overlap, tmp_search_trace_delta, &
                                                 tmp_search_trace_delta_median, &
                                                 tmp_search_trace_delta_max, &
                                                 tmp_search_trace_selected_n_bins, &
                                                 tmp_search_trace_occupancy_failed, &
                                                 tmp_search_trace_n_pooled_residuals, &
                                                 tmp_search_trace_min_bin_occupancy, &
                                                 tmp_search_trace_mean_bin_occupancy, &
                                                 tmp_search_trace_max_bin_occupancy, &
                                                 tmp_search_trace_sturges_bins, tmp_search_trace_fd_bins, &
                                                 tmp_search_trace_shared_residual_range_low, &
                                                 tmp_search_trace_shared_residual_range_high, &
                                                 tmp_candidate_confidence_interval, &
                                                 tmp_bootstrapping_top_k_jsds, &
                                                 tmp_prev_global_js_divergence, &
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
                                                 bin_count_mode, &
                                                 neighborhood_size_mode, &
                                                 delta_abs_tolerance, divergence_mode, grid_points_factor, grid_min_points, ierr)
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
        integer(int32), intent(in) :: n_points
            !! Number of reference points of the run (sizes the work arrays only, with the rerun option)
            !! DM_MIN(1_int32)
        integer(int32), intent(in) :: n_neighbors
            !! Number of neighbors per neighborhood of the run (sizes the work arrays only, with the
            !! rerun option)
            !! DM_MIN(1_int32)
            !! DM_MAX(max_n_genes_all_studies)
        integer(int32), intent(in) :: max_n_points_candidate
            !! Capacity of the per-point work arrays
            !! DM_OUTPUT_FROM(max_n_points_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(n_points)
        integer(int32), intent(in) :: max_n_neighbors_candidate
            !! Capacity of the per-neighbor work arrays
            !! DM_OUTPUT_FROM(max_n_neighbors_candidate, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
            !! DM_MIN(n_neighbors)
        integer(int32), intent(in) :: n_bootstrapping_top_k_jsds
            !! Bootstrap heap size of the rerun parameter search
            !! DM_OUTPUT_FROM(n_bootstrapping_top_k_jsds, calc_js_comp_test_permutation_bounds, tox_data_integration_js_comp_test_permutation_bounds, AUTO)
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
        integer(int32), dimension(max_n_genes_all_studies, n_studies), intent(out) :: tmp_gene_means_perms
            !! Working array of the search and the run: each study's sorting permutation of its gene means
        integer(int32), dimension(max_n_genes_all_studies*n_studies), intent(out) :: tmp_gene_means_perm_all
            !! Working array of the search and the run: sorting permutation of the pooled gene means
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_x_star
            !! Working array of the search and the run: the reference points
        integer(int32), dimension(max_n_neighbors_candidate, max_n_points_candidate, n_studies), intent(out) :: &
            tmp_neighborhood_indices
            !! Working array of the search and the run: every study's neighbor gene indices
        integer(int32), dimension(2, max_n_points_candidate, n_studies), intent(out) :: tmp_neighborhood_range
            !! Working array of the search and the run: neighborhood spans
        real(real64), dimension(max_n_reps_all_studies, max_n_neighbors_candidate, max_n_points_candidate), &
            intent(out) :: tmp_neighborhood_residuals_gathered
            !! Working array of the search and the run: one study's gathered neighborhood residuals
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
        real(real64), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals
            !! Working array of the search and the run: one reference point's pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_pooled_residuals_perm
            !! Working array of the search and the run: sorting permutation of the pooled residuals
        integer(int32), dimension(max_n_reps_all_studies*max_n_neighbors_candidate*n_studies), intent(out) :: &
            tmp_residual_to_study_perm
            !! Working array of the search and the run (Stone modes): pool indices grouped by study
        integer(int32), dimension(M_MAX_N_BINS), intent(out) :: tmp_bin_counts_search
            !! Working array of the search and the run: bin counts of the bin-count search
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_ends
            !! Working array of the search and the run: the last index of each study's pooled block
        integer(int32), dimension(n_studies), intent(out) :: tmp_study_block_sizes
            !! Working array of the search and the run (Stone modes): each study's pooled block size
        integer(int32), dimension(n_studies), intent(out) :: tmp_point_n_neighbors
            !! Working array of the search and the run: per-study neighbor counts, then sort scratch
        integer(int32), dimension(max_n_points_candidate), intent(out) :: tmp_search_n_bins_per_point
            !! Working array: the search's chosen bin counts
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_low
            !! Working array: the search's chosen lower residual-range bounds
        real(real64), dimension(max_n_points_candidate), intent(out) :: tmp_search_shared_residual_range_high
            !! Working array: the search's chosen upper residual-range bounds
        real(real64), dimension(2, n_studies), intent(out) :: tmp_search_confidence_interval
            !! Working array: the search's chosen confidence interval
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_points
            !! Working array: the search's trace
        integer(int32), dimension(M_MAX_CANDIDATE_PAIRS), intent(out) :: tmp_search_trace_n_neighbors
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
        real(real64), dimension(2, n_studies), intent(out) :: tmp_candidate_confidence_interval
            !! Working array of the search: the current candidate's confidence interval
        real(real64), dimension(n_bootstrapping_top_k_jsds, 2, n_studies), intent(out) :: tmp_bootstrapping_top_k_jsds
            !! Working array of the search: the bootstrap's bottom-k/top-k heaps
        real(real64), dimension(n_studies), intent(out) :: tmp_prev_global_js_divergence
            !! Working array of the search: the previous admissible candidate's global JSD
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
            !! `.true.`: run the fixed-k parameter search on the relabelled data first and use its
            !! `(n_points, n_neighbors)`; `.false.`: use `n_points` and `n_neighbors`
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
        integer(int32), intent(in), optional :: neighborhood_size_mode
            !! Rerun search: which neighborhood size the k-NN candidate grid pairs with each
            !! `n_points` (used only with `rerun_parameter_search`). Small neighborhoods are not the
            !! default: a higher JSCompTest value from small neighborhoods only shows that they are
            !! too small
            !!
            !! | Mode | Value |
            !! |------|-------|
            !! | Large neighborhoods, `KX_FACTORS(1)` = 0.25 (default) | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_LARGE(variable)]] |
            !! | Small neighborhoods, `KX_FACTORS(2)` = 0.5 | [[tox_data_integration_js_comp_test_fixed_k_search_impl(module):MODE_NEIGHBORHOOD_SIZE_SMALL(variable)]] |
            !! DM_DEFAULT(M_MODE_NEIGHBORHOOD_SIZE_LARGE)
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
        real(real64), intent(in), optional :: grid_points_factor
            !! Rerun search (used only with `rerun_parameter_search`): factor in the k-NN candidate
            !! grid's first reference-point count `ceiling(grid_points_factor * sqrt(max_n_genes_all_studies))`
            !! DM_MIN(above(0.0_real64))
            !! DM_DEFAULT(M_GRID_POINTS_FACTOR_DEFAULT)
        integer(int32), intent(in), optional :: grid_min_points
            !! Rerun search (used only with `rerun_parameter_search`): the k-NN candidate grid's
            !! floor, its last and coarsest `n_points` candidate
            !! DM_MIN(1_int32)
            !! DM_MAX(MAX_POINTS)
            !! DM_DEFAULT(M_GRID_MIN_POINTS_DEFAULT)
        integer(int32), intent(out) :: ierr
            !! Error code; ERR_INVALID_INPUT, or the code of the search or the run (see above)

        integer(int32) :: n_samples, i_permutation, actual_random_seed, callee_ierr, actual_n_bootstraps, actual_join_method
        integer(int32) :: used_n_points, used_n_neighbors
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

        ! The checks that do not depend on the assignment run once here, with this routine's argument
        ! positions. Every data set below runs the core, which does not repeat them.
        call check_permutation_input(n_studies, max_n_genes_all_studies, max_n_reps_all_studies, expression, &
                                     actual_rerun, max_n_points_candidate, 7_int32, tmp_sample_slots, n_samples, &
                                     max_n_neighbors_candidate=max_n_neighbors_candidate, neighbors_arg_pos=8_int32, &
                                     grid_points_factor=grid_points_factor, grid_min_points=grid_min_points, ierr=ierr)
        if (is_err(ierr)) return

        ! The observed data: the identity assignment. An error here ends the test at once.
        call init_perm(tmp_sample_assignment(1:n_samples))
        call run_js_comp_test_on_sample_assignment_core(n_studies=n_studies, &
            max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
            expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, n_points=n_points, &
            n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
            max_n_neighbors_candidate=max_n_neighbors_candidate, &
            n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, global_js_divergence=global_js_divergence_observed, &
            used_n_points=used_n_points, used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, &
            tmp_permuted_expression=tmp_permuted_expression, tmp_gene_means=tmp_gene_means, &
            tmp_residuals=tmp_residuals, tmp_gene_means_perms=tmp_gene_means_perms, &
            tmp_gene_means_perm_all=tmp_gene_means_perm_all, tmp_x_star=tmp_x_star, &
            tmp_neighborhood_indices=tmp_neighborhood_indices, tmp_neighborhood_range=tmp_neighborhood_range, &
            tmp_neighborhood_residuals_gathered=tmp_neighborhood_residuals_gathered, tmp_point_pool=tmp_point_pool, &
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
            tmp_point_n_neighbors=tmp_point_n_neighbors, tmp_search_n_bins_per_point=tmp_search_n_bins_per_point, &
            tmp_search_shared_residual_range_low=tmp_search_shared_residual_range_low, &
            tmp_search_shared_residual_range_high=tmp_search_shared_residual_range_high, &
            tmp_search_confidence_interval=tmp_search_confidence_interval, &
            tmp_search_trace_n_points=tmp_search_trace_n_points, &
            tmp_search_trace_n_neighbors=tmp_search_trace_n_neighbors, &
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
            actual_rerun=actual_rerun, actual_n_bootstraps=actual_n_bootstraps, actual_join_method=actual_join_method, &
            min_residuals_per_bin=min_residuals_per_bin, min_neighbor_overlap=min_neighbor_overlap, &
            succeeding_ci_overlap=succeeding_ci_overlap, plateau_mode=plateau_mode, &
            delta_median_threshold=delta_median_threshold, delta_max_threshold=delta_max_threshold, &
            delta_epsilon=delta_epsilon, delta_min_consecutive_transitions=delta_min_consecutive_transitions, &
            m_min=m_min, m_max=m_max, gamma_occupancy=gamma_occupancy, &
            lower_residual_range_quantile=lower_residual_range_quantile, &
            upper_residual_range_quantile=upper_residual_range_quantile, &
            two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
            random_seed=actual_random_seed, bin_count_mode=bin_count_mode, &
            neighborhood_size_mode=neighborhood_size_mode, delta_abs_tolerance=delta_abs_tolerance, &
            divergence_mode=divergence_mode, grid_points_factor=grid_points_factor, grid_min_points=grid_min_points, &
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
            call run_js_comp_test_on_sample_assignment_core(n_studies=n_studies, &
                max_n_genes_all_studies=max_n_genes_all_studies, max_n_reps_all_studies=max_n_reps_all_studies, &
                expression=expression, n_samples=n_samples, sample_assignment=tmp_sample_assignment, &
                n_points=n_points, n_neighbors=n_neighbors, max_n_points_candidate=max_n_points_candidate, &
                max_n_neighbors_candidate=max_n_neighbors_candidate, &
                n_bootstrapping_top_k_jsds=n_bootstrapping_top_k_jsds, &
                global_js_divergence=tmp_permuted_global_js_divergence, used_n_points=used_n_points, &
                used_n_neighbors=used_n_neighbors, tmp_sample_slots=tmp_sample_slots, &
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
                actual_rerun=actual_rerun, actual_n_bootstraps=actual_n_bootstraps, &
                actual_join_method=actual_join_method, min_residuals_per_bin=min_residuals_per_bin, &
                min_neighbor_overlap=min_neighbor_overlap, succeeding_ci_overlap=succeeding_ci_overlap, &
                plateau_mode=plateau_mode, delta_median_threshold=delta_median_threshold, &
                delta_max_threshold=delta_max_threshold, delta_epsilon=delta_epsilon, &
                delta_min_consecutive_transitions=delta_min_consecutive_transitions, m_min=m_min, m_max=m_max, &
                gamma_occupancy=gamma_occupancy, lower_residual_range_quantile=lower_residual_range_quantile, &
                upper_residual_range_quantile=upper_residual_range_quantile, &
                two_sided_bootstrapping_significance_level=two_sided_bootstrapping_significance_level, &
                random_seed=actual_random_seed, bin_count_mode=bin_count_mode, &
                neighborhood_size_mode=neighborhood_size_mode, delta_abs_tolerance=delta_abs_tolerance, &
                divergence_mode=divergence_mode, grid_points_factor=grid_points_factor, &
                grid_min_points=grid_min_points, ierr=callee_ierr)
            if (is_err(callee_ierr)) then
                ! Invalid input of a relabelled data set fails this permutation only (a safeguard, kept in step
                ! with the adaptive test); anything else ends the test.
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
    end subroutine run_js_comp_test_permutation_impl
end module tox_data_integration_js_comp_test_permutation_fixed_k_impl
