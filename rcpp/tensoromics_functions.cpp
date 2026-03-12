#include <Rcpp.h>
using namespace Rcpp;


 

// ===================================================================
// FORTRAN FUNCTIONS
// ===================================================================

extern "C" {
void compute_edf_c(
  const double* values,
  const int* n_values,
  double* unique_values,
  double* cdf_values,
  int* n_unique,
  int* ierr
);
void compute_edf_expert_c(
  const double* values,
  const int* n_values,
  const int* perm,
  double* unique_values,
  double* cdf_values,
  int* n_unique,
  int* ierr
);
void compute_all_contributions_c(
  const double* trajectories,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  const int* factor_indices,
  const int* n_selected_factors,
  const int* dependent_indices,
  const int* n_selected_dependents,
  const char* mode,
  double* local_contributions,
  double* total_contributions,
  double* temp_factors,
  double* temp_dependent,
  int* ierr
);
void compute_baselines_factor_dependent_c(
  const double* factor,
  const double* dependent,
  const int* n_timepoints,
  const char* mode,
  double* factor_baseline,
  double* dependent_baseline,
  int* ierr
);
void compute_contributions_c(
  const double* factor,
  const double* dependent,
  const int* n_dims,
  const char* mode,
  double* local_contributions,
  double* total_contribution,
  int* ierr
);
void perform_permutation_test_c(
  const double* trajectories,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  const int* factor_idx,
  const int* dependent_idx,
  const int* sample_idx,
  const char* mode,
  const int* n_permutations,
  double* local_contributions,
  double* total_contributions,
  double* temp_factor,
  double* temp_dependent,
  int* ierr,
  const int* random_seed
);
void compute_p_values_c(
  const double* local_contributions_observed,
  const double* total_contribution_observed,
  const double* local_contributions_perm,
  const double* total_contributions_perm,
  const int* n_timepoints,
  const int* n_permutations,
  double* local_p_values,
  double* total_p_value,
  int* ierr
);

void compute_velocity_trajectories_c(
  const double* trajectories,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  double* velocity,
  int* ierr
);
void compute_acceleration_from_velocity_c(
  const double* velocity,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  double* acceleration,
  int* ierr
);
void compute_velocity_trajectory_c(
  const double* trajectory,
  const int* n_timepoints,
  double* velocity,
  int* ierr
);
void compute_acceleration_from_velocity_trajectory_c(
  const double* velocity,
  const int* n_timepoints,
  double* acceleration,
  int* ierr
);
void compute_velocity_acceleration_contributions_c(
  const double* trajectories,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  const char* mode,
  double* factor_workspace,
  double* dependent_workspace,
  double* contributions_workspace,
  double* contrib_velocity,
  double* velocity_contribution_series,
  double* contrib_acceleration,
  double* acceleration_contribution_series,
  int* ierr
);
void compute_velocity_acceleration_contributions_alloc_c(
  const double* trajectories,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  const char* mode,
  double* contrib_velocity,
  double* velocity_contribution_series,
  double* contrib_acceleration,
  double* acceleration_contribution_series,
  int* ierr
);
void cluster_factor_trajectories_k_means_c(
  const int* n_clusters,
  const double* trajectories,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  double* centroids,
  int* labels,
  int* label_counts,
  const int* max_iterations,
  int* ierr
);

void k_means_clustering_c(
  const int* n_clusters,
  const double* data_points,
  const int* n_points,
  const int* n_dims,
  double* centroids,
  int* labels,
  int* label_counts,
  const int* max_iterations,
  int* ierr
);

void linkage_clustering_c(
  double* distances,
  const int* n_points,
  int* merge_i,
  int* merge_j,
  double* heights,
  int* cluster_sizes,
  const char* method,
  int* ierr
);

void relative_axes_changes_from_shift_vector_c(
  const double* vec,
  const int* n_axes,
  double* contributions,
  int* ierr
);

void relative_axes_expression_from_expression_vector_c(
  const double* vec,
  const int* n_axes,
  double* contributions,
  int* ierr
);

void clock_hand_angle_between_vectors_c(
  const double* v1,
  const double* v2,
  const int* n_dims,
  double* signed_angle,
  const int* selected_axes_for_signed,
  int* ierr
);

void clock_hand_angles_for_shift_vectors_c(
  const double* origins,
  const double* targets,
  const int* n_dims,
  const int* n_vecs,
  const int* vecs_selection_mask,
  const int* n_selected_vecs,
  const int* selected_axes_for_signed,
  double* signed_angles,
  int* ierr
);

void omics_vector_RAP_projection_c(
  const double* vecs,
  const int* n_axes,
  const int* n_vecs,
  const int* vecs_selection_mask,
  const int* n_selected_vecs,
  const int* axes_selection_mask,
  const int* n_selected_axes,
  double* projections,
  int* ierr
);

void omics_field_RAP_projection_c(
  const double* vecs,
  const int* n_axes,
  const int* n_vecs,
  const int* vecs_selection_mask,
  const int* n_selected_vecs,
  const int* axes_selection_mask,
  const int* n_selected_axes,
  double* projections,
  int* ierr
);

void normalize_by_std_dev_c(
  const int* n_genes,
  const int* n_tissues,
  const double* input_matrix,
  double* output_matrix,
  int* ierr
);

void quantile_normalization_c(
  const int* n_genes,
  const int* n_tissues,
  const double* input_matrix,
  double* output_matrix,
  double* temp_col,
  double* rank_means,
  int* perm,
  int* stack_left,
  int* stack_right,
  const int* max_stack,
  int* ierr
);

void log2_transformation_c(
  const int* n_genes,
  const int* n_tissues,
  const double* input_matrix,
  double* output_matrix,
  int* ierr
);

void normalize_unit_length_c(
  double* vector,
  const int* n_dims,
  int* ierr
);

void calc_tiss_avg_c(
  const int* n_genes,
  const int* n_grps,
  const int* group_s,
  const int* group_c,
  const double* input_matrix,
  double* output_matrix,
  int* ierr
);

void calc_fchange_c(
  const int* n_genes,
  const int* n_cols,
  const int* n_pairs,
  const int* control_cols,
  const int* cond_cols,
  const double* input_matrix,
  double* output_matrix,
  int* ierr
);

void fjct_compute_jsd_c(
  const int* family_idx,
  const int* gene_to_family_S1,
  const int* gene_to_family_S2,
  const int* n_genes_S1,
  const int* n_genes_S2,
  const double* neighborhood_residuals_S1,
  const double* neighborhood_residuals_S2,
  const int* neighborhood_genes_S1,
  const int* neighborhood_genes_S2,
  const int* n_reps_S1,
  const int* n_reps_S2,
  const int* n_neighbors,
  const int* n_points,
  const int* n_bins,
  const double* shared_residual_range,
  double* js_divergences,
  int* included_n_reps_S1,
  int* included_n_reps_S2,
  int* total_included_n_reps,
  double* global_js_divergence,
  double* weights,
  int* ierr
);

void fjct_compute_jsd_expert_c(
  const double* neighborhood_residuals_S1,
  const double* neighborhood_residuals_S2,
  const int* n_reps_S1,
  const int* n_reps_S2,
  const int* n_neighbors,
  const int* n_points,
  const int* neighbor_mask_S1,
  const int* neighbor_mask_S2,
  const int* n_bins,
  const double* shared_residual_range,
  double* js_divergences,
  int* included_n_reps_S1,
  int* included_n_reps_S2,
  int* total_included_n_reps,
  double* global_js_divergence,
  double* weights,
  double* pmf_S1,
  double* pmf_S2,
  int* tmp_counts,
  int* ierr
);

void fjct_compute_contribution_scores_c(
  const double* global_js_divergences,
  const int* total_included_n_reps_per_f,
  const int* k_families,
  const double* support_weights,
  double* contribution_scores,
  int* ierr
);

void gjct_permutation_test_c(
  const double* neighborhood_residuals_S1,
  const double* neighborhood_residuals_S2,
  const int* n_reps_S1,
  const int* n_reps_S2,
  const int* n_neighbors,
  const int* n_points,
  const double* global_jsd_observed,
  const int* n_bins,
  const double* shared_residual_range,
  const int* n_permutations,
  double* jsd_null,
  double* p_value,
  int* ierr,
  const int* random_seed
);

void gjct_permutation_test_filtered_c(
  const double* neighborhood_residuals_S1,
  const double* neighborhood_residuals_S2,
  const int* n_reps_S1,
  const int* n_reps_S2,
  const int* n_neighbors,
  const int* n_points,
  const double* global_jsd_observed,
  const int* n_bins,
  const double* shared_residual_range,
  const int* n_permutations,
  double* jsd_null,
  double* p_value,
  int* ierr,
  const int* random_seed,
  const int* neighbor_mask_S1,
  const int* neighbor_mask_S2
);

void compute_weighted_global_divergence_c(
  const double* js_divergences,
  const int* n_points,
  const int* included_n_residuals_S1,
  const int* included_n_residuals_S2,
  double* global_js_divergence,
  double* weights,
  int* ierr
);

void compute_divergence_per_reference_point_c(
  const double* pmf_S1,
  const double* pmf_S2,
  const int* n_points,
  const int* n_bins,
  double* js_divergences,
  int* ierr
);

void build_residual_histograms_c(
  const double* neighborhood_residuals,
  const int* n_reps,
  const int* n_neighbors,
  const int* n_points,
  const double* shared_residual_range,
  const int* n_bins,
  int* counts,
  double* pmf,
  int* included_n_residuals,
  int* ierr
);

void build_residual_histograms_filtered_c(
  const double* neighborhood_residuals,
  const int* n_reps,
  const int* n_neighbors,
  const int* n_points,
  const double* shared_residual_range,
  const int* n_bins,
  int* counts,
  double* pmf,
  int* included_n_residuals,
  int* ierr,
  const int* neighbor_mask
);

void determine_shared_residual_range_c(
  const double* neighborhood_residuals_S1,
  const double* neighborhood_residuals_S2,
  const int* n_reps_S1,
  const int* n_reps_S2,
  const int* n_neighbors,
  const int* n_points,
  const double* residual_range_quantile,
  double* shared_residual_range,
  int* ierr
);

void determine_shared_residual_range_expert_c(
  const double* residual_pool,
  int* residual_pool_perm,
  const int* n_pool,
  const double* residual_range_quantile,
  double* shared_R,
  int* ierr
);

void compute_gene_means_c(
  const int* n_genes,
  const int* n_reps,
  const double* expr,
  double* means,
  int* ierr
);

void compute_residuals_c(
  const int* n_genes,
  const int* n_reps,
  const double* expr,
  const double* means,
  double* resid,
  int* ierr
);

void pool_means_c(
  const int* n_genes_S1,
  const double* mean_S1,
  const int* n_genes_S2,
  const double* mean_S2,
  const int* n_points,
  const int* n_pool,
  double* x_star,
  int* ierr
);

void pool_means_expert_c(
  const double* pooled_means,
  int* pooled_means_perm,
  const int* pool_size,
  const int* n_points,
  const int* n_pool,
  double* x_star,
  int* ierr
);

void calc_neighborhood_size_c(
  const int* n_pool,
  const int* n_points,
  const int* n_genes_S,
  const double* mean_S,
  const int* desired_size,
  int* n_neighbors,
  int* ierr
);

void construct_neighborhoods_c(
  const int* n_points,
  const double* x_star,
  const int* n_genes_S,
  const double* mean_S,
  const int* n_reps_S,
  const double* resid_S,
  double* neighborhood_residuals,
  int* neighborhood_indices,
  const int* n_neighbors,
  int* ierr
);

    

void normalization_pipeline_c(
  const int* n_genes,
  const int* n_tissues,
  const double* input_matrix,
  double* buf_stddev,
  double* buf_quant,
  double* buf_avg,
  double* buf_log,
  double* temp_col,
  double* rank_means,
  int* perm,
  int* stack_left,
  int* stack_right,
  const int* max_stack,
  const int* group_s,
  const int* group_c,
  const int* n_grps,
  int* ierr
);

void normalize_variable_timeseries_C(
  const double* v,
  double* v_norm,
  const int* n_points,
  int* ierr,
  int* status
);

void normalize_single_trajectory_C(
  const double* trajectory,
  double* trajectory_norm,
  const int* n_factors,
  const int* n_timepoints,
  int* ierr,
  int* status
);

void normalize_all_trajectories_C(
  const double* trajectories,
  double* trajectories_norm,
  const int* n_factors,
  const int* n_samples,
  const int* n_timepoints,
  int* ierr,
  int* status
);

void compute_family_scaling_c(
  const int* n_genes,
  const int* n_families,
  const double* distances,
  const int* gene_to_fam,
  double* dscale,
  double* loess_x,
  double* loess_y,
  int* indices_used,
  int* ierr
);

void compute_family_scaling_expert_c(
  const int* n_genes,
  const int* n_families,
  const double* distances,
  const int* gene_to_fam,
  double* dscale,
  double* loess_x,
  double* loess_y,
  int* indices_used,
  int* perm_tmp,
  int* stack_left_tmp,
  int* stack_right_tmp,
  double* family_distances,
  int* ierr
);

void compute_rdi_c(
  const int* n_genes,
  const int* n_families,
  const double* distances,
  const int* gene_to_fam,
  const double* dscale,
  double* rdi,
  double* sorted_rdi,
  int* perm,
  int* stack_left,
  int* stack_right
);

void identify_outliers_c(
  const int* n_genes,
  const double* rdi,
  const double* sorted_rdi,
  int* is_outlier_int,
  double* threshold,
  const double* percentile
);

void detect_outliers_c(
  const int* n_genes,
  const int* n_families,
  const double* distances,
  const int* gene_to_fam,
  double* work_array,
  int* perm,
  int* stack_left,
  int* stack_right,
  int* is_outlier_int,
  double* loess_x,
  double* loess_y,
  int* loess_n,
  const double* percentile,
  int* ierr
);

void euclidean_distance_c(
  const double* vec1,
  const double* vec2,
  const int* d,
  double* result
);

void distance_to_centroid_c(
  const int* n_genes,
  const int* n_families,
  const double* genes,
  const double* centroids,
  const int* gene_to_fam,
  double* distances,
  const int* d
);

void compute_tissue_versatility_c(
  const int* n_axes,
  const int* n_vectors,
  const double* expression_vectors,
  const int* exp_vecs_selection_index,
  const int* n_selected_vectors,
  const int* axes_selection,
  const int* n_selected_axes,
  double* tissue_versatilities,
  double* tissue_angles_deg,
  int* ierr
);

void compute_shift_vector_field_c(
  const int* d,
  const int* n_genes,
  const int* n_families,
  const double* expression_vectors,
  const double* family_centroids,
  const int* gene_to_centroid,
  double* shift_vectors,
  int* ierr
);


void detect_neofunctionalization_c(
  const double* ancestors,
  const int* n_families,
  const double* genes,
  const int* n_axes,
  const int* gene_to_fam,
  const int* n_genes,
  const double* thresholds,
  int* neofunc,
  int* ierr
);

void mask_check_state_c(
  const int* bit_mask,
  const int* n_mask_chunks,
  const int* i_gene,
  int* state,
  int* ierr
);

void mask_chunk_count_c(
  const int* n_genes,
  int* count,
  int* ierr
);

void calc_work_arr_paralog_subsets_size(
  const int* max_subset_size,
  const int* n_genes,
  int* work_array_size,
  const int* filtered_paralogs_mask,
  const int* n_mask_chunks,
  int* ierr
);

void filter_paralogs_by_pattern_dosage_effect(
  const double* gene_angles,
  const double* threshold,
  const int* n_genes,
  const int* n_families,
  const int* gene_to_fam,
  int* masks,
  const int* n_mask_chunks,
  int* ierr
);

void filter_paralogs_by_pattern_subfunctionalization_c(
  const double* gene_angles,
  const double* threshold,
  const int* n_genes,
  const int* n_families,
  const int* gene_to_fam,
  int* masks,
  const int* n_mask_chunks,
  int* ierr
);

void detect_subfunctionalization_c(
  const double* ancestor,
  const double* genes,
  const int* n_genes,
  const int* n_dims,
  const double* rdi_threshold,
  const int* filtered_paralogs_mask,
  const int* n_mask_chunks,
  int* n_results,
  const int* max_subset_size,
  int* work_arr_paralog_subsets,
  const int* n_paralog_subsets,
  int* active_mask,
  double* temp_paralog_vector,
  const double* paralog_norms,
  const int* sorted_paralog_norms_perm,
  double* temp_work_array,
  int* ierr
);

void detect_dosage_effect_c(
  const double* ancestor,
  const double* genes,
  const int* n_genes,
  const int* n_dims,
  const int* filtered_paralogs_mask,
  const int* n_mask_chunks,
  int* n_results,
  const int* max_subset_size,
  int* work_arr_paralog_subsets,
  const int* n_paralog_subsets,
  int* active_mask,
  double* temp_paralog_vector,
  const double* max_angle,
  const double* gain_gamma,
  int* ierr
);

void mean_vector_c(
  const double* expression_vectors,
  const int* n_axes,
  const int* n_genes,
  const int* gene_indices,
  const int* n_selected_genes,
  double* centroid_col,
  int* ierr
);

void group_centroid_c(
  const double* expression_vectors,
  const int* n_axes,
  const int* n_genes,
  const int* gene_to_family,
  const int* n_families,
  double* centroid_matrix,
  const char* mode,
  const int* ortholog_set,
  int* selected_indices,
  const int* selected_indices_len,
  int* ierr
);

void build_kd_index_C(
  const double* points,
  const int* num_dimensions,
  const int* num_points,
  int* kd_indices,
  const int* dimension_order,
  int* workspace,
  double* value_buffer,
  int* permutation,
  int* left_stack,
  int* right_stack,
  int* ierr
);

void which_c(
  const int* mask,
  const int* n,
  int* idx_out,
  const int* m_max,
  int* m_out,
  int* ierr
);

void loess_smooth_2d_c(
  const int* n_total,
  const int* n_target,
  const double* x_ref,
  const double* y_ref,
  const int* indices_used,
  int* n_used,
  const double* x_query,
  const double* kernel_sigma,
  const double* kernel_cutoff,
  double* y_out,
  int* ierr
);

void deserialize_int_nd_C(
  int* arr,
  const int* arr_size,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void deserialize_real_nd_C(
  double* arr,
  const int* arr_size,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void deserialize_char_nd_C(
  char* ascii_arr,
  const int* clen,
  const int* total_array_size,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void serialize_int_nd_C(
  const void* arr,
  const int* dims,
  const int* ndim,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void serialize_real_nd_C(
  const void* arr,
  const int* dims,
  const int* ndim,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void serialize_char_nd_C(
  const char* ascii_arr,
  const int* dims,
  const int* ndim,
  const int* clen,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void deserialize_logical_nd_C(
  int* arr,
  const int* arr_size,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void serialize_logical_nd_C(
  const int* arr,
  const int* dims,
  const int* ndim,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void serialize_complex_nd_C(
  const void* arr,
  const int* dims,
  const int* ndim,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void deserialize_complex_nd_C(
  void* arr,
  const int* arr_size,
  const char* filename_ascii,
  const int* fn_len,
  int* ierr
);

void get_array_metadata_C(
  const char* filename_ascii,
  const int* fn_len,
  int* dims_out,
  const int* dims_out_capacity,
  int* ndims,
  int* ierr,
  int* clen
);

void build_bst_index_C(
  const double* values,
  const int* num_values,
  int* sorted_indices,
  int* left_stack,
  int* right_stack,
  int* ierr
);

void bst_range_query_C(
  const double* values,
  const int* sorted_indices,
  const int* num_values,
  const double* lower_bound,
  const double* upper_bound,
  int* output_indices,
  int* num_matches,
  int* ierr
);

void build_spherical_kd_C(
  const double* vectors,
  const int* num_dimensions,
  const int* num_vectors,
  int* sphere_indices,
  const int* dimension_order,
  int* workspace,
  double* value_buffer,
  int* permutation,
  int* left_stack,
  int* right_stack,
  int* ierr
);

void read_expression_vectors_tsv_C(
  const char* file_list_raw,
  const int* file_list_len,
  const int* n_files,
  const char* gene_ids_raw,
  const int* gene_ids_len,
  const int* n_genes,
  double* expression_vectors,
  const int* n_samples,
  const int* n_header_rows,
  const int* gene_col,
  const int* value_cols,
  const int* n_value_cols,
  const char* delimiter_raw,
  int* ierr
);

void read_gene_ids_from_tsv_file_C(
  const char* filename_raw,
  const int* fn_len,
  char* gene_ids_raw,
  int* gene_ids_len,
  int* n_genes,
  const int* n_header_rows,
  const int* gene_col,
  int* ierr
);

void read_orthofinder_file_C(
  const char* filename_raw,
  const int* fn_len,
  const char* gene_ids_raw,
  const int* gene_ids_len,
  const int* n_genes,
  char* family_ids_raw,
  int* family_ids_len,
  int* n_families,
  int* gene_to_fam,
  int* ierr
);

void filter_unassigned_genes_C(
  const char* gene_ids_raw,
  const int* gene_ids_len,
  const int* n_genes,
  const int* gene_to_fam,
  int* mask,
  int* n_genes_kept,
  int* ierr
);

void validate_data_structure_C(
  const int* n_genes,
  const int* n_families,
  const int* n_samples,
  const char* gene_ids_raw,
  const int* gene_ids_len,
  const char* gene_family_ids_raw,
  const int* fam_len,
  const int* gene_to_fam,
  const double* expression_vectors,
  const double* family_centroids,
  const double* shift_vectors,
  int* ierr
);

void validate_gene_to_family_mapping_C(
  const int* gene_to_fam,
  const int* n_genes,
  const int* n_families,
  int* ierr
);

void validate_expression_data_C(
  const double* expression_vectors,
  const int* n_genes,
  const int* n_samples,
  const int* check_non_negative,
  int* ierr
);

void validate_family_centroids_C(
  const double* family_centroids,
  const int* n_families,
  const int* n_samples,
  int* ierr
);

void validate_shift_vectors_C(
  const double* shift_vectors,
  const double* expression_vectors,
  const double* family_centroids,
  const int* gene_to_fam,
  const int* n_genes,
  const int* n_samples,
  const int* n_families,
  int* ierr
);

void validate_string_array_uniqueness_C(
  const char* str_arr,
  const int* str_len,
  const int* n_strings,
  int* ierr
);

void validate_all_data_C(
  const int* n_genes,
  const int* n_families,
  const int* n_samples,
  const char* gene_ids_raw,
  const int* gene_len,
  const char* gene_family_ids_raw,
  const int* fam_len,
  const int* gene_to_fam,
  const double* expression_vectors,
  const double* family_centroids,
  const double* shift_vectors,
  int* ierr
);
}
//' Calculate k-means clustering of factor trajectories
//'
//' @param trajectories Numeric vector of trajectories (factors x samples x timepoints)
//' @param centroids Numeric matrix to store cluster centroids (factors x clusters)
//' @param n_clusters Integer number of clusters
//' @param n_factors Integer number of factors
//' @param n_samples Integer number of samples
//' @param n_timepoints Integer number of timepoints
//' @param max_iterations Integer maximum number of k-means iterations
//' @return List with centroids, labels, and label counts
// [[Rcpp::export]]
List tox_cluster_factor_trajectories_k_means_rcpp(NumericVector trajectories,
                                                  NumericMatrix centroids,
                                                  int n_clusters,
                                                  int n_factors,
                                                  int n_samples,
                                                  int n_timepoints,
                                                  int max_iterations = 300) {
    int n_points = n_samples * n_timepoints;
    IntegerVector labels(n_points);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;

    cluster_factor_trajectories_k_means_c(&n_clusters,
                                          trajectories.begin(),
                                          &n_factors,
                                          &n_samples,
                                          &n_timepoints,
                                          centroids.begin(),
                                          labels.begin(),
                                          label_counts.begin(),
                                          &max_iterations,
                                          &ierr);

   

    return List::create(Named("centroids") = centroids,
                        Named("labels") = labels,
                        Named("label_counts") = label_counts,
                        Named("ierr") = ierr);
}

//' Calculate k-means clustering of data points
//'
//' @param data_points Numeric matrix of data points (points x dimensions)
//' @param centroids Numeric matrix to store cluster centroids (clusters x dimensions)
//' @param n_clusters Integer number of clusters
//' @param n_points Integer number of data points
//' @param n_dims Integer number of dimensions
//' @param max_iterations Integer maximum number of k-means iterations
//' @return List with centroids, labels, and label counts
// [[Rcpp::export]]
List tox_k_means_clustering_rcpp(int n_clusters,
                                 NumericMatrix data_points,
                                 int n_points,
                                 int n_dims,
                                 NumericMatrix centroids,
                                 int max_iterations = 300) {
    IntegerVector labels(n_points);
    IntegerVector label_counts(n_clusters);
    int ierr = 0;

    
    k_means_clustering_c(&n_clusters,
                         data_points.begin(),
                         &n_points,
                         &n_dims,
                         centroids.begin(),
                         labels.begin(),
                         label_counts.begin(),
                         &max_iterations,
                         &ierr);



    return List::create(Named("centroids") = centroids,
                        Named("labels") = labels,
                        Named("label_counts") = label_counts,
                        Named("ierr") = ierr);
}

  //' Perform hierarchical linkage clustering on a distance matrix
  //'
  //' @param distances Numeric square distance matrix (n_points x n_points)
  //' @param method Linkage method: "average", "weighted", or "ward"
  //' @return List with merge indices, heights, cluster sizes, and error code
  // [[Rcpp::export]]
  List tox_linkage_clustering_rcpp(NumericMatrix distances,
                   std::string method) {
    int n_points = distances.nrow();
    IntegerVector merge_i(n_points - 1);
    IntegerVector merge_j(n_points - 1);
    NumericVector heights(n_points - 1);
    IntegerVector cluster_sizes(n_points - 1);
    int ierr = 0;

    char method_c[8] = {' '};
    size_t len = std::min(method.size(), size_t(8));
    std::memcpy(method_c, method.c_str(), len);

    linkage_clustering_c(distances.begin(),
               &n_points,
               merge_i.begin(),
               merge_j.begin(),
               heights.begin(),
               cluster_sizes.begin(),
               method_c,
               &ierr);

    return List::create(Named("merge_i") = merge_i,
              Named("merge_j") = merge_j,
              Named("heights") = heights,
              Named("cluster_sizes") = cluster_sizes,
              Named("ierr") = ierr);
  }


//' Calculate Euclidean distance between two vectors
//'
//' @param vec1 First numeric vector
//' @param vec2 Second numeric vector
//' @return Euclidean distance (double)
// [[Rcpp::export]]
double tox_euclidean_distance_rcpp(NumericVector vec1,
                                   NumericVector vec2) {
    int d = vec1.length();
    double result = 0.0;

    euclidean_distance_c(vec1.begin(),
                         vec2.begin(),
                         &d,
                         &result);

    return result;
}


//' Calculate relative axis changes from a shift vector

//' @param vec Numeric vector of shift values
//' @return List with contributions 
// [[Rcpp::export]]
List tox_relative_axes_changes_from_shift_vector_rcpp(NumericVector vec) {

    int n_axes = vec.length();
    NumericVector contributions(n_axes);
    int ierr = 0;

    relative_axes_changes_from_shift_vector_c(vec.begin(),
                                              &n_axes,
                                              contributions.begin(),
                                              &ierr);

    return List::create(Named("contributions") = contributions,
                        Named("ierr") = ierr);
}
 
//' Calculate relative axis expression from an expression vector
//'
//' @param vec Numeric vector of expression values
//' @return List with contributions
// [[Rcpp::export]]
List tox_relative_axes_expression_from_expression_vector_rcpp(NumericVector vec) {

    int n_axes = vec.length();
    NumericVector contributions(n_axes);
    int ierr = 0;

    relative_axes_expression_from_expression_vector_c(vec.begin(),
                                                      &n_axes,
                                                      contributions.begin(),
                                                      &ierr);

    return List::create(Named("contributions") = contributions,
                        Named("ierr") = ierr);
}


//' Calculate clock hand angle between two vectors
//'
//' @param v1 First numeric vector
//' @param v2 Second numeric vector
//' @param selected_axes_for_signed Integer vector of axes to consider for signed angle calculation
//' @return List with signed angle
// [[Rcpp::export]]
List tox_clock_hand_angle_between_vectors_rcpp(NumericVector v1,
                                               NumericVector v2,
                                               IntegerVector selected_axes_for_signed) {
    int n_dims = v1.length();
    double signed_angle = 0.0;
    int ierr = 0;

    clock_hand_angle_between_vectors_c(v1.begin(),
                                       v2.begin(),
                                       &n_dims,
                                       &signed_angle,
                                       selected_axes_for_signed.begin(),
                                       &ierr);

    return List::create(Named("signed_angle") = signed_angle,
                        Named("ierr") = ierr);
}

//' Calculate clock hand angles for shift vectors
//'
//' @param origins Numeric matrix of origin vectors (axes x vectors)
//' @param targets Numeric matrix of target vectors (axes x vectors)
//' @param vecs_selection_mask Integer vector indicating selected vectors
//' @param selected_axes_for_signed Integer vector indicating selected axes for signed angle calculation
//' @return List with signed angles 
// [[Rcpp::export]]
List tox_clock_hand_angles_for_shift_vectors_rcpp(NumericMatrix origins,
                                                  NumericMatrix targets,
                                                  IntegerVector vecs_selection_mask,
                                                  IntegerVector selected_axes_for_signed) {
    int n_dims = origins.nrow();
    int n_vecs = origins.ncol();
    NumericVector signed_angles(n_vecs);
    int ierr = 0;

    // Count selected vectors
    int n_selected_vecs = 0;
    for (int i = 0; i < n_vecs; ++i) {
        if (vecs_selection_mask[i] != 0) {
            ++n_selected_vecs;
        }
    }

    clock_hand_angles_for_shift_vectors_c(origins.begin(),
                                          targets.begin(),
                                          &n_dims,
                                          &n_vecs,
                                          vecs_selection_mask.begin(),
                                          &n_selected_vecs,
                                          selected_axes_for_signed.begin(),
                                          signed_angles.begin(),
                                          &ierr);

    return List::create(Named("signed_angles") = signed_angles,
                        Named("ierr") = ierr);
}


//' Calculate distances from genes to their family centroids
//'
//' @param genes Numeric vector of gene expression values
//' @param centroids Numeric vector of family centroid values
//' @param gene_to_fam Integer vector mapping each gene to its family
//' @param d Number of dimensions
//' @return Numeric vector of distances
// [[Rcpp::export]]
NumericVector tox_distance_to_centroid_rcpp(NumericVector genes,
                                            NumericVector centroids,
                                            IntegerVector gene_to_fam,
                                            int d) {
    int n_genes = genes.length() / d;
    int n_families = centroids.length() / d;

    NumericVector distances(n_genes);

    distance_to_centroid_c(&n_genes,
                           &n_families,
                           genes.begin(),
                           centroids.begin(),
                           gene_to_fam.begin(),
                           distances.begin(),
                           &d);

    return distances;
}

//' Determine shared residual range for JSD calculation (expert method)
//'
//'@param residual_pool Numeric vector of residuals from both sets
//'@param residual_pool_perm Integer vector for permutation of residual pool
//'@param residual_range_quantile Double quantile to determine shared residual range (e.g.,95.0 for 95th percentile)
//'@return List with shared residual range
// [[Rcpp::export]]
List tox_determine_shared_residual_range_expert_rcpp(NumericVector residual_pool,
                                                     IntegerVector residual_pool_perm,
                                                     double residual_range_quantile = 95.0) {

    int pool_size = residual_pool.size();
    double shared_R = 0.0;
    int ierr = 0;

    determine_shared_residual_range_expert_c(residual_pool.begin(),
                                             residual_pool_perm.begin(),
                                             &pool_size,
                                             &residual_range_quantile,
                                             &shared_R,
                                             &ierr);

    return List::create(Named("shared_R") = shared_R,
                        Named("ierr") = ierr);
}

//' Determine shared residual range for JSD calculation
//'
//'@param neighborhood_residuals_S1 Numeric vector of neighborhood residuals for set 1 (reps x neighbors x points)
//'@param neighborhood_residuals_S2 Numeric vector of neighborhood residuals for set 2
//'@param residual_range_quantile Double quantile to determine shared residual range (e.g.,95.0 for 95th percentile)
//'@return List with shared residual range
// [[Rcpp::export]]
List tox_determine_shared_residual_range_rcpp(NumericVector neighborhood_residuals_S1,
                                              NumericVector neighborhood_residuals_S2,
                                              double residual_range_quantile = 95.0) {

    IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];

    double shared_R = 0.0;
    int ierr = 0;

    determine_shared_residual_range_c(neighborhood_residuals_S1.begin(),
                                      neighborhood_residuals_S2.begin(),
                                      &n_reps_S1,
                                      &n_reps_S2,
                                      &n_neighbors,
                                      &n_points,
                                      &residual_range_quantile,
                                      &shared_R,
                                      &ierr);

    return List::create(Named("shared_R") = shared_R,
                        Named("ierr") = ierr);
}

//' Build residual histograms for JSD calculation
//'
//' @param neighborhood_residuals Numeric vector of neighborhood residuals (reps x neighbors x points)
//' @param shared_residual_range Double maximum residual value to consider for histogram binning
//' @param n_bins Integer number of histogram bins
//' @return List with counts, pmf and included number of residuals
// [[Rcpp::export]]
List tox_build_residual_histograms_rcpp(NumericVector neighborhood_residuals,
                                        double shared_residual_range,
                                        int n_bins) {

    IntegerVector dims = neighborhood_residuals.attr("dim");
    int n_reps = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];

    IntegerMatrix counts(n_points, n_bins);
    NumericMatrix pmf(n_points, n_bins);
    IntegerVector included_n_residuals(n_points);

    int ierr = 0;

    build_residual_histograms_c(neighborhood_residuals.begin(),
                                &n_reps,
                                &n_neighbors,
                                &n_points,
                                &shared_residual_range,
                                &n_bins,
                                counts.begin(),
                                pmf.begin(),
                                included_n_residuals.begin(),
                                &ierr);

    return List::create(Named("counts") = counts,
                        Named("pmf") = pmf,
                        Named("included_n_residuals") = included_n_residuals,
                        Named("ierr") = ierr);
}

//' Build residual histograms for JSD calculation with neighbor filtering
//'
//' @param neighborhood_residuals Numeric vector of neighborhood residuals (reps x neighbors x points)
//' @param shared_residual_range Double maximum residual value to consider for histogram binning
//' @param n_bins Integer number of histogram bins
//' @param neighbor_mask Integer vector indicating which neighbors to include (length should be equal to number of neighbors)
//' @return List with counts, pmf and included number of residuals
// [[Rcpp::export]]
List tox_build_residual_histograms_filtered_rcpp(NumericVector neighborhood_residuals,
                                                 double shared_residual_range,
                                                 int n_bins,
                                                 IntegerVector neighbor_mask) {

    IntegerVector dims = neighborhood_residuals.attr("dim");
    int n_reps = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];

    IntegerMatrix counts(n_points, n_bins);
    NumericMatrix pmf(n_points, n_bins);
    IntegerVector included_n_residuals(n_points);

    int ierr = 0;

    build_residual_histograms_filtered_c(neighborhood_residuals.begin(),
                                         &n_reps,
                                         &n_neighbors,
                                         &n_points,
                                         &shared_residual_range,
                                         &n_bins,
                                         counts.begin(),
                                         pmf.begin(),
                                         included_n_residuals.begin(),
                                         &ierr,
                                         neighbor_mask.begin());

    return List::create(Named("counts") = counts,
                        Named("pmf") = pmf,
                        Named("included_n_residuals") = included_n_residuals,
                        Named("ierr") = ierr);
}
//' Compute JSD divergence per reference point
//'
//' @param pmf_S1 Numeric matrix of probability mass functions for set 1 (points x bins)
//' @param pmf_S2 Numeric matrix of probability mass functions for set 2 (points x bins)
//' @return List with JSD divergences per point
// [[Rcpp::export]]
List tox_compute_divergence_per_reference_point_rcpp(NumericMatrix pmf_S1,
                                                     NumericMatrix pmf_S2) {

    int n_points = pmf_S1.nrow();
    int n_bins = pmf_S1.ncol();

    NumericVector js_divergences(n_points);
    int ierr = 0;

    compute_divergence_per_reference_point_c(pmf_S1.begin(),
                                             pmf_S2.begin(),
                                             &n_points,
                                             &n_bins,
                                             js_divergences.begin(),
                                             &ierr);

    return List::create(Named("js_divergences") = js_divergences,
                        Named("ierr") = ierr);
}
//' Compute weighted global JSD divergence
//'
//' @param js_divergences Numeric vector of JSD divergences per point
//' @param included_n_residuals_S1 Integer vector of included residuals for set 1
//' @param included_n_residuals_S2 Integer vector of included residuals for set 2
//' @return List with global JSD divergence and weights
// [[Rcpp::export]]
List tox_compute_weighted_global_divergence_rcpp(NumericVector js_divergences,
                                                 IntegerVector included_n_residuals_S1,
                                                 IntegerVector included_n_residuals_S2) {
    int n_points = js_divergences.size();

    double global_jsd = 0.0;
    NumericVector weights(n_points);
    int ierr = 0;

    compute_weighted_global_divergence_c(js_divergences.begin(),
                                         &n_points,
                                         included_n_residuals_S1.begin(),
                                         included_n_residuals_S2.begin(),
                                         &global_jsd,
                                         weights.begin(),
                                         &ierr);

    return List::create(Named("global_js_divergence") = global_jsd,
                        Named("weights") = weights,
                        Named("ierr") = ierr);
}
//' Perform permutation test for global JSD divergence
//'
//'@param neighborhood_residuals_S1 Numeric vector of neighborhood residuals for set 1 (reps x neighbors x points)
//'@param neighborhood_residuals_S2 Numeric vector of neighborhood residuals for set 2
//'@param global_jsd_observed Double observed global JSD divergence
//'@param n_bins Integer number of histogram bins
//'@param shared_residual_range Double maximum residual value to consider for histogram binning
//'@param n_permutations Integer number of permutations to perform
//'@param random_seed Integer random seed for permutation reproducibility
//'@return List with null distribution of JSD divergences and p-value
// [[Rcpp::export]]
List tox_gjct_permutation_test_rcpp(NumericVector neighborhood_residuals_S1,
                                    NumericVector neighborhood_residuals_S2,
                                    double global_jsd_observed,
                                    int n_bins,
                                    double shared_residual_range,
                                    int n_permutations,
                                    int random_seed) {

    IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];

    NumericVector jsd_null(n_permutations);
    double p_value = 0.0;
    int ierr = 0;

    gjct_permutation_test_c(neighborhood_residuals_S1.begin(),
                            neighborhood_residuals_S2.begin(),
                            &n_reps_S1,
                            &n_reps_S2,
                            &n_neighbors,
                            &n_points,
                            &global_jsd_observed,
                            &n_bins,
                            &shared_residual_range,
                            &n_permutations,
                            jsd_null.begin(),
                            &p_value,
                            &ierr,
                            &random_seed);

    return List::create(Named("jsd_null") = jsd_null,
                        Named("p_value")  = p_value,
                        Named("ierr")     = ierr);
}
//' Perform permutation test for global JSD divergence with neighbor filtering
//'
//'@param neighborhood_residuals_S1 Numeric vector of neighborhood residuals for set 1
//'@param neighborhood_residuals_S2 Numeric vector of neighborhood residuals for set 2
//'@param global_jsd_observed Double observed global JSD divergence
//'@param n_bins Integer number of histogram bins
//'@param shared_residual_range Double maximum residual value to consider for histogram binning
//'@param n_permutations Integer number of permutations to perform
//'@param random_seed Integer random seed for permutation reproducibility
//'@param neighbor_mask_S1 Integer vector of neighbor indices for set 1
//'@param neighbor_mask_S2 Integer vector of neighbor indices for set 2
//'@return List with null distribution of JSD divergences and  p-value
// [[Rcpp::export]]
List tox_gjct_permutation_test_filtered_rcpp(NumericVector neighborhood_residuals_S1,
                                             NumericVector neighborhood_residuals_S2,
                                             double global_jsd_observed,
                                             int n_bins,
                                             double shared_residual_range,
                                             int n_permutations,
                                             IntegerVector neighbor_mask_S1,
                                             IntegerVector neighbor_mask_S2,
                                             int random_seed) {

    IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];

    NumericVector jsd_null(n_permutations);
    double p_value = 0.0;
    int ierr = 0;

    gjct_permutation_test_filtered_c(neighborhood_residuals_S1.begin(),
                                     neighborhood_residuals_S2.begin(),
                                     &n_reps_S1,
                                     &n_reps_S2,
                                     &n_neighbors,
                                     &n_points,
                                     &global_jsd_observed,
                                     &n_bins,
                                     &shared_residual_range,
                                     &n_permutations,
                                     jsd_null.begin(),
                                     &p_value,
                                     &ierr,
                                     &random_seed,
                                     neighbor_mask_S1.begin(),
                                     neighbor_mask_S2.begin());

    return List::create(Named("jsd_null") = jsd_null,
                        Named("p_value")  = p_value,
                        Named("ierr")     = ierr);
}
//' Calculate neighborhood size for JSD calculation
//'
//' @param n_pool Integer number of pooled genes to consider for neighborhood construction
//' @param n_points Integer number of reference points
//' @param n_genes_S Integer number of genes in the dataset
//' @param mean_S Numeric vector of mean expression values for each gene
//' @param desired_size Integer desired neighborhood size (if 0, the function will determine the optimal size)
//' @return List with calculated neighborhood size
// [[Rcpp::export]]
List tox_calc_neighborhood_size_rcpp(int n_pool,
                                     int n_points,
                                     int n_genes_S,
                                     NumericVector mean_S,
                                     int desired_size = 0) {
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(&n_pool,
                             &n_points,
                             &n_genes_S,
                             mean_S.begin(),
                             &desired_size,
                             &n_neighbors,
                             &ierr);

    return List::create(Named("n_neighbors") = n_neighbors,
                        Named("ierr")        = ierr);
}

//' Construct neighborhoods for JSD calculation
//'
//' @param x_star Numeric vector of reference points
//' @param n_pool Integer number of pooled genes to consider for neighborhood construction
//' @param mean_S Numeric vector of mean expression values for each gene
//' @param resid_S Numeric matrix of residuals for each gene and reference point
//' @param desired_n_neighbors Integer desired neighborhood size (if 0, the function will determine the optimal size)
//' @return List with constructed neighborhoods and error code
// [[Rcpp::export]]
List tox_construct_neighborhoods_rcpp(NumericVector x_star,
                                      int n_pool,
                                      NumericVector mean_S,
                                      NumericMatrix resid_S,
                                      int desired_n_neighbors = 0) {

    int n_points  = x_star.size();
    int n_genes_S = mean_S.size();
    int n_reps_S  = resid_S.nrow();
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(&n_pool,
                             &n_points,
                             &n_genes_S,
                             mean_S.begin(),
                             &desired_n_neighbors,
                             &n_neighbors,
                             &ierr);

    if (ierr != 0) {
        NumericVector neigh_res(0);
        IntegerMatrix neigh_idx(0, n_points);

        return List::create(Named("neighborhood_residuals") = neigh_res,
                            Named("neighborhood_indices")   = neigh_idx,
                            Named("ierr")                   = ierr);
    }

    NumericVector neigh_res(n_reps_S * n_neighbors * n_points);
    IntegerMatrix neigh_idx(n_neighbors, n_points);

    construct_neighborhoods_c(&n_points,
                              x_star.begin(),
                              &n_genes_S,
                              mean_S.begin(),
                              &n_reps_S,
                              resid_S.begin(),
                              neigh_res.begin(),
                              neigh_idx.begin(),
                              &n_neighbors,
                              &ierr);

    neigh_res.attr("dim") = IntegerVector::create(
        n_reps_S,
        n_neighbors,
        n_points
    );

    return List::create(Named("neighborhood_residuals") = neigh_res,
                        Named("neighborhood_indices")   = neigh_idx,
                        Named("ierr")                   = ierr);
}
//' Compute gene means for JSD calculation
//'
//' @param expr Numeric matrix of gene expression values (replicates x genes)
//' @return List with gene means 
// [[Rcpp::export]]
List tox_compute_gene_means_rcpp(NumericMatrix expr) {

    int n_reps  = expr.nrow();
    int n_genes = expr.ncol();

    NumericVector means(n_genes);
    int ierr = 0;

    compute_gene_means_c(&n_genes,
                         &n_reps,
                         expr.begin(),
                         means.begin(),
                         &ierr);

    return List::create(Named("means") = means,
                        Named("ierr")  = ierr);
}

//' Compute residuals for JSD calculation
//'
//' @param expr Numeric matrix of gene expression values (replicates x genes)
//' @param means Numeric vector of gene means (length equal to number of genes)
//' @return List with residuals matrix 
// [[Rcpp::export]]
List tox_compute_residuals_rcpp(NumericMatrix expr,
                                NumericVector means) {

    int n_reps  = expr.nrow();
    int n_genes = expr.ncol();

    NumericMatrix resid(n_reps, n_genes);
    int ierr = 0;

    compute_residuals_c(&n_genes,
                        &n_reps,
                        expr.begin(),
                        means.begin(),
                        resid.begin(),
                        &ierr);

    return List::create(Named("resid") = resid,
                        Named("ierr")  = ierr);
}
//' Calculate pooled means for JSD calculation
//'@param mean_S1 Numeric vector of gene means for set 1
//'@param mean_S2 Numeric vector of gene means for set 2
//'@param n_points Integer number of reference points to calculate
//'@return List with pooled means 
// [[Rcpp::export]]
List tox_pool_means_rcpp(NumericVector mean_S1,
                         NumericVector mean_S2,
                         int n_points) {

    int n_genes_S1 = mean_S1.size();
    int n_genes_S2 = mean_S2.size();

    NumericVector x_star(n_points);
    int n_pool = 0;
    int ierr = 0;

    pool_means_c(&n_genes_S1,
                 mean_S1.begin(),
                 &n_genes_S2,
                 mean_S2.begin(),
                 &n_points,
                 &n_pool,
                 x_star.begin(),
                 &ierr);

    return List::create(Named("n_pool") = n_pool,
                        Named("x_star") = x_star,
                        Named("ierr")   = ierr);
}

//' Calculate pooled means for JSD calculation using expert method
//'
//'@param pooled_means Numeric vector of pooled means (length equal to number of reference points)
//'@param pooled_perm Integer vector for permutation of pooled means
//'@param n_points Integer number of reference points to calculate
//'@return List with pooled means
// [[Rcpp::export]]
List tox_pool_means_expert_rcpp(NumericVector pooled_means,
                                IntegerVector pooled_perm,
                                int n_points) {

    NumericVector x_star(n_points);
    int pool_size = pooled_means.size();
    int n_pool = 0;
    int ierr = 0;

    pool_means_expert_c(pooled_means.begin(),
                        pooled_perm.begin(),
                        &pool_size,
                        &n_points,
                        &n_pool,
                        x_star.begin(),
                        &ierr);

    return List::create(Named("n_pool") = n_pool,
                        Named("x_star") = x_star,
                        Named("ierr")   = ierr);
}
//' Compute JSD divergences for each reference point and global JSD divergence
//'
//' @param family_idx Integer index of the gene family being analyzed
//' @param gene_to_family_S1 Integer vector mapping each gene in set 1 to its family
//' @param gene_to_family_S2 Integer vector mapping each gene in set 2 to its family
//' @param neighborhood_residuals_S1 Numeric vector of neighborhood residuals for set 1 (reps x neighbors x points)
//' @param neighborhood_residuals_S2 Numeric vector of neighborhood residuals for set 2 (reps x neighbors x points)
//' @param neighborhood_genes_S1 Integer matrix of gene indices for neighborhoods in set 1 (neighbors x points)
//' @param neighborhood_genes_S2 Integer matrix of gene indices for neighborhoods in set 2 (neighbors x points)
//' @param n_bins Integer number of histogram bins for JSD calculation
//' @param shared_residual_range Double maximum residual value to consider for histogram binning
//' @return List with JSD divergences per point, included number of replicates, total included replicates, global JSD divergence, weights, and error code
// [[Rcpp::export]]
List tox_fjct_compute_jsd_alloc_rcpp(int family_idx,
                                     IntegerVector gene_to_family_S1,
                                     IntegerVector gene_to_family_S2,
                                     NumericVector neighborhood_residuals_S1,
                                     NumericVector neighborhood_residuals_S2,
                                     IntegerMatrix neighborhood_genes_S1,
                                     IntegerMatrix neighborhood_genes_S2,
                                     int n_bins,
                                     double shared_residual_range) {

    IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];
    NumericVector jsd(n_points);
    IntegerVector inc1(n_points);
    IntegerVector inc2(n_points);
    int n_genes_S1 = gene_to_family_S1.size();
    int n_genes_S2 = gene_to_family_S2.size();
    int total_included = 0;
    double global_jsd = 0.0;
    NumericVector weights(n_points);
    int ierr = 0;

    fjct_compute_jsd_c(&family_idx,
                       gene_to_family_S1.begin(),
                       gene_to_family_S2.begin(),
                       &n_genes_S1,
                       &n_genes_S2,
                       neighborhood_residuals_S1.begin(),
                       neighborhood_residuals_S2.begin(),
                       neighborhood_genes_S1.begin(),
                       neighborhood_genes_S2.begin(),
                       &n_reps_S1,
                       &n_reps_S2,
                       &n_neighbors,
                       &n_points,
                       &n_bins,
                       &shared_residual_range,
                       jsd.begin(),
                       inc1.begin(),
                       inc2.begin(),
                       &total_included,
                       &global_jsd,
                       weights.begin(),
                       &ierr);

    return List::create(Named("js_divergences") = jsd,
                        Named("included_n_reps_S1") = inc1,
                        Named("included_n_reps_S2") = inc2,
                        Named("total_included_n_reps") = total_included,
                        Named("global_js_divergence") = global_jsd,
                        Named("weights") = weights,
                        Named("ierr") = ierr);
}

//' Compute JSD divergences for each reference point and global JSD divergence using expert method
//'
//' @param neighborhood_residuals_S1 Numeric vector of neighborhood residuals for set 1 (reps x neighbors x points)
//' @param neighborhood_residuals_S2 Numeric vector of neighborhood residuals for set 2 (reps x neighbors x points)
//' @param neighbor_mask_S1 Integer matrix of neighbor masks for set 1 (neighbors x points)
//' @param neighbor_mask_S2 Integer matrix of neighbor masks for set 2 (neighbors x points)
//' @param n_bins Integer number of histogram bins for JSD calculation
//' @param shared_residual_range Double maximum residual value to consider for histogram binning
//' @return List with JSD divergences per point, included number of replicates, total included replicates, global JSD divergence and weights
// [[Rcpp::export]]
List tox_fjct_compute_jsd_expert_rcpp(NumericVector neighborhood_residuals_S1,
                                      NumericVector neighborhood_residuals_S2,
                                      IntegerMatrix neighbor_mask_S1,
                                      IntegerMatrix neighbor_mask_S2,
                                      int n_bins,
                                      double shared_residual_range) {

    IntegerVector dims = neighborhood_residuals_S1.attr("dim");
    int n_reps_S1 = dims[0];
    int n_neighbors = dims[1];
    int n_points = dims[2];
    dims = neighborhood_residuals_S2.attr("dim");
    int n_reps_S2 = dims[0];
    NumericVector jsd(n_points);
    IntegerVector inc1(n_points);
    IntegerVector inc2(n_points);
    int total_included = 0;
    double global_jsd = 0.0;
    NumericVector weights(n_points);

    NumericMatrix pmf_S1(n_points, n_bins);
    NumericMatrix pmf_S2(n_points, n_bins);
    IntegerMatrix tmp_counts(n_points, n_bins);

    int ierr = 0;

    fjct_compute_jsd_expert_c(neighborhood_residuals_S1.begin(),
                              neighborhood_residuals_S2.begin(),
                              &n_reps_S1,
                              &n_reps_S2,
                              &n_neighbors,
                              &n_points,
                              neighbor_mask_S1.begin(),
                              neighbor_mask_S2.begin(),
                              &n_bins,
                              &shared_residual_range,
                              jsd.begin(),
                              inc1.begin(),
                              inc2.begin(),
                              &total_included,
                              &global_jsd,
                              weights.begin(),
                              pmf_S1.begin(),
                              pmf_S2.begin(),
                              tmp_counts.begin(),
                              &ierr);

    return List::create(Named("js_divergences") = jsd,
                        Named("included_n_reps_S1") = inc1,
                        Named("included_n_reps_S2") = inc2,
                        Named("total_included_n_reps") = total_included,
                        Named("global_js_divergence") = global_jsd,
                        Named("weights") = weights,
                        Named("pmf_S1") = pmf_S1,
                        Named("pmf_S2") = pmf_S2,
                        Named("tmp_counts") = tmp_counts,
                        Named("ierr") = ierr);
}

//' Compute contribution scores for each family based on global JSD divergences and included replicates
//'
//' @param global_js_divergences Numeric vector of global JSD divergences for each family
//' @param total_included_n_reps_per_f Integer vector of total included replicates for each family
//' @return List with support weights and contribution scores
// [[Rcpp::export]]
List tox_fjct_compute_contribution_scores_rcpp(NumericVector global_js_divergences,
                                               IntegerVector total_included_n_reps_per_f) {

    int k_families = global_js_divergences.size();
    NumericVector support_weights(k_families);
    NumericVector contribution_scores(k_families);
    int ierr = 0;

    fjct_compute_contribution_scores_c(global_js_divergences.begin(),
                                       total_included_n_reps_per_f.begin(),
                                       &k_families,
                                       support_weights.begin(),
                                       contribution_scores.begin(),
                                       &ierr);

    return List::create(Named("support_weights") = support_weights,
                        Named("contribution_scores") = contribution_scores,
                        Named("ierr") = ierr);
}



//' Calculate Tissue Versatility
//'
//' @param expression_vectors Numeric matrix of expression vectors (axes x vectors)
//' @param vector_selection Integer vector indicating selected vectors
//' @param axis_selection Integer vector indicating selected axes
//' @return List with tissue versatilities, angles, and selection counts
// [[Rcpp::export]]
List tox_calculate_tissue_versatility_rcpp(NumericMatrix expression_vectors,
                                           IntegerVector vector_selection,
                                           IntegerVector axis_selection) {

    int n_axes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();
    int n_selected_vectors = sum(vector_selection);
    int n_selected_axes = sum(axis_selection);

    NumericVector tissue_versatilities(n_selected_vectors);
    NumericVector tissue_angles_deg(n_selected_vectors);
    int ierr = 0;

    compute_tissue_versatility_c(&n_axes,
                                 &n_vectors,
                                 expression_vectors.begin(),
                                 vector_selection.begin(),
                                 &n_selected_vectors,
                                 axis_selection.begin(),
                                 &n_selected_axes,
                                 tissue_versatilities.begin(),
                                 tissue_angles_deg.begin(),
                                 &ierr);

    return List::create(Named("tissue_versatilities") = tissue_versatilities,
                        Named("tissue_angles_deg") = tissue_angles_deg,
                        Named("n_selected_vectors") = n_selected_vectors,
                        Named("n_selected_axes") = n_selected_axes,
                        Named("ierr") = ierr);
}


//' Normalize gene expression values by standard deviation
//'
//' @param input Numeric matrix (genes x tissues)
//' @return List with normalized output matrix 
// [[Rcpp::export]]
List tox_normalize_by_std_dev_rcpp(NumericMatrix input) {

    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    normalize_by_std_dev_c(&n_genes,
                           &n_tissues,
                           input.begin(),
                           output.begin(),
                           &ierr);

    return List::create(Named("output_vector") = output,
                        Named("ierr") = ierr);
}


//' Perform quantile normalization of gene expression values
//'
//' @param input Numeric matrix (genes x tissues)
//' @return List with quantile-normalized output matrix and intermediate results
// [[Rcpp::export]]
List tox_quantile_normalization_rcpp(NumericMatrix input) {

    int n_genes = input.nrow();
    int n_tissues = input.ncol();

    int max_stack = static_cast<int>(std::ceil(std::log2(n_genes)) + 10);

    NumericMatrix output(n_genes, n_tissues);
    NumericVector temp_col(n_genes);
    NumericVector rank_means(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(max_stack);
    IntegerVector stack_right(max_stack);
    int ierr = 0;

    quantile_normalization_c(&n_genes,
                             &n_tissues,
                             input.begin(),
                             output.begin(),
                             temp_col.begin(),
                             rank_means.begin(),
                             perm.begin(),
                             stack_left.begin(),
                             stack_right.begin(),
                             &max_stack,
                             &ierr);

    return List::create(Named("output_vector") = output,
                        Named("rank_means") = rank_means,
                        Named("perm") = perm,
                        Named("ierr") = ierr);
}

//' Apply log2(x + 1) transformation to gene expression values
//'
//' @param input Numeric matrix (genes x tissues)
//' @return List with log2-transformed output matrix
// [[Rcpp::export]]
List tox_log2_transformation_rcpp(NumericMatrix input) {

    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    NumericMatrix output(n_genes, n_tissues);
    int ierr = 0;

    log2_transformation_c(&n_genes,
                          &n_tissues,
                          input.begin(),
                          output.begin(),
                          &ierr);

    return List::create(Named("output_vector") = output,
                        Named("ierr") = ierr);
}

  //' Normalize a numeric vector to unit length
  //'
  //' @param vector Numeric vector
  //' @return List with normalized vector and error code
  // [[Rcpp::export]]
  List tox_normalize_unit_length_rcpp(NumericVector vector) {

    int n_dims = vector.size();
    int ierr = 0;

    normalize_unit_length_c(vector.begin(),
                &n_dims,
                &ierr);

    return List::create(Named("vector") = vector,
              Named("ierr") = ierr);
  }

//' Compute scalar baselines for a factor and dependent variable
//'
//' @param factor Numeric vector of factor values
//' @param dependent Numeric vector of dependent variable values
//' @param mode String indicating baseline calculation mode (e.g., 'raw', 'min', 'mean')
//' @return List with factor baseline, dependent baseline, and error code
// [[Rcpp::export]]
List tox_compute_baselines_factor_dependent_rcpp(NumericVector factor,
                                                 NumericVector dependent,
                                                 std::string mode) {
  int n_timepoints = factor.size();
  
  // Fortran expects mode as a char array (e.g., 'raw', 'min', 'mean'), pad/truncate to 8 chars
  char mode_c[8] = {' '};
  size_t len = std::min(mode.size(), size_t(8));
  std::memcpy(mode_c, mode.c_str(), len);

  double factor_baseline = 0.0;
  double dependent_baseline = 0.0;
  int ierr = 0;

  compute_baselines_factor_dependent_c(factor.begin(),
                                       dependent.begin(),
                                       &n_timepoints,
                                       mode_c,
                                       &factor_baseline,
                                       &dependent_baseline,
                                       &ierr);

  return List::create(Named("factor_baseline") = factor_baseline,
                      Named("dependent_baseline") = dependent_baseline,
                      Named("ierr") = ierr);
}

//' Calculate average expression across replicates for each tissue group
//'
//' @param input Numeric matrix (genes x samples)
//' @param group_s Integer vector: start column index for each group
//' @param group_c Integer vector: number of columns per group
//' @return List with averaged output matrix 
// [[Rcpp::export]]
List tox_calc_tiss_avg_rcpp(NumericMatrix input,
                            IntegerVector group_s,
                            IntegerVector group_c) {

    int n_gene = input.nrow();
    int n_grps = group_s.size();
    NumericMatrix output(n_gene, n_grps);
    int ierr = 0;

    calc_tiss_avg_c(&n_gene,
                    &n_grps,
                    group_s.begin(),
                    group_c.begin(),
                    input.begin(),
                    output.begin(),
                    &ierr);

    return List::create(Named("output_vector") = output,
                        Named("ierr") = ierr);
}


//' Calculate fold change between control and condition columns
//'
//' @param input Numeric matrix (genes x samples)
//' @param control_cols Integer vector of control column indices
//' @param cond_cols Integer vector of condition column indices
//' @return List with fold change output matrix 
// [[Rcpp::export]]
List tox_calc_fchange_rcpp(NumericMatrix input,
                           IntegerVector control_cols,
                           IntegerVector cond_cols) {

    int n_genes = input.nrow();
    int n_cols = input.ncol();
    int n_pairs = control_cols.size();
    NumericMatrix output(n_genes, n_pairs);
    int ierr = 0;

    calc_fchange_c(&n_genes,
                   &n_cols,
                   &n_pairs,
                   control_cols.begin(),
                   cond_cols.begin(),
                   input.begin(),
                   output.begin(),
                   &ierr);

    return List::create(Named("output_vector") = output,
                        Named("ierr") = ierr);
}


//' Complete normalization pipeline for gene expression data
//'
//' @param input Numeric matrix (genes x tissues)
//' @param group_s Integer vector: start column index for each group
//' @param group_c Integer vector: number of columns per group
//' @return List with pipeline output matrices
// [[Rcpp::export]]
List tox_normalization_pipeline_rcpp(NumericMatrix input,
                                     IntegerVector group_s,
                                     IntegerVector group_c) {

    int n_genes = input.nrow();
    int n_tissues = input.ncol();
    int n_grps = group_s.size();

    int max_stack = static_cast<int>(std::ceil(std::log2(n_genes)) + 10);

    NumericMatrix buf_stddev(n_genes, n_tissues);
    NumericMatrix buf_quant(n_genes, n_tissues);
    NumericMatrix buf_avg(n_genes, n_grps);
    NumericMatrix buf_log(n_genes, n_grps);
    NumericVector temp_col(n_genes);
    NumericVector rank_means(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(max_stack);
    IntegerVector stack_right(max_stack);
    int ierr = 0;

    normalization_pipeline_c(&n_genes,
                             &n_tissues,
                             input.begin(),
                             buf_stddev.begin(),
                             buf_quant.begin(),
                             buf_avg.begin(),
                             buf_log.begin(),
                             temp_col.begin(),
                             rank_means.begin(),
                             perm.begin(),
                             stack_left.begin(),
                             stack_right.begin(),
                             &max_stack,
                             group_s.begin(),
                             group_c.begin(),
                             &n_grps,
                             &ierr);

    return List::create(Named("buf_stddev") = buf_stddev,
                        Named("buf_quant") = buf_quant,
                        Named("buf_avg") = buf_avg,
                        Named("buf_log") = buf_log,
                        Named("rank_means") = rank_means,
                        Named("perm") = perm,
                        Named("ierr") = ierr);
}


            //' Normalize a single time series using min-max scaling
            //'
            //' @param v Numeric vector representing one time series
            //' @return List with normalized vector, status code, and error code
            // [[Rcpp::export]]
            List tox_normalize_variable_timeseries_rcpp(NumericVector v) {

              int n_points = v.size();
              NumericVector v_norm(n_points);
              int ierr = 0;
              int status = 0;

              normalize_variable_timeseries_C(v.begin(),
                                              v_norm.begin(),
                                              &n_points,
                                              &ierr,
                                              &status);

              return List::create(Named("v_norm") = v_norm,
                                  Named("status") = status,
                                  Named("ierr") = ierr);
            }


            //' Normalize all factors in a single trajectory independently
            //'
            //' @param trajectory Numeric matrix (timepoints x factors)
            //' @return List with normalized trajectory, status code, and error code
            // [[Rcpp::export]]
            List tox_normalize_single_trajectory_rcpp(NumericMatrix trajectory) {

              int n_timepoints = trajectory.nrow();
              int n_factors = trajectory.ncol();
              NumericMatrix trajectory_norm(n_timepoints, n_factors);
              int ierr = 0;
              int status = 0;

              normalize_single_trajectory_C(trajectory.begin(),
                                            trajectory_norm.begin(),
                                            &n_factors,
                                            &n_timepoints,
                                            &ierr,
                                            &status);

              return List::create(Named("traj_norm") = trajectory_norm,
                                  Named("status") = status,
                                  Named("ierr") = ierr);
            }


            //' Normalize all trajectories across factors, samples, and timepoints
            //'
            //' @param trajectories Numeric array flattened as vector (factors x samples x timepoints)
            //' @param n_factors Integer number of factors
            //' @param n_samples Integer number of samples
            //' @param n_timepoints Integer number of timepoints
            //' @return List with normalized trajectories, status code, and error code
            // [[Rcpp::export]]
            List tox_normalize_all_trajectories_rcpp(NumericVector trajectories,
                                                     int n_factors,
                                                     int n_samples,
                                                     int n_timepoints) {

              NumericVector trajectories_norm(n_factors * n_samples * n_timepoints);
              int ierr = 0;
              int status = 0;

              normalize_all_trajectories_C(trajectories.begin(),
                                           trajectories_norm.begin(),
                                           &n_factors,
                                           &n_samples,
                                           &n_timepoints,
                                           &ierr,
                                           &status);

              return List::create(Named("traj_norm") = trajectories_norm,
                                  Named("status") = status,
                                  Named("ierr") = ierr);
            }


//' Calculate shift vector field for gene expression vectors
//'
//' @param expression_vectors Numeric matrix (axes x genes)
//' @param family_centroids Numeric matrix (axes x families)
//' @param gene_to_centroid Integer vector mapping each gene to its centroid
//' @return List with shift vectors 
// [[Rcpp::export]]
List tox_compute_shift_vector_field_rcpp(NumericMatrix expression_vectors,
                                         NumericMatrix family_centroids,
                                         IntegerVector gene_to_centroid) {


    int n_axes_genes = expression_vectors.nrow();
    int n_vectors = expression_vectors.ncol();
    int n_axes_centroids = family_centroids.nrow();
    int n_families = family_centroids.ncol();

    NumericMatrix shift_vectors(2 * n_axes_genes, n_vectors);
    int ierr = 0;

    compute_shift_vector_field_c(&n_axes_genes,
                                 &n_vectors,
                                 &n_families,
                                 expression_vectors.begin(),
                                 family_centroids.begin(),
                                 gene_to_centroid.begin(),
                                 shift_vectors.begin(),
                                 &ierr);

    NumericVector flat(shift_vectors.begin(), shift_vectors.end());

    return List::create(Named("shift_vectors") = flat,
                        Named("ierr") = ierr);
}



//' Compute the element-wise mean for a given set of gene expression vectors
//'
//' @param expression_vectors Numeric matrix (axes x genes)
//' @param gene_indices Integer vector of gene indices to include in the mean
//' @return List with centroid vector 
// [[Rcpp::export]]
List tox_mean_vector_rcpp(NumericMatrix expression_vectors,
                          IntegerVector gene_indices) {


    int n_axes = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int n_selected_genes = gene_indices.length();

    NumericVector centroid_col(n_axes);
    int ierr = 0;

    mean_vector_c(expression_vectors.begin(),
                  &n_axes,
                  &n_genes,
                  gene_indices.begin(),
                  &n_selected_genes,
                  centroid_col.begin(),
                  &ierr);

    return List::create(Named("centroid_col") = centroid_col,
                        Named("ierr") = ierr);
}


//' Calculate gene family centroids
//'
//' @param expression_vectors Numeric matrix (axes x genes)
//' @param gene_to_family Integer vector mapping genes to family IDs
//' @param n_families Integer number of gene families
//' @param ortholog_set Integer vector indicating ortholog subset
//' @param mode Character string for mode ('all' or 'ortho')
//' @return List with centroid matrix
// [[Rcpp::export]]
List tox_group_centroid_rcpp(NumericMatrix expression_vectors,
                             IntegerVector gene_to_family,
                             int n_families,
                             IntegerVector ortholog_set,
                             String mode) {


    int n_axes = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();

    NumericMatrix centroid_matrix(n_axes, n_families);
    IntegerVector selected_indices(n_genes);
    int selected_indices_len = n_genes;
    int ierr = 0;

    // Convert String to std::string first
    std::string mode_str(mode);

    // Create a char array with size 10 (matching Fortran expectation)
    std::vector<char> mode_c(10, ' ');
    for (size_t i = 0; i < mode_str.size() && i < 10; ++i) {
        mode_c[i] = mode_str[i];
    }

    group_centroid_c(expression_vectors.begin(),
                     &n_axes,
                     &n_genes,
                     gene_to_family.begin(),
                     &n_families,
                     centroid_matrix.begin(),
                     mode_c.data(),
                     ortholog_set.begin(),
                     selected_indices.begin(),
                     &selected_indices_len,
                     &ierr);

    return List::create(Named("centroid_matrix") = centroid_matrix,
                        Named("ierr") = ierr);
}


//' Detect neofunctionalization for genes
//'
//' @param ancestors Numeric matrix of ancestor vectors (axes x families)
//' @param genes Numeric matrix of gene vectors (axes x genes)
//' @param gene_to_fam Integer vector mapping each gene to its family index
//' @param thresholds Numeric vector of per-axis thresholds
//' @return List with neofunctionalization matrix
// [[Rcpp::export]]
List tox_detect_neofunctionalization_rcpp(NumericMatrix ancestors,
                                          NumericMatrix genes,
                                          IntegerVector gene_to_fam,
                                          NumericVector thresholds) {


    int n_axes = ancestors.nrow();
    int n_families = ancestors.ncol();
    int n_genes = genes.ncol();

    IntegerMatrix neofunc_int(n_genes, n_axes);
    int ierr = 0;

    detect_neofunctionalization_c(ancestors.begin(),
                                  &n_families,
                                  genes.begin(),
                                  &n_axes,
                                  gene_to_fam.begin(),
                                  &n_genes,
                                  thresholds.begin(),
                                  neofunc_int.begin(),
                                  &ierr);

    LogicalMatrix neofunc(n_genes, n_axes);
    for (int j = 0; j < n_axes; ++j) {
      for (int i = 0; i < n_genes; ++i) {
        neofunc(i, j) = (neofunc_int(i, j) != 0);
      }
    }

    return List::create(Named("neofunc") = neofunc,
                        Named("ierr") = ierr);
  }


//' Check the state of a specific gene in a bit mask
//'
//' @param bit_mask Integer vector representing a chunked bit mask
//' @param i_gene Integer index of the gene to check
//' @return List with logical state 
// [[Rcpp::export]]
List tox_mask_check_state_rcpp(IntegerVector bit_mask, int i_gene) {

    int n_mask_chunks = bit_mask.size();
    int state = 0;
    int ierr = 0;

    mask_check_state_c(bit_mask.begin(),
                       &n_mask_chunks,
                       &i_gene,
                       &state,
                       &ierr);

    return List::create(Named("state") = (state != 0),
                        Named("ierr") = ierr);
}


//' Compute number of 32-bit chunks needed to encode n_genes in a bit mask
//'
//' @param n_genes Number of genes to encode
//' @return List with count of chunks
// [[Rcpp::export]]
List tox_mask_chunk_count_rcpp(int n_genes) {

    int count = 0;
    int ierr = 0;

    mask_chunk_count_c(&n_genes,
                       &count,
                       &ierr);

    return List::create(Named("count") = count,
                        Named("ierr") = ierr);
  }


//' Compute required work array size for paralog subset analysis
//'
//' @param max_subset_size Desired maximum subset size
//' @param n_genes Number of genes
//' @param filtered_paralogs_mask Chunked bit mask of filtered paralogs
//' @return List with adjusted max subset size and work array size
// [[Rcpp::export]]
List tox_calc_work_arr_paralog_subsets_size_rcpp(int max_subset_size,
                                                 int n_genes,
                                                 IntegerVector filtered_paralogs_mask) {

    int n_mask_chunks = filtered_paralogs_mask.size();
    int work_array_size = 0;
    int ierr = 0;

    calc_work_arr_paralog_subsets_size(&max_subset_size,
                                       &n_genes,
                                       &work_array_size,
                                       filtered_paralogs_mask.begin(),
                                       &n_mask_chunks,
                                       &ierr);

    return List::create(Named("actual_max_subset_size") = max_subset_size,
                        Named("work_array_size") = work_array_size,
                        Named("ierr") = ierr);
  }


//' Filter paralogs by dosage-effect pattern
//'
//' @param gene_angles Angles for all genes
//' @param threshold Filtering threshold
//' @param gene_to_fam Gene-to-family mapping
//' @param n_families Number of families
//' @return List with masks
// [[Rcpp::export]]
List tox_filter_paralogs_by_pattern_dosage_effect_rcpp(NumericVector gene_angles,
                                                       double threshold,
                                                       IntegerVector gene_to_fam,
                                                       int n_families) {

    int n_genes = gene_angles.size();
    int n_mask_chunks = (n_genes + 31) / 32;
    IntegerMatrix masks(n_mask_chunks, n_families);
    int ierr = 0;

    filter_paralogs_by_pattern_dosage_effect(gene_angles.begin(),
                                             &threshold,
                                             &n_genes,
                                             &n_families,
                                             gene_to_fam.begin(),
                                             masks.begin(),
                                             &n_mask_chunks,
                                             &ierr);

    return List::create(Named("masks") = masks,
                        Named("ierr") = ierr);
  }


//' Filter paralogs by subfunctionalization pattern
//'
//' @param gene_angles Angles for all genes
//' @param threshold Filtering threshold
//' @param gene_to_fam Gene-to-family mapping
//' @param n_families Number of families
//' @return List with masks
// [[Rcpp::export]]
List tox_filter_paralogs_by_pattern_subfunctionalization_rcpp(NumericVector gene_angles,
                                                              double threshold,
                                                              IntegerVector gene_to_fam,
                                                              int n_families) {
    int n_genes = gene_angles.size();
    int n_mask_chunks = (n_genes + 31) / 32;
    IntegerMatrix masks(n_mask_chunks, n_families);
    int ierr = 0;

    filter_paralogs_by_pattern_subfunctionalization_c(gene_angles.begin(),
                                                      &threshold,
                                                      &n_genes,
                                                      &n_families,
                                                      gene_to_fam.begin(),
                                                      masks.begin(),
                                                      &n_mask_chunks,
                                                      &ierr);

    return List::create(Named("masks") = masks,
                        Named("ierr") = ierr);
  }

  //' Detect subfunctionalization among paralogs
  //'
  //' @param ancestor Ancestor vector
  //' @param genes Gene matrix (dims x genes)
  //' @param rdi_threshold Maximum residual distance
  //' @param filtered_paralogs_mask Chunked bit mask of filtered paralogs
  //' @param max_subset_size Desired maximum subset size
  //' @param paralog_norms Norm of each paralog vector
  //' @param sorted_paralog_norms_perm Permutation of indices sorted by paralog norm
  //' @return List with n_results, results matrix and actual_max_subset_size
  // [[Rcpp::export]]
  List tox_detect_subfunctionalization_rcpp(NumericVector ancestor,
                                            NumericMatrix genes,
                                            double rdi_threshold,
                                            IntegerVector filtered_paralogs_mask,
                                            int max_subset_size,
                                            NumericVector paralog_norms,
                                            IntegerVector sorted_paralog_norms_perm) {
    int n_dims = ancestor.size();
    int n_genes = genes.ncol();
    int n_mask_chunks = filtered_paralogs_mask.size();
    int ierr = 0;

    int work_array_size = 0;
    calc_work_arr_paralog_subsets_size(&max_subset_size,
                                       &n_genes,
                                       &work_array_size,
                                       filtered_paralogs_mask.begin(),
                                       &n_mask_chunks,
                                       &ierr);

    if (work_array_size < 0) {
      work_array_size = 0;
    }

    int n_paralog_subsets = work_array_size;
    IntegerMatrix work_arr_paralog_subsets(n_mask_chunks, n_paralog_subsets);
    IntegerVector active_mask(n_mask_chunks);
    NumericVector temp_paralog_vector(n_dims);
    NumericVector temp_work_array(n_genes);
    int n_results = 0;

    detect_subfunctionalization_c(ancestor.begin(),
                                  genes.begin(),
                                  &n_genes,
                                  &n_dims,
                                  &rdi_threshold,
                                  filtered_paralogs_mask.begin(),
                                  &n_mask_chunks,
                                  &n_results,
                                  &max_subset_size,
                                  work_arr_paralog_subsets.begin(),
                                  &n_paralog_subsets,
                                  active_mask.begin(),
                                  temp_paralog_vector.begin(),
                                  paralog_norms.begin(),
                                  sorted_paralog_norms_perm.begin(),
                                  temp_work_array.begin(),
                                  &ierr);

    if (n_results < 0) {
      n_results = 0;
    }
    if (n_results > n_paralog_subsets) {
      n_results = n_paralog_subsets;
    }

    IntegerMatrix results(n_mask_chunks, n_results);
    for (int j = 0; j < n_results; ++j) {
      for (int i = 0; i < n_mask_chunks; ++i) {
        results(i, j) = work_arr_paralog_subsets(i, j);
      }
    }

    return List::create(Named("n_results") = n_results,
                        Named("results") = results,
                        Named("actual_max_subset_size") = max_subset_size,
                        Named("ierr") = ierr);
  }

  //' Detect dosage effect among paralogs
  //'
  //' @param ancestor Ancestor vector
  //' @param genes Gene matrix (dims x genes)
  //' @param filtered_paralogs_mask Chunked bit mask of filtered paralogs
  //' @param max_subset_size Desired maximum subset size
  //' @param gain_gamma Required gain threshold
  //' @param max_angle Maximum angle threshold (radians)
  //' @return List with n_results, results matrix and actual_max_subset_size
  // [[Rcpp::export]]
  List tox_detect_dosage_effect_rcpp(NumericVector ancestor,
                                     NumericMatrix genes,
                                     IntegerVector filtered_paralogs_mask,
                                     int max_subset_size,
                                     double gain_gamma,
                                     double max_angle) {
    int n_dims = ancestor.size();
    int n_genes = genes.ncol();
    int n_mask_chunks = filtered_paralogs_mask.size();
    int ierr = 0;

    int work_array_size = 0;
    calc_work_arr_paralog_subsets_size(&max_subset_size,
                                       &n_genes,
                                       &work_array_size,
                                       filtered_paralogs_mask.begin(),
                                       &n_mask_chunks,
                                       &ierr);

    if (work_array_size < 0) {
      work_array_size = 0;
    }

    int n_paralog_subsets = work_array_size;
    IntegerMatrix work_arr_paralog_subsets(n_mask_chunks, n_paralog_subsets);
    IntegerVector active_mask(n_mask_chunks);
    NumericVector temp_paralog_vector(n_dims);
    int n_results = 0;

    detect_dosage_effect_c(ancestor.begin(),
                           genes.begin(),
                           &n_genes,
                           &n_dims,
                           filtered_paralogs_mask.begin(),
                           &n_mask_chunks,
                           &n_results,
                           &max_subset_size,
                           work_arr_paralog_subsets.begin(),
                           &n_paralog_subsets,
                           active_mask.begin(),
                           temp_paralog_vector.begin(),
                           &max_angle,
                           &gain_gamma,
                           &ierr);

    if (n_results < 0) {
      n_results = 0;
    }
    if (n_results > n_paralog_subsets) {
      n_results = n_paralog_subsets;
    }

    IntegerMatrix results(n_mask_chunks, n_results);
    for (int j = 0; j < n_results; ++j) {
      for (int i = 0; i < n_mask_chunks; ++i) {
        results(i, j) = work_arr_paralog_subsets(i, j);
      }
    }

    return List::create(Named("n_results") = n_results,
                        Named("results") = results,
                        Named("actual_max_subset_size") = max_subset_size,
                        Named("ierr") = ierr);
  }


//' Compute family scaling factors using LOESS smoothing
//'
//' @param distances Numeric vector of gene distances
//' @param gene_to_fam Integer vector mapping genes to family indices
//' @param n_families Integer number of families
//' @return List with scaling factors, and loess fit
// [[Rcpp::export]]
List tox_compute_family_scaling_rcpp(NumericVector distances,
                                     IntegerVector gene_to_fam,
                                     int n_families) {
    int n_genes = distances.size();
    NumericVector dscale(n_families);
    NumericVector loess_x(n_families);
    NumericVector loess_y(n_families);
    IntegerVector indices_used(n_families);
    int ierr = 0;

    compute_family_scaling_c(&n_genes,
                             &n_families,
                             distances.begin(),
                             gene_to_fam.begin(),
                             dscale.begin(),
                             loess_x.begin(),
                             loess_y.begin(),
                             indices_used.begin(),
                             &ierr);

    return List::create(Named("dscale") = dscale,
                        Named("loess_x") = loess_x,
                        Named("loess_y") = loess_y,
                        Named("indices_used") = indices_used,
                        Named("ierr") = ierr);
}

//' Expert version of compute_family_scaling with user-provided work arrays
//'
//' @param n_families Integer number of families
//' @param distances Numeric vector of gene distances
//' @param gene_to_fam Integer vector mapping genes to family indices
//' @param perm_tmp Pre-allocated permutation array
//' @param stack_left_tmp Pre-allocated stack array (left)
//' @param stack_right_tmp Pre-allocated stack array (right)
//' @param family_distances Pre-allocated work array for family distances
//' @return List with scaling factors, loess fit, amd work array
// [[Rcpp::export]]
List tox_compute_family_scaling_expert_rcpp(int n_families,
                                            NumericVector distances,
                                            IntegerVector gene_to_fam,
                                            IntegerVector perm_tmp,
                                            IntegerVector stack_left_tmp,
                                            IntegerVector stack_right_tmp,
                                            NumericVector family_distances) {
    int n_genes = distances.size();
    NumericVector dscale(n_families);
    NumericVector loess_x(n_families);
    NumericVector loess_y(n_families);
    IntegerVector indices_used(n_families);
    int ierr = 0;

    compute_family_scaling_expert_c(&n_genes,
                                    &n_families,
                                    distances.begin(),
                                    gene_to_fam.begin(),
                                    dscale.begin(),
                                    loess_x.begin(),
                                    loess_y.begin(),
                                    indices_used.begin(),
                                    perm_tmp.begin(),
                                    stack_left_tmp.begin(),
                                    stack_right_tmp.begin(),
                                    family_distances.begin(),
                                    &ierr);

    return List::create(Named("dscale") = dscale,
                        Named("loess_x") = loess_x,
                        Named("loess_y") = loess_y,
                        Named("indices_used") = indices_used,
                        Named("perm_tmp") = perm_tmp,
                        Named("stack_left_tmp") = stack_left_tmp,
                        Named("stack_right_tmp") = stack_right_tmp,
                        Named("family_distances") = family_distances,
                        Named("ierr") = ierr);
}


//' Compute Residual Distance Index (RDI) for genes and identify outliers
//'
//'@param distances Numeric vector of gene distances
//' @param gene_to_fam Integer vector mapping genes to family indices
//' @param dscale Numeric vector of family scaling factors
//'@return List with RDI values, sorted RDI, permutation, and stack arrays
// [[Rcpp::export]]
List tox_compute_rdi_rcpp(NumericVector distances,
                          IntegerVector gene_to_fam,
                          NumericVector dscale) {

    int n_genes = distances.size();
    int n_families = dscale.size();
    NumericVector rdi(n_genes);
    NumericVector sorted_rdi(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(n_genes);
    IntegerVector stack_right(n_genes);

    compute_rdi_c(&n_genes,
                  &n_families,
                  distances.begin(),
                  gene_to_fam.begin(),
                  dscale.begin(),
                  rdi.begin(),
                  sorted_rdi.begin(),
                  perm.begin(),
                  stack_left.begin(),
                  stack_right.begin());

    return List::create(Named("rdi") = rdi,
                        Named("sorted_rdi") = sorted_rdi,
                        Named("perm") = perm,
                        Named("stack_left") = stack_left,
                        Named("stack_right") = stack_right);
}
//' Identify outliers based on RDI percentiles
//'
//' @param rdi Numeric vector of RDI values
//' @param percentile Percentile threshold
//' @return List with logical outlier vector and threshold value
// [[Rcpp::export]]
List tox_identify_outliers_rcpp(NumericVector rdi,
                                double percentile) {
    int n_genes = rdi.size();
    NumericVector sorted_rdi = clone(rdi);
    std::sort(sorted_rdi.begin(), sorted_rdi.end());
    IntegerVector is_outlier_int(n_genes);
    double threshold = 0.0;

    identify_outliers_c(&n_genes,
                        rdi.begin(),
                        sorted_rdi.begin(),
                        is_outlier_int.begin(),
                        &threshold,
                        &percentile);

    // Convert integer 0/1 flags to logical vector for R
    LogicalVector is_outlier(n_genes);
    for (int i = 0; i < n_genes; ++i) {
        is_outlier[i] = (is_outlier_int[i] != 0);
    }

    return List::create(Named("is_outlier") = is_outlier,
                        Named("threshold") = threshold);
}

//' Complete outlier detection workflow
//'
//' @param distances Numeric vector of gene distances
//' @param gene_to_fam Integer vector mapping genes to family indices
//' @param n_families Integer number of families
//' @param percentile Percentile threshold for outlier detection
//' @return List with outlier flags, and loess fit
// [[Rcpp::export]]
List tox_detect_outliers_rcpp(NumericVector distances,
                              IntegerVector gene_to_fam,
                              int n_families,
                              double percentile) {

    int n_genes = distances.size();
    NumericVector work_array(n_genes);
    IntegerVector perm(n_genes);
    IntegerVector stack_left(n_genes);
    IntegerVector stack_right(n_genes);
    IntegerVector is_outlier_int(n_genes);
    NumericVector loess_x(n_families);
    NumericVector loess_y(n_families);
    IntegerVector loess_n(n_families);
    int ierr = 0;

    detect_outliers_c(&n_genes,
                      &n_families,
                      distances.begin(),
                      gene_to_fam.begin(),
                      work_array.begin(),
                      perm.begin(),
                      stack_left.begin(),
                      stack_right.begin(),
                      is_outlier_int.begin(),
                      loess_x.begin(),
                      loess_y.begin(),
                      loess_n.begin(),
                      &percentile,
                      &ierr);

    // Convert integer flags to logical vector for R
    LogicalVector is_outlier(n_genes);
    for (int i = 0; i < n_genes; ++i) {
        is_outlier[i] = (is_outlier_int[i] != 0);
    }

    return List::create(Named("is_outlier") = is_outlier,
                        Named("loess_x") = loess_x,
                        Named("loess_y") = loess_y,
                        Named("loess_n") = loess_n,
                        Named("ierr") = ierr);
}
//' Build a k-d tree index for efficient nearest neighbor search
//'
//' @param X Numeric matrix of data points (dimensions x points)
//' @param dim_order Integer vector specifying the order of dimensions for splitting
//' @return List with k-d tree index 
// [[Rcpp::export]]
List tox_build_kd_index_rcpp(NumericMatrix X,
                             IntegerVector dim_order) {

    int d = X.nrow();
    int n = X.ncol();

    IntegerVector kd_ix(n);
    IntegerVector work(n);
    NumericVector subarray(n);
    IntegerVector perm(n);
    IntegerVector stack_left(n);
    IntegerVector stack_right(n);
    int ierr = 0;

    build_kd_index_C(X.begin(),
                     &d,
                     &n,
                     kd_ix.begin(),
                     dim_order.begin(),
                     work.begin(),
                     subarray.begin(),
                     perm.begin(),
                     stack_left.begin(),
                     stack_right.begin(),
                     &ierr);

    return List::create(Named("kd_ix") = kd_ix,
                        Named("ierr") = ierr);
}

//'Identify indices of non-zero elements in a mask vector
//'
//' @param mask Integer vector mask (0/1 values)
//' @param m_max Maximum number of indices to return
//' @return List with output indices and number of indices found
// [[Rcpp::export]]
List tox_which_rcpp(IntegerVector mask, 
                    int m_max) {

  int n = mask.size();

    IntegerVector idx_out(m_max);
    int m_out = 0;
    int ierr = 0;

    which_c(mask.begin(),
            &n,
            idx_out.begin(),
            &m_max,
            &m_out,
            &ierr);

    if (m_out > m_max) {
      m_out = m_max;
    }

    return List::create(Named("idx_out") = idx_out,
                        Named("m_out") = m_out,
                        Named("ierr") = ierr);
}

//' Perform 2D LOESS smoothing for gene family scaling
//'
//' @param n_total Total number of data points
//' @param n_target Number of target points to predict
//' @param x_ref Numeric vector of reference x values
//' @param y_ref Numeric vector of reference y values
//' @param indices_used Integer vector of indices of reference points used in the fit
//' @param n_used Number of reference points used in the fit
//' @param x_query Numeric vector of x values for target points
//' @param kernel_sigma Numeric value for kernel bandwidth
//' @param kernel_cutoff Numeric value for kernel cutoff distance
//' @return List with predicte y values
// [[Rcpp::export]]
List tox_loess_smooth_2d_rcpp(int n_total,
                             int n_target,
                             NumericVector x_ref,
                             NumericVector y_ref,
                             IntegerVector indices_used,
                             int n_used,
                             NumericVector x_query,
                             double kernel_sigma,
                             double kernel_cutoff) {

    NumericVector y_out(n_target);
    int ierr = 0;

    loess_smooth_2d_c(&n_total,
                      &n_target,
                      x_ref.begin(),
                      y_ref.begin(),
                      indices_used.begin(),
                      &n_used,
                      x_query.begin(),
                      &kernel_sigma,
                      &kernel_cutoff,
                      y_out.begin(),
                      &ierr);

    return List::create(Named("y_out") = y_out,
                        Named("ierr") = ierr);
}


// Helper function to convert string to char array with null termination
inline std::vector<char> filename_to_ascii(const std::string& filename, int& out_len) {
    // Allocate +1 for null termination like Python version
    std::vector<char> ascii(filename.size() + 1, 0);
    for (size_t i = 0; i < filename.size(); ++i) {
        ascii[i] = static_cast<char>(filename[i]);
    }
    // ascii[filename.size()] = 0;  // already initialized to 0
    out_len = static_cast<int>(ascii.size());
    return ascii;
}

//' Calculate serialization of integer arrays to binary files with specified filenames
//'
//' @param arr Integer vector (potentially with dim attribute for multi-dimensional arrays)
//' @param filename String specifying the output filename for serialization
//' @return  error code from serialization function
// [[Rcpp::export]]
int tox_serialize_int_array_rcpp(IntegerVector arr, std::string filename) {
    IntegerVector dim = arr.hasAttribute("dim")
                            ? as<IntegerVector>(arr.attr("dim"))
                            : IntegerVector::create((int)arr.size());
    int ndim = dim.size();
    std::vector<int> dims(ndim);
    for (int i = 0; i < ndim; ++i) {
        dims[i] = dim[i];
    }

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    int ierr = 0;

    serialize_int_nd_C((void*)arr.begin(),
                       dims.data(),
                       &ndim,
                       fname.data(),
                       &fn_len,
                       &ierr);

    return ierr;
}

//' Calculate serialization of real (double) arrays to binary files with specified filenames
//'
//' @param arr Numeric vector (potentially with dim attribute for multi-dimensional arrays)
//' @param filename String specifying the output filename for serialization
//' @return  error code from serialization function
// [[Rcpp::export]]
int tox_serialize_real_array_rcpp(NumericVector arr, 
                                  std::string filename) {

    IntegerVector dim = arr.hasAttribute("dim")
                            ? as<IntegerVector>(arr.attr("dim"))
                            : IntegerVector::create((int)arr.size());
    int ndim = dim.size();
    std::vector<int> dims(ndim);
    for (int i = 0; i < ndim; ++i) {
        dims[i] = dim[i];
    }

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    int ierr = 0;

    serialize_real_nd_C((void*)arr.begin(),
                        dims.data(),
                        &ndim,
                        fname.data(),
                        &fn_len,
                        &ierr);

    return ierr;
}

//' Calculate serialization of character arrays to binary files with specified filenames, including conversion to fixed-length ASCII representation
//'
//' @param carr Character vector (potentially with dim attribute for multi-dimensional arrays)
//' @param filename String specifying the output filename for serialization
//' @return  error code from serialization function
// [[Rcpp::export]]
int tox_serialize_char_array_rcpp(CharacterVector carr, 
                                  std::string filename) {

    IntegerVector dim = carr.hasAttribute("dim")
                            ? as<IntegerVector>(carr.attr("dim"))
                            : IntegerVector::create((int)carr.size());
    int ndim = dim.size();
    std::vector<int> dims(ndim);
    for (int i = 0; i < ndim; ++i) {
        dims[i] = dim[i];
    }

    int total = 1;
    for (int i = 0; i < ndim; ++i) {
        total *= dims[i];
    }

    int clen = 0;
    for (int i = 0; i < total; ++i) {
        std::string s = as<std::string>(carr[i]);
        if ((int)s.size() > clen) {
            clen = (int)s.size();
        }
    }
    if (clen == 0) {
        clen = 1;
    }

    std::vector<char> ascii_flat(total * clen, 0);
    for (int i = 0; i < total; ++i) {
        std::string s = as<std::string>(carr[i]);
        for (int j = 0; j < (int)s.size() && j < clen; ++j) {
            ascii_flat[i * clen + j] = static_cast<char>(s[j]);
        }
    }

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    int ierr = 0;

    serialize_char_nd_C(ascii_flat.data(),
                        dims.data(),
                        &ndim,
                        &clen,
                        fname.data(),
                        &fn_len,
                        &ierr);

    return ierr;
}

//' Calculate serialization of logical arrays to binary files with specified filenames, including conversion to integer representation
//'
//' @param arr Logical vector (potentially with dim attribute for multi-dimensional arrays)
//' @param filename String specifying the output filename for serialization
//' @return  error code from serialization function
// [[Rcpp::export]]
int tox_serialize_logical_array_rcpp(LogicalVector arr, 
                                     std::string filename) {

    IntegerVector dim = arr.hasAttribute("dim")
                            ? as<IntegerVector>(arr.attr("dim"))
                            : IntegerVector::create((int)arr.size());
    int ndim = dim.size();
    std::vector<int> dims(ndim);
    for (int i = 0; i < ndim; ++i) {
        dims[i] = dim[i];
    }

    // Convert logical to integer for serialization
    int total = 1;
    for (int i = 0; i < ndim; ++i) {
        total *= dims[i];
    }

    std::vector<int> int_arr(total);
    for (int i = 0; i < total; ++i) {
        int_arr[i] = arr[i] ? 1 : 0;
    }

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    int ierr = 0;

    serialize_logical_nd_C(int_arr.data(),
                           dims.data(),
                           &ndim,
                           fname.data(),
                           &fn_len,
                           &ierr);

    return ierr;
}


//' Calculate serialization of complex arrays to binary files with specified filenames, including conversion to interleaved real/imaginary representation
//'
//' @param arr Complex vector (potentially with dim attribute for multi-dimensional arrays)
//' @param filename String specifying the output filename for serialization
//' @return  error code from serialization function
// [[Rcpp::export]]
int tox_serialize_complex_array_rcpp(ComplexVector arr, 
                                     std::string filename) {

    IntegerVector dim = arr.hasAttribute("dim")
                            ? as<IntegerVector>(arr.attr("dim"))
                            : IntegerVector::create((int)arr.size());
    int ndim = dim.size();
    std::vector<int> dims(ndim);
    for (int i = 0; i < ndim; ++i) {
        dims[i] = dim[i];
    }

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    int ierr = 0;

    // Call Fortran C binding directly with complex array
    serialize_complex_nd_C((void*)arr.begin(),
                           dims.data(),
                           &ndim,
                           fname.data(),
                           &fn_len,
                           &ierr);

    return ierr;
}
//' Calculate deserialization of integer arrays from binary files with specified filenames
//'
//' @param filename String specifying the input filename for deserialization
//' @param max_dims Maximum number of dimensions to read (default 5)
//' @return List with deserialized values and dimensions
// [[Rcpp::export]]
List tox_deserialize_int_array_rcpp(std::string filename, 
                                    int max_dims = 5) {

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    std::vector<int> dims_out(max_dims);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

    get_array_metadata_C(fname.data(),
                         &fn_len,
                         dims_out.data(),
                         &max_dims,
                         &ndims,
                         &ierr,
                         &clen);

    if (ierr != 0) {
        return List::create(Named("ierr") = ierr);
    }

    int total = 1;
    for (int i = 0; i < ndims; ++i) {
        total *= dims_out[i];
    }

    IntegerVector out(total);

    deserialize_int_nd_C(out.begin(),
                         &total,
                         fname.data(),
                         &fn_len,
                         &ierr);

    return List::create(Named("values") = out,
                        Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin() + ndims),
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}


//' Calculate deserialization of real (double) arrays from binary files with specified filenames
//'
//' @param filename String specifying the input filename for deserialization
//' @param max_dims Maximum number of dimensions to read (default 5)
//' @return List with deserialized values and dimensions
// [[Rcpp::export]]
List tox_deserialize_real_array_rcpp(std::string filename, 
                                     int max_dims = 5) {

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    std::vector<int> dims_out(max_dims);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

      get_array_metadata_C(fname.data(),
                         &fn_len,
                         dims_out.data(),
                         &max_dims,
                         &ndims,
                         &ierr,
                         &clen);

    if (ierr != 0) {
        return List::create(Named("ierr") = ierr);
    }

    int total = 1;
    for (int i = 0; i < ndims; ++i) {
        total *= dims_out[i];
    }

    NumericVector out(total);

    deserialize_real_nd_C(out.begin(),
                          &total,
                          fname.data(),
                          &fn_len,
                          &ierr);

    return List::create(Named("values") = out,
                        Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin() + ndims),
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}

//' Calculate deserialization of real (double) arrays from binary files with specified filenames
//'
//' @param filename String specifying the input filename for deserialization
//' @param max_dims Maximum number of dimensions to read (default 5)
//' @return List with deserialized values and dimensions
// [[Rcpp::export]]
List tox_deserialize_char_array_rcpp(std::string filename, 
                                     int max_dims = 5) {

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    std::vector<int> dims_out(max_dims);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

      get_array_metadata_C(fname.data(),
                         &fn_len,
                         dims_out.data(),
                         &max_dims,
                         &ndims,
                         &ierr,
                         &clen);

    if (ierr != 0) {
        return List::create(Named("ierr") = ierr);
    }

    int total = 1;
    for (int i = 0; i < ndims; ++i) {
        total *= dims_out[i];
    }

    std::vector<char> ascii_out(total * clen);

    deserialize_char_nd_C(ascii_out.data(),
                          &clen,
                          &total,
                          fname.data(),
                          &fn_len,
                          &ierr);

    CharacterVector out(total);
    for (int i = 0; i < total; ++i) {
        std::string s;
        for (int j = 0; j < clen; ++j) {
            char ch = ascii_out[i * clen + j];
            if (ch == '\0') {
                break;
            }
            s.push_back(ch);
        }

        // Trim trailing whitespace
        size_t end = s.find_last_not_of(" \t\n\r\f\v");
        if (end != std::string::npos) {
            s = s.substr(0, end + 1);
        } else {
            s.clear();  // String was all whitespace
        }

        out[i] = s;
    }

    return List::create(Named("values") = out,
                        Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin() + ndims),
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}


//' Calculate deserialization of character arrays from binary files with specified filenames, including conversion from fixed-length ASCII representation
//'
//' @param filename String specifying the input filename for deserialization
//' @param max_dims Maximum number of dimensions to read (default 5)
//' @return List with deserialized values and dimensions
// [[Rcpp::export]]
List tox_deserialize_logical_array_rcpp(std::string filename, 
                                        int max_dims = 5) {

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    std::vector<int> dims_out(max_dims);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

     get_array_metadata_C(fname.data(),
                         &fn_len,
                         dims_out.data(),
                         &max_dims,
                         &ndims,
                         &ierr,
                         &clen);

    if (ierr != 0) {
        return List::create(Named("ierr") = ierr);
    }

    int total = 1;
    for (int i = 0; i < ndims; ++i) {
        total *= dims_out[i];
    }

    std::vector<int> int_out(total);

    deserialize_logical_nd_C(int_out.data(),
                             &total,
                             fname.data(),
                             &fn_len,
                             &ierr);

    // Convert integer to logical
    LogicalVector out(total);
    for (int i = 0; i < total; ++i) {
        out[i] = (int_out[i] != 0);
    }

    return List::create(Named("values") = out,
                        Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin() + ndims),
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}


//' Calculate deserialization of complex arrays from binary files with specified filenames, including conversion from interleaved real/imaginary representation
//'
//' @param filename String specifying the input filename for deserialization
//' @param max_dims Maximum number of dimensions to read (default 5)
//' @return List with deserialized values and dimensions
// [[Rcpp::export]]
List tox_deserialize_complex_array_rcpp(std::string filename, 
                                        int max_dims = 5) {

    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    std::vector<int> dims_out(max_dims);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

      get_array_metadata_C(fname.data(),
                         &fn_len,
                         dims_out.data(),
                         &max_dims,
                         &ndims,
                         &ierr,
                         &clen);

    if (ierr != 0) {
        return List::create(Named("ierr") = ierr);
    }

    int total = 1;
    for (int i = 0; i < ndims; ++i) {
        total *= dims_out[i];
    }

    ComplexVector out(total);

    deserialize_complex_nd_C((void*)out.begin(),
                             &total,
                             fname.data(),
                             &fn_len,
                             &ierr);

    return List::create(Named("values") = out,
                        Named("dims") = IntegerVector(dims_out.begin(), dims_out.begin() + ndims),
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}

//' Get metadata of arrays stored in binary files, including dimensions, number of dimensions
//'
//' @param filename String specifying the input filename for metadata retrieval
//' @param dims_out_capacity Maximum number of dimensions to retrieve (default 5)
//' @param with_clen Logical indicating whether to include character length in the output (default FALSE)
//' @return List with dimensions, number of dimensions, and optionally character length and error code
// [[Rcpp::export]]
List tox_get_array_metadata_rcpp(std::string filename,
                                int dims_out_capacity = 5,
                                bool with_clen = false) {
    int fn_len = 0;
    auto fname = filename_to_ascii(filename, fn_len);
    IntegerVector dims_res(dims_out_capacity);
    int ndims = 0;
    int ierr = 0;
    int clen = 0;

        get_array_metadata_C(fname.data(),
                             &fn_len, 
                             dims_res.begin(), 
                             &dims_out_capacity, 
                             &ndims, 
                             &ierr,
                             &clen);

    if (with_clen) {
        return List::create(Named("dims") = dims_res,
                            Named("ndim") = ndims,
                            Named("clen") = clen,
                            Named("ierr") = ierr);
    }

    return List::create(Named("dims") = dims_res,
                        Named("ndim") = ndims,
                        Named("ierr") = ierr);
}

//' Build a binary search tree index for a numeric vector, returning sorted indices and auxiliary stacks for traversal
//'
//' @param values Numeric vector to index
//' @return Integer vector of sorted indices (1-based Fortran indices preserved) for the BST index
// [[Rcpp::export]]
IntegerVector tox_build_bst_index_rcpp(NumericVector values) {

    int num_values = static_cast<int>(values.size());
    IntegerVector sorted_indices(num_values);
    IntegerVector left_stack(num_values);
    IntegerVector right_stack(num_values);
    int ierr = 0;

    build_bst_index_C(values.begin(),
                      &num_values,
                      sorted_indices.begin(),
                      left_stack.begin(),
                      right_stack.begin(),
                      &ierr);

    // Return sorted indices (1-based Fortran indices preserved)
    return sorted_indices;
}

//' Perform a range query on a binary search tree index for a numeric vector, returning indices of values within the specified bounds
//'
//' @param values Numeric vector to query
//' @param sorted_indices Integer vector of sorted indices from the BST index
//' @param lower_bound Numeric value specifying the lower bound of the range query
//' @param upper_bound Numeric value specifying the upper bound of the range query
//' @return List with output indices of values within the range and number of matches
// [[Rcpp::export]]
List tox_bst_range_query_rcpp(NumericVector values,
                          IntegerVector sorted_indices,
                          double lower_bound,
                          double upper_bound) {

    int num_values = static_cast<int>(values.size());
    IntegerVector output_indices(num_values);
    int num_matches = 0;
    int ierr = 0;

    bst_range_query_C(values.begin(),
                      sorted_indices.begin(),
                      &num_values,
                      &lower_bound,
                      &upper_bound,
                      output_indices.begin(),
                      &num_matches,
                      &ierr);

    return List::create(Named("output_indices") = output_indices,
                        Named("num_matches") = num_matches,
                        Named("ierr") = ierr);
}

//' Build a spherical k-d tree index for a numeric matrix, returning sorted indices and auxiliary stacks for traversal
//'
//' @param V Numeric matrix of data points (dimensions x points)
//' @param dim_order Integer vector specifying the order of dimensions for splitting
//' @return List with spherical k-d tree index 
// [[Rcpp::export]]
List tox_build_spherical_kd_rcpp(NumericMatrix V, 
                                 IntegerVector dim_order) {
    int d = V.nrow();
    int n = V.ncol();

    IntegerVector sphere_ix(n);
    IntegerVector work(n);
    NumericVector subarray(n);
    IntegerVector perm(n);
    IntegerVector stack_left(n);
    IntegerVector stack_right(n);
    int ierr = 0;

    build_spherical_kd_C(V.begin(),
                         &d,
                         &n,
                         sphere_ix.begin(),
                         dim_order.begin(),
                         work.begin(),
                         subarray.begin(),
                         perm.begin(),
                         stack_left.begin(),
                         stack_right.begin(),
                         &ierr);

    return List::create(Named("sphere_ix") = sphere_ix,
                        Named("ierr") = ierr);
}

//' Calculate reading of expression vectors from TSV files
//'
//' @param file_list_raw Raw matrix of filenames (one per column) for TSV files to read
//' @param gene_ids_raw Raw matrix of gene IDs to match against the TSV files
//' @param value_cols Integer vector of column indices in the TSV files that contain expression values
//' @param delimiter_raw Raw vector specifying the delimiter used in the TSV files
//' @param n_samples Integer specifying the expected number of samples (rows) in the TSV files
//' @param n_header_rows Integer specifying the number of header rows to skip in the TSV files
//' @param gene_col Integer specifying the column index in the TSV files that contains gene IDs
//' @return List with expression vector matrix and error code from reading function
// [[Rcpp::export]]
List tox_read_expression_vectors_tsv_rcpp(RawMatrix file_list_raw,
                                         RawMatrix gene_ids_raw,
                                         IntegerVector value_cols,
                                         RawVector delimiter_raw,
                                         int n_samples,
                                         int n_header_rows,
                                         int gene_col) {

    int file_list_len = file_list_raw.nrow();
    int n_files = file_list_raw.ncol();
    int gene_ids_len = gene_ids_raw.nrow();
    int n_genes = gene_ids_raw.ncol();
    int n_value_cols = value_cols.size();

    NumericMatrix expression_vectors(n_samples, n_genes);
    int ierr = 0;

    read_expression_vectors_tsv_C(reinterpret_cast<const char*>(file_list_raw.begin()),
                                  &file_list_len,
                                  &n_files,
                                  reinterpret_cast<const char*>(gene_ids_raw.begin()),
                                  &gene_ids_len,
                                  &n_genes,
                                  expression_vectors.begin(),
                                  &n_samples,
                                  &n_header_rows,
                                  &gene_col,
                                  value_cols.begin(),
                                  &n_value_cols,
                                  reinterpret_cast<const char*>(delimiter_raw.begin()),
                                  &ierr);

    return List::create(Named("expression_vectors") = expression_vectors,
                        Named("ierr") = ierr);
}

 //' Calculate reading of gene IDs from TSV files, returning a raw matrix of gene IDs and an error code
 //'
 //' @param filename_raw Raw vector specifying the filename of the TSV file to read
 //' @param n_genes Integer specifying the expected number of genes (columns) in the TSV
 //' @param gene_ids_len Integer specifying the maximum length of gene IDs (rows) to read from the TSV file
 //' @param n_header_rows Integer specifying the number of header rows to skip in the TSV file 
 //' @param gene_col Integer specifying the column index in the TSV file that contains gene IDs
 //' @return List with raw matrix of gene IDs and error code from reading function         
// [[Rcpp::export]]
List tox_read_gene_ids_from_tsv_file_rcpp(RawVector filename_raw,
                                         int n_genes,
                                         int gene_ids_len,
                                         int n_header_rows,
                                         int gene_col) {

    RawMatrix gene_ids_raw(gene_ids_len, n_genes);
    std::fill(gene_ids_raw.begin(), gene_ids_raw.end(), static_cast<Rbyte>(0));

    int ierr = 0;
    int fn_len = filename_raw.size();

    read_gene_ids_from_tsv_file_C(reinterpret_cast<const char*>(filename_raw.begin()),
                                  &fn_len,
                                  reinterpret_cast<char*>(gene_ids_raw.begin()),
                                  &gene_ids_len,
                                  &n_genes,
                                  &n_header_rows,
                                  &gene_col,
                                  &ierr);

    return List::create(Named("gene_ids_raw") = gene_ids_raw,
                        Named("ierr") = ierr);
}

//' Calculate reading of OrthoFinder output files to extract family IDs and gene-to-family mappings
//'
//' @param filename_raw Raw vector specifying the filename of the OrthoFinder output file to read
//' @param gene_ids_raw Raw matrix of gene IDs to match against the OrthoFinder output file
//' @param n_families Integer specifying the expected number of families (columns) in the OrthoFinder output file
//' @param family_ids_len Integer specifying the maximum length of family IDs (rows) to read from the OrthoFinder output file
//' @return List with raw matrix of family IDs, integer vector mapping genes to families, and error code from reading function
// [[Rcpp::export]]
List tox_read_orthofinder_file_rcpp(RawVector filename_raw,
                                   RawMatrix gene_ids_raw,
                                   int n_families,
                                   int family_ids_len) {

    int gene_ids_len = gene_ids_raw.nrow();
    int n_genes = gene_ids_raw.ncol();

    RawMatrix family_ids_raw(family_ids_len, n_families);
    std::fill(family_ids_raw.begin(), family_ids_raw.end(), static_cast<Rbyte>(0));
    IntegerVector gene_to_fam(n_genes);
    int ierr = 0;
    int fn_len = filename_raw.size();

    read_orthofinder_file_C(reinterpret_cast<const char*>(filename_raw.begin()),
                            &fn_len,
                            reinterpret_cast<const char*>(gene_ids_raw.begin()),
                            &gene_ids_len,
                            &n_genes,
                            reinterpret_cast<char*>(family_ids_raw.begin()),
                            &family_ids_len,
                            &n_families,
                            gene_to_fam.begin(),
                            &ierr);

    return List::create(Named("family_ids_raw") = family_ids_raw,
                        Named("gene_to_fam") = gene_to_fam,
                        Named("ierr") = ierr);
}

//' Calculate validation of data structure for gene IDs, family IDs, gene-to-family mappings, expression vectors, family centroids, and shift vectors
//'
//' @param gene_ids_raw Raw matrix of gene IDs to validate
//' @param gene_family_ids_raw Raw matrix of family IDs to validate
//' @param gene_to_fam Integer vector mapping genes to families to validate
//' @param expression_vectors Numeric matrix of expression vectors to validate
//' @param family_centroids Numeric matrix of family centroids to validate
//' @param shift_vectors Numeric matrix of shift vectors to validate
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_data_structure_rcpp(RawMatrix gene_ids_raw,
                                     RawMatrix gene_family_ids_raw,
                                     IntegerVector gene_to_fam,
                                     NumericMatrix expression_vectors,
                                     NumericMatrix family_centroids,
                                     NumericMatrix shift_vectors) {

    int n_genes = gene_ids_raw.ncol();
    int gene_ids_len = gene_ids_raw.nrow();
    int n_families = gene_family_ids_raw.ncol();
    int fam_len = gene_family_ids_raw.nrow();
    int n_samples = expression_vectors.nrow();
    int ierr = 0;

    validate_data_structure_C(&n_genes,
                             &n_families,
                             &n_samples,
                             reinterpret_cast<const char*>(gene_ids_raw.begin()),
                             &gene_ids_len,
                             reinterpret_cast<const char*>(gene_family_ids_raw.begin()),
                             &fam_len,
                             gene_to_fam.begin(),
                             expression_vectors.begin(),
                             family_centroids.begin(),
                             shift_vectors.begin(),
                             &ierr);

    return ierr;
}

//' Calculate filtering of unassigned genes based on gene IDs and gene-to-family mappings
//'
//' @param gene_ids_raw Raw matrix of gene IDs to filter
//' @param gene_to_fam Integer vector mapping genes to families to use for filtering
//' @return List with integer mask indicating which genes are assigned to families, number of genes kept after filtering, and error code from filtering function
// [[Rcpp::export]]

List tox_filter_unassigned_genes_rcpp(RawMatrix gene_ids_raw,
                                      IntegerVector gene_to_fam) {

    int gene_ids_len = gene_ids_raw.nrow();
    int n_genes = gene_ids_raw.ncol();
    IntegerVector mask(n_genes);
    int n_genes_kept = 0;
    int ierr = 0;

    filter_unassigned_genes_C(reinterpret_cast<const char*>(gene_ids_raw.begin()),
                              &gene_ids_len,
                              &n_genes,
                              gene_to_fam.begin(),
                              mask.begin(),
                              &n_genes_kept,
                              &ierr);

    return List::create(Named("mask") = mask,
                        Named("n_genes_kept") = n_genes_kept,
                        Named("ierr") = ierr);
}
//' Calculate validation of gene-to-family mapping based on integer vector mapping genes to families and the number of families
//'
//' @param gene_to_fam Integer vector mapping genes to families to validate
//' @param n_families Integer specifying the expected number of families for validation
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_gene_to_family_mapping_rcpp(IntegerVector gene_to_fam,
                                             int n_families) {
    int n_genes = gene_to_fam.size();
    int ierr = 0;

    validate_gene_to_family_mapping_C(gene_to_fam.begin(),
                                     &n_genes,
                                     &n_families,
                                     &ierr);

    return ierr;
}
//' Calculate validation of expression data based on expression vector matrix and checks for NaN, infinite, or optionally negative values
//'
//' @param expression_vectors Numeric matrix of expression vectors to validate
//' @param check_non_negative Logical indicating whether to check for negative values in the expression data (default FALSE)
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_expression_data_rcpp(NumericMatrix expression_vectors,
                                      bool check_non_negative) {

    int n_samples = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int flag = check_non_negative ? 1 : 0;
    int ierr = 0;

    validate_expression_data_C(expression_vectors.begin(),
                              &n_genes,
                              &n_samples,
                              &flag,
                              &ierr);

    return ierr;
}
//' Calculate validation of family centroids based on family centroid matrix and checks for NaN or infinite values
//'
//' @param family_centroids Numeric matrix of family centroids to validate
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_family_centroids_rcpp(NumericMatrix family_centroids) {

    int n_samples = family_centroids.nrow();
    int n_families = family_centroids.ncol();
    int ierr = 0;

    validate_family_centroids_C(family_centroids.begin(),
                                &n_families,
                                &n_samples,
                                &ierr);

    return ierr;
}
//' Calculate validation of shift vectors based on shift vector matrix, expression vectors, family centroids, and gene-to-family mappings
//'
//' @param shift_vectors Numeric matrix of shift vectors to validate
//' @param expression_vectors Numeric matrix of expression vectors to use for validation
//' @param family_centroids Numeric matrix of family centroids to use for validation
//' @param gene_to_fam Integer vector mapping genes to families to use for validation
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_shift_vectors_rcpp(NumericMatrix shift_vectors,
                                    NumericMatrix expression_vectors,
                                    NumericMatrix family_centroids,
                                    IntegerVector gene_to_fam) {

    int n_samples = expression_vectors.nrow();
    int n_genes = expression_vectors.ncol();
    int n_families = family_centroids.ncol();
    int ierr = 0;

    validate_shift_vectors_C(shift_vectors.begin(),
                             expression_vectors.begin(),
                             family_centroids.begin(),
                             gene_to_fam.begin(),
                             &n_genes,
                             &n_samples,
                             &n_families,
                             &ierr);

    return ierr;
}

//' Calculate validation of string array uniqueness based on raw matrix of strings, checking for duplicates
//'
//' @param string_arr_raw Raw matrix of strings to validate for uniqueness
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_string_array_uniqueness_rcpp(RawMatrix string_arr_raw) {

    int str_len = string_arr_raw.nrow();
    int n_strings = string_arr_raw.ncol();
    int ierr = 0;

    validate_string_array_uniqueness_C(reinterpret_cast<const char*>(string_arr_raw.begin()),
                                       &str_len,
                                       &n_strings,
                                       &ierr);

    return ierr;
}

//' Calculate validation of all data components together, including gene IDs, family IDs, gene-to-family mappings, expression vectors, family centroids, and shift vectors
//'
//' @param gene_ids_raw Raw matrix of gene IDs to validate
//' @param gene_family_ids_raw Raw matrix of family IDs to validate
//' @param gene_to_fam Integer vector mapping genes to families to validate
//' @param expression_vectors Numeric matrix of expression vectors to validate
//' @param family_centroids Numeric matrix of family centroids to validate
//' @param shift_vectors Numeric matrix of shift vectors to validate
//' @return  error code from validation function
// [[Rcpp::export]]
int tox_validate_all_data_rcpp(RawMatrix gene_ids_raw,
                               RawMatrix gene_family_ids_raw,
                               IntegerVector gene_to_fam,
                               NumericMatrix expression_vectors,
                               NumericMatrix family_centroids,
                               NumericMatrix shift_vectors) {

    int n_genes = gene_ids_raw.ncol();
    int gene_len = gene_ids_raw.nrow();
    int n_families = gene_family_ids_raw.ncol();
    int fam_len = gene_family_ids_raw.nrow();
    int n_samples = expression_vectors.nrow();
    int ierr = 0;

    validate_all_data_C(&n_genes,
                        &n_families,
                        &n_samples,
                        reinterpret_cast<const char*>(gene_ids_raw.begin()),
                        &gene_len,
                        reinterpret_cast<const char*>(gene_family_ids_raw.begin()),
                        &fam_len,
                        gene_to_fam.begin(),
                        expression_vectors.begin(),
                        family_centroids.begin(),
                        shift_vectors.begin(),
                        &ierr);

    return ierr;
}


//' Calculate projection of expression vectors onto family centroids using the RAP method, with selection masks for vectors and axes
//'
//' @param vecs Numeric matrix of expression vectors (dimensions x vectors)
//' @param vecs_selection_mask Integer vector indicating which vectors to include in the projection (1 for include, 0 for exclude)
//' @param axes_selection_mask Integer vector indicating which dimensions (axes) to include in the projection
//' @return List with matrix of projected values and error code from projection function
// [[Rcpp::export]]
List tox_omics_vector_RAP_projection_rcpp(NumericMatrix vecs,
                                         IntegerVector vecs_selection_mask,
                                         IntegerVector axes_selection_mask) {
    int n_axes = vecs.nrow();
    int n_vecs = vecs.ncol();

    int n_selected_vecs = 0;
    for (int j = 0; j < n_vecs; ++j) if (vecs_selection_mask[j] == 1) ++n_selected_vecs;

    int n_selected_axes = 0;
    for (int i = 0; i < n_axes; ++i) if (axes_selection_mask[i] == 1) ++n_selected_axes;

    NumericMatrix projections(n_selected_axes, n_selected_vecs);
    int ierr = 0;

    omics_vector_RAP_projection_c(vecs.begin(),
                                 &n_axes,
                                 &n_vecs,
                                 vecs_selection_mask.begin(),
                                 &n_selected_vecs,
                                 axes_selection_mask.begin(),
                                 &n_selected_axes,
                                 projections.begin(),
                                 &ierr);

    return List::create(Named("projections") = projections,
                        Named("ierr") = ierr);
}

//' Calculate projection of expression vectors onto family centroids using RAP method, with optimized C implementation for field-based projections
//'
//' @param vecs Numeric matrix of expression vectors (dimensions x vectors)
//' @param vecs_selection_mask Integer vector indicating which vectors to include in the projection (1 for include, 0 for exclude)
//' @param axes_selection_mask Integer vector indicating which dimensions (axes) to include in the projection
//' @return List with matrix of projected values and error code from projection function
// [[Rcpp::export]]
List tox_omics_field_RAP_projection_rcpp(NumericMatrix vecs,
                                        IntegerVector vecs_selection_mask,
                                        IntegerVector axes_selection_mask) {
    int n_rows = vecs.nrow();
    int n_vecs = vecs.ncol();
    int n_axes = n_rows / 2;

    int n_selected_vecs = 0;
    for (int j = 0; j < n_vecs; ++j) if (vecs_selection_mask[j] == 1) ++n_selected_vecs;

    int n_selected_axes = 0;
    for (int i = 0; i < n_axes; ++i) if (axes_selection_mask[i] == 1) ++n_selected_axes;

    NumericMatrix projections(n_selected_axes, n_selected_vecs);
    int ierr = 0;

    omics_field_RAP_projection_c(vecs.begin(),
                                &n_axes,
                                &n_vecs,
                                vecs_selection_mask.begin(),
                                &n_selected_vecs,
                                axes_selection_mask.begin(),
                                &n_selected_axes,
                                projections.begin(),
                                &ierr);

    return List::create(Named("projections") = projections,
                        Named("ierr") = ierr);
}

//' Calculate contributions of factors to dependent variables based on specified mode (e.g., raw, min, mean) using optimized C implementation
//'
//' @param factor Numeric vector of factor values
//' @param dependent Numeric vector of dependent variable values
//' @param mode Integer specifying the mode of contribution calculation (e.g., 0 for raw, 1 for min, 2 for mean)
//' @return List with local contributions for each dimension, total contribution, and error code from contribution calculation function
// [[Rcpp::export]]
List tox_compute_contributions_rcpp(NumericVector factor,
                                    NumericVector dependent,
                                    std::string mode) {

  int n_dims = factor.size();
  NumericVector local_contributions(n_dims);
  double total_contribution = 0.0;
  int ierr = 0;

  
  char mode_c[8] = {' '};
  size_t len = std::min(mode.size(), size_t(8));
  std::memcpy(mode_c, mode.c_str(), len);

  compute_contributions_c(factor.begin(),
                         dependent.begin(),
                         &n_dims,
                         mode_c,
                         local_contributions.begin(),
                         &total_contribution,
                         &ierr);

  return List::create(Named("local_contributions") = local_contributions,
                      Named("total_contribution") = total_contribution,
                      Named("ierr") = ierr);
}
//' Calculate contributions of multiple factors to multiple dependent variables across samples and timepoints based on specified mode (e.g., raw, min, mean) using optimized C implementation
//'
//' @param trajectories Numeric vector containing factor and dependent variable trajectories (dimensions: n_factors x n_samples x n_timepoints)
//' @param n_factors Integer specifying the number of factors in the trajectories
//' @param n_samples Integer specifying the number of samples in the trajectories
//' @param n_timepoints Integer specifying the number of timepoints in the trajectories
//' @param factor_indices Integer vector specifying the indices of the factors to include in the contribution calculations
//' @param n_selected_factors Integer specifying the number of selected factors to include in the contribution
//' @param dependent_indices Integer vector specifying the indices of the dependent variables to include in the contribution calculations
//' @param n_selected_dependents Integer specifying the number of selected dependent variables to include in the contribution calculations
//' @param mode String specifying the mode of contribution calculation (e.g., "raw", "min", "mean")
//' @return List with local contributions for each factor-dependent pair across samples and timepoints and  total contributions for each factor-dependent pair across samples
// [[Rcpp::export]]
List tox_compute_all_contributions_rcpp(NumericVector trajectories,
                                        int n_factors,
                                        int n_samples,
                                        int n_timepoints,
                                        IntegerVector factor_indices,
                                        int n_selected_factors,
                                        IntegerVector dependent_indices,
                                        int n_selected_dependents,
                                        std::string mode) {

  // Dimensions for output arrays
  int n_tp = n_timepoints;
  // local_contributions: (n_selected_factors x n_selected_dependents x n_samples x n_timepoints)
  NumericVector local_contributions(n_selected_factors * n_selected_dependents * n_samples * n_timepoints);
  // total_contributions: (n_selected_factors x n_selected_dependents x n_samples)
  NumericVector total_contributions(n_selected_factors * n_selected_dependents * n_samples);
  // temp arrays for Fortran workspace
  NumericVector temp_factors(n_timepoints);
  NumericVector temp_dependent(n_timepoints);
  int ierr = 0;

  // Fortran expects mode as a char array (e.g., 'raw', 'min', 'mean'), pad/truncate to 8 chars
  char mode_c[8] = {' '};
  size_t len = std::min(mode.size(), size_t(8));
  std::memcpy(mode_c, mode.c_str(), len);

  compute_all_contributions_c(trajectories.begin(),
                              &n_factors,
                              &n_samples,
                              &n_timepoints,
                              factor_indices.begin(),
                              &n_selected_factors,
                              dependent_indices.begin(),
                              &n_selected_dependents,
                              mode_c,
                              local_contributions.begin(),
                              total_contributions.begin(),
                              temp_factors.begin(),
                              temp_dependent.begin(),
                              &ierr);

  return List::create(Named("local_contributions") = local_contributions,
                      Named("total_contributions") = total_contributions,
                      Named("ierr") = ierr);
}
//' Compute permutation test contributions of factors to dependent variables based on specified mode.
//'
//' @param trajectories Numeric vector containing factor and dependent variable trajectories (dimensions: n_factors x n_samples x n_timepoints)
//' @param n_factors Integer specifying the number of factors in the trajectories
//' @param n_samples Integer specifying the number of samples in the trajectories
//' @param n_timepoints Integer specifying the number of timepoints in the trajectories
//' @param factor_idx Integer specifying the index of the factor to test contributions for
//' @param dependent_idx Integer specifying the index of the dependent variable to test contributions for 
//' @param sample_idx Integer specifying the index of the sample to test contributions for
//' @param mode String specifying the baseline mode ('raw', 'min', 'mean')
//' @param n_permutations Integer specifying the number of permutations to perform
//' @param random_seed Integer seed for random number generation in permutation test
//' @return List with local contributions for each timepoint across permutations, total contributions across permutations, and error code from permutation test function
// [[Rcpp::export]]
List tox_perform_permutation_test_rcpp(NumericVector trajectories,
                                       int n_factors,
                                       int n_samples,
                                       int n_timepoints,
                                       int factor_idx,
                                       int dependent_idx,
                                       int sample_idx,
                                       std::string mode,
                                       int n_permutations,
                                       int random_seed) {
 
                                        
  char mode_c[8] = {' '};
  size_t len = std::min(mode.size(), size_t(8));
  std::memcpy(mode_c, mode.c_str(), len);

  NumericMatrix local_contributions(n_timepoints, n_permutations);
  NumericVector total_contributions(n_permutations);
  NumericVector temp_factor(n_timepoints);
  NumericVector temp_dependent(n_timepoints);
  int ierr = 0;

  perform_permutation_test_c(trajectories.begin(),
                             &n_factors,
                             &n_samples,
                             &n_timepoints,
                             &factor_idx,
                             &dependent_idx,
                             &sample_idx,
                             mode_c,
                             &n_permutations,
                             local_contributions.begin(),
                             total_contributions.begin(),
                             temp_factor.begin(),
                             temp_dependent.begin(),
                             &ierr,
                             &random_seed);

  return List::create(Named("local_contributions") = local_contributions,
                      Named("total_contributions") = total_contributions,
                      Named("ierr") = ierr);
}


//' compute p-values for observed contributions against permutation distributions, with optimized C implementation
//'
//' @param local_contributions_observed Numeric vector of local contributions observed for each timepoint
//' @param total_contribution_observed Numeric value of total contribution observed across timepoints
//' @param local_contributions_perm Numeric matrix of local contributions from permutations (dimensions: n_timepoints x n_permutations)
//' @param total_contributions_perm Numeric vector of total contributions from permutations (length: n_permutations)
//' @param n_timepoints Integer specifying the number of timepoints for local contributions
//' @param n_permutations Integer specifying the number of permutations for the null distribution
//' @return List with local p-values for each timepoint, total p-value, and error
// [[Rcpp::export]]
List tox_compute_p_values_rcpp(NumericVector local_contributions_observed,
                              double total_contribution_observed,
                              NumericMatrix local_contributions_perm,
                              NumericVector total_contributions_perm,
                              int n_timepoints,
                              int n_permutations) {

  NumericVector local_p_values(n_timepoints);
  double total_p_value = 0.0;
  int ierr = 0;

  compute_p_values_c(local_contributions_observed.begin(),
                     &total_contribution_observed,
                     local_contributions_perm.begin(),
                     total_contributions_perm.begin(),
                     &n_timepoints,
                     &n_permutations,
                     local_p_values.begin(),
                     &total_p_value,
                     &ierr);

  return List::create(Named("local_p_values") = local_p_values,
                      Named("total_p_value") = total_p_value,
                      Named("ierr") = ierr);
}


//' Compute velocity for a single trajectory
//'
//' @param trajectory Numeric vector representing a single trajectory (length = n_timepoints)
//' @param n_timepoints Integer number of timepoints
//' @return Numeric vector of velocity values (length = n_timepoints - 1)
// [[Rcpp::export]]
NumericVector tox_compute_velocity_trajectory_rcpp(NumericVector trajectory,
                                                   int n_timepoints) {
  NumericVector velocity(n_timepoints - 1);
  int ierr = 0;
  compute_velocity_trajectory_c(trajectory.begin(),
                                &n_timepoints,
                                velocity.begin(),
                                &ierr);
  return velocity;
}

//' Compute acceleration from velocity for a single trajectory
//'
//' @param velocity Numeric vector of velocity values (length = n_timepoints - 1)
//' @param n_timepoints Integer number of timepoints
//' @return Numeric vector of acceleration values (length = n_timepoints - 2)
// [[Rcpp::export]]
NumericVector tox_compute_acceleration_from_velocity_trajectory_rcpp(NumericVector velocity,
                                                           int n_timepoints) {
  NumericVector acceleration(n_timepoints - 2);
  int ierr = 0;

  compute_acceleration_from_velocity_trajectory_c(velocity.begin(),
                                                  &n_timepoints,
                                                  acceleration.begin(),
                                                  &ierr);
  return acceleration;
}

//' Compute velocity for all trajectories 
//'
//' @param trajectories Numeric vector of flattened trajectories (factors x samples x timepoints)
//' @param n_factors Integer number of factors
//' @param n_samples Integer number of samples
//' @param n_timepoints Integer number of timepoints
//' @return List with numeric vector of velocity values (length = (n_timepoints - 1) x n_factors x n_samples) and error code from velocity calculation function
// [[Rcpp::export]]
List tox_compute_velocity_trajectories_rcpp(NumericVector trajectories,
                                            int n_factors,
                                            int n_samples,
                                            int n_timepoints) {

  NumericVector velocity((n_timepoints - 1) * n_factors * n_samples);
  int ierr = 0;

  compute_velocity_trajectories_c(trajectories.begin(),
                                  &n_factors,
                                  &n_samples,
                                  &n_timepoints,
                                  velocity.begin(),
                                  &ierr);
velocity.attr("dim") = IntegerVector::create(n_timepoints - 1, n_factors, n_samples);
return List::create(Named("velocity") = velocity,
                    Named("ierr") = ierr);
}

//' Compute acceleration from velocity for all trajectories (batch mode)
//'
//' @param velocity Numeric vector of velocity values (length = (n_timepoints - 1) x n_factors x n_samples)
//' @param n_factors Integer number of factors
//' @param n_samples Integer number of samples
//' @param n_timepoints Integer number of timepoints
//' @return List with numeric vector of acceleration values (length = (n_timepoints - 2) x n_factors x n_samples) and error code from acceleration calculation function
// [[Rcpp::export]]
List tox_compute_acceleration_from_velocity_rcpp(NumericVector velocity,
                                                 int n_factors,
                                                 int n_samples,
                                                 int n_timepoints) {

  NumericVector acceleration((n_timepoints - 2) * n_factors * n_samples);
  int ierr = 0;

  compute_acceleration_from_velocity_c(velocity.begin(),
                                       &n_factors,
                                       &n_samples,
                                       &n_timepoints,
                                       acceleration.begin(),
                                       &ierr);
  acceleration.attr("dim") = IntegerVector::create(n_timepoints - 2, n_factors, n_samples);
  return List::create(Named("acceleration") = acceleration,
                      Named("ierr") = ierr);
}

//' Compute velocity and acceleration contributions for all trajectories
//'
//' @param trajectories Numeric vector of flattened trajectories (factors x samples x timepoints)
//' @param n_factors Integer number of factors
//' @param n_samples Integer number of samples
//' @param n_timepoints Integer number of timepoints
//' @param mode String indicating calculation mode (e.g., 'raw', 'min', 'mean')
//' @return List with velocity and acceleration contributions, series, and error code
// [[Rcpp::export]]
List tox_compute_velocity_acceleration_contributions_rcpp(NumericVector trajectories,
                                                          int n_factors,
                                                          int n_samples,
                                                          int n_timepoints,
                                                          std::string mode) {

  // Allocate workspace arrays
  NumericVector factor_workspace((n_timepoints - 1) * n_factors);
  NumericVector dependent_workspace(n_timepoints - 1);
  NumericVector contributions_workspace(n_timepoints - 1);

  // Allocate output arrays with correct shapes
  NumericVector contrib_velocity(n_factors * n_factors * n_samples);
  NumericVector velocity_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
  NumericVector contrib_acceleration(n_factors * n_factors * n_samples);
  NumericVector acceleration_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
  int ierr = 0;
  char mode_c[8] = {' '};
  size_t len = std::min(mode.size(), size_t(8));
  std::memcpy(mode_c, mode.c_str(), len);

  compute_velocity_acceleration_contributions_c(
    trajectories.begin(),
    &n_factors,
    &n_samples,
    &n_timepoints,
    mode_c,
    factor_workspace.begin(),
    dependent_workspace.begin(),
    contributions_workspace.begin(),
    contrib_velocity.begin(),
    velocity_contribution_series.begin(),
    contrib_acceleration.begin(),
    acceleration_contribution_series.begin(),
    &ierr
  );

  // Set dimensions for output arrays
  contrib_velocity.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
  velocity_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);
  contrib_acceleration.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
  acceleration_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);

  return List::create(Named("contrib_velocity") = contrib_velocity,
                      Named("velocity_contribution_series") = velocity_contribution_series,
                      Named("contrib_acceleration") = contrib_acceleration,
                      Named("acceleration_contribution_series") = acceleration_contribution_series,
                      Named("ierr") = ierr
  );
}

//' Compute velocity and acceleration contributions for all trajectories (alloc version)
//'
//' @param trajectories Numeric vector of flattened trajectories (factors x samples x timepoints)
//' @param n_factors Integer number of factors
//' @param n_samples Integer number of samples
//' @param n_timepoints Integer number of timepoints
//' @param mode String indicating calculation mode (e.g., 'raw', 'min', 'mean')
//' @return List with velocity and acceleration contributions, series, and error code
// [[Rcpp::export]]
List tox_compute_velocity_acceleration_contributions_alloc_rcpp(NumericVector trajectories,
                                                               int n_factors,
                                                               int n_samples,
                                                               int n_timepoints,
                                                               std::string mode) {

  NumericVector contrib_velocity(n_factors * n_factors * n_samples);
  NumericVector velocity_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
  NumericVector contrib_acceleration(n_factors * n_factors * n_samples);
  NumericVector acceleration_contribution_series(n_timepoints * n_factors * n_factors * n_samples);
  int ierr = 0;
  char mode_c[8] = {' '};
  size_t len = std::min(mode.size(), size_t(8));
  std::memcpy(mode_c, mode.c_str(), len);

  compute_velocity_acceleration_contributions_alloc_c(trajectories.begin(),
                                                     &n_factors,
                                                     &n_samples,
                                                     &n_timepoints,
                                                     mode_c,
                                                     contrib_velocity.begin(),
                                                     velocity_contribution_series.begin(),
                                                     contrib_acceleration.begin(),
                                                     acceleration_contribution_series.begin(),
                                                     &ierr);
  // Set dimensions for output arrays
  contrib_velocity.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
  velocity_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);
  contrib_acceleration.attr("dim") = IntegerVector::create(n_factors, n_factors, n_samples);
  acceleration_contribution_series.attr("dim") = IntegerVector::create(n_timepoints, n_factors, n_factors, n_samples);

  return List::create(Named("contrib_velocity") = contrib_velocity,
                      Named("velocity_contribution_series") = velocity_contribution_series,
                      Named("contrib_acceleration") = contrib_acceleration,
                      Named("acceleration_contribution_series") = acceleration_contribution_series,
                      Named("ierr") = ierr);
}

//' Compute empirical distribution function (EDF) for a numeric vector of values
//'
//' @param values Numeric vector of values to compute EDF for
//' @return List with unique values, CDF values, number of unique values, and error code from EDF computation function
// [[Rcpp::export]]
List tox_compute_edf_rcpp(NumericVector values) {

  int n_values = values.size();
  NumericVector unique_values(n_values);
  NumericVector cdf_values(n_values);
  int n_unique = 0;
  int ierr = 0;

  compute_edf_c(values.begin(),
                &n_values,
                unique_values.begin(),
                cdf_values.begin(),
                &n_unique,
                &ierr);

     return List::create(Named("unique_values") = unique_values,
                         Named("cdf_values") = cdf_values,
                         Named("n_unique") = n_unique,
                         Named("ierr") = ierr);
}

//' Compute empirical distribution function (EDF) with a pre-sorted permutation
//'
//' @param values Numeric vector of values to compute EDF for
//' @param perm Integer vector of permutation indices sorted by values[perm]
//' @return List with unique values, CDF values, number of unique values, and error code from EDF expert computation
// [[Rcpp::export]]
List tox_compute_edf_expert_rcpp(NumericVector values,
                                 IntegerVector perm) {

  int n_values = values.size();
  NumericVector unique_values(n_values);
  NumericVector cdf_values(n_values);
  int n_unique = 0;
  int ierr = 0;

  compute_edf_expert_c(values.begin(),
                       &n_values,
                       perm.begin(),
                       unique_values.begin(),
                       cdf_values.begin(),
                       &n_unique,
                       &ierr);
 
     return List::create(Named("unique_values") = unique_values,
                         Named("cdf_values") = cdf_values,
                         Named("n_unique") = n_unique,
                         Named("ierr") = ierr);
}