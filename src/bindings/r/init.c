// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
// weak so the one library still loads into a non-R host (Python/ctypes);
// see helper/codegen/emit/c_call.py
#pragma weak COMPLEX
#pragma weak INTEGER
#pragma weak LENGTH
#pragma weak LOGICAL
#pragma weak REAL
#pragma weak STRING_ELT
#pragma weak TYPEOF
#pragma weak XLENGTH
#pragma weak SET_STRING_ELT
#pragma weak SET_VECTOR_ELT
#pragma weak R_CHAR
#pragma weak R_alloc
#pragma weak R_DimSymbol
#pragma weak R_NamesSymbol
#pragma weak R_NilValue
#pragma weak R_NaString
#pragma weak R_registerRoutines
#pragma weak R_useDynamicSymbols
#pragma weak Rf_allocVector
#pragma weak Rf_asInteger
#pragma weak Rf_asLogical
#pragma weak Rf_asReal
#pragma weak Rf_coerceVector
#pragma weak Rf_duplicate
#pragma weak Rf_getAttrib
#pragma weak Rf_length
#pragma weak Rf_mkChar
#pragma weak Rf_mkCharLen
#pragma weak Rf_protect
#pragma weak Rf_setAttrib
#pragma weak Rf_unprotect
#pragma weak Rf_ScalarComplex
#pragma weak Rf_ScalarInteger
#pragma weak Rf_ScalarLogical
#pragma weak Rf_ScalarReal

// forward declarations of the .Call entry points
SEXP build_bst_index_call(SEXP);
SEXP bst_range_query_call(SEXP, SEXP, SEXP, SEXP);
SEXP build_kd_index_call(SEXP, SEXP);
SEXP build_spherical_kd_call(SEXP, SEXP);
SEXP deserialize_char_helper_call(SEXP, SEXP, SEXP);
SEXP deserialize_complex_helper_call(SEXP, SEXP);
SEXP deserialize_int_helper_call(SEXP, SEXP);
SEXP deserialize_logical_helper_call(SEXP, SEXP);
SEXP deserialize_real_helper_call(SEXP, SEXP);
SEXP serialize_char_helper_call(SEXP, SEXP);
SEXP serialize_complex_helper_call(SEXP, SEXP);
SEXP serialize_int_helper_call(SEXP, SEXP);
SEXP serialize_logical_helper_call(SEXP, SEXP);
SEXP serialize_real_helper_call(SEXP, SEXP);
SEXP get_array_metadata_call(SEXP, SEXP);
SEXP loess_smooth_2d_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP compute_edf_expert_call(SEXP, SEXP);
SEXP compute_edf_call(SEXP);
SEXP compute_scaled_distance_quantile_call(SEXP, SEXP, SEXP, SEXP);
SEXP cluster_factor_trajectories_k_means_call(SEXP, SEXP, SEXP);
SEXP k_means_clustering_call(SEXP, SEXP, SEXP);
SEXP linkage_clustering_call(SEXP, SEXP);
SEXP create_zip_archive_call(SEXP, SEXP, SEXP);
SEXP save_tox_data_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP get_tox_data_dims_call(SEXP);
SEXP read_tox_data_into_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP compute_gene_means_call(SEXP);
SEXP compute_residuals_call(SEXP, SEXP);
SEXP pool_means_call(SEXP, SEXP, SEXP);
SEXP pool_means_expert_call(SEXP, SEXP, SEXP);
SEXP calc_neighborhood_size_call(SEXP, SEXP, SEXP, SEXP);
SEXP construct_neighborhoods_call(SEXP, SEXP, SEXP, SEXP);
SEXP gjct_permutation_test_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP gjct_permutation_test_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP determine_shared_residual_range_expert_call(SEXP, SEXP, SEXP);
SEXP determine_shared_residual_range_call(SEXP, SEXP, SEXP);
SEXP build_residual_histograms_call(SEXP, SEXP, SEXP, SEXP);
SEXP compute_divergence_per_reference_point_call(SEXP, SEXP);
SEXP compute_weighted_global_divergence_call(SEXP, SEXP, SEXP);
SEXP fjct_compute_jsd_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP fjct_compute_jsd_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP fjct_compute_contribution_scores_call(SEXP, SEXP);
SEXP read_expression_vectors_tsv_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP read_gene_ids_from_tsv_file_call(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP read_orthofinder_file_call(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP get_unassigned_mask_call(SEXP);
SEXP validate_data_structure_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP validate_gene_to_family_mapping_call(SEXP, SEXP);
SEXP validate_expression_data_call(SEXP, SEXP);
SEXP validate_family_centroids_call(SEXP);
SEXP validate_shift_vectors_call(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP validate_string_array_uniqueness_call(SEXP);
SEXP validate_all_data_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP euclidean_distance_call(SEXP, SEXP);
SEXP distance_to_centroid_call(SEXP, SEXP, SEXP);
SEXP mean_vector_call(SEXP, SEXP);
SEXP group_centroid_orthologs_expert_call(SEXP, SEXP, SEXP, SEXP);
SEXP group_centroid_orthologs_call(SEXP, SEXP, SEXP, SEXP);
SEXP group_centroid_all_expert_call(SEXP, SEXP, SEXP);
SEXP group_centroid_all_call(SEXP, SEXP, SEXP);
SEXP compute_family_scaling_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP compute_family_scaling_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP compute_rdi_expert_call(SEXP, SEXP, SEXP);
SEXP compute_rdi_call(SEXP, SEXP, SEXP);
SEXP identify_outliers_call(SEXP, SEXP, SEXP, SEXP);
SEXP detect_outliers_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP detect_outliers_call(SEXP, SEXP, SEXP, SEXP);
SEXP loess_fit_plain_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP loess_fit_plain_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP loess_fit_robust_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP loess_fit_robust_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP tox_loess_required_workspace_call(SEXP, SEXP, SEXP);
SEXP normalize_unit_length_call(SEXP);
SEXP normalization_pipeline_call(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP normalize_by_std_dev_call(SEXP, SEXP, SEXP);
SEXP root_mean_sq_normalization_call(SEXP);
SEXP quantile_normalization_expert_call(SEXP);
SEXP quantile_normalization_call(SEXP);
SEXP log2_transformation_call(SEXP);
SEXP calc_tiss_avg_call(SEXP, SEXP);
SEXP calc_fchange_call(SEXP, SEXP, SEXP);
SEXP detect_neofunctionalization_call(SEXP, SEXP, SEXP, SEXP);
SEXP detect_dosage_effect_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP detect_dosage_effect_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP detect_subfunctionalization_expert_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP detect_subfunctionalization_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP filter_paralogs_by_pattern_dosage_effect_call(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP filter_paralogs_by_pattern_subfunctionalization_call(SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP mask_check_state_call(SEXP, SEXP);
SEXP mask_chunk_count_call(SEXP);
SEXP calc_work_arr_paralog_subsets_size_call(SEXP, SEXP, SEXP);
SEXP omics_vector_RAP_projection_call(SEXP, SEXP, SEXP);
SEXP omics_field_RAP_projection_call(SEXP, SEXP, SEXP);
SEXP clock_hand_angle_between_vectors_call(SEXP, SEXP, SEXP);
SEXP clock_hand_angles_for_shift_vectors_call(SEXP, SEXP, SEXP);
SEXP relative_axes_changes_from_shift_vector_call(SEXP);
SEXP relative_axes_expression_from_expression_vector_call(SEXP);
SEXP compute_shift_vector_field_call(SEXP, SEXP, SEXP);
SEXP compute_tissue_versatility_call(SEXP, SEXP, SEXP);
SEXP perform_permutation_test_call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
SEXP compute_p_values_call(SEXP, SEXP, SEXP, SEXP);
SEXP compute_contributions_call(SEXP, SEXP, SEXP);
SEXP compute_all_contributions_call(SEXP, SEXP, SEXP, SEXP);
SEXP compute_baselines_factor_dependent_call(SEXP, SEXP, SEXP);
SEXP compute_velocity_trajectory_call(SEXP);
SEXP compute_acceleration_from_velocity_trajectory_call(SEXP, SEXP);
SEXP compute_velocity_trajectories_call(SEXP);
SEXP compute_acceleration_from_velocity_call(SEXP, SEXP);
SEXP compute_velocity_acceleration_contributions_expert_call(SEXP, SEXP);
SEXP compute_velocity_acceleration_contributions_call(SEXP, SEXP);
SEXP normalize_variable_timeseries_call(SEXP);
SEXP normalize_single_trajectory_call(SEXP);
SEXP normalize_all_trajectories_expert_call(SEXP);
SEXP normalize_all_trajectories_call(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"build_bst_index_call", (DL_FUNC) &build_bst_index_call, 1},
    {"bst_range_query_call", (DL_FUNC) &bst_range_query_call, 4},
    {"build_kd_index_call", (DL_FUNC) &build_kd_index_call, 2},
    {"build_spherical_kd_call", (DL_FUNC) &build_spherical_kd_call, 2},
    {"deserialize_char_helper_call", (DL_FUNC) &deserialize_char_helper_call, 3},
    {"deserialize_complex_helper_call", (DL_FUNC) &deserialize_complex_helper_call, 2},
    {"deserialize_int_helper_call", (DL_FUNC) &deserialize_int_helper_call, 2},
    {"deserialize_logical_helper_call", (DL_FUNC) &deserialize_logical_helper_call, 2},
    {"deserialize_real_helper_call", (DL_FUNC) &deserialize_real_helper_call, 2},
    {"serialize_char_helper_call", (DL_FUNC) &serialize_char_helper_call, 2},
    {"serialize_complex_helper_call", (DL_FUNC) &serialize_complex_helper_call, 2},
    {"serialize_int_helper_call", (DL_FUNC) &serialize_int_helper_call, 2},
    {"serialize_logical_helper_call", (DL_FUNC) &serialize_logical_helper_call, 2},
    {"serialize_real_helper_call", (DL_FUNC) &serialize_real_helper_call, 2},
    {"get_array_metadata_call", (DL_FUNC) &get_array_metadata_call, 2},
    {"loess_smooth_2d_call", (DL_FUNC) &loess_smooth_2d_call, 6},
    {"compute_edf_expert_call", (DL_FUNC) &compute_edf_expert_call, 2},
    {"compute_edf_call", (DL_FUNC) &compute_edf_call, 1},
    {"compute_scaled_distance_quantile_call", (DL_FUNC) &compute_scaled_distance_quantile_call, 4},
    {"cluster_factor_trajectories_k_means_call", (DL_FUNC) &cluster_factor_trajectories_k_means_call, 3},
    {"k_means_clustering_call", (DL_FUNC) &k_means_clustering_call, 3},
    {"linkage_clustering_call", (DL_FUNC) &linkage_clustering_call, 2},
    {"create_zip_archive_call", (DL_FUNC) &create_zip_archive_call, 3},
    {"save_tox_data_call", (DL_FUNC) &save_tox_data_call, 13},
    {"get_tox_data_dims_call", (DL_FUNC) &get_tox_data_dims_call, 1},
    {"read_tox_data_into_call", (DL_FUNC) &read_tox_data_into_call, 12},
    {"compute_gene_means_call", (DL_FUNC) &compute_gene_means_call, 1},
    {"compute_residuals_call", (DL_FUNC) &compute_residuals_call, 2},
    {"pool_means_call", (DL_FUNC) &pool_means_call, 3},
    {"pool_means_expert_call", (DL_FUNC) &pool_means_expert_call, 3},
    {"calc_neighborhood_size_call", (DL_FUNC) &calc_neighborhood_size_call, 4},
    {"construct_neighborhoods_call", (DL_FUNC) &construct_neighborhoods_call, 4},
    {"gjct_permutation_test_call", (DL_FUNC) &gjct_permutation_test_call, 9},
    {"gjct_permutation_test_expert_call", (DL_FUNC) &gjct_permutation_test_expert_call, 9},
    {"determine_shared_residual_range_expert_call", (DL_FUNC) &determine_shared_residual_range_expert_call, 3},
    {"determine_shared_residual_range_call", (DL_FUNC) &determine_shared_residual_range_call, 3},
    {"build_residual_histograms_call", (DL_FUNC) &build_residual_histograms_call, 4},
    {"compute_divergence_per_reference_point_call", (DL_FUNC) &compute_divergence_per_reference_point_call, 2},
    {"compute_weighted_global_divergence_call", (DL_FUNC) &compute_weighted_global_divergence_call, 3},
    {"fjct_compute_jsd_call", (DL_FUNC) &fjct_compute_jsd_call, 9},
    {"fjct_compute_jsd_expert_call", (DL_FUNC) &fjct_compute_jsd_expert_call, 6},
    {"fjct_compute_contribution_scores_call", (DL_FUNC) &fjct_compute_contribution_scores_call, 2},
    {"read_expression_vectors_tsv_call", (DL_FUNC) &read_expression_vectors_tsv_call, 8},
    {"read_gene_ids_from_tsv_file_call", (DL_FUNC) &read_gene_ids_from_tsv_file_call, 5},
    {"read_orthofinder_file_call", (DL_FUNC) &read_orthofinder_file_call, 5},
    {"get_unassigned_mask_call", (DL_FUNC) &get_unassigned_mask_call, 1},
    {"validate_data_structure_call", (DL_FUNC) &validate_data_structure_call, 9},
    {"validate_gene_to_family_mapping_call", (DL_FUNC) &validate_gene_to_family_mapping_call, 2},
    {"validate_expression_data_call", (DL_FUNC) &validate_expression_data_call, 2},
    {"validate_family_centroids_call", (DL_FUNC) &validate_family_centroids_call, 1},
    {"validate_shift_vectors_call", (DL_FUNC) &validate_shift_vectors_call, 5},
    {"validate_string_array_uniqueness_call", (DL_FUNC) &validate_string_array_uniqueness_call, 1},
    {"validate_all_data_call", (DL_FUNC) &validate_all_data_call, 11},
    {"euclidean_distance_call", (DL_FUNC) &euclidean_distance_call, 2},
    {"distance_to_centroid_call", (DL_FUNC) &distance_to_centroid_call, 3},
    {"mean_vector_call", (DL_FUNC) &mean_vector_call, 2},
    {"group_centroid_orthologs_expert_call", (DL_FUNC) &group_centroid_orthologs_expert_call, 4},
    {"group_centroid_orthologs_call", (DL_FUNC) &group_centroid_orthologs_call, 4},
    {"group_centroid_all_expert_call", (DL_FUNC) &group_centroid_all_expert_call, 3},
    {"group_centroid_all_call", (DL_FUNC) &group_centroid_all_call, 3},
    {"compute_family_scaling_expert_call", (DL_FUNC) &compute_family_scaling_expert_call, 9},
    {"compute_family_scaling_call", (DL_FUNC) &compute_family_scaling_call, 7},
    {"compute_rdi_expert_call", (DL_FUNC) &compute_rdi_expert_call, 3},
    {"compute_rdi_call", (DL_FUNC) &compute_rdi_call, 3},
    {"identify_outliers_call", (DL_FUNC) &identify_outliers_call, 4},
    {"detect_outliers_expert_call", (DL_FUNC) &detect_outliers_expert_call, 6},
    {"detect_outliers_call", (DL_FUNC) &detect_outliers_call, 4},
    {"loess_fit_plain_expert_call", (DL_FUNC) &loess_fit_plain_expert_call, 11},
    {"loess_fit_plain_call", (DL_FUNC) &loess_fit_plain_call, 9},
    {"loess_fit_robust_expert_call", (DL_FUNC) &loess_fit_robust_expert_call, 12},
    {"loess_fit_robust_call", (DL_FUNC) &loess_fit_robust_call, 10},
    {"tox_loess_required_workspace_call", (DL_FUNC) &tox_loess_required_workspace_call, 3},
    {"normalize_unit_length_call", (DL_FUNC) &normalize_unit_length_call, 1},
    {"normalization_pipeline_call", (DL_FUNC) &normalization_pipeline_call, 5},
    {"normalize_by_std_dev_call", (DL_FUNC) &normalize_by_std_dev_call, 3},
    {"root_mean_sq_normalization_call", (DL_FUNC) &root_mean_sq_normalization_call, 1},
    {"quantile_normalization_expert_call", (DL_FUNC) &quantile_normalization_expert_call, 1},
    {"quantile_normalization_call", (DL_FUNC) &quantile_normalization_call, 1},
    {"log2_transformation_call", (DL_FUNC) &log2_transformation_call, 1},
    {"calc_tiss_avg_call", (DL_FUNC) &calc_tiss_avg_call, 2},
    {"calc_fchange_call", (DL_FUNC) &calc_fchange_call, 3},
    {"detect_neofunctionalization_call", (DL_FUNC) &detect_neofunctionalization_call, 4},
    {"detect_dosage_effect_expert_call", (DL_FUNC) &detect_dosage_effect_expert_call, 7},
    {"detect_dosage_effect_call", (DL_FUNC) &detect_dosage_effect_call, 7},
    {"detect_subfunctionalization_expert_call", (DL_FUNC) &detect_subfunctionalization_expert_call, 8},
    {"detect_subfunctionalization_call", (DL_FUNC) &detect_subfunctionalization_call, 8},
    {"filter_paralogs_by_pattern_dosage_effect_call", (DL_FUNC) &filter_paralogs_by_pattern_dosage_effect_call, 5},
    {"filter_paralogs_by_pattern_subfunctionalization_call", (DL_FUNC) &filter_paralogs_by_pattern_subfunctionalization_call, 5},
    {"mask_check_state_call", (DL_FUNC) &mask_check_state_call, 2},
    {"mask_chunk_count_call", (DL_FUNC) &mask_chunk_count_call, 1},
    {"calc_work_arr_paralog_subsets_size_call", (DL_FUNC) &calc_work_arr_paralog_subsets_size_call, 3},
    {"omics_vector_RAP_projection_call", (DL_FUNC) &omics_vector_RAP_projection_call, 3},
    {"omics_field_RAP_projection_call", (DL_FUNC) &omics_field_RAP_projection_call, 3},
    {"clock_hand_angle_between_vectors_call", (DL_FUNC) &clock_hand_angle_between_vectors_call, 3},
    {"clock_hand_angles_for_shift_vectors_call", (DL_FUNC) &clock_hand_angles_for_shift_vectors_call, 3},
    {"relative_axes_changes_from_shift_vector_call", (DL_FUNC) &relative_axes_changes_from_shift_vector_call, 1},
    {"relative_axes_expression_from_expression_vector_call", (DL_FUNC) &relative_axes_expression_from_expression_vector_call, 1},
    {"compute_shift_vector_field_call", (DL_FUNC) &compute_shift_vector_field_call, 3},
    {"compute_tissue_versatility_call", (DL_FUNC) &compute_tissue_versatility_call, 3},
    {"perform_permutation_test_call", (DL_FUNC) &perform_permutation_test_call, 7},
    {"compute_p_values_call", (DL_FUNC) &compute_p_values_call, 4},
    {"compute_contributions_call", (DL_FUNC) &compute_contributions_call, 3},
    {"compute_all_contributions_call", (DL_FUNC) &compute_all_contributions_call, 4},
    {"compute_baselines_factor_dependent_call", (DL_FUNC) &compute_baselines_factor_dependent_call, 3},
    {"compute_velocity_trajectory_call", (DL_FUNC) &compute_velocity_trajectory_call, 1},
    {"compute_acceleration_from_velocity_trajectory_call", (DL_FUNC) &compute_acceleration_from_velocity_trajectory_call, 2},
    {"compute_velocity_trajectories_call", (DL_FUNC) &compute_velocity_trajectories_call, 1},
    {"compute_acceleration_from_velocity_call", (DL_FUNC) &compute_acceleration_from_velocity_call, 2},
    {"compute_velocity_acceleration_contributions_expert_call", (DL_FUNC) &compute_velocity_acceleration_contributions_expert_call, 2},
    {"compute_velocity_acceleration_contributions_call", (DL_FUNC) &compute_velocity_acceleration_contributions_call, 2},
    {"normalize_variable_timeseries_call", (DL_FUNC) &normalize_variable_timeseries_call, 1},
    {"normalize_single_trajectory_call", (DL_FUNC) &normalize_single_trajectory_call, 1},
    {"normalize_all_trajectories_expert_call", (DL_FUNC) &normalize_all_trajectories_expert_call, 1},
    {"normalize_all_trajectories_call", (DL_FUNC) &normalize_all_trajectories_call, 1},
    {NULL, NULL, 0}
};

void R_init_tensoromics(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
#endif  // R binding
