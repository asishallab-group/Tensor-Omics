// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void run_js_comp_test_on_sample_assignment_c(const int*, const int*, const int*, const double*, const int*, const int*, const int*, const int*, double*, int*, int*, const unsigned char*, const int*, const char*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const char*, const char*, const double*, const char*, const char*, const double*, const int*, const char*, int*);
void run_js_comp_test_adaptive_on_sample_assignment_c(const int*, const int*, const int*, const double*, const int*, const int*, const int*, const int*, const int*, double*, int*, int*, int*, int*, const unsigned char*, const int*, const char*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const double*, const double*, const int*, const int*, const char*, const double*, const char*, const char*, const char*, const char*, int*);
void run_js_comp_test_label_permutation_test_c(const int*, const int*, const int*, const double*, const int*, const int*, const int*, double*, double*, double*, int*, const unsigned char*, const int*, const char*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const char*, const char*, const double*, const char*, const char*, const double*, const int*, const char*, int*);
void run_js_comp_test_adaptive_label_permutation_test_c(const int*, const int*, const int*, const double*, const int*, const int*, const int*, const int*, double*, double*, double*, int*, const unsigned char*, const int*, const char*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const double*, const double*, const int*, const int*, const char*, const double*, const char*, const char*, const char*, const char*, int*);

SEXP run_js_comp_test_on_sample_assignment_call(SEXP expression, SEXP sample_assignment, SEXP n_points, SEXP n_neighbors, SEXP rerun_parameter_search, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP bin_count_mode, SEXP neighborhood_size_mode, SEXP delta_abs_tolerance, SEXP bootstrap_mode, SEXP grid_sample_size_mode, SEXP grid_points_factor, SEXP grid_min_points, SEXP grid_floor_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[2];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[1];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[0];
    int n_samples = (int) Rf_length(sample_assignment);

    // scalar inputs, pulled from their length-1 vectors
    int n_points_v = Rf_asInteger(n_points);
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    unsigned char rerun_parameter_search_v = (Rf_asLogical(rerun_parameter_search) == TRUE) ? 1 : 0;
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double min_neighbor_overlap_v = Rf_asReal(min_neighbor_overlap);
    double succeeding_ci_overlap_v = Rf_asReal(succeeding_ci_overlap);
    double delta_median_threshold_v = Rf_asReal(delta_median_threshold);
    double delta_max_threshold_v = Rf_asReal(delta_max_threshold);
    double delta_epsilon_v = Rf_asReal(delta_epsilon);
    int delta_min_consecutive_transitions_v = Rf_asInteger(delta_min_consecutive_transitions);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);
    int random_seed_v = Rf_asInteger(random_seed);
    double delta_abs_tolerance_v = Rf_asReal(delta_abs_tolerance);
    double grid_points_factor_v = Rf_asReal(grid_points_factor);
    int grid_min_points_v = Rf_asInteger(grid_min_points);

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);
    char* plateau_mode_c = tox_char_in(plateau_mode, 19);
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);
    char* neighborhood_size_mode_c = tox_char_in(neighborhood_size_mode, 23);
    char* bootstrap_mode_c = tox_char_in(bootstrap_mode, 19);
    char* grid_sample_size_mode_c = tox_char_in(grid_sample_size_mode, 26);
    char* grid_floor_mode_c = tox_char_in(grid_floor_mode, 16);

    // outputs and work space
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    int used_n_points = 0;
    int used_n_neighbors = 0;
    int ierr = 0;

    run_js_comp_test_on_sample_assignment_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(expression),
        &n_samples,
        INTEGER(sample_assignment),
        &n_points_v,
        &n_neighbors_v,
        REAL(global_js_divergence),
        &used_n_points,
        &used_n_neighbors,
        &rerun_parameter_search_v,
        &n_bootstraps_v,
        join_method_c,
        &min_residuals_per_bin_v,
        &min_neighbor_overlap_v,
        &succeeding_ci_overlap_v,
        plateau_mode_c,
        &delta_median_threshold_v,
        &delta_max_threshold_v,
        &delta_epsilon_v,
        &delta_min_consecutive_transitions_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        &two_sided_bootstrapping_significance_level_v,
        &random_seed_v,
        bin_count_mode_c,
        neighborhood_size_mode_c,
        &delta_abs_tolerance_v,
        bootstrap_mode_c,
        grid_sample_size_mode_c,
        &grid_points_factor_v,
        &grid_min_points_v,
        grid_floor_mode_c,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, global_js_divergence);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(used_n_points));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(used_n_neighbors));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("used_n_points"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("used_n_neighbors"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_adaptive_on_sample_assignment_call(SEXP expression, SEXP sample_assignment, SEXP k_start, SEXP k_step, SEXP k_max, SEXP rerun_parameter_search, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP tau, SEXP mad_distance_factor, SEXP max_pooled_residuals, SEXP min_study_neighbors, SEXP bin_count_mode, SEXP delta_abs_tolerance, SEXP bootstrap_mode, SEXP tau_reference_mode, SEXP overlap_extension_mode, SEXP seed_step_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[2];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[1];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[0];
    int n_samples = (int) Rf_length(sample_assignment);

    // scalar inputs, pulled from their length-1 vectors
    int k_start_v = Rf_asInteger(k_start);
    int k_step_v = Rf_asInteger(k_step);
    int k_max_v = Rf_asInteger(k_max);
    unsigned char rerun_parameter_search_v = (Rf_asLogical(rerun_parameter_search) == TRUE) ? 1 : 0;
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double min_neighbor_overlap_v = Rf_asReal(min_neighbor_overlap);
    double succeeding_ci_overlap_v = Rf_asReal(succeeding_ci_overlap);
    double delta_median_threshold_v = Rf_asReal(delta_median_threshold);
    double delta_max_threshold_v = Rf_asReal(delta_max_threshold);
    double delta_epsilon_v = Rf_asReal(delta_epsilon);
    int delta_min_consecutive_transitions_v = Rf_asInteger(delta_min_consecutive_transitions);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);
    int random_seed_v = Rf_asInteger(random_seed);
    double tau_v = Rf_asReal(tau);
    double mad_distance_factor_v = Rf_asReal(mad_distance_factor);
    int max_pooled_residuals_v = Rf_asInteger(max_pooled_residuals);
    int min_study_neighbors_v = Rf_asInteger(min_study_neighbors);
    double delta_abs_tolerance_v = Rf_asReal(delta_abs_tolerance);

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);
    char* plateau_mode_c = tox_char_in(plateau_mode, 19);
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);
    char* bootstrap_mode_c = tox_char_in(bootstrap_mode, 19);
    char* tau_reference_mode_c = tox_char_in(tau_reference_mode, 28);
    char* overlap_extension_mode_c = tox_char_in(overlap_extension_mode, 26);
    char* seed_step_mode_c = tox_char_in(seed_step_mode, 20);

    // outputs and work space
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    int used_k_start = 0;
    int used_k_step = 0;
    int used_k_max = 0;
    int n_points = 0;
    int ierr = 0;

    run_js_comp_test_adaptive_on_sample_assignment_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(expression),
        &n_samples,
        INTEGER(sample_assignment),
        &k_start_v,
        &k_step_v,
        &k_max_v,
        REAL(global_js_divergence),
        &used_k_start,
        &used_k_step,
        &used_k_max,
        &n_points,
        &rerun_parameter_search_v,
        &n_bootstraps_v,
        join_method_c,
        &min_residuals_per_bin_v,
        &min_neighbor_overlap_v,
        &succeeding_ci_overlap_v,
        plateau_mode_c,
        &delta_median_threshold_v,
        &delta_max_threshold_v,
        &delta_epsilon_v,
        &delta_min_consecutive_transitions_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        &two_sided_bootstrapping_significance_level_v,
        &random_seed_v,
        &tau_v,
        &mad_distance_factor_v,
        &max_pooled_residuals_v,
        &min_study_neighbors_v,
        bin_count_mode_c,
        &delta_abs_tolerance_v,
        bootstrap_mode_c,
        tau_reference_mode_c,
        overlap_extension_mode_c,
        seed_step_mode_c,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
    SET_VECTOR_ELT(_out, 0, global_js_divergence);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(used_k_start));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(used_k_step));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(used_k_max));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(n_points));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("used_k_start"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("used_k_step"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("used_k_max"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_points"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_label_permutation_test_call(SEXP expression, SEXP n_points, SEXP n_neighbors, SEXP n_permutations, SEXP rerun_parameter_search, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP bin_count_mode, SEXP neighborhood_size_mode, SEXP delta_abs_tolerance, SEXP bootstrap_mode, SEXP grid_sample_size_mode, SEXP grid_points_factor, SEXP grid_min_points, SEXP grid_floor_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[2];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[1];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int n_points_v = Rf_asInteger(n_points);
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    int n_permutations_v = Rf_asInteger(n_permutations);
    unsigned char rerun_parameter_search_v = (Rf_asLogical(rerun_parameter_search) == TRUE) ? 1 : 0;
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double min_neighbor_overlap_v = Rf_asReal(min_neighbor_overlap);
    double succeeding_ci_overlap_v = Rf_asReal(succeeding_ci_overlap);
    double delta_median_threshold_v = Rf_asReal(delta_median_threshold);
    double delta_max_threshold_v = Rf_asReal(delta_max_threshold);
    double delta_epsilon_v = Rf_asReal(delta_epsilon);
    int delta_min_consecutive_transitions_v = Rf_asInteger(delta_min_consecutive_transitions);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);
    int random_seed_v = Rf_asInteger(random_seed);
    double delta_abs_tolerance_v = Rf_asReal(delta_abs_tolerance);
    double grid_points_factor_v = Rf_asReal(grid_points_factor);
    int grid_min_points_v = Rf_asInteger(grid_min_points);

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);
    char* plateau_mode_c = tox_char_in(plateau_mode, 19);
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);
    char* neighborhood_size_mode_c = tox_char_in(neighborhood_size_mode, 23);
    char* bootstrap_mode_c = tox_char_in(bootstrap_mode, 19);
    char* grid_sample_size_mode_c = tox_char_in(grid_sample_size_mode, 26);
    char* grid_floor_mode_c = tox_char_in(grid_floor_mode, 16);

    // outputs and work space
    SEXP global_js_divergence_observed = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP jsd_null = PROTECT(Rf_allocVector(REALSXP, n_permutations_v * n_studies)); nprot++;
    { SEXP jsd_null_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(jsd_null_dim)[0] = n_permutations_v; INTEGER(jsd_null_dim)[1] = n_studies; Rf_setAttrib(jsd_null, R_DimSymbol, jsd_null_dim); UNPROTECT(1); }
    int n_failed_permutations = 0;
    int ierr = 0;

    run_js_comp_test_label_permutation_test_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(expression),
        &n_points_v,
        &n_neighbors_v,
        &n_permutations_v,
        REAL(global_js_divergence_observed),
        REAL(p_values),
        REAL(jsd_null),
        &n_failed_permutations,
        &rerun_parameter_search_v,
        &n_bootstraps_v,
        join_method_c,
        &min_residuals_per_bin_v,
        &min_neighbor_overlap_v,
        &succeeding_ci_overlap_v,
        plateau_mode_c,
        &delta_median_threshold_v,
        &delta_max_threshold_v,
        &delta_epsilon_v,
        &delta_min_consecutive_transitions_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        &two_sided_bootstrapping_significance_level_v,
        &random_seed_v,
        bin_count_mode_c,
        neighborhood_size_mode_c,
        &delta_abs_tolerance_v,
        bootstrap_mode_c,
        grid_sample_size_mode_c,
        &grid_points_factor_v,
        &grid_min_points_v,
        grid_floor_mode_c,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, global_js_divergence_observed);
    SET_VECTOR_ELT(_out, 1, p_values);
    SET_VECTOR_ELT(_out, 2, jsd_null);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_failed_permutations));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence_observed"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("jsd_null"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_failed_permutations"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_adaptive_label_permutation_test_call(SEXP expression, SEXP k_start, SEXP k_step, SEXP k_max, SEXP n_permutations, SEXP rerun_parameter_search, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP tau, SEXP mad_distance_factor, SEXP max_pooled_residuals, SEXP min_study_neighbors, SEXP bin_count_mode, SEXP delta_abs_tolerance, SEXP bootstrap_mode, SEXP tau_reference_mode, SEXP overlap_extension_mode, SEXP seed_step_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[2];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[1];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int k_start_v = Rf_asInteger(k_start);
    int k_step_v = Rf_asInteger(k_step);
    int k_max_v = Rf_asInteger(k_max);
    int n_permutations_v = Rf_asInteger(n_permutations);
    unsigned char rerun_parameter_search_v = (Rf_asLogical(rerun_parameter_search) == TRUE) ? 1 : 0;
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double min_neighbor_overlap_v = Rf_asReal(min_neighbor_overlap);
    double succeeding_ci_overlap_v = Rf_asReal(succeeding_ci_overlap);
    double delta_median_threshold_v = Rf_asReal(delta_median_threshold);
    double delta_max_threshold_v = Rf_asReal(delta_max_threshold);
    double delta_epsilon_v = Rf_asReal(delta_epsilon);
    int delta_min_consecutive_transitions_v = Rf_asInteger(delta_min_consecutive_transitions);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);
    double lower_residual_range_quantile_v = Rf_asReal(lower_residual_range_quantile);
    double upper_residual_range_quantile_v = Rf_asReal(upper_residual_range_quantile);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);
    int random_seed_v = Rf_asInteger(random_seed);
    double tau_v = Rf_asReal(tau);
    double mad_distance_factor_v = Rf_asReal(mad_distance_factor);
    int max_pooled_residuals_v = Rf_asInteger(max_pooled_residuals);
    int min_study_neighbors_v = Rf_asInteger(min_study_neighbors);
    double delta_abs_tolerance_v = Rf_asReal(delta_abs_tolerance);

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);
    char* plateau_mode_c = tox_char_in(plateau_mode, 19);
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);
    char* bootstrap_mode_c = tox_char_in(bootstrap_mode, 19);
    char* tau_reference_mode_c = tox_char_in(tau_reference_mode, 28);
    char* overlap_extension_mode_c = tox_char_in(overlap_extension_mode, 26);
    char* seed_step_mode_c = tox_char_in(seed_step_mode, 20);

    // outputs and work space
    SEXP global_js_divergence_observed = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP jsd_null = PROTECT(Rf_allocVector(REALSXP, n_permutations_v * n_studies)); nprot++;
    { SEXP jsd_null_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(jsd_null_dim)[0] = n_permutations_v; INTEGER(jsd_null_dim)[1] = n_studies; Rf_setAttrib(jsd_null, R_DimSymbol, jsd_null_dim); UNPROTECT(1); }
    int n_failed_permutations = 0;
    int ierr = 0;

    run_js_comp_test_adaptive_label_permutation_test_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(expression),
        &k_start_v,
        &k_step_v,
        &k_max_v,
        &n_permutations_v,
        REAL(global_js_divergence_observed),
        REAL(p_values),
        REAL(jsd_null),
        &n_failed_permutations,
        &rerun_parameter_search_v,
        &n_bootstraps_v,
        join_method_c,
        &min_residuals_per_bin_v,
        &min_neighbor_overlap_v,
        &succeeding_ci_overlap_v,
        plateau_mode_c,
        &delta_median_threshold_v,
        &delta_max_threshold_v,
        &delta_epsilon_v,
        &delta_min_consecutive_transitions_v,
        &m_min_v,
        &m_max_v,
        &gamma_occupancy_v,
        &lower_residual_range_quantile_v,
        &upper_residual_range_quantile_v,
        &two_sided_bootstrapping_significance_level_v,
        &random_seed_v,
        &tau_v,
        &mad_distance_factor_v,
        &max_pooled_residuals_v,
        &min_study_neighbors_v,
        bin_count_mode_c,
        &delta_abs_tolerance_v,
        bootstrap_mode_c,
        tau_reference_mode_c,
        overlap_extension_mode_c,
        seed_step_mode_c,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, global_js_divergence_observed);
    SET_VECTOR_ELT(_out, 1, p_values);
    SET_VECTOR_ELT(_out, 2, jsd_null);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_failed_permutations));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence_observed"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("jsd_null"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_failed_permutations"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
