// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void generate_adaptive_js_comp_test_candidates_c(const int*, const int*, int*, int*, int*);
void run_js_comp_test_adaptive_parameter_search_c(const int*, const int*, const int*, const double*, const double*, const int*, const char*, const int*, int*, int*, int*, int*, int*, double*, double*, double*, unsigned char*, int*, int*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, double*, int*, unsigned char*, int*, int*, double*, int*, int*, int*, double*, double*, int*, int*, int*, int*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const double*, const double*, const int*, const int*, const char*, int*);

SEXP generate_adaptive_js_comp_test_candidates_call(SEXP max_n_genes_all_studies, SEXP n_studies) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int max_n_genes_all_studies_v = Rf_asInteger(max_n_genes_all_studies);
    int n_studies_v = Rf_asInteger(n_studies);

    // outputs and work space
    SEXP candidates_k_start_k_step_k_max = PROTECT(Rf_allocVector(INTSXP, 3 * 16)); nprot++;
    { SEXP candidates_k_start_k_step_k_max_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(candidates_k_start_k_step_k_max_dim)[0] = 3; INTEGER(candidates_k_start_k_step_k_max_dim)[1] = 16; Rf_setAttrib(candidates_k_start_k_step_k_max, R_DimSymbol, candidates_k_start_k_step_k_max_dim); UNPROTECT(1); }
    int n_candidates = 0;
    int ierr = 0;

    generate_adaptive_js_comp_test_candidates_c(
        &max_n_genes_all_studies_v,
        &n_studies_v,
        INTEGER(candidates_k_start_k_step_k_max),
        &n_candidates,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, candidates_k_start_k_step_k_max);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_candidates));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("candidates_k_start_k_step_k_max"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_candidates"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_adaptive_parameter_search_call(SEXP gene_means, SEXP residuals, SEXP n_bootstraps, SEXP join_method, SEXP max_n_points_candidate, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP tau, SEXP mad_distance_factor, SEXP max_pooled_residuals, SEXP min_study_neighbors, SEXP bin_count_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[1];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[0];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(residuals, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    int max_n_points_candidate_v = Rf_asInteger(max_n_points_candidate);
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

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);
    char* plateau_mode_c = tox_char_in(plateau_mode, 19);
    char* bin_count_mode_c = tox_char_in(bin_count_mode, 29);

    // outputs and work space
    int k_start = 0;
    int k_step = 0;
    int k_max = 0;
    int n_points = 0;
    SEXP n_bins_per_point = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v)); nprot++;
    SEXP shared_residual_range_low = PROTECT(Rf_allocVector(REALSXP, max_n_points_candidate_v)); nprot++;
    SEXP shared_residual_range_high = PROTECT(Rf_allocVector(REALSXP, max_n_points_candidate_v)); nprot++;
    SEXP best_candidate_confidence_interval = PROTECT(Rf_allocVector(REALSXP, 2 * n_studies)); nprot++;
    { SEXP best_candidate_confidence_interval_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(best_candidate_confidence_interval_dim)[0] = 2; INTEGER(best_candidate_confidence_interval_dim)[1] = n_studies; Rf_setAttrib(best_candidate_confidence_interval, R_DimSymbol, best_candidate_confidence_interval_dim); UNPROTECT(1); }
    unsigned char plateau_established = 0;
    int n_admissible_evaluated = 0;
    SEXP trace_k_start = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP trace_k_step = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP trace_k_max = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP trace_n_points = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP trace_global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_global_js_divergence_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_global_js_divergence_dim)[0] = n_studies; INTEGER(trace_global_js_divergence_dim)[1] = 16; Rf_setAttrib(trace_global_js_divergence, R_DimSymbol, trace_global_js_divergence_dim); UNPROTECT(1); }
    SEXP trace_ci_lower = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_ci_lower_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_ci_lower_dim)[0] = n_studies; INTEGER(trace_ci_lower_dim)[1] = 16; Rf_setAttrib(trace_ci_lower, R_DimSymbol, trace_ci_lower_dim); UNPROTECT(1); }
    SEXP trace_ci_upper = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_ci_upper_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_ci_upper_dim)[0] = n_studies; INTEGER(trace_ci_upper_dim)[1] = 16; Rf_setAttrib(trace_ci_upper, R_DimSymbol, trace_ci_upper_dim); UNPROTECT(1); }
    SEXP trace_ci_width = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_ci_width_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_ci_width_dim)[0] = n_studies; INTEGER(trace_ci_width_dim)[1] = 16; Rf_setAttrib(trace_ci_width, R_DimSymbol, trace_ci_width_dim); UNPROTECT(1); }
    SEXP trace_ci_width_relative = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_ci_width_relative_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_ci_width_relative_dim)[0] = n_studies; INTEGER(trace_ci_width_relative_dim)[1] = 16; Rf_setAttrib(trace_ci_width_relative, R_DimSymbol, trace_ci_width_relative_dim); UNPROTECT(1); }
    SEXP trace_ci_overlap = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_ci_overlap_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_ci_overlap_dim)[0] = n_studies; INTEGER(trace_ci_overlap_dim)[1] = 16; Rf_setAttrib(trace_ci_overlap, R_DimSymbol, trace_ci_overlap_dim); UNPROTECT(1); }
    SEXP trace_delta = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_delta_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_delta_dim)[0] = n_studies; INTEGER(trace_delta_dim)[1] = 16; Rf_setAttrib(trace_delta, R_DimSymbol, trace_delta_dim); UNPROTECT(1); }
    SEXP trace_delta_median = PROTECT(Rf_allocVector(REALSXP, 16)); nprot++;
    SEXP trace_delta_max = PROTECT(Rf_allocVector(REALSXP, 16)); nprot++;
    SEXP trace_selected_n_bins = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_selected_n_bins_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_selected_n_bins_dim)[0] = max_n_points_candidate_v; INTEGER(trace_selected_n_bins_dim)[1] = 16; Rf_setAttrib(trace_selected_n_bins, R_DimSymbol, trace_selected_n_bins_dim); UNPROTECT(1); }
    unsigned char* trace_occupancy_failed_c = tox_bool_alloc(max_n_points_candidate_v * 16);
    SEXP trace_n_pooled_residuals = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_n_pooled_residuals_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_n_pooled_residuals_dim)[0] = max_n_points_candidate_v; INTEGER(trace_n_pooled_residuals_dim)[1] = 16; Rf_setAttrib(trace_n_pooled_residuals, R_DimSymbol, trace_n_pooled_residuals_dim); UNPROTECT(1); }
    SEXP trace_min_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_min_bin_occupancy_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_min_bin_occupancy_dim)[0] = max_n_points_candidate_v; INTEGER(trace_min_bin_occupancy_dim)[1] = 16; Rf_setAttrib(trace_min_bin_occupancy, R_DimSymbol, trace_min_bin_occupancy_dim); UNPROTECT(1); }
    SEXP trace_mean_bin_occupancy = PROTECT(Rf_allocVector(REALSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_mean_bin_occupancy_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_mean_bin_occupancy_dim)[0] = max_n_points_candidate_v; INTEGER(trace_mean_bin_occupancy_dim)[1] = 16; Rf_setAttrib(trace_mean_bin_occupancy, R_DimSymbol, trace_mean_bin_occupancy_dim); UNPROTECT(1); }
    SEXP trace_max_bin_occupancy = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_max_bin_occupancy_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_max_bin_occupancy_dim)[0] = max_n_points_candidate_v; INTEGER(trace_max_bin_occupancy_dim)[1] = 16; Rf_setAttrib(trace_max_bin_occupancy, R_DimSymbol, trace_max_bin_occupancy_dim); UNPROTECT(1); }
    SEXP trace_sturges_bins = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_sturges_bins_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_sturges_bins_dim)[0] = max_n_points_candidate_v; INTEGER(trace_sturges_bins_dim)[1] = 16; Rf_setAttrib(trace_sturges_bins, R_DimSymbol, trace_sturges_bins_dim); UNPROTECT(1); }
    SEXP trace_fd_bins = PROTECT(Rf_allocVector(INTSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_fd_bins_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_fd_bins_dim)[0] = max_n_points_candidate_v; INTEGER(trace_fd_bins_dim)[1] = 16; Rf_setAttrib(trace_fd_bins, R_DimSymbol, trace_fd_bins_dim); UNPROTECT(1); }
    SEXP trace_shared_residual_range_low = PROTECT(Rf_allocVector(REALSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_shared_residual_range_low_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_shared_residual_range_low_dim)[0] = max_n_points_candidate_v; INTEGER(trace_shared_residual_range_low_dim)[1] = 16; Rf_setAttrib(trace_shared_residual_range_low, R_DimSymbol, trace_shared_residual_range_low_dim); UNPROTECT(1); }
    SEXP trace_shared_residual_range_high = PROTECT(Rf_allocVector(REALSXP, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_shared_residual_range_high_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_shared_residual_range_high_dim)[0] = max_n_points_candidate_v; INTEGER(trace_shared_residual_range_high_dim)[1] = 16; Rf_setAttrib(trace_shared_residual_range_high, R_DimSymbol, trace_shared_residual_range_high_dim); UNPROTECT(1); }
    int n_candidates_tried = 0;
    SEXP candidate_k_start = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP candidate_n_points = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP candidate_status = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    int ierr = 0;

    run_js_comp_test_adaptive_parameter_search_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(gene_means),
        REAL(residuals),
        &n_bootstraps_v,
        join_method_c,
        &max_n_points_candidate_v,
        &k_start,
        &k_step,
        &k_max,
        &n_points,
        INTEGER(n_bins_per_point),
        REAL(shared_residual_range_low),
        REAL(shared_residual_range_high),
        REAL(best_candidate_confidence_interval),
        &plateau_established,
        &n_admissible_evaluated,
        INTEGER(trace_k_start),
        INTEGER(trace_k_step),
        INTEGER(trace_k_max),
        INTEGER(trace_n_points),
        REAL(trace_global_js_divergence),
        REAL(trace_ci_lower),
        REAL(trace_ci_upper),
        REAL(trace_ci_width),
        REAL(trace_ci_width_relative),
        REAL(trace_ci_overlap),
        REAL(trace_delta),
        REAL(trace_delta_median),
        REAL(trace_delta_max),
        INTEGER(trace_selected_n_bins),
        trace_occupancy_failed_c,
        INTEGER(trace_n_pooled_residuals),
        INTEGER(trace_min_bin_occupancy),
        REAL(trace_mean_bin_occupancy),
        INTEGER(trace_max_bin_occupancy),
        INTEGER(trace_sturges_bins),
        INTEGER(trace_fd_bins),
        REAL(trace_shared_residual_range_low),
        REAL(trace_shared_residual_range_high),
        &n_candidates_tried,
        INTEGER(candidate_k_start),
        INTEGER(candidate_n_points),
        INTEGER(candidate_status),
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
        &ierr
    );

    // convert the outputs back
    SEXP trace_occupancy_failed = PROTECT(tox_bool_out(trace_occupancy_failed_c, max_n_points_candidate_v * 16)); nprot++;
    { SEXP trace_occupancy_failed_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_occupancy_failed_dim)[0] = max_n_points_candidate_v; INTEGER(trace_occupancy_failed_dim)[1] = 16; Rf_setAttrib(trace_occupancy_failed, R_DimSymbol, trace_occupancy_failed_dim); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 38)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(k_start));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(k_step));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(k_max));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_points));
    SET_VECTOR_ELT(_out, 4, n_bins_per_point);
    SET_VECTOR_ELT(_out, 5, shared_residual_range_low);
    SET_VECTOR_ELT(_out, 6, shared_residual_range_high);
    SET_VECTOR_ELT(_out, 7, best_candidate_confidence_interval);
    SET_VECTOR_ELT(_out, 8, Rf_ScalarLogical(plateau_established != 0));
    SET_VECTOR_ELT(_out, 9, Rf_ScalarInteger(n_admissible_evaluated));
    SET_VECTOR_ELT(_out, 10, trace_k_start);
    SET_VECTOR_ELT(_out, 11, trace_k_step);
    SET_VECTOR_ELT(_out, 12, trace_k_max);
    SET_VECTOR_ELT(_out, 13, trace_n_points);
    SET_VECTOR_ELT(_out, 14, trace_global_js_divergence);
    SET_VECTOR_ELT(_out, 15, trace_ci_lower);
    SET_VECTOR_ELT(_out, 16, trace_ci_upper);
    SET_VECTOR_ELT(_out, 17, trace_ci_width);
    SET_VECTOR_ELT(_out, 18, trace_ci_width_relative);
    SET_VECTOR_ELT(_out, 19, trace_ci_overlap);
    SET_VECTOR_ELT(_out, 20, trace_delta);
    SET_VECTOR_ELT(_out, 21, trace_delta_median);
    SET_VECTOR_ELT(_out, 22, trace_delta_max);
    SET_VECTOR_ELT(_out, 23, trace_selected_n_bins);
    SET_VECTOR_ELT(_out, 24, trace_occupancy_failed);
    SET_VECTOR_ELT(_out, 25, trace_n_pooled_residuals);
    SET_VECTOR_ELT(_out, 26, trace_min_bin_occupancy);
    SET_VECTOR_ELT(_out, 27, trace_mean_bin_occupancy);
    SET_VECTOR_ELT(_out, 28, trace_max_bin_occupancy);
    SET_VECTOR_ELT(_out, 29, trace_sturges_bins);
    SET_VECTOR_ELT(_out, 30, trace_fd_bins);
    SET_VECTOR_ELT(_out, 31, trace_shared_residual_range_low);
    SET_VECTOR_ELT(_out, 32, trace_shared_residual_range_high);
    SET_VECTOR_ELT(_out, 33, Rf_ScalarInteger(n_candidates_tried));
    SET_VECTOR_ELT(_out, 34, candidate_k_start);
    SET_VECTOR_ELT(_out, 35, candidate_n_points);
    SET_VECTOR_ELT(_out, 36, candidate_status);
    SET_VECTOR_ELT(_out, 37, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 38)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("k_start"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("k_step"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("k_max"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_points"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_bins_per_point"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("best_candidate_confidence_interval"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("plateau_established"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("n_admissible_evaluated"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("trace_k_start"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("trace_k_step"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("trace_k_max"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("trace_n_points"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("trace_global_js_divergence"));
    SET_STRING_ELT(_nms, 15, Rf_mkChar("trace_ci_lower"));
    SET_STRING_ELT(_nms, 16, Rf_mkChar("trace_ci_upper"));
    SET_STRING_ELT(_nms, 17, Rf_mkChar("trace_ci_width"));
    SET_STRING_ELT(_nms, 18, Rf_mkChar("trace_ci_width_relative"));
    SET_STRING_ELT(_nms, 19, Rf_mkChar("trace_ci_overlap"));
    SET_STRING_ELT(_nms, 20, Rf_mkChar("trace_delta"));
    SET_STRING_ELT(_nms, 21, Rf_mkChar("trace_delta_median"));
    SET_STRING_ELT(_nms, 22, Rf_mkChar("trace_delta_max"));
    SET_STRING_ELT(_nms, 23, Rf_mkChar("trace_selected_n_bins"));
    SET_STRING_ELT(_nms, 24, Rf_mkChar("trace_occupancy_failed"));
    SET_STRING_ELT(_nms, 25, Rf_mkChar("trace_n_pooled_residuals"));
    SET_STRING_ELT(_nms, 26, Rf_mkChar("trace_min_bin_occupancy"));
    SET_STRING_ELT(_nms, 27, Rf_mkChar("trace_mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 28, Rf_mkChar("trace_max_bin_occupancy"));
    SET_STRING_ELT(_nms, 29, Rf_mkChar("trace_sturges_bins"));
    SET_STRING_ELT(_nms, 30, Rf_mkChar("trace_fd_bins"));
    SET_STRING_ELT(_nms, 31, Rf_mkChar("trace_shared_residual_range_low"));
    SET_STRING_ELT(_nms, 32, Rf_mkChar("trace_shared_residual_range_high"));
    SET_STRING_ELT(_nms, 33, Rf_mkChar("n_candidates_tried"));
    SET_STRING_ELT(_nms, 34, Rf_mkChar("candidate_k_start"));
    SET_STRING_ELT(_nms, 35, Rf_mkChar("candidate_n_points"));
    SET_STRING_ELT(_nms, 36, Rf_mkChar("candidate_status"));
    SET_STRING_ELT(_nms, 37, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
