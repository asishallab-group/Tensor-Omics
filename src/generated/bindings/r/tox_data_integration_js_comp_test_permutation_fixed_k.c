// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void run_js_comp_test_on_sample_assignment_c(const int*, const int*, const int*, const double*, const int*, const int*, const int*, const int*, double*, int*, int*, const unsigned char*, const int*, const char*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const char*, const char*, const double*, const char*, const double*, const int*, int*);
void run_js_comp_test_permutation_c(const int*, const int*, const int*, const double*, const int*, const int*, const int*, double*, double*, double*, int*, const unsigned char*, const int*, const char*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const int*, const int*, const double*, const double*, const double*, const double*, const int*, const char*, const char*, const double*, const char*, const double*, const int*, int*);

SEXP run_js_comp_test_on_sample_assignment_call(SEXP expression, SEXP sample_assignment, SEXP n_points, SEXP n_neighbors, SEXP rerun_parameter_search, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP bin_count_mode, SEXP neighborhood_size_mode, SEXP delta_abs_tolerance, SEXP divergence_mode, SEXP grid_points_factor, SEXP grid_min_points) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[2];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[1];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[0];
    int n_samples = (int) Rf_length(sample_assignment);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_points = Rf_asInteger(n_points);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    unsigned char _v_rerun_parameter_search = (Rf_asLogical(rerun_parameter_search) == TRUE) ? 1 : 0;
    int _v_n_bootstraps = Rf_asInteger(n_bootstraps);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);
    double _v_min_neighbor_overlap = Rf_asReal(min_neighbor_overlap);
    double _v_succeeding_ci_overlap = Rf_asReal(succeeding_ci_overlap);
    double _v_delta_median_threshold = Rf_asReal(delta_median_threshold);
    double _v_delta_max_threshold = Rf_asReal(delta_max_threshold);
    double _v_delta_epsilon = Rf_asReal(delta_epsilon);
    int _v_delta_min_consecutive_transitions = Rf_asInteger(delta_min_consecutive_transitions);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    double _v_gamma_occupancy = Rf_asReal(gamma_occupancy);
    double _v_lower_residual_range_quantile = Rf_asReal(lower_residual_range_quantile);
    double _v_upper_residual_range_quantile = Rf_asReal(upper_residual_range_quantile);
    double _v_two_sided_bootstrapping_significance_level = Rf_asReal(two_sided_bootstrapping_significance_level);
    int _v_random_seed = Rf_asInteger(random_seed);
    double _v_delta_abs_tolerance = Rf_asReal(delta_abs_tolerance);
    double _v_grid_points_factor = Rf_asReal(grid_points_factor);
    int _v_grid_min_points = Rf_asInteger(grid_min_points);

    // convert what Fortran cannot take from R directly
    const int* _c_sample_assignment = tox_int_in(sample_assignment, "sample_assignment");
    char* _c_join_method = tox_char_in(join_method, 11, "join_method");
    char* _c_plateau_mode = tox_char_in(plateau_mode, 19, "plateau_mode");
    char* _c_bin_count_mode = tox_char_in(bin_count_mode, 29, "bin_count_mode");
    char* _c_neighborhood_size_mode = tox_char_in(neighborhood_size_mode, 23, "neighborhood_size_mode");
    char* _c_divergence_mode = tox_char_in(divergence_mode, 22, "divergence_mode");

    // outputs and work space
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    int used_n_points = 0;
    int used_n_neighbors = 0;
    int ierr = 0;

    run_js_comp_test_on_sample_assignment_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(expression),
        &n_samples,
        _c_sample_assignment,
        &_v_n_points,
        &_v_n_neighbors,
        REAL(global_js_divergence),
        &used_n_points,
        &used_n_neighbors,
        &_v_rerun_parameter_search,
        &_v_n_bootstraps,
        _c_join_method,
        &_v_min_residuals_per_bin,
        &_v_min_neighbor_overlap,
        &_v_succeeding_ci_overlap,
        _c_plateau_mode,
        &_v_delta_median_threshold,
        &_v_delta_max_threshold,
        &_v_delta_epsilon,
        &_v_delta_min_consecutive_transitions,
        &_v_m_min,
        &_v_m_max,
        &_v_gamma_occupancy,
        &_v_lower_residual_range_quantile,
        &_v_upper_residual_range_quantile,
        &_v_two_sided_bootstrapping_significance_level,
        &_v_random_seed,
        _c_bin_count_mode,
        _c_neighborhood_size_mode,
        &_v_delta_abs_tolerance,
        _c_divergence_mode,
        &_v_grid_points_factor,
        &_v_grid_min_points,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, global_js_divergence);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(used_n_points));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(used_n_neighbors));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("used_n_points"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("used_n_neighbors"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP run_js_comp_test_permutation_call(SEXP expression, SEXP n_points, SEXP n_neighbors, SEXP n_permutations, SEXP rerun_parameter_search, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP m_min, SEXP m_max, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed, SEXP bin_count_mode, SEXP neighborhood_size_mode, SEXP delta_abs_tolerance, SEXP divergence_mode, SEXP grid_points_factor, SEXP grid_min_points) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[2];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[1];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(expression, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_points = Rf_asInteger(n_points);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    int _v_n_permutations = Rf_asInteger(n_permutations);
    unsigned char _v_rerun_parameter_search = (Rf_asLogical(rerun_parameter_search) == TRUE) ? 1 : 0;
    int _v_n_bootstraps = Rf_asInteger(n_bootstraps);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);
    double _v_min_neighbor_overlap = Rf_asReal(min_neighbor_overlap);
    double _v_succeeding_ci_overlap = Rf_asReal(succeeding_ci_overlap);
    double _v_delta_median_threshold = Rf_asReal(delta_median_threshold);
    double _v_delta_max_threshold = Rf_asReal(delta_max_threshold);
    double _v_delta_epsilon = Rf_asReal(delta_epsilon);
    int _v_delta_min_consecutive_transitions = Rf_asInteger(delta_min_consecutive_transitions);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    double _v_gamma_occupancy = Rf_asReal(gamma_occupancy);
    double _v_lower_residual_range_quantile = Rf_asReal(lower_residual_range_quantile);
    double _v_upper_residual_range_quantile = Rf_asReal(upper_residual_range_quantile);
    double _v_two_sided_bootstrapping_significance_level = Rf_asReal(two_sided_bootstrapping_significance_level);
    int _v_random_seed = Rf_asInteger(random_seed);
    double _v_delta_abs_tolerance = Rf_asReal(delta_abs_tolerance);
    double _v_grid_points_factor = Rf_asReal(grid_points_factor);
    int _v_grid_min_points = Rf_asInteger(grid_min_points);

    // convert what Fortran cannot take from R directly
    char* _c_join_method = tox_char_in(join_method, 11, "join_method");
    char* _c_plateau_mode = tox_char_in(plateau_mode, 19, "plateau_mode");
    char* _c_bin_count_mode = tox_char_in(bin_count_mode, 29, "bin_count_mode");
    char* _c_neighborhood_size_mode = tox_char_in(neighborhood_size_mode, 23, "neighborhood_size_mode");
    char* _c_divergence_mode = tox_char_in(divergence_mode, 22, "divergence_mode");

    // outputs and work space
    SEXP global_js_divergence_observed = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    SEXP jsd_null = PROTECT(Rf_allocVector(REALSXP, _v_n_permutations * n_studies)); _nprot++;
    { SEXP _dim_jsd_null = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_jsd_null)[0] = _v_n_permutations; INTEGER(_dim_jsd_null)[1] = n_studies; Rf_setAttrib(jsd_null, R_DimSymbol, _dim_jsd_null); UNPROTECT(1); }
    int n_failed_permutations = 0;
    int ierr = 0;

    run_js_comp_test_permutation_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(expression),
        &_v_n_points,
        &_v_n_neighbors,
        &_v_n_permutations,
        REAL(global_js_divergence_observed),
        REAL(p_values),
        REAL(jsd_null),
        &n_failed_permutations,
        &_v_rerun_parameter_search,
        &_v_n_bootstraps,
        _c_join_method,
        &_v_min_residuals_per_bin,
        &_v_min_neighbor_overlap,
        &_v_succeeding_ci_overlap,
        _c_plateau_mode,
        &_v_delta_median_threshold,
        &_v_delta_max_threshold,
        &_v_delta_epsilon,
        &_v_delta_min_consecutive_transitions,
        &_v_m_min,
        &_v_m_max,
        &_v_gamma_occupancy,
        &_v_lower_residual_range_quantile,
        &_v_upper_residual_range_quantile,
        &_v_two_sided_bootstrapping_significance_level,
        &_v_random_seed,
        _c_bin_count_mode,
        _c_neighborhood_size_mode,
        &_v_delta_abs_tolerance,
        _c_divergence_mode,
        &_v_grid_points_factor,
        &_v_grid_min_points,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); _nprot++;
    SET_VECTOR_ELT(_out, 0, global_js_divergence_observed);
    SET_VECTOR_ELT(_out, 1, p_values);
    SET_VECTOR_ELT(_out, 2, jsd_null);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_failed_permutations));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence_observed"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("jsd_null"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_failed_permutations"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
