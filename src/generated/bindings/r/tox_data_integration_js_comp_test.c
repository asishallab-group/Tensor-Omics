// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void estimate_bin_count_c(const double*, const int*, const int*, const int*, const double*, int*, int*, int*, int*);
void estimate_bin_count_expert_c(const double*, const int*, const int*, const int*, const int*, const double*, int*, int*, int*, int*);
void determine_bin_count_occupancy_c(const double*, const int*, const int*, const int*, const double*, int*, unsigned char*, int*, int*, double*, int*, int*, int*, const int*, const int*, const int*, const double*, int*);
void determine_bin_count_occupancy_expert_c(const double*, const int*, const int*, const int*, const int*, const double*, int*, unsigned char*, int*, int*, double*, int*, int*, int*, int*, const int*, const int*, const int*, const double*, int*);
void generate_js_comp_test_candidates_c(const int*, const double*, const int*, const int*, const double*, int*, int*, int*, int*);
void generate_js_comp_test_candidates_expert_c(const int*, const double*, const int*, const int*, const int*, const double*, int*, int*, int*, int*);
void check_neighborhood_overlaps_c(const int*, const int*, const double*, unsigned char*, int*);
void check_mean_pmf_min_counts_c(const int*, const int*, const int*, const int*, const int*, unsigned char*, int*);
void check_plateau_condition_c(const double*, double*, const int*, int*, int*, const int*, const char*, const double*, unsigned char*, int*);
void check_effect_size_plateau_condition_c(const double*, const double*, const int*, const unsigned char*, const double*, const double*, const double*, const int*, int*, double*, double*, double*, unsigned char*, int*);
void create_mean_pmf_c(const double*, const int*, const int*, const int*, const int*, const int*, double*, int*, int*, int*);
void create_mean_pmf_only_c(const double*, const int*, const int*, const int*, double*, int*);
void bootstrap_histogram_c(const int*, const int*, const int*, const int*, const int*, const int*, const int*, double*, const double*, const int*, int*);
void run_js_comp_test_c(const int*, const int*, const int*, const int*, const int*, const int*, const double*, const double*, const int*, const double*, const double*, int*, int*, double*, int*, int*, double*, int*, int*, double*, double*, double*, double*, const int*, const int*, int*);
void run_js_comp_test_parameter_search_c(const int*, const int*, const int*, const double*, const double*, const double*, const int*, const char*, int*, int*, int*, double*, unsigned char*, int*, int*, int*, double*, double*, double*, double*, double*, double*, double*, double*, const int*, const double*, const double*, const char*, const double*, const double*, const double*, const int*, const double*, const int*, int*);

SEXP estimate_bin_count_call(SEXP residuals, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(residuals);

    // scalar inputs, pulled from their length-1 vectors
    int max_n_reps_all_studies_v = Rf_asInteger(max_n_reps_all_studies);
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // outputs and work space
    int n_bins = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    estimate_bin_count_c(
        REAL(residuals),
        &n_residuals,
        &max_n_reps_all_studies_v,
        &n_neighbors_v,
        &shared_residual_range_v,
        &n_bins,
        &sturges_bins,
        &fd_bins,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP estimate_bin_count_expert_call(SEXP residuals, SEXP residuals_perm, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(residuals);

    // scalar inputs, pulled from their length-1 vectors
    int max_n_reps_all_studies_v = Rf_asInteger(max_n_reps_all_studies);
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // outputs and work space
    int n_bins = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    estimate_bin_count_expert_c(
        REAL(residuals),
        INTEGER(residuals_perm),
        &n_residuals,
        &max_n_reps_all_studies_v,
        &n_neighbors_v,
        &shared_residual_range_v,
        &n_bins,
        &sturges_bins,
        &fd_bins,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP determine_bin_count_occupancy_call(SEXP pooled_residuals, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP shared_residual_range, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin, SEXP gamma_occupancy) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int max_n_reps_all_studies_v = Rf_asInteger(max_n_reps_all_studies);
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char occupancy_failed = 0;
    int n_pooled_residuals = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    determine_bin_count_occupancy_c(
        REAL(pooled_residuals),
        &n_residuals,
        &max_n_reps_all_studies_v,
        &n_neighbors_v,
        &shared_residual_range_v,
        &selected_n_bins,
        &occupancy_failed,
        &n_pooled_residuals,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &sturges_bins,
        &fd_bins,
        &m_min_v,
        &m_max_v,
        &min_residuals_per_bin_v,
        &gamma_occupancy_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(occupancy_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_pooled_residuals));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP determine_bin_count_occupancy_expert_call(SEXP pooled_residuals, SEXP pooled_residuals_perm, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP shared_residual_range, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin, SEXP gamma_occupancy) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int max_n_reps_all_studies_v = Rf_asInteger(max_n_reps_all_studies);
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int m_min_v = Rf_asInteger(m_min);
    int m_max_v = Rf_asInteger(m_max);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double gamma_occupancy_v = Rf_asReal(gamma_occupancy);

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char occupancy_failed = 0;
    int n_pooled_residuals = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int* tmp_bin_counts = (int*) R_alloc(256, sizeof(int));
    int ierr = 0;

    determine_bin_count_occupancy_expert_c(
        REAL(pooled_residuals),
        INTEGER(pooled_residuals_perm),
        &n_residuals,
        &max_n_reps_all_studies_v,
        &n_neighbors_v,
        &shared_residual_range_v,
        &selected_n_bins,
        &occupancy_failed,
        &n_pooled_residuals,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &sturges_bins,
        &fd_bins,
        tmp_bin_counts,
        &m_min_v,
        &m_max_v,
        &min_residuals_per_bin_v,
        &gamma_occupancy_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(occupancy_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_pooled_residuals));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP generate_js_comp_test_candidates_call(SEXP max_n_genes_all_studies, SEXP residuals, SEXP max_n_reps_all_studies, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(residuals);

    // scalar inputs, pulled from their length-1 vectors
    int max_n_genes_all_studies_v = Rf_asInteger(max_n_genes_all_studies);
    int max_n_reps_all_studies_v = Rf_asInteger(max_n_reps_all_studies);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // outputs and work space
    SEXP candidates_n_points_n_neighbors = PROTECT(Rf_allocVector(INTSXP, 2 * 16)); nprot++;
    { SEXP candidates_n_points_n_neighbors_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(candidates_n_points_n_neighbors_dim)[0] = 2; INTEGER(candidates_n_points_n_neighbors_dim)[1] = 16; Rf_setAttrib(candidates_n_points_n_neighbors, R_DimSymbol, candidates_n_points_n_neighbors_dim); UNPROTECT(1); }
    SEXP n_bins_candidates = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    int n_candidates = 0;
    int ierr = 0;

    generate_js_comp_test_candidates_c(
        &max_n_genes_all_studies_v,
        REAL(residuals),
        &n_residuals,
        &max_n_reps_all_studies_v,
        &shared_residual_range_v,
        INTEGER(candidates_n_points_n_neighbors),
        INTEGER(n_bins_candidates),
        &n_candidates,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, candidates_n_points_n_neighbors);
    SET_VECTOR_ELT(_out, 1, n_bins_candidates);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_candidates));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("candidates_n_points_n_neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_bins_candidates"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_candidates"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP generate_js_comp_test_candidates_expert_call(SEXP max_n_genes_all_studies, SEXP residuals, SEXP residuals_perm, SEXP max_n_reps_all_studies, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(residuals);

    // scalar inputs, pulled from their length-1 vectors
    int max_n_genes_all_studies_v = Rf_asInteger(max_n_genes_all_studies);
    int max_n_reps_all_studies_v = Rf_asInteger(max_n_reps_all_studies);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // outputs and work space
    SEXP candidates_n_points_n_neighbors = PROTECT(Rf_allocVector(INTSXP, 2 * 16)); nprot++;
    { SEXP candidates_n_points_n_neighbors_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(candidates_n_points_n_neighbors_dim)[0] = 2; INTEGER(candidates_n_points_n_neighbors_dim)[1] = 16; Rf_setAttrib(candidates_n_points_n_neighbors, R_DimSymbol, candidates_n_points_n_neighbors_dim); UNPROTECT(1); }
    SEXP n_bins_candidates = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    int n_candidates = 0;
    int ierr = 0;

    generate_js_comp_test_candidates_expert_c(
        &max_n_genes_all_studies_v,
        REAL(residuals),
        INTEGER(residuals_perm),
        &n_residuals,
        &max_n_reps_all_studies_v,
        &shared_residual_range_v,
        INTEGER(candidates_n_points_n_neighbors),
        INTEGER(n_bins_candidates),
        &n_candidates,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, candidates_n_points_n_neighbors);
    SET_VECTOR_ELT(_out, 1, n_bins_candidates);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_candidates));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("candidates_n_points_n_neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_bins_candidates"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_candidates"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP check_neighborhood_overlaps_call(SEXP neighborhood_range, SEXP min_neighbor_overlap) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(neighborhood_range, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double min_neighbor_overlap_v = Rf_asReal(min_neighbor_overlap);

    // outputs and work space
    unsigned char all_have_min_neighbor_overlap = 0;
    int ierr = 0;

    check_neighborhood_overlaps_c(
        INTEGER(neighborhood_range),
        &n_points,
        &min_neighbor_overlap_v,
        &all_have_min_neighbor_overlap,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(all_have_min_neighbor_overlap != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("all_have_min_neighbor_overlap"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP check_mean_pmf_min_counts_call(SEXP mean_pmf_counts, SEXP n_bins_per_point, SEXP min_count) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int min_count_v = Rf_asInteger(min_count);

    // outputs and work space
    unsigned char all_bins_have_min_count = 0;
    int ierr = 0;

    check_mean_pmf_min_counts_c(
        INTEGER(mean_pmf_counts),
        &n_bins,
        INTEGER(n_bins_per_point),
        &n_points,
        &min_count_v,
        &all_bins_have_min_count,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(all_bins_have_min_count != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("all_bins_have_min_count"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP check_plateau_condition_call(SEXP confidence_interval, SEXP best_candidate_pair_confidence_interval, SEXP best_candidate_index, SEXP best_exceeded_ci_overlap_count, SEXP candidate_index, SEXP join_method, SEXP succeeding_ci_overlap) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(confidence_interval, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int best_candidate_index_v = Rf_asInteger(best_candidate_index);
    int best_exceeded_ci_overlap_count_v = Rf_asInteger(best_exceeded_ci_overlap_count);
    int candidate_index_v = Rf_asInteger(candidate_index);
    double succeeding_ci_overlap_v = Rf_asReal(succeeding_ci_overlap);

    // copy what is modified in place, so the caller's stays intact
    SEXP best_candidate_pair_confidence_interval_out = PROTECT(Rf_duplicate(best_candidate_pair_confidence_interval)); nprot++;

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);

    // outputs and work space
    unsigned char plateau_found = 0;
    int ierr = 0;

    check_plateau_condition_c(
        REAL(confidence_interval),
        REAL(best_candidate_pair_confidence_interval_out),
        &n_studies,
        &best_candidate_index_v,
        &best_exceeded_ci_overlap_count_v,
        &candidate_index_v,
        join_method_c,
        &succeeding_ci_overlap_v,
        &plateau_found,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, best_candidate_pair_confidence_interval_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(best_candidate_index_v));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(best_exceeded_ci_overlap_count_v));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarLogical(plateau_found != 0));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("best_candidate_pair_confidence_interval"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("best_candidate_index"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("best_exceeded_ci_overlap_count"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("plateau_found"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP check_effect_size_plateau_condition_call(SEXP global_js_divergence, SEXP prev_global_js_divergence, SEXP has_previous, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP n_consecutive_ok) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = (int) Rf_length(global_js_divergence);

    // scalar inputs, pulled from their length-1 vectors
    unsigned char has_previous_v = (Rf_asLogical(has_previous) == TRUE) ? 1 : 0;
    double delta_median_threshold_v = Rf_asReal(delta_median_threshold);
    double delta_max_threshold_v = Rf_asReal(delta_max_threshold);
    double delta_epsilon_v = Rf_asReal(delta_epsilon);
    int delta_min_consecutive_transitions_v = Rf_asInteger(delta_min_consecutive_transitions);
    int n_consecutive_ok_v = Rf_asInteger(n_consecutive_ok);

    // outputs and work space
    SEXP delta = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    double delta_median = 0;
    double delta_max = 0;
    unsigned char plateau_found = 0;
    int ierr = 0;

    check_effect_size_plateau_condition_c(
        REAL(global_js_divergence),
        REAL(prev_global_js_divergence),
        &n_studies,
        &has_previous_v,
        &delta_median_threshold_v,
        &delta_max_threshold_v,
        &delta_epsilon_v,
        &delta_min_consecutive_transitions_v,
        &n_consecutive_ok_v,
        REAL(delta),
        &delta_median,
        &delta_max,
        &plateau_found,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_consecutive_ok_v));
    SET_VECTOR_ELT(_out, 1, delta);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(delta_median));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(delta_max));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarLogical(plateau_found != 0));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_consecutive_ok"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("delta"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("delta_median"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("delta_max"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("plateau_found"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP create_mean_pmf_call(SEXP pmfs, SEXP counts, SEXP included_n_reps) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[1];
    int n_studies = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[2];

    // outputs and work space
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_bins * n_points)); nprot++;
    { SEXP mean_pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_dim)[0] = n_bins; INTEGER(mean_pmf_dim)[1] = n_points; Rf_setAttrib(mean_pmf, R_DimSymbol, mean_pmf_dim); UNPROTECT(1); }
    SEXP mean_pmf_included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, n_bins * n_points)); nprot++;
    { SEXP mean_pmf_counts_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_counts_dim)[0] = n_bins; INTEGER(mean_pmf_counts_dim)[1] = n_points; Rf_setAttrib(mean_pmf_counts, R_DimSymbol, mean_pmf_counts_dim); UNPROTECT(1); }
    int ierr = 0;

    create_mean_pmf_c(
        REAL(pmfs),
        INTEGER(counts),
        &n_bins,
        &n_points,
        &n_studies,
        INTEGER(included_n_reps),
        REAL(mean_pmf),
        INTEGER(mean_pmf_included_n_reps),
        INTEGER(mean_pmf_counts),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, mean_pmf);
    SET_VECTOR_ELT(_out, 1, mean_pmf_included_n_reps);
    SET_VECTOR_ELT(_out, 2, mean_pmf_counts);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("mean_pmf_included_n_reps"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("mean_pmf_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP create_mean_pmf_only_call(SEXP pmfs) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[1];
    int n_studies = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[2];

    // outputs and work space
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_bins * n_points)); nprot++;
    { SEXP mean_pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_dim)[0] = n_bins; INTEGER(mean_pmf_dim)[1] = n_points; Rf_setAttrib(mean_pmf, R_DimSymbol, mean_pmf_dim); UNPROTECT(1); }
    int ierr = 0;

    create_mean_pmf_only_c(
        REAL(pmfs),
        &n_bins,
        &n_points,
        &n_studies,
        REAL(mean_pmf),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, mean_pmf);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bootstrap_histogram_call(SEXP n_bootstraps, SEXP mean_pmf_counts, SEXP mean_pmf_included_n_reps, SEXP included_n_reps, SEXP confidence_interval, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];
    int n_studies = INTEGER(Rf_getAttrib(included_n_reps, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);
    int random_seed_v = Rf_asInteger(random_seed);

    // copy what is modified in place, so the caller's stays intact
    SEXP confidence_interval_out = PROTECT(Rf_duplicate(confidence_interval)); nprot++;

    // outputs and work space
    int ierr = 0;

    bootstrap_histogram_c(
        &n_bootstraps_v,
        &n_bins,
        &n_points,
        &n_studies,
        INTEGER(mean_pmf_counts),
        INTEGER(mean_pmf_included_n_reps),
        INTEGER(included_n_reps),
        REAL(confidence_interval_out),
        &two_sided_bootstrapping_significance_level_v,
        &random_seed_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, confidence_interval_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("confidence_interval"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_call(SEXP n_neighbors, SEXP n_bins, SEXP shared_residual_range, SEXP gene_means, SEXP gene_means_perms, SEXP residuals, SEXP x_star, SEXP n_permutations, SEXP random_seed) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[1];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[0];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(residuals, R_DimSymbol))[0];
    int n_points = (int) Rf_length(x_star);

    // scalar inputs, pulled from their length-1 vectors
    int n_neighbors_v = Rf_asInteger(n_neighbors);
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int n_permutations_v = Rf_asInteger(n_permutations);
    int random_seed_v = Rf_asInteger(random_seed);

    // outputs and work space
    SEXP neighborhood_indices = PROTECT(Rf_allocVector(INTSXP, n_neighbors_v * n_points * n_studies)); nprot++;
    { SEXP neighborhood_indices_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_indices_dim)[0] = n_neighbors_v; INTEGER(neighborhood_indices_dim)[1] = n_points; INTEGER(neighborhood_indices_dim)[2] = n_studies; Rf_setAttrib(neighborhood_indices, R_DimSymbol, neighborhood_indices_dim); UNPROTECT(1); }
    SEXP neighborhood_range = PROTECT(Rf_allocVector(INTSXP, 2 * n_points * n_studies)); nprot++;
    { SEXP neighborhood_range_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_range_dim)[0] = 2; INTEGER(neighborhood_range_dim)[1] = n_points; INTEGER(neighborhood_range_dim)[2] = n_studies; Rf_setAttrib(neighborhood_range, R_DimSymbol, neighborhood_range_dim); UNPROTECT(1); }
    SEXP pmfs = PROTECT(Rf_allocVector(REALSXP, n_bins_v * n_points * n_studies)); nprot++;
    { SEXP pmfs_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(pmfs_dim)[0] = n_bins_v; INTEGER(pmfs_dim)[1] = n_points; INTEGER(pmfs_dim)[2] = n_studies; Rf_setAttrib(pmfs, R_DimSymbol, pmfs_dim); UNPROTECT(1); }
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, n_bins_v * n_points * n_studies)); nprot++;
    { SEXP counts_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(counts_dim)[0] = n_bins_v; INTEGER(counts_dim)[1] = n_points; INTEGER(counts_dim)[2] = n_studies; Rf_setAttrib(counts, R_DimSymbol, counts_dim); UNPROTECT(1); }
    SEXP included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points * n_studies)); nprot++;
    { SEXP included_n_reps_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(included_n_reps_dim)[0] = n_points; INTEGER(included_n_reps_dim)[1] = n_studies; Rf_setAttrib(included_n_reps, R_DimSymbol, included_n_reps_dim); UNPROTECT(1); }
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_bins_v * n_points)); nprot++;
    { SEXP mean_pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_dim)[0] = n_bins_v; INTEGER(mean_pmf_dim)[1] = n_points; Rf_setAttrib(mean_pmf, R_DimSymbol, mean_pmf_dim); UNPROTECT(1); }
    SEXP mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, n_bins_v * n_points)); nprot++;
    { SEXP mean_pmf_counts_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(mean_pmf_counts_dim)[0] = n_bins_v; INTEGER(mean_pmf_counts_dim)[1] = n_points; Rf_setAttrib(mean_pmf_counts, R_DimSymbol, mean_pmf_counts_dim); UNPROTECT(1); }
    SEXP mean_pmf_included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP js_divergences_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(js_divergences_dim)[0] = n_points; INTEGER(js_divergences_dim)[1] = n_studies; Rf_setAttrib(js_divergences, R_DimSymbol, js_divergences_dim); UNPROTECT(1); }
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points * n_studies)); nprot++;
    { SEXP weights_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(weights_dim)[0] = n_points; INTEGER(weights_dim)[1] = n_studies; Rf_setAttrib(weights, R_DimSymbol, weights_dim); UNPROTECT(1); }
    SEXP global_js_divergence = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    int ierr = 0;

    run_js_comp_test_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        &n_points,
        &n_neighbors_v,
        &n_bins_v,
        &shared_residual_range_v,
        REAL(gene_means),
        INTEGER(gene_means_perms),
        REAL(residuals),
        REAL(x_star),
        INTEGER(neighborhood_indices),
        INTEGER(neighborhood_range),
        REAL(pmfs),
        INTEGER(counts),
        INTEGER(included_n_reps),
        REAL(mean_pmf),
        INTEGER(mean_pmf_counts),
        INTEGER(mean_pmf_included_n_reps),
        REAL(js_divergences),
        REAL(weights),
        REAL(global_js_divergence),
        REAL(p_values),
        &n_permutations_v,
        &random_seed_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 13)); nprot++;
    SET_VECTOR_ELT(_out, 0, neighborhood_indices);
    SET_VECTOR_ELT(_out, 1, neighborhood_range);
    SET_VECTOR_ELT(_out, 2, pmfs);
    SET_VECTOR_ELT(_out, 3, counts);
    SET_VECTOR_ELT(_out, 4, included_n_reps);
    SET_VECTOR_ELT(_out, 5, mean_pmf);
    SET_VECTOR_ELT(_out, 6, mean_pmf_counts);
    SET_VECTOR_ELT(_out, 7, mean_pmf_included_n_reps);
    SET_VECTOR_ELT(_out, 8, js_divergences);
    SET_VECTOR_ELT(_out, 9, weights);
    SET_VECTOR_ELT(_out, 10, global_js_divergence);
    SET_VECTOR_ELT(_out, 11, p_values);
    SET_VECTOR_ELT(_out, 12, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 13)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighborhood_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("neighborhood_range"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("pmfs"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("counts"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("included_n_reps"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("mean_pmf_counts"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("mean_pmf_included_n_reps"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP run_js_comp_test_parameter_search_call(SEXP gene_means, SEXP residuals, SEXP shared_residual_range, SEXP n_bootstraps, SEXP join_method, SEXP min_residuals_per_bin, SEXP min_neighbor_overlap, SEXP succeeding_ci_overlap, SEXP plateau_mode, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[1];
    int max_n_genes_all_studies = INTEGER(Rf_getAttrib(gene_means, R_DimSymbol))[0];
    int max_n_reps_all_studies = INTEGER(Rf_getAttrib(residuals, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    int min_residuals_per_bin_v = Rf_asInteger(min_residuals_per_bin);
    double min_neighbor_overlap_v = Rf_asReal(min_neighbor_overlap);
    double succeeding_ci_overlap_v = Rf_asReal(succeeding_ci_overlap);
    double delta_median_threshold_v = Rf_asReal(delta_median_threshold);
    double delta_max_threshold_v = Rf_asReal(delta_max_threshold);
    double delta_epsilon_v = Rf_asReal(delta_epsilon);
    int delta_min_consecutive_transitions_v = Rf_asInteger(delta_min_consecutive_transitions);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);
    int random_seed_v = Rf_asInteger(random_seed);

    // convert what Fortran cannot take from R directly
    char* join_method_c = tox_char_in(join_method, 11);
    char* plateau_mode_c = tox_char_in(plateau_mode, 19);

    // outputs and work space
    int n_points = 0;
    int n_neighbors = 0;
    int n_bins = 0;
    SEXP best_candidate_pair_confidence_interval = PROTECT(Rf_allocVector(REALSXP, 2 * n_studies)); nprot++;
    { SEXP best_candidate_pair_confidence_interval_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(best_candidate_pair_confidence_interval_dim)[0] = 2; INTEGER(best_candidate_pair_confidence_interval_dim)[1] = n_studies; Rf_setAttrib(best_candidate_pair_confidence_interval, R_DimSymbol, best_candidate_pair_confidence_interval_dim); UNPROTECT(1); }
    unsigned char plateau_established = 0;
    int n_admissible_evaluated = 0;
    SEXP trace_n_points = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
    SEXP trace_n_neighbors = PROTECT(Rf_allocVector(INTSXP, 16)); nprot++;
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
    SEXP trace_delta = PROTECT(Rf_allocVector(REALSXP, n_studies * 16)); nprot++;
    { SEXP trace_delta_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(trace_delta_dim)[0] = n_studies; INTEGER(trace_delta_dim)[1] = 16; Rf_setAttrib(trace_delta, R_DimSymbol, trace_delta_dim); UNPROTECT(1); }
    SEXP trace_delta_median = PROTECT(Rf_allocVector(REALSXP, 16)); nprot++;
    SEXP trace_delta_max = PROTECT(Rf_allocVector(REALSXP, 16)); nprot++;
    int ierr = 0;

    run_js_comp_test_parameter_search_c(
        &n_studies,
        &max_n_genes_all_studies,
        &max_n_reps_all_studies,
        REAL(gene_means),
        REAL(residuals),
        &shared_residual_range_v,
        &n_bootstraps_v,
        join_method_c,
        &n_points,
        &n_neighbors,
        &n_bins,
        REAL(best_candidate_pair_confidence_interval),
        &plateau_established,
        &n_admissible_evaluated,
        INTEGER(trace_n_points),
        INTEGER(trace_n_neighbors),
        REAL(trace_global_js_divergence),
        REAL(trace_ci_lower),
        REAL(trace_ci_upper),
        REAL(trace_ci_width),
        REAL(trace_ci_width_relative),
        REAL(trace_delta),
        REAL(trace_delta_median),
        REAL(trace_delta_max),
        &min_residuals_per_bin_v,
        &min_neighbor_overlap_v,
        &succeeding_ci_overlap_v,
        plateau_mode_c,
        &delta_median_threshold_v,
        &delta_max_threshold_v,
        &delta_epsilon_v,
        &delta_min_consecutive_transitions_v,
        &two_sided_bootstrapping_significance_level_v,
        &random_seed_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 17)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_points));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_neighbors));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_bins));
    SET_VECTOR_ELT(_out, 3, best_candidate_pair_confidence_interval);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarLogical(plateau_established != 0));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(n_admissible_evaluated));
    SET_VECTOR_ELT(_out, 6, trace_n_points);
    SET_VECTOR_ELT(_out, 7, trace_n_neighbors);
    SET_VECTOR_ELT(_out, 8, trace_global_js_divergence);
    SET_VECTOR_ELT(_out, 9, trace_ci_lower);
    SET_VECTOR_ELT(_out, 10, trace_ci_upper);
    SET_VECTOR_ELT(_out, 11, trace_ci_width);
    SET_VECTOR_ELT(_out, 12, trace_ci_width_relative);
    SET_VECTOR_ELT(_out, 13, trace_delta);
    SET_VECTOR_ELT(_out, 14, trace_delta_median);
    SET_VECTOR_ELT(_out, 15, trace_delta_max);
    SET_VECTOR_ELT(_out, 16, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 17)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_points"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_neighbors"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_bins"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("best_candidate_pair_confidence_interval"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("plateau_established"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("n_admissible_evaluated"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("trace_n_points"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("trace_n_neighbors"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("trace_global_js_divergence"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("trace_ci_lower"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("trace_ci_upper"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("trace_ci_width"));
    SET_STRING_ELT(_nms, 12, Rf_mkChar("trace_ci_width_relative"));
    SET_STRING_ELT(_nms, 13, Rf_mkChar("trace_delta"));
    SET_STRING_ELT(_nms, 14, Rf_mkChar("trace_delta_median"));
    SET_STRING_ELT(_nms, 15, Rf_mkChar("trace_delta_max"));
    SET_STRING_ELT(_nms, 16, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
