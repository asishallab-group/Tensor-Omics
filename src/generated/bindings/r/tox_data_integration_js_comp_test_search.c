// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void check_neighborhood_overlaps_c(const int*, const int*, const double*, unsigned char*, int*);
void check_mean_pmf_min_counts_c(const int*, const int*, const int*, const int*, const int*, unsigned char*, int*);
void check_plateau_condition_c(const double*, double*, const int*, int*, int*, const int*, const char*, const double*, double*, unsigned char*, int*);
void check_effect_size_plateau_condition_c(const double*, const double*, const int*, const unsigned char*, const double*, const double*, const double*, const int*, int*, double*, double*, double*, unsigned char*, int*);

SEXP check_neighborhood_overlaps_call(SEXP neighborhood_range, SEXP min_neighbor_overlap) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(neighborhood_range, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_min_neighbor_overlap = Rf_asReal(min_neighbor_overlap);

    // convert what Fortran cannot take from R directly
    const int* _c_neighborhood_range = tox_int_in(neighborhood_range, "neighborhood_range");

    // outputs and work space
    unsigned char all_have_min_neighbor_overlap = 0;
    int ierr = 0;

    check_neighborhood_overlaps_c(
        _c_neighborhood_range,
        &n_points,
        &_v_min_neighbor_overlap,
        &all_have_min_neighbor_overlap,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(all_have_min_neighbor_overlap != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("all_have_min_neighbor_overlap"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP check_mean_pmf_min_counts_call(SEXP mean_pmf_counts, SEXP n_bins_per_point, SEXP min_count) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_min_count = Rf_asInteger(min_count);

    // convert what Fortran cannot take from R directly
    const int* _c_mean_pmf_counts = tox_int_in(mean_pmf_counts, "mean_pmf_counts");
    const int* _c_n_bins_per_point = tox_int_in(n_bins_per_point, "n_bins_per_point");

    // outputs and work space
    unsigned char all_bins_have_min_count = 0;
    int ierr = 0;

    check_mean_pmf_min_counts_c(
        _c_mean_pmf_counts,
        &n_bins,
        _c_n_bins_per_point,
        &n_points,
        &_v_min_count,
        &all_bins_have_min_count,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(all_bins_have_min_count != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("all_bins_have_min_count"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP check_plateau_condition_call(SEXP confidence_interval, SEXP best_candidate_pair_confidence_interval, SEXP best_candidate_index, SEXP best_exceeded_ci_overlap_count, SEXP candidate_index, SEXP join_method, SEXP succeeding_ci_overlap) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = INTEGER(Rf_getAttrib(confidence_interval, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_best_candidate_index = Rf_asInteger(best_candidate_index);
    int _v_best_exceeded_ci_overlap_count = Rf_asInteger(best_exceeded_ci_overlap_count);
    int _v_candidate_index = Rf_asInteger(candidate_index);
    double _v_succeeding_ci_overlap = Rf_asReal(succeeding_ci_overlap);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_best_candidate_pair_confidence_interval = PROTECT(Rf_duplicate(best_candidate_pair_confidence_interval)); _nprot++;

    // convert what Fortran cannot take from R directly
    char* _c_join_method = tox_char_in(join_method, 11, "join_method");

    // outputs and work space
    SEXP ci_overlap = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    unsigned char plateau_found = 0;
    int ierr = 0;

    check_plateau_condition_c(
        REAL(confidence_interval),
        REAL(_out_best_candidate_pair_confidence_interval),
        &n_studies,
        &_v_best_candidate_index,
        &_v_best_exceeded_ci_overlap_count,
        &_v_candidate_index,
        _c_join_method,
        &_v_succeeding_ci_overlap,
        REAL(ci_overlap),
        &plateau_found,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_best_candidate_pair_confidence_interval);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(_v_best_candidate_index));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(_v_best_exceeded_ci_overlap_count));
    SET_VECTOR_ELT(_out, 3, ci_overlap);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarLogical(plateau_found != 0));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("best_candidate_pair_confidence_interval"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("best_candidate_index"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("best_exceeded_ci_overlap_count"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ci_overlap"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("plateau_found"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP check_effect_size_plateau_condition_call(SEXP global_js_divergence, SEXP prev_global_js_divergence, SEXP has_previous, SEXP delta_median_threshold, SEXP delta_max_threshold, SEXP delta_epsilon, SEXP delta_min_consecutive_transitions, SEXP n_consecutive_ok) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_studies = (int) Rf_length(global_js_divergence);

    // scalar inputs, pulled from their length-1 vectors
    unsigned char _v_has_previous = (Rf_asLogical(has_previous) == TRUE) ? 1 : 0;
    double _v_delta_median_threshold = Rf_asReal(delta_median_threshold);
    double _v_delta_max_threshold = Rf_asReal(delta_max_threshold);
    double _v_delta_epsilon = Rf_asReal(delta_epsilon);
    int _v_delta_min_consecutive_transitions = Rf_asInteger(delta_min_consecutive_transitions);
    int _v_n_consecutive_ok = Rf_asInteger(n_consecutive_ok);

    // outputs and work space
    SEXP delta = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    double delta_median = 0;
    double delta_max = 0;
    unsigned char plateau_found = 0;
    int ierr = 0;

    check_effect_size_plateau_condition_c(
        REAL(global_js_divergence),
        REAL(prev_global_js_divergence),
        &n_studies,
        &_v_has_previous,
        &_v_delta_median_threshold,
        &_v_delta_max_threshold,
        &_v_delta_epsilon,
        &_v_delta_min_consecutive_transitions,
        &_v_n_consecutive_ok,
        REAL(delta),
        &delta_median,
        &delta_max,
        &plateau_found,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(_v_n_consecutive_ok));
    SET_VECTOR_ELT(_out, 1, delta);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(delta_median));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(delta_max));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarLogical(plateau_found != 0));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_consecutive_ok"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("delta"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("delta_median"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("delta_max"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("plateau_found"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
