// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void check_neighborhood_overlaps_c(const int*, const int*, const double*, unsigned char*, int*);
void check_mean_pmf_min_counts_c(const int*, const int*, const int*, const int*, const int*, unsigned char*, int*);
void check_plateau_condition_c(const double*, double*, const int*, int*, int*, const int*, const char*, const double*, double*, unsigned char*, int*);
void check_effect_size_plateau_condition_c(const double*, const double*, const int*, const unsigned char*, const double*, const double*, const double*, const int*, int*, double*, double*, double*, unsigned char*, int*);

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
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];

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
    SEXP ci_overlap = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
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
        REAL(ci_overlap),
        &plateau_found,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
    SET_VECTOR_ELT(_out, 0, best_candidate_pair_confidence_interval_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(best_candidate_index_v));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(best_exceeded_ci_overlap_count_v));
    SET_VECTOR_ELT(_out, 3, ci_overlap);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarLogical(plateau_found != 0));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("best_candidate_pair_confidence_interval"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("best_candidate_index"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("best_exceeded_ci_overlap_count"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ci_overlap"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("plateau_found"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
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

#endif  // R binding
