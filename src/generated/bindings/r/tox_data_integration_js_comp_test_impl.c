// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void calc_js_comp_test_candidate_bounds_c(const int*, int*, int*, int*);
void calc_js_comp_test_n_top_k_jsds_c(const int*, const double*, int*, int*);

SEXP calc_js_comp_test_candidate_bounds_call(SEXP max_n_genes_all_studies) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int max_n_genes_all_studies_v = Rf_asInteger(max_n_genes_all_studies);

    // outputs and work space
    int max_n_points_candidate = 0;
    int max_n_neighbors_candidate = 0;
    int ierr = 0;

    calc_js_comp_test_candidate_bounds_c(
        &max_n_genes_all_studies_v,
        &max_n_points_candidate,
        &max_n_neighbors_candidate,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(max_n_points_candidate));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(max_n_neighbors_candidate));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("max_n_points_candidate"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("max_n_neighbors_candidate"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP calc_js_comp_test_n_top_k_jsds_call(SEXP n_bootstraps, SEXP two_sided_bootstrapping_significance_level) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_bootstraps_v = Rf_asInteger(n_bootstraps);
    double two_sided_bootstrapping_significance_level_v = Rf_asReal(two_sided_bootstrapping_significance_level);

    // outputs and work space
    int n_top_k = 0;
    int ierr = 0;

    calc_js_comp_test_n_top_k_jsds_c(
        &n_bootstraps_v,
        &two_sided_bootstrapping_significance_level_v,
        &n_top_k,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_top_k));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_top_k"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
