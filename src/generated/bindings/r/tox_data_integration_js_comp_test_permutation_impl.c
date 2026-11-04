// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void calc_js_comp_test_permutation_bounds_c(const int*, const int*, const int*, const int*, const double*, int*, int*, int*, int*);
void calc_adaptive_js_comp_test_permutation_bounds_c(const int*, const int*, const int*, const int*, const double*, int*, int*, int*);

SEXP calc_js_comp_test_permutation_bounds_call(SEXP max_n_genes_all_studies, SEXP n_points, SEXP n_neighbors, SEXP n_bootstraps, SEXP two_sided_bootstrapping_significance_level) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_genes_all_studies = Rf_asInteger(max_n_genes_all_studies);
    int _v_n_points = Rf_asInteger(n_points);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    int _v_n_bootstraps = Rf_asInteger(n_bootstraps);
    double _v_two_sided_bootstrapping_significance_level = Rf_asReal(two_sided_bootstrapping_significance_level);

    // outputs and work space
    int max_n_points_candidate = 0;
    int max_n_neighbors_candidate = 0;
    int n_bootstrapping_top_k_jsds = 0;
    int ierr = 0;

    calc_js_comp_test_permutation_bounds_c(
        &_v_max_n_genes_all_studies,
        &_v_n_points,
        &_v_n_neighbors,
        &_v_n_bootstraps,
        &_v_two_sided_bootstrapping_significance_level,
        &max_n_points_candidate,
        &max_n_neighbors_candidate,
        &n_bootstrapping_top_k_jsds,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(max_n_points_candidate));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(max_n_neighbors_candidate));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_bootstrapping_top_k_jsds));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("max_n_points_candidate"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("max_n_neighbors_candidate"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_bootstrapping_top_k_jsds"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_adaptive_js_comp_test_permutation_bounds_call(SEXP max_n_genes_all_studies, SEXP n_studies, SEXP k_start, SEXP n_bootstraps, SEXP two_sided_bootstrapping_significance_level) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_genes_all_studies = Rf_asInteger(max_n_genes_all_studies);
    int _v_n_studies = Rf_asInteger(n_studies);
    int _v_k_start = Rf_asInteger(k_start);
    int _v_n_bootstraps = Rf_asInteger(n_bootstraps);
    double _v_two_sided_bootstrapping_significance_level = Rf_asReal(two_sided_bootstrapping_significance_level);

    // outputs and work space
    int max_n_points_candidate = 0;
    int n_bootstrapping_top_k_jsds = 0;
    int ierr = 0;

    calc_adaptive_js_comp_test_permutation_bounds_c(
        &_v_max_n_genes_all_studies,
        &_v_n_studies,
        &_v_k_start,
        &_v_n_bootstraps,
        &_v_two_sided_bootstrapping_significance_level,
        &max_n_points_candidate,
        &n_bootstrapping_top_k_jsds,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(max_n_points_candidate));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_bootstrapping_top_k_jsds));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("max_n_points_candidate"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_bootstrapping_top_k_jsds"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
