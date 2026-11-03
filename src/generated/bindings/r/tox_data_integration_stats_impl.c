// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void calc_js_comp_test_n_top_k_jsds_c(const int*, const double*, int*, int*);

SEXP calc_js_comp_test_n_top_k_jsds_call(SEXP n_bootstraps, SEXP two_sided_bootstrapping_significance_level) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bootstraps = Rf_asInteger(n_bootstraps);
    double _v_two_sided_bootstrapping_significance_level = Rf_asReal(two_sided_bootstrapping_significance_level);

    // outputs and work space
    int n_top_k = 0;
    int ierr = 0;

    calc_js_comp_test_n_top_k_jsds_c(
        &_v_n_bootstraps,
        &_v_two_sided_bootstrapping_significance_level,
        &n_top_k,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_top_k));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_top_k"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
