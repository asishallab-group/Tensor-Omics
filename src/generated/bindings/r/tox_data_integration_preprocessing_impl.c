// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void calc_neighborhood_size_c(const int*, const int*, const int*, const double*, const int*, int*, int*);

SEXP calc_neighborhood_size_call(SEXP n_pool, SEXP n_points, SEXP mean_S, SEXP desired_size) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S = (int) Rf_length(mean_S);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_pool = Rf_asInteger(n_pool);
    int _v_n_points = Rf_asInteger(n_points);
    int _v_desired_size = Rf_asInteger(desired_size);

    // outputs and work space
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(
        &_v_n_pool,
        &_v_n_points,
        &n_genes_S,
        REAL(mean_S),
        &_v_desired_size,
        &n_neighbors,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_neighbors));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
