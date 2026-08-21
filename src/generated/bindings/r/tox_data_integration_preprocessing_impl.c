// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void calc_neighborhood_size_c(const int*, const int*, const int*, const double*, const int*, int*, int*);

SEXP calc_neighborhood_size_call(SEXP n_pool, SEXP n_points, SEXP mean_S, SEXP desired_size) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S = (int) Rf_length(mean_S);

    // scalar inputs, pulled from their length-1 vectors
    int n_pool_v = Rf_asInteger(n_pool);
    int n_points_v = Rf_asInteger(n_points);
    int desired_size_v = Rf_asInteger(desired_size);

    // outputs and work space
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(
        &n_pool_v,
        &n_points_v,
        &n_genes_S,
        REAL(mean_S),
        &desired_size_v,
        &n_neighbors,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_neighbors));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
