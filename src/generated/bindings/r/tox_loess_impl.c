// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void tox_loess_required_workspace_c(const int*, const int*, int*, int*, const unsigned char*, int*);

SEXP tox_loess_required_workspace_call(SEXP n_dim, SEXP max_neighborhood_size, SEXP save_factorization) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_dim_v = Rf_asInteger(n_dim);
    int max_neighborhood_size_v = Rf_asInteger(max_neighborhood_size);
    unsigned char save_factorization_v = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;

    // outputs and work space
    int int_workspace_size = 0;
    int real_workspace_size = 0;
    int ierr = 0;

    tox_loess_required_workspace_c(
        &n_dim_v,
        &max_neighborhood_size_v,
        &int_workspace_size,
        &real_workspace_size,
        &save_factorization_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(int_workspace_size));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(real_workspace_size));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("int_workspace_size"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("real_workspace_size"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
