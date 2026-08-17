// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void tox_stc_accept_ensemble_svd_workspace_c(const int*, int*, int*);

SEXP tox_stc_accept_ensemble_svd_workspace_call(SEXP n_dimensions) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_dimensions_v = Rf_asInteger(n_dimensions);

    // outputs and work space
    int lwork = 0;
    int ierr = 0;

    tox_stc_accept_ensemble_svd_workspace_c(
        &n_dimensions_v,
        &lwork,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(lwork));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("lwork"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
