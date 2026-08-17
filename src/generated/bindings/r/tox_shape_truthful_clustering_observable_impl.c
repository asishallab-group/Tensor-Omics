// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void tox_stc_observable_svd_workspace_c(const int*, const int*, int*, int*, int*);

SEXP tox_stc_observable_svd_workspace_call(SEXP n_dimensions, SEXP n_selected_member) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_dimensions_v = Rf_asInteger(n_dimensions);
    int n_selected_member_v = Rf_asInteger(n_selected_member);

    // outputs and work space
    int lwork = 0;
    int iwork_size = 0;
    int ierr = 0;

    tox_stc_observable_svd_workspace_c(
        &n_dimensions_v,
        &n_selected_member_v,
        &lwork,
        &iwork_size,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(lwork));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(iwork_size));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("lwork"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("iwork_size"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
