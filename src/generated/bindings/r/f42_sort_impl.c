// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void sort_real_get_perm_c(const double*, const int*, int*, int*);

SEXP sort_real_get_perm_call(SEXP array) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n = (int) Rf_length(array);

    // outputs and work space
    SEXP perm = PROTECT(Rf_allocVector(INTSXP, n)); nprot++;
    int ierr = 0;

    sort_real_get_perm_c(
        REAL(array),
        &n,
        INTEGER(perm),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, perm);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("perm"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
