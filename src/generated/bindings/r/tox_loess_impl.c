// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void tox_loess_required_workspace_c(const int*, const int*, int*, int*, const unsigned char*, int*);

SEXP tox_loess_required_workspace_call(SEXP n_dim, SEXP max_neighborhood_size, SEXP save_factorization) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_n_dim = Rf_asInteger(n_dim);
    int _v_max_neighborhood_size = Rf_asInteger(max_neighborhood_size);
    unsigned char _v_save_factorization = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;

    // outputs and work space
    int int_workspace_size = 0;
    int real_workspace_size = 0;
    int ierr = 0;

    tox_loess_required_workspace_c(
        &_v_n_dim,
        &_v_max_neighborhood_size,
        &int_workspace_size,
        &real_workspace_size,
        &_v_save_factorization,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(int_workspace_size));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(real_workspace_size));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("int_workspace_size"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("real_workspace_size"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
