// Generated. Do not edit.
#if !defined(NO_R_INTERFACE) && !defined(NO_C_INTERFACE)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void get_array_metadata_c(const char*, const int*, int*, const int*, int*, int*, int*);

SEXP get_array_metadata_call(SEXP filename, SEXP dims_out_capacity) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);

    // scalar inputs, pulled from their length-1 vectors
    int dims_out_capacity_v = Rf_asInteger(dims_out_capacity);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);

    // outputs and work space
    SEXP dims_out = PROTECT(Rf_allocVector(INTSXP, dims_out_capacity_v)); nprot++;
    int ndims = 0;
    int type_code = 0;
    int ierr = 0;

    get_array_metadata_c(
        filename_c,
        &filename_strlen,
        INTEGER(dims_out),
        &dims_out_capacity_v,
        &ndims,
        &type_code,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, dims_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ndims));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(type_code));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("dims_out"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ndims"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("type_code"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R interface
