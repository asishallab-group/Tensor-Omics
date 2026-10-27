// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void deserialize_complex_helper_c(double _Complex*, const int*, const int*, const int*, const char*, const int*, int*);

SEXP deserialize_complex_helper_call(SEXP arr_shape, SEXP filename) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_elements = tox_prod(arr_shape);
    int n_arr_shape_elements = (int) Rf_length(arr_shape);
    int filename_strlen = tox_max_strlen(filename);

    // convert what Fortran cannot take from R directly
    char* _c_filename = tox_char_in(filename, filename_strlen, "filename");

    // outputs and work space
    SEXP arr = PROTECT(Rf_allocVector(CPLXSXP, tox_prod(arr_shape))); _nprot++;
    Rf_setAttrib(arr, R_DimSymbol, arr_shape);
    int ierr = 0;

    deserialize_complex_helper_c(
        (double _Complex*) COMPLEX(arr),
        &n_elements,
        INTEGER(arr_shape),
        &n_arr_shape_elements,
        _c_filename,
        &filename_strlen,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, arr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("arr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
