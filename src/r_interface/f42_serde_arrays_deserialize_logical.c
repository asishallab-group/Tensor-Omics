// Generated. Do not edit.
#if !defined(NO_R_INTERFACE) && !defined(NO_C_INTERFACE)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void deserialize_logical_helper_c(unsigned char*, const int*, const int*, const int*, const char*, const int*, int*);

SEXP deserialize_logical_helper_call(SEXP arr_shape, SEXP filename) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_elements = tox_prod(arr_shape);
    int n_arr_shape_elements = (int) Rf_length(arr_shape);
    int filename_strlen = tox_max_strlen(filename);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);

    // outputs and work space
    unsigned char* arr_c = tox_bool_alloc(tox_prod(arr_shape));
    int ierr = 0;

    deserialize_logical_helper_c(
        arr_c,
        &n_elements,
        INTEGER(arr_shape),
        &n_arr_shape_elements,
        filename_c,
        &filename_strlen,
        &ierr
    );

    // convert the outputs back
    SEXP arr = PROTECT(tox_bool_out(arr_c, tox_prod(arr_shape))); nprot++;
    Rf_setAttrib(arr, R_DimSymbol, arr_shape);

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, arr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("arr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R interface
