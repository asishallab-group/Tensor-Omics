// Generated. Do not edit.
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void serialize_char_helper_c(const char*, const int*, const int*, const int*, const int*, const char*, const int*, int*);

SEXP serialize_char_helper_call(SEXP arr, SEXP filename) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    SEXP arr_shape = PROTECT(tox_shape_of(arr)); nprot++;
    int arr_strlen = tox_max_strlen(arr);
    int n_strings = (int) Rf_length(arr);
    int n_arr_shape_elements = (int) Rf_length(arr_shape);
    int filename_strlen = tox_max_strlen(filename);

    // convert what Fortran cannot take from R directly
    char* arr_c = tox_char_in(arr, arr_strlen);
    char* filename_c = tox_char_in(filename, filename_strlen);

    // outputs and work space
    int ierr = 0;

    serialize_char_helper_c(
        arr_c,
        &arr_strlen,
        &n_strings,
        INTEGER(arr_shape),
        &n_arr_shape_elements,
        filename_c,
        &filename_strlen,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}
