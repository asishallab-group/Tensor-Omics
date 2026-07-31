// Generated. Do not edit.
#if !defined(NO_R_INTERFACE) && !defined(NO_C_INTERFACE)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void build_bst_index_c(const double*, const int*, int*, int*);
void bst_range_query_c(const double*, const int*, const int*, const double*, const double*, int*, int*, int*);

SEXP build_bst_index_call(SEXP values) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_values = (int) Rf_length(values);

    // outputs and work space
    SEXP sorted_indices = PROTECT(Rf_allocVector(INTSXP, n_values)); nprot++;
    int ierr = 0;

    build_bst_index_c(
        REAL(values),
        &n_values,
        INTEGER(sorted_indices),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, sorted_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("sorted_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bst_range_query_call(SEXP values, SEXP sorted_indices, SEXP lower_bound, SEXP upper_bound) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_values = (int) Rf_length(values);

    // scalar inputs, pulled from their length-1 vectors
    double lower_bound_v = Rf_asReal(lower_bound);
    double upper_bound_v = Rf_asReal(upper_bound);

    // outputs and work space
    SEXP output_indices = PROTECT(Rf_allocVector(INTSXP, n_values)); nprot++;
    int n_matches = 0;
    int ierr = 0;

    bst_range_query_c(
        REAL(values),
        INTEGER(sorted_indices),
        &n_values,
        &lower_bound_v,
        &upper_bound_v,
        INTEGER(output_indices),
        &n_matches,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, output_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_matches));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("output_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_matches"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R interface
