// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void compute_edf_c(const double*, const int*, double*, double*, int*, int*);
void compute_edf_expert_c(const double*, const int*, const int*, double*, double*, int*, int*);
void calc_percentile_c(const double*, const int*, const double*, double*, const int*, int*);
void calc_percentile_expert_c(const double*, const int*, const int*, const double*, double*, const int*, int*);

SEXP compute_edf_call(SEXP values) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_values = (int) Rf_length(values);

    // outputs and work space
    SEXP unique_values = PROTECT(Rf_allocVector(REALSXP, n_values)); nprot++;
    SEXP cdf_values = PROTECT(Rf_allocVector(REALSXP, n_values)); nprot++;
    int n_unique = 0;
    int ierr = 0;

    compute_edf_c(
        REAL(values),
        &n_values,
        REAL(unique_values),
        REAL(cdf_values),
        &n_unique,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, unique_values);
    SET_VECTOR_ELT(_out, 1, cdf_values);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_unique));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("unique_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("cdf_values"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_unique"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_edf_expert_call(SEXP values, SEXP values_perm) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_values = (int) Rf_length(values);

    // outputs and work space
    SEXP unique_values = PROTECT(Rf_allocVector(REALSXP, n_values)); nprot++;
    SEXP cdf_values = PROTECT(Rf_allocVector(REALSXP, n_values)); nprot++;
    int n_unique = 0;
    int ierr = 0;

    compute_edf_expert_c(
        REAL(values),
        &n_values,
        INTEGER(values_perm),
        REAL(unique_values),
        REAL(cdf_values),
        &n_unique,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, unique_values);
    SET_VECTOR_ELT(_out, 1, cdf_values);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_unique));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("unique_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("cdf_values"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_unique"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP calc_percentile_call(SEXP array, SEXP percentile, SEXP n_considered) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_array = (int) Rf_length(array);

    // scalar inputs, pulled from their length-1 vectors
    double percentile_v = Rf_asReal(percentile);
    int n_considered_v = Rf_asInteger(n_considered);

    // outputs and work space
    double value = 0;
    int ierr = 0;

    calc_percentile_c(
        REAL(array),
        &n_array,
        &percentile_v,
        &value,
        &n_considered_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(value));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("value"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP calc_percentile_expert_call(SEXP array, SEXP array_perm, SEXP percentile, SEXP n_considered) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_array = (int) Rf_length(array);

    // scalar inputs, pulled from their length-1 vectors
    double percentile_v = Rf_asReal(percentile);
    int n_considered_v = Rf_asInteger(n_considered);

    // outputs and work space
    double value = 0;
    int ierr = 0;

    calc_percentile_expert_c(
        REAL(array),
        &n_array,
        INTEGER(array_perm),
        &percentile_v,
        &value,
        &n_considered_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(value));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("value"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
