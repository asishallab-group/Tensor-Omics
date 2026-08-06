// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void loess_fit_plain_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, double*, int*);
void loess_fit_robust_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, double*, int*);

SEXP loess_fit_plain_call(SEXP x, SEXP y, SEXP weights, SEXP eval_points, SEXP span, SEXP degree, SEXP max_neighborhood_size, SEXP compute_influence, SEXP save_factorization) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n = (int) Rf_length(x);

    // scalar inputs, pulled from their length-1 vectors
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    int max_neighborhood_size_v = Rf_asInteger(max_neighborhood_size);
    unsigned char compute_influence_v = (Rf_asLogical(compute_influence) == TRUE) ? 1 : 0;
    unsigned char save_factorization_v = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;

    // outputs and work space
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, n)); nprot++;
    int ierr = 0;

    loess_fit_plain_c(
        &n,
        REAL(x),
        REAL(y),
        REAL(weights),
        REAL(eval_points),
        &span_v,
        &degree_v,
        &max_neighborhood_size_v,
        &compute_influence_v,
        &save_factorization_v,
        REAL(fitted_values),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, fitted_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("fitted_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP loess_fit_robust_call(SEXP x, SEXP y, SEXP weights, SEXP eval_points, SEXP span, SEXP degree, SEXP max_neighborhood_size, SEXP compute_influence, SEXP save_factorization, SEXP n_iters) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n = (int) Rf_length(x);

    // scalar inputs, pulled from their length-1 vectors
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    int max_neighborhood_size_v = Rf_asInteger(max_neighborhood_size);
    unsigned char compute_influence_v = (Rf_asLogical(compute_influence) == TRUE) ? 1 : 0;
    unsigned char save_factorization_v = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;
    int n_iters_v = Rf_asInteger(n_iters);

    // outputs and work space
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, n)); nprot++;
    int ierr = 0;

    loess_fit_robust_c(
        &n,
        REAL(x),
        REAL(y),
        REAL(weights),
        REAL(eval_points),
        &span_v,
        &degree_v,
        &max_neighborhood_size_v,
        &compute_influence_v,
        &save_factorization_v,
        &n_iters_v,
        REAL(fitted_values),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, fitted_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("fitted_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
