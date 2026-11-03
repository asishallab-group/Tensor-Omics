// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void loess_fit_plain_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, double*, int*);
void loess_fit_robust_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, double*, int*);

SEXP loess_fit_plain_call(SEXP x, SEXP y, SEXP weights, SEXP eval_points, SEXP span, SEXP degree, SEXP max_neighborhood_size, SEXP compute_influence, SEXP save_factorization) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n = (int) Rf_length(x);

    // scalar inputs, pulled from their length-1 vectors
    double _v_span = Rf_asReal(span);
    int _v_degree = Rf_asInteger(degree);
    int _v_max_neighborhood_size = Rf_asInteger(max_neighborhood_size);
    unsigned char _v_compute_influence = (Rf_asLogical(compute_influence) == TRUE) ? 1 : 0;
    unsigned char _v_save_factorization = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;

    // outputs and work space
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, n)); _nprot++;
    int ierr = 0;

    loess_fit_plain_c(
        &n,
        REAL(x),
        REAL(y),
        REAL(weights),
        REAL(eval_points),
        &_v_span,
        &_v_degree,
        &_v_max_neighborhood_size,
        &_v_compute_influence,
        &_v_save_factorization,
        REAL(fitted_values),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, fitted_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("fitted_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP loess_fit_robust_call(SEXP x, SEXP y, SEXP weights, SEXP eval_points, SEXP span, SEXP degree, SEXP max_neighborhood_size, SEXP compute_influence, SEXP save_factorization, SEXP n_iters) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n = (int) Rf_length(x);

    // scalar inputs, pulled from their length-1 vectors
    double _v_span = Rf_asReal(span);
    int _v_degree = Rf_asInteger(degree);
    int _v_max_neighborhood_size = Rf_asInteger(max_neighborhood_size);
    unsigned char _v_compute_influence = (Rf_asLogical(compute_influence) == TRUE) ? 1 : 0;
    unsigned char _v_save_factorization = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;
    int _v_n_iters = Rf_asInteger(n_iters);

    // outputs and work space
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, n)); _nprot++;
    int ierr = 0;

    loess_fit_robust_c(
        &n,
        REAL(x),
        REAL(y),
        REAL(weights),
        REAL(eval_points),
        &_v_span,
        &_v_degree,
        &_v_max_neighborhood_size,
        &_v_compute_influence,
        &_v_save_factorization,
        &_v_n_iters,
        REAL(fitted_values),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, fitted_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("fitted_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
