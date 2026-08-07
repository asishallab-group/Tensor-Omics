// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void loess_smooth_2d_c(const int*, const int*, const double*, const double*, const int*, const int*, const double*, const double*, const double*, double*, int*);
void compute_scaled_distance_quantile_c(const int*, const double*, const double*, const int*, double*, const double*, int*);

SEXP loess_smooth_2d_call(SEXP x_ref, SEXP y_ref, SEXP indices_used, SEXP x_query, SEXP kernel_sigma, SEXP kernel_cutoff) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_total = (int) Rf_length(x_ref);
    int n_target = (int) Rf_length(x_query);
    int n_used = (int) Rf_length(indices_used);

    // scalar inputs, pulled from their length-1 vectors
    double kernel_sigma_v = Rf_asReal(kernel_sigma);
    double kernel_cutoff_v = Rf_asReal(kernel_cutoff);

    // outputs and work space
    SEXP y_out = PROTECT(Rf_allocVector(REALSXP, n_target)); nprot++;
    int ierr = 0;

    loess_smooth_2d_c(
        &n_total,
        &n_target,
        REAL(x_ref),
        REAL(y_ref),
        INTEGER(indices_used),
        &n_used,
        REAL(x_query),
        &kernel_sigma_v,
        &kernel_cutoff_v,
        REAL(y_out),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, y_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("y_out"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_scaled_distance_quantile_call(SEXP rdi, SEXP sorted_rdi, SEXP perm, SEXP c_const) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(rdi);

    // scalar inputs, pulled from their length-1 vectors
    double c_const_v = Rf_asReal(c_const);

    // outputs and work space
    SEXP quantile = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_scaled_distance_quantile_c(
        &n_genes,
        REAL(rdi),
        REAL(sorted_rdi),
        INTEGER(perm),
        REAL(quantile),
        &c_const_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, quantile);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("quantile"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
