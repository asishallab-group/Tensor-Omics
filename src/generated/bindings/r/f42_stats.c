// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void loess_smooth_2d_c(const int*, const int*, const double*, const double*, const int*, const int*, const double*, const double*, const double*, double*, int*);
void compute_edf_c(const double*, const int*, double*, double*, int*, int*);
void compute_edf_expert_c(const double*, const int*, const int*, double*, double*, int*, int*);
void calc_quantile_c(const double*, const int*, const double*, double*, const int*, int*);
void calc_quantile_expert_c(const double*, const int*, const int*, const double*, double*, const int*, int*);
void compute_scaled_distance_tail_probability_c(const int*, const double*, const double*, double*, const double*, int*);
void compute_scaled_distance_tail_probability_expert_c(const int*, const double*, const double*, const int*, double*, const double*, int*);

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

SEXP calc_quantile_call(SEXP array, SEXP level, SEXP n_considered) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_array = (int) Rf_length(array);

    // scalar inputs, pulled from their length-1 vectors
    double level_v = Rf_asReal(level);
    int n_considered_v = Rf_asInteger(n_considered);

    // outputs and work space
    double value = 0;
    int ierr = 0;

    calc_quantile_c(
        REAL(array),
        &n_array,
        &level_v,
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

SEXP calc_quantile_expert_call(SEXP array, SEXP array_perm, SEXP level, SEXP n_considered) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_array = (int) Rf_length(array);

    // scalar inputs, pulled from their length-1 vectors
    double level_v = Rf_asReal(level);
    int n_considered_v = Rf_asInteger(n_considered);

    // outputs and work space
    double value = 0;
    int ierr = 0;

    calc_quantile_expert_c(
        REAL(array),
        &n_array,
        INTEGER(array_perm),
        &level_v,
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

SEXP compute_scaled_distance_tail_probability_call(SEXP rdi, SEXP sorted_rdi, SEXP c_const) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(rdi);

    // scalar inputs, pulled from their length-1 vectors
    double c_const_v = Rf_asReal(c_const);

    // outputs and work space
    SEXP tail_probability = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_scaled_distance_tail_probability_c(
        &n_genes,
        REAL(rdi),
        REAL(sorted_rdi),
        REAL(tail_probability),
        &c_const_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, tail_probability);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("tail_probability"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_scaled_distance_tail_probability_expert_call(SEXP rdi, SEXP sorted_rdi, SEXP sorted_rdi_perm, SEXP c_const) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(rdi);

    // scalar inputs, pulled from their length-1 vectors
    double c_const_v = Rf_asReal(c_const);

    // outputs and work space
    SEXP tail_probability = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_scaled_distance_tail_probability_expert_c(
        &n_genes,
        REAL(rdi),
        REAL(sorted_rdi),
        INTEGER(sorted_rdi_perm),
        REAL(tail_probability),
        &c_const_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, tail_probability);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("tail_probability"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
