// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void loess_fit_plain_expert_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, int*, const int*, double*, const int*, double*, double*, int*);
void loess_fit_plain_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, double*, int*);
void loess_fit_robust_expert_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, int*, const int*, double*, const int*, double*, double*, double*, double*, int*, double*, int*);
void loess_fit_robust_c(const int*, const double*, const double*, const double*, const double*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, double*, int*);
void loess_c(const double*, const int*, const double*, const int*, const double*, const int*, double*, const char*, const int*, int*);

SEXP loess_fit_plain_expert_call(SEXP x, SEXP y, SEXP weights, SEXP eval_points, SEXP span, SEXP degree, SEXP max_neighborhood_size, SEXP compute_influence, SEXP save_factorization, SEXP int_workspace_size, SEXP real_workspace_size) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n = (int) Rf_length(x);

    // scalar inputs, pulled from their length-1 vectors
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    int max_neighborhood_size_v = Rf_asInteger(max_neighborhood_size);
    unsigned char compute_influence_v = (Rf_asLogical(compute_influence) == TRUE) ? 1 : 0;
    unsigned char save_factorization_v = (Rf_asLogical(save_factorization) == TRUE) ? 1 : 0;
    int int_workspace_size_v = Rf_asInteger(int_workspace_size);
    int real_workspace_size_v = Rf_asInteger(real_workspace_size);

    // outputs and work space
    int* tmp_int_workspace = (int*) R_alloc(int_workspace_size_v, sizeof(int));
    double* tmp_real_workspace = (double*) R_alloc(real_workspace_size_v, sizeof(double));
    double* tmp_hat_diag = (double*) R_alloc(n, sizeof(double));
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, n)); nprot++;
    int ierr = 0;

    loess_fit_plain_expert_c(
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
        tmp_int_workspace,
        &int_workspace_size_v,
        tmp_real_workspace,
        &real_workspace_size_v,
        tmp_hat_diag,
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

SEXP loess_fit_robust_expert_call(SEXP x, SEXP y, SEXP weights, SEXP eval_points, SEXP span, SEXP degree, SEXP max_neighborhood_size, SEXP compute_influence, SEXP save_factorization, SEXP n_iters, SEXP int_workspace_size, SEXP real_workspace_size) {
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
    int int_workspace_size_v = Rf_asInteger(int_workspace_size);
    int real_workspace_size_v = Rf_asInteger(real_workspace_size);

    // outputs and work space
    int* tmp_int_workspace = (int*) R_alloc(int_workspace_size_v, sizeof(int));
    double* tmp_real_workspace = (double*) R_alloc(real_workspace_size_v, sizeof(double));
    double* tmp_hat_diag = (double*) R_alloc(n, sizeof(double));
    double* tmp_robust_weights = (double*) R_alloc(n, sizeof(double));
    double* tmp_combined_weights = (double*) R_alloc(n, sizeof(double));
    double* tmp_residuals = (double*) R_alloc(n, sizeof(double));
    int* tmp_permutation_indices = (int*) R_alloc(n, sizeof(int));
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, n)); nprot++;
    int ierr = 0;

    loess_fit_robust_expert_c(
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
        tmp_int_workspace,
        &int_workspace_size_v,
        tmp_real_workspace,
        &real_workspace_size_v,
        tmp_hat_diag,
        tmp_robust_weights,
        tmp_combined_weights,
        tmp_residuals,
        tmp_permutation_indices,
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

SEXP loess_call(SEXP x, SEXP y, SEXP span, SEXP degree, SEXP mode, SEXP n_iters) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_x_elements = (int) Rf_length(x);
    int n_y_elements = (int) Rf_length(y);

    // scalar inputs, pulled from their length-1 vectors
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    int n_iters_v = Rf_asInteger(n_iters);

    // convert what Fortran cannot take from R directly
    char* mode_c = tox_char_in(mode, 6);

    // outputs and work space
    SEXP fitted_values = PROTECT(Rf_allocVector(REALSXP, ((int) Rf_length(y)))); nprot++;
    int ierr = 0;

    loess_c(
        REAL(x),
        &n_x_elements,
        REAL(y),
        &n_y_elements,
        &span_v,
        &degree_v,
        REAL(fitted_values),
        mode_c,
        &n_iters_v,
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
