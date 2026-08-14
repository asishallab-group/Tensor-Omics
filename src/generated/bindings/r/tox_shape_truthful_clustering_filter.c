// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void filter_ensembles_by_stop_condition_c(const int*, const int*, const unsigned char*, unsigned char*, int*);
void filter_ensembles_by_dimension_c(const int*, const int*, const int*, const unsigned char*, const int*, const int*, unsigned char*, int*);
void filter_ensembles_by_var_explained_c(const int*, const int*, const double*, const int*, const int*, const unsigned char*, const double*, unsigned char*, int*);
void filter_ensembles_c(const int*, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const int*, const unsigned char*, const int*, const unsigned char*, const int*, const int*, const double*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, int*);

SEXP filter_ensembles_by_stop_condition_call(SEXP ensemble_stop_reason, SEXP allowed_stop_reasons) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int allowed_stop_reasons_size = 0;
    if (allowed_stop_reasons != R_NilValue) {
        allowed_stop_reasons_size = (int) Rf_length(allowed_stop_reasons);
    }

    // derived from the inputs, not asked of the caller
    int n_ensembles = (int) Rf_length(ensemble_stop_reason);

    // convert what Fortran cannot take from R directly
    unsigned char* allowed_stop_reasons_c = tox_bool_in(allowed_stop_reasons);

    // outputs and work space
    unsigned char* eligible_c = tox_bool_alloc(n_ensembles);
    int ierr = 0;

    filter_ensembles_by_stop_condition_c(
        &n_ensembles,
        INTEGER(ensemble_stop_reason),
        allowed_stop_reasons != R_NilValue ? allowed_stop_reasons_c : NULL,
        eligible_c,
        &ierr
    );

    // convert the outputs back
    SEXP eligible = PROTECT(tox_bool_out(eligible_c, n_ensembles)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, eligible);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("eligible"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP filter_ensembles_by_dimension_call(SEXP n_dimensions, SEXP ensemble_d_final, SEXP ensemble_has_final, SEXP d_min, SEXP d_max) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* d_min_p = NULL;
    int d_min_size = 0;
    if (d_min != R_NilValue) {
        d_min_size = (int) Rf_length(d_min);
        d_min_p = INTEGER(d_min);
    }
    const int* d_max_p = NULL;
    int d_max_size = 0;
    if (d_max != R_NilValue) {
        d_max_size = (int) Rf_length(d_max);
        d_max_p = INTEGER(d_max);
    }

    // derived from the inputs, not asked of the caller
    int n_ensembles = (int) Rf_length(ensemble_d_final);

    // scalar inputs, pulled from their length-1 vectors
    int n_dimensions_v = Rf_asInteger(n_dimensions);

    // convert what Fortran cannot take from R directly
    unsigned char* ensemble_has_final_c = tox_bool_in(ensemble_has_final);

    // outputs and work space
    unsigned char* eligible_c = tox_bool_alloc(n_ensembles);
    int ierr = 0;

    filter_ensembles_by_dimension_c(
        &n_dimensions_v,
        &n_ensembles,
        INTEGER(ensemble_d_final),
        ensemble_has_final_c,
        d_min_p,
        d_max_p,
        eligible_c,
        &ierr
    );

    // convert the outputs back
    SEXP eligible = PROTECT(tox_bool_out(eligible_c, n_ensembles)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, eligible);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("eligible"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP filter_ensembles_by_var_explained_call(SEXP ensemble_S_final, SEXP ensemble_d_final, SEXP ensemble_k_final, SEXP ensemble_has_final, SEXP var_explained_min) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const double* var_explained_min_p = NULL;
    int var_explained_min_size = 0;
    if (var_explained_min != R_NilValue) {
        var_explained_min_size = (int) Rf_length(var_explained_min);
        var_explained_min_p = REAL(var_explained_min);
    }

    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(ensemble_S_final, R_DimSymbol))[0];
    int n_ensembles = INTEGER(Rf_getAttrib(ensemble_S_final, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    unsigned char* ensemble_has_final_c = tox_bool_in(ensemble_has_final);

    // outputs and work space
    unsigned char* eligible_c = tox_bool_alloc(n_ensembles);
    int ierr = 0;

    filter_ensembles_by_var_explained_c(
        &n_dimensions,
        &n_ensembles,
        REAL(ensemble_S_final),
        INTEGER(ensemble_d_final),
        INTEGER(ensemble_k_final),
        ensemble_has_final_c,
        var_explained_min_p,
        eligible_c,
        &ierr
    );

    // convert the outputs back
    SEXP eligible = PROTECT(tox_bool_out(eligible_c, n_ensembles)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, eligible);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("eligible"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP filter_ensembles_call(SEXP ensemble_U_history, SEXP ensemble_d_history, SEXP ensemble_S_history, SEXP ensemble_mu_history, SEXP ensemble_G_history, SEXP ensemble_k_history, SEXP ensemble_accepted_history, SEXP ensemble_stop_reason, SEXP allowed_stop_reasons, SEXP d_min, SEXP d_max, SEXP var_explained_min) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int allowed_stop_reasons_size = 0;
    if (allowed_stop_reasons != R_NilValue) {
        allowed_stop_reasons_size = (int) Rf_length(allowed_stop_reasons);
    }
    const int* d_min_p = NULL;
    int d_min_size = 0;
    if (d_min != R_NilValue) {
        d_min_size = (int) Rf_length(d_min);
        d_min_p = INTEGER(d_min);
    }
    const int* d_max_p = NULL;
    int d_max_size = 0;
    if (d_max != R_NilValue) {
        d_max_size = (int) Rf_length(d_max);
        d_max_p = INTEGER(d_max);
    }
    const double* var_explained_min_p = NULL;
    int var_explained_min_size = 0;
    if (var_explained_min != R_NilValue) {
        var_explained_min_size = (int) Rf_length(var_explained_min);
        var_explained_min_p = REAL(var_explained_min);
    }

    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[0];
    int o = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[2];
    int n_ensembles = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[3];

    // convert what Fortran cannot take from R directly
    unsigned char* ensemble_accepted_history_c = tox_bool_in(ensemble_accepted_history);
    unsigned char* allowed_stop_reasons_c = tox_bool_in(allowed_stop_reasons);

    // outputs and work space
    unsigned char* eligible_c = tox_bool_alloc(n_ensembles);
    unsigned char* eligible_by_stop_condition_c = tox_bool_alloc(n_ensembles);
    unsigned char* eligible_by_dimension_c = tox_bool_alloc(n_ensembles);
    unsigned char* eligible_by_var_explained_c = tox_bool_alloc(n_ensembles);
    int ierr = 0;

    filter_ensembles_c(
        &n_dimensions,
        &o,
        &n_ensembles,
        REAL(ensemble_U_history),
        INTEGER(ensemble_d_history),
        REAL(ensemble_S_history),
        REAL(ensemble_mu_history),
        REAL(ensemble_G_history),
        INTEGER(ensemble_k_history),
        ensemble_accepted_history_c,
        INTEGER(ensemble_stop_reason),
        allowed_stop_reasons != R_NilValue ? allowed_stop_reasons_c : NULL,
        d_min_p,
        d_max_p,
        var_explained_min_p,
        eligible_c,
        eligible_by_stop_condition_c,
        eligible_by_dimension_c,
        eligible_by_var_explained_c,
        &ierr
    );

    // convert the outputs back
    SEXP eligible = PROTECT(tox_bool_out(eligible_c, n_ensembles)); nprot++;
    SEXP eligible_by_stop_condition = PROTECT(tox_bool_out(eligible_by_stop_condition_c, n_ensembles)); nprot++;
    SEXP eligible_by_dimension = PROTECT(tox_bool_out(eligible_by_dimension_c, n_ensembles)); nprot++;
    SEXP eligible_by_var_explained = PROTECT(tox_bool_out(eligible_by_var_explained_c, n_ensembles)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, eligible);
    SET_VECTOR_ELT(_out, 1, eligible_by_stop_condition);
    SET_VECTOR_ELT(_out, 2, eligible_by_dimension);
    SET_VECTOR_ELT(_out, 3, eligible_by_var_explained);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("eligible"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("eligible_by_stop_condition"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("eligible_by_dimension"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("eligible_by_var_explained"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
