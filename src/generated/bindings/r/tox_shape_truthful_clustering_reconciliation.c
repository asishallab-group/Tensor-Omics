// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void ensemble_reconciliation_c(const unsigned char*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const double*, const double*, const int*, const unsigned char*, const int*, const char*, const double*, const unsigned char*, const unsigned char*, const int*, const int*, const double*, const int*, int*, int*, double*, unsigned char*, unsigned char*, unsigned char*, unsigned char*, int*);
void merge_to_super_ensembles_c(const unsigned char*, const unsigned char*, const int*, const int*, const char*, const double*, const unsigned char*, const int*, int*, int*, double*, int*);

SEXP ensemble_reconciliation_call(SEXP ensemble_masks, SEXP ensemble_stop_reason, SEXP ensemble_U_history, SEXP ensemble_d_history, SEXP ensemble_S_history, SEXP ensemble_mu_history, SEXP ensemble_G_history, SEXP ensemble_k_history, SEXP ensemble_accepted_history, SEXP mode, SEXP min_overlap_coefficient, SEXP report_overlap_coefficient, SEXP allowed_stop_reasons, SEXP filter_dim_min, SEXP filter_dim_max, SEXP var_explained_min, SEXP max_group_size) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int allowed_stop_reasons_size = 0;
    if (allowed_stop_reasons != R_NilValue) {
        allowed_stop_reasons_size = (int) Rf_length(allowed_stop_reasons);
    }
    const int* filter_dim_min_p = NULL;
    int filter_dim_min_size = 0;
    if (filter_dim_min != R_NilValue) {
        filter_dim_min_size = (int) Rf_length(filter_dim_min);
        filter_dim_min_p = INTEGER(filter_dim_min);
    }
    const int* filter_dim_max_p = NULL;
    int filter_dim_max_size = 0;
    if (filter_dim_max != R_NilValue) {
        filter_dim_max_size = (int) Rf_length(filter_dim_max);
        filter_dim_max_p = INTEGER(filter_dim_max);
    }
    const double* var_explained_min_p = NULL;
    int var_explained_min_size = 0;
    if (var_explained_min != R_NilValue) {
        var_explained_min_size = (int) Rf_length(var_explained_min);
        var_explained_min_p = REAL(var_explained_min);
    }

    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[0];
    int n_ensembles = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[1];
    int o = INTEGER(Rf_getAttrib(ensemble_U_history, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double min_overlap_coefficient_v = Rf_asReal(min_overlap_coefficient);
    unsigned char report_overlap_coefficient_v = (Rf_asLogical(report_overlap_coefficient) == TRUE) ? 1 : 0;
    int max_group_size_v = Rf_asInteger(max_group_size);

    // convert what Fortran cannot take from R directly
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);
    unsigned char* ensemble_accepted_history_c = tox_bool_in(ensemble_accepted_history);
    char* mode_c = tox_char_in(mode, 25);
    unsigned char* allowed_stop_reasons_c = tox_bool_in(allowed_stop_reasons);

    // outputs and work space
    SEXP super_ensembles = PROTECT(Rf_allocVector(INTSXP, max_group_size_v * (n_ensembles*(n_ensembles-1)))); nprot++;
    { SEXP super_ensembles_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(super_ensembles_dim)[0] = max_group_size_v; INTEGER(super_ensembles_dim)[1] = n_ensembles*(n_ensembles-1); Rf_setAttrib(super_ensembles, R_DimSymbol, super_ensembles_dim); UNPROTECT(1); }
    int n_super_ensembles = 0;
    SEXP super_ensembles_overlap_coefficient = PROTECT(Rf_allocVector(REALSXP, (max_group_size_v-1) * (n_ensembles*(n_ensembles-1)))); nprot++;
    { SEXP super_ensembles_overlap_coefficient_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(super_ensembles_overlap_coefficient_dim)[0] = max_group_size_v-1; INTEGER(super_ensembles_overlap_coefficient_dim)[1] = n_ensembles*(n_ensembles-1); Rf_setAttrib(super_ensembles_overlap_coefficient, R_DimSymbol, super_ensembles_overlap_coefficient_dim); UNPROTECT(1); }
    unsigned char* eligible_c = tox_bool_alloc(n_ensembles);
    unsigned char* eligible_by_stop_condition_c = tox_bool_alloc(n_ensembles);
    unsigned char* eligible_by_dimension_c = tox_bool_alloc(n_ensembles);
    unsigned char* eligible_by_var_explained_c = tox_bool_alloc(n_ensembles);
    int ierr = 0;

    ensemble_reconciliation_c(
        ensemble_masks_c,
        INTEGER(ensemble_stop_reason),
        &n_dimensions,
        &n_vectors,
        &n_ensembles,
        REAL(ensemble_U_history),
        INTEGER(ensemble_d_history),
        REAL(ensemble_S_history),
        REAL(ensemble_mu_history),
        REAL(ensemble_G_history),
        INTEGER(ensemble_k_history),
        ensemble_accepted_history_c,
        &o,
        mode_c,
        &min_overlap_coefficient_v,
        &report_overlap_coefficient_v,
        allowed_stop_reasons != R_NilValue ? allowed_stop_reasons_c : NULL,
        filter_dim_min_p,
        filter_dim_max_p,
        var_explained_min_p,
        &max_group_size_v,
        INTEGER(super_ensembles),
        &n_super_ensembles,
        REAL(super_ensembles_overlap_coefficient),
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

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 8)); nprot++;
    SET_VECTOR_ELT(_out, 0, super_ensembles);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_super_ensembles));
    SET_VECTOR_ELT(_out, 2, super_ensembles_overlap_coefficient);
    SET_VECTOR_ELT(_out, 3, eligible);
    SET_VECTOR_ELT(_out, 4, eligible_by_stop_condition);
    SET_VECTOR_ELT(_out, 5, eligible_by_dimension);
    SET_VECTOR_ELT(_out, 6, eligible_by_var_explained);
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 8)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("super_ensembles"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_super_ensembles"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("super_ensembles_overlap_coefficient"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("eligible"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("eligible_by_stop_condition"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("eligible_by_dimension"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("eligible_by_var_explained"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP merge_to_super_ensembles_call(SEXP ensemble_masks, SEXP eligible, SEXP mode, SEXP min_overlap_coefficient, SEXP report_overlap_coefficient, SEXP max_group_size) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_vectors = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[0];
    int n_ensembles = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double min_overlap_coefficient_v = Rf_asReal(min_overlap_coefficient);
    unsigned char report_overlap_coefficient_v = (Rf_asLogical(report_overlap_coefficient) == TRUE) ? 1 : 0;
    int max_group_size_v = Rf_asInteger(max_group_size);

    // convert what Fortran cannot take from R directly
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);
    unsigned char* eligible_c = tox_bool_in(eligible);
    char* mode_c = tox_char_in(mode, 25);

    // outputs and work space
    SEXP super_ensembles = PROTECT(Rf_allocVector(INTSXP, max_group_size_v * (n_ensembles*(n_ensembles-1)))); nprot++;
    { SEXP super_ensembles_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(super_ensembles_dim)[0] = max_group_size_v; INTEGER(super_ensembles_dim)[1] = n_ensembles*(n_ensembles-1); Rf_setAttrib(super_ensembles, R_DimSymbol, super_ensembles_dim); UNPROTECT(1); }
    int n_super_ensembles = 0;
    SEXP super_ensembles_overlap_coefficient = PROTECT(Rf_allocVector(REALSXP, (max_group_size_v-1) * (n_ensembles*(n_ensembles-1)))); nprot++;
    { SEXP super_ensembles_overlap_coefficient_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(super_ensembles_overlap_coefficient_dim)[0] = max_group_size_v-1; INTEGER(super_ensembles_overlap_coefficient_dim)[1] = n_ensembles*(n_ensembles-1); Rf_setAttrib(super_ensembles_overlap_coefficient, R_DimSymbol, super_ensembles_overlap_coefficient_dim); UNPROTECT(1); }
    int ierr = 0;

    merge_to_super_ensembles_c(
        ensemble_masks_c,
        eligible_c,
        &n_vectors,
        &n_ensembles,
        mode_c,
        &min_overlap_coefficient_v,
        &report_overlap_coefficient_v,
        &max_group_size_v,
        INTEGER(super_ensembles),
        &n_super_ensembles,
        REAL(super_ensembles_overlap_coefficient),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, super_ensembles);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_super_ensembles));
    SET_VECTOR_ELT(_out, 2, super_ensembles_overlap_coefficient);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("super_ensembles"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_super_ensembles"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("super_ensembles_overlap_coefficient"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
