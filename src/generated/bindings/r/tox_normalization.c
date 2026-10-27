// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void normalize_unit_length_c(double*, const int*, int*);
void normalization_pipeline_c(const int*, const int*, const double*, double*, const int*, const int*, const double*, const int*, const unsigned char*, int*);
void normalize_by_std_dev_c(const int*, const int*, const double*, double*, const double*, const int*, int*);
void root_mean_sq_normalization_c(const int*, const int*, const double*, double*, int*);
void quantile_normalization_c(const int*, const int*, const double*, double*, double*, int*);
void log2_transformation_c(const int*, const int*, const double*, double*, int*);
void calc_tiss_avg_c(const int*, const int*, const int*, const double*, const int*, double*, int*);
void calc_fchange_c(const int*, const int*, const int*, const int*, const int*, const double*, double*, int*);

SEXP normalize_unit_length_call(SEXP vector) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = (int) Rf_length(vector);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_vector = PROTECT(Rf_duplicate(vector)); _nprot++;

    // outputs and work space
    int ierr = 0;

    normalize_unit_length_c(
        REAL(_out_vector),
        &n_dims,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_vector);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("vector"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP normalization_pipeline_call(SEXP expr, SEXP reps_per_tissue, SEXP span, SEXP degree, SEXP use_quantile) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];
    int n_tissues = (int) Rf_length(reps_per_tissue);

    // scalar inputs, pulled from their length-1 vectors
    double _v_span = Rf_asReal(span);
    int _v_degree = Rf_asInteger(degree);
    unsigned char _v_use_quantile = (Rf_asLogical(use_quantile) == TRUE) ? 1 : 0;

    // convert what Fortran cannot take from R directly
    const int* _c_reps_per_tissue = tox_int_in(reps_per_tissue, "reps_per_tissue");

    // outputs and work space
    SEXP log_transformed_expr = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); _nprot++;
    { SEXP _dim_log_transformed_expr = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_log_transformed_expr)[0] = n_tissues; INTEGER(_dim_log_transformed_expr)[1] = n_genes; Rf_setAttrib(log_transformed_expr, R_DimSymbol, _dim_log_transformed_expr); UNPROTECT(1); }
    int ierr = 0;

    normalization_pipeline_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(log_transformed_expr),
        _c_reps_per_tissue,
        &n_tissues,
        &_v_span,
        &_v_degree,
        &_v_use_quantile,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, log_transformed_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("log_transformed_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP normalize_by_std_dev_call(SEXP expr, SEXP span, SEXP degree) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    double _v_span = Rf_asReal(span);
    int _v_degree = Rf_asInteger(degree);

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); _nprot++;
    { SEXP _dim_normalized_expr = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_normalized_expr)[0] = n_replicates; INTEGER(_dim_normalized_expr)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, _dim_normalized_expr); UNPROTECT(1); }
    int ierr = 0;

    normalize_by_std_dev_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        &_v_span,
        &_v_degree,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP root_mean_sq_normalization_call(SEXP expr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); _nprot++;
    { SEXP _dim_normalized_expr = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_normalized_expr)[0] = n_replicates; INTEGER(_dim_normalized_expr)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, _dim_normalized_expr); UNPROTECT(1); }
    int ierr = 0;

    root_mean_sq_normalization_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP quantile_normalization_call(SEXP expr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); _nprot++;
    { SEXP _dim_normalized_expr = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_normalized_expr)[0] = n_replicates; INTEGER(_dim_normalized_expr)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, _dim_normalized_expr); UNPROTECT(1); }
    SEXP rank_means = PROTECT(Rf_allocVector(REALSXP, n_genes)); _nprot++;
    int ierr = 0;

    quantile_normalization_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        REAL(rank_means),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, rank_means);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("rank_means"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP log2_transformation_call(SEXP expr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_tissues = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP transformed_expr = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); _nprot++;
    { SEXP _dim_transformed_expr = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_transformed_expr)[0] = n_tissues; INTEGER(_dim_transformed_expr)[1] = n_genes; Rf_setAttrib(transformed_expr, R_DimSymbol, _dim_transformed_expr); UNPROTECT(1); }
    int ierr = 0;

    log2_transformation_c(
        &n_genes,
        &n_tissues,
        REAL(expr),
        REAL(transformed_expr),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, transformed_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("transformed_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_tiss_avg_call(SEXP reps_per_tissue, SEXP expr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_tissues = (int) Rf_length(reps_per_tissue);
    int n_expr_elements_dim_1 = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // convert what Fortran cannot take from R directly
    const int* _c_reps_per_tissue = tox_int_in(reps_per_tissue, "reps_per_tissue");

    // outputs and work space
    SEXP tissue_averages = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); _nprot++;
    { SEXP _dim_tissue_averages = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_tissue_averages)[0] = n_tissues; INTEGER(_dim_tissue_averages)[1] = n_genes; Rf_setAttrib(tissue_averages, R_DimSymbol, _dim_tissue_averages); UNPROTECT(1); }
    int ierr = 0;

    calc_tiss_avg_c(
        &n_genes,
        &n_tissues,
        _c_reps_per_tissue,
        REAL(expr),
        &n_expr_elements_dim_1,
        REAL(tissue_averages),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, tissue_averages);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("tissue_averages"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_fchange_call(SEXP control_tissues, SEXP condition_tissues, SEXP expr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_tissues = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];
    int n_pairs = (int) Rf_length(control_tissues);

    // convert what Fortran cannot take from R directly
    const int* _c_control_tissues = tox_int_in(control_tissues, "control_tissues");
    const int* _c_condition_tissues = tox_int_in(condition_tissues, "condition_tissues");

    // outputs and work space
    SEXP fold_changes = PROTECT(Rf_allocVector(REALSXP, n_pairs * n_genes)); _nprot++;
    { SEXP _dim_fold_changes = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_fold_changes)[0] = n_pairs; INTEGER(_dim_fold_changes)[1] = n_genes; Rf_setAttrib(fold_changes, R_DimSymbol, _dim_fold_changes); UNPROTECT(1); }
    int ierr = 0;

    calc_fchange_c(
        &n_genes,
        &n_tissues,
        &n_pairs,
        _c_control_tissues,
        _c_condition_tissues,
        REAL(expr),
        REAL(fold_changes),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, fold_changes);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("fold_changes"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
