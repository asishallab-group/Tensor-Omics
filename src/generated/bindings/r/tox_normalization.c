// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void normalize_unit_length_c(double*, const int*, int*);
void normalization_pipeline_expert_c(const int*, const int*, const double*, double*, const int*, const int*, double*, double*, int*, double*, int*, const int*, double*, const int*, double*, double*, double*, double*, double*, double*, int*, const double*, const int*, const unsigned char*, int*);
void normalization_pipeline_c(const int*, const int*, const double*, double*, const int*, const int*, const double*, const int*, const unsigned char*, int*);
void normalize_by_std_dev_expert_c(const int*, const int*, const double*, double*, double*, double*, int*, double*, int*, const int*, double*, const int*, double*, double*, double*, double*, double*, double*, int*, const double*, const int*, int*);
void normalize_by_std_dev_c(const int*, const int*, const double*, double*, const double*, const int*, int*);
void root_mean_sq_normalization_c(const int*, const int*, const double*, double*, int*);
void quantile_normalization_expert_c(const int*, const int*, const double*, double*, double*, double*, int*, int*);
void quantile_normalization_c(const int*, const int*, const double*, double*, double*, int*);
void log2_transformation_c(const int*, const int*, const double*, double*, int*);
void calc_tiss_avg_c(const int*, const int*, const int*, const double*, double*, int*);
void calc_fchange_c(const int*, const int*, const int*, const int*, const int*, const double*, double*, int*);

SEXP normalize_unit_length_call(SEXP vector) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = (int) Rf_length(vector);

    // copy what is modified in place, so the caller's stays intact
    SEXP vector_out = PROTECT(Rf_duplicate(vector)); nprot++;

    // outputs and work space
    int ierr = 0;

    normalize_unit_length_c(
        REAL(vector_out),
        &n_dims,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, vector_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("vector"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP normalization_pipeline_expert_call(SEXP expr, SEXP reps_per_tissue, SEXP int_workspace_size, SEXP real_workspace_size, SEXP span, SEXP degree, SEXP use_quantile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];
    int n_tissues = (int) Rf_length(reps_per_tissue);

    // scalar inputs, pulled from their length-1 vectors
    int int_workspace_size_v = Rf_asInteger(int_workspace_size);
    int real_workspace_size_v = Rf_asInteger(real_workspace_size);
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    unsigned char use_quantile_v = (Rf_asLogical(use_quantile) == TRUE) ? 1 : 0;

    // outputs and work space
    SEXP log_transformed_expr = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); nprot++;
    { SEXP log_transformed_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(log_transformed_expr_dim)[0] = n_tissues; INTEGER(log_transformed_expr_dim)[1] = n_genes; Rf_setAttrib(log_transformed_expr, R_DimSymbol, log_transformed_expr_dim); UNPROTECT(1); }
    double* tmp_expr_copy = (double*) R_alloc(n_replicates * n_genes, sizeof(double));
    double* tmp_loess_y = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_indices_used = (int*) R_alloc(n_genes, sizeof(int));
    double* tmp_yhat_global = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_int_workspace = (int*) R_alloc(int_workspace_size_v, sizeof(int));
    double* tmp_real_workspace = (double*) R_alloc(real_workspace_size_v, sizeof(double));
    double* tmp_hat_diag = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_loess_weights = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_eval_points = (double*) R_alloc(n_genes * 1, sizeof(double));
    double* tmp_robust_weights = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_combined_weights = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_residuals = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_permutation_indices = (int*) R_alloc(n_genes, sizeof(int));
    int ierr = 0;

    normalization_pipeline_expert_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(log_transformed_expr),
        INTEGER(reps_per_tissue),
        &n_tissues,
        tmp_expr_copy,
        tmp_loess_y,
        tmp_indices_used,
        tmp_yhat_global,
        tmp_int_workspace,
        &int_workspace_size_v,
        tmp_real_workspace,
        &real_workspace_size_v,
        tmp_hat_diag,
        tmp_loess_weights,
        tmp_eval_points,
        tmp_robust_weights,
        tmp_combined_weights,
        tmp_residuals,
        tmp_permutation_indices,
        &span_v,
        &degree_v,
        &use_quantile_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, log_transformed_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("log_transformed_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP normalization_pipeline_call(SEXP expr, SEXP reps_per_tissue, SEXP span, SEXP degree, SEXP use_quantile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];
    int n_tissues = (int) Rf_length(reps_per_tissue);

    // scalar inputs, pulled from their length-1 vectors
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    unsigned char use_quantile_v = (Rf_asLogical(use_quantile) == TRUE) ? 1 : 0;

    // outputs and work space
    SEXP log_transformed_expr = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); nprot++;
    { SEXP log_transformed_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(log_transformed_expr_dim)[0] = n_tissues; INTEGER(log_transformed_expr_dim)[1] = n_genes; Rf_setAttrib(log_transformed_expr, R_DimSymbol, log_transformed_expr_dim); UNPROTECT(1); }
    int ierr = 0;

    normalization_pipeline_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(log_transformed_expr),
        INTEGER(reps_per_tissue),
        &n_tissues,
        &span_v,
        &degree_v,
        &use_quantile_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, log_transformed_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("log_transformed_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP normalize_by_std_dev_expert_call(SEXP expr, SEXP int_workspace_size, SEXP real_workspace_size, SEXP span, SEXP degree) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int int_workspace_size_v = Rf_asInteger(int_workspace_size);
    int real_workspace_size_v = Rf_asInteger(real_workspace_size);
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); nprot++;
    { SEXP normalized_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(normalized_expr_dim)[0] = n_replicates; INTEGER(normalized_expr_dim)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, normalized_expr_dim); UNPROTECT(1); }
    double* tmp_loess_x = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_loess_y = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_indices_used = (int*) R_alloc(n_genes, sizeof(int));
    double* tmp_yhat_global = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_int_workspace = (int*) R_alloc(int_workspace_size_v, sizeof(int));
    double* tmp_real_workspace = (double*) R_alloc(real_workspace_size_v, sizeof(double));
    double* tmp_hat_diag = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_loess_weights = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_eval_points = (double*) R_alloc(n_genes * 1, sizeof(double));
    double* tmp_robust_weights = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_combined_weights = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_residuals = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_permutation_indices = (int*) R_alloc(n_genes, sizeof(int));
    int ierr = 0;

    normalize_by_std_dev_expert_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        tmp_loess_x,
        tmp_loess_y,
        tmp_indices_used,
        tmp_yhat_global,
        tmp_int_workspace,
        &int_workspace_size_v,
        tmp_real_workspace,
        &real_workspace_size_v,
        tmp_hat_diag,
        tmp_loess_weights,
        tmp_eval_points,
        tmp_robust_weights,
        tmp_combined_weights,
        tmp_residuals,
        tmp_permutation_indices,
        &span_v,
        &degree_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP normalize_by_std_dev_call(SEXP expr, SEXP span, SEXP degree) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); nprot++;
    { SEXP normalized_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(normalized_expr_dim)[0] = n_replicates; INTEGER(normalized_expr_dim)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, normalized_expr_dim); UNPROTECT(1); }
    int ierr = 0;

    normalize_by_std_dev_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        &span_v,
        &degree_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP root_mean_sq_normalization_call(SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); nprot++;
    { SEXP normalized_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(normalized_expr_dim)[0] = n_replicates; INTEGER(normalized_expr_dim)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, normalized_expr_dim); UNPROTECT(1); }
    int ierr = 0;

    root_mean_sq_normalization_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP quantile_normalization_expert_call(SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); nprot++;
    { SEXP normalized_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(normalized_expr_dim)[0] = n_replicates; INTEGER(normalized_expr_dim)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, normalized_expr_dim); UNPROTECT(1); }
    SEXP rank_means = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    double* tmp_genes_row = (double*) R_alloc(n_genes, sizeof(double));
    int* tmp_perm = (int*) R_alloc(n_genes, sizeof(int));
    int ierr = 0;

    quantile_normalization_expert_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        REAL(rank_means),
        tmp_genes_row,
        tmp_perm,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, rank_means);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("rank_means"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP quantile_normalization_call(SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_replicates = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP normalized_expr = PROTECT(Rf_allocVector(REALSXP, n_replicates * n_genes)); nprot++;
    { SEXP normalized_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(normalized_expr_dim)[0] = n_replicates; INTEGER(normalized_expr_dim)[1] = n_genes; Rf_setAttrib(normalized_expr, R_DimSymbol, normalized_expr_dim); UNPROTECT(1); }
    SEXP rank_means = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    quantile_normalization_c(
        &n_genes,
        &n_replicates,
        REAL(expr),
        REAL(normalized_expr),
        REAL(rank_means),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, normalized_expr);
    SET_VECTOR_ELT(_out, 1, rank_means);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normalized_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("rank_means"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP log2_transformation_call(SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_tissues = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP transformed_expr = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); nprot++;
    { SEXP transformed_expr_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(transformed_expr_dim)[0] = n_tissues; INTEGER(transformed_expr_dim)[1] = n_genes; Rf_setAttrib(transformed_expr, R_DimSymbol, transformed_expr_dim); UNPROTECT(1); }
    int ierr = 0;

    log2_transformation_c(
        &n_genes,
        &n_tissues,
        REAL(expr),
        REAL(transformed_expr),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, transformed_expr);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("transformed_expr"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP calc_tiss_avg_call(SEXP reps_per_tissue, SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_tissues = (int) Rf_length(reps_per_tissue);

    // outputs and work space
    SEXP tissue_averages = PROTECT(Rf_allocVector(REALSXP, n_tissues * n_genes)); nprot++;
    { SEXP tissue_averages_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(tissue_averages_dim)[0] = n_tissues; INTEGER(tissue_averages_dim)[1] = n_genes; Rf_setAttrib(tissue_averages, R_DimSymbol, tissue_averages_dim); UNPROTECT(1); }
    int ierr = 0;

    calc_tiss_avg_c(
        &n_genes,
        &n_tissues,
        INTEGER(reps_per_tissue),
        REAL(expr),
        REAL(tissue_averages),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, tissue_averages);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("tissue_averages"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP calc_fchange_call(SEXP control_tissues, SEXP condition_tissues, SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_tissues = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];
    int n_pairs = (int) Rf_length(control_tissues);

    // outputs and work space
    SEXP fold_changes = PROTECT(Rf_allocVector(REALSXP, n_pairs * n_genes)); nprot++;
    { SEXP fold_changes_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(fold_changes_dim)[0] = n_pairs; INTEGER(fold_changes_dim)[1] = n_genes; Rf_setAttrib(fold_changes, R_DimSymbol, fold_changes_dim); UNPROTECT(1); }
    int ierr = 0;

    calc_fchange_c(
        &n_genes,
        &n_tissues,
        &n_pairs,
        INTEGER(control_tissues),
        INTEGER(condition_tissues),
        REAL(expr),
        REAL(fold_changes),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, fold_changes);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("fold_changes"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
