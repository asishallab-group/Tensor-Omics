// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void compute_family_scaling_expert_c(const int*, const int*, const double*, const int*, double*, double*, double*, int*, int*, int*, int*, int*, const int*, double*, const int*, double*, double*, double*, double*, double*, double*, int*, double*, const double*, const int*, const char*, const int*, double*, int*, double*, int*);
void compute_family_scaling_c(const int*, const int*, const double*, const int*, double*, double*, double*, int*, const double*, const int*, const char*, const int*, double*, int*, int*);
void compute_rdi_expert_c(const int*, const double*, const int*, const double*, const int*, double*, double*, int*, int*, int*, int*);
void compute_rdi_c(const int*, const double*, const int*, const double*, const int*, double*, double*, int*, int*);
void identify_outliers_c(const int*, const double*, const double*, const int*, unsigned char*, double*, double*, const double*, int*);
void detect_outliers_expert_c(const int*, const int*, const double*, const int*, int*, int*, int*, int*, const int*, double*, const int*, double*, double*, double*, double*, double*, double*, int*, double*, double*, double*, int*, double*, double*, double*, double*, unsigned char*, double*, double*, int*, double*, int*, const double*);
void detect_outliers_c(const int*, const int*, const double*, const int*, unsigned char*, double*, double*, int*, double*, int*, const double*);

SEXP compute_family_scaling_expert_call(SEXP n_families, SEXP distances, SEXP gene_to_fam, SEXP int_workspace_size, SEXP real_workspace_size, SEXP span, SEXP degree, SEXP mode, SEXP n_iters) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    int int_workspace_size_v = Rf_asInteger(int_workspace_size);
    int real_workspace_size_v = Rf_asInteger(real_workspace_size);
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    int n_iters_v = Rf_asInteger(n_iters);

    // convert what Fortran cannot take from R directly
    char* mode_c = tox_char_in(mode, 6);

    // outputs and work space
    SEXP dscale = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_x = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_y = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP indices_used = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    int* tmp_perm = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_stack_left = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_stack_right = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_int_workspace = (int*) R_alloc(int_workspace_size_v, sizeof(int));
    double* tmp_real_workspace = (double*) R_alloc(real_workspace_size_v, sizeof(double));
    double* tmp_diagl = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_weights = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_eval_points = (double*) R_alloc(n_families_v * 1, sizeof(double));
    double* tmp_robust_weights = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_combined_weights = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_residuals = (double*) R_alloc(n_families_v, sizeof(double));
    int* tmp_permutation_indices = (int*) R_alloc(n_families_v, sizeof(int));
    double* tmp_fitted_values = (double*) R_alloc(n_families_v, sizeof(double));
    double low_sd_cutoff = 0;
    SEXP excluded_low_sd = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    double* tmp_means_aux = (double*) R_alloc(n_families_v, sizeof(double));
    int ierr = 0;

    compute_family_scaling_expert_c(
        &n_genes,
        &n_families_v,
        REAL(distances),
        INTEGER(gene_to_fam),
        REAL(dscale),
        REAL(loess_x),
        REAL(loess_y),
        INTEGER(indices_used),
        tmp_perm,
        tmp_stack_left,
        tmp_stack_right,
        tmp_int_workspace,
        &int_workspace_size_v,
        tmp_real_workspace,
        &real_workspace_size_v,
        tmp_diagl,
        tmp_weights,
        tmp_eval_points,
        tmp_robust_weights,
        tmp_combined_weights,
        tmp_residuals,
        tmp_permutation_indices,
        tmp_fitted_values,
        &span_v,
        &degree_v,
        mode_c,
        &n_iters_v,
        &low_sd_cutoff,
        INTEGER(excluded_low_sd),
        tmp_means_aux,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); nprot++;
    SET_VECTOR_ELT(_out, 0, dscale);
    SET_VECTOR_ELT(_out, 1, loess_x);
    SET_VECTOR_ELT(_out, 2, loess_y);
    SET_VECTOR_ELT(_out, 3, indices_used);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(low_sd_cutoff));
    SET_VECTOR_ELT(_out, 5, excluded_low_sd);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("dscale"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("loess_x"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("loess_y"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("indices_used"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("low_sd_cutoff"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("excluded_low_sd"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_family_scaling_call(SEXP n_families, SEXP distances, SEXP gene_to_fam, SEXP span, SEXP degree, SEXP mode, SEXP n_iters) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    double span_v = Rf_asReal(span);
    int degree_v = Rf_asInteger(degree);
    int n_iters_v = Rf_asInteger(n_iters);

    // convert what Fortran cannot take from R directly
    char* mode_c = tox_char_in(mode, 6);

    // outputs and work space
    SEXP dscale = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_x = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_y = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP indices_used = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    double low_sd_cutoff = 0;
    SEXP excluded_low_sd = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    int ierr = 0;

    compute_family_scaling_c(
        &n_genes,
        &n_families_v,
        REAL(distances),
        INTEGER(gene_to_fam),
        REAL(dscale),
        REAL(loess_x),
        REAL(loess_y),
        INTEGER(indices_used),
        &span_v,
        &degree_v,
        mode_c,
        &n_iters_v,
        &low_sd_cutoff,
        INTEGER(excluded_low_sd),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); nprot++;
    SET_VECTOR_ELT(_out, 0, dscale);
    SET_VECTOR_ELT(_out, 1, loess_x);
    SET_VECTOR_ELT(_out, 2, loess_y);
    SET_VECTOR_ELT(_out, 3, indices_used);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(low_sd_cutoff));
    SET_VECTOR_ELT(_out, 5, excluded_low_sd);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("dscale"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("loess_x"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("loess_y"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("indices_used"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("low_sd_cutoff"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("excluded_low_sd"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_rdi_expert_call(SEXP distances, SEXP gene_to_fam, SEXP dscale) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);
    int n_dscale_elements = (int) Rf_length(dscale);

    // outputs and work space
    SEXP rdi = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    SEXP sorted_rdi = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    SEXP perm = PROTECT(Rf_allocVector(INTSXP, n_genes)); nprot++;
    int* tmp_stack_left = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_stack_right = (int*) R_alloc(n_genes, sizeof(int));
    int ierr = 0;

    compute_rdi_expert_c(
        &n_genes,
        REAL(distances),
        INTEGER(gene_to_fam),
        REAL(dscale),
        &n_dscale_elements,
        REAL(rdi),
        REAL(sorted_rdi),
        INTEGER(perm),
        tmp_stack_left,
        tmp_stack_right,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, rdi);
    SET_VECTOR_ELT(_out, 1, sorted_rdi);
    SET_VECTOR_ELT(_out, 2, perm);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("rdi"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sorted_rdi"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("perm"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_rdi_call(SEXP distances, SEXP gene_to_fam, SEXP dscale) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);
    int n_dscale_elements = (int) Rf_length(dscale);

    // outputs and work space
    SEXP rdi = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    SEXP sorted_rdi = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    SEXP perm = PROTECT(Rf_allocVector(INTSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_rdi_c(
        &n_genes,
        REAL(distances),
        INTEGER(gene_to_fam),
        REAL(dscale),
        &n_dscale_elements,
        REAL(rdi),
        REAL(sorted_rdi),
        INTEGER(perm),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, rdi);
    SET_VECTOR_ELT(_out, 1, sorted_rdi);
    SET_VECTOR_ELT(_out, 2, perm);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("rdi"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sorted_rdi"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("perm"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP identify_outliers_call(SEXP rdi, SEXP sorted_rdi, SEXP perm, SEXP percentile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(rdi);

    // scalar inputs, pulled from their length-1 vectors
    double percentile_v = Rf_asReal(percentile);

    // outputs and work space
    unsigned char* is_outlier_c = tox_bool_alloc(n_genes);
    double threshold = 0;
    SEXP quantile = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    identify_outliers_c(
        &n_genes,
        REAL(rdi),
        REAL(sorted_rdi),
        INTEGER(perm),
        is_outlier_c,
        &threshold,
        REAL(quantile),
        &percentile_v,
        &ierr
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(is_outlier_c, n_genes)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, is_outlier);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(threshold));
    SET_VECTOR_ELT(_out, 2, quantile);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("threshold"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("quantile"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_outliers_expert_call(SEXP n_families, SEXP distances, SEXP gene_to_fam, SEXP int_workspace_size, SEXP real_workspace_size, SEXP percentile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    int int_workspace_size_v = Rf_asInteger(int_workspace_size);
    int real_workspace_size_v = Rf_asInteger(real_workspace_size);
    double percentile_v = Rf_asReal(percentile);

    // outputs and work space
    int* tmp_perm = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_stack_left = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_stack_right = (int*) R_alloc(n_genes, sizeof(int));
    int* tmp_int_workspace = (int*) R_alloc(int_workspace_size_v, sizeof(int));
    double* tmp_real_workspace = (double*) R_alloc(real_workspace_size_v, sizeof(double));
    double* tmp_diagl = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_weights = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_eval_points = (double*) R_alloc(n_families_v * 1, sizeof(double));
    double* tmp_robust_weights = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_combined_weights = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_residuals = (double*) R_alloc(n_families_v, sizeof(double));
    int* tmp_permutation_indices = (int*) R_alloc(n_families_v, sizeof(int));
    double* tmp_fitted_values = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_means_aux = (double*) R_alloc(n_families_v, sizeof(double));
    double* tmp_dscale = (double*) R_alloc(n_families_v, sizeof(double));
    int* tmp_excluded_low_sd = (int*) R_alloc(n_families_v, sizeof(int));
    double tmp_low_sd_cutoff = 0;
    double* tmp_rdi = (double*) R_alloc(n_genes, sizeof(double));
    double* tmp_sorted_rdi = (double*) R_alloc(n_genes, sizeof(double));
    double tmp_threshold = 0;
    unsigned char* is_outlier_c = tox_bool_alloc(n_genes);
    SEXP loess_x = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_y = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_n = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP quantile = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    detect_outliers_expert_c(
        &n_genes,
        &n_families_v,
        REAL(distances),
        INTEGER(gene_to_fam),
        tmp_perm,
        tmp_stack_left,
        tmp_stack_right,
        tmp_int_workspace,
        &int_workspace_size_v,
        tmp_real_workspace,
        &real_workspace_size_v,
        tmp_diagl,
        tmp_weights,
        tmp_eval_points,
        tmp_robust_weights,
        tmp_combined_weights,
        tmp_residuals,
        tmp_permutation_indices,
        tmp_fitted_values,
        tmp_means_aux,
        tmp_dscale,
        tmp_excluded_low_sd,
        &tmp_low_sd_cutoff,
        tmp_rdi,
        tmp_sorted_rdi,
        &tmp_threshold,
        is_outlier_c,
        REAL(loess_x),
        REAL(loess_y),
        INTEGER(loess_n),
        REAL(quantile),
        &ierr,
        &percentile_v
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(is_outlier_c, n_genes)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
    SET_VECTOR_ELT(_out, 0, is_outlier);
    SET_VECTOR_ELT(_out, 1, loess_x);
    SET_VECTOR_ELT(_out, 2, loess_y);
    SET_VECTOR_ELT(_out, 3, loess_n);
    SET_VECTOR_ELT(_out, 4, quantile);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("loess_x"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("loess_y"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("loess_n"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("quantile"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_outliers_call(SEXP n_families, SEXP distances, SEXP gene_to_fam, SEXP percentile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    double percentile_v = Rf_asReal(percentile);

    // outputs and work space
    unsigned char* is_outlier_c = tox_bool_alloc(n_genes);
    SEXP loess_x = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_y = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP loess_n = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP quantile = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    detect_outliers_c(
        &n_genes,
        &n_families_v,
        REAL(distances),
        INTEGER(gene_to_fam),
        is_outlier_c,
        REAL(loess_x),
        REAL(loess_y),
        INTEGER(loess_n),
        REAL(quantile),
        &ierr,
        &percentile_v
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(is_outlier_c, n_genes)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
    SET_VECTOR_ELT(_out, 0, is_outlier);
    SET_VECTOR_ELT(_out, 1, loess_x);
    SET_VECTOR_ELT(_out, 2, loess_y);
    SET_VECTOR_ELT(_out, 3, loess_n);
    SET_VECTOR_ELT(_out, 4, quantile);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("loess_x"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("loess_y"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("loess_n"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("quantile"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
