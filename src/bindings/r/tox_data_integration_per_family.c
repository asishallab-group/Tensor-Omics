// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void fjct_compute_jsd_expert_c(const int*, const int*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*, int*, int*, double*, double*, unsigned char*, unsigned char*, double*, double*, int*, int*);
void fjct_compute_jsd_c(const int*, const int*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*, int*, int*, double*, double*, int*);
void fjct_compute_masked_jsd_expert_c(const double*, const double*, const int*, const int*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, const double*, double*, int*, int*, int*, double*, double*, double*, double*, int*, int*);
void fjct_compute_masked_jsd_c(const double*, const double*, const int*, const int*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, const double*, double*, int*, int*, int*, double*, double*, double*, double*, int*);
void fjct_compute_contribution_scores_c(const double*, const int*, const int*, double*, double*, int*);

SEXP fjct_compute_jsd_expert_call(SEXP family_idx, SEXP gene_to_family_S1, SEXP gene_to_family_S2, SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighborhood_genes_S1, SEXP neighborhood_genes_S2, SEXP n_bins, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) Rf_length(gene_to_family_S1);
    int n_genes_S2 = (int) Rf_length(gene_to_family_S2);
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int family_idx_v = Rf_asInteger(family_idx);
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    unsigned char* tmp_neighbor_mask_S1_c = tox_bool_alloc(n_neighbors * n_points);
    unsigned char* tmp_neighbor_mask_S2_c = tox_bool_alloc(n_neighbors * n_points);
    double* tmp_pmf_S1 = (double*) R_alloc(n_points * n_bins_v, sizeof(double));
    double* tmp_pmf_S2 = (double*) R_alloc(n_points * n_bins_v, sizeof(double));
    int* tmp_counts = (int*) R_alloc(n_points * n_bins_v, sizeof(int));
    int ierr = 0;

    fjct_compute_jsd_expert_c(
        &family_idx_v,
        INTEGER(gene_to_family_S1),
        INTEGER(gene_to_family_S2),
        &n_genes_S1,
        &n_genes_S2,
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        INTEGER(neighborhood_genes_S1),
        INTEGER(neighborhood_genes_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &n_bins_v,
        &shared_residual_range_v,
        REAL(js_divergences),
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &total_included_n_reps,
        &global_js_divergence,
        REAL(weights),
        tmp_neighbor_mask_S1_c,
        tmp_neighbor_mask_S2_c,
        tmp_pmf_S1,
        tmp_pmf_S2,
        tmp_counts,
        &ierr
    );

    // convert the outputs back
    SEXP tmp_neighbor_mask_S1 = PROTECT(tox_bool_out(tmp_neighbor_mask_S1_c, n_neighbors * n_points)); nprot++;
    { SEXP tmp_neighbor_mask_S1_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(tmp_neighbor_mask_S1_dim)[0] = n_neighbors; INTEGER(tmp_neighbor_mask_S1_dim)[1] = n_points; Rf_setAttrib(tmp_neighbor_mask_S1, R_DimSymbol, tmp_neighbor_mask_S1_dim); UNPROTECT(1); }
    SEXP tmp_neighbor_mask_S2 = PROTECT(tox_bool_out(tmp_neighbor_mask_S2_c, n_neighbors * n_points)); nprot++;
    { SEXP tmp_neighbor_mask_S2_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(tmp_neighbor_mask_S2_dim)[0] = n_neighbors; INTEGER(tmp_neighbor_mask_S2_dim)[1] = n_points; Rf_setAttrib(tmp_neighbor_mask_S2, R_DimSymbol, tmp_neighbor_mask_S2_dim); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(total_included_n_reps));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 5, weights);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("included_n_reps_S1"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps_S2"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("total_included_n_reps"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP fjct_compute_jsd_call(SEXP family_idx, SEXP gene_to_family_S1, SEXP gene_to_family_S2, SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighborhood_genes_S1, SEXP neighborhood_genes_S2, SEXP n_bins, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) Rf_length(gene_to_family_S1);
    int n_genes_S2 = (int) Rf_length(gene_to_family_S2);
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int family_idx_v = Rf_asInteger(family_idx);
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    int ierr = 0;

    fjct_compute_jsd_c(
        &family_idx_v,
        INTEGER(gene_to_family_S1),
        INTEGER(gene_to_family_S2),
        &n_genes_S1,
        &n_genes_S2,
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        INTEGER(neighborhood_genes_S1),
        INTEGER(neighborhood_genes_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &n_bins_v,
        &shared_residual_range_v,
        REAL(js_divergences),
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &total_included_n_reps,
        &global_js_divergence,
        REAL(weights),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(total_included_n_reps));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 5, weights);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("included_n_reps_S1"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps_S2"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("total_included_n_reps"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP fjct_compute_masked_jsd_expert_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2, SEXP n_bins, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // convert what Fortran cannot take from R directly
    unsigned char* neighbor_mask_S1_c = tox_bool_in(neighbor_mask_S1);
    unsigned char* neighbor_mask_S2_c = tox_bool_in(neighbor_mask_S2);

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP pmf_S1 = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins_v)); nprot++;
    { SEXP pmf_S1_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(pmf_S1_dim)[0] = n_points; INTEGER(pmf_S1_dim)[1] = n_bins_v; Rf_setAttrib(pmf_S1, R_DimSymbol, pmf_S1_dim); UNPROTECT(1); }
    SEXP pmf_S2 = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins_v)); nprot++;
    { SEXP pmf_S2_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(pmf_S2_dim)[0] = n_points; INTEGER(pmf_S2_dim)[1] = n_bins_v; Rf_setAttrib(pmf_S2, R_DimSymbol, pmf_S2_dim); UNPROTECT(1); }
    int* tmp_counts = (int*) R_alloc(n_points * n_bins_v, sizeof(int));
    int ierr = 0;

    fjct_compute_masked_jsd_expert_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        neighbor_mask_S1_c,
        neighbor_mask_S2_c,
        &n_bins_v,
        &shared_residual_range_v,
        REAL(js_divergences),
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &total_included_n_reps,
        &global_js_divergence,
        REAL(weights),
        REAL(pmf_S1),
        REAL(pmf_S2),
        tmp_counts,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(total_included_n_reps));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 5, weights);
    SET_VECTOR_ELT(_out, 6, pmf_S1);
    SET_VECTOR_ELT(_out, 7, pmf_S2);
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("included_n_reps_S1"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps_S2"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("total_included_n_reps"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("pmf_S1"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("pmf_S2"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP fjct_compute_masked_jsd_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2, SEXP n_bins, SEXP shared_residual_range) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);

    // convert what Fortran cannot take from R directly
    unsigned char* neighbor_mask_S1_c = tox_bool_in(neighbor_mask_S1);
    unsigned char* neighbor_mask_S2_c = tox_bool_in(neighbor_mask_S2);

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    SEXP pmf_S1 = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins_v)); nprot++;
    { SEXP pmf_S1_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(pmf_S1_dim)[0] = n_points; INTEGER(pmf_S1_dim)[1] = n_bins_v; Rf_setAttrib(pmf_S1, R_DimSymbol, pmf_S1_dim); UNPROTECT(1); }
    SEXP pmf_S2 = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins_v)); nprot++;
    { SEXP pmf_S2_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(pmf_S2_dim)[0] = n_points; INTEGER(pmf_S2_dim)[1] = n_bins_v; Rf_setAttrib(pmf_S2, R_DimSymbol, pmf_S2_dim); UNPROTECT(1); }
    int ierr = 0;

    fjct_compute_masked_jsd_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        neighbor_mask_S1_c,
        neighbor_mask_S2_c,
        &n_bins_v,
        &shared_residual_range_v,
        REAL(js_divergences),
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &total_included_n_reps,
        &global_js_divergence,
        REAL(weights),
        REAL(pmf_S1),
        REAL(pmf_S2),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(total_included_n_reps));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 5, weights);
    SET_VECTOR_ELT(_out, 6, pmf_S1);
    SET_VECTOR_ELT(_out, 7, pmf_S2);
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("included_n_reps_S1"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps_S2"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("total_included_n_reps"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("pmf_S1"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("pmf_S2"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP fjct_compute_contribution_scores_call(SEXP global_js_divergences, SEXP total_included_n_reps_per_f) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int k_families = (int) Rf_length(global_js_divergences);

    // outputs and work space
    SEXP support_weights = PROTECT(Rf_allocVector(REALSXP, k_families)); nprot++;
    SEXP contribution_scores = PROTECT(Rf_allocVector(REALSXP, k_families)); nprot++;
    int ierr = 0;

    fjct_compute_contribution_scores_c(
        REAL(global_js_divergences),
        INTEGER(total_included_n_reps_per_f),
        &k_families,
        REAL(support_weights),
        REAL(contribution_scores),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, support_weights);
    SET_VECTOR_ELT(_out, 1, contribution_scores);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("support_weights"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("contribution_scores"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
