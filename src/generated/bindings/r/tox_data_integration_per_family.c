// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void fjct_compute_jsd_c(const int*, const int*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*, int*, int*, double*, double*, int*);
void fjct_compute_masked_jsd_c(const double*, const double*, const int*, const int*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, const double*, double*, int*, int*, int*, double*, double*, double*, double*, int*);
void fjct_compute_contribution_scores_c(const double*, const int*, const int*, double*, double*, int*);

SEXP fjct_compute_jsd_call(SEXP family_idx, SEXP gene_to_family_S1, SEXP gene_to_family_S2, SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighborhood_genes_S1, SEXP neighborhood_genes_S2, SEXP n_bins, SEXP shared_residual_range) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) Rf_length(gene_to_family_S1);
    int n_genes_S2 = (int) Rf_length(gene_to_family_S2);
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int _v_family_idx = Rf_asInteger(family_idx);
    int _v_n_bins = Rf_asInteger(n_bins);
    double _v_shared_residual_range = Rf_asReal(shared_residual_range);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_family_S1 = tox_int_in(gene_to_family_S1, "gene_to_family_S1");
    const int* _c_gene_to_family_S2 = tox_int_in(gene_to_family_S2, "gene_to_family_S2");
    const int* _c_neighborhood_genes_S1 = tox_int_in(neighborhood_genes_S1, "neighborhood_genes_S1");
    const int* _c_neighborhood_genes_S2 = tox_int_in(neighborhood_genes_S2, "neighborhood_genes_S2");

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    int ierr = 0;

    fjct_compute_jsd_c(
        &_v_family_idx,
        _c_gene_to_family_S1,
        _c_gene_to_family_S2,
        &n_genes_S1,
        &n_genes_S2,
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        _c_neighborhood_genes_S1,
        _c_neighborhood_genes_S2,
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &_v_n_bins,
        &_v_shared_residual_range,
        REAL(js_divergences),
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &total_included_n_reps,
        &global_js_divergence,
        REAL(weights),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); _nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(total_included_n_reps));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 5, weights);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("included_n_reps_S1"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps_S2"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("total_included_n_reps"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP fjct_compute_masked_jsd_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2, SEXP n_bins, SEXP shared_residual_range) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bins = Rf_asInteger(n_bins);
    double _v_shared_residual_range = Rf_asReal(shared_residual_range);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_neighbor_mask_S1 = tox_bool_in(neighbor_mask_S1, "neighbor_mask_S1");
    unsigned char* _c_neighbor_mask_S2 = tox_bool_in(neighbor_mask_S2, "neighbor_mask_S2");

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    int total_included_n_reps = 0;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    SEXP pmf_S1 = PROTECT(Rf_allocVector(REALSXP, n_points * _v_n_bins)); _nprot++;
    { SEXP _dim_pmf_S1 = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_pmf_S1)[0] = n_points; INTEGER(_dim_pmf_S1)[1] = _v_n_bins; Rf_setAttrib(pmf_S1, R_DimSymbol, _dim_pmf_S1); UNPROTECT(1); }
    SEXP pmf_S2 = PROTECT(Rf_allocVector(REALSXP, n_points * _v_n_bins)); _nprot++;
    { SEXP _dim_pmf_S2 = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_pmf_S2)[0] = n_points; INTEGER(_dim_pmf_S2)[1] = _v_n_bins; Rf_setAttrib(pmf_S2, R_DimSymbol, _dim_pmf_S2); UNPROTECT(1); }
    int ierr = 0;

    fjct_compute_masked_jsd_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        _c_neighbor_mask_S1,
        _c_neighbor_mask_S2,
        &_v_n_bins,
        &_v_shared_residual_range,
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

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); _nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(total_included_n_reps));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 5, weights);
    SET_VECTOR_ELT(_out, 6, pmf_S1);
    SET_VECTOR_ELT(_out, 7, pmf_S2);
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); _nprot++;
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
    UNPROTECT(_nprot);
    return _out;
}

SEXP fjct_compute_contribution_scores_call(SEXP global_js_divergences, SEXP total_included_n_reps_per_f) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int k_families = (int) Rf_length(global_js_divergences);

    // convert what Fortran cannot take from R directly
    const int* _c_total_included_n_reps_per_f = tox_int_in(total_included_n_reps_per_f, "total_included_n_reps_per_f");

    // outputs and work space
    SEXP support_weights = PROTECT(Rf_allocVector(REALSXP, k_families)); _nprot++;
    SEXP contribution_scores = PROTECT(Rf_allocVector(REALSXP, k_families)); _nprot++;
    int ierr = 0;

    fjct_compute_contribution_scores_c(
        REAL(global_js_divergences),
        _c_total_included_n_reps_per_f,
        &k_families,
        REAL(support_weights),
        REAL(contribution_scores),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, support_weights);
    SET_VECTOR_ELT(_out, 1, contribution_scores);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("support_weights"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("contribution_scores"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
