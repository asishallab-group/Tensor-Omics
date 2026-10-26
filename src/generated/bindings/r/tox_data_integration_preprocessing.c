// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void compute_gene_means_c(const int*, const int*, const double*, double*, int*);
void compute_residuals_c(const int*, const int*, const double*, const double*, double*, int*);
void pool_means_c(const double*, const int*, const int*, int*, double*, int*);
void pool_means_expert_c(const double*, const int*, const int*, const int*, int*, double*, int*);
void pool_study_means_c(const int*, const double*, const int*, const double*, const int*, int*, double*, int*);
void construct_neighborhoods_c(const int*, const double*, const int*, const double*, const int*, const double*, double*, int*, const int*, int*);

SEXP compute_gene_means_call(SEXP expr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_reps = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP means = PROTECT(Rf_allocVector(REALSXP, n_genes)); _nprot++;
    int ierr = 0;

    compute_gene_means_c(
        &n_genes,
        &n_reps,
        REAL(expr),
        REAL(means),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, means);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("means"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_residuals_call(SEXP expr, SEXP means) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_reps = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP resid = PROTECT(Rf_allocVector(REALSXP, n_reps * n_genes)); _nprot++;
    { SEXP _dim_resid = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_resid)[0] = n_reps; INTEGER(_dim_resid)[1] = n_genes; Rf_setAttrib(resid, R_DimSymbol, _dim_resid); UNPROTECT(1); }
    int ierr = 0;

    compute_residuals_c(
        &n_genes,
        &n_reps,
        REAL(expr),
        REAL(means),
        REAL(resid),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, resid);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("resid"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP pool_means_call(SEXP pooled_means, SEXP n_points) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(pooled_means);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_points = Rf_asInteger(n_points);

    // outputs and work space
    int n_pool = 0;
    SEXP x_star = PROTECT(Rf_allocVector(REALSXP, _v_n_points)); _nprot++;
    int ierr = 0;

    pool_means_c(
        REAL(pooled_means),
        &pool_size,
        &_v_n_points,
        &n_pool,
        REAL(x_star),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_pool));
    SET_VECTOR_ELT(_out, 1, x_star);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_pool"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("x_star"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP pool_means_expert_call(SEXP pooled_means, SEXP pooled_means_perm, SEXP n_points) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(pooled_means);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_points = Rf_asInteger(n_points);

    // outputs and work space
    int n_pool = 0;
    SEXP x_star = PROTECT(Rf_allocVector(REALSXP, _v_n_points)); _nprot++;
    int ierr = 0;

    pool_means_expert_c(
        REAL(pooled_means),
        INTEGER(pooled_means_perm),
        &pool_size,
        &_v_n_points,
        &n_pool,
        REAL(x_star),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_pool));
    SET_VECTOR_ELT(_out, 1, x_star);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_pool"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("x_star"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP pool_study_means_call(SEXP mean_S1, SEXP mean_S2, SEXP n_points) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) Rf_length(mean_S1);
    int n_genes_S2 = (int) Rf_length(mean_S2);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_points = Rf_asInteger(n_points);

    // outputs and work space
    int n_pool = 0;
    SEXP x_star = PROTECT(Rf_allocVector(REALSXP, _v_n_points)); _nprot++;
    int ierr = 0;

    pool_study_means_c(
        &n_genes_S1,
        REAL(mean_S1),
        &n_genes_S2,
        REAL(mean_S2),
        &_v_n_points,
        &n_pool,
        REAL(x_star),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_pool));
    SET_VECTOR_ELT(_out, 1, x_star);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_pool"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("x_star"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP construct_neighborhoods_call(SEXP x_star, SEXP mean_S, SEXP resid_S, SEXP n_neighbors) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = (int) Rf_length(x_star);
    int n_genes_S = (int) Rf_length(mean_S);
    int n_reps_S = INTEGER(Rf_getAttrib(resid_S, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_neighbors = Rf_asInteger(n_neighbors);

    // outputs and work space
    SEXP neighborhood_residuals = PROTECT(Rf_allocVector(REALSXP, n_reps_S * _v_n_neighbors * n_points)); _nprot++;
    { SEXP _dim_neighborhood_residuals = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_neighborhood_residuals)[0] = n_reps_S; INTEGER(_dim_neighborhood_residuals)[1] = _v_n_neighbors; INTEGER(_dim_neighborhood_residuals)[2] = n_points; Rf_setAttrib(neighborhood_residuals, R_DimSymbol, _dim_neighborhood_residuals); UNPROTECT(1); }
    SEXP neighborhood_indices = PROTECT(Rf_allocVector(INTSXP, _v_n_neighbors * n_points)); _nprot++;
    { SEXP _dim_neighborhood_indices = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_neighborhood_indices)[0] = _v_n_neighbors; INTEGER(_dim_neighborhood_indices)[1] = n_points; Rf_setAttrib(neighborhood_indices, R_DimSymbol, _dim_neighborhood_indices); UNPROTECT(1); }
    int ierr = 0;

    construct_neighborhoods_c(
        &n_points,
        REAL(x_star),
        &n_genes_S,
        REAL(mean_S),
        &n_reps_S,
        REAL(resid_S),
        REAL(neighborhood_residuals),
        INTEGER(neighborhood_indices),
        &_v_n_neighbors,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, neighborhood_residuals);
    SET_VECTOR_ELT(_out, 1, neighborhood_indices);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighborhood_residuals"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("neighborhood_indices"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
