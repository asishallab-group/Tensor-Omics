// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void compute_gene_means_c(const int*, const int*, const double*, double*, int*);
void compute_residuals_c(const int*, const int*, const double*, const double*, double*, int*);
void pool_means_c(const int*, const double*, const int*, const double*, const int*, int*, double*, int*);
void pool_means_expert_c(const double*, const int*, const int*, const int*, int*, double*, int*);
void calc_neighborhood_size_c(const int*, const int*, const int*, const double*, const int*, int*, int*);
void construct_neighborhoods_c(const int*, const double*, const int*, const double*, const int*, const double*, double*, int*, const int*, int*);
void gjct_permutation_test_c(const double*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, int*, const int*, const unsigned char*, const unsigned char*);
void gjct_permutation_test_expert_c(double*, double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, double*, double*, double*, int*, int*, int*, double*, double*, int*, const int*, const unsigned char*, const unsigned char*);
void determine_shared_residual_range_expert_c(const double*, const int*, const int*, double*, int*, const double*);
void determine_shared_residual_range_c(const double*, const double*, const int*, const int*, const int*, const int*, double*, int*, const double*);
void build_residual_histograms_c(const double*, const int*, const int*, const int*, const double*, const int*, int*, double*, int*, int*, const unsigned char*);
void compute_divergence_per_reference_point_c(const double*, const double*, const int*, const int*, double*, int*);
void compute_weighted_global_divergence_c(const double*, const int*, const int*, const int*, double*, double*, int*);
void fjct_compute_jsd_c(const int*, const int*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*, int*, int*, double*, double*, int*);
void fjct_compute_jsd_expert_c(const double*, const double*, const int*, const int*, const int*, const int*, const unsigned char*, const unsigned char*, const int*, const double*, double*, int*, int*, int*, double*, double*, double*, double*, int*, int*);
void fjct_compute_contribution_scores_c(const double*, const int*, const int*, double*, double*, int*);

SEXP compute_gene_means_call(SEXP expr) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_reps = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP means = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_gene_means_c(
        &n_genes,
        &n_reps,
        REAL(expr),
        REAL(means),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, means);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("means"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_residuals_call(SEXP expr, SEXP means) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[1];
    int n_reps = INTEGER(Rf_getAttrib(expr, R_DimSymbol))[0];

    // outputs and work space
    SEXP resid = PROTECT(Rf_allocVector(REALSXP, n_reps * n_genes)); nprot++;
    { SEXP resid_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(resid_dim)[0] = n_reps; INTEGER(resid_dim)[1] = n_genes; Rf_setAttrib(resid, R_DimSymbol, resid_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_residuals_c(
        &n_genes,
        &n_reps,
        REAL(expr),
        REAL(means),
        REAL(resid),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, resid);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("resid"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP pool_means_call(SEXP mean_S1, SEXP mean_S2, SEXP n_points) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S1 = (int) Rf_length(mean_S1);
    int n_genes_S2 = (int) Rf_length(mean_S2);

    // scalar inputs, pulled from their length-1 vectors
    int n_points_v = Rf_asInteger(n_points);

    // outputs and work space
    int n_pool = 0;
    SEXP x_star = PROTECT(Rf_allocVector(REALSXP, n_points_v)); nprot++;
    int ierr = 0;

    pool_means_c(
        &n_genes_S1,
        REAL(mean_S1),
        &n_genes_S2,
        REAL(mean_S2),
        &n_points_v,
        &n_pool,
        REAL(x_star),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_pool));
    SET_VECTOR_ELT(_out, 1, x_star);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_pool"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("x_star"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP pool_means_expert_call(SEXP pooled_means, SEXP pooled_means_perm, SEXP n_points) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(pooled_means);

    // scalar inputs, pulled from their length-1 vectors
    int n_points_v = Rf_asInteger(n_points);

    // outputs and work space
    int n_pool = 0;
    SEXP x_star = PROTECT(Rf_allocVector(REALSXP, n_points_v)); nprot++;
    int ierr = 0;

    pool_means_expert_c(
        REAL(pooled_means),
        INTEGER(pooled_means_perm),
        &pool_size,
        &n_points_v,
        &n_pool,
        REAL(x_star),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_pool));
    SET_VECTOR_ELT(_out, 1, x_star);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_pool"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("x_star"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP calc_neighborhood_size_call(SEXP n_pool, SEXP n_points, SEXP mean_S, SEXP desired_size) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes_S = (int) Rf_length(mean_S);

    // scalar inputs, pulled from their length-1 vectors
    int n_pool_v = Rf_asInteger(n_pool);
    int n_points_v = Rf_asInteger(n_points);
    int desired_size_v = Rf_asInteger(desired_size);

    // outputs and work space
    int n_neighbors = 0;
    int ierr = 0;

    calc_neighborhood_size_c(
        &n_pool_v,
        &n_points_v,
        &n_genes_S,
        REAL(mean_S),
        &desired_size_v,
        &n_neighbors,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_neighbors));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP construct_neighborhoods_call(SEXP x_star, SEXP mean_S, SEXP resid_S, SEXP n_neighbors) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = (int) Rf_length(x_star);
    int n_genes_S = (int) Rf_length(mean_S);
    int n_reps_S = INTEGER(Rf_getAttrib(resid_S, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int n_neighbors_v = Rf_asInteger(n_neighbors);

    // outputs and work space
    SEXP neighborhood_residuals = PROTECT(Rf_allocVector(REALSXP, n_reps_S * n_neighbors_v * n_points)); nprot++;
    { SEXP neighborhood_residuals_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(neighborhood_residuals_dim)[0] = n_reps_S; INTEGER(neighborhood_residuals_dim)[1] = n_neighbors_v; INTEGER(neighborhood_residuals_dim)[2] = n_points; Rf_setAttrib(neighborhood_residuals, R_DimSymbol, neighborhood_residuals_dim); UNPROTECT(1); }
    SEXP neighborhood_indices = PROTECT(Rf_allocVector(INTSXP, n_neighbors_v * n_points)); nprot++;
    { SEXP neighborhood_indices_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(neighborhood_indices_dim)[0] = n_neighbors_v; INTEGER(neighborhood_indices_dim)[1] = n_points; Rf_setAttrib(neighborhood_indices, R_DimSymbol, neighborhood_indices_dim); UNPROTECT(1); }
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
        &n_neighbors_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, neighborhood_residuals);
    SET_VECTOR_ELT(_out, 1, neighborhood_indices);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighborhood_residuals"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("neighborhood_indices"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP gjct_permutation_test_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP global_jsd_observed, SEXP n_bins, SEXP shared_residual_range, SEXP n_permutations, SEXP random_seed, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = NULL;
    int random_seed_size = 0;
    if (random_seed != R_NilValue) {
        random_seed_size = (int) Rf_length(random_seed);
        random_seed_p = INTEGER(random_seed);
    }
    int neighbor_mask_S1_size = 0;
    int neighbor_mask_S1_dim[2]; for (int i = 0; i < 2; ++i) neighbor_mask_S1_dim[i] = 0;
    if (neighbor_mask_S1 != R_NilValue) {
        neighbor_mask_S1_size = (int) Rf_length(neighbor_mask_S1);
        SEXP neighbor_mask_S1_d = Rf_getAttrib(neighbor_mask_S1, R_DimSymbol);
        if (neighbor_mask_S1_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(neighbor_mask_S1_d); ++i) neighbor_mask_S1_dim[i] = INTEGER(neighbor_mask_S1_d)[i];
        }
    }
    int neighbor_mask_S2_size = 0;
    int neighbor_mask_S2_dim[2]; for (int i = 0; i < 2; ++i) neighbor_mask_S2_dim[i] = 0;
    if (neighbor_mask_S2 != R_NilValue) {
        neighbor_mask_S2_size = (int) Rf_length(neighbor_mask_S2);
        SEXP neighbor_mask_S2_d = Rf_getAttrib(neighbor_mask_S2, R_DimSymbol);
        if (neighbor_mask_S2_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(neighbor_mask_S2_d); ++i) neighbor_mask_S2_dim[i] = INTEGER(neighbor_mask_S2_d)[i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double global_jsd_observed_v = Rf_asReal(global_jsd_observed);
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int n_permutations_v = Rf_asInteger(n_permutations);

    // convert what Fortran cannot take from R directly
    unsigned char* neighbor_mask_S1_c = tox_bool_in(neighbor_mask_S1);
    unsigned char* neighbor_mask_S2_c = tox_bool_in(neighbor_mask_S2);

    // outputs and work space
    SEXP jsd_null = PROTECT(Rf_allocVector(REALSXP, n_permutations_v)); nprot++;
    double p_value = 0;
    int ierr = 0;

    gjct_permutation_test_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &global_jsd_observed_v,
        &n_bins_v,
        &shared_residual_range_v,
        &n_permutations_v,
        REAL(jsd_null),
        &p_value,
        &ierr,
        random_seed_p,
        neighbor_mask_S1 != R_NilValue ? neighbor_mask_S1_c : NULL,
        neighbor_mask_S2 != R_NilValue ? neighbor_mask_S2_c : NULL
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, jsd_null);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(p_value));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("jsd_null"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("p_value"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP gjct_permutation_test_expert_call(SEXP neighborhood_residuals_S1_copy, SEXP neighborhood_residuals_S2_copy, SEXP global_jsd_observed, SEXP n_bins, SEXP shared_residual_range, SEXP n_permutations, SEXP random_seed, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = NULL;
    int random_seed_size = 0;
    if (random_seed != R_NilValue) {
        random_seed_size = (int) Rf_length(random_seed);
        random_seed_p = INTEGER(random_seed);
    }
    int neighbor_mask_S1_size = 0;
    int neighbor_mask_S1_dim[2]; for (int i = 0; i < 2; ++i) neighbor_mask_S1_dim[i] = 0;
    if (neighbor_mask_S1 != R_NilValue) {
        neighbor_mask_S1_size = (int) Rf_length(neighbor_mask_S1);
        SEXP neighbor_mask_S1_d = Rf_getAttrib(neighbor_mask_S1, R_DimSymbol);
        if (neighbor_mask_S1_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(neighbor_mask_S1_d); ++i) neighbor_mask_S1_dim[i] = INTEGER(neighbor_mask_S1_d)[i];
        }
    }
    int neighbor_mask_S2_size = 0;
    int neighbor_mask_S2_dim[2]; for (int i = 0; i < 2; ++i) neighbor_mask_S2_dim[i] = 0;
    if (neighbor_mask_S2 != R_NilValue) {
        neighbor_mask_S2_size = (int) Rf_length(neighbor_mask_S2);
        SEXP neighbor_mask_S2_d = Rf_getAttrib(neighbor_mask_S2, R_DimSymbol);
        if (neighbor_mask_S2_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(neighbor_mask_S2_d); ++i) neighbor_mask_S2_dim[i] = INTEGER(neighbor_mask_S2_d)[i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1_copy, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2_copy, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1_copy, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1_copy, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double global_jsd_observed_v = Rf_asReal(global_jsd_observed);
    int n_bins_v = Rf_asInteger(n_bins);
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int n_permutations_v = Rf_asInteger(n_permutations);

    // copy what is modified in place, so the caller's stays intact
    SEXP neighborhood_residuals_S1_copy_out = PROTECT(Rf_duplicate(neighborhood_residuals_S1_copy)); nprot++;
    SEXP neighborhood_residuals_S2_copy_out = PROTECT(Rf_duplicate(neighborhood_residuals_S2_copy)); nprot++;

    // convert what Fortran cannot take from R directly
    unsigned char* neighbor_mask_S1_c = tox_bool_in(neighbor_mask_S1);
    unsigned char* neighbor_mask_S2_c = tox_bool_in(neighbor_mask_S2);

    // outputs and work space
    SEXP jsd_null = PROTECT(Rf_allocVector(REALSXP, n_permutations_v)); nprot++;
    double p_value = 0;
    double* tmp_pool = (double*) R_alloc((n_reps_S1 + n_reps_S2) * n_neighbors, sizeof(double));
    double* tmp_pmf_S1 = (double*) R_alloc(n_points * n_bins_v, sizeof(double));
    double* tmp_pmf_S2 = (double*) R_alloc(n_points * n_bins_v, sizeof(double));
    int* tmp_counts = (int*) R_alloc(n_points * n_bins_v, sizeof(int));
    int* tmp_included_n_reps_S1 = (int*) R_alloc(n_points, sizeof(int));
    int* tmp_included_n_reps_S2 = (int*) R_alloc(n_points, sizeof(int));
    double* tmp_js_divergences = (double*) R_alloc(n_points, sizeof(double));
    double* tmp_weights = (double*) R_alloc(n_points, sizeof(double));
    int ierr = 0;

    gjct_permutation_test_expert_c(
        REAL(neighborhood_residuals_S1_copy_out),
        REAL(neighborhood_residuals_S2_copy_out),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &global_jsd_observed_v,
        &n_bins_v,
        &shared_residual_range_v,
        &n_permutations_v,
        REAL(jsd_null),
        &p_value,
        tmp_pool,
        tmp_pmf_S1,
        tmp_pmf_S2,
        tmp_counts,
        tmp_included_n_reps_S1,
        tmp_included_n_reps_S2,
        tmp_js_divergences,
        tmp_weights,
        &ierr,
        random_seed_p,
        neighbor_mask_S1 != R_NilValue ? neighbor_mask_S1_c : NULL,
        neighbor_mask_S2 != R_NilValue ? neighbor_mask_S2_c : NULL
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, neighborhood_residuals_S1_copy_out);
    SET_VECTOR_ELT(_out, 1, neighborhood_residuals_S2_copy_out);
    SET_VECTOR_ELT(_out, 2, jsd_null);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(p_value));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighborhood_residuals_S1_copy"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("neighborhood_residuals_S2_copy"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("jsd_null"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("p_value"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP determine_shared_residual_range_expert_call(SEXP abs_residual_pool, SEXP abs_residual_pool_perm, SEXP residual_range_quantile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(abs_residual_pool);

    // scalar inputs, pulled from their length-1 vectors
    double residual_range_quantile_v = Rf_asReal(residual_range_quantile);

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_expert_c(
        REAL(abs_residual_pool),
        INTEGER(abs_residual_pool_perm),
        &pool_size,
        &shared_residual_range,
        &ierr,
        &residual_range_quantile_v
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(shared_residual_range));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("shared_residual_range"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP determine_shared_residual_range_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP residual_range_quantile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double residual_range_quantile_v = Rf_asReal(residual_range_quantile);

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &shared_residual_range,
        &ierr,
        &residual_range_quantile_v
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(shared_residual_range));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("shared_residual_range"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP build_residual_histograms_call(SEXP neighborhood_residuals, SEXP shared_residual_range, SEXP n_bins, SEXP neighbor_mask) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int neighbor_mask_size = 0;
    int neighbor_mask_dim[2]; for (int i = 0; i < 2; ++i) neighbor_mask_dim[i] = 0;
    if (neighbor_mask != R_NilValue) {
        neighbor_mask_size = (int) Rf_length(neighbor_mask);
        SEXP neighbor_mask_d = Rf_getAttrib(neighbor_mask, R_DimSymbol);
        if (neighbor_mask_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(neighbor_mask_d); ++i) neighbor_mask_dim[i] = INTEGER(neighbor_mask_d)[i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps = INTEGER(Rf_getAttrib(neighborhood_residuals, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double shared_residual_range_v = Rf_asReal(shared_residual_range);
    int n_bins_v = Rf_asInteger(n_bins);

    // convert what Fortran cannot take from R directly
    unsigned char* neighbor_mask_c = tox_bool_in(neighbor_mask);

    // outputs and work space
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, n_points * n_bins_v)); nprot++;
    { SEXP counts_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(counts_dim)[0] = n_points; INTEGER(counts_dim)[1] = n_bins_v; Rf_setAttrib(counts, R_DimSymbol, counts_dim); UNPROTECT(1); }
    SEXP pmf = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins_v)); nprot++;
    { SEXP pmf_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(pmf_dim)[0] = n_points; INTEGER(pmf_dim)[1] = n_bins_v; Rf_setAttrib(pmf, R_DimSymbol, pmf_dim); UNPROTECT(1); }
    SEXP included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int ierr = 0;

    build_residual_histograms_c(
        REAL(neighborhood_residuals),
        &n_reps,
        &n_neighbors,
        &n_points,
        &shared_residual_range_v,
        &n_bins_v,
        INTEGER(counts),
        REAL(pmf),
        INTEGER(included_n_reps),
        &ierr,
        neighbor_mask != R_NilValue ? neighbor_mask_c : NULL
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, counts);
    SET_VECTOR_ELT(_out, 1, pmf);
    SET_VECTOR_ELT(_out, 2, included_n_reps);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("counts"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("pmf"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_divergence_per_reference_point_call(SEXP pmf_S1, SEXP pmf_S2) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(pmf_S1, R_DimSymbol))[0];
    int n_bins = INTEGER(Rf_getAttrib(pmf_S1, R_DimSymbol))[1];

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    int ierr = 0;

    compute_divergence_per_reference_point_c(
        REAL(pmf_S1),
        REAL(pmf_S2),
        &n_points,
        &n_bins,
        REAL(js_divergences),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_weighted_global_divergence_call(SEXP js_divergences, SEXP included_n_reps_S1, SEXP included_n_reps_S2) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = (int) Rf_length(js_divergences);

    // outputs and work space
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); nprot++;
    int ierr = 0;

    compute_weighted_global_divergence_c(
        REAL(js_divergences),
        &n_points,
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &global_js_divergence,
        REAL(weights),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 1, weights);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
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

SEXP fjct_compute_jsd_expert_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2, SEXP n_bins, SEXP shared_residual_range) {
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

    fjct_compute_jsd_expert_c(
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
