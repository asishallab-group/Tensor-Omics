// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void gjct_permutation_test_expert_c(const double*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, double*, double*, double*, double*, double*, int*, int*, int*, double*, double*, const int*, const unsigned char*, const unsigned char*, int*);
void gjct_permutation_test_c(const double*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, const int*, const unsigned char*, const unsigned char*, int*);

SEXP gjct_permutation_test_expert_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP global_jsd_observed, SEXP n_bins, SEXP shared_residual_range, SEXP n_permutations, SEXP random_seed, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2) {
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
    double* tmp_residuals_S1 = (double*) R_alloc(n_reps_S1 * n_neighbors * n_points, sizeof(double));
    double* tmp_residuals_S2 = (double*) R_alloc(n_reps_S2 * n_neighbors * n_points, sizeof(double));
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
        tmp_residuals_S1,
        tmp_residuals_S2,
        tmp_pool,
        tmp_pmf_S1,
        tmp_pmf_S2,
        tmp_counts,
        tmp_included_n_reps_S1,
        tmp_included_n_reps_S2,
        tmp_js_divergences,
        tmp_weights,
        random_seed_p,
        neighbor_mask_S1 != R_NilValue ? neighbor_mask_S1_c : NULL,
        neighbor_mask_S2 != R_NilValue ? neighbor_mask_S2_c : NULL,
        &ierr
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
        random_seed_p,
        neighbor_mask_S1 != R_NilValue ? neighbor_mask_S1_c : NULL,
        neighbor_mask_S2 != R_NilValue ? neighbor_mask_S2_c : NULL,
        &ierr
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

#endif  // R binding
