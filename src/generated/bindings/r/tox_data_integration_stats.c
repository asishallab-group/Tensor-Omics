// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void gjct_permutation_test_c(const double*, const double*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const int*, double*, double*, const int*, const unsigned char*, const unsigned char*, int*);

SEXP gjct_permutation_test_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP global_jsd_observed, SEXP n_bins, SEXP shared_residual_range, SEXP n_permutations, SEXP random_seed, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2) {
    int _nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* _p_random_seed = NULL;
    int _size_random_seed = 0;
    if (random_seed != R_NilValue) {
        _size_random_seed = (int) Rf_length(random_seed);
        _p_random_seed = INTEGER(random_seed);
    }
    int _size_neighbor_mask_S1 = 0;
    int _dim_neighbor_mask_S1[2]; for (int _i = 0; _i < 2; ++_i) _dim_neighbor_mask_S1[_i] = 0;
    if (neighbor_mask_S1 != R_NilValue) {
        _size_neighbor_mask_S1 = (int) Rf_length(neighbor_mask_S1);
        SEXP _d_neighbor_mask_S1 = Rf_getAttrib(neighbor_mask_S1, R_DimSymbol);
        if (_d_neighbor_mask_S1 != R_NilValue) {
            for (int _i = 0; _i < 2 && _i < (int) Rf_length(_d_neighbor_mask_S1); ++_i) _dim_neighbor_mask_S1[_i] = INTEGER(_d_neighbor_mask_S1)[_i];
        }
    }
    int _size_neighbor_mask_S2 = 0;
    int _dim_neighbor_mask_S2[2]; for (int _i = 0; _i < 2; ++_i) _dim_neighbor_mask_S2[_i] = 0;
    if (neighbor_mask_S2 != R_NilValue) {
        _size_neighbor_mask_S2 = (int) Rf_length(neighbor_mask_S2);
        SEXP _d_neighbor_mask_S2 = Rf_getAttrib(neighbor_mask_S2, R_DimSymbol);
        if (_d_neighbor_mask_S2 != R_NilValue) {
            for (int _i = 0; _i < 2 && _i < (int) Rf_length(_d_neighbor_mask_S2); ++_i) _dim_neighbor_mask_S2[_i] = INTEGER(_d_neighbor_mask_S2)[_i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double _v_global_jsd_observed = Rf_asReal(global_jsd_observed);
    int _v_n_bins = Rf_asInteger(n_bins);
    double _v_shared_residual_range = Rf_asReal(shared_residual_range);
    int _v_n_permutations = Rf_asInteger(n_permutations);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_neighbor_mask_S1 = tox_bool_in(neighbor_mask_S1, "neighbor_mask_S1");
    unsigned char* _c_neighbor_mask_S2 = tox_bool_in(neighbor_mask_S2, "neighbor_mask_S2");

    // outputs and work space
    SEXP jsd_null = PROTECT(Rf_allocVector(REALSXP, _v_n_permutations)); _nprot++;
    double p_value = 0;
    int ierr = 0;

    gjct_permutation_test_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &_v_global_jsd_observed,
        &_v_n_bins,
        &_v_shared_residual_range,
        &_v_n_permutations,
        REAL(jsd_null),
        &p_value,
        _p_random_seed,
        neighbor_mask_S1 != R_NilValue ? _c_neighbor_mask_S1 : NULL,
        neighbor_mask_S2 != R_NilValue ? _c_neighbor_mask_S2 : NULL,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, jsd_null);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(p_value));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("jsd_null"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("p_value"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
