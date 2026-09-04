// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void determine_shared_residual_range_c(const double*, const int*, double*, const double*, int*);
void determine_shared_residual_range_expert_c(const double*, const int*, const int*, double*, const double*, int*);
void determine_study_shared_residual_range_c(const double*, const double*, const int*, const int*, const int*, const int*, double*, const double*, int*);
void build_residual_histograms_c(const double*, const int*, const int*, const int*, const double*, const int*, int*, double*, int*, const unsigned char*, int*);
void compute_divergence_per_reference_point_c(const double*, const double*, const int*, const int*, double*, int*);
void compute_weighted_global_divergence_c(const double*, const int*, const int*, const int*, double*, double*, int*);

SEXP determine_shared_residual_range_call(SEXP abs_residual_pool, SEXP residual_range_quantile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(abs_residual_pool);

    // scalar inputs, pulled from their length-1 vectors
    double residual_range_quantile_v = Rf_asReal(residual_range_quantile);

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_c(
        REAL(abs_residual_pool),
        &pool_size,
        &shared_residual_range,
        &residual_range_quantile_v,
        &ierr
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
        &residual_range_quantile_v,
        &ierr
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

SEXP determine_study_shared_residual_range_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP residual_range_quantile) {
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

    determine_study_shared_residual_range_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &shared_residual_range,
        &residual_range_quantile_v,
        &ierr
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
        neighbor_mask != R_NilValue ? neighbor_mask_c : NULL,
        &ierr
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

#endif  // R binding
