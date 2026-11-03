// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void determine_shared_residual_range_c(const double*, const int*, double*, const double*, int*);
void determine_shared_residual_range_expert_c(const double*, const int*, const int*, double*, const double*, int*);
void determine_study_shared_residual_range_c(const double*, const double*, const int*, const int*, const int*, const int*, double*, const double*, int*);
void build_residual_histograms_c(const double*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, int*, double*, int*, const unsigned char*, int*);
void calc_pmf_c(const int*, const int*, const int*, const int*, double*, int*);
void compute_divergence_per_reference_point_c(const double*, const double*, const int*, const int*, double*, int*);
void compute_weighted_global_divergence_c(const double*, const int*, const int*, const int*, double*, double*, int*);
void jct_compute_jsd_pipeline_c(const double*, const double*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*, int*, double*, double*, double*, double*, const unsigned char*, const unsigned char*, int*);
void create_mean_pmf_c(const double*, const int*, const int*, const int*, const int*, const int*, double*, int*, int*, int*);
void create_mean_pmf_only_c(const double*, const int*, const int*, const int*, double*, int*);

SEXP determine_shared_residual_range_call(SEXP abs_residual_pool, SEXP residual_range_quantile_level) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(abs_residual_pool);

    // scalar inputs, pulled from their length-1 vectors
    double _v_residual_range_quantile_level = Rf_asReal(residual_range_quantile_level);

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_c(
        REAL(abs_residual_pool),
        &pool_size,
        &shared_residual_range,
        &_v_residual_range_quantile_level,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(shared_residual_range));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("shared_residual_range"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_shared_residual_range_expert_call(SEXP abs_residual_pool, SEXP abs_residual_pool_perm, SEXP residual_range_quantile_level) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int pool_size = (int) Rf_length(abs_residual_pool);

    // scalar inputs, pulled from their length-1 vectors
    double _v_residual_range_quantile_level = Rf_asReal(residual_range_quantile_level);

    // convert what Fortran cannot take from R directly
    const int* _c_abs_residual_pool_perm = tox_int_in(abs_residual_pool_perm, "abs_residual_pool_perm");

    // outputs and work space
    double shared_residual_range = 0;
    int ierr = 0;

    determine_shared_residual_range_expert_c(
        REAL(abs_residual_pool),
        _c_abs_residual_pool_perm,
        &pool_size,
        &shared_residual_range,
        &_v_residual_range_quantile_level,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(shared_residual_range));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("shared_residual_range"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_study_shared_residual_range_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP residual_range_quantile_level) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_reps_S1 = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[0];
    int n_reps_S2 = INTEGER(Rf_getAttrib(neighborhood_residuals_S2, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals_S1, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    double _v_residual_range_quantile_level = Rf_asReal(residual_range_quantile_level);

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
        &_v_residual_range_quantile_level,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(shared_residual_range));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("shared_residual_range"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP build_residual_histograms_call(SEXP neighborhood_residuals, SEXP shared_residual_range_low, SEXP shared_residual_range_high, SEXP max_n_bins, SEXP n_bins_per_point, SEXP neighbor_mask) {
    int _nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int _size_neighbor_mask = 0;
    int _dim_neighbor_mask[2]; for (int _i = 0; _i < 2; ++_i) _dim_neighbor_mask[_i] = 0;
    if (neighbor_mask != R_NilValue) {
        _size_neighbor_mask = (int) Rf_length(neighbor_mask);
        SEXP _d_neighbor_mask = Rf_getAttrib(neighbor_mask, R_DimSymbol);
        if (_d_neighbor_mask != R_NilValue) {
            for (int _i = 0; _i < 2 && _i < (int) Rf_length(_d_neighbor_mask); ++_i) _dim_neighbor_mask[_i] = INTEGER(_d_neighbor_mask)[_i];
        }
    }

    // derived from the inputs, not asked of the caller
    int n_reps = INTEGER(Rf_getAttrib(neighborhood_residuals, R_DimSymbol))[0];
    int n_neighbors = INTEGER(Rf_getAttrib(neighborhood_residuals, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(neighborhood_residuals, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_bins = Rf_asInteger(max_n_bins);

    // convert what Fortran cannot take from R directly
    const int* _c_n_bins_per_point = tox_int_in(n_bins_per_point, "n_bins_per_point");
    unsigned char* _c_neighbor_mask = tox_bool_in(neighbor_mask, "neighbor_mask");

    // outputs and work space
    SEXP counts = PROTECT(Rf_allocVector(INTSXP, n_points * _v_max_n_bins)); _nprot++;
    { SEXP _dim_counts = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_counts)[0] = n_points; INTEGER(_dim_counts)[1] = _v_max_n_bins; Rf_setAttrib(counts, R_DimSymbol, _dim_counts); UNPROTECT(1); }
    SEXP pmf = PROTECT(Rf_allocVector(REALSXP, n_points * _v_max_n_bins)); _nprot++;
    { SEXP _dim_pmf = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_pmf)[0] = n_points; INTEGER(_dim_pmf)[1] = _v_max_n_bins; Rf_setAttrib(pmf, R_DimSymbol, _dim_pmf); UNPROTECT(1); }
    SEXP included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    int ierr = 0;

    build_residual_histograms_c(
        REAL(neighborhood_residuals),
        &n_reps,
        &n_neighbors,
        &n_points,
        REAL(shared_residual_range_low),
        REAL(shared_residual_range_high),
        &_v_max_n_bins,
        _c_n_bins_per_point,
        INTEGER(counts),
        REAL(pmf),
        INTEGER(included_n_reps),
        neighbor_mask != R_NilValue ? _c_neighbor_mask : NULL,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, counts);
    SET_VECTOR_ELT(_out, 1, pmf);
    SET_VECTOR_ELT(_out, 2, included_n_reps);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("counts"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("pmf"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_pmf_call(SEXP counts, SEXP included_n_reps) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(counts, R_DimSymbol))[0];
    int n_bins = INTEGER(Rf_getAttrib(counts, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    const int* _c_counts = tox_int_in(counts, "counts");
    const int* _c_included_n_reps = tox_int_in(included_n_reps, "included_n_reps");

    // outputs and work space
    SEXP pmf = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins)); _nprot++;
    { SEXP _dim_pmf = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_pmf)[0] = n_points; INTEGER(_dim_pmf)[1] = n_bins; Rf_setAttrib(pmf, R_DimSymbol, _dim_pmf); UNPROTECT(1); }
    int ierr = 0;

    calc_pmf_c(
        _c_counts,
        _c_included_n_reps,
        &n_points,
        &n_bins,
        REAL(pmf),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, pmf);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("pmf"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_divergence_per_reference_point_call(SEXP pmf_S1, SEXP pmf_S2) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(pmf_S1, R_DimSymbol))[0];
    int n_bins = INTEGER(Rf_getAttrib(pmf_S1, R_DimSymbol))[1];

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    int ierr = 0;

    compute_divergence_per_reference_point_c(
        REAL(pmf_S1),
        REAL(pmf_S2),
        &n_points,
        &n_bins,
        REAL(js_divergences),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_weighted_global_divergence_call(SEXP js_divergences, SEXP included_n_reps_S1, SEXP included_n_reps_S2) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = (int) Rf_length(js_divergences);

    // convert what Fortran cannot take from R directly
    const int* _c_included_n_reps_S1 = tox_int_in(included_n_reps_S1, "included_n_reps_S1");
    const int* _c_included_n_reps_S2 = tox_int_in(included_n_reps_S2, "included_n_reps_S2");

    // outputs and work space
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    int ierr = 0;

    compute_weighted_global_divergence_c(
        REAL(js_divergences),
        &n_points,
        _c_included_n_reps_S1,
        _c_included_n_reps_S2,
        &global_js_divergence,
        REAL(weights),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 1, weights);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP jct_compute_jsd_pipeline_call(SEXP neighborhood_residuals_S1, SEXP neighborhood_residuals_S2, SEXP n_bins, SEXP shared_residual_range, SEXP neighbor_mask_S1, SEXP neighbor_mask_S2) {
    int _nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
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
    int _v_n_bins = Rf_asInteger(n_bins);
    double _v_shared_residual_range = Rf_asReal(shared_residual_range);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_neighbor_mask_S1 = tox_bool_in(neighbor_mask_S1, "neighbor_mask_S1");
    unsigned char* _c_neighbor_mask_S2 = tox_bool_in(neighbor_mask_S2, "neighbor_mask_S2");

    // outputs and work space
    SEXP js_divergences = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    SEXP included_n_reps_S1 = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    SEXP included_n_reps_S2 = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    double global_js_divergence = 0;
    SEXP weights = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    SEXP pmf_S1 = PROTECT(Rf_allocVector(REALSXP, n_points * _v_n_bins)); _nprot++;
    { SEXP _dim_pmf_S1 = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_pmf_S1)[0] = n_points; INTEGER(_dim_pmf_S1)[1] = _v_n_bins; Rf_setAttrib(pmf_S1, R_DimSymbol, _dim_pmf_S1); UNPROTECT(1); }
    SEXP pmf_S2 = PROTECT(Rf_allocVector(REALSXP, n_points * _v_n_bins)); _nprot++;
    { SEXP _dim_pmf_S2 = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_pmf_S2)[0] = n_points; INTEGER(_dim_pmf_S2)[1] = _v_n_bins; Rf_setAttrib(pmf_S2, R_DimSymbol, _dim_pmf_S2); UNPROTECT(1); }
    int ierr = 0;

    jct_compute_jsd_pipeline_c(
        REAL(neighborhood_residuals_S1),
        REAL(neighborhood_residuals_S2),
        &n_reps_S1,
        &n_reps_S2,
        &n_neighbors,
        &n_points,
        &_v_n_bins,
        &_v_shared_residual_range,
        REAL(js_divergences),
        INTEGER(included_n_reps_S1),
        INTEGER(included_n_reps_S2),
        &global_js_divergence,
        REAL(weights),
        REAL(pmf_S1),
        REAL(pmf_S2),
        neighbor_mask_S1 != R_NilValue ? _c_neighbor_mask_S1 : NULL,
        neighbor_mask_S2 != R_NilValue ? _c_neighbor_mask_S2 : NULL,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 8)); _nprot++;
    SET_VECTOR_ELT(_out, 0, js_divergences);
    SET_VECTOR_ELT(_out, 1, included_n_reps_S1);
    SET_VECTOR_ELT(_out, 2, included_n_reps_S2);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(global_js_divergence));
    SET_VECTOR_ELT(_out, 4, weights);
    SET_VECTOR_ELT(_out, 5, pmf_S1);
    SET_VECTOR_ELT(_out, 6, pmf_S2);
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 8)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("js_divergences"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("included_n_reps_S1"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("included_n_reps_S2"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("global_js_divergence"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("weights"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("pmf_S1"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("pmf_S2"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP create_mean_pmf_call(SEXP pmfs, SEXP counts, SEXP included_n_reps) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[0];
    int n_studies = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[2];

    // convert what Fortran cannot take from R directly
    const int* _c_counts = tox_int_in(counts, "counts");
    const int* _c_included_n_reps = tox_int_in(included_n_reps, "included_n_reps");

    // outputs and work space
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins)); _nprot++;
    { SEXP _dim_mean_pmf = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_mean_pmf)[0] = n_points; INTEGER(_dim_mean_pmf)[1] = n_bins; Rf_setAttrib(mean_pmf, R_DimSymbol, _dim_mean_pmf); UNPROTECT(1); }
    SEXP mean_pmf_included_n_reps = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    SEXP mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, n_points * n_bins)); _nprot++;
    { SEXP _dim_mean_pmf_counts = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_mean_pmf_counts)[0] = n_points; INTEGER(_dim_mean_pmf_counts)[1] = n_bins; Rf_setAttrib(mean_pmf_counts, R_DimSymbol, _dim_mean_pmf_counts); UNPROTECT(1); }
    int ierr = 0;

    create_mean_pmf_c(
        REAL(pmfs),
        _c_counts,
        &n_bins,
        &n_points,
        &n_studies,
        _c_included_n_reps,
        REAL(mean_pmf),
        INTEGER(mean_pmf_included_n_reps),
        INTEGER(mean_pmf_counts),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, mean_pmf);
    SET_VECTOR_ELT(_out, 1, mean_pmf_included_n_reps);
    SET_VECTOR_ELT(_out, 2, mean_pmf_counts);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("mean_pmf_included_n_reps"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("mean_pmf_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP create_mean_pmf_only_call(SEXP pmfs) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[0];
    int n_studies = INTEGER(Rf_getAttrib(pmfs, R_DimSymbol))[2];

    // outputs and work space
    SEXP mean_pmf = PROTECT(Rf_allocVector(REALSXP, n_points * n_bins)); _nprot++;
    { SEXP _dim_mean_pmf = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_mean_pmf)[0] = n_points; INTEGER(_dim_mean_pmf)[1] = n_bins; Rf_setAttrib(mean_pmf, R_DimSymbol, _dim_mean_pmf); UNPROTECT(1); }
    int ierr = 0;

    create_mean_pmf_only_c(
        REAL(pmfs),
        &n_bins,
        &n_points,
        &n_studies,
        REAL(mean_pmf),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, mean_pmf);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("mean_pmf"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
