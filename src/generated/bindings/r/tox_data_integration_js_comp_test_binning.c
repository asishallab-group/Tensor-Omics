// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void estimate_bin_count_c(const double*, const int*, const int*, const int*, const double*, int*, int*, int*, int*);
void estimate_bin_count_expert_c(const double*, const int*, const int*, const int*, const int*, const double*, int*, int*, int*, int*);
void determine_bin_count_occupancy_c(const double*, const int*, const int*, const int*, int*, unsigned char*, double*, double*, int*, int*, double*, int*, int*, int*, const int*, const int*, const int*, const double*, const double*, const double*, int*);
void determine_bin_count_occupancy_expert_c(const double*, const int*, const int*, const int*, const int*, int*, unsigned char*, double*, double*, int*, int*, double*, int*, int*, int*, int*, const int*, const int*, const int*, const double*, const double*, const double*, int*);
void determine_bin_count_occupancy_exhaustive_c(const double*, const int*, const int*, const double*, const double*, int*, unsigned char*, int*, double*, int*, const int*, const int*, const int*, int*);
void determine_bin_count_occupancy_exhaustive_expert_c(const double*, const int*, const int*, const int*, const double*, const double*, int*, unsigned char*, int*, double*, int*, const int*, const int*, const int*, int*);
void determine_bin_count_stone_c(const double*, const int*, const int*, const int*, const int*, const int*, int*, unsigned char*, double*, double*, int*, int*, double*, int*, int*, int*, const char*, const int*, const int*, const int*, const double*, const double*, int*);
void determine_bin_count_stone_expert_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, int*, unsigned char*, double*, double*, int*, int*, double*, int*, int*, int*, int*, int*, int*, const char*, const int*, const int*, const int*, const double*, const double*, int*);

SEXP estimate_bin_count_call(SEXP residuals, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP shared_residual_range) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_reps_all_studies = Rf_asInteger(max_n_reps_all_studies);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    double _v_shared_residual_range = Rf_asReal(shared_residual_range);

    // outputs and work space
    int n_bins = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    estimate_bin_count_c(
        REAL(residuals),
        &n_residuals,
        &_v_max_n_reps_all_studies,
        &_v_n_neighbors,
        &_v_shared_residual_range,
        &n_bins,
        &sturges_bins,
        &fd_bins,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP estimate_bin_count_expert_call(SEXP residuals, SEXP residuals_perm, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP shared_residual_range) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_reps_all_studies = Rf_asInteger(max_n_reps_all_studies);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    double _v_shared_residual_range = Rf_asReal(shared_residual_range);

    // convert what Fortran cannot take from R directly
    const int* _c_residuals_perm = tox_int_in(residuals_perm, "residuals_perm");

    // outputs and work space
    int n_bins = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    estimate_bin_count_expert_c(
        REAL(residuals),
        _c_residuals_perm,
        &n_residuals,
        &_v_max_n_reps_all_studies,
        &_v_n_neighbors,
        &_v_shared_residual_range,
        &n_bins,
        &sturges_bins,
        &fd_bins,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_bin_count_occupancy_call(SEXP pooled_residuals, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_reps_all_studies = Rf_asInteger(max_n_reps_all_studies);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);
    double _v_gamma_occupancy = Rf_asReal(gamma_occupancy);
    double _v_lower_residual_range_quantile = Rf_asReal(lower_residual_range_quantile);
    double _v_upper_residual_range_quantile = Rf_asReal(upper_residual_range_quantile);

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char occupancy_failed = 0;
    double shared_residual_range_low = 0;
    double shared_residual_range_high = 0;
    int n_pooled_residuals = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    determine_bin_count_occupancy_c(
        REAL(pooled_residuals),
        &n_residuals,
        &_v_max_n_reps_all_studies,
        &_v_n_neighbors,
        &selected_n_bins,
        &occupancy_failed,
        &shared_residual_range_low,
        &shared_residual_range_high,
        &n_pooled_residuals,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &sturges_bins,
        &fd_bins,
        &_v_m_min,
        &_v_m_max,
        &_v_min_residuals_per_bin,
        &_v_gamma_occupancy,
        &_v_lower_residual_range_quantile,
        &_v_upper_residual_range_quantile,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 11)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(occupancy_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(shared_residual_range_low));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(shared_residual_range_high));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(n_pooled_residuals));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 9, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 10, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 11)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_bin_count_occupancy_expert_call(SEXP pooled_residuals, SEXP pooled_residuals_perm, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin, SEXP gamma_occupancy, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_n_reps_all_studies = Rf_asInteger(max_n_reps_all_studies);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);
    double _v_gamma_occupancy = Rf_asReal(gamma_occupancy);
    double _v_lower_residual_range_quantile = Rf_asReal(lower_residual_range_quantile);
    double _v_upper_residual_range_quantile = Rf_asReal(upper_residual_range_quantile);

    // convert what Fortran cannot take from R directly
    const int* _c_pooled_residuals_perm = tox_int_in(pooled_residuals_perm, "pooled_residuals_perm");

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char occupancy_failed = 0;
    double shared_residual_range_low = 0;
    double shared_residual_range_high = 0;
    int n_pooled_residuals = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int* tmp_bin_counts = (int*) R_alloc(256, sizeof(int));
    int ierr = 0;

    determine_bin_count_occupancy_expert_c(
        REAL(pooled_residuals),
        _c_pooled_residuals_perm,
        &n_residuals,
        &_v_max_n_reps_all_studies,
        &_v_n_neighbors,
        &selected_n_bins,
        &occupancy_failed,
        &shared_residual_range_low,
        &shared_residual_range_high,
        &n_pooled_residuals,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &sturges_bins,
        &fd_bins,
        tmp_bin_counts,
        &_v_m_min,
        &_v_m_max,
        &_v_min_residuals_per_bin,
        &_v_gamma_occupancy,
        &_v_lower_residual_range_quantile,
        &_v_upper_residual_range_quantile,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 11)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(occupancy_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(shared_residual_range_low));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(shared_residual_range_high));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(n_pooled_residuals));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 9, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 10, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 11)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_bin_count_occupancy_exhaustive_call(SEXP pooled_residuals, SEXP n_pooled_residuals, SEXP shared_residual_range_low, SEXP shared_residual_range_high, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_pooled_residuals = Rf_asInteger(n_pooled_residuals);
    double _v_shared_residual_range_low = Rf_asReal(shared_residual_range_low);
    double _v_shared_residual_range_high = Rf_asReal(shared_residual_range_high);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char occupancy_failed = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int ierr = 0;

    determine_bin_count_occupancy_exhaustive_c(
        REAL(pooled_residuals),
        &n_residuals,
        &_v_n_pooled_residuals,
        &_v_shared_residual_range_low,
        &_v_shared_residual_range_high,
        &selected_n_bins,
        &occupancy_failed,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &_v_m_min,
        &_v_m_max,
        &_v_min_residuals_per_bin,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(occupancy_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_bin_count_occupancy_exhaustive_expert_call(SEXP pooled_residuals, SEXP pooled_residuals_perm, SEXP n_pooled_residuals, SEXP shared_residual_range_low, SEXP shared_residual_range_high, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_pooled_residuals = Rf_asInteger(n_pooled_residuals);
    double _v_shared_residual_range_low = Rf_asReal(shared_residual_range_low);
    double _v_shared_residual_range_high = Rf_asReal(shared_residual_range_high);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);

    // convert what Fortran cannot take from R directly
    const int* _c_pooled_residuals_perm = tox_int_in(pooled_residuals_perm, "pooled_residuals_perm");

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char occupancy_failed = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int ierr = 0;

    determine_bin_count_occupancy_exhaustive_expert_c(
        REAL(pooled_residuals),
        _c_pooled_residuals_perm,
        &n_residuals,
        &_v_n_pooled_residuals,
        &_v_shared_residual_range_low,
        &_v_shared_residual_range_high,
        &selected_n_bins,
        &occupancy_failed,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &_v_m_min,
        &_v_m_max,
        &_v_min_residuals_per_bin,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(occupancy_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("occupancy_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_bin_count_stone_call(SEXP pooled_residuals, SEXP residual_to_study, SEXP n_studies, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP bin_count_mode, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_studies = Rf_asInteger(n_studies);
    int _v_max_n_reps_all_studies = Rf_asInteger(max_n_reps_all_studies);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);
    double _v_lower_residual_range_quantile = Rf_asReal(lower_residual_range_quantile);
    double _v_upper_residual_range_quantile = Rf_asReal(upper_residual_range_quantile);

    // convert what Fortran cannot take from R directly
    const int* _c_residual_to_study = tox_int_in(residual_to_study, "residual_to_study");
    char* _c_bin_count_mode = tox_char_in(bin_count_mode, 29, "bin_count_mode");

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char stone_failed = 0;
    double shared_residual_range_low = 0;
    double shared_residual_range_high = 0;
    int n_pooled_residuals = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int ierr = 0;

    determine_bin_count_stone_c(
        REAL(pooled_residuals),
        &n_residuals,
        _c_residual_to_study,
        &_v_n_studies,
        &_v_max_n_reps_all_studies,
        &_v_n_neighbors,
        &selected_n_bins,
        &stone_failed,
        &shared_residual_range_low,
        &shared_residual_range_high,
        &n_pooled_residuals,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &sturges_bins,
        &fd_bins,
        _c_bin_count_mode,
        &_v_m_min,
        &_v_m_max,
        &_v_min_residuals_per_bin,
        &_v_lower_residual_range_quantile,
        &_v_upper_residual_range_quantile,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 11)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(stone_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(shared_residual_range_low));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(shared_residual_range_high));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(n_pooled_residuals));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 9, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 10, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 11)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("stone_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP determine_bin_count_stone_expert_call(SEXP pooled_residuals, SEXP pooled_residuals_perm, SEXP residual_to_study, SEXP n_studies, SEXP max_n_reps_all_studies, SEXP n_neighbors, SEXP bin_count_mode, SEXP m_min, SEXP m_max, SEXP min_residuals_per_bin, SEXP lower_residual_range_quantile, SEXP upper_residual_range_quantile) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_residuals = (int) Rf_length(pooled_residuals);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_studies = Rf_asInteger(n_studies);
    int _v_max_n_reps_all_studies = Rf_asInteger(max_n_reps_all_studies);
    int _v_n_neighbors = Rf_asInteger(n_neighbors);
    int _v_m_min = Rf_asInteger(m_min);
    int _v_m_max = Rf_asInteger(m_max);
    int _v_min_residuals_per_bin = Rf_asInteger(min_residuals_per_bin);
    double _v_lower_residual_range_quantile = Rf_asReal(lower_residual_range_quantile);
    double _v_upper_residual_range_quantile = Rf_asReal(upper_residual_range_quantile);

    // convert what Fortran cannot take from R directly
    const int* _c_pooled_residuals_perm = tox_int_in(pooled_residuals_perm, "pooled_residuals_perm");
    const int* _c_residual_to_study = tox_int_in(residual_to_study, "residual_to_study");
    char* _c_bin_count_mode = tox_char_in(bin_count_mode, 29, "bin_count_mode");

    // outputs and work space
    int selected_n_bins = 0;
    unsigned char stone_failed = 0;
    double shared_residual_range_low = 0;
    double shared_residual_range_high = 0;
    int n_pooled_residuals = 0;
    int min_bin_occupancy = 0;
    double mean_bin_occupancy = 0;
    int max_bin_occupancy = 0;
    int sturges_bins = 0;
    int fd_bins = 0;
    int* tmp_bin_counts = (int*) R_alloc(256, sizeof(int));
    int* tmp_residual_to_study_perm = (int*) R_alloc(n_residuals, sizeof(int));
    int* tmp_study_block_sizes = (int*) R_alloc(_v_n_studies, sizeof(int));
    int ierr = 0;

    determine_bin_count_stone_expert_c(
        REAL(pooled_residuals),
        _c_pooled_residuals_perm,
        &n_residuals,
        _c_residual_to_study,
        &_v_n_studies,
        &_v_max_n_reps_all_studies,
        &_v_n_neighbors,
        &selected_n_bins,
        &stone_failed,
        &shared_residual_range_low,
        &shared_residual_range_high,
        &n_pooled_residuals,
        &min_bin_occupancy,
        &mean_bin_occupancy,
        &max_bin_occupancy,
        &sturges_bins,
        &fd_bins,
        tmp_bin_counts,
        tmp_residual_to_study_perm,
        tmp_study_block_sizes,
        _c_bin_count_mode,
        &_v_m_min,
        &_v_m_max,
        &_v_min_residuals_per_bin,
        &_v_lower_residual_range_quantile,
        &_v_upper_residual_range_quantile,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 11)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(selected_n_bins));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarLogical(stone_failed != 0));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(shared_residual_range_low));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(shared_residual_range_high));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(n_pooled_residuals));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(min_bin_occupancy));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarReal(mean_bin_occupancy));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(max_bin_occupancy));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(sturges_bins));
    SET_VECTOR_ELT(_out, 9, Rf_ScalarInteger(fd_bins));
    SET_VECTOR_ELT(_out, 10, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 11)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("selected_n_bins"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("stone_failed"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("shared_residual_range_low"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("shared_residual_range_high"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_pooled_residuals"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("min_bin_occupancy"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("mean_bin_occupancy"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("max_bin_occupancy"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("sturges_bins"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("fd_bins"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
