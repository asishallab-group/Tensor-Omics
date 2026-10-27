// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void normalize_variable_timeseries_c(const double*, double*, const int*, int*, int*);
void normalize_single_trajectory_c(const double*, double*, const int*, const int*, int*, int*);
void normalize_all_trajectories_c(const double*, double*, const int*, const int*, const int*, int*, int*);

SEXP normalize_variable_timeseries_call(SEXP v) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = (int) Rf_length(v);

    // outputs and work space
    SEXP v_norm = PROTECT(Rf_allocVector(REALSXP, n_points)); _nprot++;
    int status = 0;
    int ierr = 0;

    normalize_variable_timeseries_c(
        REAL(v),
        REAL(v_norm),
        &n_points,
        &status,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, v_norm);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(status));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("v_norm"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP normalize_single_trajectory_call(SEXP trajectory) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectory, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectory, R_DimSymbol))[0];

    // outputs and work space
    SEXP trajectory_norm = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors)); _nprot++;
    { SEXP _dim_trajectory_norm = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_trajectory_norm)[0] = n_timepoints; INTEGER(_dim_trajectory_norm)[1] = n_factors; Rf_setAttrib(trajectory_norm, R_DimSymbol, _dim_trajectory_norm); UNPROTECT(1); }
    SEXP status = PROTECT(Rf_allocVector(INTSXP, n_factors)); _nprot++;
    int ierr = 0;

    normalize_single_trajectory_c(
        REAL(trajectory),
        REAL(trajectory_norm),
        &n_factors,
        &n_timepoints,
        INTEGER(status),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, trajectory_norm);
    SET_VECTOR_ELT(_out, 1, status);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("trajectory_norm"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP normalize_all_trajectories_call(SEXP trajectories) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // outputs and work space
    SEXP trajectories_norm = PROTECT(Rf_allocVector(REALSXP, n_factors * n_samples * n_timepoints)); _nprot++;
    { SEXP _dim_trajectories_norm = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_trajectories_norm)[0] = n_factors; INTEGER(_dim_trajectories_norm)[1] = n_samples; INTEGER(_dim_trajectories_norm)[2] = n_timepoints; Rf_setAttrib(trajectories_norm, R_DimSymbol, _dim_trajectories_norm); UNPROTECT(1); }
    SEXP status = PROTECT(Rf_allocVector(INTSXP, n_factors * n_samples)); _nprot++;
    { SEXP _dim_status = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_status)[0] = n_factors; INTEGER(_dim_status)[1] = n_samples; Rf_setAttrib(status, R_DimSymbol, _dim_status); UNPROTECT(1); }
    int ierr = 0;

    normalize_all_trajectories_c(
        REAL(trajectories),
        REAL(trajectories_norm),
        &n_factors,
        &n_samples,
        &n_timepoints,
        INTEGER(status),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, trajectories_norm);
    SET_VECTOR_ELT(_out, 1, status);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("trajectories_norm"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
