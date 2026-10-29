// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void perform_permutation_test_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, const int*, double*, double*, const int*, int*);
void compute_p_values_c(const double*, const double*, const double*, const double*, const int*, const int*, double*, double*, int*);
void compute_contributions_c(const double*, const double*, const int*, const char*, double*, double*, int*);
void compute_all_contributions_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, double*, double*, int*);
void compute_baselines_factor_dependent_c(const int*, const double*, const double*, const char*, double*, double*, int*);
void compute_velocity_trajectory_c(const double*, double*, const int*, int*);
void compute_acceleration_from_velocity_trajectory_c(const double*, const int*, double*, const int*, int*);
void compute_velocity_trajectories_c(const double*, double*, const int*, const int*, const int*, int*);
void compute_acceleration_from_velocity_c(const double*, const int*, double*, const int*, const int*, const int*, int*);
void compute_velocity_acceleration_contributions_c(const double*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, int*);

SEXP perform_permutation_test_call(SEXP trajectories, SEXP factor_idx, SEXP dependent_idx, SEXP sample_idx, SEXP baseline_mode, SEXP n_permutations, SEXP random_seed) {
    int _nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* _p_random_seed = NULL;
    int _size_random_seed = 0;
    if (random_seed != R_NilValue) {
        _size_random_seed = (int) Rf_length(random_seed);
        _p_random_seed = INTEGER(random_seed);
    }

    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int _v_factor_idx = Rf_asInteger(factor_idx);
    int _v_dependent_idx = Rf_asInteger(dependent_idx);
    int _v_sample_idx = Rf_asInteger(sample_idx);
    int _v_n_permutations = Rf_asInteger(n_permutations);

    // convert what Fortran cannot take from R directly
    char* _c_baseline_mode = tox_char_in(baseline_mode, 13, "baseline_mode");

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_timepoints * _v_n_permutations)); _nprot++;
    { SEXP _dim_local_contributions = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_local_contributions)[0] = n_timepoints; INTEGER(_dim_local_contributions)[1] = _v_n_permutations; Rf_setAttrib(local_contributions, R_DimSymbol, _dim_local_contributions); UNPROTECT(1); }
    SEXP total_contributions = PROTECT(Rf_allocVector(REALSXP, _v_n_permutations)); _nprot++;
    int ierr = 0;

    perform_permutation_test_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &_v_factor_idx,
        &_v_dependent_idx,
        &_v_sample_idx,
        _c_baseline_mode,
        &_v_n_permutations,
        REAL(local_contributions),
        REAL(total_contributions),
        _p_random_seed,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, total_contributions);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contributions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_p_values_call(SEXP local_contributions_observed, SEXP total_contribution_observed, SEXP local_contributions_perm_test, SEXP total_contributions_perm_test) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) Rf_length(local_contributions_observed);
    int n_permutations = INTEGER(Rf_getAttrib(local_contributions_perm_test, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_total_contribution_observed = Rf_asReal(total_contribution_observed);

    // outputs and work space
    SEXP local_p_values = PROTECT(Rf_allocVector(REALSXP, n_timepoints)); _nprot++;
    double total_p_value = 0;
    int ierr = 0;

    compute_p_values_c(
        REAL(local_contributions_observed),
        &_v_total_contribution_observed,
        REAL(local_contributions_perm_test),
        REAL(total_contributions_perm_test),
        &n_timepoints,
        &n_permutations,
        REAL(local_p_values),
        &total_p_value,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, local_p_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(total_p_value));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_p_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_p_value"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_contributions_call(SEXP factor, SEXP dependent, SEXP baseline_mode) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = (int) Rf_length(factor);

    // convert what Fortran cannot take from R directly
    char* _c_baseline_mode = tox_char_in(baseline_mode, 13, "baseline_mode");

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_dims)); _nprot++;
    double total_contribution = 0;
    int ierr = 0;

    compute_contributions_c(
        REAL(factor),
        REAL(dependent),
        &n_dims,
        _c_baseline_mode,
        REAL(local_contributions),
        &total_contribution,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(total_contribution));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contribution"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_all_contributions_call(SEXP trajectories, SEXP factor_indices, SEXP dependent_indices, SEXP baseline_mode) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];
    int n_selected_factors = (int) Rf_length(factor_indices);
    int n_selected_dependents = (int) Rf_length(dependent_indices);

    // convert what Fortran cannot take from R directly
    const int* _c_factor_indices = tox_int_in(factor_indices, "factor_indices");
    const int* _c_dependent_indices = tox_int_in(dependent_indices, "dependent_indices");
    char* _c_baseline_mode = tox_char_in(baseline_mode, 13, "baseline_mode");

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_selected_factors * n_selected_dependents * n_samples)); _nprot++;
    { SEXP _dim_local_contributions = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(_dim_local_contributions)[0] = n_timepoints; INTEGER(_dim_local_contributions)[1] = n_selected_factors; INTEGER(_dim_local_contributions)[2] = n_selected_dependents; INTEGER(_dim_local_contributions)[3] = n_samples; Rf_setAttrib(local_contributions, R_DimSymbol, _dim_local_contributions); UNPROTECT(1); }
    SEXP total_contributions = PROTECT(Rf_allocVector(REALSXP, n_selected_factors * n_selected_dependents * n_samples)); _nprot++;
    { SEXP _dim_total_contributions = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_total_contributions)[0] = n_selected_factors; INTEGER(_dim_total_contributions)[1] = n_selected_dependents; INTEGER(_dim_total_contributions)[2] = n_samples; Rf_setAttrib(total_contributions, R_DimSymbol, _dim_total_contributions); UNPROTECT(1); }
    int ierr = 0;

    compute_all_contributions_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        _c_factor_indices,
        &n_selected_factors,
        _c_dependent_indices,
        &n_selected_dependents,
        _c_baseline_mode,
        REAL(local_contributions),
        REAL(total_contributions),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, total_contributions);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contributions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_baselines_factor_dependent_call(SEXP factor, SEXP dependent, SEXP baseline_mode) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) Rf_length(factor);

    // convert what Fortran cannot take from R directly
    char* _c_baseline_mode = tox_char_in(baseline_mode, 13, "baseline_mode");

    // outputs and work space
    double factor_baseline = 0;
    double dependent_baseline = 0;
    int ierr = 0;

    compute_baselines_factor_dependent_c(
        &n_timepoints,
        REAL(factor),
        REAL(dependent),
        _c_baseline_mode,
        &factor_baseline,
        &dependent_baseline,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(factor_baseline));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(dependent_baseline));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("factor_baseline"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("dependent_baseline"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_velocity_trajectory_call(SEXP trajectory) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) Rf_length(trajectory);

    // outputs and work space
    SEXP velocity = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, n_timepoints - 1)))); _nprot++;
    int ierr = 0;

    compute_velocity_trajectory_c(
        REAL(trajectory),
        REAL(velocity),
        &n_timepoints,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, velocity);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_acceleration_from_velocity_trajectory_call(SEXP velocity, SEXP n_timepoints) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_velocity_elements = (int) Rf_length(velocity);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_timepoints = Rf_asInteger(n_timepoints);

    // outputs and work space
    SEXP acceleration = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, _v_n_timepoints - 2)))); _nprot++;
    int ierr = 0;

    compute_acceleration_from_velocity_trajectory_c(
        REAL(velocity),
        &n_velocity_elements,
        REAL(acceleration),
        &_v_n_timepoints,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, acceleration);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("acceleration"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_velocity_trajectories_call(SEXP trajectories) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // outputs and work space
    SEXP velocity = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, n_timepoints - 1)) * n_factors * n_samples)); _nprot++;
    { SEXP _dim_velocity = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_velocity)[0] = tox_imax(0, n_timepoints - 1); INTEGER(_dim_velocity)[1] = n_factors; INTEGER(_dim_velocity)[2] = n_samples; Rf_setAttrib(velocity, R_DimSymbol, _dim_velocity); UNPROTECT(1); }
    int ierr = 0;

    compute_velocity_trajectories_c(
        REAL(trajectories),
        REAL(velocity),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, velocity);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_acceleration_from_velocity_call(SEXP velocity, SEXP n_timepoints) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_velocity_elements_dim_1 = INTEGER(Rf_getAttrib(velocity, R_DimSymbol))[0];
    int n_factors = INTEGER(Rf_getAttrib(velocity, R_DimSymbol))[1];
    int n_samples = INTEGER(Rf_getAttrib(velocity, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_timepoints = Rf_asInteger(n_timepoints);

    // outputs and work space
    SEXP acceleration = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, _v_n_timepoints - 2)) * n_factors * n_samples)); _nprot++;
    { SEXP _dim_acceleration = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_acceleration)[0] = tox_imax(0, _v_n_timepoints - 2); INTEGER(_dim_acceleration)[1] = n_factors; INTEGER(_dim_acceleration)[2] = n_samples; Rf_setAttrib(acceleration, R_DimSymbol, _dim_acceleration); UNPROTECT(1); }
    int ierr = 0;

    compute_acceleration_from_velocity_c(
        REAL(velocity),
        &n_velocity_elements_dim_1,
        REAL(acceleration),
        &n_factors,
        &n_samples,
        &_v_n_timepoints,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, acceleration);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("acceleration"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_velocity_acceleration_contributions_call(SEXP trajectories, SEXP baseline_mode) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // convert what Fortran cannot take from R directly
    char* _c_baseline_mode = tox_char_in(baseline_mode, 13, "baseline_mode");

    // outputs and work space
    SEXP contrib_velocity = PROTECT(Rf_allocVector(REALSXP, n_factors * n_factors * n_samples)); _nprot++;
    { SEXP _dim_contrib_velocity = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_contrib_velocity)[0] = n_factors; INTEGER(_dim_contrib_velocity)[1] = n_factors; INTEGER(_dim_contrib_velocity)[2] = n_samples; Rf_setAttrib(contrib_velocity, R_DimSymbol, _dim_contrib_velocity); UNPROTECT(1); }
    SEXP velocity_contribution_series = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors * n_factors * n_samples)); _nprot++;
    { SEXP _dim_velocity_contribution_series = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(_dim_velocity_contribution_series)[0] = n_timepoints; INTEGER(_dim_velocity_contribution_series)[1] = n_factors; INTEGER(_dim_velocity_contribution_series)[2] = n_factors; INTEGER(_dim_velocity_contribution_series)[3] = n_samples; Rf_setAttrib(velocity_contribution_series, R_DimSymbol, _dim_velocity_contribution_series); UNPROTECT(1); }
    SEXP contrib_acceleration = PROTECT(Rf_allocVector(REALSXP, n_factors * n_factors * n_samples)); _nprot++;
    { SEXP _dim_contrib_acceleration = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(_dim_contrib_acceleration)[0] = n_factors; INTEGER(_dim_contrib_acceleration)[1] = n_factors; INTEGER(_dim_contrib_acceleration)[2] = n_samples; Rf_setAttrib(contrib_acceleration, R_DimSymbol, _dim_contrib_acceleration); UNPROTECT(1); }
    SEXP acceleration_contribution_series = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors * n_factors * n_samples)); _nprot++;
    { SEXP _dim_acceleration_contribution_series = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(_dim_acceleration_contribution_series)[0] = n_timepoints; INTEGER(_dim_acceleration_contribution_series)[1] = n_factors; INTEGER(_dim_acceleration_contribution_series)[2] = n_factors; INTEGER(_dim_acceleration_contribution_series)[3] = n_samples; Rf_setAttrib(acceleration_contribution_series, R_DimSymbol, _dim_acceleration_contribution_series); UNPROTECT(1); }
    int ierr = 0;

    compute_velocity_acceleration_contributions_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        _c_baseline_mode,
        REAL(contrib_velocity),
        REAL(velocity_contribution_series),
        REAL(contrib_acceleration),
        REAL(acceleration_contribution_series),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); _nprot++;
    SET_VECTOR_ELT(_out, 0, contrib_velocity);
    SET_VECTOR_ELT(_out, 1, velocity_contribution_series);
    SET_VECTOR_ELT(_out, 2, contrib_acceleration);
    SET_VECTOR_ELT(_out, 3, acceleration_contribution_series);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("contrib_velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("velocity_contribution_series"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("contrib_acceleration"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("acceleration_contribution_series"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
