// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void perform_permutation_test_expert_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, const int*, double*, double*, double*, double*, const int*, int*);
void perform_permutation_test_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, const int*, double*, double*, const int*, int*);
void compute_p_values_c(const double*, const double*, const double*, const double*, const int*, const int*, double*, double*, int*);
void compute_contributions_c(const double*, const double*, const int*, const char*, double*, double*, int*);
void compute_all_contributions_expert_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, int*);
void compute_all_contributions_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, double*, double*, int*);
void compute_baselines_factor_dependent_c(const int*, const double*, const double*, const char*, double*, double*, int*);
void compute_velocity_trajectory_c(const double*, double*, const int*, int*);
void compute_acceleration_from_velocity_trajectory_c(const double*, double*, const int*, int*);
void compute_velocity_trajectories_c(const double*, double*, const int*, const int*, const int*, int*);
void compute_acceleration_from_velocity_c(const double*, double*, const int*, const int*, const int*, int*);
void compute_velocity_acceleration_contributions_expert_c(const double*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, double*, double*, double*, int*);
void compute_velocity_acceleration_contributions_c(const double*, const int*, const int*, const int*, const char*, double*, double*, double*, double*, int*);

SEXP perform_permutation_test_expert_call(SEXP trajectories, SEXP factor_idx, SEXP dependent_idx, SEXP sample_idx, SEXP baseline_mode, SEXP n_permutations, SEXP random_seed) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = NULL;
    int random_seed_size = 0;
    if (random_seed != R_NilValue) {
        random_seed_size = (int) Rf_length(random_seed);
        random_seed_p = INTEGER(random_seed);
    }

    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int factor_idx_v = Rf_asInteger(factor_idx);
    int dependent_idx_v = Rf_asInteger(dependent_idx);
    int sample_idx_v = Rf_asInteger(sample_idx);
    int n_permutations_v = Rf_asInteger(n_permutations);

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_permutations_v)); nprot++;
    { SEXP local_contributions_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(local_contributions_dim)[0] = n_timepoints; INTEGER(local_contributions_dim)[1] = n_permutations_v; Rf_setAttrib(local_contributions, R_DimSymbol, local_contributions_dim); UNPROTECT(1); }
    SEXP total_contributions = PROTECT(Rf_allocVector(REALSXP, n_permutations_v)); nprot++;
    double* tmp_factor = (double*) R_alloc(n_timepoints, sizeof(double));
    double* tmp_dependent = (double*) R_alloc(n_timepoints, sizeof(double));
    int ierr = 0;

    perform_permutation_test_expert_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &factor_idx_v,
        &dependent_idx_v,
        &sample_idx_v,
        baseline_mode_c,
        &n_permutations_v,
        REAL(local_contributions),
        REAL(total_contributions),
        tmp_factor,
        tmp_dependent,
        random_seed_p,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, total_contributions);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contributions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP perform_permutation_test_call(SEXP trajectories, SEXP factor_idx, SEXP dependent_idx, SEXP sample_idx, SEXP baseline_mode, SEXP n_permutations, SEXP random_seed) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* random_seed_p = NULL;
    int random_seed_size = 0;
    if (random_seed != R_NilValue) {
        random_seed_size = (int) Rf_length(random_seed);
        random_seed_p = INTEGER(random_seed);
    }

    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int factor_idx_v = Rf_asInteger(factor_idx);
    int dependent_idx_v = Rf_asInteger(dependent_idx);
    int sample_idx_v = Rf_asInteger(sample_idx);
    int n_permutations_v = Rf_asInteger(n_permutations);

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_permutations_v)); nprot++;
    { SEXP local_contributions_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(local_contributions_dim)[0] = n_timepoints; INTEGER(local_contributions_dim)[1] = n_permutations_v; Rf_setAttrib(local_contributions, R_DimSymbol, local_contributions_dim); UNPROTECT(1); }
    SEXP total_contributions = PROTECT(Rf_allocVector(REALSXP, n_permutations_v)); nprot++;
    int ierr = 0;

    perform_permutation_test_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &factor_idx_v,
        &dependent_idx_v,
        &sample_idx_v,
        baseline_mode_c,
        &n_permutations_v,
        REAL(local_contributions),
        REAL(total_contributions),
        random_seed_p,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, total_contributions);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contributions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_p_values_call(SEXP local_contributions_observed, SEXP total_contribution_observed, SEXP local_contributions_perm_test, SEXP total_contributions_perm_test) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) Rf_length(local_contributions_observed);
    int n_permutations = INTEGER(Rf_getAttrib(local_contributions_perm_test, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double total_contribution_observed_v = Rf_asReal(total_contribution_observed);

    // outputs and work space
    SEXP local_p_values = PROTECT(Rf_allocVector(REALSXP, n_timepoints)); nprot++;
    double total_p_value = 0;
    int ierr = 0;

    compute_p_values_c(
        REAL(local_contributions_observed),
        &total_contribution_observed_v,
        REAL(local_contributions_perm_test),
        REAL(total_contributions_perm_test),
        &n_timepoints,
        &n_permutations,
        REAL(local_p_values),
        &total_p_value,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, local_p_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(total_p_value));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_p_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_p_value"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_contributions_call(SEXP factor, SEXP dependent, SEXP baseline_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = (int) Rf_length(factor);

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_dims)); nprot++;
    double total_contribution = 0;
    int ierr = 0;

    compute_contributions_c(
        REAL(factor),
        REAL(dependent),
        &n_dims,
        baseline_mode_c,
        REAL(local_contributions),
        &total_contribution,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(total_contribution));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contribution"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_all_contributions_expert_call(SEXP trajectories, SEXP factor_indices, SEXP dependent_indices, SEXP baseline_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];
    int n_selected_factors = (int) Rf_length(factor_indices);
    int n_selected_dependents = (int) Rf_length(dependent_indices);

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_selected_factors * n_selected_dependents * n_samples)); nprot++;
    { SEXP local_contributions_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(local_contributions_dim)[0] = n_timepoints; INTEGER(local_contributions_dim)[1] = n_selected_factors; INTEGER(local_contributions_dim)[2] = n_selected_dependents; INTEGER(local_contributions_dim)[3] = n_samples; Rf_setAttrib(local_contributions, R_DimSymbol, local_contributions_dim); UNPROTECT(1); }
    SEXP total_contributions = PROTECT(Rf_allocVector(REALSXP, n_selected_factors * n_selected_dependents * n_samples)); nprot++;
    { SEXP total_contributions_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(total_contributions_dim)[0] = n_selected_factors; INTEGER(total_contributions_dim)[1] = n_selected_dependents; INTEGER(total_contributions_dim)[2] = n_samples; Rf_setAttrib(total_contributions, R_DimSymbol, total_contributions_dim); UNPROTECT(1); }
    double* tmp_factors = (double*) R_alloc(n_timepoints * n_selected_factors, sizeof(double));
    double* tmp_dependent = (double*) R_alloc(n_timepoints, sizeof(double));
    int ierr = 0;

    compute_all_contributions_expert_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        INTEGER(factor_indices),
        &n_selected_factors,
        INTEGER(dependent_indices),
        &n_selected_dependents,
        baseline_mode_c,
        REAL(local_contributions),
        REAL(total_contributions),
        tmp_factors,
        tmp_dependent,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, total_contributions);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contributions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_all_contributions_call(SEXP trajectories, SEXP factor_indices, SEXP dependent_indices, SEXP baseline_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];
    int n_selected_factors = (int) Rf_length(factor_indices);
    int n_selected_dependents = (int) Rf_length(dependent_indices);

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    SEXP local_contributions = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_selected_factors * n_selected_dependents * n_samples)); nprot++;
    { SEXP local_contributions_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(local_contributions_dim)[0] = n_timepoints; INTEGER(local_contributions_dim)[1] = n_selected_factors; INTEGER(local_contributions_dim)[2] = n_selected_dependents; INTEGER(local_contributions_dim)[3] = n_samples; Rf_setAttrib(local_contributions, R_DimSymbol, local_contributions_dim); UNPROTECT(1); }
    SEXP total_contributions = PROTECT(Rf_allocVector(REALSXP, n_selected_factors * n_selected_dependents * n_samples)); nprot++;
    { SEXP total_contributions_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(total_contributions_dim)[0] = n_selected_factors; INTEGER(total_contributions_dim)[1] = n_selected_dependents; INTEGER(total_contributions_dim)[2] = n_samples; Rf_setAttrib(total_contributions, R_DimSymbol, total_contributions_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_all_contributions_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        INTEGER(factor_indices),
        &n_selected_factors,
        INTEGER(dependent_indices),
        &n_selected_dependents,
        baseline_mode_c,
        REAL(local_contributions),
        REAL(total_contributions),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, local_contributions);
    SET_VECTOR_ELT(_out, 1, total_contributions);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("local_contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("total_contributions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_baselines_factor_dependent_call(SEXP factor, SEXP dependent, SEXP baseline_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) Rf_length(factor);

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    double factor_baseline = 0;
    double dependent_baseline = 0;
    int ierr = 0;

    compute_baselines_factor_dependent_c(
        &n_timepoints,
        REAL(factor),
        REAL(dependent),
        baseline_mode_c,
        &factor_baseline,
        &dependent_baseline,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(factor_baseline));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(dependent_baseline));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("factor_baseline"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("dependent_baseline"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_velocity_trajectory_call(SEXP trajectory) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_timepoints = (int) Rf_length(trajectory);

    // outputs and work space
    SEXP velocity = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, n_timepoints - 1)))); nprot++;
    int ierr = 0;

    compute_velocity_trajectory_c(
        REAL(trajectory),
        REAL(velocity),
        &n_timepoints,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, velocity);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_acceleration_from_velocity_trajectory_call(SEXP velocity, SEXP n_timepoints) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_timepoints_v = Rf_asInteger(n_timepoints);

    // outputs and work space
    SEXP acceleration = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, n_timepoints_v - 2)))); nprot++;
    int ierr = 0;

    compute_acceleration_from_velocity_trajectory_c(
        REAL(velocity),
        REAL(acceleration),
        &n_timepoints_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, acceleration);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("acceleration"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_velocity_trajectories_call(SEXP trajectories) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // outputs and work space
    SEXP velocity = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, n_timepoints - 1)) * n_factors * n_samples)); nprot++;
    { SEXP velocity_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(velocity_dim)[0] = tox_imax(0, n_timepoints - 1); INTEGER(velocity_dim)[1] = n_factors; INTEGER(velocity_dim)[2] = n_samples; Rf_setAttrib(velocity, R_DimSymbol, velocity_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_velocity_trajectories_c(
        REAL(trajectories),
        REAL(velocity),
        &n_factors,
        &n_samples,
        &n_timepoints,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, velocity);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_acceleration_from_velocity_call(SEXP velocity, SEXP n_timepoints) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(velocity, R_DimSymbol))[1];
    int n_samples = INTEGER(Rf_getAttrib(velocity, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int n_timepoints_v = Rf_asInteger(n_timepoints);

    // outputs and work space
    SEXP acceleration = PROTECT(Rf_allocVector(REALSXP, (tox_imax(0, n_timepoints_v - 2)) * n_factors * n_samples)); nprot++;
    { SEXP acceleration_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(acceleration_dim)[0] = tox_imax(0, n_timepoints_v - 2); INTEGER(acceleration_dim)[1] = n_factors; INTEGER(acceleration_dim)[2] = n_samples; Rf_setAttrib(acceleration, R_DimSymbol, acceleration_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_acceleration_from_velocity_c(
        REAL(velocity),
        REAL(acceleration),
        &n_factors,
        &n_samples,
        &n_timepoints_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, acceleration);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("acceleration"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_velocity_acceleration_contributions_expert_call(SEXP trajectories, SEXP baseline_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    double* tmp_factors = (double*) R_alloc((n_timepoints - 1) * n_factors, sizeof(double));
    double* tmp_dependent = (double*) R_alloc((n_timepoints - 1), sizeof(double));
    double* tmp_contributions = (double*) R_alloc((n_timepoints - 1), sizeof(double));
    SEXP contrib_velocity = PROTECT(Rf_allocVector(REALSXP, n_factors * n_factors * n_samples)); nprot++;
    { SEXP contrib_velocity_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(contrib_velocity_dim)[0] = n_factors; INTEGER(contrib_velocity_dim)[1] = n_factors; INTEGER(contrib_velocity_dim)[2] = n_samples; Rf_setAttrib(contrib_velocity, R_DimSymbol, contrib_velocity_dim); UNPROTECT(1); }
    SEXP velocity_contribution_series = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors * n_factors * n_samples)); nprot++;
    { SEXP velocity_contribution_series_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(velocity_contribution_series_dim)[0] = n_timepoints; INTEGER(velocity_contribution_series_dim)[1] = n_factors; INTEGER(velocity_contribution_series_dim)[2] = n_factors; INTEGER(velocity_contribution_series_dim)[3] = n_samples; Rf_setAttrib(velocity_contribution_series, R_DimSymbol, velocity_contribution_series_dim); UNPROTECT(1); }
    SEXP contrib_acceleration = PROTECT(Rf_allocVector(REALSXP, n_factors * n_factors * n_samples)); nprot++;
    { SEXP contrib_acceleration_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(contrib_acceleration_dim)[0] = n_factors; INTEGER(contrib_acceleration_dim)[1] = n_factors; INTEGER(contrib_acceleration_dim)[2] = n_samples; Rf_setAttrib(contrib_acceleration, R_DimSymbol, contrib_acceleration_dim); UNPROTECT(1); }
    SEXP acceleration_contribution_series = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors * n_factors * n_samples)); nprot++;
    { SEXP acceleration_contribution_series_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(acceleration_contribution_series_dim)[0] = n_timepoints; INTEGER(acceleration_contribution_series_dim)[1] = n_factors; INTEGER(acceleration_contribution_series_dim)[2] = n_factors; INTEGER(acceleration_contribution_series_dim)[3] = n_samples; Rf_setAttrib(acceleration_contribution_series, R_DimSymbol, acceleration_contribution_series_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_velocity_acceleration_contributions_expert_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        baseline_mode_c,
        tmp_factors,
        tmp_dependent,
        tmp_contributions,
        REAL(contrib_velocity),
        REAL(velocity_contribution_series),
        REAL(contrib_acceleration),
        REAL(acceleration_contribution_series),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, contrib_velocity);
    SET_VECTOR_ELT(_out, 1, velocity_contribution_series);
    SET_VECTOR_ELT(_out, 2, contrib_acceleration);
    SET_VECTOR_ELT(_out, 3, acceleration_contribution_series);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("contrib_velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("velocity_contribution_series"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("contrib_acceleration"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("acceleration_contribution_series"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_velocity_acceleration_contributions_call(SEXP trajectories, SEXP baseline_mode) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // convert what Fortran cannot take from R directly
    char* baseline_mode_c = tox_char_in(baseline_mode, 4);

    // outputs and work space
    SEXP contrib_velocity = PROTECT(Rf_allocVector(REALSXP, n_factors * n_factors * n_samples)); nprot++;
    { SEXP contrib_velocity_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(contrib_velocity_dim)[0] = n_factors; INTEGER(contrib_velocity_dim)[1] = n_factors; INTEGER(contrib_velocity_dim)[2] = n_samples; Rf_setAttrib(contrib_velocity, R_DimSymbol, contrib_velocity_dim); UNPROTECT(1); }
    SEXP velocity_contribution_series = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors * n_factors * n_samples)); nprot++;
    { SEXP velocity_contribution_series_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(velocity_contribution_series_dim)[0] = n_timepoints; INTEGER(velocity_contribution_series_dim)[1] = n_factors; INTEGER(velocity_contribution_series_dim)[2] = n_factors; INTEGER(velocity_contribution_series_dim)[3] = n_samples; Rf_setAttrib(velocity_contribution_series, R_DimSymbol, velocity_contribution_series_dim); UNPROTECT(1); }
    SEXP contrib_acceleration = PROTECT(Rf_allocVector(REALSXP, n_factors * n_factors * n_samples)); nprot++;
    { SEXP contrib_acceleration_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(contrib_acceleration_dim)[0] = n_factors; INTEGER(contrib_acceleration_dim)[1] = n_factors; INTEGER(contrib_acceleration_dim)[2] = n_samples; Rf_setAttrib(contrib_acceleration, R_DimSymbol, contrib_acceleration_dim); UNPROTECT(1); }
    SEXP acceleration_contribution_series = PROTECT(Rf_allocVector(REALSXP, n_timepoints * n_factors * n_factors * n_samples)); nprot++;
    { SEXP acceleration_contribution_series_dim = PROTECT(Rf_allocVector(INTSXP, 4)); INTEGER(acceleration_contribution_series_dim)[0] = n_timepoints; INTEGER(acceleration_contribution_series_dim)[1] = n_factors; INTEGER(acceleration_contribution_series_dim)[2] = n_factors; INTEGER(acceleration_contribution_series_dim)[3] = n_samples; Rf_setAttrib(acceleration_contribution_series, R_DimSymbol, acceleration_contribution_series_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_velocity_acceleration_contributions_c(
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        baseline_mode_c,
        REAL(contrib_velocity),
        REAL(velocity_contribution_series),
        REAL(contrib_acceleration),
        REAL(acceleration_contribution_series),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, contrib_velocity);
    SET_VECTOR_ELT(_out, 1, velocity_contribution_series);
    SET_VECTOR_ELT(_out, 2, contrib_acceleration);
    SET_VECTOR_ELT(_out, 3, acceleration_contribution_series);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("contrib_velocity"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("velocity_contribution_series"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("contrib_acceleration"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("acceleration_contribution_series"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
