// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void gjct_permutation_test_c(const int*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, double*, double*, const int*, int*);
void bootstrap_histogram_c(const int*, const int*, const int*, const int*, const int*, const int*, const int*, double*, const double*, const int*, int*);

SEXP gjct_permutation_test_call(SEXP n_permutations, SEXP mean_pmf_counts, SEXP mean_pmf, SEXP mean_pmf_included_n_reps, SEXP included_n_reps, SEXP global_jsd_observed, SEXP random_seed) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];
    int n_studies = INTEGER(Rf_getAttrib(included_n_reps, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_permutations = Rf_asInteger(n_permutations);
    int _v_random_seed = Rf_asInteger(random_seed);

    // convert what Fortran cannot take from R directly
    const int* _c_mean_pmf_counts = tox_int_in(mean_pmf_counts, "mean_pmf_counts");
    const int* _c_mean_pmf_included_n_reps = tox_int_in(mean_pmf_included_n_reps, "mean_pmf_included_n_reps");
    const int* _c_included_n_reps = tox_int_in(included_n_reps, "included_n_reps");

    // outputs and work space
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    SEXP p_values_observed_consensus = PROTECT(Rf_allocVector(REALSXP, n_studies)); _nprot++;
    int ierr = 0;

    gjct_permutation_test_c(
        &_v_n_permutations,
        &n_bins,
        &n_points,
        &n_studies,
        _c_mean_pmf_counts,
        REAL(mean_pmf),
        _c_mean_pmf_included_n_reps,
        _c_included_n_reps,
        REAL(global_jsd_observed),
        REAL(p_values),
        REAL(p_values_observed_consensus),
        &_v_random_seed,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, p_values);
    SET_VECTOR_ELT(_out, 1, p_values_observed_consensus);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("p_values_observed_consensus"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bootstrap_histogram_call(SEXP n_bootstraps, SEXP mean_pmf_counts, SEXP mean_pmf_included_n_reps, SEXP included_n_reps, SEXP confidence_interval, SEXP two_sided_bootstrapping_significance_level, SEXP random_seed) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];
    int n_studies = INTEGER(Rf_getAttrib(included_n_reps, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bootstraps = Rf_asInteger(n_bootstraps);
    double _v_two_sided_bootstrapping_significance_level = Rf_asReal(two_sided_bootstrapping_significance_level);
    int _v_random_seed = Rf_asInteger(random_seed);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_confidence_interval = PROTECT(Rf_duplicate(confidence_interval)); _nprot++;

    // convert what Fortran cannot take from R directly
    const int* _c_mean_pmf_counts = tox_int_in(mean_pmf_counts, "mean_pmf_counts");
    const int* _c_mean_pmf_included_n_reps = tox_int_in(mean_pmf_included_n_reps, "mean_pmf_included_n_reps");
    const int* _c_included_n_reps = tox_int_in(included_n_reps, "included_n_reps");

    // outputs and work space
    int ierr = 0;

    bootstrap_histogram_c(
        &_v_n_bootstraps,
        &n_bins,
        &n_points,
        &n_studies,
        _c_mean_pmf_counts,
        _c_mean_pmf_included_n_reps,
        _c_included_n_reps,
        REAL(_out_confidence_interval),
        &_v_two_sided_bootstrapping_significance_level,
        &_v_random_seed,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_confidence_interval);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("confidence_interval"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
