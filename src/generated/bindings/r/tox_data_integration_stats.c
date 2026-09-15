// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void gjct_permutation_test_c(const int*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, double*, const int*, int*);

SEXP gjct_permutation_test_call(SEXP n_permutations, SEXP mean_pmf_counts, SEXP mean_pmf, SEXP mean_pmf_included_n_reps, SEXP included_n_reps, SEXP global_jsd_observed, SEXP random_seed) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bins = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(mean_pmf_counts, R_DimSymbol))[1];
    int n_studies = INTEGER(Rf_getAttrib(included_n_reps, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_permutations_v = Rf_asInteger(n_permutations);
    int random_seed_v = Rf_asInteger(random_seed);

    // outputs and work space
    SEXP p_values = PROTECT(Rf_allocVector(REALSXP, n_studies)); nprot++;
    int ierr = 0;

    gjct_permutation_test_c(
        &n_permutations_v,
        &n_bins,
        &n_points,
        &n_studies,
        INTEGER(mean_pmf_counts),
        REAL(mean_pmf),
        INTEGER(mean_pmf_included_n_reps),
        INTEGER(included_n_reps),
        REAL(global_jsd_observed),
        REAL(p_values),
        &random_seed_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, p_values);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("p_values"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
