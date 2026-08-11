// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void ensemble_reconciliation_c(const unsigned char*, const int*, const int*, const char*, const double*, const unsigned char*, const int*, int*, int*, double*, int*);

SEXP ensemble_reconciliation_call(SEXP ensemble_masks, SEXP mode, SEXP min_overlap_coefficient, SEXP report_overlap_coefficient, SEXP max_group_size) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_vectors = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[0];
    int n_ensembles = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double min_overlap_coefficient_v = Rf_asReal(min_overlap_coefficient);
    unsigned char report_overlap_coefficient_v = (Rf_asLogical(report_overlap_coefficient) == TRUE) ? 1 : 0;
    int max_group_size_v = Rf_asInteger(max_group_size);

    // convert what Fortran cannot take from R directly
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);
    char* mode_c = tox_char_in(mode, 25);

    // outputs and work space
    SEXP super_ensembles = PROTECT(Rf_allocVector(INTSXP, max_group_size_v * (n_ensembles*(n_ensembles-1)))); nprot++;
    { SEXP super_ensembles_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(super_ensembles_dim)[0] = max_group_size_v; INTEGER(super_ensembles_dim)[1] = n_ensembles*(n_ensembles-1); Rf_setAttrib(super_ensembles, R_DimSymbol, super_ensembles_dim); UNPROTECT(1); }
    int n_super_ensembles = 0;
    SEXP super_ensembles_overlap_coefficient = PROTECT(Rf_allocVector(REALSXP, (max_group_size_v-1) * (n_ensembles*(n_ensembles-1)))); nprot++;
    { SEXP super_ensembles_overlap_coefficient_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(super_ensembles_overlap_coefficient_dim)[0] = max_group_size_v-1; INTEGER(super_ensembles_overlap_coefficient_dim)[1] = n_ensembles*(n_ensembles-1); Rf_setAttrib(super_ensembles_overlap_coefficient, R_DimSymbol, super_ensembles_overlap_coefficient_dim); UNPROTECT(1); }
    int ierr = 0;

    ensemble_reconciliation_c(
        ensemble_masks_c,
        &n_vectors,
        &n_ensembles,
        mode_c,
        &min_overlap_coefficient_v,
        &report_overlap_coefficient_v,
        &max_group_size_v,
        INTEGER(super_ensembles),
        &n_super_ensembles,
        REAL(super_ensembles_overlap_coefficient),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, super_ensembles);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_super_ensembles));
    SET_VECTOR_ELT(_out, 2, super_ensembles_overlap_coefficient);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("super_ensembles"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_super_ensembles"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("super_ensembles_overlap_coefficient"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
