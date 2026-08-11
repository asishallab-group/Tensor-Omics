// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void serialize_stc_points_as_csv_c(const char*, const int*, const int*, const int*, const int*, const int*, const unsigned char*, const unsigned char*, const unsigned char*, const int*, int*);
void serialize_stc_ensemble_overlap_as_csv_c(const char*, const int*, const int*, const int*, const unsigned char*, int*);
void serialize_stc_super_ensembles_as_tsv_c(const char*, const int*, const int*, const int*, const int*, int*);

SEXP serialize_stc_points_as_csv_call(SEXP filename, SEXP n_super_ensembles, SEXP seed_selection_mask, SEXP ensemble_masks, SEXP ensemble_low_confidence_masks, SEXP super_ensembles) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int n_vectors = (int) Rf_length(seed_selection_mask);
    int n_selected_seed = tox_sum_true(seed_selection_mask);
    int max_group_size = INTEGER(Rf_getAttrib(super_ensembles, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int n_super_ensembles_v = Rf_asInteger(n_super_ensembles);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);
    unsigned char* seed_selection_mask_c = tox_bool_in(seed_selection_mask);
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);
    unsigned char* ensemble_low_confidence_masks_c = tox_bool_in(ensemble_low_confidence_masks);

    // outputs and work space
    int ierr = 0;

    serialize_stc_points_as_csv_c(
        filename_c,
        &filename_strlen,
        &n_vectors,
        &n_selected_seed,
        &max_group_size,
        &n_super_ensembles_v,
        seed_selection_mask_c,
        ensemble_masks_c,
        ensemble_low_confidence_masks_c,
        INTEGER(super_ensembles),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP serialize_stc_ensemble_overlap_as_csv_call(SEXP filename, SEXP ensemble_masks) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int n_vectors = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[0];
    int n_selected_seed = INTEGER(Rf_getAttrib(ensemble_masks, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);
    unsigned char* ensemble_masks_c = tox_bool_in(ensemble_masks);

    // outputs and work space
    int ierr = 0;

    serialize_stc_ensemble_overlap_as_csv_c(
        filename_c,
        &filename_strlen,
        &n_vectors,
        &n_selected_seed,
        ensemble_masks_c,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP serialize_stc_super_ensembles_as_tsv_call(SEXP filename, SEXP super_ensembles) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int max_group_size = INTEGER(Rf_getAttrib(super_ensembles, R_DimSymbol))[0];
    int n_super_ensembles = INTEGER(Rf_getAttrib(super_ensembles, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);

    // outputs and work space
    int ierr = 0;

    serialize_stc_super_ensembles_as_tsv_c(
        filename_c,
        &filename_strlen,
        &max_group_size,
        &n_super_ensembles,
        INTEGER(super_ensembles),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
