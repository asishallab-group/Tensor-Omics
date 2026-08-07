// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void density_labels_c(const double*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*);
void seeds_c(const double*, const int*, const int*, const int*, const int*, const int*, const double*, unsigned char*, int*);

SEXP density_labels_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP k_density, SEXP bandwidth_percentile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int k_density_v = Rf_asInteger(k_density);
    double bandwidth_percentile_v = Rf_asReal(bandwidth_percentile);

    // outputs and work space
    SEXP labels = PROTECT(Rf_allocVector(REALSXP, n_vectors)); nprot++;
    int ierr = 0;

    density_labels_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        &k_density_v,
        &bandwidth_percentile_v,
        REAL(labels),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, labels);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("labels"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP seeds_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP k_density, SEXP bandwidth_percentile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int k_density_v = Rf_asInteger(k_density);
    double bandwidth_percentile_v = Rf_asReal(bandwidth_percentile);

    // outputs and work space
    unsigned char* is_seed_mask_c = tox_bool_alloc(n_vectors);
    int ierr = 0;

    seeds_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        &k_density_v,
        &bandwidth_percentile_v,
        is_seed_mask_c,
        &ierr
    );

    // convert the outputs back
    SEXP is_seed_mask = PROTECT(tox_bool_out(is_seed_mask_c, n_vectors)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, is_seed_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_seed_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
