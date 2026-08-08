// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void calc_ensemble_growth_radius_c(const double*, const int*, const int*, const int*, const int*, const int*, const int*, const double*, double*, int*);
void grow_ensemble_c(const double*, const int*, const int*, const int*, const int*, const unsigned char*, const double*, unsigned char*, int*);

SEXP calc_ensemble_growth_radius_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP seed_index, SEXP k_min, SEXP radius_percentile) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int seed_index_v = Rf_asInteger(seed_index);
    int k_min_v = Rf_asInteger(k_min);
    double radius_percentile_v = Rf_asReal(radius_percentile);

    // outputs and work space
    double growth_radius = 0;
    int ierr = 0;

    calc_ensemble_growth_radius_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        &seed_index_v,
        &k_min_v,
        &radius_percentile_v,
        &growth_radius,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(growth_radius));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("growth_radius"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP grow_ensemble_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP is_member_mask, SEXP growth_radius) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double growth_radius_v = Rf_asReal(growth_radius);

    // convert what Fortran cannot take from R directly
    unsigned char* is_member_mask_c = tox_bool_in(is_member_mask);

    // outputs and work space
    unsigned char* is_member_mask_next_c = tox_bool_alloc(n_vectors);
    int ierr = 0;

    grow_ensemble_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        is_member_mask_c,
        &growth_radius_v,
        is_member_mask_next_c,
        &ierr
    );

    // convert the outputs back
    SEXP is_member_mask_next = PROTECT(tox_bool_out(is_member_mask_next_c, n_vectors)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, is_member_mask_next);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_member_mask_next"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
