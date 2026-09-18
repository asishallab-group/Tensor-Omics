// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void build_kd_index_c(const double*, const int*, const int*, int*, const int*, int*);
void build_spherical_kd_c(const double*, const int*, const int*, int*, const int*, int*);
void vicinity_vectors_c(const double*, const double*, const int*, const int*, const double*, const int*, const int*, unsigned char*, int*);
void vicinity_vectors_count_c(const double*, const double*, const int*, const int*, const double*, const int*, const int*, int*, int*);

SEXP build_kd_index_call(SEXP points, SEXP dimension_order) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // outputs and work space
    SEXP kd_indices = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int ierr = 0;

    build_kd_index_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, kd_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("kd_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP build_spherical_kd_call(SEXP points, SEXP dimension_order) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // outputs and work space
    SEXP kd_indices = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    int ierr = 0;

    build_spherical_kd_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, kd_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("kd_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP vicinity_vectors_call(SEXP query_point, SEXP points, SEXP r, SEXP dimension_order, SEXP kd_indices) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = (int) Rf_length(query_point);
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double r_v = Rf_asReal(r);

    // outputs and work space
    unsigned char* vicinity_mask_c = tox_bool_alloc(n_points);
    int ierr = 0;

    vicinity_vectors_c(
        REAL(query_point),
        REAL(points),
        &n_dimensions,
        &n_points,
        &r_v,
        INTEGER(dimension_order),
        INTEGER(kd_indices),
        vicinity_mask_c,
        &ierr
    );

    // convert the outputs back
    SEXP vicinity_mask = PROTECT(tox_bool_out(vicinity_mask_c, n_points)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, vicinity_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("vicinity_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP vicinity_vectors_count_call(SEXP query_point, SEXP points, SEXP r, SEXP dimension_order, SEXP kd_indices) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = (int) Rf_length(query_point);
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double r_v = Rf_asReal(r);

    // outputs and work space
    int n_neighbors = 0;
    int ierr = 0;

    vicinity_vectors_count_c(
        REAL(query_point),
        REAL(points),
        &n_dimensions,
        &n_points,
        &r_v,
        INTEGER(dimension_order),
        INTEGER(kd_indices),
        &n_neighbors,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_neighbors));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
