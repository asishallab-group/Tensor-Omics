// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void build_kd_index_c(const double*, const int*, const int*, int*, const int*, int*);
void build_spherical_kd_c(const double*, const int*, const int*, int*, const int*, int*);
void kd_knn_query_c(const double*, const int*, const int*, const int*, const int*, const double*, const int*, int*, double*, int*);
void kd_range_query_mask_c(const double*, const int*, const int*, const int*, const int*, const double*, const double*, unsigned char*, int*);
void kd_range_query_list_c(const double*, const int*, const int*, const int*, const int*, const double*, const double*, int*, int*, int*);
void kd_range_query_count_c(const double*, const int*, const int*, const int*, const int*, const double*, const double*, int*, int*);

SEXP build_kd_index_call(SEXP points, SEXP dimension_order) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");

    // outputs and work space
    SEXP kd_indices = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    int ierr = 0;

    build_kd_index_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        INTEGER(kd_indices),
        _c_dimension_order,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, kd_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("kd_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP build_spherical_kd_call(SEXP points, SEXP dimension_order) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");

    // outputs and work space
    SEXP kd_indices = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    int ierr = 0;

    build_spherical_kd_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        INTEGER(kd_indices),
        _c_dimension_order,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, kd_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("kd_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP kd_knn_query_call(SEXP points, SEXP kd_indices, SEXP dimension_order, SEXP query_point, SEXP k_neighbors) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_k_neighbors = Rf_asInteger(k_neighbors);

    // convert what Fortran cannot take from R directly
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");

    // outputs and work space
    SEXP neighbors = PROTECT(Rf_allocVector(INTSXP, _v_k_neighbors)); _nprot++;
    SEXP distances = PROTECT(Rf_allocVector(REALSXP, _v_k_neighbors)); _nprot++;
    int ierr = 0;

    kd_knn_query_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        _c_kd_indices,
        _c_dimension_order,
        REAL(query_point),
        &_v_k_neighbors,
        INTEGER(neighbors),
        REAL(distances),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, neighbors);
    SET_VECTOR_ELT(_out, 1, distances);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("distances"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP kd_range_query_mask_call(SEXP points, SEXP kd_indices, SEXP dimension_order, SEXP query_point, SEXP radius) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_radius = Rf_asReal(radius);

    // convert what Fortran cannot take from R directly
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");

    // outputs and work space
    unsigned char* _c_in_radius_mask = tox_bool_alloc(n_points);
    int ierr = 0;

    kd_range_query_mask_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        _c_kd_indices,
        _c_dimension_order,
        REAL(query_point),
        &_v_radius,
        _c_in_radius_mask,
        &ierr
    );

    // convert the outputs back
    SEXP in_radius_mask = PROTECT(tox_bool_out(_c_in_radius_mask, n_points)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, in_radius_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("in_radius_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP kd_range_query_list_call(SEXP points, SEXP kd_indices, SEXP dimension_order, SEXP query_point, SEXP radius) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_radius = Rf_asReal(radius);

    // convert what Fortran cannot take from R directly
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");

    // outputs and work space
    SEXP neighbors = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    int n_found = 0;
    int ierr = 0;

    kd_range_query_list_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        _c_kd_indices,
        _c_dimension_order,
        REAL(query_point),
        &_v_radius,
        INTEGER(neighbors),
        &n_found,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, neighbors);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_found));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighbors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_found"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP kd_range_query_count_call(SEXP points, SEXP kd_indices, SEXP dimension_order, SEXP query_point, SEXP radius) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(points, R_DimSymbol))[0];
    int n_points = INTEGER(Rf_getAttrib(points, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    double _v_radius = Rf_asReal(radius);

    // convert what Fortran cannot take from R directly
    const int* _c_kd_indices = tox_int_in(kd_indices, "kd_indices");
    const int* _c_dimension_order = tox_int_in(dimension_order, "dimension_order");

    // outputs and work space
    int neighbor_count = 0;
    int ierr = 0;

    kd_range_query_count_c(
        REAL(points),
        &n_dimensions,
        &n_points,
        _c_kd_indices,
        _c_dimension_order,
        REAL(query_point),
        &_v_radius,
        &neighbor_count,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(neighbor_count));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neighbor_count"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
