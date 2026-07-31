// Generated. Do not edit.
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void build_kd_index_c(const double*, const int*, const int*, int*, const int*, int*);
void build_spherical_kd_c(const double*, const int*, const int*, int*, const int*, int*);

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
