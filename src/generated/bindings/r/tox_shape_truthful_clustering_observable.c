// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void normal_error_c(const int*, const double*, const int*, double*, int*);
void tangent_scales_c(const int*, const double*, const int*, double*, int*);
void observable_c(const double*, const int*, const int*, const unsigned char*, const int*, double*, double*, double*, int*, double*, double*, double*, int*);

SEXP normal_error_call(SEXP d, SEXP eigenvalues) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = (int) Rf_length(eigenvalues);

    // scalar inputs, pulled from their length-1 vectors
    int d_v = Rf_asInteger(d);

    // outputs and work space
    double normal_error_value = 0;
    int ierr = 0;

    normal_error_c(
        &d_v,
        REAL(eigenvalues),
        &n_dimensions,
        &normal_error_value,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(normal_error_value));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("normal_error_value"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP tangent_scales_call(SEXP d, SEXP eigenvalues) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = (int) Rf_length(eigenvalues);

    // scalar inputs, pulled from their length-1 vectors
    int d_v = Rf_asInteger(d);

    // outputs and work space
    SEXP tangent_scales_value = PROTECT(Rf_allocVector(REALSXP, d_v)); nprot++;
    int ierr = 0;

    tangent_scales_c(
        &d_v,
        REAL(eigenvalues),
        &n_dimensions,
        REAL(tangent_scales_value),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, tangent_scales_value);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("tangent_scales_value"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP observable_call(SEXP vectors, SEXP member_selection_mask) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_selected_member = tox_sum_true(member_selection_mask);

    // convert what Fortran cannot take from R directly
    unsigned char* member_selection_mask_c = tox_bool_in(member_selection_mask);

    // outputs and work space
    SEXP U = PROTECT(Rf_allocVector(REALSXP, n_dimensions * n_dimensions)); nprot++;
    { SEXP U_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(U_dim)[0] = n_dimensions; INTEGER(U_dim)[1] = n_dimensions; Rf_setAttrib(U, R_DimSymbol, U_dim); UNPROTECT(1); }
    SEXP eigenvalues = PROTECT(Rf_allocVector(REALSXP, n_dimensions)); nprot++;
    SEXP mu = PROTECT(Rf_allocVector(REALSXP, n_dimensions)); nprot++;
    int d = 0;
    double G = 0;
    double normal_error_value = 0;
    SEXP tangent_scales_value = PROTECT(Rf_allocVector(REALSXP, n_dimensions)); nprot++;
    int ierr = 0;

    observable_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        member_selection_mask_c,
        &n_selected_member,
        REAL(U),
        REAL(eigenvalues),
        REAL(mu),
        &d,
        &G,
        &normal_error_value,
        REAL(tangent_scales_value),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 8)); nprot++;
    SET_VECTOR_ELT(_out, 0, U);
    SET_VECTOR_ELT(_out, 1, eigenvalues);
    SET_VECTOR_ELT(_out, 2, mu);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(d));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(G));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarReal(normal_error_value));
    SET_VECTOR_ELT(_out, 6, tangent_scales_value);
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 8)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("U"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("eigenvalues"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("mu"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("d"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("G"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("normal_error_value"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("tangent_scales_value"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
