// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void omics_vector_RAP_projection_c(const double*, const int*, const int*, const unsigned char*, const int*, const unsigned char*, const int*, double*, int*);
void omics_field_RAP_projection_c(const double*, const int*, const int*, const unsigned char*, const int*, const unsigned char*, const int*, double*, int*);
void clock_hand_angle_between_vectors_c(const double*, const double*, const int*, double*, const int*, int*);
void clock_hand_angles_for_shift_vectors_c(const double*, const int*, const int*, const unsigned char*, const int*, const int*, double*, int*);
void compute_relative_axis_contributions_c(const double*, const int*, double*, int*);
void relative_axes_changes_from_shift_vector_c(const double*, const int*, double*, int*);
void relative_axes_expression_from_expression_vector_c(const double*, const int*, double*, int*);

SEXP omics_vector_RAP_projection_call(SEXP vecs, SEXP vecs_selection_mask, SEXP axes_selection_mask) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(vecs, R_DimSymbol))[0];
    int n_vecs = INTEGER(Rf_getAttrib(vecs, R_DimSymbol))[1];
    int n_selected_vecs = tox_sum_true(vecs_selection_mask);
    int n_selected_axes = tox_sum_true(axes_selection_mask);

    // convert what Fortran cannot take from R directly
    unsigned char* vecs_selection_mask_c = tox_bool_in(vecs_selection_mask);
    unsigned char* axes_selection_mask_c = tox_bool_in(axes_selection_mask);

    // outputs and work space
    SEXP projections = PROTECT(Rf_allocVector(REALSXP, n_selected_axes * n_selected_vecs)); nprot++;
    { SEXP projections_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(projections_dim)[0] = n_selected_axes; INTEGER(projections_dim)[1] = n_selected_vecs; Rf_setAttrib(projections, R_DimSymbol, projections_dim); UNPROTECT(1); }
    int ierr = 0;

    omics_vector_RAP_projection_c(
        REAL(vecs),
        &n_axes,
        &n_vecs,
        vecs_selection_mask_c,
        &n_selected_vecs,
        axes_selection_mask_c,
        &n_selected_axes,
        REAL(projections),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, projections);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("projections"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP omics_field_RAP_projection_call(SEXP fields, SEXP fields_selection_mask, SEXP axes_selection_mask) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(fields, R_DimSymbol))[0];
    int n_fields = INTEGER(Rf_getAttrib(fields, R_DimSymbol))[2];
    int n_selected_fields = tox_sum_true(fields_selection_mask);
    int n_selected_axes = tox_sum_true(axes_selection_mask);

    // convert what Fortran cannot take from R directly
    unsigned char* fields_selection_mask_c = tox_bool_in(fields_selection_mask);
    unsigned char* axes_selection_mask_c = tox_bool_in(axes_selection_mask);

    // outputs and work space
    SEXP projections = PROTECT(Rf_allocVector(REALSXP, n_selected_axes * n_selected_fields)); nprot++;
    { SEXP projections_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(projections_dim)[0] = n_selected_axes; INTEGER(projections_dim)[1] = n_selected_fields; Rf_setAttrib(projections, R_DimSymbol, projections_dim); UNPROTECT(1); }
    int ierr = 0;

    omics_field_RAP_projection_c(
        REAL(fields),
        &n_axes,
        &n_fields,
        fields_selection_mask_c,
        &n_selected_fields,
        axes_selection_mask_c,
        &n_selected_axes,
        REAL(projections),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, projections);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("projections"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP clock_hand_angle_between_vectors_call(SEXP v1, SEXP v2, SEXP selected_axes_for_signed) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = (int) Rf_length(v1);

    // outputs and work space
    double signed_angle = 0;
    int ierr = 0;

    clock_hand_angle_between_vectors_c(
        REAL(v1),
        REAL(v2),
        &n_dims,
        &signed_angle,
        INTEGER(selected_axes_for_signed),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(signed_angle));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("signed_angle"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP clock_hand_angles_for_shift_vectors_call(SEXP fields, SEXP fields_selection_mask, SEXP selected_axes_for_signed) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = INTEGER(Rf_getAttrib(fields, R_DimSymbol))[0];
    int n_fields = INTEGER(Rf_getAttrib(fields, R_DimSymbol))[2];
    int n_selected_fields = tox_sum_true(fields_selection_mask);

    // convert what Fortran cannot take from R directly
    unsigned char* fields_selection_mask_c = tox_bool_in(fields_selection_mask);

    // outputs and work space
    SEXP signed_angles = PROTECT(Rf_allocVector(REALSXP, n_selected_fields)); nprot++;
    int ierr = 0;

    clock_hand_angles_for_shift_vectors_c(
        REAL(fields),
        &n_dims,
        &n_fields,
        fields_selection_mask_c,
        &n_selected_fields,
        INTEGER(selected_axes_for_signed),
        REAL(signed_angles),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, signed_angles);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("signed_angles"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_relative_axis_contributions_call(SEXP vec) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = (int) Rf_length(vec);

    // outputs and work space
    SEXP contributions = PROTECT(Rf_allocVector(REALSXP, n_axes)); nprot++;
    int ierr = 0;

    compute_relative_axis_contributions_c(
        REAL(vec),
        &n_axes,
        REAL(contributions),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, contributions);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP relative_axes_changes_from_shift_vector_call(SEXP vec) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = (int) Rf_length(vec);

    // outputs and work space
    SEXP contributions = PROTECT(Rf_allocVector(REALSXP, n_axes)); nprot++;
    int ierr = 0;

    relative_axes_changes_from_shift_vector_c(
        REAL(vec),
        &n_axes,
        REAL(contributions),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, contributions);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP relative_axes_expression_from_expression_vector_call(SEXP vec) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = (int) Rf_length(vec);

    // outputs and work space
    SEXP contributions = PROTECT(Rf_allocVector(REALSXP, n_axes)); nprot++;
    int ierr = 0;

    relative_axes_expression_from_expression_vector_c(
        REAL(vec),
        &n_axes,
        REAL(contributions),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, contributions);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("contributions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
