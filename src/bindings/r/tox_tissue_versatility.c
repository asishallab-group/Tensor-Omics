// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void compute_tissue_versatility_c(const int*, const int*, const double*, const unsigned char*, const int*, const unsigned char*, const int*, double*, double*, int*);

SEXP compute_tissue_versatility_call(SEXP expression_vectors, SEXP exp_vecs_selection_index, SEXP n_selected_vectors, SEXP axes_selection, SEXP n_selected_axes) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_selected_vectors_v = Rf_asInteger(n_selected_vectors);
    int n_selected_axes_v = Rf_asInteger(n_selected_axes);

    // convert what Fortran cannot take from R directly
    unsigned char* exp_vecs_selection_index_c = tox_bool_in(exp_vecs_selection_index);
    unsigned char* axes_selection_c = tox_bool_in(axes_selection);

    // outputs and work space
    SEXP tissue_versatilities = PROTECT(Rf_allocVector(REALSXP, n_selected_vectors_v)); nprot++;
    SEXP tissue_angles_deg = PROTECT(Rf_allocVector(REALSXP, n_selected_vectors_v)); nprot++;
    int ierr = 0;

    compute_tissue_versatility_c(
        &n_axes,
        &n_vectors,
        REAL(expression_vectors),
        exp_vecs_selection_index_c,
        &n_selected_vectors_v,
        axes_selection_c,
        &n_selected_axes_v,
        REAL(tissue_versatilities),
        REAL(tissue_angles_deg),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, tissue_versatilities);
    SET_VECTOR_ELT(_out, 1, tissue_angles_deg);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("tissue_versatilities"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("tissue_angles_deg"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
