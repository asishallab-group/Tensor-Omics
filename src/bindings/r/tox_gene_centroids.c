// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void mean_vector_c(const double*, const int*, const int*, const int*, const int*, double*, int*);
void group_centroid_c(const double*, const int*, const int*, const int*, const int*, double*, const char*, int*, int*, const unsigned char*);

SEXP mean_vector_call(SEXP expression_vectors, SEXP gene_indices) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_selected_genes = (int) Rf_length(gene_indices);

    // outputs and work space
    SEXP centroid = PROTECT(Rf_allocVector(REALSXP, n_axes)); nprot++;
    int ierr = 0;

    mean_vector_c(
        REAL(expression_vectors),
        &n_axes,
        &n_genes,
        INTEGER(gene_indices),
        &n_selected_genes,
        REAL(centroid),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, centroid);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroid"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP group_centroid_call(SEXP expression_vectors, SEXP gene_to_family, SEXP n_families, SEXP mode, SEXP ortholog_set) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int ortholog_set_size = 0;
    if (ortholog_set != R_NilValue) {
        ortholog_set_size = (int) Rf_length(ortholog_set);
    }

    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);

    // convert what Fortran cannot take from R directly
    char* mode_c = tox_char_in(mode, 15);
    unsigned char* ortholog_set_c = tox_bool_in(ortholog_set);

    // outputs and work space
    SEXP centroid_matrix = PROTECT(Rf_allocVector(REALSXP, n_axes * n_families_v)); nprot++;
    { SEXP centroid_matrix_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(centroid_matrix_dim)[0] = n_axes; INTEGER(centroid_matrix_dim)[1] = n_families_v; Rf_setAttrib(centroid_matrix, R_DimSymbol, centroid_matrix_dim); UNPROTECT(1); }
    int* tmp_group_indices = (int*) R_alloc(n_genes, sizeof(int));
    int ierr = 0;

    group_centroid_c(
        REAL(expression_vectors),
        &n_axes,
        &n_genes,
        INTEGER(gene_to_family),
        &n_families_v,
        REAL(centroid_matrix),
        mode_c,
        tmp_group_indices,
        &ierr,
        ortholog_set != R_NilValue ? ortholog_set_c : NULL
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, centroid_matrix);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroid_matrix"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
