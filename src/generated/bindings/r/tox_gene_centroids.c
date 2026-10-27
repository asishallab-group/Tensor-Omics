// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void mean_vector_c(const double*, const int*, const int*, const int*, const int*, double*, int*);
void group_centroid_orthologs_c(const double*, const int*, const int*, const int*, const int*, double*, const unsigned char*, int*);
void group_centroid_all_c(const double*, const int*, const int*, const int*, const int*, double*, int*);

SEXP mean_vector_call(SEXP expression_vectors, SEXP gene_indices) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_selected_genes = (int) Rf_length(gene_indices);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_indices = tox_int_in(gene_indices, "gene_indices");

    // outputs and work space
    SEXP centroid = PROTECT(Rf_allocVector(REALSXP, n_axes)); _nprot++;
    int ierr = 0;

    mean_vector_c(
        REAL(expression_vectors),
        &n_axes,
        &n_genes,
        _c_gene_indices,
        &n_selected_genes,
        REAL(centroid),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, centroid);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroid"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP group_centroid_orthologs_call(SEXP expression_vectors, SEXP gene_to_family, SEXP n_families, SEXP ortholog_set) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_families = Rf_asInteger(n_families);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_family = tox_int_in(gene_to_family, "gene_to_family");
    unsigned char* _c_ortholog_set = tox_bool_in(ortholog_set, "ortholog_set");

    // outputs and work space
    SEXP centroid_matrix = PROTECT(Rf_allocVector(REALSXP, n_axes * _v_n_families)); _nprot++;
    { SEXP _dim_centroid_matrix = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_centroid_matrix)[0] = n_axes; INTEGER(_dim_centroid_matrix)[1] = _v_n_families; Rf_setAttrib(centroid_matrix, R_DimSymbol, _dim_centroid_matrix); UNPROTECT(1); }
    int ierr = 0;

    group_centroid_orthologs_c(
        REAL(expression_vectors),
        &n_axes,
        &n_genes,
        _c_gene_to_family,
        &_v_n_families,
        REAL(centroid_matrix),
        _c_ortholog_set,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, centroid_matrix);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroid_matrix"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP group_centroid_all_call(SEXP expression_vectors, SEXP gene_to_family, SEXP n_families) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_families = Rf_asInteger(n_families);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_family = tox_int_in(gene_to_family, "gene_to_family");

    // outputs and work space
    SEXP centroid_matrix = PROTECT(Rf_allocVector(REALSXP, n_axes * _v_n_families)); _nprot++;
    { SEXP _dim_centroid_matrix = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_centroid_matrix)[0] = n_axes; INTEGER(_dim_centroid_matrix)[1] = _v_n_families; Rf_setAttrib(centroid_matrix, R_DimSymbol, _dim_centroid_matrix); UNPROTECT(1); }
    int ierr = 0;

    group_centroid_all_c(
        REAL(expression_vectors),
        &n_axes,
        &n_genes,
        _c_gene_to_family,
        &_v_n_families,
        REAL(centroid_matrix),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, centroid_matrix);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroid_matrix"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
