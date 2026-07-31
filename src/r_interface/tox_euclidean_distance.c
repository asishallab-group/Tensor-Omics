// Generated. Do not edit.
#if !defined(NO_R_INTERFACE) && !defined(NO_C_INTERFACE)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void euclidean_distance_c(const double*, const double*, const int*, double*, int*);
void distance_to_centroid_c(const int*, const int*, const double*, const double*, const int*, double*, const int*, int*);

SEXP euclidean_distance_call(SEXP vec1, SEXP vec2) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_elements = (int) Rf_length(vec1);

    // outputs and work space
    double result = 0;
    int ierr = 0;

    euclidean_distance_c(
        REAL(vec1),
        REAL(vec2),
        &n_elements,
        &result,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(result));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("result"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP distance_to_centroid_call(SEXP genes, SEXP centroids, SEXP gene_to_fam) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_families = INTEGER(Rf_getAttrib(centroids, R_DimSymbol))[1];
    int n_tissues = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[0];

    // outputs and work space
    SEXP distances = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    distance_to_centroid_c(
        &n_genes,
        &n_families,
        REAL(genes),
        REAL(centroids),
        INTEGER(gene_to_fam),
        REAL(distances),
        &n_tissues,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, distances);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("distances"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R interface
