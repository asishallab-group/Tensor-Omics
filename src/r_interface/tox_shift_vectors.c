// Generated. Do not edit.
#if !defined(NO_R_INTERFACE) && !defined(NO_C_INTERFACE)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void compute_shift_vector_field_c(const int*, const int*, const int*, const double*, const double*, const int*, double*, int*);

SEXP compute_shift_vector_field_call(SEXP expression_vectors, SEXP family_centroids, SEXP gene_to_fam) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_tissues = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_families = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[1];

    // outputs and work space
    SEXP shift_vectors = PROTECT(Rf_allocVector(REALSXP, n_tissues * 2 * n_genes)); nprot++;
    { SEXP shift_vectors_dim = PROTECT(Rf_allocVector(INTSXP, 3)); INTEGER(shift_vectors_dim)[0] = n_tissues; INTEGER(shift_vectors_dim)[1] = 2; INTEGER(shift_vectors_dim)[2] = n_genes; Rf_setAttrib(shift_vectors, R_DimSymbol, shift_vectors_dim); UNPROTECT(1); }
    int ierr = 0;

    compute_shift_vector_field_c(
        &n_tissues,
        &n_genes,
        &n_families,
        REAL(expression_vectors),
        REAL(family_centroids),
        INTEGER(gene_to_fam),
        REAL(shift_vectors),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, shift_vectors);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("shift_vectors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R interface
