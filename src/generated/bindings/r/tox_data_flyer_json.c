// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void save_flyer_json_c(const int*, const int*, const int*, const double*, const double*, const int*, const unsigned char*, const char*, const int*, const char*, const int*, const char*, const int*, const char*, const int*, const char*, const int*, const char*, const int*, int*);

SEXP save_flyer_json_call(SEXP expression_vectors, SEXP family_centroids, SEXP gene_to_fam, SEXP is_outlier, SEXP axis_labels, SEXP family_ids, SEXP gene_ids, SEXP gene_species, SEXP gene_types, SEXP filename) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_families = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[1];
    int axis_labels_strlen = tox_max_strlen(axis_labels);
    int family_ids_strlen = tox_max_strlen(family_ids);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int gene_species_strlen = tox_max_strlen(gene_species);
    int gene_types_strlen = tox_max_strlen(gene_types);
    int filename_strlen = tox_max_strlen(filename);

    // convert what Fortran cannot take from R directly
    unsigned char* is_outlier_c = tox_bool_in(is_outlier);
    char* axis_labels_c = tox_char_in(axis_labels, axis_labels_strlen);
    char* family_ids_c = tox_char_in(family_ids, family_ids_strlen);
    char* gene_ids_c = tox_char_in(gene_ids, gene_ids_strlen);
    char* gene_species_c = tox_char_in(gene_species, gene_species_strlen);
    char* gene_types_c = tox_char_in(gene_types, gene_types_strlen);
    char* filename_c = tox_char_in(filename, filename_strlen);

    // outputs and work space
    int ierr = 0;

    save_flyer_json_c(
        &n_axes,
        &n_genes,
        &n_families,
        REAL(expression_vectors),
        REAL(family_centroids),
        INTEGER(gene_to_fam),
        is_outlier_c,
        axis_labels_c,
        &axis_labels_strlen,
        family_ids_c,
        &family_ids_strlen,
        gene_ids_c,
        &gene_ids_strlen,
        gene_species_c,
        &gene_species_strlen,
        gene_types_c,
        &gene_types_strlen,
        filename_c,
        &filename_strlen,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
