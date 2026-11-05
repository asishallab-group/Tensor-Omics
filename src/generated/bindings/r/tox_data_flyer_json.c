// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void save_flyer_json_c(const int*, const int*, const int*, const double*, const double*, const int*, const unsigned char*, const char*, const int*, const char*, const int*, const char*, const int*, const char*, const int*, const char*, const int*, const char*, const int*, int*);

SEXP save_flyer_json_call(SEXP expression_vectors, SEXP family_centroids, SEXP gene_to_fam, SEXP is_outlier, SEXP axis_labels, SEXP family_ids, SEXP gene_ids, SEXP gene_species, SEXP gene_types, SEXP filename) {
    int _nprot = 0;
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
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");
    unsigned char* _c_is_outlier = tox_bool_in(is_outlier, "is_outlier");
    char* _c_axis_labels = tox_char_in(axis_labels, axis_labels_strlen, "axis_labels");
    char* _c_family_ids = tox_char_in(family_ids, family_ids_strlen, "family_ids");
    char* _c_gene_ids = tox_char_in(gene_ids, gene_ids_strlen, "gene_ids");
    char* _c_gene_species = tox_char_in(gene_species, gene_species_strlen, "gene_species");
    char* _c_gene_types = tox_char_in(gene_types, gene_types_strlen, "gene_types");
    char* _c_filename = tox_char_in(filename, filename_strlen, "filename");

    // outputs and work space
    int ierr = 0;

    save_flyer_json_c(
        &n_axes,
        &n_genes,
        &n_families,
        REAL(expression_vectors),
        REAL(family_centroids),
        _c_gene_to_fam,
        _c_is_outlier,
        _c_axis_labels,
        &axis_labels_strlen,
        _c_family_ids,
        &family_ids_strlen,
        _c_gene_ids,
        &gene_ids_strlen,
        _c_gene_species,
        &gene_species_strlen,
        _c_gene_types,
        &gene_types_strlen,
        _c_filename,
        &filename_strlen,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
