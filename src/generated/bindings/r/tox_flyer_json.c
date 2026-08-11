// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void serialize_tox_data_as_flyer_json_c(const char*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, const double*, const char*, const int*, const int*, const double*, const int*, const int*, const unsigned char*, const char*, const int*, const char*, const int*, int*);

SEXP serialize_tox_data_as_flyer_json_call(SEXP filename, SEXP tissues, SEXP family_ids, SEXP centroids, SEXP gene_ids, SEXP genes, SEXP gene_to_fam, SEXP sorted_gene_to_fam_perm, SEXP gene_outliers, SEXP gene_species, SEXP gene_types) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int tissues_strlen = tox_max_strlen(tissues);
    int n_tissues = (int) Rf_length(tissues);
    int family_ids_strlen = tox_max_strlen(family_ids);
    int n_families = (int) Rf_length(family_ids);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_genes = (int) Rf_length(gene_ids);
    int gene_species_strlen = tox_max_strlen(gene_species);
    int gene_types_strlen = tox_max_strlen(gene_types);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);
    char* tissues_c = tox_char_in(tissues, tissues_strlen);
    char* family_ids_c = tox_char_in(family_ids, family_ids_strlen);
    char* gene_ids_c = tox_char_in(gene_ids, gene_ids_strlen);
    unsigned char* gene_outliers_c = tox_bool_in(gene_outliers);
    char* gene_species_c = tox_char_in(gene_species, gene_species_strlen);
    char* gene_types_c = tox_char_in(gene_types, gene_types_strlen);

    // outputs and work space
    int ierr = 0;

    serialize_tox_data_as_flyer_json_c(
        filename_c,
        &filename_strlen,
        tissues_c,
        &tissues_strlen,
        &n_tissues,
        family_ids_c,
        &family_ids_strlen,
        &n_families,
        REAL(centroids),
        gene_ids_c,
        &gene_ids_strlen,
        &n_genes,
        REAL(genes),
        INTEGER(gene_to_fam),
        INTEGER(sorted_gene_to_fam_perm),
        gene_outliers_c,
        gene_species_c,
        &gene_species_strlen,
        gene_types_c,
        &gene_types_strlen,
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
