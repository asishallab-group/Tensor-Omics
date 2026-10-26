// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void validate_data_structure_c(const int*, const int*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, int*);
void validate_gene_to_family_mapping_c(const int*, const int*, const int*, int*);
void validate_expression_data_c(const double*, const int*, const int*, const unsigned char*, int*);
void validate_family_centroids_c(const double*, const int*, const int*, int*);
void validate_shift_vectors_c(const double*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const int*, const int*, const int*, int*);
void validate_string_array_uniqueness_c(const char*, const int*, const int*, int*);
void validate_all_data_c(const int*, const int*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const double*, const int*, const int*, const unsigned char*, const unsigned char*, int*);

SEXP validate_data_structure_call(SEXP n_genes, SEXP n_families, SEXP n_samples, SEXP gene_ids, SEXP gene_family_ids, SEXP gene_to_fam, SEXP expression_vectors, SEXP family_centroids, SEXP shift_vectors) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = (int) Rf_length(gene_ids);
    int gene_family_ids_strlen = tox_max_strlen(gene_family_ids);
    int n_gene_family_ids_elements = (int) Rf_length(gene_family_ids);
    int n_gene_to_fam_elements = (int) Rf_length(gene_to_fam);
    int n_expression_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_expression_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_family_centroids_elements_dim_1 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[0];
    int n_family_centroids_elements_dim_2 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[1];
    int n_shift_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(shift_vectors, R_DimSymbol))[0];
    int n_shift_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(shift_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_genes = Rf_asInteger(n_genes);
    int _v_n_families = Rf_asInteger(n_families);
    int _v_n_samples = Rf_asInteger(n_samples);

    // convert what Fortran cannot take from R directly
    char* _c_gene_ids = tox_char_in(gene_ids, gene_ids_strlen);
    char* _c_gene_family_ids = tox_char_in(gene_family_ids, gene_family_ids_strlen);

    // outputs and work space
    int ierr = 0;

    validate_data_structure_c(
        &_v_n_genes,
        &_v_n_families,
        &_v_n_samples,
        _c_gene_ids,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        _c_gene_family_ids,
        &gene_family_ids_strlen,
        &n_gene_family_ids_elements,
        INTEGER(gene_to_fam),
        &n_gene_to_fam_elements,
        REAL(expression_vectors),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        REAL(family_centroids),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        REAL(shift_vectors),
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
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

SEXP validate_gene_to_family_mapping_call(SEXP gene_to_fam, SEXP n_families) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_gene_to_fam_elements = (int) Rf_length(gene_to_fam);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_families = Rf_asInteger(n_families);

    // outputs and work space
    int ierr = 0;

    validate_gene_to_family_mapping_c(
        INTEGER(gene_to_fam),
        &n_gene_to_fam_elements,
        &_v_n_families,
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

SEXP validate_expression_data_call(SEXP expression_vectors, SEXP check_non_negative) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_expression_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_expression_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    unsigned char _v_check_non_negative = (Rf_asLogical(check_non_negative) == TRUE) ? 1 : 0;

    // outputs and work space
    int ierr = 0;

    validate_expression_data_c(
        REAL(expression_vectors),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        &_v_check_non_negative,
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

SEXP validate_family_centroids_call(SEXP family_centroids) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_family_centroids_elements_dim_1 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[0];
    int n_family_centroids_elements_dim_2 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[1];

    // outputs and work space
    int ierr = 0;

    validate_family_centroids_c(
        REAL(family_centroids),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
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

SEXP validate_shift_vectors_call(SEXP shift_vectors, SEXP expression_vectors, SEXP family_centroids, SEXP gene_to_fam, SEXP n_samples) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_shift_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(shift_vectors, R_DimSymbol))[0];
    int n_shift_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(shift_vectors, R_DimSymbol))[1];
    int n_expression_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_expression_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_family_centroids_elements_dim_1 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[0];
    int n_family_centroids_elements_dim_2 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[1];
    int n_gene_to_fam_elements = (int) Rf_length(gene_to_fam);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_samples = Rf_asInteger(n_samples);

    // outputs and work space
    int ierr = 0;

    validate_shift_vectors_c(
        REAL(shift_vectors),
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        REAL(expression_vectors),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        REAL(family_centroids),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        INTEGER(gene_to_fam),
        &n_gene_to_fam_elements,
        &_v_n_samples,
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

SEXP validate_string_array_uniqueness_call(SEXP str_arr) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int str_arr_strlen = tox_max_strlen(str_arr);
    int n_str_arr_elements = (int) Rf_length(str_arr);

    // convert what Fortran cannot take from R directly
    char* _c_str_arr = tox_char_in(str_arr, str_arr_strlen);

    // outputs and work space
    int ierr = 0;

    validate_string_array_uniqueness_c(
        _c_str_arr,
        &str_arr_strlen,
        &n_str_arr_elements,
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

SEXP validate_all_data_call(SEXP n_genes, SEXP n_families, SEXP n_samples, SEXP gene_ids, SEXP gene_family_ids, SEXP gene_to_fam, SEXP expression_vectors, SEXP family_centroids, SEXP shift_vectors, SEXP check_uniqueness, SEXP check_shift_consistency) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = (int) Rf_length(gene_ids);
    int gene_family_ids_strlen = tox_max_strlen(gene_family_ids);
    int n_gene_family_ids_elements = (int) Rf_length(gene_family_ids);
    int n_gene_to_fam_elements = (int) Rf_length(gene_to_fam);
    int n_expression_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_expression_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_family_centroids_elements_dim_1 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[0];
    int n_family_centroids_elements_dim_2 = INTEGER(Rf_getAttrib(family_centroids, R_DimSymbol))[1];
    int n_shift_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(shift_vectors, R_DimSymbol))[0];
    int n_shift_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(shift_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_genes = Rf_asInteger(n_genes);
    int _v_n_families = Rf_asInteger(n_families);
    int _v_n_samples = Rf_asInteger(n_samples);
    unsigned char _v_check_uniqueness = (Rf_asLogical(check_uniqueness) == TRUE) ? 1 : 0;
    unsigned char _v_check_shift_consistency = (Rf_asLogical(check_shift_consistency) == TRUE) ? 1 : 0;

    // convert what Fortran cannot take from R directly
    char* _c_gene_ids = tox_char_in(gene_ids, gene_ids_strlen);
    char* _c_gene_family_ids = tox_char_in(gene_family_ids, gene_family_ids_strlen);

    // outputs and work space
    int ierr = 0;

    validate_all_data_c(
        &_v_n_genes,
        &_v_n_families,
        &_v_n_samples,
        _c_gene_ids,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        _c_gene_family_ids,
        &gene_family_ids_strlen,
        &n_gene_family_ids_elements,
        INTEGER(gene_to_fam),
        &n_gene_to_fam_elements,
        REAL(expression_vectors),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        REAL(family_centroids),
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        REAL(shift_vectors),
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        &_v_check_uniqueness,
        &_v_check_shift_consistency,
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
