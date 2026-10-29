// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void read_expression_vectors_tsv_c(const char*, const int*, const int*, const char*, const int*, const int*, double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, const char*, int*);
void read_gene_ids_from_tsv_file_c(const char*, const int*, char*, const int*, const int*, const int*, const int*, int*);
void read_orthofinder_file_c(const char*, const int*, const char*, const int*, const int*, char*, const int*, const int*, int*, const int*, int*);
void get_unassigned_mask_c(const int*, const int*, unsigned char*, int*, int*);

SEXP read_expression_vectors_tsv_call(SEXP file_list, SEXP gene_ids, SEXP expression_vectors, SEXP n_header_rows, SEXP gene_col, SEXP value_cols, SEXP start_row, SEXP delimiter) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int file_list_strlen = tox_max_strlen(file_list);
    int n_file_list_elements = (int) Rf_length(file_list);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = (int) Rf_length(gene_ids);
    int n_expression_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_expression_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_value_cols_elements = (int) Rf_length(value_cols);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_header_rows = Rf_asInteger(n_header_rows);
    int _v_gene_col = Rf_asInteger(gene_col);
    int _v_start_row = Rf_asInteger(start_row);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_expression_vectors = PROTECT(Rf_duplicate(expression_vectors)); _nprot++;

    // convert what Fortran cannot take from R directly
    char* _c_file_list = tox_char_in(file_list, file_list_strlen, "file_list");
    char* _c_gene_ids = tox_char_in(gene_ids, gene_ids_strlen, "gene_ids");
    const int* _c_value_cols = tox_int_in(value_cols, "value_cols");
    char* _c_delimiter = tox_char_in(delimiter, 1, "delimiter");

    // outputs and work space
    int ierr = 0;

    read_expression_vectors_tsv_c(
        _c_file_list,
        &file_list_strlen,
        &n_file_list_elements,
        _c_gene_ids,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        REAL(_out_expression_vectors),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        &_v_n_header_rows,
        &_v_gene_col,
        _c_value_cols,
        &n_value_cols_elements,
        &_v_start_row,
        _c_delimiter,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_expression_vectors);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("expression_vectors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP read_gene_ids_from_tsv_file_call(SEXP filename, SEXP gene_ids_strlen, SEXP n_gene_ids_elements, SEXP n_header_rows, SEXP gene_col) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);

    // scalar inputs, pulled from their length-1 vectors
    int _v_gene_ids_strlen = Rf_asInteger(gene_ids_strlen);
    int _v_n_gene_ids_elements = Rf_asInteger(n_gene_ids_elements);
    int _v_n_header_rows = Rf_asInteger(n_header_rows);
    int _v_gene_col = Rf_asInteger(gene_col);

    // convert what Fortran cannot take from R directly
    char* _c_filename = tox_char_in(filename, filename_strlen, "filename");

    // outputs and work space
    char* _c_gene_ids = tox_char_alloc(_v_gene_ids_strlen, _v_n_gene_ids_elements);
    int ierr = 0;

    read_gene_ids_from_tsv_file_c(
        _c_filename,
        &filename_strlen,
        _c_gene_ids,
        &_v_gene_ids_strlen,
        &_v_n_gene_ids_elements,
        &_v_n_header_rows,
        &_v_gene_col,
        &ierr
    );

    // convert the outputs back
    SEXP gene_ids = PROTECT(tox_char_out(_c_gene_ids, _v_gene_ids_strlen, _v_n_gene_ids_elements)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, gene_ids);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("gene_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP read_orthofinder_file_call(SEXP filename, SEXP gene_ids, SEXP family_ids_strlen, SEXP n_family_ids_elements, SEXP n_gene_to_fam_elements) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = (int) Rf_length(gene_ids);

    // scalar inputs, pulled from their length-1 vectors
    int _v_family_ids_strlen = Rf_asInteger(family_ids_strlen);
    int _v_n_family_ids_elements = Rf_asInteger(n_family_ids_elements);
    int _v_n_gene_to_fam_elements = Rf_asInteger(n_gene_to_fam_elements);

    // convert what Fortran cannot take from R directly
    char* _c_filename = tox_char_in(filename, filename_strlen, "filename");
    char* _c_gene_ids = tox_char_in(gene_ids, gene_ids_strlen, "gene_ids");

    // outputs and work space
    char* _c_family_ids = tox_char_alloc(_v_family_ids_strlen, _v_n_family_ids_elements);
    SEXP gene_to_fam = PROTECT(Rf_allocVector(INTSXP, _v_n_gene_to_fam_elements)); _nprot++;
    int ierr = 0;

    read_orthofinder_file_c(
        _c_filename,
        &filename_strlen,
        _c_gene_ids,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        _c_family_ids,
        &_v_family_ids_strlen,
        &_v_n_family_ids_elements,
        INTEGER(gene_to_fam),
        &_v_n_gene_to_fam_elements,
        &ierr
    );

    // convert the outputs back
    SEXP family_ids = PROTECT(tox_char_out(_c_family_ids, _v_family_ids_strlen, _v_n_family_ids_elements)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, family_ids);
    SET_VECTOR_ELT(_out, 1, gene_to_fam);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("gene_to_fam"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP get_unassigned_mask_call(SEXP gene_to_fam) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_gene_to_fam_elements = (int) Rf_length(gene_to_fam);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    unsigned char* _c_mask = tox_bool_alloc(((int) Rf_length(gene_to_fam)));
    int n_genes_kept = 0;
    int ierr = 0;

    get_unassigned_mask_c(
        _c_gene_to_fam,
        &n_gene_to_fam_elements,
        _c_mask,
        &n_genes_kept,
        &ierr
    );

    // convert the outputs back
    SEXP mask = PROTECT(tox_bool_out(_c_mask, ((int) Rf_length(gene_to_fam)))); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_genes_kept));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_genes_kept"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
