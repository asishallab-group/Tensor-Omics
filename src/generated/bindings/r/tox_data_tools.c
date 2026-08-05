// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void read_expression_vectors_tsv_c(const char*, const int*, const int*, const char*, const int*, const int*, double*, const int*, const int*, const int*, const int*, const int*, const int*, const int*, int*, const char*);
void read_gene_ids_from_tsv_file_c(const char*, const int*, char*, const int*, const int*, const int*, const int*, int*);
void read_orthofinder_file_c(const char*, const int*, const char*, const int*, const int*, char*, const int*, const int*, int*, const int*, int*);
void get_unassigned_mask_c(const int*, const int*, unsigned char*, int*, int*);

SEXP read_expression_vectors_tsv_call(SEXP file_list, SEXP gene_ids, SEXP expression_vectors, SEXP n_header_rows, SEXP gene_col, SEXP value_cols, SEXP start_row, SEXP delimiter) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int file_list_strlen = tox_max_strlen(file_list);
    int n_file_list_elements = (int) Rf_length(file_list);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = (int) Rf_length(gene_ids);
    int n_expression_vectors_elements_dim_1 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_expression_vectors_elements_dim_2 = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_value_cols_elements = (int) Rf_length(value_cols);

    // scalar inputs, pulled from their length-1 vectors
    int n_header_rows_v = Rf_asInteger(n_header_rows);
    int gene_col_v = Rf_asInteger(gene_col);
    int start_row_v = Rf_asInteger(start_row);

    // copy what is modified in place, so the caller's stays intact
    SEXP expression_vectors_out = PROTECT(Rf_duplicate(expression_vectors)); nprot++;

    // convert what Fortran cannot take from R directly
    char* file_list_c = tox_char_in(file_list, file_list_strlen);
    char* gene_ids_c = tox_char_in(gene_ids, gene_ids_strlen);
    char* delimiter_c = tox_char_in(delimiter, 1);

    // outputs and work space
    int ierr = 0;

    read_expression_vectors_tsv_c(
        file_list_c,
        &file_list_strlen,
        &n_file_list_elements,
        gene_ids_c,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        REAL(expression_vectors_out),
        &n_expression_vectors_elements_dim_1,
        &n_expression_vectors_elements_dim_2,
        &n_header_rows_v,
        &gene_col_v,
        INTEGER(value_cols),
        &n_value_cols_elements,
        &start_row_v,
        &ierr,
        delimiter_c
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, expression_vectors_out);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("expression_vectors"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP read_gene_ids_from_tsv_file_call(SEXP filename, SEXP gene_ids_strlen, SEXP n_gene_ids_elements, SEXP n_header_rows, SEXP gene_col) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);

    // scalar inputs, pulled from their length-1 vectors
    int gene_ids_strlen_v = Rf_asInteger(gene_ids_strlen);
    int n_gene_ids_elements_v = Rf_asInteger(n_gene_ids_elements);
    int n_header_rows_v = Rf_asInteger(n_header_rows);
    int gene_col_v = Rf_asInteger(gene_col);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);

    // outputs and work space
    char* gene_ids_c = tox_char_alloc(gene_ids_strlen_v, n_gene_ids_elements_v);
    int ierr = 0;

    read_gene_ids_from_tsv_file_c(
        filename_c,
        &filename_strlen,
        gene_ids_c,
        &gene_ids_strlen_v,
        &n_gene_ids_elements_v,
        &n_header_rows_v,
        &gene_col_v,
        &ierr
    );

    // convert the outputs back
    SEXP gene_ids = PROTECT(tox_char_out(gene_ids_c, gene_ids_strlen_v, n_gene_ids_elements_v)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, gene_ids);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("gene_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP read_orthofinder_file_call(SEXP filename, SEXP gene_ids, SEXP family_ids_strlen, SEXP n_family_ids_elements, SEXP n_gene_to_fam_elements) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int filename_strlen = tox_max_strlen(filename);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = (int) Rf_length(gene_ids);

    // scalar inputs, pulled from their length-1 vectors
    int family_ids_strlen_v = Rf_asInteger(family_ids_strlen);
    int n_family_ids_elements_v = Rf_asInteger(n_family_ids_elements);
    int n_gene_to_fam_elements_v = Rf_asInteger(n_gene_to_fam_elements);

    // convert what Fortran cannot take from R directly
    char* filename_c = tox_char_in(filename, filename_strlen);
    char* gene_ids_c = tox_char_in(gene_ids, gene_ids_strlen);

    // outputs and work space
    char* family_ids_c = tox_char_alloc(family_ids_strlen_v, n_family_ids_elements_v);
    SEXP gene_to_fam = PROTECT(Rf_allocVector(INTSXP, n_gene_to_fam_elements_v)); nprot++;
    int ierr = 0;

    read_orthofinder_file_c(
        filename_c,
        &filename_strlen,
        gene_ids_c,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        family_ids_c,
        &family_ids_strlen_v,
        &n_family_ids_elements_v,
        INTEGER(gene_to_fam),
        &n_gene_to_fam_elements_v,
        &ierr
    );

    // convert the outputs back
    SEXP family_ids = PROTECT(tox_char_out(family_ids_c, family_ids_strlen_v, n_family_ids_elements_v)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, family_ids);
    SET_VECTOR_ELT(_out, 1, gene_to_fam);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("gene_to_fam"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP get_unassigned_mask_call(SEXP gene_to_fam) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_gene_to_fam_elements = (int) Rf_length(gene_to_fam);

    // outputs and work space
    unsigned char* mask_c = tox_bool_alloc(((int) Rf_length(gene_to_fam)));
    int n_genes_kept = 0;
    int ierr = 0;

    get_unassigned_mask_c(
        INTEGER(gene_to_fam),
        &n_gene_to_fam_elements,
        mask_c,
        &n_genes_kept,
        &ierr
    );

    // convert the outputs back
    SEXP mask = PROTECT(tox_bool_out(mask_c, ((int) Rf_length(gene_to_fam)))); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(n_genes_kept));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("n_genes_kept"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
