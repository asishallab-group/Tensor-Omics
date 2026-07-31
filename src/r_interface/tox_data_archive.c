// Generated. Do not edit.
#if !defined(NO_R_INTERFACE) && !defined(NO_C_INTERFACE)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void create_zip_archive_c(const char*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, int*);
void save_tox_data_c(const char*, const int*, int*, const char*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, const int*, const int*, const char*, const int*, const char*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*);
void get_tox_data_dims_c(const char*, const int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
void read_tox_data_into_c(const char*, const int*, const int*, const int*, char*, const int*, const int*, double*, const int*, int*, const int*, const int*, char*, const int*, const int*, double*, const int*, const int*, double*, int*);

SEXP create_zip_archive_call(SEXP zip_filename, SEXP keys, SEXP filenames) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);
    int keys_strlen = tox_max_strlen(keys);
    int n_keys_elements = (int) Rf_length(keys);
    int filenames_strlen = tox_max_strlen(filenames);
    int n_filenames_elements = (int) Rf_length(filenames);

    // convert what Fortran cannot take from R directly
    char* zip_filename_c = tox_char_in(zip_filename, zip_filename_strlen);
    char* keys_c = tox_char_in(keys, keys_strlen);
    char* filenames_c = tox_char_in(filenames, filenames_strlen);

    // outputs and work space
    int ierr = 0;

    create_zip_archive_c(
        zip_filename_c,
        &zip_filename_strlen,
        keys_c,
        &keys_strlen,
        &n_keys_elements,
        filenames_c,
        &filenames_strlen,
        &n_filenames_elements,
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

SEXP save_tox_data_call(SEXP zip_filename, SEXP gene_ids, SEXP gene_ids_file, SEXP expression, SEXP expression_file, SEXP gene_to_family, SEXP gene_to_family_file, SEXP family_ids, SEXP family_ids_file, SEXP family_centroids, SEXP family_centroids_file, SEXP shift_vectors, SEXP shift_vectors_file) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int gene_ids_size = 0;
    if (gene_ids != R_NilValue) {
        gene_ids_size = (int) Rf_length(gene_ids);
    }
    int gene_ids_file_size = 0;
    if (gene_ids_file != R_NilValue) {
        gene_ids_file_size = (int) Rf_length(gene_ids_file);
    }
    const double* expression_p = NULL;
    int expression_size = 0;
    int expression_dim[2]; for (int i = 0; i < 2; ++i) expression_dim[i] = 0;
    if (expression != R_NilValue) {
        expression_size = (int) Rf_length(expression);
        expression_p = REAL(expression);
        SEXP expression_d = Rf_getAttrib(expression, R_DimSymbol);
        if (expression_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(expression_d); ++i) expression_dim[i] = INTEGER(expression_d)[i];
        }
    }
    int expression_file_size = 0;
    if (expression_file != R_NilValue) {
        expression_file_size = (int) Rf_length(expression_file);
    }
    const int* gene_to_family_p = NULL;
    int gene_to_family_size = 0;
    if (gene_to_family != R_NilValue) {
        gene_to_family_size = (int) Rf_length(gene_to_family);
        gene_to_family_p = INTEGER(gene_to_family);
    }
    int gene_to_family_file_size = 0;
    if (gene_to_family_file != R_NilValue) {
        gene_to_family_file_size = (int) Rf_length(gene_to_family_file);
    }
    int family_ids_size = 0;
    if (family_ids != R_NilValue) {
        family_ids_size = (int) Rf_length(family_ids);
    }
    int family_ids_file_size = 0;
    if (family_ids_file != R_NilValue) {
        family_ids_file_size = (int) Rf_length(family_ids_file);
    }
    const double* family_centroids_p = NULL;
    int family_centroids_size = 0;
    int family_centroids_dim[2]; for (int i = 0; i < 2; ++i) family_centroids_dim[i] = 0;
    if (family_centroids != R_NilValue) {
        family_centroids_size = (int) Rf_length(family_centroids);
        family_centroids_p = REAL(family_centroids);
        SEXP family_centroids_d = Rf_getAttrib(family_centroids, R_DimSymbol);
        if (family_centroids_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(family_centroids_d); ++i) family_centroids_dim[i] = INTEGER(family_centroids_d)[i];
        }
    }
    int family_centroids_file_size = 0;
    if (family_centroids_file != R_NilValue) {
        family_centroids_file_size = (int) Rf_length(family_centroids_file);
    }
    const double* shift_vectors_p = NULL;
    int shift_vectors_size = 0;
    int shift_vectors_dim[2]; for (int i = 0; i < 2; ++i) shift_vectors_dim[i] = 0;
    if (shift_vectors != R_NilValue) {
        shift_vectors_size = (int) Rf_length(shift_vectors);
        shift_vectors_p = REAL(shift_vectors);
        SEXP shift_vectors_d = Rf_getAttrib(shift_vectors, R_DimSymbol);
        if (shift_vectors_d != R_NilValue) {
            for (int i = 0; i < 2 && i < (int) Rf_length(shift_vectors_d); ++i) shift_vectors_dim[i] = INTEGER(shift_vectors_d)[i];
        }
    }
    int shift_vectors_file_size = 0;
    if (shift_vectors_file != R_NilValue) {
        shift_vectors_file_size = (int) Rf_length(shift_vectors_file);
    }

    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = gene_ids_size;
    int gene_ids_file_strlen = tox_max_strlen(gene_ids_file);
    int n_expression_elements_dim_1 = expression_dim[0];
    int n_expression_elements_dim_2 = expression_dim[1];
    int expression_file_strlen = tox_max_strlen(expression_file);
    int n_gene_to_family_elements = gene_to_family_size;
    int gene_to_family_file_strlen = tox_max_strlen(gene_to_family_file);
    int family_ids_strlen = tox_max_strlen(family_ids);
    int n_family_ids_elements = family_ids_size;
    int family_ids_file_strlen = tox_max_strlen(family_ids_file);
    int n_family_centroids_elements_dim_1 = family_centroids_dim[0];
    int n_family_centroids_elements_dim_2 = family_centroids_dim[1];
    int family_centroids_file_strlen = tox_max_strlen(family_centroids_file);
    int n_shift_vectors_elements_dim_1 = shift_vectors_dim[0];
    int n_shift_vectors_elements_dim_2 = shift_vectors_dim[1];
    int shift_vectors_file_strlen = tox_max_strlen(shift_vectors_file);

    // convert what Fortran cannot take from R directly
    char* zip_filename_c = tox_char_in(zip_filename, zip_filename_strlen);
    char* gene_ids_c = tox_char_in(gene_ids, gene_ids_strlen);
    char* gene_ids_file_c = tox_char_in(gene_ids_file, gene_ids_file_strlen);
    char* expression_file_c = tox_char_in(expression_file, expression_file_strlen);
    char* gene_to_family_file_c = tox_char_in(gene_to_family_file, gene_to_family_file_strlen);
    char* family_ids_c = tox_char_in(family_ids, family_ids_strlen);
    char* family_ids_file_c = tox_char_in(family_ids_file, family_ids_file_strlen);
    char* family_centroids_file_c = tox_char_in(family_centroids_file, family_centroids_file_strlen);
    char* shift_vectors_file_c = tox_char_in(shift_vectors_file, shift_vectors_file_strlen);

    // outputs and work space
    int ierr = 0;

    save_tox_data_c(
        zip_filename_c,
        &zip_filename_strlen,
        &ierr,
        gene_ids != R_NilValue ? gene_ids_c : NULL,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        gene_ids_file != R_NilValue ? gene_ids_file_c : NULL,
        &gene_ids_file_strlen,
        expression_p,
        &n_expression_elements_dim_1,
        &n_expression_elements_dim_2,
        expression_file != R_NilValue ? expression_file_c : NULL,
        &expression_file_strlen,
        gene_to_family_p,
        &n_gene_to_family_elements,
        gene_to_family_file != R_NilValue ? gene_to_family_file_c : NULL,
        &gene_to_family_file_strlen,
        family_ids != R_NilValue ? family_ids_c : NULL,
        &family_ids_strlen,
        &n_family_ids_elements,
        family_ids_file != R_NilValue ? family_ids_file_c : NULL,
        &family_ids_file_strlen,
        family_centroids_p,
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        family_centroids_file != R_NilValue ? family_centroids_file_c : NULL,
        &family_centroids_file_strlen,
        shift_vectors_p,
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        shift_vectors_file != R_NilValue ? shift_vectors_file_c : NULL,
        &shift_vectors_file_strlen
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 1)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 1)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP get_tox_data_dims_call(SEXP zip_filename) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);

    // convert what Fortran cannot take from R directly
    char* zip_filename_c = tox_char_in(zip_filename, zip_filename_strlen);

    // outputs and work space
    int n_gene_ids = 0;
    int gene_id_len = 0;
    int n_expression_rows = 0;
    int n_expression_cols = 0;
    int n_gene_to_family = 0;
    int n_family_ids = 0;
    int family_id_len = 0;
    int n_family_centroids_rows = 0;
    int n_family_centroids_cols = 0;
    int n_shift_vectors_rows = 0;
    int n_shift_vectors_cols = 0;
    int ierr = 0;

    get_tox_data_dims_c(
        zip_filename_c,
        &zip_filename_strlen,
        &n_gene_ids,
        &gene_id_len,
        &n_expression_rows,
        &n_expression_cols,
        &n_gene_to_family,
        &n_family_ids,
        &family_id_len,
        &n_family_centroids_rows,
        &n_family_centroids_cols,
        &n_shift_vectors_rows,
        &n_shift_vectors_cols,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 12)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_gene_ids));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(gene_id_len));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(n_expression_rows));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(n_expression_cols));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(n_gene_to_family));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(n_family_ids));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(family_id_len));
    SET_VECTOR_ELT(_out, 7, Rf_ScalarInteger(n_family_centroids_rows));
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(n_family_centroids_cols));
    SET_VECTOR_ELT(_out, 9, Rf_ScalarInteger(n_shift_vectors_rows));
    SET_VECTOR_ELT(_out, 10, Rf_ScalarInteger(n_shift_vectors_cols));
    SET_VECTOR_ELT(_out, 11, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 12)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_gene_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("gene_id_len"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("n_expression_rows"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("n_expression_cols"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("n_gene_to_family"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("n_family_ids"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("family_id_len"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("n_family_centroids_rows"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("n_family_centroids_cols"));
    SET_STRING_ELT(_nms, 9, Rf_mkChar("n_shift_vectors_rows"));
    SET_STRING_ELT(_nms, 10, Rf_mkChar("n_shift_vectors_cols"));
    SET_STRING_ELT(_nms, 11, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP read_tox_data_into_call(SEXP zip_filename, SEXP n_gene_ids, SEXP gene_id_len, SEXP n_expression_rows, SEXP n_expression_cols, SEXP n_gene_to_family, SEXP n_family_ids, SEXP family_id_len, SEXP n_family_centroids_rows, SEXP n_family_centroids_cols, SEXP n_shift_vectors_rows, SEXP n_shift_vectors_cols) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);

    // scalar inputs, pulled from their length-1 vectors
    int n_gene_ids_v = Rf_asInteger(n_gene_ids);
    int gene_id_len_v = Rf_asInteger(gene_id_len);
    int n_expression_rows_v = Rf_asInteger(n_expression_rows);
    int n_expression_cols_v = Rf_asInteger(n_expression_cols);
    int n_gene_to_family_v = Rf_asInteger(n_gene_to_family);
    int n_family_ids_v = Rf_asInteger(n_family_ids);
    int family_id_len_v = Rf_asInteger(family_id_len);
    int n_family_centroids_rows_v = Rf_asInteger(n_family_centroids_rows);
    int n_family_centroids_cols_v = Rf_asInteger(n_family_centroids_cols);
    int n_shift_vectors_rows_v = Rf_asInteger(n_shift_vectors_rows);
    int n_shift_vectors_cols_v = Rf_asInteger(n_shift_vectors_cols);

    // convert what Fortran cannot take from R directly
    char* zip_filename_c = tox_char_in(zip_filename, zip_filename_strlen);

    // outputs and work space
    char* gene_ids_c = tox_char_alloc(gene_id_len_v, n_gene_ids_v);
    SEXP expression = PROTECT(Rf_allocVector(REALSXP, n_expression_rows_v * n_expression_cols_v)); nprot++;
    { SEXP expression_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(expression_dim)[0] = n_expression_rows_v; INTEGER(expression_dim)[1] = n_expression_cols_v; Rf_setAttrib(expression, R_DimSymbol, expression_dim); UNPROTECT(1); }
    SEXP gene_to_family = PROTECT(Rf_allocVector(INTSXP, n_gene_to_family_v)); nprot++;
    char* family_ids_c = tox_char_alloc(family_id_len_v, n_family_ids_v);
    SEXP family_centroids = PROTECT(Rf_allocVector(REALSXP, n_family_centroids_rows_v * n_family_centroids_cols_v)); nprot++;
    { SEXP family_centroids_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(family_centroids_dim)[0] = n_family_centroids_rows_v; INTEGER(family_centroids_dim)[1] = n_family_centroids_cols_v; Rf_setAttrib(family_centroids, R_DimSymbol, family_centroids_dim); UNPROTECT(1); }
    SEXP shift_vectors = PROTECT(Rf_allocVector(REALSXP, n_shift_vectors_rows_v * n_shift_vectors_cols_v)); nprot++;
    { SEXP shift_vectors_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(shift_vectors_dim)[0] = n_shift_vectors_rows_v; INTEGER(shift_vectors_dim)[1] = n_shift_vectors_cols_v; Rf_setAttrib(shift_vectors, R_DimSymbol, shift_vectors_dim); UNPROTECT(1); }
    int ierr = 0;

    read_tox_data_into_c(
        zip_filename_c,
        &zip_filename_strlen,
        &n_gene_ids_v,
        &gene_id_len_v,
        gene_ids_c,
        &n_expression_rows_v,
        &n_expression_cols_v,
        REAL(expression),
        &n_gene_to_family_v,
        INTEGER(gene_to_family),
        &n_family_ids_v,
        &family_id_len_v,
        family_ids_c,
        &n_family_centroids_rows_v,
        &n_family_centroids_cols_v,
        REAL(family_centroids),
        &n_shift_vectors_rows_v,
        &n_shift_vectors_cols_v,
        REAL(shift_vectors),
        &ierr
    );

    // convert the outputs back
    SEXP gene_ids = PROTECT(tox_char_out(gene_ids_c, gene_id_len_v, n_gene_ids_v)); nprot++;
    SEXP family_ids = PROTECT(tox_char_out(family_ids_c, family_id_len_v, n_family_ids_v)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); nprot++;
    SET_VECTOR_ELT(_out, 0, gene_ids);
    SET_VECTOR_ELT(_out, 1, expression);
    SET_VECTOR_ELT(_out, 2, gene_to_family);
    SET_VECTOR_ELT(_out, 3, family_ids);
    SET_VECTOR_ELT(_out, 4, family_centroids);
    SET_VECTOR_ELT(_out, 5, shift_vectors);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("gene_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("expression"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("gene_to_family"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("family_ids"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("family_centroids"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("shift_vectors"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R interface
