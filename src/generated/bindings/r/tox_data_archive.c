// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void create_zip_archive_c(const char*, const int*, const char*, const int*, const int*, const char*, const int*, const int*, int*);
void save_tox_data_c(const char*, const int*, const char*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, const int*, const int*, const char*, const int*, const char*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, const double*, const int*, const int*, const char*, const int*, int*);
void get_tox_data_dims_c(const char*, const int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*, int*);
void read_tox_data_into_c(const char*, const int*, const int*, const int*, char*, const int*, const int*, double*, const int*, int*, const int*, const int*, char*, const int*, const int*, double*, const int*, const int*, double*, int*);

SEXP create_zip_archive_call(SEXP zip_filename, SEXP keys, SEXP filenames) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);
    int keys_strlen = tox_max_strlen(keys);
    int n_keys_elements = (int) Rf_length(keys);
    int filenames_strlen = tox_max_strlen(filenames);
    int n_filenames_elements = (int) Rf_length(filenames);

    // convert what Fortran cannot take from R directly
    char* _c_zip_filename = tox_char_in(zip_filename, zip_filename_strlen, "zip_filename");
    char* _c_keys = tox_char_in(keys, keys_strlen, "keys");
    char* _c_filenames = tox_char_in(filenames, filenames_strlen, "filenames");

    // outputs and work space
    int ierr = 0;

    create_zip_archive_c(
        _c_zip_filename,
        &zip_filename_strlen,
        _c_keys,
        &keys_strlen,
        &n_keys_elements,
        _c_filenames,
        &filenames_strlen,
        &n_filenames_elements,
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

SEXP save_tox_data_call(SEXP zip_filename, SEXP gene_ids, SEXP gene_ids_file, SEXP expression, SEXP expression_file, SEXP gene_to_family, SEXP gene_to_family_file, SEXP family_ids, SEXP family_ids_file, SEXP family_centroids, SEXP family_centroids_file, SEXP shift_vectors, SEXP shift_vectors_file) {
    int _nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    int _size_gene_ids = 0;
    if (gene_ids != R_NilValue) {
        _size_gene_ids = (int) Rf_length(gene_ids);
    }
    int _size_gene_ids_file = 0;
    if (gene_ids_file != R_NilValue) {
        _size_gene_ids_file = (int) Rf_length(gene_ids_file);
    }
    const double* _p_expression = NULL;
    int _size_expression = 0;
    int _dim_expression[2]; for (int _i = 0; _i < 2; ++_i) _dim_expression[_i] = 0;
    if (expression != R_NilValue) {
        _size_expression = (int) Rf_length(expression);
        _p_expression = REAL(expression);
        SEXP _d_expression = Rf_getAttrib(expression, R_DimSymbol);
        if (_d_expression != R_NilValue) {
            for (int _i = 0; _i < 2 && _i < (int) Rf_length(_d_expression); ++_i) _dim_expression[_i] = INTEGER(_d_expression)[_i];
        }
    }
    int _size_expression_file = 0;
    if (expression_file != R_NilValue) {
        _size_expression_file = (int) Rf_length(expression_file);
    }
    int _size_gene_to_family = 0;
    if (gene_to_family != R_NilValue) {
        _size_gene_to_family = (int) Rf_length(gene_to_family);
    }
    int _size_gene_to_family_file = 0;
    if (gene_to_family_file != R_NilValue) {
        _size_gene_to_family_file = (int) Rf_length(gene_to_family_file);
    }
    int _size_family_ids = 0;
    if (family_ids != R_NilValue) {
        _size_family_ids = (int) Rf_length(family_ids);
    }
    int _size_family_ids_file = 0;
    if (family_ids_file != R_NilValue) {
        _size_family_ids_file = (int) Rf_length(family_ids_file);
    }
    const double* _p_family_centroids = NULL;
    int _size_family_centroids = 0;
    int _dim_family_centroids[2]; for (int _i = 0; _i < 2; ++_i) _dim_family_centroids[_i] = 0;
    if (family_centroids != R_NilValue) {
        _size_family_centroids = (int) Rf_length(family_centroids);
        _p_family_centroids = REAL(family_centroids);
        SEXP _d_family_centroids = Rf_getAttrib(family_centroids, R_DimSymbol);
        if (_d_family_centroids != R_NilValue) {
            for (int _i = 0; _i < 2 && _i < (int) Rf_length(_d_family_centroids); ++_i) _dim_family_centroids[_i] = INTEGER(_d_family_centroids)[_i];
        }
    }
    int _size_family_centroids_file = 0;
    if (family_centroids_file != R_NilValue) {
        _size_family_centroids_file = (int) Rf_length(family_centroids_file);
    }
    const double* _p_shift_vectors = NULL;
    int _size_shift_vectors = 0;
    int _dim_shift_vectors[2]; for (int _i = 0; _i < 2; ++_i) _dim_shift_vectors[_i] = 0;
    if (shift_vectors != R_NilValue) {
        _size_shift_vectors = (int) Rf_length(shift_vectors);
        _p_shift_vectors = REAL(shift_vectors);
        SEXP _d_shift_vectors = Rf_getAttrib(shift_vectors, R_DimSymbol);
        if (_d_shift_vectors != R_NilValue) {
            for (int _i = 0; _i < 2 && _i < (int) Rf_length(_d_shift_vectors); ++_i) _dim_shift_vectors[_i] = INTEGER(_d_shift_vectors)[_i];
        }
    }
    int _size_shift_vectors_file = 0;
    if (shift_vectors_file != R_NilValue) {
        _size_shift_vectors_file = (int) Rf_length(shift_vectors_file);
    }

    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);
    int gene_ids_strlen = tox_max_strlen(gene_ids);
    int n_gene_ids_elements = _size_gene_ids;
    int gene_ids_file_strlen = tox_max_strlen(gene_ids_file);
    int n_expression_elements_dim_1 = _dim_expression[0];
    int n_expression_elements_dim_2 = _dim_expression[1];
    int expression_file_strlen = tox_max_strlen(expression_file);
    int n_gene_to_family_elements = _size_gene_to_family;
    int gene_to_family_file_strlen = tox_max_strlen(gene_to_family_file);
    int family_ids_strlen = tox_max_strlen(family_ids);
    int n_family_ids_elements = _size_family_ids;
    int family_ids_file_strlen = tox_max_strlen(family_ids_file);
    int n_family_centroids_elements_dim_1 = _dim_family_centroids[0];
    int n_family_centroids_elements_dim_2 = _dim_family_centroids[1];
    int family_centroids_file_strlen = tox_max_strlen(family_centroids_file);
    int n_shift_vectors_elements_dim_1 = _dim_shift_vectors[0];
    int n_shift_vectors_elements_dim_2 = _dim_shift_vectors[1];
    int shift_vectors_file_strlen = tox_max_strlen(shift_vectors_file);

    // convert what Fortran cannot take from R directly
    char* _c_zip_filename = tox_char_in(zip_filename, zip_filename_strlen, "zip_filename");
    char* _c_gene_ids = tox_char_in(gene_ids, gene_ids_strlen, "gene_ids");
    char* _c_gene_ids_file = tox_char_in(gene_ids_file, gene_ids_file_strlen, "gene_ids_file");
    char* _c_expression_file = tox_char_in(expression_file, expression_file_strlen, "expression_file");
    const int* _c_gene_to_family = tox_int_in(gene_to_family, "gene_to_family");
    char* _c_gene_to_family_file = tox_char_in(gene_to_family_file, gene_to_family_file_strlen, "gene_to_family_file");
    char* _c_family_ids = tox_char_in(family_ids, family_ids_strlen, "family_ids");
    char* _c_family_ids_file = tox_char_in(family_ids_file, family_ids_file_strlen, "family_ids_file");
    char* _c_family_centroids_file = tox_char_in(family_centroids_file, family_centroids_file_strlen, "family_centroids_file");
    char* _c_shift_vectors_file = tox_char_in(shift_vectors_file, shift_vectors_file_strlen, "shift_vectors_file");

    // outputs and work space
    int ierr = 0;

    save_tox_data_c(
        _c_zip_filename,
        &zip_filename_strlen,
        gene_ids != R_NilValue ? _c_gene_ids : NULL,
        &gene_ids_strlen,
        &n_gene_ids_elements,
        gene_ids_file != R_NilValue ? _c_gene_ids_file : NULL,
        &gene_ids_file_strlen,
        _p_expression,
        &n_expression_elements_dim_1,
        &n_expression_elements_dim_2,
        expression_file != R_NilValue ? _c_expression_file : NULL,
        &expression_file_strlen,
        gene_to_family != R_NilValue ? _c_gene_to_family : NULL,
        &n_gene_to_family_elements,
        gene_to_family_file != R_NilValue ? _c_gene_to_family_file : NULL,
        &gene_to_family_file_strlen,
        family_ids != R_NilValue ? _c_family_ids : NULL,
        &family_ids_strlen,
        &n_family_ids_elements,
        family_ids_file != R_NilValue ? _c_family_ids_file : NULL,
        &family_ids_file_strlen,
        _p_family_centroids,
        &n_family_centroids_elements_dim_1,
        &n_family_centroids_elements_dim_2,
        family_centroids_file != R_NilValue ? _c_family_centroids_file : NULL,
        &family_centroids_file_strlen,
        _p_shift_vectors,
        &n_shift_vectors_elements_dim_1,
        &n_shift_vectors_elements_dim_2,
        shift_vectors_file != R_NilValue ? _c_shift_vectors_file : NULL,
        &shift_vectors_file_strlen,
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

SEXP get_tox_data_dims_call(SEXP zip_filename) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);

    // convert what Fortran cannot take from R directly
    char* _c_zip_filename = tox_char_in(zip_filename, zip_filename_strlen, "zip_filename");

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
        _c_zip_filename,
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

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 12)); _nprot++;
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
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 12)); _nprot++;
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
    UNPROTECT(_nprot);
    return _out;
}

SEXP read_tox_data_into_call(SEXP zip_filename, SEXP n_gene_ids, SEXP gene_id_len, SEXP n_expression_rows, SEXP n_expression_cols, SEXP n_gene_to_family, SEXP n_family_ids, SEXP family_id_len, SEXP n_family_centroids_rows, SEXP n_family_centroids_cols, SEXP n_shift_vectors_rows, SEXP n_shift_vectors_cols) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int zip_filename_strlen = tox_max_strlen(zip_filename);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_gene_ids = Rf_asInteger(n_gene_ids);
    int _v_gene_id_len = Rf_asInteger(gene_id_len);
    int _v_n_expression_rows = Rf_asInteger(n_expression_rows);
    int _v_n_expression_cols = Rf_asInteger(n_expression_cols);
    int _v_n_gene_to_family = Rf_asInteger(n_gene_to_family);
    int _v_n_family_ids = Rf_asInteger(n_family_ids);
    int _v_family_id_len = Rf_asInteger(family_id_len);
    int _v_n_family_centroids_rows = Rf_asInteger(n_family_centroids_rows);
    int _v_n_family_centroids_cols = Rf_asInteger(n_family_centroids_cols);
    int _v_n_shift_vectors_rows = Rf_asInteger(n_shift_vectors_rows);
    int _v_n_shift_vectors_cols = Rf_asInteger(n_shift_vectors_cols);

    // convert what Fortran cannot take from R directly
    char* _c_zip_filename = tox_char_in(zip_filename, zip_filename_strlen, "zip_filename");

    // outputs and work space
    char* _c_gene_ids = tox_char_alloc(_v_gene_id_len, _v_n_gene_ids);
    SEXP expression = PROTECT(Rf_allocVector(REALSXP, _v_n_expression_rows * _v_n_expression_cols)); _nprot++;
    { SEXP _dim_expression = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_expression)[0] = _v_n_expression_rows; INTEGER(_dim_expression)[1] = _v_n_expression_cols; Rf_setAttrib(expression, R_DimSymbol, _dim_expression); UNPROTECT(1); }
    SEXP gene_to_family = PROTECT(Rf_allocVector(INTSXP, _v_n_gene_to_family)); _nprot++;
    char* _c_family_ids = tox_char_alloc(_v_family_id_len, _v_n_family_ids);
    SEXP family_centroids = PROTECT(Rf_allocVector(REALSXP, _v_n_family_centroids_rows * _v_n_family_centroids_cols)); _nprot++;
    { SEXP _dim_family_centroids = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_family_centroids)[0] = _v_n_family_centroids_rows; INTEGER(_dim_family_centroids)[1] = _v_n_family_centroids_cols; Rf_setAttrib(family_centroids, R_DimSymbol, _dim_family_centroids); UNPROTECT(1); }
    SEXP shift_vectors = PROTECT(Rf_allocVector(REALSXP, _v_n_shift_vectors_rows * _v_n_shift_vectors_cols)); _nprot++;
    { SEXP _dim_shift_vectors = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_shift_vectors)[0] = _v_n_shift_vectors_rows; INTEGER(_dim_shift_vectors)[1] = _v_n_shift_vectors_cols; Rf_setAttrib(shift_vectors, R_DimSymbol, _dim_shift_vectors); UNPROTECT(1); }
    int ierr = 0;

    read_tox_data_into_c(
        _c_zip_filename,
        &zip_filename_strlen,
        &_v_n_gene_ids,
        &_v_gene_id_len,
        _c_gene_ids,
        &_v_n_expression_rows,
        &_v_n_expression_cols,
        REAL(expression),
        &_v_n_gene_to_family,
        INTEGER(gene_to_family),
        &_v_n_family_ids,
        &_v_family_id_len,
        _c_family_ids,
        &_v_n_family_centroids_rows,
        &_v_n_family_centroids_cols,
        REAL(family_centroids),
        &_v_n_shift_vectors_rows,
        &_v_n_shift_vectors_cols,
        REAL(shift_vectors),
        &ierr
    );

    // convert the outputs back
    SEXP gene_ids = PROTECT(tox_char_out(_c_gene_ids, _v_gene_id_len, _v_n_gene_ids)); _nprot++;
    SEXP family_ids = PROTECT(tox_char_out(_c_family_ids, _v_family_id_len, _v_n_family_ids)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); _nprot++;
    SET_VECTOR_ELT(_out, 0, gene_ids);
    SET_VECTOR_ELT(_out, 1, expression);
    SET_VECTOR_ELT(_out, 2, gene_to_family);
    SET_VECTOR_ELT(_out, 3, family_ids);
    SET_VECTOR_ELT(_out, 4, family_centroids);
    SET_VECTOR_ELT(_out, 5, shift_vectors);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("gene_ids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("expression"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("gene_to_family"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("family_ids"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("family_centroids"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("shift_vectors"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
