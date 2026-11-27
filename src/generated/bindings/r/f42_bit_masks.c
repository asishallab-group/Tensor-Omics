// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h f7ca108707d9ae14 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void bitmask_n_words_c(const int*, int*, int*);
void bitmask_test_c(const int*, const int*, const int*, const int*, unsigned char*, int*);
void bitmask_from_logical_c(const int*, const int*, const unsigned char*, int*, int*);
void bitmask_to_logical_c(const int*, const int*, const int*, unsigned char*, int*);
void bitmasks_from_logical_2D_c(const int*, const int*, const int*, const unsigned char*, int*, int*);
void bitmasks_to_logical_2D_c(const int*, const int*, const int*, const int*, unsigned char*, int*);

SEXP bitmask_n_words_call(SEXP n_bits) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);

    // outputs and work space
    int n_words = 0;
    int ierr = 0;

    bitmask_n_words_c(
        &_v_n_bits,
        &n_words,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_words));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_words"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bitmask_test_call(SEXP n_bits, SEXP bitmask, SEXP i_bit) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_words = (int) Rf_length(bitmask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);
    int _v_i_bit = Rf_asInteger(i_bit);

    // convert what Fortran cannot take from R directly
    const int* _c_bitmask = tox_words_in(bitmask, "bitmask");

    // outputs and work space
    unsigned char is_set = 0;
    int ierr = 0;

    bitmask_test_c(
        &_v_n_bits,
        &n_words,
        _c_bitmask,
        &_v_i_bit,
        &is_set,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(is_set != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_set"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bitmask_from_logical_call(SEXP n_words, SEXP flags) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bits = (int) Rf_length(flags);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_words = Rf_asInteger(n_words);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_flags = tox_bool_in(flags, "flags");

    // outputs and work space
    R_xlen_t _n_bitmask = tox_count("bitmask", (double) _v_n_words);
    SEXP bitmask = PROTECT(Rf_allocVector(INTSXP, _n_bitmask)); _nprot++;
    int ierr = 0;

    bitmask_from_logical_c(
        &n_bits,
        &_v_n_words,
        _c_flags,
        INTEGER(bitmask),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, bitmask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("bitmask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bitmask_to_logical_call(SEXP n_bits, SEXP bitmask) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_words = (int) Rf_length(bitmask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);

    // convert what Fortran cannot take from R directly
    const int* _c_bitmask = tox_words_in(bitmask, "bitmask");

    // outputs and work space
    R_xlen_t _n_flags = tox_count("flags", (double) _v_n_bits);
    unsigned char* _c_flags = tox_bool_alloc(_n_flags);
    int ierr = 0;

    bitmask_to_logical_c(
        &_v_n_bits,
        &n_words,
        _c_bitmask,
        _c_flags,
        &ierr
    );

    // convert the outputs back, unless the library reports an error
    SEXP flags = R_NilValue;
    if (!tox_failed(ierr)) {
        flags = PROTECT(tox_bool_out(_c_flags, _n_flags)); _nprot++;
    }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, flags);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("flags"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bitmasks_from_logical_2D_call(SEXP n_words, SEXP flags) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bits = INTEGER(Rf_getAttrib(flags, R_DimSymbol))[0];
    int n_masks = INTEGER(Rf_getAttrib(flags, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_words = Rf_asInteger(n_words);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_flags = tox_bool_in(flags, "flags");

    // outputs and work space
    R_xlen_t _n_bitmasks = tox_count("bitmasks", (double) _v_n_words * (double) n_masks);
    SEXP bitmasks = PROTECT(Rf_allocVector(INTSXP, _n_bitmasks)); _nprot++;
    { SEXP _dim_bitmasks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_bitmasks)[0] = _v_n_words; INTEGER(_dim_bitmasks)[1] = n_masks; Rf_setAttrib(bitmasks, R_DimSymbol, _dim_bitmasks); UNPROTECT(1); }
    int ierr = 0;

    bitmasks_from_logical_2D_c(
        &n_bits,
        &n_masks,
        &_v_n_words,
        _c_flags,
        INTEGER(bitmasks),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, bitmasks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("bitmasks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bitmasks_to_logical_2D_call(SEXP n_bits, SEXP bitmasks) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_masks = INTEGER(Rf_getAttrib(bitmasks, R_DimSymbol))[1];
    int n_words = INTEGER(Rf_getAttrib(bitmasks, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);

    // convert what Fortran cannot take from R directly
    const int* _c_bitmasks = tox_words_in(bitmasks, "bitmasks");

    // outputs and work space
    R_xlen_t _n_flags = tox_count("flags", (double) _v_n_bits * (double) n_masks);
    unsigned char* _c_flags = tox_bool_alloc(_n_flags);
    int ierr = 0;

    bitmasks_to_logical_2D_c(
        &_v_n_bits,
        &n_masks,
        &n_words,
        _c_bitmasks,
        _c_flags,
        &ierr
    );

    // convert the outputs back, unless the library reports an error
    SEXP flags = R_NilValue;
    if (!tox_failed(ierr)) {
        flags = PROTECT(tox_bool_out(_c_flags, _n_flags)); _nprot++;
        { SEXP _dim_flags = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_flags)[0] = _v_n_bits; INTEGER(_dim_flags)[1] = n_masks; Rf_setAttrib(flags, R_DimSymbol, _dim_flags); UNPROTECT(1); }
    }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, flags);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("flags"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
