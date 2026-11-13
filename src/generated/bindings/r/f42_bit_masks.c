// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e5d25b7db337126 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void bit_mask_n_words_c(const int*, int*, int*);
void bit_mask_test_c(const int*, const int*, const int*, const int*, unsigned char*, int*);
void bit_mask_from_logical_c(const int*, const int*, const unsigned char*, int*, int*);
void bit_mask_to_logical_c(const int*, const int*, const int*, unsigned char*, int*);
void bit_masks_from_logical_2D_c(const int*, const int*, const int*, const unsigned char*, int*, int*);
void bit_masks_to_logical_2D_c(const int*, const int*, const int*, const int*, unsigned char*, int*);

SEXP bit_mask_n_words_call(SEXP n_bits) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);

    // outputs and work space
    int n_words = 0;
    int ierr = 0;

    bit_mask_n_words_c(
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

SEXP bit_mask_test_call(SEXP n_bits, SEXP bit_mask, SEXP i_bit) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_words = (int) Rf_length(bit_mask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);
    int _v_i_bit = Rf_asInteger(i_bit);

    // convert what Fortran cannot take from R directly
    const int* _c_bit_mask = tox_int_in(bit_mask, "bit_mask");

    // outputs and work space
    unsigned char is_set = 0;
    int ierr = 0;

    bit_mask_test_c(
        &_v_n_bits,
        &n_words,
        _c_bit_mask,
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

SEXP bit_mask_from_logical_call(SEXP n_words, SEXP flags) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bits = (int) Rf_length(flags);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_words = Rf_asInteger(n_words);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_flags = tox_bool_in(flags, "flags");

    // outputs and work space
    SEXP bit_mask = PROTECT(Rf_allocVector(INTSXP, _v_n_words)); _nprot++;
    int ierr = 0;

    bit_mask_from_logical_c(
        &n_bits,
        &_v_n_words,
        _c_flags,
        INTEGER(bit_mask),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, bit_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("bit_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bit_mask_to_logical_call(SEXP n_bits, SEXP bit_mask) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_words = (int) Rf_length(bit_mask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);

    // convert what Fortran cannot take from R directly
    const int* _c_bit_mask = tox_int_in(bit_mask, "bit_mask");

    // outputs and work space
    unsigned char* _c_flags = tox_bool_alloc(_v_n_bits);
    int ierr = 0;

    bit_mask_to_logical_c(
        &_v_n_bits,
        &n_words,
        _c_bit_mask,
        _c_flags,
        &ierr
    );

    // convert the outputs back
    SEXP flags = PROTECT(tox_bool_out(_c_flags, _v_n_bits)); _nprot++;

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

SEXP bit_masks_from_logical_2D_call(SEXP n_words, SEXP flags) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bits = INTEGER(Rf_getAttrib(flags, R_DimSymbol))[0];
    int n_masks = INTEGER(Rf_getAttrib(flags, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_words = Rf_asInteger(n_words);

    // convert what Fortran cannot take from R directly
    unsigned char* _c_flags = tox_bool_in(flags, "flags");

    // outputs and work space
    SEXP bit_masks = PROTECT(Rf_allocVector(INTSXP, _v_n_words * n_masks)); _nprot++;
    { SEXP _dim_bit_masks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_bit_masks)[0] = _v_n_words; INTEGER(_dim_bit_masks)[1] = n_masks; Rf_setAttrib(bit_masks, R_DimSymbol, _dim_bit_masks); UNPROTECT(1); }
    int ierr = 0;

    bit_masks_from_logical_2D_c(
        &n_bits,
        &n_masks,
        &_v_n_words,
        _c_flags,
        INTEGER(bit_masks),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, bit_masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("bit_masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP bit_masks_to_logical_2D_call(SEXP n_bits, SEXP bit_masks) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_masks = INTEGER(Rf_getAttrib(bit_masks, R_DimSymbol))[1];
    int n_words = INTEGER(Rf_getAttrib(bit_masks, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_bits = Rf_asInteger(n_bits);

    // convert what Fortran cannot take from R directly
    const int* _c_bit_masks = tox_int_in(bit_masks, "bit_masks");

    // outputs and work space
    unsigned char* _c_flags = tox_bool_alloc(_v_n_bits * n_masks);
    int ierr = 0;

    bit_masks_to_logical_2D_c(
        &_v_n_bits,
        &n_masks,
        &n_words,
        _c_bit_masks,
        _c_flags,
        &ierr
    );

    // convert the outputs back
    SEXP flags = PROTECT(tox_bool_out(_c_flags, _v_n_bits * n_masks)); _nprot++;
    { SEXP _dim_flags = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_flags)[0] = _v_n_bits; INTEGER(_dim_flags)[1] = n_masks; Rf_setAttrib(flags, R_DimSymbol, _dim_flags); UNPROTECT(1); }

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
