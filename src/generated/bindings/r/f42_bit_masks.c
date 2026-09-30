// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void bit_mask_n_words_c(const int*, int*, int*);
void bit_mask_test_c(const int*, const int*, const int*, const int*, unsigned char*, int*);
void bit_mask_from_logical_c(const int*, const int*, const unsigned char*, int*, int*);
void bit_mask_to_logical_c(const int*, const int*, const int*, unsigned char*, int*);
void bit_masks_from_logical_2D_c(const int*, const int*, const int*, const unsigned char*, int*, int*);
void bit_masks_to_logical_2D_c(const int*, const int*, const int*, const int*, unsigned char*, int*);

SEXP bit_mask_n_words_call(SEXP n_bits) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_bits_v = Rf_asInteger(n_bits);

    // outputs and work space
    int n_words = 0;
    int ierr = 0;

    bit_mask_n_words_c(
        &n_bits_v,
        &n_words,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_words));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_words"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bit_mask_test_call(SEXP n_bits, SEXP bit_mask, SEXP i_bit) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_words = (int) Rf_length(bit_mask);

    // scalar inputs, pulled from their length-1 vectors
    int n_bits_v = Rf_asInteger(n_bits);
    int i_bit_v = Rf_asInteger(i_bit);

    // outputs and work space
    unsigned char is_set = 0;
    int ierr = 0;

    bit_mask_test_c(
        &n_bits_v,
        &n_words,
        INTEGER(bit_mask),
        &i_bit_v,
        &is_set,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(is_set != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_set"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bit_mask_from_logical_call(SEXP n_words, SEXP flags) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bits = (int) Rf_length(flags);

    // scalar inputs, pulled from their length-1 vectors
    int n_words_v = Rf_asInteger(n_words);

    // convert what Fortran cannot take from R directly
    unsigned char* flags_c = tox_bool_in(flags);

    // outputs and work space
    SEXP bit_mask = PROTECT(Rf_allocVector(INTSXP, n_words_v)); nprot++;
    int ierr = 0;

    bit_mask_from_logical_c(
        &n_bits,
        &n_words_v,
        flags_c,
        INTEGER(bit_mask),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, bit_mask);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("bit_mask"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bit_mask_to_logical_call(SEXP n_bits, SEXP bit_mask) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_words = (int) Rf_length(bit_mask);

    // scalar inputs, pulled from their length-1 vectors
    int n_bits_v = Rf_asInteger(n_bits);

    // outputs and work space
    unsigned char* flags_c = tox_bool_alloc(n_bits_v);
    int ierr = 0;

    bit_mask_to_logical_c(
        &n_bits_v,
        &n_words,
        INTEGER(bit_mask),
        flags_c,
        &ierr
    );

    // convert the outputs back
    SEXP flags = PROTECT(tox_bool_out(flags_c, n_bits_v)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, flags);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("flags"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bit_masks_from_logical_2D_call(SEXP n_words, SEXP flags) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bits = INTEGER(Rf_getAttrib(flags, R_DimSymbol))[0];
    int n_masks = INTEGER(Rf_getAttrib(flags, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_words_v = Rf_asInteger(n_words);

    // convert what Fortran cannot take from R directly
    unsigned char* flags_c = tox_bool_in(flags);

    // outputs and work space
    SEXP bit_masks = PROTECT(Rf_allocVector(INTSXP, n_words_v * n_masks)); nprot++;
    { SEXP bit_masks_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(bit_masks_dim)[0] = n_words_v; INTEGER(bit_masks_dim)[1] = n_masks; Rf_setAttrib(bit_masks, R_DimSymbol, bit_masks_dim); UNPROTECT(1); }
    int ierr = 0;

    bit_masks_from_logical_2D_c(
        &n_bits,
        &n_masks,
        &n_words_v,
        flags_c,
        INTEGER(bit_masks),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, bit_masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("bit_masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP bit_masks_to_logical_2D_call(SEXP n_bits, SEXP bit_masks) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_masks = INTEGER(Rf_getAttrib(bit_masks, R_DimSymbol))[1];
    int n_words = INTEGER(Rf_getAttrib(bit_masks, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int n_bits_v = Rf_asInteger(n_bits);

    // outputs and work space
    unsigned char* flags_c = tox_bool_alloc(n_bits_v * n_masks);
    int ierr = 0;

    bit_masks_to_logical_2D_c(
        &n_bits_v,
        &n_masks,
        &n_words,
        INTEGER(bit_masks),
        flags_c,
        &ierr
    );

    // convert the outputs back
    SEXP flags = PROTECT(tox_bool_out(flags_c, n_bits_v * n_masks)); nprot++;
    { SEXP flags_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(flags_dim)[0] = n_bits_v; INTEGER(flags_dim)[1] = n_masks; Rf_setAttrib(flags, R_DimSymbol, flags_dim); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, flags);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("flags"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
