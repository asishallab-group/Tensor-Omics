// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void mask_check_state_c(const int*, const int*, const int*, unsigned char*, int*);
void mask_chunk_count_c(const int*, int*, int*);
void calc_work_arr_paralog_subsets_size_c(int*, const int*, int*, const int*, const int*, int*);

SEXP mask_check_state_call(SEXP bit_mask, SEXP i_gene) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bit_mask_elements = (int) Rf_length(bit_mask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_i_gene = Rf_asInteger(i_gene);

    // convert what Fortran cannot take from R directly
    const int* _c_bit_mask = tox_int_in(bit_mask, "bit_mask");

    // outputs and work space
    unsigned char state = 0;
    int ierr = 0;

    mask_check_state_c(
        _c_bit_mask,
        &n_bit_mask_elements,
        &_v_i_gene,
        &state,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(state != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("state"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP mask_chunk_count_call(SEXP n_genes) {
    int _nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int _v_n_genes = Rf_asInteger(n_genes);

    // outputs and work space
    int count = 0;
    int ierr = 0;

    mask_chunk_count_c(
        &_v_n_genes,
        &count,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(count));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("count"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_work_arr_paralog_subsets_size_call(SEXP max_subset_size, SEXP n_genes, SEXP filtered_paralogs_mask) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_mask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_genes = Rf_asInteger(n_genes);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_mask = tox_int_in(filtered_paralogs_mask, "filtered_paralogs_mask");

    // outputs and work space
    int work_array_size = 0;
    int ierr = 0;

    calc_work_arr_paralog_subsets_size_c(
        &_v_max_subset_size,
        &_v_n_genes,
        &work_array_size,
        _c_filtered_paralogs_mask,
        &n_mask_chunks,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(_v_max_subset_size));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(work_array_size));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("max_subset_size"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("work_array_size"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
