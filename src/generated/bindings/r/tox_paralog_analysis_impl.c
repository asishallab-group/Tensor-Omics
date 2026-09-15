// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void calc_work_arr_paralog_subsets_size_c(int*, const int*, int*, const int*, const int*, int*);

SEXP calc_work_arr_paralog_subsets_size_call(SEXP max_subset_size, SEXP n_genes, SEXP filtered_paralogs_mask) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_mask);

    // scalar inputs, pulled from their length-1 vectors
    int max_subset_size_v = Rf_asInteger(max_subset_size);
    int n_genes_v = Rf_asInteger(n_genes);

    // outputs and work space
    int work_array_size = 0;
    int ierr = 0;

    calc_work_arr_paralog_subsets_size_c(
        &max_subset_size_v,
        &n_genes_v,
        &work_array_size,
        INTEGER(filtered_paralogs_mask),
        &n_mask_chunks,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(max_subset_size_v));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(work_array_size));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("max_subset_size"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("work_array_size"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
