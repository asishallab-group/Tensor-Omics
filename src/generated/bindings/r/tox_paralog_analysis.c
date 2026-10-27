// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void detect_neofunctionalization_c(const double*, const int*, const double*, const int*, const int*, const int*, const double*, unsigned char*, int*);
void detect_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, const double*, const double*, int*);
void detect_subfunctionalization_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, const double*, const double*, const int*, int*);
void filter_paralogs_by_pattern_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
void filter_paralogs_by_pattern_subfunctionalization_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);

SEXP detect_neofunctionalization_call(SEXP ancestors, SEXP genes, SEXP gene_to_fam, SEXP thresholds) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_families = INTEGER(Rf_getAttrib(ancestors, R_DimSymbol))[1];
    int n_axes = INTEGER(Rf_getAttrib(ancestors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    unsigned char* _c_neofunc = tox_bool_alloc(n_genes * n_axes);
    int ierr = 0;

    detect_neofunctionalization_c(
        REAL(ancestors),
        &n_families,
        REAL(genes),
        &n_axes,
        _c_gene_to_fam,
        &n_genes,
        REAL(thresholds),
        _c_neofunc,
        &ierr
    );

    // convert the outputs back
    SEXP neofunc = PROTECT(tox_bool_out(_c_neofunc, n_genes * n_axes)); _nprot++;
    { SEXP _dim_neofunc = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_neofunc)[0] = n_genes; INTEGER(_dim_neofunc)[1] = n_axes; Rf_setAttrib(neofunc, R_DimSymbol, _dim_neofunc); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, neofunc);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neofunc"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP detect_dosage_effect_call(SEXP ancestor, SEXP genes, SEXP filtered_paralogs_mask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP max_angle, SEXP gain_gamma) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_mask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_paralog_subsets = Rf_asInteger(n_paralog_subsets);
    double _v_max_angle = Rf_asReal(max_angle);
    double _v_gain_gamma = Rf_asReal(gain_gamma);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_mask = tox_int_in(filtered_paralogs_mask, "filtered_paralogs_mask");

    // outputs and work space
    int n_results = 0;
    SEXP work_arr_paralog_subsets = PROTECT(Rf_allocVector(INTSXP, n_mask_chunks * _v_n_paralog_subsets)); _nprot++;
    { SEXP _dim_work_arr_paralog_subsets = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_work_arr_paralog_subsets)[0] = n_mask_chunks; INTEGER(_dim_work_arr_paralog_subsets)[1] = _v_n_paralog_subsets; Rf_setAttrib(work_arr_paralog_subsets, R_DimSymbol, _dim_work_arr_paralog_subsets); UNPROTECT(1); }
    int ierr = 0;

    detect_dosage_effect_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        _c_filtered_paralogs_mask,
        &n_mask_chunks,
        &n_results,
        &_v_max_subset_size,
        INTEGER(work_arr_paralog_subsets),
        &_v_n_paralog_subsets,
        &_v_max_angle,
        &_v_gain_gamma,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, work_arr_paralog_subsets);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("work_arr_paralog_subsets"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP detect_subfunctionalization_call(SEXP ancestor, SEXP genes, SEXP filtered_paralogs_mask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP rdi_threshold, SEXP paralog_norms, SEXP sorted_paralog_norms_perm) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_mask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_paralog_subsets = Rf_asInteger(n_paralog_subsets);
    double _v_rdi_threshold = Rf_asReal(rdi_threshold);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_mask = tox_int_in(filtered_paralogs_mask, "filtered_paralogs_mask");
    const int* _c_sorted_paralog_norms_perm = tox_int_in(sorted_paralog_norms_perm, "sorted_paralog_norms_perm");

    // outputs and work space
    int n_results = 0;
    SEXP work_arr_paralog_subsets = PROTECT(Rf_allocVector(INTSXP, n_mask_chunks * _v_n_paralog_subsets)); _nprot++;
    { SEXP _dim_work_arr_paralog_subsets = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_work_arr_paralog_subsets)[0] = n_mask_chunks; INTEGER(_dim_work_arr_paralog_subsets)[1] = _v_n_paralog_subsets; Rf_setAttrib(work_arr_paralog_subsets, R_DimSymbol, _dim_work_arr_paralog_subsets); UNPROTECT(1); }
    int ierr = 0;

    detect_subfunctionalization_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        _c_filtered_paralogs_mask,
        &n_mask_chunks,
        &n_results,
        &_v_max_subset_size,
        INTEGER(work_arr_paralog_subsets),
        &_v_n_paralog_subsets,
        &_v_rdi_threshold,
        REAL(paralog_norms),
        _c_sorted_paralog_norms_perm,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, work_arr_paralog_subsets);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("work_arr_paralog_subsets"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP filter_paralogs_by_pattern_dosage_effect_call(SEXP gene_angles, SEXP threshold, SEXP n_families, SEXP gene_to_fam, SEXP n_mask_chunks) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(gene_angles);

    // scalar inputs, pulled from their length-1 vectors
    double _v_threshold = Rf_asReal(threshold);
    int _v_n_families = Rf_asInteger(n_families);
    int _v_n_mask_chunks = Rf_asInteger(n_mask_chunks);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    SEXP masks = PROTECT(Rf_allocVector(INTSXP, _v_n_mask_chunks * _v_n_families)); _nprot++;
    { SEXP _dim_masks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_masks)[0] = _v_n_mask_chunks; INTEGER(_dim_masks)[1] = _v_n_families; Rf_setAttrib(masks, R_DimSymbol, _dim_masks); UNPROTECT(1); }
    int ierr = 0;

    filter_paralogs_by_pattern_dosage_effect_c(
        REAL(gene_angles),
        &_v_threshold,
        &n_genes,
        &_v_n_families,
        _c_gene_to_fam,
        INTEGER(masks),
        &_v_n_mask_chunks,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP filter_paralogs_by_pattern_subfunctionalization_call(SEXP gene_angles, SEXP threshold, SEXP n_families, SEXP gene_to_fam, SEXP n_mask_chunks) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(gene_angles);

    // scalar inputs, pulled from their length-1 vectors
    double _v_threshold = Rf_asReal(threshold);
    int _v_n_families = Rf_asInteger(n_families);
    int _v_n_mask_chunks = Rf_asInteger(n_mask_chunks);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    SEXP masks = PROTECT(Rf_allocVector(INTSXP, _v_n_mask_chunks * _v_n_families)); _nprot++;
    { SEXP _dim_masks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_masks)[0] = _v_n_mask_chunks; INTEGER(_dim_masks)[1] = _v_n_families; Rf_setAttrib(masks, R_DimSymbol, _dim_masks); UNPROTECT(1); }
    int ierr = 0;

    filter_paralogs_by_pattern_subfunctionalization_c(
        REAL(gene_angles),
        &_v_threshold,
        &n_genes,
        &_v_n_families,
        _c_gene_to_fam,
        INTEGER(masks),
        &_v_n_mask_chunks,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
