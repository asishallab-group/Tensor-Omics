// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 217910824620d3a1 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void detect_neofunctionalization_c(const double*, const int*, const double*, const int*, const int*, const int*, const double*, unsigned char*, int*);
void detect_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, const double*, const double*, int*);
void detect_subfunctionalization_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, const double*, int*);
void detect_subfunctionalization_expert_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, int*, double*, const double*, const double*, const int*, double*, int*);
void filter_paralogs_by_pattern_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
void filter_paralogs_by_pattern_subfunctionalization_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
void calc_paralog_norms_c(const double*, const int*, const int*, double*, int*);
void calc_paralog_subsets_bitmasks_size_c(int*, const int*, int*, const int*, const int*, int*);

SEXP detect_neofunctionalization_call(SEXP ancestors, SEXP genes, SEXP gene_to_fam, SEXP thresholds) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_families = INTEGER(Rf_getAttrib(ancestors, R_DimSymbol))[1];
    int n_axes = INTEGER(Rf_getAttrib(ancestors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    R_xlen_t _n_neofunc = tox_count("neofunc", (double) n_genes * (double) n_axes);
    unsigned char* _c_neofunc = tox_bool_alloc(_n_neofunc);
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

    // convert the outputs back, unless the library reports an error
    SEXP neofunc = R_NilValue;
    if (!tox_failed(ierr)) {
        neofunc = PROTECT(tox_bool_out(_c_neofunc, _n_neofunc)); _nprot++;
        { SEXP _dim_neofunc = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_neofunc)[0] = n_genes; INTEGER(_dim_neofunc)[1] = n_axes; Rf_setAttrib(neofunc, R_DimSymbol, _dim_neofunc); UNPROTECT(1); }
    }

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

SEXP detect_dosage_effect_call(SEXP ancestor, SEXP genes, SEXP filtered_paralogs_bitmask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP max_angle, SEXP gain_gamma) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_bitmask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_paralog_subsets = Rf_asInteger(n_paralog_subsets);
    double _v_max_angle = Rf_asReal(max_angle);
    double _v_gain_gamma = Rf_asReal(gain_gamma);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_bitmask = tox_words_in(filtered_paralogs_bitmask, "filtered_paralogs_bitmask");

    // outputs and work space
    int n_results = 0;
    R_xlen_t _n_paralog_subsets_bitmasks = tox_count("paralog_subsets_bitmasks", (double) n_mask_chunks * (double) _v_n_paralog_subsets);
    SEXP paralog_subsets_bitmasks = PROTECT(Rf_allocVector(INTSXP, _n_paralog_subsets_bitmasks)); _nprot++;
    { SEXP _dim_paralog_subsets_bitmasks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_paralog_subsets_bitmasks)[0] = n_mask_chunks; INTEGER(_dim_paralog_subsets_bitmasks)[1] = _v_n_paralog_subsets; Rf_setAttrib(paralog_subsets_bitmasks, R_DimSymbol, _dim_paralog_subsets_bitmasks); UNPROTECT(1); }
    int ierr = 0;

    detect_dosage_effect_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        _c_filtered_paralogs_bitmask,
        &n_mask_chunks,
        &n_results,
        &_v_max_subset_size,
        INTEGER(paralog_subsets_bitmasks),
        &_v_n_paralog_subsets,
        &_v_max_angle,
        &_v_gain_gamma,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, paralog_subsets_bitmasks);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("paralog_subsets_bitmasks"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP detect_subfunctionalization_call(SEXP ancestor, SEXP genes, SEXP filtered_paralogs_bitmask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP rdi_threshold) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_bitmask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_paralog_subsets = Rf_asInteger(n_paralog_subsets);
    double _v_rdi_threshold = Rf_asReal(rdi_threshold);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_bitmask = tox_words_in(filtered_paralogs_bitmask, "filtered_paralogs_bitmask");

    // outputs and work space
    int n_results = 0;
    R_xlen_t _n_paralog_subsets_bitmasks = tox_count("paralog_subsets_bitmasks", (double) n_mask_chunks * (double) _v_n_paralog_subsets);
    SEXP paralog_subsets_bitmasks = PROTECT(Rf_allocVector(INTSXP, _n_paralog_subsets_bitmasks)); _nprot++;
    { SEXP _dim_paralog_subsets_bitmasks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_paralog_subsets_bitmasks)[0] = n_mask_chunks; INTEGER(_dim_paralog_subsets_bitmasks)[1] = _v_n_paralog_subsets; Rf_setAttrib(paralog_subsets_bitmasks, R_DimSymbol, _dim_paralog_subsets_bitmasks); UNPROTECT(1); }
    int ierr = 0;

    detect_subfunctionalization_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        _c_filtered_paralogs_bitmask,
        &n_mask_chunks,
        &n_results,
        &_v_max_subset_size,
        INTEGER(paralog_subsets_bitmasks),
        &_v_n_paralog_subsets,
        &_v_rdi_threshold,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, paralog_subsets_bitmasks);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("paralog_subsets_bitmasks"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP detect_subfunctionalization_expert_call(SEXP ancestor, SEXP genes, SEXP filtered_paralogs_bitmask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP rdi_threshold, SEXP paralog_norms, SEXP paralog_norms_perm) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_bitmask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_paralog_subsets = Rf_asInteger(n_paralog_subsets);
    double _v_rdi_threshold = Rf_asReal(rdi_threshold);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_bitmask = tox_words_in(filtered_paralogs_bitmask, "filtered_paralogs_bitmask");
    const int* _c_paralog_norms_perm = tox_int_in(paralog_norms_perm, "paralog_norms_perm");

    // outputs and work space
    int n_results = 0;
    R_xlen_t _n_paralog_subsets_bitmasks = tox_count("paralog_subsets_bitmasks", (double) n_mask_chunks * (double) _v_n_paralog_subsets);
    SEXP paralog_subsets_bitmasks = PROTECT(Rf_allocVector(INTSXP, _n_paralog_subsets_bitmasks)); _nprot++;
    { SEXP _dim_paralog_subsets_bitmasks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_paralog_subsets_bitmasks)[0] = n_mask_chunks; INTEGER(_dim_paralog_subsets_bitmasks)[1] = _v_n_paralog_subsets; Rf_setAttrib(paralog_subsets_bitmasks, R_DimSymbol, _dim_paralog_subsets_bitmasks); UNPROTECT(1); }
    R_xlen_t _n_tmp_active_bitmask = tox_count("tmp_active_bitmask", (double) n_mask_chunks);
    int* tmp_active_bitmask = (int*) R_alloc((size_t) _n_tmp_active_bitmask, sizeof(int));
    R_xlen_t _n_tmp_paralog_vector = tox_count("tmp_paralog_vector", (double) n_dims);
    double* tmp_paralog_vector = (double*) R_alloc((size_t) _n_tmp_paralog_vector, sizeof(double));
    R_xlen_t _n_tmp_work_array = tox_count("tmp_work_array", (double) n_genes);
    double* tmp_work_array = (double*) R_alloc((size_t) _n_tmp_work_array, sizeof(double));
    int ierr = 0;

    detect_subfunctionalization_expert_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        _c_filtered_paralogs_bitmask,
        &n_mask_chunks,
        &n_results,
        &_v_max_subset_size,
        INTEGER(paralog_subsets_bitmasks),
        &_v_n_paralog_subsets,
        tmp_active_bitmask,
        tmp_paralog_vector,
        &_v_rdi_threshold,
        REAL(paralog_norms),
        _c_paralog_norms_perm,
        tmp_work_array,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); _nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, paralog_subsets_bitmasks);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("paralog_subsets_bitmasks"));
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
    R_xlen_t _n_family_bitmasks = tox_count("family_bitmasks", (double) _v_n_mask_chunks * (double) _v_n_families);
    SEXP family_bitmasks = PROTECT(Rf_allocVector(INTSXP, _n_family_bitmasks)); _nprot++;
    { SEXP _dim_family_bitmasks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_family_bitmasks)[0] = _v_n_mask_chunks; INTEGER(_dim_family_bitmasks)[1] = _v_n_families; Rf_setAttrib(family_bitmasks, R_DimSymbol, _dim_family_bitmasks); UNPROTECT(1); }
    int ierr = 0;

    filter_paralogs_by_pattern_dosage_effect_c(
        REAL(gene_angles),
        &_v_threshold,
        &n_genes,
        &_v_n_families,
        _c_gene_to_fam,
        INTEGER(family_bitmasks),
        &_v_n_mask_chunks,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, family_bitmasks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_bitmasks"));
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
    R_xlen_t _n_family_bitmasks = tox_count("family_bitmasks", (double) _v_n_mask_chunks * (double) _v_n_families);
    SEXP family_bitmasks = PROTECT(Rf_allocVector(INTSXP, _n_family_bitmasks)); _nprot++;
    { SEXP _dim_family_bitmasks = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(_dim_family_bitmasks)[0] = _v_n_mask_chunks; INTEGER(_dim_family_bitmasks)[1] = _v_n_families; Rf_setAttrib(family_bitmasks, R_DimSymbol, _dim_family_bitmasks); UNPROTECT(1); }
    int ierr = 0;

    filter_paralogs_by_pattern_subfunctionalization_c(
        REAL(gene_angles),
        &_v_threshold,
        &n_genes,
        &_v_n_families,
        _c_gene_to_fam,
        INTEGER(family_bitmasks),
        &_v_n_mask_chunks,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, family_bitmasks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_bitmasks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_paralog_norms_call(SEXP genes) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dims = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];

    // outputs and work space
    R_xlen_t _n_paralog_norms = tox_count("paralog_norms", (double) n_genes);
    SEXP paralog_norms = PROTECT(Rf_allocVector(REALSXP, _n_paralog_norms)); _nprot++;
    int ierr = 0;

    calc_paralog_norms_c(
        REAL(genes),
        &n_dims,
        &n_genes,
        REAL(paralog_norms),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); _nprot++;
    SET_VECTOR_ELT(_out, 0, paralog_norms);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("paralog_norms"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP calc_paralog_subsets_bitmasks_size_call(SEXP max_subset_size, SEXP n_genes, SEXP filtered_paralogs_bitmask) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_bitmask);

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_subset_size = Rf_asInteger(max_subset_size);
    int _v_n_genes = Rf_asInteger(n_genes);

    // convert what Fortran cannot take from R directly
    const int* _c_filtered_paralogs_bitmask = tox_words_in(filtered_paralogs_bitmask, "filtered_paralogs_bitmask");

    // outputs and work space
    int work_array_size = 0;
    int ierr = 0;

    calc_paralog_subsets_bitmasks_size_c(
        &_v_max_subset_size,
        &_v_n_genes,
        &work_array_size,
        _c_filtered_paralogs_bitmask,
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
