// Generated. Do not edit.
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void mask_check_state_c(const int*, const int*, const int*, unsigned char*, int*);
void detect_neofunctionalization_c(const double*, const int*, const double*, const int*, const int*, const int*, const double*, unsigned char*, int*);
void detect_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, const int*, int*, const int*, int*, const int*, int*, double*, int*, const double*, const double*);
void detect_subfunctionalization_c(const double*, const double*, const int*, const int*, const double*, const int*, const int*, int*, const int*, int*, const int*, int*, double*, const double*, const int*, double*, int*);
void mask_chunk_count_c(const int*, int*, int*);
void filter_paralogs_by_pattern_subfunctionalization_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
void filter_paralogs_by_pattern_dosage_effect_c(const double*, const double*, const int*, const int*, const int*, int*, const int*, int*);
void calc_work_arr_paralog_subsets_size_c(int*, const int*, int*, const int*, const int*, int*);

SEXP mask_check_state_call(SEXP bit_mask, SEXP i_gene) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_bit_mask_elements = (int) Rf_length(bit_mask);

    // scalar inputs, pulled from their length-1 vectors
    int i_gene_v = Rf_asInteger(i_gene);

    // outputs and work space
    unsigned char state = 0;
    int ierr = 0;

    mask_check_state_c(
        INTEGER(bit_mask),
        &n_bit_mask_elements,
        &i_gene_v,
        &state,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarLogical(state != 0));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("state"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_neofunctionalization_call(SEXP ancestors, SEXP genes, SEXP gene_to_fam, SEXP thresholds) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_families = INTEGER(Rf_getAttrib(ancestors, R_DimSymbol))[1];
    int n_axes = INTEGER(Rf_getAttrib(ancestors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];

    // outputs and work space
    unsigned char* neofunc_c = tox_bool_alloc(n_genes * n_axes);
    int ierr = 0;

    detect_neofunctionalization_c(
        REAL(ancestors),
        &n_families,
        REAL(genes),
        &n_axes,
        INTEGER(gene_to_fam),
        &n_genes,
        REAL(thresholds),
        neofunc_c,
        &ierr
    );

    // convert the outputs back
    SEXP neofunc = PROTECT(tox_bool_out(neofunc_c, n_genes * n_axes)); nprot++;
    { SEXP neofunc_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(neofunc_dim)[0] = n_genes; INTEGER(neofunc_dim)[1] = n_axes; Rf_setAttrib(neofunc, R_DimSymbol, neofunc_dim); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, neofunc);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("neofunc"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_dosage_effect_call(SEXP ancestor, SEXP genes, SEXP filtered_paralogs_mask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP max_angle, SEXP gain_gamma) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_mask);

    // scalar inputs, pulled from their length-1 vectors
    int max_subset_size_v = Rf_asInteger(max_subset_size);
    int n_paralog_subsets_v = Rf_asInteger(n_paralog_subsets);
    double max_angle_v = Rf_asReal(max_angle);
    double gain_gamma_v = Rf_asReal(gain_gamma);

    // outputs and work space
    int n_results = 0;
    SEXP work_arr_paralog_subsets = PROTECT(Rf_allocVector(INTSXP, n_mask_chunks * n_paralog_subsets_v)); nprot++;
    { SEXP work_arr_paralog_subsets_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(work_arr_paralog_subsets_dim)[0] = n_mask_chunks; INTEGER(work_arr_paralog_subsets_dim)[1] = n_paralog_subsets_v; Rf_setAttrib(work_arr_paralog_subsets, R_DimSymbol, work_arr_paralog_subsets_dim); UNPROTECT(1); }
    int* tmp_active_mask = (int*) R_alloc(n_mask_chunks, sizeof(int));
    double* tmp_paralog_vector = (double*) R_alloc(n_dims, sizeof(double));
    int ierr = 0;

    detect_dosage_effect_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        INTEGER(filtered_paralogs_mask),
        &n_mask_chunks,
        &n_results,
        &max_subset_size_v,
        INTEGER(work_arr_paralog_subsets),
        &n_paralog_subsets_v,
        tmp_active_mask,
        tmp_paralog_vector,
        &ierr,
        &max_angle_v,
        &gain_gamma_v
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, work_arr_paralog_subsets);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("work_arr_paralog_subsets"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_subfunctionalization_call(SEXP ancestor, SEXP genes, SEXP rdi_threshold, SEXP filtered_paralogs_mask, SEXP max_subset_size, SEXP n_paralog_subsets, SEXP paralog_norms, SEXP sorted_paralog_norms_perm) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = INTEGER(Rf_getAttrib(genes, R_DimSymbol))[1];
    int n_dims = (int) Rf_length(ancestor);
    int n_mask_chunks = (int) Rf_length(filtered_paralogs_mask);

    // scalar inputs, pulled from their length-1 vectors
    double rdi_threshold_v = Rf_asReal(rdi_threshold);
    int max_subset_size_v = Rf_asInteger(max_subset_size);
    int n_paralog_subsets_v = Rf_asInteger(n_paralog_subsets);

    // outputs and work space
    int n_results = 0;
    SEXP work_arr_paralog_subsets = PROTECT(Rf_allocVector(INTSXP, n_mask_chunks * n_paralog_subsets_v)); nprot++;
    { SEXP work_arr_paralog_subsets_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(work_arr_paralog_subsets_dim)[0] = n_mask_chunks; INTEGER(work_arr_paralog_subsets_dim)[1] = n_paralog_subsets_v; Rf_setAttrib(work_arr_paralog_subsets, R_DimSymbol, work_arr_paralog_subsets_dim); UNPROTECT(1); }
    int* tmp_active_mask = (int*) R_alloc(n_mask_chunks, sizeof(int));
    double* tmp_paralog_vector = (double*) R_alloc(n_dims, sizeof(double));
    double* tmp_work_array = (double*) R_alloc(n_genes, sizeof(double));
    int ierr = 0;

    detect_subfunctionalization_c(
        REAL(ancestor),
        REAL(genes),
        &n_genes,
        &n_dims,
        &rdi_threshold_v,
        INTEGER(filtered_paralogs_mask),
        &n_mask_chunks,
        &n_results,
        &max_subset_size_v,
        INTEGER(work_arr_paralog_subsets),
        &n_paralog_subsets_v,
        tmp_active_mask,
        tmp_paralog_vector,
        REAL(paralog_norms),
        INTEGER(sorted_paralog_norms_perm),
        tmp_work_array,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(n_results));
    SET_VECTOR_ELT(_out, 1, work_arr_paralog_subsets);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("n_results"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("work_arr_paralog_subsets"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP mask_chunk_count_call(SEXP n_genes) {
    int nprot = 0;
    // scalar inputs, pulled from their length-1 vectors
    int n_genes_v = Rf_asInteger(n_genes);

    // outputs and work space
    int count = 0;
    int ierr = 0;

    mask_chunk_count_c(
        &n_genes_v,
        &count,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarInteger(count));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("count"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP filter_paralogs_by_pattern_subfunctionalization_call(SEXP gene_angles, SEXP threshold, SEXP n_families, SEXP gene_to_fam, SEXP n_mask_chunks) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(gene_angles);

    // scalar inputs, pulled from their length-1 vectors
    double threshold_v = Rf_asReal(threshold);
    int n_families_v = Rf_asInteger(n_families);
    int n_mask_chunks_v = Rf_asInteger(n_mask_chunks);

    // outputs and work space
    SEXP masks = PROTECT(Rf_allocVector(INTSXP, n_mask_chunks_v * n_families_v)); nprot++;
    { SEXP masks_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(masks_dim)[0] = n_mask_chunks_v; INTEGER(masks_dim)[1] = n_families_v; Rf_setAttrib(masks, R_DimSymbol, masks_dim); UNPROTECT(1); }
    int ierr = 0;

    filter_paralogs_by_pattern_subfunctionalization_c(
        REAL(gene_angles),
        &threshold_v,
        &n_genes,
        &n_families_v,
        INTEGER(gene_to_fam),
        INTEGER(masks),
        &n_mask_chunks_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP filter_paralogs_by_pattern_dosage_effect_call(SEXP gene_angles, SEXP threshold, SEXP n_families, SEXP gene_to_fam, SEXP n_mask_chunks) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(gene_angles);

    // scalar inputs, pulled from their length-1 vectors
    double threshold_v = Rf_asReal(threshold);
    int n_families_v = Rf_asInteger(n_families);
    int n_mask_chunks_v = Rf_asInteger(n_mask_chunks);

    // outputs and work space
    SEXP masks = PROTECT(Rf_allocVector(INTSXP, n_mask_chunks_v * n_families_v)); nprot++;
    { SEXP masks_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(masks_dim)[0] = n_mask_chunks_v; INTEGER(masks_dim)[1] = n_families_v; Rf_setAttrib(masks, R_DimSymbol, masks_dim); UNPROTECT(1); }
    int ierr = 0;

    filter_paralogs_by_pattern_dosage_effect_c(
        REAL(gene_angles),
        &threshold_v,
        &n_genes,
        &n_families_v,
        INTEGER(gene_to_fam),
        INTEGER(masks),
        &n_mask_chunks_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, masks);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

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
