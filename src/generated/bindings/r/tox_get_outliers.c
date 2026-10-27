// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void compute_family_scaling_c(const int*, const int*, const double*, const int*, double*, double*, double*, int*, const double*, const int*, const char*, const int*, double*, int*, int*);
void compute_rdi_c(const int*, const double*, const int*, const double*, const int*, double*, double*, int*, int*);
void identify_outliers_c(const int*, const double*, const double*, const int*, unsigned char*, double*, double*, const double*, int*);
void detect_outliers_c(const int*, const int*, const double*, const int*, unsigned char*, double*, double*, int*, double*, const double*, int*);

SEXP compute_family_scaling_call(SEXP n_families, SEXP distances, SEXP gene_to_fam, SEXP span, SEXP degree, SEXP mode, SEXP n_iters) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_families = Rf_asInteger(n_families);
    double _v_span = Rf_asReal(span);
    int _v_degree = Rf_asInteger(degree);
    int _v_n_iters = Rf_asInteger(n_iters);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");
    char* _c_mode = tox_char_in(mode, 6, "mode");

    // outputs and work space
    SEXP dscale = PROTECT(Rf_allocVector(REALSXP, _v_n_families)); _nprot++;
    SEXP loess_x = PROTECT(Rf_allocVector(REALSXP, _v_n_families)); _nprot++;
    SEXP loess_y = PROTECT(Rf_allocVector(REALSXP, _v_n_families)); _nprot++;
    SEXP indices_used = PROTECT(Rf_allocVector(INTSXP, _v_n_families)); _nprot++;
    double low_sd_cutoff = 0;
    SEXP excluded_low_sd = PROTECT(Rf_allocVector(INTSXP, _v_n_families)); _nprot++;
    int ierr = 0;

    compute_family_scaling_c(
        &n_genes,
        &_v_n_families,
        REAL(distances),
        _c_gene_to_fam,
        REAL(dscale),
        REAL(loess_x),
        REAL(loess_y),
        INTEGER(indices_used),
        &_v_span,
        &_v_degree,
        _c_mode,
        &_v_n_iters,
        &low_sd_cutoff,
        INTEGER(excluded_low_sd),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); _nprot++;
    SET_VECTOR_ELT(_out, 0, dscale);
    SET_VECTOR_ELT(_out, 1, loess_x);
    SET_VECTOR_ELT(_out, 2, loess_y);
    SET_VECTOR_ELT(_out, 3, indices_used);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(low_sd_cutoff));
    SET_VECTOR_ELT(_out, 5, excluded_low_sd);
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("dscale"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("loess_x"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("loess_y"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("indices_used"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("low_sd_cutoff"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("excluded_low_sd"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP compute_rdi_call(SEXP distances, SEXP gene_to_fam, SEXP dscale) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);
    int n_dscale_elements = (int) Rf_length(dscale);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    SEXP rdi = PROTECT(Rf_allocVector(REALSXP, n_genes)); _nprot++;
    SEXP sorted_rdi = PROTECT(Rf_allocVector(REALSXP, n_genes)); _nprot++;
    SEXP perm = PROTECT(Rf_allocVector(INTSXP, n_genes)); _nprot++;
    int ierr = 0;

    compute_rdi_c(
        &n_genes,
        REAL(distances),
        _c_gene_to_fam,
        REAL(dscale),
        &n_dscale_elements,
        REAL(rdi),
        REAL(sorted_rdi),
        INTEGER(perm),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, rdi);
    SET_VECTOR_ELT(_out, 1, sorted_rdi);
    SET_VECTOR_ELT(_out, 2, perm);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("rdi"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("sorted_rdi"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("perm"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP identify_outliers_call(SEXP rdi, SEXP sorted_rdi, SEXP perm, SEXP quantile_level) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(rdi);

    // scalar inputs, pulled from their length-1 vectors
    double _v_quantile_level = Rf_asReal(quantile_level);

    // convert what Fortran cannot take from R directly
    const int* _c_perm = tox_int_in(perm, "perm");

    // outputs and work space
    unsigned char* _c_is_outlier = tox_bool_alloc(n_genes);
    double threshold = 0;
    SEXP tail_probability = PROTECT(Rf_allocVector(REALSXP, n_genes)); _nprot++;
    int ierr = 0;

    identify_outliers_c(
        &n_genes,
        REAL(rdi),
        REAL(sorted_rdi),
        _c_perm,
        _c_is_outlier,
        &threshold,
        REAL(tail_probability),
        &_v_quantile_level,
        &ierr
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(_c_is_outlier, n_genes)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, is_outlier);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(threshold));
    SET_VECTOR_ELT(_out, 2, tail_probability);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("threshold"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("tail_probability"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP detect_outliers_call(SEXP n_families, SEXP distances, SEXP gene_to_fam, SEXP quantile_level) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(distances);

    // scalar inputs, pulled from their length-1 vectors
    int _v_n_families = Rf_asInteger(n_families);
    double _v_quantile_level = Rf_asReal(quantile_level);

    // convert what Fortran cannot take from R directly
    const int* _c_gene_to_fam = tox_int_in(gene_to_fam, "gene_to_fam");

    // outputs and work space
    unsigned char* _c_is_outlier = tox_bool_alloc(n_genes);
    SEXP loess_x = PROTECT(Rf_allocVector(REALSXP, _v_n_families)); _nprot++;
    SEXP loess_y = PROTECT(Rf_allocVector(REALSXP, _v_n_families)); _nprot++;
    SEXP loess_n = PROTECT(Rf_allocVector(INTSXP, _v_n_families)); _nprot++;
    SEXP tail_probability = PROTECT(Rf_allocVector(REALSXP, n_genes)); _nprot++;
    int ierr = 0;

    detect_outliers_c(
        &n_genes,
        &_v_n_families,
        REAL(distances),
        _c_gene_to_fam,
        _c_is_outlier,
        REAL(loess_x),
        REAL(loess_y),
        INTEGER(loess_n),
        REAL(tail_probability),
        &_v_quantile_level,
        &ierr
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(_c_is_outlier, n_genes)); _nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); _nprot++;
    SET_VECTOR_ELT(_out, 0, is_outlier);
    SET_VECTOR_ELT(_out, 1, loess_x);
    SET_VECTOR_ELT(_out, 2, loess_y);
    SET_VECTOR_ELT(_out, 3, loess_n);
    SET_VECTOR_ELT(_out, 4, tail_probability);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("loess_x"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("loess_y"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("loess_n"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("tail_probability"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
