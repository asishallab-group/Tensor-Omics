// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 0e1e7c507a726932 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void compute_family_direction_c(const int*, const int*, const int*, const double*, const int*, double*, double*, int*, int*, const double*, const double*, int*);
void compute_angular_deviations_c(const int*, const int*, const int*, const double*, const double*, const int*, double*, int*);
void compute_family_direction_rap_c(const int*, const int*, const double*, const int*, double*, double*, int*, int*, const double*, const double*, int*);
void compute_angular_deviations_rap_c(const int*, const int*, const double*, const double*, const int*, double*, int*);
void compute_relative_angular_deviations_c(const int*, const int*, const double*, const double*, const int*, double*, int*);
void compute_angle_outlier_threshold_c(const int*, const double*, double*, const double*, int*);
void flag_angle_outliers_c(const int*, const double*, const double*, unsigned char*, int*);
void detect_angle_outliers_c(const int*, const int*, const int*, const double*, const int*, double*, double*, int*, int*, double*, double*, unsigned char*, int*, const double*, const double*, const double*, int*);
void detect_angle_outliers_rap_c(const int*, const int*, const double*, const int*, double*, double*, int*, int*, double*, double*, unsigned char*, int*, const double*, const double*, const double*, int*);

SEXP compute_family_direction_call(SEXP n_families, SEXP expression_vectors, SEXP gene_to_fam, SEXP min_angular_dispersion, SEXP max_angular_dispersion) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    double min_angular_dispersion_v = Rf_asReal(min_angular_dispersion);
    double max_angular_dispersion_v = Rf_asReal(max_angular_dispersion);

    // outputs and work space
    SEXP family_directions = PROTECT(Rf_allocVector(REALSXP, n_axes * n_families_v)); nprot++;
    { SEXP family_directions_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(family_directions_dim)[0] = n_axes; INTEGER(family_directions_dim)[1] = n_families_v; Rf_setAttrib(family_directions, R_DimSymbol, family_directions_dim); UNPROTECT(1); }
    SEXP angular_dispersions = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP member_counts = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP status = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    int ierr = 0;

    compute_family_direction_c(
        &n_axes,
        &n_genes,
        &n_families_v,
        REAL(expression_vectors),
        INTEGER(gene_to_fam),
        REAL(family_directions),
        REAL(angular_dispersions),
        INTEGER(member_counts),
        INTEGER(status),
        &min_angular_dispersion_v,
        &max_angular_dispersion_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, family_directions);
    SET_VECTOR_ELT(_out, 1, angular_dispersions);
    SET_VECTOR_ELT(_out, 2, member_counts);
    SET_VECTOR_ELT(_out, 3, status);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_directions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("angular_dispersions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("member_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_angular_deviations_call(SEXP expression_vectors, SEXP family_directions, SEXP gene_to_fam) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];
    int n_families = INTEGER(Rf_getAttrib(family_directions, R_DimSymbol))[1];

    // outputs and work space
    SEXP angular_deviations = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_angular_deviations_c(
        &n_axes,
        &n_genes,
        &n_families,
        REAL(expression_vectors),
        REAL(family_directions),
        INTEGER(gene_to_fam),
        REAL(angular_deviations),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, angular_deviations);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("angular_deviations"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_family_direction_rap_call(SEXP n_families, SEXP signed_angles, SEXP gene_to_fam, SEXP min_angular_dispersion, SEXP max_angular_dispersion) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(signed_angles);

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    double min_angular_dispersion_v = Rf_asReal(min_angular_dispersion);
    double max_angular_dispersion_v = Rf_asReal(max_angular_dispersion);

    // outputs and work space
    SEXP family_mean_angles = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP angular_dispersions = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP member_counts = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP status = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    int ierr = 0;

    compute_family_direction_rap_c(
        &n_genes,
        &n_families_v,
        REAL(signed_angles),
        INTEGER(gene_to_fam),
        REAL(family_mean_angles),
        REAL(angular_dispersions),
        INTEGER(member_counts),
        INTEGER(status),
        &min_angular_dispersion_v,
        &max_angular_dispersion_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 5)); nprot++;
    SET_VECTOR_ELT(_out, 0, family_mean_angles);
    SET_VECTOR_ELT(_out, 1, angular_dispersions);
    SET_VECTOR_ELT(_out, 2, member_counts);
    SET_VECTOR_ELT(_out, 3, status);
    SET_VECTOR_ELT(_out, 4, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 5)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_mean_angles"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("angular_dispersions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("member_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_angular_deviations_rap_call(SEXP signed_angles, SEXP family_mean_angles, SEXP gene_to_fam) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(signed_angles);
    int n_families = (int) Rf_length(family_mean_angles);

    // outputs and work space
    SEXP angular_deviations = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_angular_deviations_rap_c(
        &n_genes,
        &n_families,
        REAL(signed_angles),
        REAL(family_mean_angles),
        INTEGER(gene_to_fam),
        REAL(angular_deviations),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, angular_deviations);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("angular_deviations"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_relative_angular_deviations_call(SEXP angular_deviations, SEXP angular_dispersions, SEXP gene_to_fam) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(angular_deviations);
    int n_families = (int) Rf_length(angular_dispersions);

    // outputs and work space
    SEXP relative_angular_deviations = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    int ierr = 0;

    compute_relative_angular_deviations_c(
        &n_genes,
        &n_families,
        REAL(angular_deviations),
        REAL(angular_dispersions),
        INTEGER(gene_to_fam),
        REAL(relative_angular_deviations),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, relative_angular_deviations);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("relative_angular_deviations"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP compute_angle_outlier_threshold_call(SEXP relative_angular_deviations, SEXP quantile_level) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(relative_angular_deviations);

    // scalar inputs, pulled from their length-1 vectors
    double quantile_level_v = Rf_asReal(quantile_level);

    // outputs and work space
    double threshold = 0;
    int ierr = 0;

    compute_angle_outlier_threshold_c(
        &n_genes,
        REAL(relative_angular_deviations),
        &threshold,
        &quantile_level_v,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(threshold));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("threshold"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP flag_angle_outliers_call(SEXP relative_angular_deviations, SEXP threshold) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(relative_angular_deviations);

    // scalar inputs, pulled from their length-1 vectors
    double threshold_v = Rf_asReal(threshold);

    // outputs and work space
    unsigned char* is_outlier_c = tox_bool_alloc(n_genes);
    int ierr = 0;

    flag_angle_outliers_c(
        &n_genes,
        REAL(relative_angular_deviations),
        &threshold_v,
        is_outlier_c,
        &ierr
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(is_outlier_c, n_genes)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, is_outlier);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_angle_outliers_call(SEXP n_families, SEXP expression_vectors, SEXP gene_to_fam, SEXP quantile_level, SEXP min_angular_dispersion, SEXP max_angular_dispersion) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_axes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[0];
    int n_genes = INTEGER(Rf_getAttrib(expression_vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    double quantile_level_v = Rf_asReal(quantile_level);
    double min_angular_dispersion_v = Rf_asReal(min_angular_dispersion);
    double max_angular_dispersion_v = Rf_asReal(max_angular_dispersion);

    // outputs and work space
    SEXP family_directions = PROTECT(Rf_allocVector(REALSXP, n_axes * n_families_v)); nprot++;
    { SEXP family_directions_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(family_directions_dim)[0] = n_axes; INTEGER(family_directions_dim)[1] = n_families_v; Rf_setAttrib(family_directions, R_DimSymbol, family_directions_dim); UNPROTECT(1); }
    SEXP angular_dispersions = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP member_counts = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP status = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP relative_angular_deviations = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    double threshold = 0;
    unsigned char* is_outlier_c = tox_bool_alloc(n_genes);
    SEXP gene_status = PROTECT(Rf_allocVector(INTSXP, n_genes)); nprot++;
    int ierr = 0;

    detect_angle_outliers_c(
        &n_axes,
        &n_genes,
        &n_families_v,
        REAL(expression_vectors),
        INTEGER(gene_to_fam),
        REAL(family_directions),
        REAL(angular_dispersions),
        INTEGER(member_counts),
        INTEGER(status),
        REAL(relative_angular_deviations),
        &threshold,
        is_outlier_c,
        INTEGER(gene_status),
        &quantile_level_v,
        &min_angular_dispersion_v,
        &max_angular_dispersion_v,
        &ierr
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(is_outlier_c, n_genes)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); nprot++;
    SET_VECTOR_ELT(_out, 0, family_directions);
    SET_VECTOR_ELT(_out, 1, angular_dispersions);
    SET_VECTOR_ELT(_out, 2, member_counts);
    SET_VECTOR_ELT(_out, 3, status);
    SET_VECTOR_ELT(_out, 4, relative_angular_deviations);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarReal(threshold));
    SET_VECTOR_ELT(_out, 6, is_outlier);
    SET_VECTOR_ELT(_out, 7, gene_status);
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_directions"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("angular_dispersions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("member_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("relative_angular_deviations"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("threshold"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("gene_status"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP detect_angle_outliers_rap_call(SEXP n_families, SEXP signed_angles, SEXP gene_to_fam, SEXP quantile_level, SEXP min_angular_dispersion, SEXP max_angular_dispersion) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_genes = (int) Rf_length(signed_angles);

    // scalar inputs, pulled from their length-1 vectors
    int n_families_v = Rf_asInteger(n_families);
    double quantile_level_v = Rf_asReal(quantile_level);
    double min_angular_dispersion_v = Rf_asReal(min_angular_dispersion);
    double max_angular_dispersion_v = Rf_asReal(max_angular_dispersion);

    // outputs and work space
    SEXP family_mean_angles = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP angular_dispersions = PROTECT(Rf_allocVector(REALSXP, n_families_v)); nprot++;
    SEXP member_counts = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP status = PROTECT(Rf_allocVector(INTSXP, n_families_v)); nprot++;
    SEXP relative_angular_deviations = PROTECT(Rf_allocVector(REALSXP, n_genes)); nprot++;
    double threshold = 0;
    unsigned char* is_outlier_c = tox_bool_alloc(n_genes);
    SEXP gene_status = PROTECT(Rf_allocVector(INTSXP, n_genes)); nprot++;
    int ierr = 0;

    detect_angle_outliers_rap_c(
        &n_genes,
        &n_families_v,
        REAL(signed_angles),
        INTEGER(gene_to_fam),
        REAL(family_mean_angles),
        REAL(angular_dispersions),
        INTEGER(member_counts),
        INTEGER(status),
        REAL(relative_angular_deviations),
        &threshold,
        is_outlier_c,
        INTEGER(gene_status),
        &quantile_level_v,
        &min_angular_dispersion_v,
        &max_angular_dispersion_v,
        &ierr
    );

    // convert the outputs back
    SEXP is_outlier = PROTECT(tox_bool_out(is_outlier_c, n_genes)); nprot++;

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 9)); nprot++;
    SET_VECTOR_ELT(_out, 0, family_mean_angles);
    SET_VECTOR_ELT(_out, 1, angular_dispersions);
    SET_VECTOR_ELT(_out, 2, member_counts);
    SET_VECTOR_ELT(_out, 3, status);
    SET_VECTOR_ELT(_out, 4, relative_angular_deviations);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarReal(threshold));
    SET_VECTOR_ELT(_out, 6, is_outlier);
    SET_VECTOR_ELT(_out, 7, gene_status);
    SET_VECTOR_ELT(_out, 8, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 9)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("family_mean_angles"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("angular_dispersions"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("member_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("status"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("relative_angular_deviations"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("threshold"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("is_outlier"));
    SET_STRING_ELT(_nms, 7, Rf_mkChar("gene_status"));
    SET_STRING_ELT(_nms, 8, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
