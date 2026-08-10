// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void sample_estimator_anchors_c(const double*, const int*, const int*, int*, int*);
void grow_estimator_anchor_clouds_c(const double*, const int*, const int*, const int*, const int*, const double*, unsigned char*, int*, int*);
void estimate_stc_parameters_c(const double*, const int*, const int*, const int*, const int*, const int*, const double*, const int*, const double*, const double*, double*, double*, double*, double*, double*, double*, int*);

SEXP sample_estimator_anchors_call(SEXP density_labels, SEXP n_anchors) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_vectors = (int) Rf_length(density_labels);

    // scalar inputs, pulled from their length-1 vectors
    int n_anchors_v = Rf_asInteger(n_anchors);

    // outputs and work space
    SEXP anchor_indices = PROTECT(Rf_allocVector(INTSXP, n_anchors_v)); nprot++;
    int ierr = 0;

    sample_estimator_anchors_c(
        REAL(density_labels),
        &n_vectors,
        &n_anchors_v,
        INTEGER(anchor_indices),
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 2)); nprot++;
    SET_VECTOR_ELT(_out, 0, anchor_indices);
    SET_VECTOR_ELT(_out, 1, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 2)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("anchor_indices"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP grow_estimator_anchor_clouds_call(SEXP vectors, SEXP anchor_indices, SEXP seed_max_set_size) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];
    int n_anchors = (int) Rf_length(anchor_indices);

    // scalar inputs, pulled from their length-1 vectors
    double seed_max_set_size_v = Rf_asReal(seed_max_set_size);

    // outputs and work space
    unsigned char* cloud_masks_c = tox_bool_alloc(n_vectors * n_anchors);
    SEXP cloud_sizes = PROTECT(Rf_allocVector(INTSXP, n_anchors)); nprot++;
    int ierr = 0;

    grow_estimator_anchor_clouds_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(anchor_indices),
        &n_anchors,
        &seed_max_set_size_v,
        cloud_masks_c,
        INTEGER(cloud_sizes),
        &ierr
    );

    // convert the outputs back
    SEXP cloud_masks = PROTECT(tox_bool_out(cloud_masks_c, n_vectors * n_anchors)); nprot++;
    { SEXP cloud_masks_dim = PROTECT(Rf_allocVector(INTSXP, 2)); INTEGER(cloud_masks_dim)[0] = n_vectors; INTEGER(cloud_masks_dim)[1] = n_anchors; Rf_setAttrib(cloud_masks, R_DimSymbol, cloud_masks_dim); UNPROTECT(1); }

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 3)); nprot++;
    SET_VECTOR_ELT(_out, 0, cloud_masks);
    SET_VECTOR_ELT(_out, 1, cloud_sizes);
    SET_VECTOR_ELT(_out, 2, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 3)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("cloud_masks"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("cloud_sizes"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP estimate_stc_parameters_call(SEXP vectors, SEXP kd_indices, SEXP dimension_order, SEXP k_density, SEXP bandwidth_percentile, SEXP n_anchors, SEXP seed_max_set_size, SEXP first_quartile_percentile) {
    int nprot = 0;
    // optionals: a null pointer and size 0 when the caller omits them
    const int* k_density_p = NULL;
    int k_density_size = 0;
    if (k_density != R_NilValue) {
        k_density_size = (int) Rf_length(k_density);
        k_density_p = INTEGER(k_density);
    }
    const double* bandwidth_percentile_p = NULL;
    int bandwidth_percentile_size = 0;
    if (bandwidth_percentile != R_NilValue) {
        bandwidth_percentile_size = (int) Rf_length(bandwidth_percentile);
        bandwidth_percentile_p = REAL(bandwidth_percentile);
    }

    // derived from the inputs, not asked of the caller
    int n_dimensions = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[0];
    int n_vectors = INTEGER(Rf_getAttrib(vectors, R_DimSymbol))[1];

    // scalar inputs, pulled from their length-1 vectors
    int n_anchors_v = Rf_asInteger(n_anchors);
    double seed_max_set_size_v = Rf_asReal(seed_max_set_size);
    double first_quartile_percentile_v = Rf_asReal(first_quartile_percentile);

    // outputs and work space
    double estimated_k_min = 0;
    double estimated_k_density = 0;
    double estimated_density_quantile = 0;
    double estimated_chordal_dist_max_as_prcnt_of_range = 0;
    double estimated_G_max = 0;
    double estimated_d_max = 0;
    int ierr = 0;

    estimate_stc_parameters_c(
        REAL(vectors),
        &n_dimensions,
        &n_vectors,
        INTEGER(kd_indices),
        INTEGER(dimension_order),
        k_density_p,
        bandwidth_percentile_p,
        &n_anchors_v,
        &seed_max_set_size_v,
        &first_quartile_percentile_v,
        &estimated_k_min,
        &estimated_k_density,
        &estimated_density_quantile,
        &estimated_chordal_dist_max_as_prcnt_of_range,
        &estimated_G_max,
        &estimated_d_max,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 7)); nprot++;
    SET_VECTOR_ELT(_out, 0, Rf_ScalarReal(estimated_k_min));
    SET_VECTOR_ELT(_out, 1, Rf_ScalarReal(estimated_k_density));
    SET_VECTOR_ELT(_out, 2, Rf_ScalarReal(estimated_density_quantile));
    SET_VECTOR_ELT(_out, 3, Rf_ScalarReal(estimated_chordal_dist_max_as_prcnt_of_range));
    SET_VECTOR_ELT(_out, 4, Rf_ScalarReal(estimated_G_max));
    SET_VECTOR_ELT(_out, 5, Rf_ScalarReal(estimated_d_max));
    SET_VECTOR_ELT(_out, 6, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 7)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("estimated_k_min"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("estimated_k_density"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("estimated_density_quantile"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("estimated_chordal_dist_max_as_prcnt_of_range"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("estimated_G_max"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("estimated_d_max"));
    SET_STRING_ELT(_nms, 6, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
