// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"
// tox_marshal.h 7902a51790babb50 -- its hash, so that fpm, which only hashes this file, recompiles it when the header changes

// the Fortran C-ABI symbols this module calls
void cluster_factor_trajectories_k_means_c(const int*, const double*, const int*, const int*, const int*, double*, int*, int*, const int*, int*);
void k_means_clustering_c(const int*, const double*, const int*, const int*, double*, int*, int*, const int*, int*);
void linkage_clustering_c(double*, const int*, int*, int*, double*, int*, const char*, int*);

SEXP cluster_factor_trajectories_k_means_call(SEXP trajectories, SEXP centroids, SEXP max_iterations) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_clusters = INTEGER(Rf_getAttrib(centroids, R_DimSymbol))[1];
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_iterations = Rf_asInteger(max_iterations);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_centroids = PROTECT(Rf_duplicate(centroids)); _nprot++;

    // outputs and work space
    SEXP labels = PROTECT(Rf_allocVector(INTSXP, (n_samples * n_timepoints))); _nprot++;
    SEXP label_counts = PROTECT(Rf_allocVector(INTSXP, n_clusters)); _nprot++;
    int ierr = 0;

    cluster_factor_trajectories_k_means_c(
        &n_clusters,
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        REAL(_out_centroids),
        INTEGER(labels),
        INTEGER(label_counts),
        &_v_max_iterations,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_centroids);
    SET_VECTOR_ELT(_out, 1, labels);
    SET_VECTOR_ELT(_out, 2, label_counts);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("labels"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("label_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP k_means_clustering_call(SEXP data_points, SEXP centroids, SEXP max_iterations) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_clusters = INTEGER(Rf_getAttrib(centroids, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(data_points, R_DimSymbol))[1];
    int n_dims = INTEGER(Rf_getAttrib(data_points, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int _v_max_iterations = Rf_asInteger(max_iterations);

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_centroids = PROTECT(Rf_duplicate(centroids)); _nprot++;

    // outputs and work space
    SEXP labels = PROTECT(Rf_allocVector(INTSXP, n_points)); _nprot++;
    SEXP label_counts = PROTECT(Rf_allocVector(INTSXP, n_clusters)); _nprot++;
    int ierr = 0;

    k_means_clustering_c(
        &n_clusters,
        REAL(data_points),
        &n_points,
        &n_dims,
        REAL(_out_centroids),
        INTEGER(labels),
        INTEGER(label_counts),
        &_v_max_iterations,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_centroids);
    SET_VECTOR_ELT(_out, 1, labels);
    SET_VECTOR_ELT(_out, 2, label_counts);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("labels"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("label_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

SEXP linkage_clustering_call(SEXP distances, SEXP method) {
    int _nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(distances, R_DimSymbol))[0];

    // copy what is modified in place, so the caller's stays intact
    SEXP _out_distances = PROTECT(Rf_duplicate(distances)); _nprot++;

    // convert what Fortran cannot take from R directly
    char* _c_method = tox_char_in(method, 8, "method");

    // outputs and work space
    SEXP merge_i = PROTECT(Rf_allocVector(INTSXP, (n_points - 1))); _nprot++;
    SEXP merge_j = PROTECT(Rf_allocVector(INTSXP, (n_points - 1))); _nprot++;
    SEXP heights = PROTECT(Rf_allocVector(REALSXP, (n_points - 1))); _nprot++;
    SEXP cluster_sizes = PROTECT(Rf_allocVector(INTSXP, (n_points - 1))); _nprot++;
    int ierr = 0;

    linkage_clustering_c(
        REAL(_out_distances),
        &n_points,
        INTEGER(merge_i),
        INTEGER(merge_j),
        REAL(heights),
        INTEGER(cluster_sizes),
        _c_method,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); _nprot++;
    SET_VECTOR_ELT(_out, 0, _out_distances);
    SET_VECTOR_ELT(_out, 1, merge_i);
    SET_VECTOR_ELT(_out, 2, merge_j);
    SET_VECTOR_ELT(_out, 3, heights);
    SET_VECTOR_ELT(_out, 4, cluster_sizes);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); _nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("distances"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("merge_i"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("merge_j"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("heights"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("cluster_sizes"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(_nprot);
    return _out;
}

#endif  // R binding
