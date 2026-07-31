// Generated. Do not edit.
#if !defined(NO_R_BINDING) && !defined(NO_C_BINDING)
#include <R.h>
#include <Rinternals.h>
#include "tox_marshal.h"

// the Fortran C-ABI symbols this module calls
void cluster_factor_trajectories_k_means_c(const int*, const double*, const int*, const int*, const int*, double*, int*, int*, int*, const int*);
void k_means_clustering_c(const int*, const double*, const int*, const int*, double*, int*, int*, int*, const int*);
void linkage_clustering_c(double*, const int*, int*, int*, double*, int*, const char*, int*);

SEXP cluster_factor_trajectories_k_means_call(SEXP trajectories, SEXP centroids, SEXP max_iterations) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_clusters = INTEGER(Rf_getAttrib(centroids, R_DimSymbol))[1];
    int n_factors = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[0];
    int n_samples = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[1];
    int n_timepoints = INTEGER(Rf_getAttrib(trajectories, R_DimSymbol))[2];

    // scalar inputs, pulled from their length-1 vectors
    int max_iterations_v = Rf_asInteger(max_iterations);

    // copy what is modified in place, so the caller's stays intact
    SEXP centroids_out = PROTECT(Rf_duplicate(centroids)); nprot++;

    // outputs and work space
    SEXP labels = PROTECT(Rf_allocVector(INTSXP, (n_samples*n_timepoints))); nprot++;
    SEXP label_counts = PROTECT(Rf_allocVector(INTSXP, n_clusters)); nprot++;
    int ierr = 0;

    cluster_factor_trajectories_k_means_c(
        &n_clusters,
        REAL(trajectories),
        &n_factors,
        &n_samples,
        &n_timepoints,
        REAL(centroids_out),
        INTEGER(labels),
        INTEGER(label_counts),
        &ierr,
        &max_iterations_v
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, centroids_out);
    SET_VECTOR_ELT(_out, 1, labels);
    SET_VECTOR_ELT(_out, 2, label_counts);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("labels"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("label_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP k_means_clustering_call(SEXP data_points, SEXP centroids, SEXP max_iterations) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_clusters = INTEGER(Rf_getAttrib(centroids, R_DimSymbol))[1];
    int n_points = INTEGER(Rf_getAttrib(data_points, R_DimSymbol))[1];
    int n_dims = INTEGER(Rf_getAttrib(data_points, R_DimSymbol))[0];

    // scalar inputs, pulled from their length-1 vectors
    int max_iterations_v = Rf_asInteger(max_iterations);

    // copy what is modified in place, so the caller's stays intact
    SEXP centroids_out = PROTECT(Rf_duplicate(centroids)); nprot++;

    // outputs and work space
    SEXP labels = PROTECT(Rf_allocVector(INTSXP, n_points)); nprot++;
    SEXP label_counts = PROTECT(Rf_allocVector(INTSXP, n_clusters)); nprot++;
    int ierr = 0;

    k_means_clustering_c(
        &n_clusters,
        REAL(data_points),
        &n_points,
        &n_dims,
        REAL(centroids_out),
        INTEGER(labels),
        INTEGER(label_counts),
        &ierr,
        &max_iterations_v
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 4)); nprot++;
    SET_VECTOR_ELT(_out, 0, centroids_out);
    SET_VECTOR_ELT(_out, 1, labels);
    SET_VECTOR_ELT(_out, 2, label_counts);
    SET_VECTOR_ELT(_out, 3, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 4)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("centroids"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("labels"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("label_counts"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

SEXP linkage_clustering_call(SEXP distances, SEXP method) {
    int nprot = 0;
    // derived from the inputs, not asked of the caller
    int n_points = INTEGER(Rf_getAttrib(distances, R_DimSymbol))[0];

    // copy what is modified in place, so the caller's stays intact
    SEXP distances_out = PROTECT(Rf_duplicate(distances)); nprot++;

    // convert what Fortran cannot take from R directly
    char* method_c = tox_char_in(method, 8);

    // outputs and work space
    SEXP merge_i = PROTECT(Rf_allocVector(INTSXP, (n_points - 1))); nprot++;
    SEXP merge_j = PROTECT(Rf_allocVector(INTSXP, (n_points - 1))); nprot++;
    SEXP heights = PROTECT(Rf_allocVector(REALSXP, (n_points - 1))); nprot++;
    SEXP cluster_sizes = PROTECT(Rf_allocVector(INTSXP, (n_points - 1))); nprot++;
    int ierr = 0;

    linkage_clustering_c(
        REAL(distances_out),
        &n_points,
        INTEGER(merge_i),
        INTEGER(merge_j),
        REAL(heights),
        INTEGER(cluster_sizes),
        method_c,
        &ierr
    );

    SEXP _out = PROTECT(Rf_allocVector(VECSXP, 6)); nprot++;
    SET_VECTOR_ELT(_out, 0, distances_out);
    SET_VECTOR_ELT(_out, 1, merge_i);
    SET_VECTOR_ELT(_out, 2, merge_j);
    SET_VECTOR_ELT(_out, 3, heights);
    SET_VECTOR_ELT(_out, 4, cluster_sizes);
    SET_VECTOR_ELT(_out, 5, Rf_ScalarInteger(ierr));
    SEXP _nms = PROTECT(Rf_allocVector(STRSXP, 6)); nprot++;
    SET_STRING_ELT(_nms, 0, Rf_mkChar("distances"));
    SET_STRING_ELT(_nms, 1, Rf_mkChar("merge_i"));
    SET_STRING_ELT(_nms, 2, Rf_mkChar("merge_j"));
    SET_STRING_ELT(_nms, 3, Rf_mkChar("heights"));
    SET_STRING_ELT(_nms, 4, Rf_mkChar("cluster_sizes"));
    SET_STRING_ELT(_nms, 5, Rf_mkChar("ierr"));
    Rf_setAttrib(_out, R_NamesSymbol, _nms);
    UNPROTECT(nprot);
    return _out;
}

#endif  // R binding
