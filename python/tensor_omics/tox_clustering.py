"""Python binding to Generated from the kernel; do not edit -- regenerate instead.

Generated from tox_clustering. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.cluster_factor_trajectories_k_means_c.restype = None
_lib.cluster_factor_trajectories_k_means_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENTS = ("n_clusters", "trajectories", "n_factors", "n_samples", "n_timepoints", "centroids", "labels", "label_counts", "ierr", "max_iterations",)

_lib.k_means_clustering_c.restype = None
_lib.k_means_clustering_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_K_MEANS_CLUSTERING_ARGUMENTS = ("n_clusters", "data_points", "n_points", "n_dims", "centroids", "labels", "label_counts", "ierr", "max_iterations",)

_lib.linkage_clustering_c.restype = None
_lib.linkage_clustering_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LINKAGE_CLUSTERING_ARGUMENTS = ("distances", "n_points", "merge_i", "merge_j", "heights", "cluster_sizes", "method", "ierr",)

def cluster_factor_trajectories_k_means(
        trajectories,
        centroids,
        max_iterations,
):
    r"""Performs k-means clustering on factor trajectories, so factor evolution over time

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_factors, n_clusters,), column-major (order='F'), modified in place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.

        The final values will be the final centroids of the clusters
    max_iterations : int
        number of maximum iterations of the clustering

    Returns
    -------
    dict
        with keys:

        labels : np.ndarray[np.int32] of shape (n_samples*n_timepoints,)
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.

            each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        label_counts : np.ndarray[np.int32] of shape (n_clusters,)
            holds the number of points having the respective label assigned

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_clustering::cluster_factor_trajectories_k_means`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'trajectories' must be an array of np.float64: {error}") from None
    if trajectories.ndim != 3:
        raise ValueError(f"'trajectories' must have 3 dimensions, but has {trajectories.ndim}")
    if not isinstance(centroids, np.ndarray) or centroids.dtype != np.float64:
        raise TypeError("'centroids' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if centroids.ndim != 2:
        raise ValueError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}")
    if not centroids.flags.f_contiguous:
        raise ValueError("'centroids' is modified in place, so it must already be column-major (order='F')")

    # what the inputs already say, rather than asking for it again
    n_clusters = centroids.shape[1]
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if centroids.shape[0] != n_factors:
        raise ValueError(f"'centroids' has {centroids.shape[0]} along axis 0, but "
            f"'trajectories' implies n_factors == {n_factors}"
        )

    # outputs and work arrays, which the caller never sees
    labels = np.empty((n_samples*n_timepoints,), dtype=np.int32, order='C')
    label_counts = np.empty((n_clusters,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.cluster_factor_trajectories_k_means_c(
        ctypes.byref(ctypes.c_int(n_clusters)),
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        centroids,
        labels,
        label_counts,
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_int(max_iterations)),
    )

    check_err_code(ierr.value, _CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENTS)

    return {
        "labels": labels,
        "label_counts": label_counts,
    }

def k_means_clustering(
        data_points,
        centroids,
        max_iterations=300,
):
    r"""k-means clustering algorithm

    Parameters
    ----------
    data_points : np.ndarray[np.float64] of shape (n_dims, n_points,), column-major (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_dims, n_clusters,), column-major (order='F'), modified in place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.

        The final values will be the final centroids of the clusters
    max_iterations : int, optional, default 300
        number of maximum iterations of the clustering.
        The default value is `300_int32`.

    Returns
    -------
    dict
        with keys:

        labels : np.ndarray[np.int32] of shape (n_points,)
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.

            each label is the index of its related cluster -> `1<=label<=n_clusters=k`
        label_counts : np.ndarray[np.int32] of shape (n_clusters,)
            holds the number of points having the respective label assigned

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_clustering::k_means_clustering`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        data_points = np.asfortranarray(data_points, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'data_points' must be an array of np.float64: {error}") from None
    if data_points.ndim != 2:
        raise ValueError(f"'data_points' must have 2 dimensions, but has {data_points.ndim}")
    if not isinstance(centroids, np.ndarray) or centroids.dtype != np.float64:
        raise TypeError("'centroids' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if centroids.ndim != 2:
        raise ValueError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}")
    if not centroids.flags.f_contiguous:
        raise ValueError("'centroids' is modified in place, so it must already be column-major (order='F')")

    # what the inputs already say, rather than asking for it again
    n_clusters = centroids.shape[1]
    n_points = data_points.shape[1]
    n_dims = data_points.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if centroids.shape[0] != n_dims:
        raise ValueError(f"'centroids' has {centroids.shape[0]} along axis 0, but "
            f"'data_points' implies n_dims == {n_dims}"
        )

    # outputs and work arrays, which the caller never sees
    labels = np.empty((n_points,), dtype=np.int32, order='C')
    label_counts = np.empty((n_clusters,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.k_means_clustering_c(
        ctypes.byref(ctypes.c_int(n_clusters)),
        data_points,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_dims)),
        centroids,
        labels,
        label_counts,
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_int(max_iterations)),
    )

    check_err_code(ierr.value, _K_MEANS_CLUSTERING_ARGUMENTS)

    return {
        "labels": labels,
        "label_counts": label_counts,
    }

def linkage_clustering(
        distances,
        method,
):
    r"""Perform linkage clustering on a distance matrix.

    Parameters
    ----------
    distances : np.ndarray[np.float64] of shape (n_points, n_points,), column-major (order='F'), modified in place
        symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.

        @note
        This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
        So there is no need to copy an existing distance matrix, just pass the original.
        @endnote

        The distance-matrix structure (symmetry, non-negativity, zero diagonal) is checked below,
        not by the finiteness contract, so this argument opts out of that.
        NaN is permitted for this value.
        Infinite values are permitted for this value.
    method : str, one of 'average' | 'weighted' | 'ward'
        used algorithm
        The minimum valid value is `0_int32`.
        The maximum valid value is `2_int32`.

    Returns
    -------
    dict
        with keys:

        merge_i : np.ndarray[np.int32] of shape (n_points - 1,)
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        merge_j : np.ndarray[np.int32] of shape (n_points - 1,)
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes
        heights : np.ndarray[np.float64] of shape (n_points - 1,)
            height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter
        cluster_sizes : np.ndarray[np.int32] of shape (n_points - 1,)
            size of cluster at iteration k

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_clustering::linkage_clustering`.
    """
    # accept anything array-like, converting only when C needs it
    if not isinstance(distances, np.ndarray) or distances.dtype != np.float64:
        raise TypeError("'distances' is modified in place, so it must already be a numpy array of {}".format(np.float64))
    if distances.ndim != 2:
        raise ValueError(f"'distances' must have 2 dimensions, but has {distances.ndim}")
    if not distances.flags.f_contiguous:
        raise ValueError("'distances' is modified in place, so it must already be column-major (order='F')")
    method = np.array([str(method).lower().encode()], dtype="S8")

    # what the inputs already say, rather than asking for it again
    n_points = distances.shape[0]

    # outputs and work arrays, which the caller never sees
    merge_i = np.empty((n_points - 1,), dtype=np.int32, order='C')
    merge_j = np.empty((n_points - 1,), dtype=np.int32, order='C')
    heights = np.empty((n_points - 1,), dtype=np.float64, order='C')
    cluster_sizes = np.empty((n_points - 1,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.linkage_clustering_c(
        distances,
        ctypes.byref(ctypes.c_int(n_points)),
        merge_i,
        merge_j,
        heights,
        cluster_sizes,
        method,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _LINKAGE_CLUSTERING_ARGUMENTS)

    return {
        "merge_i": merge_i,
        "merge_j": merge_j,
        "heights": heights,
        "cluster_sizes": cluster_sizes,
    }
