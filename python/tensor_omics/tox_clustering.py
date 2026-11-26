r"""tox_clustering

Clustering for tensor omics: k-means over factors and trajectories, and hierarchical
(agglomerative) linkage clustering over a precomputed distance matrix.

`k_means_clustering` partitions points in n dimensions; `cluster_factor_trajectories_k_means`
applies the same to whole time series, treating each sample's trajectory, over all its factors
and timepoints, as one point.
`linkage_clustering` takes the distances already computed and merges under the linkage
criterion asked for, so the same matrix can be re-clustered without recomputing it.

Python binding, generated from tox_clustering. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.cluster_factor_trajectories_k_means_c.restype = None
_lib.cluster_factor_trajectories_k_means_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENTS = ("n_clusters", "trajectories", "n_factors", "n_samples", "n_timepoints", "centroids", "labels", "label_counts", "max_iterations", "ierr",)
#: For a derived argument, the one the caller passed it in
_CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENT_SOURCES = ("centroids", None, "trajectories", "trajectories", "trajectories", None, None, None, None, None,)
#: For an array the C wrapper compares with its declaration, that shape
_CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENT_SHAPES = (None, None, None, None, None, "n_factors*n_timepoints, n_clusters", None, None, None, None,)

_lib.k_means_clustering_c.restype = None
_lib.k_means_clustering_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_K_MEANS_CLUSTERING_ARGUMENTS = ("n_clusters", "data_points", "n_points", "n_dims", "centroids", "labels", "label_counts", "max_iterations", "ierr",)
#: For a derived argument, the one the caller passed it in
_K_MEANS_CLUSTERING_ARGUMENT_SOURCES = ("centroids", None, "data_points", "data_points", None, None, None, None, None,)

_lib.linkage_clustering_c.restype = None
_lib.linkage_clustering_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_LINKAGE_CLUSTERING_ARGUMENTS = ("distances", "n_points", "merge_i", "merge_j", "heights", "cluster_sizes", "method", "ierr",)
#: For a derived argument, the one the caller passed it in
_LINKAGE_CLUSTERING_ARGUMENT_SOURCES = (None, "distances", None, None, None, None, None, None,)

def cluster_factor_trajectories_k_means(
        trajectories,
        centroids,
        max_iterations,
):
    r"""Performs k-means clustering on whole factor trajectories, so which samples evolve alike

    One sample is one point: its trajectory is flattened over factors and time into a single
    vector, `factor1_t1, factor2_t1, ..., factor1_t2, ...`, so the clustering answers which
    samples follow a similar course. Clustering the individual `(sample, timepoint)` states
    instead is :func:`tensor_omics.k_means_clustering` on the same
    array read as `(n_factors, n_samples*n_timepoints)`, which is a different question. The
    distances and means are scaled as that procedure documents, with the same limits.

    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints,), column-major (order='F')
        matrix with the trajectories to cluster
    centroids : np.ndarray[np.float64] of shape (n_factors*n_timepoints, n_clusters,), column-major (order='F'), modified in place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.

        One centroid is a whole flattened trajectory, so it has `n_factors*n_timepoints` entries.

        The final values will be the final centroids of the clusters. A cluster that ends up
        without samples (its `label_counts` is 0) keeps the centroid it had when it lost its
        last sample, or its starting centroid if it never had one.
    max_iterations : int
        number of maximum iterations of the clustering. Unlike
        :func:`tensor_omics.k_means_clustering`, which defaults to
        300, this one has to be given.
        The minimum valid value is `1`.

    Returns
    -------
    dict
        with keys:

        labels : np.ndarray[np.int32] of shape (n_samples,), read-only
            array of labels, one per sample, so the first label is the first sample's.

            each label is the index of its related cluster -> `1<=label<=n_clusters=k`
            A result is a value; call `.copy()` to obtain a modifiable array.
        label_counts : np.ndarray[np.int32] of shape (n_clusters,), read-only
            holds the number of samples having the respective label assigned
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_clustering::cluster_factor_trajectories_k_means`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    trajectories = _contract.reals(trajectories, "trajectories", ndim=3, order="F")
    if not isinstance(centroids, _np.ndarray) or centroids.dtype != _np.float64:
        raise _ToxArgumentError("'centroids' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="centroids")
    if centroids.ndim != 2:
        raise _ToxArgumentError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}", argument="centroids")
    if not centroids.flags.f_contiguous:
        raise _ToxArgumentError("'centroids' is modified in place, so it must already be column-major (order='F')", argument="centroids")
    max_iterations = _contract.integer(max_iterations, "max_iterations")

    # what the inputs already say, rather than asking for it again
    n_clusters = centroids.shape[1]
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]
    n_centroids_elements_dim_1 = centroids.shape[0]

    # outputs and work arrays, which the caller never sees
    labels = _np.empty((n_samples,), dtype=_np.int32, order='C')
    label_counts = _np.empty((n_clusters,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.cluster_factor_trajectories_k_means_c(
        _ctypes.byref(_ctypes.c_int(n_clusters)),
        trajectories,
        _ctypes.byref(_ctypes.c_int(n_factors)),
        _ctypes.byref(_ctypes.c_int(n_samples)),
        _ctypes.byref(_ctypes.c_int(n_timepoints)),
        centroids,
        _ctypes.byref(_ctypes.c_int(n_centroids_elements_dim_1)),
        labels,
        label_counts,
        _ctypes.byref(_ctypes.c_int(max_iterations)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENTS, _CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENT_SOURCES, shapes=_CLUSTER_FACTOR_TRAJECTORIES_K_MEANS_ARGUMENT_SHAPES)

    # a result is a value: modify a copy, not this
    labels.flags.writeable = False
    label_counts.flags.writeable = False

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

    1. Assigns each data point to one of `k` clusters whose centroid is closest
    2. Recalculates the centroids using the mean of its assigned points
    3. repeat 1-2 until assignment remains unchanged

    A point exactly as close to several centroids goes to the last of them, the one with the
    highest cluster index.

    Each point's squared distances are compared scaled by a power of two picked for that point,
    by f42's rule for sums of squares, from its largest coordinate difference to the nearest
    centroid in that measure, and each coordinate of a centroid is the mean of its cluster's
    values in that dimension, summed scaled by the power of two picked for that set. So for
    finite input the squared distances that decide a label neither overflow nor vanish beside
    larger values, and no centroid's sum overflows; only the distances to centroids far from a
    point may overflow, and those centroids lose. One limit stays: where the build or the host
    flushes subnormal numbers to zero, coordinates closer than `2**-1022` read as equal.

    Ordinary data needs no scaling and gets the arithmetic of the plain formulas: its results
    are bit-identical to those of the code before this scaling in builds that do not contract a
    multiplication and an addition into one fused operation (FMA); where a build does, a near
    tie can fall the other way.

    Parameters
    ----------
    data_points : np.ndarray[np.float64] of shape (n_dims, n_points,), column-major (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_dims, n_clusters,), column-major (order='F'), modified in place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.

        The final values will be the final centroids of the clusters. A cluster that ends up
        without points (its `label_counts` is 0) keeps the centroid it had when it lost its
        last point, or its starting centroid if it never had one.
    max_iterations : int, optional, default 300
        number of maximum iterations of the clustering. Defaults to 300, while
        :func:`tensor_omics.cluster_factor_trajectories_k_means`
        requires it.
        The minimum valid value is `1`.
        The default value is `300`.

    Returns
    -------
    dict
        with keys:

        labels : np.ndarray[np.int32] of shape (n_points,), read-only
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.

            each label is the index of its related cluster -> `1<=label<=n_clusters=k`
            A result is a value; call `.copy()` to obtain a modifiable array.
        label_counts : np.ndarray[np.int32] of shape (n_clusters,), read-only
            holds the number of points having the respective label assigned
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_clustering::k_means_clustering`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    data_points = _contract.reals(data_points, "data_points", ndim=2, order="F")
    if not isinstance(centroids, _np.ndarray) or centroids.dtype != _np.float64:
        raise _ToxArgumentError("'centroids' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="centroids")
    if centroids.ndim != 2:
        raise _ToxArgumentError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}", argument="centroids")
    if not centroids.flags.f_contiguous:
        raise _ToxArgumentError("'centroids' is modified in place, so it must already be column-major (order='F')", argument="centroids")
    max_iterations = _contract.integer(max_iterations, "max_iterations")

    # what the inputs already say, rather than asking for it again
    n_clusters = centroids.shape[1]
    n_points = data_points.shape[1]
    n_dims = data_points.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if centroids.shape[0] != n_dims:
        raise _ToxArgumentError(f"'centroids' has {centroids.shape[0]} along axis 0, but 'data_points' implies n_dims == {n_dims}; 'centroids' is declared with shape (n_dims, n_clusters)", argument="centroids")

    # outputs and work arrays, which the caller never sees
    labels = _np.empty((n_points,), dtype=_np.int32, order='C')
    label_counts = _np.empty((n_clusters,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.k_means_clustering_c(
        _ctypes.byref(_ctypes.c_int(n_clusters)),
        data_points,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_dims)),
        centroids,
        labels,
        label_counts,
        _ctypes.byref(_ctypes.c_int(max_iterations)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _K_MEANS_CLUSTERING_ARGUMENTS, _K_MEANS_CLUSTERING_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    labels.flags.writeable = False
    label_counts.flags.writeable = False

    return {
        "labels": labels,
        "label_counts": label_counts,
    }

def linkage_clustering(
        distances,
        method,
):
    r"""Perform linkage clustering on a distance matrix.

    The bottom triangle is used as scratch and restored from the top triangle before
    returning, on success or on error, so the matrix comes back unchanged. There is no
    need to copy it before calling.

    UPGMA and WPGMA average each pair of distances scaled by the power of two f42's rule picks
    for the larger, and Ward combines its three distances scaled for the largest of them, so for
    finite input no update overflows, and small distances beside large ones keep their value.
    Ward's distances can still exceed every input distance, as they grow with the clusters'
    sizes: a merge whose true height passes `huge` stops the run with `ERR_NAN_INF`. When an
    input exceeds `2**1007`, Ward keeps its scratch distances scaled by `2**-16`, so that one
    passing `huge` on the way is still compared right; inputs below `2**-1006` are then
    subnormal in that scratch, and read as 0 where subnormals are flushed.

    Every update is clamped where rounding can leave the range it lies in mathematically: an
    average into the pair it averages, a Ward update to at least the merge's height. So the
    heights never decrease, and equal distances average to themselves. Distances that need no
    scaling get the plain formulas: their results are bit-identical to those of the code before
    this scaling in builds that do not contract a multiplication and an addition into one fused
    operation (FMA), except where the clamps apply.

    Parameters
    ----------
    distances : np.ndarray[np.float64] of shape (n_points, n_points,), column-major (order='F'), modified in place
        symmetric distance matrix, holding the non-negative distances between points. Distance of X->X is always zero.

        The bottom triangle is used as scratch and restored from the top triangle before
        returning, on success or on error, so the matrix comes back unchanged. There is no
        need to copy it before calling.

        Its structure (symmetry, non-negativity, zero diagonal) is validated by the
        distance-matrix naming convention in the generated wrapper.
    method : str, one of 'average' | 'weighted' | 'ward'
        used algorithm
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

    Returns
    -------
    dict
        with keys:

        merge_i : np.ndarray[np.int32] of shape (n_points - 1,), read-only
            the first of the two clusters merged at merge k: a positive label is a point, its
            1-based index; a negative label `-j` is the cluster made by merge j. Of the two, this
            is the one whose lowest point index is the higher.

            The signs are the reverse of R's `hclust`, whose `merge` holds points as negative
            and merged clusters as positive numbers.
            A result is a value; call `.copy()` to obtain a modifiable array.
        merge_j : np.ndarray[np.int32] of shape (n_points - 1,), read-only
            the second of the two clusters merged at merge k, labelled as in `merge_i`: the one
            whose lowest point index is the lower.

            Two pairs at the same distance merge in the order of their lowest point indices: the
            pair whose lower cluster has the lowest point index first, then the one whose higher
            cluster has.
            A result is a value; call `.copy()` to obtain a modifiable array.
        heights : np.ndarray[np.float64] of shape (n_points - 1,), read-only
            the linkage distance at which merge k joins its two clusters, by the method's
            criterion. The heights never decrease from one merge to the next.
            A result is a value; call `.copy()` to obtain a modifiable array.
        cluster_sizes : np.ndarray[np.int32] of shape (n_points - 1,), read-only
            number of points in the cluster made by merge k
            A result is a value; call `.copy()` to obtain a modifiable array.

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_clustering::linkage_clustering`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    if not isinstance(distances, _np.ndarray) or distances.dtype != _np.float64:
        raise _ToxArgumentError("'distances' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="distances")
    if distances.ndim != 2:
        raise _ToxArgumentError(f"'distances' must have 2 dimensions, but has {distances.ndim}", argument="distances")
    if not distances.flags.f_contiguous:
        raise _ToxArgumentError("'distances' is modified in place, so it must already be column-major (order='F')", argument="distances")
    method = _np.array([_contract.mode(method, "method", ("average", "weighted", "ward",)).ljust(8)], dtype="S8")

    # what the inputs already say, rather than asking for it again
    n_points = distances.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if distances.shape[1] != n_points:
        raise _ToxArgumentError(f"'distances' has {distances.shape[1]} along axis 1, but 'distances' implies n_points == {n_points}; 'distances' is declared with shape (n_points, n_points)", argument="distances")

    # outputs and work arrays, which the caller never sees
    merge_i = _np.empty((max(0, n_points - 1),), dtype=_np.int32, order='C')
    merge_j = _np.empty((max(0, n_points - 1),), dtype=_np.int32, order='C')
    heights = _np.empty((max(0, n_points - 1),), dtype=_np.float64, order='C')
    cluster_sizes = _np.empty((max(0, n_points - 1),), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.linkage_clustering_c(
        distances,
        _ctypes.byref(_ctypes.c_int(n_points)),
        merge_i,
        merge_j,
        heights,
        cluster_sizes,
        method,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _LINKAGE_CLUSTERING_ARGUMENTS, _LINKAGE_CLUSTERING_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    merge_i.flags.writeable = False
    merge_j.flags.writeable = False
    heights.flags.writeable = False
    cluster_sizes.flags.writeable = False

    return {
        "merge_i": merge_i,
        "merge_j": merge_j,
        "heights": heights,
        "cluster_sizes": cluster_sizes,
    }
