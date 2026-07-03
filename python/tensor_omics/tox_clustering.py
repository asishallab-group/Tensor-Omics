from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def cluster_factor_trajectories_k_means(
        trajectories,
        centroids,
        max_iterations=300
        ):
    """
    Parameters
    ----------
    trajectories : np.ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_factors, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.The centroids should be unique. This is not checked in this routine.The final values will be the final centroids of the clusters
    max_iterations : int, optional
        number of maximum iterations of the clusteringThe default value is `300_int32`.

    Returns
    -------
    results : dict
        labels : np.ndarray[np.int32] of shape (n_samples * n_timepoints,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : np.ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned


    Notes
    -----
    Performs k-means clustering on factor trajectories, so factor evolution over time
    """

    # ensure all array inputs are numpy arrays
    trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    assert type(centroids) is np.ndarray and centroids.flags.f_contiguous and centroids.dtype == np.float64, "'centroids' must be column-major numpy array (order='F')"

    # extract dimension arguments
    n_clusters = centroids.shape[1]
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]


    # Create temporaries and/or outputs
    labels = np.empty((n_samples * n_timepoints,), dtype=np.int32, order='F')
    label_counts = np.empty((n_clusters,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)
    max_iterations = ctypes.c_int(max_iterations)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.cluster_factor_trajectories_k_means_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=3, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.cluster_factor_trajectories_k_means_c.restype = None

    tox.cluster_factor_trajectories_k_means_c(
        ctypes.byref(ctypes.c_int(n_clusters)),
        trajectories,
        ctypes.byref(ctypes.c_int(n_factors)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ctypes.c_int(n_timepoints)),
        centroids,
        labels,
        label_counts,
        ctypes.byref(ierr),
        ctypes.byref(max_iterations)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroids.setflags(write=False)
    labels.setflags(write=False)
    label_counts.setflags(write=False)

    return {
        "labels": labels,
        "label_counts": label_counts
    }


def k_means_clustering(
        data_points,
        centroids,
        max_iterations=300
        ):
    """
    performs k-means clustering on `data_points`


    Parameters
    ----------
    data_points : np.ndarray[np.float64] of shape (n_dims, n_points) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : np.ndarray[np.float64] of shape (n_dims, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.The centroids should be unique. This is not checked in this routine.The final values will be the final centroids of the clusters
    max_iterations : int, optional
        number of maximum iterations of the clusteringThe default value is `300_int32`.

    Returns
    -------
    results : dict
        labels : np.ndarray[np.int32] of shape (n_points,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : np.ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned


    Notes
    -----
    k-means clustering algorithm:1. Assigns each data point to one of `k` clusters whose centroid is clostest2. Recalculates the centroids using the mean of its assigned points3. repeat 1-2 until assignment remains unchanged
    """

    # ensure all array inputs are numpy arrays
    data_points = np.asfortranarray(data_points, dtype=np.float64)
    assert type(centroids) is np.ndarray and centroids.flags.f_contiguous and centroids.dtype == np.float64, "'centroids' must be column-major numpy array (order='F')"

    # extract dimension arguments
    n_clusters = centroids.shape[1]
    n_points = data_points.shape[1]
    n_dims = data_points.shape[0]


    # Create temporaries and/or outputs
    labels = np.empty((n_points,), dtype=np.int32, order='F')
    label_counts = np.empty((n_clusters,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)
    max_iterations = ctypes.c_int(max_iterations)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.k_means_clustering_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.k_means_clustering_c.restype = None

    tox.k_means_clustering_c(
        ctypes.byref(ctypes.c_int(n_clusters)),
        data_points,
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_dims)),
        centroids,
        labels,
        label_counts,
        ctypes.byref(ierr),
        ctypes.byref(max_iterations)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroids.setflags(write=False)
    labels.setflags(write=False)
    label_counts.setflags(write=False)

    return {
        "labels": labels,
        "label_counts": label_counts
    }


def linkage_clustering(
        distances,
        method
        ):
    """
    Perform linkage clustering on a distance matrix.


    Parameters
    ----------
    distances : np.ndarray[np.float64] of shape (n_points, n_points) in column-major layout (order='F'), modified in-place
        symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.@noteThis subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.So there is no need to copy an existing distance matrix, just pass the original.@endnote
    method : str
        used algorithm|      Method      |   Value    |
        |------------------|------------|
        | Average / UPGMA  | "average"  |
        | Weighted / WPGMA | "weighted" |
        |       Ward       |   "ward"   |

    Returns
    -------
    results : dict
        merge_i : np.ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes,
        merge_j : np.ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes,
        heights : np.ndarray[np.float64] of shape (n_points - 1,) in column-major layout (order='F')
            height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter,
        cluster_sizes : np.ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            size of cluster at iteration k


    Notes
    -----
    @noteThis subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.So there is no need to copy an existing distance matrix, just pass the original.@endnote
    """

    # ensure all array inputs are numpy arrays
    assert type(distances) is np.ndarray and distances.flags.f_contiguous and distances.dtype == np.float64, "'distances' must be column-major numpy array (order='F')"
    method = np.asarray(method)

    # extract dimension arguments
    n_points = distances.shape[0]


    # Create temporaries and/or outputs
    merge_i = np.empty((n_points - 1,), dtype=np.int32, order='F')
    merge_j = np.empty((n_points - 1,), dtype=np.int32, order='F')
    heights = np.empty((n_points - 1,), dtype=np.float64, order='F')
    cluster_sizes = np.empty((n_points - 1,), dtype=np.int32, order='F')
    method = method.astype(f"S{8}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.linkage_clustering_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{8}"),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.linkage_clustering_c.restype = None

    tox.linkage_clustering_c(
        distances,
        ctypes.byref(ctypes.c_int(n_points)),
        merge_i,
        merge_j,
        heights,
        cluster_sizes,
        method,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    distances.setflags(write=False)
    merge_i.setflags(write=False)
    merge_j.setflags(write=False)
    heights.setflags(write=False)
    cluster_sizes.setflags(write=False)

    return {
        "merge_i": merge_i,
        "merge_j": merge_j,
        "heights": heights,
        "cluster_sizes": cluster_sizes
    }
