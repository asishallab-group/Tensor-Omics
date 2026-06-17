from error_handling import check_err_code

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
    trajectories : ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : ndarray[np.float64] of shape (n_factors, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.
        The final values will be the final centroids of the clusters
    max_iterations : int, optional
        number of maximum iterations of the clustering
        The default value is `300_int32`.

    Returns
    -------
    results : dict
        labels : ndarray[np.int32] of shape (n_samples * n_timepoints,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned

    Notes
    -----
    Performs k-means clustering on factor trajectories, so factor evolution over time
    """

    # ensure all array inputs are numpy arrays
    trajectories = np.asfortranarray(trajectories, dtype=np.float64)
    centroids = np.asfortranarray(centroids, dtype=np.float64)

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
        ctypes.byref(ctypes.c_int(max_iterations))
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
    data_points : ndarray[np.float64] of shape (n_dims, n_points) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : ndarray[np.float64] of shape (n_dims, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.
        The final values will be the final centroids of the clusters
    max_iterations : int, optional
        number of maximum iterations of the clustering
        The default value is `300_int32`.

    Returns
    -------
    results : dict
        labels : ndarray[np.int32] of shape (n_points,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned

    Notes
    -----
    k-means clustering algorithm:
    1. Assigns each data point to one of `k` clusters whose centroid is clostest
    2. Recalculates the centroids using the mean of its assigned points
    3. repeat 1-2 until assignment remains unchanged
    """

    # ensure all array inputs are numpy arrays
    data_points = np.asfortranarray(data_points, dtype=np.float64)
    centroids = np.asfortranarray(centroids, dtype=np.float64)

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
        ctypes.byref(ctypes.c_int(max_iterations))
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
    distances : ndarray[np.float64] of shape (n_points, n_points) in column-major layout (order='F'), modified in-place
        symmetric distance matrix, holding the positive distances between points. Distance of X->X is always zero.
        @note
        This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
        So there is no need to copy an existing distance matrix, just pass the original.
        @endnote
    method : str
        used algorithm
        |      Method      |                      Value                     |
        |------------------|------------------------------------------------|
        | Average / UPGMA  |    "average"      |
        | Weighted / WPGMA |    "weighted"     |
        |      Ward        |    "ward"         |

    Returns
    -------
    results : dict
        merge_i : ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes,
        merge_j : ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            holds cluster labels of the merged node pair at iteration k -> positives relate to leafs/data point indices, negatives to inner nodes,
        heights : ndarray[np.float64] of shape (n_points - 1,) in column-major layout (order='F')
            height of the shorter branch of the merge, e.g. if (A,B)+(C) merges to ((A,B),C), the branch to (A,B) is shorter,
        cluster_sizes : ndarray[np.int32] of shape (n_points - 1,) in column-major layout (order='F')
            size of cluster at iteration k

    Notes
    -----
    @note
    This subroutine operates in-place in the bottom triangle of the distance matrix and recovers it using the top triangle once done or on error.
    So there is no need to copy an existing distance matrix, just pass the original.
    @endnote
    """

    # ensure all array inputs are numpy arrays
    distances = np.asfortranarray(distances, dtype=np.float64)
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


def mean_vector(
        expression_vectors,
        gene_indices
        ):
    """
    Parameters
    ----------
    expression_vectors : ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_indices : ndarray[np.int32] of shape (n_selected_genes,) in column-major layout (order='F')
        An array containing the column indices of the selected genes in 'expression_vectors'.

    Returns
    -------
    centroid : ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
        The output vector representing the computed centroid.

    Notes
    -----
    Computes the element-wise mean for a given set of vectors.
    """

    # ensure all array inputs are numpy arrays
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    gene_indices = np.ascontiguousarray(gene_indices, dtype=np.int32)

    # extract dimension arguments
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_selected_genes = gene_indices.shape[0]

    # Create temporaries and/or outputs
    centroid = np.empty((n_axes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mean_vector_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mean_vector_c.restype = None

    tox.mean_vector_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_indices,
        ctypes.byref(ctypes.c_int(n_selected_genes)),
        centroid,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroid.setflags(write=False)

    return centroid


def group_centroid(
        expression_vectors,
        gene_to_family,
        n_families,
        mode,
        ortholog_set=None
        ):
    """
    Parameters
    ----------
    expression_vectors : ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        An array mapping each gene (by index) to a family ID.
    mode : str
        used mode for grouping
        |       Mode       |                             Value                               |
        |------------------|-----------------------------------------------------------------|
        | Group Orthologs  |   "group_orthologs"    |
        |    Group all     |      "group_all"       |
    ortholog_set : ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F'), optional
        A logical array indicating if a gene is part of a specific subset (e.g., orthologs).
        This optional argument needs to be passed if used mode (`mode`) is [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]].

    Returns
    -------
    results : dict
        centroid_matrix : ndarray[np.float64] of shape (n_axes, n_families) in column-major layout (order='F')
            The output matrix (n_axes x n_families) to store the computed centroids.,
        tmp_selected_indices : ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
            An output array for storing indices.

    Notes
    -----
    Iterates over families, filters gene indices, and computes centroids.
    """

    # ensure all array inputs are numpy arrays
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    mode = np.asarray(mode)
    if ortholog_set is not None:
        ortholog_set = np.ascontiguousarray(ortholog_set, dtype=np.int32)

    # extract dimension arguments
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Create temporaries and/or outputs
    centroid_matrix = np.empty((n_axes, n_families), dtype=np.float64, order='F')
    mode = mode.astype(f"S{15}", order="F")
    tmp_selected_indices = np.empty((n_genes,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.group_centroid_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{15}"),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        nullable(np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32))
    )
    tox.group_centroid_c.restype = None

    tox.group_centroid_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_families)),
        centroid_matrix,
        mode,
        tmp_selected_indices,
        ctypes.byref(ierr),
        ortholog_set
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    centroid_matrix.setflags(write=False)

    return centroid_matrix


def mask_get_first_successor_idx(
        bit_mask
        ):
    """
    Parameters
    ----------
    bit_mask : ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active genes

    Returns
    -------
    idx : int
        index of last active gene

    Notes
    -----
    Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.
    """

    # ensure all array inputs are numpy arrays
    bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)

    # extract dimension arguments
    n_bit_mask_elements = bit_mask.shape[0]

    # Create temporaries and/or outputs
    idx = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.mask_get_first_successor_idx_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.mask_get_first_successor_idx_c.restype = None

    tox.mask_get_first_successor_idx_c(
        bit_mask,
        ctypes.byref(ctypes.c_int(n_bit_mask_elements)),
        ctypes.byref(idx),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return idx


def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    Parameters
    ----------
    gene_angles : ndarray[np.float64] of shape (n_genes,) in column-major layout (order='F')
        vector, holding the angles between ancestor and genes (0<=angle<=Pi)
    threshold : float
        filter threshold
    gene_to_fam : ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        a mapping of gene index to family index, so gene i is related to `gene_angles(i)` and part of family `j=gene_to_fam(i)`.

    Returns
    -------
    masks : ndarray[np.int32] of shape (n_mask_chunks, n_families) in column-major layout (order='F')
        bit mask that will have indices of genes kept by pattern set to 1, else 0

    Notes
    -----
    This subroutine prefilters the genes for dosage effect,
    as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.
    """

    # ensure all array inputs are numpy arrays
    gene_angles = np.ascontiguousarray(gene_angles, dtype=np.float64)
    gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)

    # extract dimension arguments
    n_genes = gene_angles.shape[0]

    # Create temporaries and/or outputs
    threshold = ctypes.c_double(threshold)
    masks = np.empty((n_mask_chunks, n_families), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.filter_paralogs_by_pattern_dosage_effect_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.filter_paralogs_by_pattern_dosage_effect_c.restype = None

    tox.filter_paralogs_by_pattern_dosage_effect_c(
        gene_angles,
        ctypes.byref(ctypes.c_double(threshold)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        gene_to_fam,
        masks,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    masks.setflags(write=False)

    return masks


def calc_work_arr_paralog_subsets_size(
        max_subset_size,
        n_genes,
        filtered_paralogs_mask
        ):
    """
    Parameters
    ----------
    max_subset_size : int, modified in-place
        maximum size that a subset must not exceed.
        @warning
        If the desired size is too large and leads to an integer overflow, `max_subset_size` will be set to the maximum valid size.
        Also, size will be set to number of genes in `filtered_paralogs_mask` if larger.
        @endwarning
    n_genes : int
        number of genes
    filtered_paralogs_mask : ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        Output mask with all genes disabled that did not pass the filter
        M_DM_FROM_JUST_INFO to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.

    Returns
    -------
    work_array_size : int
        The calculated needed work array size in absolute worst case scenario. Look into source for details.

    Notes
    -----
    The `detect_*` subroutines need a work array for the to be tested subsets.
    In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
    This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
    This subroutine calculates the needed size for the work array.
    """

    # ensure all array inputs are numpy arrays
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)

    # extract dimension arguments
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # Create temporaries and/or outputs
    max_subset_size = ctypes.c_int(max_subset_size)
    n_genes = ctypes.c_int(n_genes)
    work_array_size = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.calc_work_arr_paralog_subsets_size_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.calc_work_arr_paralog_subsets_size_c.restype = None

    tox.calc_work_arr_paralog_subsets_size_c(
        ctypes.byref(ctypes.c_int(max_subset_size)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(work_array_size),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return work_array_size


def omics_vector_RAP_projection(
        vecs,
        vecs_selection_mask,
        n_selected_vecs,
        axes_selection_mask,
        n_selected_axes
        ):
    """
    Parameters
    ----------
    vecs : ndarray[np.float64] of shape (n_axes, n_vecs) in column-major layout (order='F')
        matrix with expression vectors
    vecs_selection_mask : ndarray[np.int32] of shape (n_vecs,) in column-major layout (order='F')
        `.true.` for vectors where projection is to be computed
    axes_selection_mask : ndarray[np.int32] of shape (n_axes,) in column-major layout (order='F')
        `.true.` for axes to be included in RAP

    Returns
    -------
    projections : ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs) in column-major layout (order='F')
        projected vectors

    Notes
    -----
    Project selected vectors (e.g. expression vectors) onto the RAP constructed from a selected set of axes.
    """

    # ensure all array inputs are numpy arrays
    vecs = np.asfortranarray(vecs, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.int32)

    # extract dimension arguments
    n_axes = vecs.shape[0]
    n_vecs = vecs.shape[1]

    # Create temporaries and/or outputs
    projections = np.empty((n_selected_axes, n_selected_vecs), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.omics_vector_RAP_projection_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.omics_vector_RAP_projection_c.restype = None

    tox.omics_vector_RAP_projection_c(
        vecs,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        projections,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    projections.setflags(write=False)

    return projections


def omics_field_RAP_projection(
        vecs,
        vecs_selection_mask,
        n_selected_vecs,
        axes_selection_mask,
        n_selected_axes
        ):
    """
    Parameters
    ----------
    vecs : ndarray[np.float64] of shape (2 * n_axes, n_vecs) in column-major layout (order='F')
        matrix with vector fields, first n rows mean vector origin, last n rows vector targets
    vecs_selection_mask : ndarray[np.int32] of shape (n_vecs,) in column-major layout (order='F')
        `.true.` for vectors where projection is to be computed
    axes_selection_mask : ndarray[np.int32] of shape (n_axes,) in column-major layout (order='F')
        `.true.` for axes to be included in RAP

    Returns
    -------
    projections : ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs) in column-major layout (order='F')
        projected vectors

    Notes
    -----
    Project selected vector fields (e.g. shift vectors) onto the RAP constructed from a selected set of axes.
    """

    # ensure all array inputs are numpy arrays
    vecs = np.asfortranarray(vecs, dtype=np.float64)
    vecs_selection_mask = np.ascontiguousarray(vecs_selection_mask, dtype=np.int32)
    axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.int32)

    # extract dimension arguments
    n_axes = axes_selection_mask.shape[0]
    n_vecs = vecs.shape[1]

    # Create temporaries and/or outputs
    projections = np.empty((n_selected_axes, n_selected_vecs), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.omics_field_RAP_projection_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.omics_field_RAP_projection_c.restype = None

    tox.omics_field_RAP_projection_c(
        vecs,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_vecs)),
        vecs_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        projections,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    projections.setflags(write=False)

    return projections


def project_selected_vecs_onto_rap(
        selected_vecs
        ):
    """
    Parameters
    ----------
    selected_vecs : ndarray[np.float64] of shape (n_selected_axes, n_selected_vecs) in column-major layout (order='F'), modified in-place
        matrix with vectors for selected axes

    Returns
    -------
    None

    Notes
    -----
    Projects selected vectors onto its RAP
    """

    # ensure all array inputs are numpy arrays
    selected_vecs = np.asfortranarray(selected_vecs, dtype=np.float64)

    # extract dimension arguments
    n_selected_axes = selected_vecs.shape[0]
    n_selected_vecs = selected_vecs.shape[1]

    # Create temporaries and/or outputs
    ierr = ctypes.c_int(0)

    # define ctypes interface
    def nullable(ty):
        @classmethod
        def from_param(cls, obj):
            if obj is not None:
                return ty.from_param(obj)
        return type(ty.__name__, (ty,), {'from_param': from_param})

    tox.project_selected_vecs_onto_rap_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.project_selected_vecs_onto_rap_c.restype = None

    tox.project_selected_vecs_onto_rap_c(
        selected_vecs,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        ctypes.byref(ctypes.c_int(n_selected_vecs)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    selected_vecs.setflags(write=False)

    return None
