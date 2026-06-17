from error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def save_tox_data(
        zip_filename,
        gene_ids,
        gene_ids_file,
        expression,
        expression_file,
        gene_to_family,
        gene_to_family_file,
        family_ids,
        family_ids_file,
        family_centroids,
        family_centroids_file,
        shift_vectors,
        shift_vectors_file
        ):
    """
    Save standard tox data

    Parameters
    ----------
    zip_filename : str
        Zip filename
    gene_ids : ndarray[f"S{gene_ids_strlen}"] of shape (n_gene_ids_elements,) in column-major layout (order='F'), optional
        Gene ids array, will be saved if provided
        M_DOC_NO_DEFAULT
    gene_ids_file : str, optional
        Name of the gene ids file
    expression : ndarray[np.float64] of shape (n_expression_elements_dim_1, n_expression_elements_dim_2) in column-major layout (order='F'), optional
        Expression vectors array, will be saved if provided
    expression_file : str, optional
        Name of the expression file
    gene_to_family : ndarray[np.int32] of shape (n_gene_to_family_elements,) in column-major layout (order='F'), optional
        Gene to family mapping array, will be saved if provided
    gene_to_family_file : str, optional
        Name of the gene to family mapping file
    family_ids : ndarray[f"S{family_ids_strlen}"] of shape (n_family_ids_elements,) in column-major layout (order='F'), optional
        Family ids array, will be saved if provided
    family_ids_file : str, optional
        Name of the family ids file
    family_centroids : ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2) in column-major layout (order='F'), optional
        Family centroids array, will be saved if provided
    family_centroids_file : str, optional
        Name of the family centroids file
    shift_vectors : ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2) in column-major layout (order='F'), optional
        Shift vectors array, will be saved if provided
    shift_vectors_file : str, optional
        Name of the shift vectors file

    Returns
    -------
    None
    """

    # ensure all array inputs are numpy arrays
    zip_filename = np.asarray(zip_filename)
    gene_ids = np.asarray(gene_ids)
    gene_ids_file = np.asarray(gene_ids_file)
    expression = np.asfortranarray(expression, dtype=np.float64)
    expression_file = np.asarray(expression_file)
    gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    gene_to_family_file = np.asarray(gene_to_family_file)
    family_ids = np.asarray(family_ids)
    family_ids_file = np.asarray(family_ids_file)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    family_centroids_file = np.asarray(family_centroids_file)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    shift_vectors_file = np.asarray(shift_vectors_file)

    # extract dimension arguments
    zip_filename_strlen = zip_filename.dtype.itemsize // zip_filename.dtype.alignment
    gene_ids_strlen = gene_ids.dtype.itemsize // gene_ids.dtype.alignment
    n_gene_ids_elements = gene_ids.shape[0]
    gene_ids_file_strlen = gene_ids_file.dtype.itemsize // gene_ids_file.dtype.alignment
    n_expression_elements_dim_1 = expression.shape[0]
    n_expression_elements_dim_2 = expression.shape[1]
    expression_file_strlen = expression_file.dtype.itemsize // expression_file.dtype.alignment
    n_gene_to_family_elements = gene_to_family.shape[0]
    gene_to_family_file_strlen = gene_to_family_file.dtype.itemsize // gene_to_family_file.dtype.alignment
    family_ids_strlen = family_ids.dtype.itemsize // family_ids.dtype.alignment
    n_family_ids_elements = family_ids.shape[0]
    family_ids_file_strlen = family_ids_file.dtype.itemsize // family_ids_file.dtype.alignment
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    family_centroids_file_strlen = family_centroids_file.dtype.itemsize // family_centroids_file.dtype.alignment
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]
    shift_vectors_file_strlen = shift_vectors_file.dtype.itemsize // shift_vectors_file.dtype.alignment

    # Create temporaries and/or outputs
    zip_filename = zip_filename.astype(f"S{zip_filename_strlen}", order="F")
    ierr = ctypes.c_int(0)
    gene_ids = gene_ids.astype(f"S{gene_ids_strlen}", order="F")
    gene_ids_file = gene_ids_file.astype(f"S{gene_ids_file_strlen}", order="F")
    expression_file = expression_file.astype(f"S{expression_file_strlen}", order="F")
    gene_to_family_file = gene_to_family_file.astype(f"S{gene_to_family_file_strlen}", order="F")
    family_ids = family_ids.astype(f"S{family_ids_strlen}", order="F")
    family_ids_file = family_ids_file.astype(f"S{family_ids_file_strlen}", order="F")
    family_centroids_file = family_centroids_file.astype(f"S{family_centroids_file_strlen}", order="F")
    shift_vectors_file = shift_vectors_file.astype(f"S{shift_vectors_file_strlen}", order="F")

    # define ctypes interface
    tox.save_tox_data_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{zip_filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{gene_ids_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{gene_ids_file_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{expression_file_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{gene_to_family_file_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{family_ids_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{family_ids_file_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{family_centroids_file_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{shift_vectors_file_strlen}"),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.save_tox_data_c.restype = None

    tox.save_tox_data_c(
        zip_filename,
        ctypes.byref(ctypes.c_int(zip_filename_strlen)),
        ctypes.byref(ierr),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        gene_ids_file,
        ctypes.byref(ctypes.c_int(gene_ids_file_strlen)),
        expression,
        ctypes.byref(ctypes.c_int(n_expression_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_elements_dim_2)),
        expression_file,
        ctypes.byref(ctypes.c_int(expression_file_strlen)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_gene_to_family_elements)),
        gene_to_family_file,
        ctypes.byref(ctypes.c_int(gene_to_family_file_strlen)),
        family_ids,
        ctypes.byref(ctypes.c_int(family_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_family_ids_elements)),
        family_ids_file,
        ctypes.byref(ctypes.c_int(family_ids_file_strlen)),
        family_centroids,
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_2)),
        family_centroids_file,
        ctypes.byref(ctypes.c_int(family_centroids_file_strlen)),
        shift_vectors,
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_2)),
        shift_vectors_file,
        ctypes.byref(ctypes.c_int(shift_vectors_file_strlen))
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return None


def cluster_factor_trajectories_k_means(
        trajectories,
        centroids,
        max_iterations
        ):
    """
    Performs k-means clustering on factor trajectories, so factor evolution over time

    Parameters
    ----------
    trajectories : ndarray[np.float64] of shape (n_factors, n_samples, n_timepoints) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : ndarray[np.float64] of shape (n_factors, n_clusters) in column-major layout (order='F'), modified in-place
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.
        The final values will be the final centroids of the clusters
    max_iterations : int
        number of maximum iterations of the clustering

    Returns
    -------
    results : dict
        labels : ndarray[np.int32] of shape (n_samples * n_timepoints,) in column-major layout (order='F')
            array of labels, each index corresponds to the respective point's index, so first label is first point's label.
            each label is the index of its related cluster -> `1<=label<=n_clusters=k`,
        label_counts : ndarray[np.int32] of shape (n_clusters,) in column-major layout (order='F')
            holds the number of points having the respective label assigned
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


def read_expression_vectors_tsv(
        file_list,
        gene_ids,
        expression_vectors,
        n_header_rows,
        gene_col,
        value_cols,
        start_row,
        delimiter
        ):
    """
    Read expression vectors from csv/tsv files

    Parameters
    ----------
    file_list : ndarray[f"S{file_list_strlen}"] of shape (n_file_list_elements,) in column-major layout (order='F')
        List of files to read from
    gene_ids : ndarray[f"S{gene_ids_strlen}"] of shape (n_gene_ids_elements,) in column-major layout (order='F')
        Array of gene IDS
    expression_vectors : ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2) in column-major layout (order='F'), modified in-place
        Array of expression vectors
    n_header_rows : int
        Number of header rows to skip
    gene_col : int
        Index of column with gene_ids
    value_cols : ndarray[np.int32] of shape (n_value_cols_elements,) in column-major layout (order='F')
        Indicies of columns containing values
    start_row : int
        Row in the expression vectors to start in
    delimiter : str, optional
        optional delimiter
        M_DOC_DEFAULT('\t')

    Returns
    -------
    None
    """

    # ensure all array inputs are numpy arrays
    file_list = np.asarray(file_list)
    gene_ids = np.asarray(gene_ids)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    value_cols = np.ascontiguousarray(value_cols, dtype=np.int32)
    delimiter = np.asarray(delimiter)

    # extract dimension arguments
    file_list_strlen = file_list.dtype.itemsize // file_list.dtype.alignment
    n_file_list_elements = file_list.shape[0]
    gene_ids_strlen = gene_ids.dtype.itemsize // gene_ids.dtype.alignment
    n_gene_ids_elements = gene_ids.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_value_cols_elements = value_cols.shape[0]

    # Create temporaries and/or outputs
    file_list = file_list.astype(f"S{file_list_strlen}", order="F")
    gene_ids = gene_ids.astype(f"S{gene_ids_strlen}", order="F")
    n_header_rows = ctypes.c_int(n_header_rows)
    gene_col = ctypes.c_int(gene_col)
    start_row = ctypes.c_int(start_row)
    ierr = ctypes.c_int(0)
    delimiter = delimiter.astype(f"S{1}", order="F")

    # define ctypes interface
    tox.read_expression_vectors_tsv_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{file_list_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{gene_ids_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{1}")
    )
    tox.read_expression_vectors_tsv_c.restype = None

    tox.read_expression_vectors_tsv_c(
        file_list,
        ctypes.byref(ctypes.c_int(file_list_strlen)),
        ctypes.byref(ctypes.c_int(n_file_list_elements)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_2)),
        ctypes.byref(ctypes.c_int(n_header_rows)),
        ctypes.byref(ctypes.c_int(gene_col)),
        value_cols,
        ctypes.byref(ctypes.c_int(n_value_cols_elements)),
        ctypes.byref(ctypes.c_int(start_row)),
        ctypes.byref(ierr),
        delimiter
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    expression_vectors.setflags(write=False)

    return None


def mean_vector(
        expression_vectors,
        gene_indices
        ):
    """
    Computes the element-wise mean for a given set of vectors.

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
    Iterates over families, filters gene indices, and computes centroids.

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

    Returns
    -------
    results : dict
        centroid_matrix : ndarray[np.float64] of shape (n_axes, n_families) in column-major layout (order='F')
            The output matrix (n_axes x n_families) to store the computed centroids.,
        tmp_selected_indices : ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
            An output array for storing indices.
    """

    # ensure all array inputs are numpy arrays
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    mode = np.asarray(mode)
    if ortholog_set is None:
        ortholog_set = [0] * expression_vectors.shape[1]
    ortholog_set = np.ascontiguousarray(ortholog_set, dtype=np.int32)
    print(ortholog_set)

    # extract dimension arguments
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Create temporaries and/or outputs
    centroid_matrix = np.empty((n_axes, n_families), dtype=np.float64, order='F')
    mode = mode.astype(f"S{15}", order="F")
    tmp_selected_indices = np.empty((n_genes,), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
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
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32)
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
    tmp_selected_indices.setflags(write=False)

    return centroid_matrix


def mask_get_first_successor_idx(
        bit_mask
        ):
    """
    Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.

    Parameters
    ----------
    bit_mask : ndarray[np.int32] of shape (n_bit_mask_elements,) in column-major layout (order='F')
        chunked mask to mark active genes

    Returns
    -------
    idx : int
        index of last active gene
    """

    # ensure all array inputs are numpy arrays
    bit_mask = np.ascontiguousarray(bit_mask, dtype=np.int32)

    # extract dimension arguments
    n_bit_mask_elements = bit_mask.shape[0]

    # Create temporaries and/or outputs
    idx = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    # define ctypes interface
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


def detect_dosage_effect(
        ancestor,
        genes,
        filtered_paralogs_mask,
        max_subset_size,
        n_paralog_subsets,
        max_angle,
        gain_gamma
        ):
    """
    Identifies subsets of paralogs with small angle to the `ancestor` (max_angle) and sum to a magnitude significantly exceeding `norm(ancestor)` (gain)

    Parameters
    ----------
    ancestor : ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : ndarray[np.float64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    filtered_paralogs_mask : ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
        It is recommended to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.
    max_subset_size : int
        maximum subset size of checked gene subsets.
        It is *VERY IMPORTANT* to compute this argument using [[tox_paralog_analysis(module):calc_work_arr_paralog_subsets_size]]'s output `max_subset_size`.
    max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
    gain_gamma : float, optional
        positive magnitude gain for dosage effect, default 0.1

    Returns
    -------
    results : dict
        n_results : int
            number of resulting subsets. They are stored as the first `n_results` elements of `work_arr_paralog_subsets`,
        work_arr_paralog_subsets : ndarray[np.int32] of shape (n_mask_chunks, n_paralog_subsets) in column-major layout (order='F')
            working array to hold bitmask encoded subsets for detection.
            @note
            Each bitmask is built of 32 bit chunks. `(n_genes + 31) / 32` is equivalent to `ceil(n_genes / 32)` and represents the number of chunks
            @endnote,
        active_mask : ndarray[np.int32] of shape (n_mask_chunks,) in column-major layout (order='F')
            working array to hold the extended subsets,
        temp_paralog_vector : ndarray[np.float64] of shape (n_dims,) in column-major layout (order='F')
            vector used for pruning subsets
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.float64)
    genes = np.asfortranarray(genes, dtype=np.float64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]

    # Create temporaries and/or outputs
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks,), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)
    max_angle = ctypes.c_double(max_angle)
    gain_gamma = ctypes.c_double(gain_gamma)

    # define ctypes interface
    tox.detect_dosage_effect_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.float64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double)
    )
    tox.detect_dosage_effect_c.restype = None

    tox.detect_dosage_effect_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(ctypes.c_int(max_subset_size)),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_double(max_angle)),
        ctypes.byref(ctypes.c_double(gain_gamma))
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)

    return {
        "n_results": n_results,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector
    }


def filter_paralogs_by_pattern_dosage_effect(
        gene_angles,
        threshold,
        n_families,
        gene_to_fam,
        n_mask_chunks
        ):
    """
    This subroutine prefilters the genes for dosage effect,
    as genes that are already too distant in angle to the ancestor don't match the pattern and don't need to be tried as subset extensions.

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
    The `detect_*` subroutines need a work array for the to be tested subsets.
    In worst case, all need to be tried and subsets that cannot be extended will be kept as results.
    This is the reason why the work array holds the results as well, as all subsets that are stored in the array can be results as well.
    This subroutine calculates the needed size for the work array.

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
        It is recommended to compute this argument using [[tox_paralog_analysis(module):filter_paralogs_by_pattern]]'s output `masks(:, family_idx)`.

    Returns
    -------
    work_array_size : int
        The calculated needed work array size in absolute worst case scenario. Look into source for details.
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
