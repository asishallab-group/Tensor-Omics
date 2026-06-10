from error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def deserialize_char_nd(
        filename
        ):
    """
    Subroutine to deserialize a flat character array from a file

    Parameters
    ----------
    filename : str
        Name of the file to read
    """

    # ensure all array inputs are numpy arrays
    flat = np.asarray(flat)
    filename = np.asarray(filename)

    # extract dimension arguments
    flat_strlen = flat.dtype.itemsize // flat.dtype.alignment
    n_flat_elements = flat.shape[0]
    filename_strlen = filename.dtype.itemsize // filename.dtype.alignment

    # Create temporaries and/or outputs
    flat = np.zeros((n_flat_elements), dtype=f"S{flat_strlen}", order='F')
    filename = filename.astype(f"S{filename_strlen}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    tox.deserialize_char_nd_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{flat_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.deserialize_char_nd_c.restype = None

    tox.deserialize_char_nd_c(
        flat,
        ctypes.byref(ctypes.c_int(flat_strlen)),
        ctypes.byref(ctypes.c_int(n_flat_elements)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    flat.setflags(write=False)

    return flat


def serialize_char_nd(
        flat,
        filename
        ):
    """
    Serialize a character array of arbitrary dimensions to a binary file.
    The file will contain a magic number, type code, dimension, shape, character length, and the array data.
    @note This routine is only called by R and serializes only flat character arrays to the memory

    Parameters
    ----------
    flat : ndarray[f"S{flat_strlen}"] of shape (flat_strlen, *) in column-major layout (order='F')
        flat array to save
    flat_shape : ndarray[np.int32] of shape (n_flat_shape_elements) in column-major layout (order='F')
        dimensions of the array
    filename : str
        output filename
    """

    # ensure all array inputs are numpy arrays
    flat = np.asarray(flat)
    flat_shape = np.ascontiguousarray(flat.shape, dtype=np.int32)
    filename = np.asarray(filename)

    # extract dimension arguments
    flat_strlen = flat.dtype.itemsize // flat.dtype.alignment
    n_flat_shape_elements = flat_shape.shape[0]
    filename_strlen = filename.dtype.itemsize // filename.dtype.alignment

    # Create temporaries and/or outputs
    flat = flat.astype(f"S{flat_strlen}", order="F")
    filename = filename.astype(f"S{filename_strlen}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    tox.serialize_char_nd_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='F_CONTIGUOUS', dtype=f"S{flat_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.serialize_char_nd_c.restype = None

    tox.serialize_char_nd_c(
        flat,
        ctypes.byref(ctypes.c_int(flat_strlen)),
        flat_shape,
        ctypes.byref(ctypes.c_int(n_flat_shape_elements)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return None


def serialize_int_3d(
        arr,
        filename
        ):
    """
    Serialize a 3D integer(int32) array to a binary file.
    The file will contain a magic number, type code, dimension, shape, and the array data.

    Parameters
    ----------
    arr : ndarray[np.int32] of shape (n_arr_elements_dim_1, n_arr_elements_dim_2, n_arr_elements_dim_3) in column-major layout (order='F')
        array to save
    filename : str
        output filename
    """

    # ensure all array inputs are numpy arrays
    arr = np.asfortranarray(arr, dtype=np.int32)
    filename = np.asarray(filename)

    # extract dimension arguments
    n_arr_elements_dim_1 = arr.shape[0]
    n_arr_elements_dim_2 = arr.shape[1]
    n_arr_elements_dim_3 = arr.shape[2]
    filename_strlen = filename.dtype.itemsize // filename.dtype.alignment

    # Create temporaries and/or outputs
    filename = filename.astype(f"S{filename_strlen}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    tox.serialize_int_3d_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=3, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.serialize_int_3d_c.restype = None

    tox.serialize_int_3d_c(
        arr,
        ctypes.byref(ctypes.c_int(n_arr_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_arr_elements_dim_2)),
        ctypes.byref(ctypes.c_int(n_arr_elements_dim_3)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return None


def serialize_logical_3d(
        arr,
        filename
        ):
    """
    Serialize a 3D logical array to a binary file.
    The file will contain a magic number, type code, dimension, shape, and the array data.

    Parameters
    ----------
    arr : ndarray[np.int32] of shape (n_arr_elements_dim_1, n_arr_elements_dim_2, n_arr_elements_dim_3) in column-major layout (order='F')
        array to save
    filename : str
        output filename
    """

    # ensure all array inputs are numpy arrays
    arr = np.asfortranarray(arr, dtype=np.int32)
    filename = np.asarray(filename)

    # extract dimension arguments
    n_arr_elements_dim_1 = arr.shape[0]
    n_arr_elements_dim_2 = arr.shape[1]
    n_arr_elements_dim_3 = arr.shape[2]
    filename_strlen = filename.dtype.itemsize // filename.dtype.alignment

    # Create temporaries and/or outputs
    filename = filename.astype(f"S{filename_strlen}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    tox.serialize_logical_3d_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=3, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.serialize_logical_3d_c.restype = None

    tox.serialize_logical_3d_c(
        arr,
        ctypes.byref(ctypes.c_int(n_arr_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_arr_elements_dim_2)),
        ctypes.byref(ctypes.c_int(n_arr_elements_dim_3)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return None


def serialize_real_nd(
        arr,
        dims,
        ndim,
        filename
        ):
    """
    Writes serialized real array from R to file with metdata.

    Parameters
    ----------
    arr : ndarray[np.real64] of shape (n_arr_elements) in column-major layout (order='F')
        array to save
    dims : ndarray[np.int32] of shape (n_dims_elements) in column-major layout (order='F')
        Dimensions of the array
    ndim : int
        Number of dimensions
    filename : str
        filename
    """

    # ensure all array inputs are numpy arrays
    arr = np.ascontiguousarray(arr, dtype=np.real64)
    dims = np.ascontiguousarray(dims, dtype=np.int32)
    filename = np.asarray(filename)

    # extract dimension arguments
    n_arr_elements = arr.shape[0]
    n_dims_elements = dims.shape[0]
    filename_strlen = filename.dtype.itemsize // filename.dtype.alignment

    # Create temporaries and/or outputs
    ndim = ctypes.c_int(ndim)
    filename = filename.astype(f"S{filename_strlen}", order="F")
    ierr = ctypes.c_int(0)

    # define ctypes interface
    tox.serialize_real_nd_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.real64),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{filename_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.serialize_real_nd_c.restype = None

    tox.serialize_real_nd_c(
        arr,
        ctypes.byref(ctypes.c_int(n_arr_elements)),
        dims,
        ctypes.byref(ctypes.c_int(n_dims_elements)),
        ctypes.byref(ctypes.c_int(ndim)),
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only


    return None


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
    gene_ids : ndarray[f"S{gene_ids_strlen}"] of shape (gene_ids_strlen, n_gene_ids_elements) in column-major layout (order='F'), optional
        Gene ids array, will be saved if provided
        This argument will be ignored if not present.
    gene_ids_file : str, optional
        Name of the gene ids file
    expression : ndarray[np.real64] of shape (n_expression_elements_dim_1, n_expression_elements_dim_2) in column-major layout (order='F'), optional
        Expression vectors array, will be saved if provided
    expression_file : str, optional
        Name of the expression file
    gene_to_family : ndarray[np.int32] of shape (n_gene_to_family_elements) in column-major layout (order='F'), optional
        Gene to family mapping array, will be saved if provided
    gene_to_family_file : str, optional
        Name of the gene to family mapping file
    family_ids : ndarray[f"S{family_ids_strlen}"] of shape (family_ids_strlen, n_family_ids_elements) in column-major layout (order='F'), optional
        Family ids array, will be saved if provided
    family_ids_file : str, optional
        Name of the family ids file
    family_centroids : ndarray[np.real64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2) in column-major layout (order='F'), optional
        Family centroids array, will be saved if provided
    family_centroids_file : str, optional
        Name of the family centroids file
    shift_vectors : ndarray[np.real64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2) in column-major layout (order='F'), optional
        Shift vectors array, will be saved if provided
    shift_vectors_file : str, optional
        Name of the shift vectors file
    """

    # ensure all array inputs are numpy arrays
    zip_filename = np.asarray(zip_filename)
    gene_ids = np.asarray(gene_ids)
    gene_ids_file = np.asarray(gene_ids_file)
    expression = np.asfortranarray(expression, dtype=np.real64)
    expression_file = np.asarray(expression_file)
    gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    gene_to_family_file = np.asarray(gene_to_family_file)
    family_ids = np.asarray(family_ids)
    family_ids_file = np.asarray(family_ids_file)
    family_centroids = np.asfortranarray(family_centroids, dtype=np.real64)
    family_centroids_file = np.asarray(family_centroids_file)
    shift_vectors = np.asfortranarray(shift_vectors, dtype=np.real64)
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
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.real64),
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
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.real64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=0, flags='C_CONTIGUOUS', dtype=f"S{family_centroids_file_strlen}"),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.real64),
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
    trajectories : ndarray[np.real64] of shape (n_factors, n_samples, n_timepoints) in column-major layout (order='F')
        matrix with data points to cluster
    centroids : ndarray[np.real64] of shape (n_factors, n_clusters) in column-major layout (order='F')
        matrix with initial centroids of the clusters, could be random data or actual points or unassigned garbage.
        The centroids should be unique. This is not checked in this routine.
        The final values will be the final centroids of the clusters
    max_iterations : int
        number of maximum iterations of the clustering
    """

    # ensure all array inputs are numpy arrays
    trajectories = np.asfortranarray(trajectories, dtype=np.real64)
    centroids = np.asfortranarray(centroids, dtype=np.real64)
    labels = np.ascontiguousarray(labels, dtype=np.int32)
    label_counts = np.ascontiguousarray(label_counts, dtype=np.int32)

    # extract dimension arguments
    n_clusters = centroids.shape[1]
    n_factors = trajectories.shape[0]
    n_samples = trajectories.shape[1]
    n_timepoints = trajectories.shape[2]

    # Create temporaries and/or outputs
    labels = np.empty((n_samples * n_timepoints), dtype=np.int32, order='F')
    label_counts = np.empty((n_clusters), dtype=np.int32, order='F')
    ierr = ctypes.c_int(0)
    max_iterations = ctypes.c_int(max_iterations)

    # define ctypes interface
    tox.cluster_factor_trajectories_k_means_c.argtypes = (
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=3, flags='F_CONTIGUOUS', dtype=np.real64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.real64),
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
        "centroids": centroids,
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
        delimiter='\t'
        ):
    """
    Read expression vectors from csv/tsv files

    Parameters
    ----------
    file_list : ndarray[f"S{file_list_strlen}"] of shape (file_list_strlen, n_file_list_elements) in column-major layout (order='F')
        List of files to read from
    gene_ids : ndarray[f"S{gene_ids_strlen}"] of shape (gene_ids_strlen, n_gene_ids_elements) in column-major layout (order='F')
        Array of gene IDS
    expression_vectors : ndarray[np.real64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2) in column-major layout (order='F')
        Array of expression vectors
    n_header_rows : int
        Number of header rows to skip
    gene_col : int
        Index of column with gene_ids
    value_cols : ndarray[np.int32] of shape (n_value_cols_elements) in column-major layout (order='F')
        Indicies of columns containing values
    start_row : int
        Row in the expression vectors to start in
    delimiter : str, optional
        optional delimiter
        The default value is `'\t'`.
    """

    # ensure all array inputs are numpy arrays
    file_list = np.asarray(file_list)
    gene_ids = np.asarray(gene_ids)
    expression_vectors = np.asfortranarray(expression_vectors, dtype=np.real64)
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
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.real64),
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

    return expression_vectors


def mask_get_first_successor_idx(
        bit_mask
        ):
    """
    Helper function that returns the index after the last active gene in `bit_mask`, so the first succeeding gene.

    Parameters
    ----------
    bit_mask : ndarray[np.int32] of shape (n_bit_mask_elements) in column-major layout (order='F')
        chunked mask to mark active genes
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


def detect_patterns(
        ancestor,
        genes,
        pattern,
        filtered_paralogs_mask,
        max_subset_size,
        dosage_max_angle,
        dosage_gain_gamma,
        subfunc_rdi_threshold,
        subfunc_paralog_norms,
        subfunc_sorted_paralog_norms_perm
        ):
    """
    Identifies subsets of paralogs where dosage effect or subfunctionalization applies, depending on `pattern`

    Parameters
    ----------
    ancestor : ndarray[np.real64] of shape (n_dims) in column-major layout (order='F')
        expression vector of ancestral ortholog
    genes : ndarray[np.real64] of shape (n_dims, n_genes) in column-major layout (order='F')
        expression vectors of genes
    pattern : int
        used pattern for detection
        |       Pattern        | Value |
        |----------------------|-------|
        |    Dosage Effect     |   0   |
        | Subfunctionalization |   1   |
    filtered_paralogs_mask : ndarray[np.int32] of shape (n_mask_chunks) in column-major layout (order='F')
        bit mask with genes' indices kept by pattern set to 1, else 0. Use `filter_paralogs_by_pattern` for its calculation
    max_subset_size : int
        maximum subset size of checked gene subsets. ***USE `calc_work_arr_paralog_subsets_size` TO DETERMINE THIS NUMBER***
    dosage_max_angle : float, optional
        in dosage mode maximum angle in radians `0<=angle<=Pi` that a subset candidate must not exceed, otherwise pruned, default is Pi
    dosage_gain_gamma : float, optional
        in dosage mode required positive magnitude gain for dosage, default 0.1
    subfunc_rdi_threshold : float, optional
        max allowed residual distance from `ancestor`
    subfunc_paralog_norms : ndarray[np.real64] of shape (n_genes) in column-major layout (order='F'), optional
        in subfunctionalization mode needed for subset pruning, holds the euclidean norms of genes (you can use the `norm` from `f42_utils` function for this)
    subfunc_sorted_paralog_norms_perm : ndarray[np.int32] of shape (n_genes) in column-major layout (order='F'), optional
        in subfunctionalization mode needed for subset pruning, as the minimum norm of the genes that could extend a subset should not be lower than the subset angle to the ancestor
        This optional argument needs to be passed if used mode is [[tox_paralog_analysis(module):MODE_SUBFUNC_PATTERN(variable)]].
    """

    # ensure all array inputs are numpy arrays
    ancestor = np.ascontiguousarray(ancestor, dtype=np.real64)
    genes = np.asfortranarray(genes, dtype=np.real64)
    filtered_paralogs_mask = np.ascontiguousarray(filtered_paralogs_mask, dtype=np.int32)
    work_arr_paralog_subsets = np.asfortranarray(work_arr_paralog_subsets, dtype=np.int32)
    active_mask = np.ascontiguousarray(active_mask, dtype=np.int32)
    temp_paralog_vector = np.ascontiguousarray(temp_paralog_vector, dtype=np.real64)
    subfunc_paralog_norms = np.ascontiguousarray(subfunc_paralog_norms, dtype=np.real64)
    subfunc_sorted_paralog_norms_perm = np.ascontiguousarray(subfunc_sorted_paralog_norms_perm, dtype=np.int32)
    subfunc_temp_work_array = np.ascontiguousarray(subfunc_temp_work_array, dtype=np.real64)

    # extract dimension arguments
    n_genes = genes.shape[1]
    n_dims = ancestor.shape[0]
    n_mask_chunks = filtered_paralogs_mask.shape[0]
    n_paralog_subsets = work_arr_paralog_subsets.shape[1]

    # Create temporaries and/or outputs
    pattern = ctypes.c_int(pattern)
    n_results = ctypes.c_int(0)
    max_subset_size = ctypes.c_int(max_subset_size)
    work_arr_paralog_subsets = np.empty((n_mask_chunks, n_paralog_subsets), dtype=np.int32, order='F')
    active_mask = np.empty((n_mask_chunks), dtype=np.int32, order='F')
    temp_paralog_vector = np.empty((n_dims), dtype=np.real64, order='F')
    dosage_max_angle = ctypes.c_double(dosage_max_angle)
    dosage_gain_gamma = ctypes.c_double(dosage_gain_gamma)
    subfunc_rdi_threshold = ctypes.c_double(subfunc_rdi_threshold)
    subfunc_temp_work_array = np.empty((n_genes), dtype=np.real64, order='F')
    ierr = ctypes.c_int(0)

    # define ctypes interface
    tox.detect_patterns_c.argtypes = (
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.real64),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.real64),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=2, flags='F_CONTIGUOUS', dtype=np.int32),
        ctypes.POINTER(ctypes.c_int),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.real64),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        ctypes.POINTER(ctypes.c_double),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.real64),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.int32),
        np.ctypeslib.ndpointer(ndim=1, flags='C_CONTIGUOUS', dtype=np.real64),
        ctypes.POINTER(ctypes.c_int)
    )
    tox.detect_patterns_c.restype = None

    tox.detect_patterns_c(
        ancestor,
        genes,
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_dims)),
        ctypes.byref(ctypes.c_int(pattern)),
        filtered_paralogs_mask,
        ctypes.byref(ctypes.c_int(n_mask_chunks)),
        ctypes.byref(n_results),
        ctypes.byref(ctypes.c_int(max_subset_size)),
        work_arr_paralog_subsets,
        ctypes.byref(ctypes.c_int(n_paralog_subsets)),
        active_mask,
        temp_paralog_vector,
        ctypes.byref(ctypes.c_double(dosage_max_angle)),
        ctypes.byref(ctypes.c_double(dosage_gain_gamma)),
        ctypes.byref(ctypes.c_double(subfunc_rdi_threshold)),
        subfunc_paralog_norms,
        subfunc_sorted_paralog_norms_perm,
        subfunc_temp_work_array,
        ctypes.byref(ierr)
    )

    # throw error on error
    check_err_code(ierr.value)

    # Mark all arrays as read-only
    work_arr_paralog_subsets.setflags(write=False)
    active_mask.setflags(write=False)
    temp_paralog_vector.setflags(write=False)
    subfunc_temp_work_array.setflags(write=False)

    return {
        "n_results": n_results,
        "work_arr_paralog_subsets": work_arr_paralog_subsets,
        "active_mask": active_mask,
        "temp_paralog_vector": temp_paralog_vector,
        "subfunc_temp_work_array": subfunc_temp_work_array
    }
