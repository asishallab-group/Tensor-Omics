from .error_handling import check_err_code

import numpy as np
import ctypes
import os

# Load library
dll_path = os.path.abspath("build/libtensor-omics.so")
tox = ctypes.CDLL(dll_path)


def mean_vector(
        expression_vectors,
        gene_indices
        ):
    """
    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_indices : np.ndarray[np.int32] of shape (n_selected_genes,) in column-major layout (order='F')
        An array containing the column indices of the selected genes in 'expression_vectors'.

    Returns
    -------
    centroid : np.ndarray[np.float64] of shape (n_axes,) in column-major layout (order='F')
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
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes) in column-major layout (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
        An array mapping each gene (by index) to a family ID.
    mode : str
        used mode for grouping|       Mode      |       Value       |
        |-----------------|-------------------|
        | Group Orthologs | "group_orthologs" |
        |    Group all    |    "group_all"    |
    ortholog_set : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F'), optional
        A logical array indicating if a gene is part of a specific subset (e.g., orthologs).This optional argument needs to be passed if used mode (`mode`) is [[tox_gene_centroids(module):MODE_GROUP_ORTHOLOGS(variable)]].

    Returns
    -------
    results : dict
        centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families) in column-major layout (order='F')
            The output matrix (n_axes x n_families) to store the computed centroids.,
        tmp_selected_indices : np.ndarray[np.int32] of shape (n_genes,) in column-major layout (order='F')
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
