"""Python interface to Module with Euclidean distance computation routines for tensor omics.

Generated from tox_euclidean_distance. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.euclidean_distance_c.restype = None
_lib.euclidean_distance_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_double),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_EUCLIDEAN_DISTANCE_ARGUMENTS = ("vec1", "vec2", "n_elements", "result", "ierr",)

_lib.distance_to_centroid_c.restype = None
_lib.distance_to_centroid_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_DISTANCE_TO_CENTROID_ARGUMENTS = ("n_genes", "n_families", "genes", "centroids", "gene_to_fam", "distances", "n_tissues", "ierr",)

def euclidean_distance(
        vec1,
        vec2,
):
    r"""Compute the Euclidean distance between two vectors.

    Parameters
    ----------
    vec1 : np.ndarray[np.float64] of shape (n_elements,)
        First expression vector
    vec2 : np.ndarray[np.float64] of shape (n_elements,)
        Second expression vector

    Returns
    -------
    result : float
        Output scalar distance

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_euclidean_distance::euclidean_distance`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        vec1 = np.ascontiguousarray(vec1, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vec1' must be an array of np.float64: {error}") from None
    if vec1.ndim != 1:
        raise ValueError(f"'vec1' must have 1 dimension, but has {vec1.ndim}")
    try:
        vec2 = np.ascontiguousarray(vec2, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vec2' must be an array of np.float64: {error}") from None
    if vec2.ndim != 1:
        raise ValueError(f"'vec2' must have 1 dimension, but has {vec2.ndim}")

    # what the inputs already say, rather than asking for it again
    n_elements = vec1.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if vec2.shape[0] != n_elements:
        raise ValueError(f"'vec2' has {vec2.shape[0]} along axis 0, but "
            f"'vec1' implies n_elements == {n_elements}"
        )

    # outputs and work arrays, which the caller never sees
    result = ctypes.c_double(0)
    ierr = ctypes.c_int(0)

    _lib.euclidean_distance_c(
        vec1,
        vec2,
        ctypes.byref(ctypes.c_int(n_elements)),
        ctypes.byref(result),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _EUCLIDEAN_DISTANCE_ARGUMENTS)

    return result.value

def distance_to_centroid(
        genes,
        centroids,
        gene_to_fam,
):
    r"""Compute distance from each gene to its corresponding family centroid.

    Parameters
    ----------
    genes : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene expression matrix (n_tissues × n_genes), column-major
    centroids : np.ndarray[np.float64] of shape (n_tissues, n_families,), column-major (order='F')
        Family centroid matrix (n_tissues × n_families), column-major, `-1.0_real64` for unassigned genes
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `genes`, using `0_int32` for unassigned genes

    Returns
    -------
    distances : np.ndarray[np.float64] of shape (n_genes,)
        Output distances array

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_euclidean_distance::distance_to_centroid`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        genes = np.asfortranarray(genes, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'genes' must be an array of np.float64: {error}") from None
    if genes.ndim != 2:
        raise ValueError(f"'genes' must have 2 dimensions, but has {genes.ndim}")
    try:
        centroids = np.asfortranarray(centroids, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'centroids' must be an array of np.float64: {error}") from None
    if centroids.ndim != 2:
        raise ValueError(f"'centroids' must have 2 dimensions, but has {centroids.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes = genes.shape[1]
    n_families = centroids.shape[1]
    n_tissues = genes.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'genes' implies n_genes == {n_genes}"
        )
    if centroids.shape[0] != n_tissues:
        raise ValueError(f"'centroids' has {centroids.shape[0]} along axis 0, but "
            f"'genes' implies n_tissues == {n_tissues}"
        )

    # outputs and work arrays, which the caller never sees
    distances = np.empty((n_genes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.distance_to_centroid_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        genes,
        centroids,
        gene_to_fam,
        distances,
        ctypes.byref(ctypes.c_int(n_tissues)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _DISTANCE_TO_CENTROID_ARGUMENTS)

    return distances
