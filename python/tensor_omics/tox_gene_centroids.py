"""Python binding to Generated from the kernel; do not edit -- regenerate instead.

Generated from tox_gene_centroids. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.mean_vector_c.restype = None
_lib.mean_vector_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MEAN_VECTOR_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_indices", "n_selected_genes", "centroid", "ierr",)

_lib.group_centroid_orthologs_expert_c.restype = None
_lib.group_centroid_orthologs_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ORTHOLOGS_EXPERT_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "tmp_group_indices", "ortholog_set", "ierr",)

_lib.group_centroid_orthologs_c.restype = None
_lib.group_centroid_orthologs_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ORTHOLOGS_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "ortholog_set", "ierr",)

_lib.group_centroid_all_expert_c.restype = None
_lib.group_centroid_all_expert_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ALL_EXPERT_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "tmp_group_indices", "ierr",)

_lib.group_centroid_all_c.restype = None
_lib.group_centroid_all_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ALL_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "ierr",)

def mean_vector(
        expression_vectors,
        gene_indices,
):
    r"""Computes the element-wise mean for a given set of vectors.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_indices : np.ndarray[np.int32] of shape (n_selected_genes,)
        An array containing the column indices of the selected genes in 'expression_vectors'.
        The minimum valid value is `1`.
        The maximum valid value is `n_genes`.

    Returns
    -------
    centroid : np.ndarray[np.float64] of shape (n_axes,)
        The output vector representing the computed centroid.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_gene_centroids::mean_vector`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_indices = np.ascontiguousarray(gene_indices, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_indices' must be an array of np.int32: {error}") from None
    if gene_indices.ndim != 1:
        raise ValueError(f"'gene_indices' must have 1 dimension, but has {gene_indices.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_selected_genes = gene_indices.shape[0]

    # outputs and work arrays, which the caller never sees
    centroid = np.empty((n_axes,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.mean_vector_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_indices,
        ctypes.byref(ctypes.c_int(n_selected_genes)),
        centroid,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _MEAN_VECTOR_ARGUMENTS)

    return centroid

def group_centroid_orthologs_expert(
        expression_vectors,
        gene_to_family,
        n_families,
        ortholog_set,
):
    r"""Iterates over families, filters gene indices, and computes centroids.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    n_families : int
        Total number of gene families to compute centroids for.
    ortholog_set : np.ndarray[np.bool_] of shape (n_genes,)
        A logical array indicating if a gene is part of a specific subset (e.g., orthologs).

    Returns
    -------
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        The output matrix (n_axes x n_families) to store the computed centroids.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_orthologs`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_family' must be an array of np.int32: {error}") from None
    if gene_to_family.ndim != 1:
        raise ValueError(f"'gene_to_family' must have 1 dimension, but has {gene_to_family.ndim}")
    try:
        ortholog_set = np.ascontiguousarray(ortholog_set, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ortholog_set' must be an array of np.bool_: {error}") from None
    if ortholog_set.ndim != 1:
        raise ValueError(f"'ortholog_set' must have 1 dimension, but has {ortholog_set.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_family.shape[0] != n_genes:
        raise ValueError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if ortholog_set.shape[0] != n_genes:
        raise ValueError(f"'ortholog_set' has {ortholog_set.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    centroid_matrix = np.empty((n_axes, n_families,), dtype=np.float64, order='F')
    tmp_group_indices = np.empty((n_genes,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.group_centroid_orthologs_expert_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_families)),
        centroid_matrix,
        tmp_group_indices,
        ortholog_set,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GROUP_CENTROID_ORTHOLOGS_EXPERT_ARGUMENTS)

    return centroid_matrix

def group_centroid_orthologs(
        expression_vectors,
        gene_to_family,
        n_families,
        ortholog_set,
):
    r"""Iterates over families, filters gene indices, and computes centroids.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    n_families : int
        Total number of gene families to compute centroids for.
    ortholog_set : np.ndarray[np.bool_] of shape (n_genes,)
        A logical array indicating if a gene is part of a specific subset (e.g., orthologs).

    Returns
    -------
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        The output matrix (n_axes x n_families) to store the computed centroids.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_orthologs_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_family' must be an array of np.int32: {error}") from None
    if gene_to_family.ndim != 1:
        raise ValueError(f"'gene_to_family' must have 1 dimension, but has {gene_to_family.ndim}")
    try:
        ortholog_set = np.ascontiguousarray(ortholog_set, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ortholog_set' must be an array of np.bool_: {error}") from None
    if ortholog_set.ndim != 1:
        raise ValueError(f"'ortholog_set' must have 1 dimension, but has {ortholog_set.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_family.shape[0] != n_genes:
        raise ValueError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )
    if ortholog_set.shape[0] != n_genes:
        raise ValueError(f"'ortholog_set' has {ortholog_set.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    centroid_matrix = np.empty((n_axes, n_families,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.group_centroid_orthologs_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_families)),
        centroid_matrix,
        ortholog_set,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GROUP_CENTROID_ORTHOLOGS_ARGUMENTS)

    return centroid_matrix

def group_centroid_all_expert(
        expression_vectors,
        gene_to_family,
        n_families,
):
    r"""Iterates over families, filters gene indices, and computes centroids.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    n_families : int
        Total number of gene families to compute centroids for.

    Returns
    -------
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        The output matrix (n_axes x n_families) to store the computed centroids.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_all`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_family' must be an array of np.int32: {error}") from None
    if gene_to_family.ndim != 1:
        raise ValueError(f"'gene_to_family' must have 1 dimension, but has {gene_to_family.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_family.shape[0] != n_genes:
        raise ValueError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    centroid_matrix = np.empty((n_axes, n_families,), dtype=np.float64, order='F')
    tmp_group_indices = np.empty((n_genes,), dtype=np.int32, order='C')
    ierr = ctypes.c_int(0)

    _lib.group_centroid_all_expert_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_families)),
        centroid_matrix,
        tmp_group_indices,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GROUP_CENTROID_ALL_EXPERT_ARGUMENTS)

    return centroid_matrix

def group_centroid_all(
        expression_vectors,
        gene_to_family,
        n_families,
):
    r"""Iterates over families, filters gene indices, and computes centroids.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_genes,), column-major (order='F')
        The input matrix of all gene expression vectors (n_axes x n_genes).
    gene_to_family : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.
    n_families : int
        Total number of gene families to compute centroids for.

    Returns
    -------
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F')
        The output matrix (n_axes x n_families) to store the computed centroids.

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_all_alloc`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        gene_to_family = np.ascontiguousarray(gene_to_family, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_family' must be an array of np.int32: {error}") from None
    if gene_to_family.ndim != 1:
        raise ValueError(f"'gene_to_family' must have 1 dimension, but has {gene_to_family.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if gene_to_family.shape[0] != n_genes:
        raise ValueError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    centroid_matrix = np.empty((n_axes, n_families,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.group_centroid_all_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_genes)),
        gene_to_family,
        ctypes.byref(ctypes.c_int(n_families)),
        centroid_matrix,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _GROUP_CENTROID_ALL_ARGUMENTS)

    return centroid_matrix
