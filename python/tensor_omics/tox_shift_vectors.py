"""tox_shift_vectors

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_shift_vectors. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.compute_shift_vector_field_c.restype = None
_lib.compute_shift_vector_field_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=3, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_SHIFT_VECTOR_FIELD_ARGUMENTS = ("n_tissues", "n_genes", "n_families", "expression_vectors", "family_centroids", "gene_to_fam", "shift_vectors", "ierr",)

def compute_shift_vector_field(
        expression_vectors,
        family_centroids,
        gene_to_fam,
):
    r"""Compute the shift vector field for all genes.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_tissues, n_genes,), column-major (order='F')
        Gene expression matrix
    family_centroids : np.ndarray[np.float64] of shape (n_tissues, n_families,), column-major (order='F')
        Family centroid matrix
    gene_to_fam : np.ndarray[np.int32] of shape (n_genes,)
        Index mapping -> each index `i` holds the family index for the corresponding gene in `expression_vectors`, using `0` for unassigned genes
        The minimum valid value is `1`.
        The maximum valid value is `n_families`.
        The value `0` is additionally accepted.

    Returns
    -------
    shift_vectors : np.ndarray[np.float64] of shape (n_tissues, 2, n_genes,), column-major (order='F')
        Output, real matrix array, stores the centroid of the gene's family in `shift_vectors(:, 1, i_gene)` (zero vector if no family assigned) and the shift vectors in `shift_vectors(:, 2, i_gene)`

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_shift_vectors::compute_shift_vector_field`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'family_centroids' must be an array of np.float64: {error}") from None
    if family_centroids.ndim != 2:
        raise ValueError(f"'family_centroids' must have 2 dimensions, but has {family_centroids.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_tissues = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_families = family_centroids.shape[1]

    # Fortran cannot check that shared extents agree; this can
    if family_centroids.shape[0] != n_tissues:
        raise ValueError(f"'family_centroids' has {family_centroids.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_tissues == {n_tissues}"
        )
    if gene_to_fam.shape[0] != n_genes:
        raise ValueError(f"'gene_to_fam' has {gene_to_fam.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_genes == {n_genes}"
        )

    # outputs and work arrays, which the caller never sees
    shift_vectors = np.empty((n_tissues, 2, n_genes,), dtype=np.float64, order='F')
    ierr = ctypes.c_int(0)

    _lib.compute_shift_vector_field_c(
        ctypes.byref(ctypes.c_int(n_tissues)),
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        expression_vectors,
        family_centroids,
        gene_to_fam,
        shift_vectors,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_SHIFT_VECTOR_FIELD_ARGUMENTS)

    return shift_vectors
