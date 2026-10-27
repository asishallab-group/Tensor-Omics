r"""tox_gene_centroids

Expression centroids of gene families.

`mean_vector` is the centroid of a set of expression vectors. `group_centroid_orthologs`
and `group_centroid_all` take the centroid of a family: over its orthologs only, or over
every gene in it -- two routines rather than one taking a flag, so which set a result is
over is visible at the call site.

Python binding, generated from tox_gene_centroids. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.mean_vector_c.restype = None
_lib.mean_vector_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_MEAN_VECTOR_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_indices", "n_selected_genes", "centroid", "ierr",)
#: For a derived argument, the one the caller passed it in
_MEAN_VECTOR_ARGUMENT_SOURCES = (None, "expression_vectors", "expression_vectors", None, "gene_indices", None, None,)

_lib.group_centroid_orthologs_c.restype = None
_lib.group_centroid_orthologs_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ORTHOLOGS_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "ortholog_set", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROUP_CENTROID_ORTHOLOGS_ARGUMENT_SOURCES = (None, "expression_vectors", "expression_vectors", None, "centroid_matrix", None, None, None,)

_lib.group_centroid_all_c.restype = None
_lib.group_centroid_all_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GROUP_CENTROID_ALL_ARGUMENTS = ("expression_vectors", "n_axes", "n_genes", "gene_to_family", "n_families", "centroid_matrix", "ierr",)
#: For a derived argument, the one the caller passed it in
_GROUP_CENTROID_ALL_ARGUMENT_SOURCES = (None, "expression_vectors", "expression_vectors", None, "centroid_matrix", None, None,)

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
    centroid : np.ndarray[np.float64] of shape (n_axes,), read-only
        The output vector representing the computed centroid.
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
    Generated from the Fortran procedure `tox_gene_centroids::mean_vector`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    gene_indices = _contract.integers(gene_indices, "gene_indices", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]
    n_selected_genes = gene_indices.shape[0]

    # outputs and work arrays, which the caller never sees
    centroid = _np.empty((n_axes,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.mean_vector_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        gene_indices,
        _ctypes.byref(_ctypes.c_int(n_selected_genes)),
        centroid,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _MEAN_VECTOR_ARGUMENTS, _MEAN_VECTOR_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    centroid.flags.writeable = False

    return centroid

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
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F'), read-only
        The output matrix (n_axes x n_families) to store the computed centroids.
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
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_orthologs`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    gene_to_family = _contract.integers(gene_to_family, "gene_to_family", ndim=1)
    n_families = _contract.integer(n_families, "n_families")
    ortholog_set = _contract.logicals(ortholog_set, "ortholog_set", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_family.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_to_family' is declared with shape (n_genes)", argument="gene_to_family")
    if ortholog_set.shape[0] != n_genes:
        raise _ToxArgumentError(f"'ortholog_set' has {ortholog_set.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'ortholog_set' is declared with shape (n_genes)", argument="ortholog_set")

    # outputs and work arrays, which the caller never sees
    centroid_matrix = _np.empty((n_axes, n_families,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.group_centroid_orthologs_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        gene_to_family,
        _ctypes.byref(_ctypes.c_int(n_families)),
        centroid_matrix,
        ortholog_set,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GROUP_CENTROID_ORTHOLOGS_ARGUMENTS, _GROUP_CENTROID_ORTHOLOGS_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    centroid_matrix.flags.writeable = False

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
    centroid_matrix : np.ndarray[np.float64] of shape (n_axes, n_families,), column-major (order='F'), read-only
        The output matrix (n_axes x n_families) to store the computed centroids.
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
    Generated from the Fortran procedure `tox_gene_centroids::group_centroid_all`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    expression_vectors = _contract.reals(expression_vectors, "expression_vectors", ndim=2, order="F")
    gene_to_family = _contract.integers(gene_to_family, "gene_to_family", ndim=1)
    n_families = _contract.integer(n_families, "n_families")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_genes = expression_vectors.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if gene_to_family.shape[0] != n_genes:
        raise _ToxArgumentError(f"'gene_to_family' has {gene_to_family.shape[0]} along axis 0, but 'expression_vectors' implies n_genes == {n_genes}; 'gene_to_family' is declared with shape (n_genes)", argument="gene_to_family")

    # outputs and work arrays, which the caller never sees
    centroid_matrix = _np.empty((n_axes, n_families,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.group_centroid_all_c(
        expression_vectors,
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_genes)),
        gene_to_family,
        _ctypes.byref(_ctypes.c_int(n_families)),
        centroid_matrix,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GROUP_CENTROID_ALL_ARGUMENTS, _GROUP_CENTROID_ALL_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    centroid_matrix.flags.writeable = False

    return centroid_matrix
