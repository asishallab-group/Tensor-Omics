"""Python interface to Semantic validation of TensorOmics data sets (dimensions, ID uniqueness, value ranges,

Generated from tox_data_validation. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.validate_data_structure_c.restype = None
_lib.validate_data_structure_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_DATA_STRUCTURE_ARGUMENTS = ("n_genes", "n_families", "n_samples", "gene_ids", "gene_family_ids", "gene_to_fam", "expression_vectors", "family_centroids", "shift_vectors", "ierr",)

_lib.validate_gene_to_family_mapping_c.restype = None
_lib.validate_gene_to_family_mapping_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_GENE_TO_FAMILY_MAPPING_ARGUMENTS = ("gene_to_fam", "n_families", "ierr",)

_lib.validate_expression_data_c.restype = None
_lib.validate_expression_data_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_EXPRESSION_DATA_ARGUMENTS = ("expression_vectors", "check_non_negative", "ierr",)

_lib.validate_family_centroids_c.restype = None
_lib.validate_family_centroids_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_FAMILY_CENTROIDS_ARGUMENTS = ("family_centroids", "ierr",)

_lib.validate_shift_vectors_c.restype = None
_lib.validate_shift_vectors_c.argtypes = (
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_SHIFT_VECTORS_ARGUMENTS = ("shift_vectors", "expression_vectors", "family_centroids", "gene_to_fam", "n_samples", "ierr",)

_lib.validate_string_array_uniqueness_c.restype = None
_lib.validate_string_array_uniqueness_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_STRING_ARRAY_UNIQUENESS_ARGUMENTS = ("str_arr", "ierr",)

_lib.validate_all_data_c.restype = None
_lib.validate_all_data_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_bool),
    ctypes.POINTER(ctypes.c_bool),
)

#: The wrapped procedure's arguments, so an error can name one
_VALIDATE_ALL_DATA_ARGUMENTS = ("n_genes", "n_families", "n_samples", "gene_ids", "gene_family_ids", "gene_to_fam", "expression_vectors", "family_centroids", "shift_vectors", "ierr", "check_uniqueness", "check_shift_consistency",)

def validate_data_structure(
        n_genes,
        n_families,
        n_samples,
        gene_ids,
        gene_family_ids,
        gene_to_fam,
        expression_vectors,
        family_centroids,
        shift_vectors,
):
    r"""Validate full data structure

    Parameters
    ----------
    n_genes : int
        Expected number of genes
    n_families : int
        Expected number of families
    n_samples : int
        Expected number of samples
    gene_ids : sequence of str, of length n_gene_ids_elements
        Gene ids
    gene_family_ids : sequence of str, of length n_gene_family_ids_elements
        Gene family ids
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        Expression vectors
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        Family centroids
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F')
        Shift vectors

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_data_structure`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_ids = np.asarray([str(_s).encode() for _s in gene_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_ids' must be a sequence of strings: {error}") from None
    if gene_ids.ndim != 1:
        raise ValueError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}")
    try:
        gene_family_ids = np.asarray([str(_s).encode() for _s in gene_family_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_family_ids' must be a sequence of strings: {error}") from None
    if gene_family_ids.ndim != 1:
        raise ValueError(f"'gene_family_ids' must have 1 dimension, but has {gene_family_ids.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
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
        shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'shift_vectors' must be an array of np.float64: {error}") from None
    if shift_vectors.ndim != 2:
        raise ValueError(f"'shift_vectors' must have 2 dimensions, but has {shift_vectors.ndim}")

    # what the inputs already say, rather than asking for it again
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]
    gene_family_ids_strlen = gene_family_ids.itemsize
    n_gene_family_ids_elements = gene_family_ids.shape[0]
    n_gene_to_fam_elements = gene_to_fam.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_data_structure_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ctypes.c_int(n_samples)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        gene_family_ids,
        ctypes.byref(ctypes.c_int(gene_family_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_family_ids_elements)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_gene_to_fam_elements)),
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_2)),
        family_centroids,
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_2)),
        shift_vectors,
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_2)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VALIDATE_DATA_STRUCTURE_ARGUMENTS)

    return None

def validate_gene_to_family_mapping(
        gene_to_fam,
        n_families,
):
    r"""Validate gene to family mapping

    Parameters
    ----------
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    n_families : int
        number of families

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_gene_to_family_mapping`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")

    # what the inputs already say, rather than asking for it again
    n_gene_to_fam_elements = gene_to_fam.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_gene_to_family_mapping_c(
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_gene_to_fam_elements)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VALIDATE_GENE_TO_FAMILY_MAPPING_ARGUMENTS)

    return None

def validate_expression_data(
        expression_vectors,
        check_non_negative,
):
    r"""Validate expresssion data

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        Expression vectors
    check_non_negative : bool
        Defines if non negative should be checked

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_expression_data`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")

    # what the inputs already say, rather than asking for it again
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_expression_data_c(
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_2)),
        ctypes.byref(ctypes.c_bool(check_non_negative)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VALIDATE_EXPRESSION_DATA_ARGUMENTS)

    return None

def validate_family_centroids(
        family_centroids,
):
    r"""Validate the family centroids

    Parameters
    ----------
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        Family centroids array

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_family_centroids`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        family_centroids = np.asfortranarray(family_centroids, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'family_centroids' must be an array of np.float64: {error}") from None
    if family_centroids.ndim != 2:
        raise ValueError(f"'family_centroids' must have 2 dimensions, but has {family_centroids.ndim}")

    # what the inputs already say, rather than asking for it again
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_family_centroids_c(
        family_centroids,
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_2)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VALIDATE_FAMILY_CENTROIDS_ARGUMENTS)

    return None

def validate_shift_vectors(
        shift_vectors,
        expression_vectors,
        family_centroids,
        gene_to_fam,
        n_samples,
):
    r"""Validates shift vectors

    Parameters
    ----------
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F')
        shift vectors
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        expression vectors
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        family centroids
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    n_samples : int
        Number of samples

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_shift_vectors`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'shift_vectors' must be an array of np.float64: {error}") from None
    if shift_vectors.ndim != 2:
        raise ValueError(f"'shift_vectors' must have 2 dimensions, but has {shift_vectors.ndim}")
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
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    n_gene_to_fam_elements = gene_to_fam.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_shift_vectors_c(
        shift_vectors,
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_2)),
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_2)),
        family_centroids,
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_2)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_gene_to_fam_elements)),
        ctypes.byref(ctypes.c_int(n_samples)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VALIDATE_SHIFT_VECTORS_ARGUMENTS)

    return None

def validate_string_array_uniqueness(
        str_arr,
):
    r"""Validate that no string appears more than once

    Parameters
    ----------
    str_arr : sequence of str, of length n_str_arr_elements
        string array

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_string_array_uniqueness`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        str_arr = np.asarray([str(_s).encode() for _s in str_arr], dtype="S")
    except TypeError as error:
        raise TypeError(f"'str_arr' must be a sequence of strings: {error}") from None
    if str_arr.ndim != 1:
        raise ValueError(f"'str_arr' must have 1 dimension, but has {str_arr.ndim}")

    # what the inputs already say, rather than asking for it again
    str_arr_strlen = str_arr.itemsize
    n_str_arr_elements = str_arr.shape[0]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_string_array_uniqueness_c(
        str_arr,
        ctypes.byref(ctypes.c_int(str_arr_strlen)),
        ctypes.byref(ctypes.c_int(n_str_arr_elements)),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _VALIDATE_STRING_ARRAY_UNIQUENESS_ARGUMENTS)

    return None

def validate_all_data(
        n_genes,
        n_families,
        n_samples,
        gene_ids,
        gene_family_ids,
        gene_to_fam,
        expression_vectors,
        family_centroids,
        shift_vectors,
        check_uniqueness=True,
        check_shift_consistency=True,
):
    r"""Comprehensive validation routine, combining all checks

    Parameters
    ----------
    n_genes : int
        Number of genes
    n_families : int
        Number of families
    n_samples : int
        Number of samples
    gene_ids : sequence of str, of length n_gene_ids_elements
        Gene ids array
    gene_family_ids : sequence of str, of length n_gene_family_ids_elements
        gene family ids
    gene_to_fam : np.ndarray[np.int32] of shape (n_gene_to_fam_elements,)
        gene to family mapping
    expression_vectors : np.ndarray[np.float64] of shape (n_expression_vectors_elements_dim_1, n_expression_vectors_elements_dim_2,), column-major (order='F')
        Expression vectors
    family_centroids : np.ndarray[np.float64] of shape (n_family_centroids_elements_dim_1, n_family_centroids_elements_dim_2,), column-major (order='F')
        family centroids
    shift_vectors : np.ndarray[np.float64] of shape (n_shift_vectors_elements_dim_1, n_shift_vectors_elements_dim_2,), column-major (order='F')
        shift vectors
    check_uniqueness : bool, optional, default True
        Check ID arrays for uniqueness.
        The default value is `.true.`.
    check_shift_consistency : bool, optional, default True
        Check consitency of shift array.
        The default value is `.true.`.

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_validation::validate_all_data`.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_ids = np.asarray([str(_s).encode() for _s in gene_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_ids' must be a sequence of strings: {error}") from None
    if gene_ids.ndim != 1:
        raise ValueError(f"'gene_ids' must have 1 dimension, but has {gene_ids.ndim}")
    try:
        gene_family_ids = np.asarray([str(_s).encode() for _s in gene_family_ids], dtype="S")
    except TypeError as error:
        raise TypeError(f"'gene_family_ids' must be a sequence of strings: {error}") from None
    if gene_family_ids.ndim != 1:
        raise ValueError(f"'gene_family_ids' must have 1 dimension, but has {gene_family_ids.ndim}")
    try:
        gene_to_fam = np.ascontiguousarray(gene_to_fam, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'gene_to_fam' must be an array of np.int32: {error}") from None
    if gene_to_fam.ndim != 1:
        raise ValueError(f"'gene_to_fam' must have 1 dimension, but has {gene_to_fam.ndim}")
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
        shift_vectors = np.asfortranarray(shift_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'shift_vectors' must be an array of np.float64: {error}") from None
    if shift_vectors.ndim != 2:
        raise ValueError(f"'shift_vectors' must have 2 dimensions, but has {shift_vectors.ndim}")

    # what the inputs already say, rather than asking for it again
    gene_ids_strlen = gene_ids.itemsize
    n_gene_ids_elements = gene_ids.shape[0]
    gene_family_ids_strlen = gene_family_ids.itemsize
    n_gene_family_ids_elements = gene_family_ids.shape[0]
    n_gene_to_fam_elements = gene_to_fam.shape[0]
    n_expression_vectors_elements_dim_1 = expression_vectors.shape[0]
    n_expression_vectors_elements_dim_2 = expression_vectors.shape[1]
    n_family_centroids_elements_dim_1 = family_centroids.shape[0]
    n_family_centroids_elements_dim_2 = family_centroids.shape[1]
    n_shift_vectors_elements_dim_1 = shift_vectors.shape[0]
    n_shift_vectors_elements_dim_2 = shift_vectors.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.validate_all_data_c(
        ctypes.byref(ctypes.c_int(n_genes)),
        ctypes.byref(ctypes.c_int(n_families)),
        ctypes.byref(ctypes.c_int(n_samples)),
        gene_ids,
        ctypes.byref(ctypes.c_int(gene_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_ids_elements)),
        gene_family_ids,
        ctypes.byref(ctypes.c_int(gene_family_ids_strlen)),
        ctypes.byref(ctypes.c_int(n_gene_family_ids_elements)),
        gene_to_fam,
        ctypes.byref(ctypes.c_int(n_gene_to_fam_elements)),
        expression_vectors,
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_expression_vectors_elements_dim_2)),
        family_centroids,
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_family_centroids_elements_dim_2)),
        shift_vectors,
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_1)),
        ctypes.byref(ctypes.c_int(n_shift_vectors_elements_dim_2)),
        ctypes.byref(ierr),
        ctypes.byref(ctypes.c_bool(check_uniqueness)),
        ctypes.byref(ctypes.c_bool(check_shift_consistency)),
    )

    check_err_code(ierr.value, _VALIDATE_ALL_DATA_ARGUMENTS)

    return None
