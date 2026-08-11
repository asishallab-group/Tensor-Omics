"""tox_stc_csv

Plain-text (CSV/TSV) companions to `tox_stc_json`'s JSON/HTML report, for the data-science

Python binding, generated from tox_stc_csv. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.serialize_stc_points_as_csv_c.restype = None
_lib.serialize_stc_points_as_csv_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_STC_POINTS_AS_CSV_ARGUMENTS = ("filename", "n_vectors", "n_selected_seed", "max_group_size", "n_super_ensembles", "seed_selection_mask", "ensemble_masks", "ensemble_low_confidence_masks", "super_ensembles", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_STC_POINTS_AS_CSV_ARGUMENT_SOURCES = (None, "seed_selection_mask", "ensemble_masks", "super_ensembles", None, None, None, None, None, None,)

_lib.serialize_stc_ensemble_overlap_as_csv_c.restype = None
_lib.serialize_stc_ensemble_overlap_as_csv_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_STC_ENSEMBLE_OVERLAP_AS_CSV_ARGUMENTS = ("filename", "n_vectors", "n_selected_seed", "ensemble_masks", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_STC_ENSEMBLE_OVERLAP_AS_CSV_ARGUMENT_SOURCES = (None, "ensemble_masks", "ensemble_masks", None, None,)

_lib.serialize_stc_super_ensembles_as_tsv_c.restype = None
_lib.serialize_stc_super_ensembles_as_tsv_c.argtypes = (
    np.ctypeslib.ndpointer(ndim=1),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='F_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_SERIALIZE_STC_SUPER_ENSEMBLES_AS_TSV_ARGUMENTS = ("filename", "max_group_size", "n_super_ensembles", "super_ensembles", "ierr",)
#: For a derived argument, the one the caller passed it in
_SERIALIZE_STC_SUPER_ENSEMBLES_AS_TSV_ARGUMENT_SOURCES = (None, "super_ensembles", "super_ensembles", None, None,)

def serialize_stc_points_as_csv(
        filename,
        n_super_ensembles,
        seed_selection_mask,
        ensemble_masks,
        ensemble_low_confidence_masks,
        super_ensembles,
):
    r"""Serializes each input vector's ensemble/super-ensemble membership as CSV

    Parameters
    ----------
    filename : str
        Name of the CSV file to write
    n_super_ensembles : int
        Number of leading columns of `super_ensembles` actually filled
    seed_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Seed selection, see `seeds`
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble accepted membership, one column per seed
    ensemble_low_confidence_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble iteration-1 fallback membership
    super_ensembles : np.ndarray[np.int32] of shape (max_group_size, n_selected_seed*(n_selected_seed-1),), column-major (order='F')
        One super-ensemble per column, 0-padded, see `ensemble_reconciliation`

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc_csv::serialize_stc_points_as_csv`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        seed_selection_mask = np.ascontiguousarray(seed_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'seed_selection_mask' must be an array of np.bool_: {error}") from None
    if seed_selection_mask.ndim != 1:
        raise ValueError(f"'seed_selection_mask' must have 1 dimension, but has {seed_selection_mask.ndim}")
    try:
        ensemble_masks = np.asfortranarray(ensemble_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_masks' must be an array of np.bool_: {error}") from None
    if ensemble_masks.ndim != 2:
        raise ValueError(f"'ensemble_masks' must have 2 dimensions, but has {ensemble_masks.ndim}")
    try:
        ensemble_low_confidence_masks = np.asfortranarray(ensemble_low_confidence_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_low_confidence_masks' must be an array of np.bool_: {error}") from None
    if ensemble_low_confidence_masks.ndim != 2:
        raise ValueError(f"'ensemble_low_confidence_masks' must have 2 dimensions, but has {ensemble_low_confidence_masks.ndim}")
    try:
        super_ensembles = np.asfortranarray(super_ensembles, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'super_ensembles' must be an array of np.int32: {error}") from None
    if super_ensembles.ndim != 2:
        raise ValueError(f"'super_ensembles' must have 2 dimensions, but has {super_ensembles.ndim}")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    n_vectors = seed_selection_mask.shape[0]
    n_selected_seed = int(seed_selection_mask.sum())
    max_group_size = super_ensembles.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if ensemble_masks.shape[0] != n_vectors:
        raise ValueError(f"'ensemble_masks' has {ensemble_masks.shape[0]} along axis 0, but "
            f"'seed_selection_mask' implies n_vectors == {n_vectors}"
        )
    if ensemble_low_confidence_masks.shape[0] != n_vectors:
        raise ValueError(f"'ensemble_low_confidence_masks' has {ensemble_low_confidence_masks.shape[0]} along axis 0, but "
            f"'seed_selection_mask' implies n_vectors == {n_vectors}"
        )
    if ensemble_low_confidence_masks.shape[1] != n_selected_seed:
        raise ValueError(f"'ensemble_low_confidence_masks' has {ensemble_low_confidence_masks.shape[1]} along axis 1, but "
            f"'ensemble_masks' implies n_selected_seed == {n_selected_seed}"
        )

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_stc_points_as_csv_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_selected_seed)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        ctypes.byref(ctypes.c_int(n_super_ensembles)),
        seed_selection_mask,
        ensemble_masks,
        ensemble_low_confidence_masks,
        super_ensembles,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_STC_POINTS_AS_CSV_ARGUMENTS, _SERIALIZE_STC_POINTS_AS_CSV_ARGUMENT_SOURCES)

    return None

def serialize_stc_ensemble_overlap_as_csv(
        filename,
        ensemble_masks,
):
    r"""Serializes the full pairwise ensemble Overlap Coefficient matrix as CSV

    Parameters
    ----------
    filename : str
        Name of the CSV file to write
    ensemble_masks : np.ndarray[np.bool_] of shape (n_vectors, n_selected_seed,), column-major (order='F')
        Per-ensemble accepted membership, one column per seed

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc_csv::serialize_stc_ensemble_overlap_as_csv`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        ensemble_masks = np.asfortranarray(ensemble_masks, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'ensemble_masks' must be an array of np.bool_: {error}") from None
    if ensemble_masks.ndim != 2:
        raise ValueError(f"'ensemble_masks' must have 2 dimensions, but has {ensemble_masks.ndim}")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    n_vectors = ensemble_masks.shape[0]
    n_selected_seed = ensemble_masks.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_stc_ensemble_overlap_as_csv_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        ctypes.byref(ctypes.c_int(n_selected_seed)),
        ensemble_masks,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_STC_ENSEMBLE_OVERLAP_AS_CSV_ARGUMENTS, _SERIALIZE_STC_ENSEMBLE_OVERLAP_AS_CSV_ARGUMENT_SOURCES)

    return None

def serialize_stc_super_ensembles_as_tsv(
        filename,
        super_ensembles,
):
    r"""Serializes the super-ensembles as a gene-family-file-style TSV

    Parameters
    ----------
    filename : str
        Name of the TSV file to write
    super_ensembles : np.ndarray[np.int32] of shape (max_group_size, n_super_ensembles,), column-major (order='F')
        One super-ensemble per column, 0-padded, see `ensemble_reconciliation`

    Returns
    -------
    None

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_stc_csv::serialize_stc_super_ensembles_as_tsv`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    filename = np.array([str(filename).encode()], dtype="S")
    try:
        super_ensembles = np.asfortranarray(super_ensembles, dtype=np.int32)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'super_ensembles' must be an array of np.int32: {error}") from None
    if super_ensembles.ndim != 2:
        raise ValueError(f"'super_ensembles' must have 2 dimensions, but has {super_ensembles.ndim}")

    # what the inputs already say, rather than asking for it again
    filename_strlen = filename.itemsize
    max_group_size = super_ensembles.shape[0]
    n_super_ensembles = super_ensembles.shape[1]

    # outputs and work arrays, which the caller never sees
    ierr = ctypes.c_int(0)

    _lib.serialize_stc_super_ensembles_as_tsv_c(
        filename,
        ctypes.byref(ctypes.c_int(filename_strlen)),
        ctypes.byref(ctypes.c_int(max_group_size)),
        ctypes.byref(ctypes.c_int(n_super_ensembles)),
        super_ensembles,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _SERIALIZE_STC_SUPER_ENSEMBLES_AS_TSV_ARGUMENTS, _SERIALIZE_STC_SUPER_ENSEMBLES_AS_TSV_ARGUMENT_SOURCES)

    return None
