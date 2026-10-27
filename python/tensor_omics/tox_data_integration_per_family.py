r"""tox_data_integration_per_family

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) JSD Calculation per family

The JSD value for a sub-neighborhood -- typically the genes of one family -- obtained by
driving the same pipeline over a masked set of neighbors. Answers whether two studies are
compatible *for this family*, which the global figure can hide either way.

Python binding, generated from tox_data_integration_per_family. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.fjct_compute_jsd_c.restype = None
_lib.fjct_compute_jsd_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FJCT_COMPUTE_JSD_ARGUMENTS = ("family_idx", "gene_to_family_S1", "gene_to_family_S2", "n_genes_S1", "n_genes_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_genes_S1", "neighborhood_genes_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "ierr",)
#: For a derived argument, the one the caller passed it in
_FJCT_COMPUTE_JSD_ARGUMENT_SOURCES = (None, None, None, "gene_to_family_S1", "gene_to_family_S2", None, None, None, None, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", None, None, None, None, None, None, None, None, None,)

_lib.fjct_compute_masked_jsd_c.restype = None
_lib.fjct_compute_masked_jsd_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FJCT_COMPUTE_MASKED_JSD_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "neighbor_mask_S1", "neighbor_mask_S2", "n_bins", "shared_residual_range", "js_divergences", "included_n_reps_S1", "included_n_reps_S2", "total_included_n_reps", "global_js_divergence", "weights", "pmf_S1", "pmf_S2", "ierr",)
#: For a derived argument, the one the caller passed it in
_FJCT_COMPUTE_MASKED_JSD_ARGUMENT_SOURCES = (None, None, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", None, None, "pmf_S1", None, None, None, None, None, None, None, None, None, None,)

_lib.fjct_compute_contribution_scores_c.restype = None
_lib.fjct_compute_contribution_scores_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_FJCT_COMPUTE_CONTRIBUTION_SCORES_ARGUMENTS = ("global_js_divergences", "total_included_n_reps_per_f", "k_families", "support_weights", "contribution_scores", "ierr",)
#: For a derived argument, the one the caller passed it in
_FJCT_COMPUTE_CONTRIBUTION_SCORES_ARGUMENT_SOURCES = (None, None, "global_js_divergences", None, None, None,)

def fjct_compute_jsd(
        family_idx,
        gene_to_family_S1,
        gene_to_family_S2,
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighborhood_genes_S1,
        neighborhood_genes_S2,
        n_bins,
        shared_residual_range,
):
    r"""Compute the family-level compatibility score between two studies for a single gene family

    Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts
    the residual samples to the genes belonging to the family `family_idx`.

    Parameters
    ----------
    family_idx : int
        Index of the family that should be analyzed
        The minimum valid value is `1`.
    gene_to_family_S1 : np.ndarray[np.int32] of shape (n_genes_S1,)
        Mapping for study 1: Each index (gene) holds the index of its family
        The minimum valid value is `0`.
    gene_to_family_S2 : np.ndarray[np.int32] of shape (n_genes_S2,)
        Mapping for study 2: Each index (gene) holds the index of its family
        The minimum valid value is `0`.
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_genes_S1 : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F')
        Indices of the selected neighborhood genes of study 1
        The minimum valid value is `1`.
        The maximum valid value is `n_genes_S1`.
    neighborhood_genes_S2 : np.ndarray[np.int32] of shape (n_neighbors, n_points,), column-major (order='F')
        Indices of the selected neighborhood genes of study 2
        The minimum valid value is `1`.
        The maximum valid value is `n_genes_S2`.
    n_bins : int
        Number of equally sized histogram bins used for the studies
    shared_residual_range : float
        Computed residual range for both studies
        The minimum valid value is `0.0`.

    Returns
    -------
    dict
        with keys:

        js_divergences : np.ndarray[np.float64] of shape (n_points,), read-only
            Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals (included ones) in study 1
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals (included ones) in study 2
            A result is a value; call `.copy()` to obtain a modifiable array.
        total_included_n_reps : int
            Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,), read-only
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
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
    Generated from the Fortran procedure `tox_data_integration_per_family::fjct_compute_jsd`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        gene_to_family_S1 = _np.ascontiguousarray(gene_to_family_S1, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_to_family_S1' must be an array of np.int32: {_error}", argument="gene_to_family_S1") from None
    if gene_to_family_S1.ndim != 1:
        raise _ToxArgumentError(f"'gene_to_family_S1' must have 1 dimension, but has {gene_to_family_S1.ndim}", argument="gene_to_family_S1")
    try:
        gene_to_family_S2 = _np.ascontiguousarray(gene_to_family_S2, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'gene_to_family_S2' must be an array of np.int32: {_error}", argument="gene_to_family_S2") from None
    if gene_to_family_S2.ndim != 1:
        raise _ToxArgumentError(f"'gene_to_family_S2' must have 1 dimension, but has {gene_to_family_S2.ndim}", argument="gene_to_family_S2")
    try:
        neighborhood_residuals_S1 = _np.asfortranarray(neighborhood_residuals_S1, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals_S1' must be an array of np.float64: {_error}", argument="neighborhood_residuals_S1") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}", argument="neighborhood_residuals_S1")
    try:
        neighborhood_residuals_S2 = _np.asfortranarray(neighborhood_residuals_S2, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' must be an array of np.float64: {_error}", argument="neighborhood_residuals_S2") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}", argument="neighborhood_residuals_S2")
    try:
        neighborhood_genes_S1 = _np.asfortranarray(neighborhood_genes_S1, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_genes_S1' must be an array of np.int32: {_error}", argument="neighborhood_genes_S1") from None
    if neighborhood_genes_S1.ndim != 2:
        raise _ToxArgumentError(f"'neighborhood_genes_S1' must have 2 dimensions, but has {neighborhood_genes_S1.ndim}", argument="neighborhood_genes_S1")
    try:
        neighborhood_genes_S2 = _np.asfortranarray(neighborhood_genes_S2, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_genes_S2' must be an array of np.int32: {_error}", argument="neighborhood_genes_S2") from None
    if neighborhood_genes_S2.ndim != 2:
        raise _ToxArgumentError(f"'neighborhood_genes_S2' must have 2 dimensions, but has {neighborhood_genes_S2.ndim}", argument="neighborhood_genes_S2")

    # what the inputs already say, rather than asking for it again
    n_genes_S1 = gene_to_family_S1.shape[0]
    n_genes_S2 = gene_to_family_S2.shape[0]
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}", argument="neighborhood_residuals_S2")
    if neighborhood_genes_S1.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_genes_S1' has {neighborhood_genes_S1.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}", argument="neighborhood_genes_S1")
    if neighborhood_genes_S2.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_genes_S2' has {neighborhood_genes_S2.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}", argument="neighborhood_genes_S2")
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but 'neighborhood_residuals_S1' implies n_points == {n_points}", argument="neighborhood_residuals_S2")
    if neighborhood_genes_S1.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighborhood_genes_S1' has {neighborhood_genes_S1.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}", argument="neighborhood_genes_S1")
    if neighborhood_genes_S2.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighborhood_genes_S2' has {neighborhood_genes_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}", argument="neighborhood_genes_S2")

    # outputs and work arrays, which the caller never sees
    js_divergences = _np.empty((n_points,), dtype=_np.float64, order='C')
    included_n_reps_S1 = _np.empty((n_points,), dtype=_np.int32, order='C')
    included_n_reps_S2 = _np.empty((n_points,), dtype=_np.int32, order='C')
    total_included_n_reps = _ctypes.c_int(0)
    global_js_divergence = _ctypes.c_double(0)
    weights = _np.empty((n_points,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.fjct_compute_jsd_c(
        _ctypes.byref(_ctypes.c_int(family_idx)),
        gene_to_family_S1,
        gene_to_family_S2,
        _ctypes.byref(_ctypes.c_int(n_genes_S1)),
        _ctypes.byref(_ctypes.c_int(n_genes_S2)),
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighborhood_genes_S1,
        neighborhood_genes_S2,
        _ctypes.byref(_ctypes.c_int(n_reps_S1)),
        _ctypes.byref(_ctypes.c_int(n_reps_S2)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_double(shared_residual_range)),
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
        _ctypes.byref(total_included_n_reps),
        _ctypes.byref(global_js_divergence),
        weights,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _FJCT_COMPUTE_JSD_ARGUMENTS, _FJCT_COMPUTE_JSD_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    js_divergences.flags.writeable = False
    included_n_reps_S1.flags.writeable = False
    included_n_reps_S2.flags.writeable = False
    weights.flags.writeable = False

    return {
        "js_divergences": js_divergences,
        "included_n_reps_S1": included_n_reps_S1,
        "included_n_reps_S2": included_n_reps_S2,
        "total_included_n_reps": total_included_n_reps.value,
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
    }

def fjct_compute_masked_jsd(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        neighbor_mask_S1,
        neighbor_mask_S2,
        n_bins,
        shared_residual_range,
):
    r"""Compute the compatibility score between two studies for a single masked sub-neighborhood

    Reuses the same conditioning-on-mean-expression pipeline as the global gJCT, but restricts the
    residual samples to the neighbors selected by `neighbor_mask_S1`/`neighbor_mask_S2`. Typically
    those are all neighbors belonging to one gene family, which is what `fjct_compute_jsd` builds
    the masks for from a family index.

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F')
        Mask selecting the neighbors of study 1 to include
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F')
        Mask selecting the neighbors of study 2 to include
    n_bins : int
        Number of equally sized histogram bins used for the studies
    shared_residual_range : float
        Computed residual range for both studies
        The minimum valid value is `0.0`.

    Returns
    -------
    dict
        with keys:

        js_divergences : np.ndarray[np.float64] of shape (n_points,), read-only
            Jensen-Shannon divergence per reference point, computed for studies S1 and S2
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps_S1 : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals (included ones) in study 1
            A result is a value; call `.copy()` to obtain a modifiable array.
        included_n_reps_S2 : np.ndarray[np.int32] of shape (n_points,), read-only
            Count of non-NaN residuals (included ones) in study 2
            A result is a value; call `.copy()` to obtain a modifiable array.
        total_included_n_reps : int
            Total number of included replicates from both studies (\( \text{sum}(included\_n\_reps\_S1) + \text{sum}(included\_n\_reps\_S2) \))
        global_js_divergence : float
            Weighted global Jensen-Shannon divergence
        weights : np.ndarray[np.float64] of shape (n_points,), read-only
            Weights used for calculating the global weighted Jensen-Shannon divergence `global_js_divergence`
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmf_S1 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Normalized histogram counts for study 1
            A result is a value; call `.copy()` to obtain a modifiable array.
        pmf_S2 : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F'), read-only
            Normalized histogram counts for study 2
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
    Generated from the Fortran procedure `tox_data_integration_per_family::fjct_compute_masked_jsd`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        neighborhood_residuals_S1 = _np.asfortranarray(neighborhood_residuals_S1, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals_S1' must be an array of np.float64: {_error}", argument="neighborhood_residuals_S1") from None
    if neighborhood_residuals_S1.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals_S1' must have 3 dimensions, but has {neighborhood_residuals_S1.ndim}", argument="neighborhood_residuals_S1")
    try:
        neighborhood_residuals_S2 = _np.asfortranarray(neighborhood_residuals_S2, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' must be an array of np.float64: {_error}", argument="neighborhood_residuals_S2") from None
    if neighborhood_residuals_S2.ndim != 3:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' must have 3 dimensions, but has {neighborhood_residuals_S2.ndim}", argument="neighborhood_residuals_S2")
    try:
        neighbor_mask_S1 = _np.asfortranarray(neighbor_mask_S1, dtype=_np.bool_)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighbor_mask_S1' must be an array of np.bool_: {_error}", argument="neighbor_mask_S1") from None
    if neighbor_mask_S1.ndim != 2:
        raise _ToxArgumentError(f"'neighbor_mask_S1' must have 2 dimensions, but has {neighbor_mask_S1.ndim}", argument="neighbor_mask_S1")
    try:
        neighbor_mask_S2 = _np.asfortranarray(neighbor_mask_S2, dtype=_np.bool_)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'neighbor_mask_S2' must be an array of np.bool_: {_error}", argument="neighbor_mask_S2") from None
    if neighbor_mask_S2.ndim != 2:
        raise _ToxArgumentError(f"'neighbor_mask_S2' must have 2 dimensions, but has {neighbor_mask_S2.ndim}", argument="neighbor_mask_S2")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that shared extents agree; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}", argument="neighborhood_residuals_S2")
    if neighbor_mask_S1.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}", argument="neighbor_mask_S1")
    if neighbor_mask_S2.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}", argument="neighbor_mask_S2")
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but 'neighborhood_residuals_S1' implies n_points == {n_points}", argument="neighborhood_residuals_S2")
    if neighbor_mask_S1.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}", argument="neighbor_mask_S1")
    if neighbor_mask_S2.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}", argument="neighbor_mask_S2")

    # outputs and work arrays, which the caller never sees
    js_divergences = _np.empty((n_points,), dtype=_np.float64, order='C')
    included_n_reps_S1 = _np.empty((n_points,), dtype=_np.int32, order='C')
    included_n_reps_S2 = _np.empty((n_points,), dtype=_np.int32, order='C')
    total_included_n_reps = _ctypes.c_int(0)
    global_js_divergence = _ctypes.c_double(0)
    weights = _np.empty((n_points,), dtype=_np.float64, order='C')
    pmf_S1 = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    pmf_S2 = _np.empty((n_points, n_bins,), dtype=_np.float64, order='F')
    ierr = _ctypes.c_int(0)

    _lib.fjct_compute_masked_jsd_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        _ctypes.byref(_ctypes.c_int(n_reps_S1)),
        _ctypes.byref(_ctypes.c_int(n_reps_S2)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        neighbor_mask_S1,
        neighbor_mask_S2,
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_double(shared_residual_range)),
        js_divergences,
        included_n_reps_S1,
        included_n_reps_S2,
        _ctypes.byref(total_included_n_reps),
        _ctypes.byref(global_js_divergence),
        weights,
        pmf_S1,
        pmf_S2,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _FJCT_COMPUTE_MASKED_JSD_ARGUMENTS, _FJCT_COMPUTE_MASKED_JSD_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    js_divergences.flags.writeable = False
    included_n_reps_S1.flags.writeable = False
    included_n_reps_S2.flags.writeable = False
    weights.flags.writeable = False
    pmf_S1.flags.writeable = False
    pmf_S2.flags.writeable = False

    return {
        "js_divergences": js_divergences,
        "included_n_reps_S1": included_n_reps_S1,
        "included_n_reps_S2": included_n_reps_S2,
        "total_included_n_reps": total_included_n_reps.value,
        "global_js_divergence": global_js_divergence.value,
        "weights": weights,
        "pmf_S1": pmf_S1,
        "pmf_S2": pmf_S2,
    }

def fjct_compute_contribution_scores(
        global_js_divergences,
        total_included_n_reps_per_f,
):
    r"""Compute the per-sub-neighborhood contribution score

    Combines

    1. how divergent the family is between the studies (`global_js_divergences`), and
    2. how much residual support the family has overall (`total_included_n_reps_per_f`),

    using the outputs of `fjct_compute_jsd`, collected for the analyzed sub-neighborhoods.

    Parameters
    ----------
    global_js_divergences : np.ndarray[np.float64] of shape (k_families,)
        Per-sub-neighborhood weighted global JSD
        The minimum valid value is `0.0`.
    total_included_n_reps_per_f : np.ndarray[np.int32] of shape (k_families,)
        Per-sub-neighborhood `total_included_n_reps`
        The minimum valid value is `0`.

    Returns
    -------
    dict
        with keys:

        support_weights : np.ndarray[np.float64] of shape (k_families,), read-only
            Per-sub-neighborhood calculated support weight (ratio between its `total_included_n_reps` and `sum(total_included_n_reps_per_f)`, zero if there were no replicates included at all)
            A result is a value; call `.copy()` to obtain a modifiable array.
        contribution_scores : np.ndarray[np.float64] of shape (k_families,), read-only
            Per-sub-neighborhood calculated contribution ( \( support\_weights_i * global\_js\_divergences_i \) )
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
    Generated from the Fortran procedure `tox_data_integration_per_family::fjct_compute_contribution_scores`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        global_js_divergences = _np.ascontiguousarray(global_js_divergences, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'global_js_divergences' must be an array of np.float64: {_error}", argument="global_js_divergences") from None
    if global_js_divergences.ndim != 1:
        raise _ToxArgumentError(f"'global_js_divergences' must have 1 dimension, but has {global_js_divergences.ndim}", argument="global_js_divergences")
    try:
        total_included_n_reps_per_f = _np.ascontiguousarray(total_included_n_reps_per_f, dtype=_np.int32)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'total_included_n_reps_per_f' must be an array of np.int32: {_error}", argument="total_included_n_reps_per_f") from None
    if total_included_n_reps_per_f.ndim != 1:
        raise _ToxArgumentError(f"'total_included_n_reps_per_f' must have 1 dimension, but has {total_included_n_reps_per_f.ndim}", argument="total_included_n_reps_per_f")

    # what the inputs already say, rather than asking for it again
    k_families = global_js_divergences.shape[0]

    # Fortran cannot check that shared extents agree; this can
    if total_included_n_reps_per_f.shape[0] != k_families:
        raise _ToxArgumentError(f"'total_included_n_reps_per_f' has {total_included_n_reps_per_f.shape[0]} along axis 0, but 'global_js_divergences' implies k_families == {k_families}", argument="total_included_n_reps_per_f")

    # outputs and work arrays, which the caller never sees
    support_weights = _np.empty((k_families,), dtype=_np.float64, order='C')
    contribution_scores = _np.empty((k_families,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.fjct_compute_contribution_scores_c(
        global_js_divergences,
        total_included_n_reps_per_f,
        _ctypes.byref(_ctypes.c_int(k_families)),
        support_weights,
        contribution_scores,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _FJCT_COMPUTE_CONTRIBUTION_SCORES_ARGUMENTS, _FJCT_COMPUTE_CONTRIBUTION_SCORES_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    support_weights.flags.writeable = False
    contribution_scores.flags.writeable = False

    return {
        "support_weights": support_weights,
        "contribution_scores": contribution_scores,
    }
