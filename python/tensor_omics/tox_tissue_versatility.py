r"""tox_tissue_versatility

Normalized tissue (axis) versatility: how uniformly a gene is expressed across tissues.

An angle-based metric. A gene expressed equally across every selected axis points along the
diagonal of that subspace and scores maximally versatile; one confined to a single tissue
points along that axis and scores minimally. Normalized, so scores over different numbers of
selected axes are comparable.

Python binding, generated from tox_tissue_versatility. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.compute_tissue_versatility_c.restype = None
_lib.compute_tissue_versatility_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_COMPUTE_TISSUE_VERSATILITY_ARGUMENTS = ("n_axes", "n_vectors", "expression_vectors", "vectors_selection_mask", "n_selected_vectors", "axes_selection_mask", "n_selected_axes", "tissue_versatilities", "tissue_angles_deg", "ierr",)
#: For a derived argument, the one the caller passed it in
_COMPUTE_TISSUE_VERSATILITY_ARGUMENT_SOURCES = ("expression_vectors", "expression_vectors", None, None, "tissue_versatilities", None, "axes_selection_mask", None, None, None,)

def compute_tissue_versatility(
        expression_vectors,
        vectors_selection_mask,
        axes_selection_mask,
):
    r"""Computes normalized tissue versatility for selected expression vectors.

    The metric is based on the angle between each gene expression vector and the space diagonal.
    Versatility is normalized to [0, 1], where 0 means uniform expression and 1 means expression in only one axis.

    The masks follow the `n_selected_` convention, so the generated wrapper validates that each
    selection count matches its mask; `n_selected_axes` (not an array extent) carries its own floor.

    Parameters
    ----------
    expression_vectors : np.ndarray[np.float64] of shape (n_axes, n_vectors,), column-major (order='F')
        2D array (n_axes, n_vectors), each column is a gene expression vector
    vectors_selection_mask : np.ndarray[np.bool_] of shape (n_vectors,)
        Logical array (n_vectors), True for vectors to process
    axes_selection_mask : np.ndarray[np.bool_] of shape (n_axes,)
        Logical array (n_axes), True for axes to include in calculation

    Returns
    -------
    dict
        with keys:

        tissue_versatilities : np.ndarray[np.float64] of shape (n_selected_vectors,), read-only
            Output, real array, length = n_selected_vectors, stores the calculated tissue versatilities
            A result is a value; call `.copy()` to obtain a modifiable array.
        tissue_angles_deg : np.ndarray[np.float64] of shape (n_selected_vectors,), read-only
            Output, real array, length = n_selected_vectors, stores the calculated angles in degrees
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
    Generated from the Fortran procedure `tox_tissue_versatility::compute_tissue_versatility`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = _np.asfortranarray(expression_vectors, dtype=_np.float64)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'expression_vectors' must be an array of np.float64: {_error}", argument="expression_vectors") from None
    if expression_vectors.ndim != 2:
        raise _ToxArgumentError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}", argument="expression_vectors")
    try:
        vectors_selection_mask = _np.ascontiguousarray(vectors_selection_mask, dtype=_np.bool_)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'vectors_selection_mask' must be an array of np.bool_: {_error}", argument="vectors_selection_mask") from None
    if vectors_selection_mask.ndim != 1:
        raise _ToxArgumentError(f"'vectors_selection_mask' must have 1 dimension, but has {vectors_selection_mask.ndim}", argument="vectors_selection_mask")
    try:
        axes_selection_mask = _np.ascontiguousarray(axes_selection_mask, dtype=_np.bool_)
    except (TypeError, ValueError) as _error:
        raise _ToxArgumentError(f"'axes_selection_mask' must be an array of np.bool_: {_error}", argument="axes_selection_mask") from None
    if axes_selection_mask.ndim != 1:
        raise _ToxArgumentError(f"'axes_selection_mask' must have 1 dimension, but has {axes_selection_mask.ndim}", argument="axes_selection_mask")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_vectors = expression_vectors.shape[1]
    n_selected_vectors = int(vectors_selection_mask.sum())
    n_selected_axes = int(axes_selection_mask.sum())

    # Fortran cannot check that an array has the shape it is declared with; this can
    if vectors_selection_mask.shape[0] != n_vectors:
        raise _ToxArgumentError(f"'vectors_selection_mask' has {vectors_selection_mask.shape[0]} along axis 0, but 'expression_vectors' implies n_vectors == {n_vectors}; 'vectors_selection_mask' has shape (n_vectors)", argument="vectors_selection_mask")
    if axes_selection_mask.shape[0] != n_axes:
        raise _ToxArgumentError(f"'axes_selection_mask' has {axes_selection_mask.shape[0]} along axis 0, but 'expression_vectors' implies n_axes == {n_axes}; 'axes_selection_mask' has shape (n_axes)", argument="axes_selection_mask")

    # outputs and work arrays, which the caller never sees
    tissue_versatilities = _np.empty((n_selected_vectors,), dtype=_np.float64, order='C')
    tissue_angles_deg = _np.empty((n_selected_vectors,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.compute_tissue_versatility_c(
        _ctypes.byref(_ctypes.c_int(n_axes)),
        _ctypes.byref(_ctypes.c_int(n_vectors)),
        expression_vectors,
        vectors_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_vectors)),
        axes_selection_mask,
        _ctypes.byref(_ctypes.c_int(n_selected_axes)),
        tissue_versatilities,
        tissue_angles_deg,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _COMPUTE_TISSUE_VERSATILITY_ARGUMENTS, _COMPUTE_TISSUE_VERSATILITY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    tissue_versatilities.flags.writeable = False
    tissue_angles_deg.flags.writeable = False

    return {
        "tissue_versatilities": tissue_versatilities,
        "tissue_angles_deg": tissue_angles_deg,
    }
