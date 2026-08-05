"""tox_tissue_versatility

Generated from the kernel; do not edit -- regenerate instead.

Python binding, generated from tox_tissue_versatility. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.compute_tissue_versatility_c.restype = None
_lib.compute_tissue_versatility_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.bool_, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
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
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_tissue_versatility::compute_tissue_versatility`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        expression_vectors = np.asfortranarray(expression_vectors, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'expression_vectors' must be an array of np.float64: {error}") from None
    if expression_vectors.ndim != 2:
        raise ValueError(f"'expression_vectors' must have 2 dimensions, but has {expression_vectors.ndim}")
    try:
        vectors_selection_mask = np.ascontiguousarray(vectors_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'vectors_selection_mask' must be an array of np.bool_: {error}") from None
    if vectors_selection_mask.ndim != 1:
        raise ValueError(f"'vectors_selection_mask' must have 1 dimension, but has {vectors_selection_mask.ndim}")
    try:
        axes_selection_mask = np.ascontiguousarray(axes_selection_mask, dtype=np.bool_)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'axes_selection_mask' must be an array of np.bool_: {error}") from None
    if axes_selection_mask.ndim != 1:
        raise ValueError(f"'axes_selection_mask' must have 1 dimension, but has {axes_selection_mask.ndim}")

    # what the inputs already say, rather than asking for it again
    n_axes = expression_vectors.shape[0]
    n_vectors = expression_vectors.shape[1]
    n_selected_vectors = int(vectors_selection_mask.sum())
    n_selected_axes = int(axes_selection_mask.sum())

    # Fortran cannot check that shared extents agree; this can
    if axes_selection_mask.shape[0] != n_axes:
        raise ValueError(f"'axes_selection_mask' has {axes_selection_mask.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_axes == {n_axes}"
        )
    if vectors_selection_mask.shape[0] != n_vectors:
        raise ValueError(f"'vectors_selection_mask' has {vectors_selection_mask.shape[0]} along axis 0, but "
            f"'expression_vectors' implies n_vectors == {n_vectors}"
        )

    # outputs and work arrays, which the caller never sees
    tissue_versatilities = np.empty((n_selected_vectors,), dtype=np.float64, order='C')
    tissue_angles_deg = np.empty((n_selected_vectors,), dtype=np.float64, order='C')
    ierr = ctypes.c_int(0)

    _lib.compute_tissue_versatility_c(
        ctypes.byref(ctypes.c_int(n_axes)),
        ctypes.byref(ctypes.c_int(n_vectors)),
        expression_vectors,
        vectors_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_vectors)),
        axes_selection_mask,
        ctypes.byref(ctypes.c_int(n_selected_axes)),
        tissue_versatilities,
        tissue_angles_deg,
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _COMPUTE_TISSUE_VERSATILITY_ARGUMENTS, _COMPUTE_TISSUE_VERSATILITY_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    tissue_versatilities.flags.writeable = False
    tissue_angles_deg.flags.writeable = False

    return {
        "tissue_versatilities": tissue_versatilities,
        "tissue_angles_deg": tissue_angles_deg,
    }
