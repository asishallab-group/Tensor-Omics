r"""f42_kd_tree

k-d tree spatial index over fixed-dimensional point sets.
Builds a k-d tree by recursively partitioning `kd_indices` around the median point along a
caller-supplied, cycling dimension order, using a stack-based (non-recursive) traversal so it
is safe to call from `pure` procedures. The tree is stored implicitly as an in-place-permuted
index array rather than as linked nodes.

Python binding, generated from f42_kd_tree. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.build_kd_index_c.restype = None
_lib.build_kd_index_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_KD_INDEX_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_KD_INDEX_ARGUMENT_SOURCES = (None, "points", "points", None, None, None,)

_lib.build_spherical_kd_c.restype = None
_lib.build_spherical_kd_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BUILD_SPHERICAL_KD_ARGUMENTS = ("points", "n_dimensions", "n_points", "kd_indices", "dimension_order", "ierr",)
#: For a derived argument, the one the caller passed it in
_BUILD_SPHERICAL_KD_ARGUMENT_SOURCES = (None, "points", "points", None, None, None,)

def build_kd_index(
        points,
        dimension_order,
):
    r"""Build a k-d tree index using a stack-based, non-recursive approach

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Data points
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order (by variance)
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.

    Returns
    -------
    kd_indices : np.ndarray[np.int32] of shape (n_points,), read-only
        Output index array (k-d tree order)
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
    Generated from the Fortran procedure `f42_kd_tree::build_kd_index`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")

    # outputs and work arrays, which the caller never sees
    kd_indices = _np.empty((n_points,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_kd_index_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BUILD_KD_INDEX_ARGUMENTS, _BUILD_KD_INDEX_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    kd_indices.flags.writeable = False

    return kd_indices

def build_spherical_kd(
        points,
        dimension_order,
):
    r"""Build a k-d tree index over points on the unit sphere (unit vectors)

    This is a thin, semantically-named wrapper: partitioning is identical to
    :func:`tensor_omics.build_kd_index` (plain per-axis median splits);
    callers are responsible for ensuring `points` are actually unit-normalized beforehand.

    Parameters
    ----------
    points : np.ndarray[np.float64] of shape (n_dimensions, n_points,), column-major (order='F')
        Data points
    dimension_order : np.ndarray[np.int32] of shape (n_dimensions,)
        Dimension order (by variance)
        The minimum valid value is `1`.
        The maximum valid value is `n_dimensions`.

    Returns
    -------
    kd_indices : np.ndarray[np.int32] of shape (n_points,), read-only
        Output index array (k-d tree order)
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
    Generated from the Fortran procedure `f42_kd_tree::build_spherical_kd`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    points = _contract.reals(points, "points", ndim=2, order="F")
    dimension_order = _contract.integers(dimension_order, "dimension_order", ndim=1)

    # what the inputs already say, rather than asking for it again
    n_dimensions = points.shape[0]
    n_points = points.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if dimension_order.shape[0] != n_dimensions:
        raise _ToxArgumentError(f"'dimension_order' has {dimension_order.shape[0]} along axis 0, but 'points' implies n_dimensions == {n_dimensions}; 'dimension_order' is declared with shape (n_dimensions)", argument="dimension_order")

    # outputs and work arrays, which the caller never sees
    kd_indices = _np.empty((n_points,), dtype=_np.int32, order='C')
    ierr = _ctypes.c_int(0)

    _lib.build_spherical_kd_c(
        points,
        _ctypes.byref(_ctypes.c_int(n_dimensions)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        kd_indices,
        dimension_order,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BUILD_SPHERICAL_KD_ARGUMENTS, _BUILD_SPHERICAL_KD_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    kd_indices.flags.writeable = False

    return kd_indices
