r"""tox_data_integration_preprocessing_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing

The step that turns expression vectors into the neighborhood residuals the rest of the test
consumes: gene-wise means, the signed deviation of each replicate from them, and the
neighborhoods of reference points those residuals are grouped into so the comparison is
conditioned on expression level rather than pooled across it.

`calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own.

Python binding, generated from tox_data_integration_preprocessing_impl. Do not edit.
"""

import ctypes
import os

import numpy as np

from .error_handling import check_err_code
from .library import load_library, nullable

_lib = load_library()

_lib.calc_neighborhood_size_c.restype = None
_lib.calc_neighborhood_size_c.argtypes = (
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='C_CONTIGUOUS'),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
    ctypes.POINTER(ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CALC_NEIGHBORHOOD_SIZE_ARGUMENTS = ("n_pool", "n_points", "n_genes_S", "mean_S", "desired_size",)
#: For a derived argument, the one the caller passed it in
_CALC_NEIGHBORHOOD_SIZE_ARGUMENT_SOURCES = (None, None, "mean_S", None, None,)

def calc_neighborhood_size(
        n_pool,
        n_points,
        mean_S,
        desired_size=1000,
):
    r"""Calculate the number of neighbors to be used for constructing neighborhoods

    Parameters
    ----------
    n_pool : int
        Total number of pooled mean-expression values across both studies
    n_points : int
        Number of reference points
    mean_S : np.ndarray[np.float64] of shape (n_genes_S,)
        Per-gene mean expression values
        NaN is permitted for this value.
    desired_size : int, optional, default 1000
        Optional desired neighborhood size
        The default value is `1000`.

    Returns
    -------
    n_neighbors : int
        Calculated neighborhood size

    Raises
    ------
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing_impl::calc_neighborhood_size`, whose argument names are
    the ones an error message reports.
    """
    # accept anything array-like, converting only when C needs it
    try:
        mean_S = np.ascontiguousarray(mean_S, dtype=np.float64)
    except (TypeError, ValueError) as error:
        raise TypeError(f"'mean_S' must be an array of np.float64: {error}") from None
    if mean_S.ndim != 1:
        raise ValueError(f"'mean_S' must have 1 dimension, but has {mean_S.ndim}")

    # what the inputs already say, rather than asking for it again
    n_genes_S = mean_S.shape[0]

    # outputs and work arrays, which the caller never sees
    n_neighbors = ctypes.c_int(0)
    ierr = ctypes.c_int(0)

    _lib.calc_neighborhood_size_c(
        ctypes.byref(ctypes.c_int(n_pool)),
        ctypes.byref(ctypes.c_int(n_points)),
        ctypes.byref(ctypes.c_int(n_genes_S)),
        mean_S,
        ctypes.byref(ctypes.c_int(desired_size)),
        ctypes.byref(n_neighbors),
        ctypes.byref(ierr),
    )

    check_err_code(ierr.value, _CALC_NEIGHBORHOOD_SIZE_ARGUMENTS, _CALC_NEIGHBORHOOD_SIZE_ARGUMENT_SOURCES)

    return n_neighbors.value
