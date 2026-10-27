r"""tox_data_integration_preprocessing_impl

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Preprocessing

The step that turns expression vectors into the neighborhood residuals the rest of the test
consumes: gene-wise means, the signed deviation of each replicate from them, and the
neighborhoods of reference points those residuals are grouped into so the comparison is
conditioned on expression level rather than pooled across it.

`calc_neighborhood_size` sizes a neighborhood for a caller that allocates its own.

Python binding, generated from tox_data_integration_preprocessing_impl. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.calc_neighborhood_size_c.restype = None
_lib.calc_neighborhood_size_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
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

    The `desired_size` works as upper limit, as the actual neighborhood size might be lower
    due to few genes with non-NaN mean.

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
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_preprocessing_impl::calc_neighborhood_size`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_pool = _contract.integer(n_pool, "n_pool")
    n_points = _contract.integer(n_points, "n_points")
    mean_S = _contract.reals(mean_S, "mean_S", ndim=1)
    desired_size = _contract.integer(desired_size, "desired_size")

    # what the inputs already say, rather than asking for it again
    n_genes_S = mean_S.shape[0]

    # outputs and work arrays, which the caller never sees
    n_neighbors = _ctypes.c_int(0)
    ierr = _ctypes.c_int(0)

    _lib.calc_neighborhood_size_c(
        _ctypes.byref(_ctypes.c_int(n_pool)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_genes_S)),
        mean_S,
        _ctypes.byref(_ctypes.c_int(desired_size)),
        _ctypes.byref(n_neighbors),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CALC_NEIGHBORHOOD_SIZE_ARGUMENTS, _CALC_NEIGHBORHOOD_SIZE_ARGUMENT_SOURCES)

    return n_neighbors.value
