r"""tox_data_integration_stats

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test

A permutation test estimating an empirical p-value for the weighted global JSD. Under the null
hypothesis that both studies are exchangeable, S1/S2 residuals are repeatedly shuffled within
each reference point's pooled neighborhood and the JSD is recomputed, giving a null
distribution against which the observed value is compared.

Work copies are shuffled, never the caller's arrays.

Python binding, generated from tox_data_integration_stats. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.gjct_permutation_test_c.restype = None
_lib.gjct_permutation_test_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=3, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _nullable(_ctypes.POINTER(_ctypes.c_int)),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    _nullable(_np.ctypeslib.ndpointer(dtype=_np.bool_, ndim=2, flags='F_CONTIGUOUS')),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_ARGUMENTS = ("neighborhood_residuals_S1", "neighborhood_residuals_S2", "n_reps_S1", "n_reps_S2", "n_neighbors", "n_points", "global_jsd_observed", "n_bins", "shared_residual_range", "n_permutations", "jsd_null", "p_value", "random_seed", "neighbor_mask_S1", "neighbor_mask_S2", "ierr",)
#: For a derived argument, the one the caller passed it in
_GJCT_PERMUTATION_TEST_ARGUMENT_SOURCES = (None, None, "neighborhood_residuals_S1", "neighborhood_residuals_S2", "neighborhood_residuals_S1", "neighborhood_residuals_S1", None, None, None, "jsd_null", None, None, None, None, None, None,)

def gjct_permutation_test(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        global_jsd_observed,
        n_bins,
        shared_residual_range,
        n_permutations,
        random_seed=None,
        neighbor_mask_S1=None,
        neighbor_mask_S2=None,
):
    r"""Estimate how likely the observed divergence is to occur by chance

    Tests the null hypothesis that both studies are exchangeable. The residuals are shuffled in
    the work copies, so the caller's own arrays are left untouched.

    Parameters
    ----------
    neighborhood_residuals_S1 : np.ndarray[np.float64] of shape (n_reps_S1, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 1, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    neighborhood_residuals_S2 : np.ndarray[np.float64] of shape (n_reps_S2, n_neighbors, n_points,), column-major (order='F')
        Computed neighborhood residuals for study 2, NaN is explicitly allowed for missing values
        NaN is permitted for this value.
    global_jsd_observed : float
        Observed global JSD value for both studies
    n_bins : int
        Number of equally sized histogram bins used for the studies
    shared_residual_range : float
        Computed residual range for both studies
        The minimum valid value is `0.0`.
    n_permutations : int
        Number of permutations to perform
    random_seed : int, optional
        Seed to use for shuffling
    neighbor_mask_S1 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 1 (e.g. for family-wise analysis)
    neighbor_mask_S2 : np.ndarray[np.bool_] of shape (n_neighbors, n_points,), column-major (order='F'), optional
        Optional mask to exclude specific neighbors from study 2 (e.g. for family-wise analysis)

    Returns
    -------
    dict
        with keys:

        jsd_null : np.ndarray[np.float64] of shape (n_permutations,), read-only
            Vector of global divergence values obtained under the null hypothesis
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_value : float
            Empirical p-value of the permutation test: \( \frac{\text{count}(jsd\_null \ge global\_jsd\_observed) + 1}{n\_permutations + 1} \)

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats::gjct_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    neighborhood_residuals_S1 = _contract.reals(neighborhood_residuals_S1, "neighborhood_residuals_S1", ndim=3, order="F")
    neighborhood_residuals_S2 = _contract.reals(neighborhood_residuals_S2, "neighborhood_residuals_S2", ndim=3, order="F")
    global_jsd_observed = _contract.real(global_jsd_observed, "global_jsd_observed")
    n_bins = _contract.integer(n_bins, "n_bins")
    shared_residual_range = _contract.real(shared_residual_range, "shared_residual_range")
    n_permutations = _contract.integer(n_permutations, "n_permutations")
    if random_seed is not None:
        random_seed = _contract.integer(random_seed, "random_seed")
    if neighbor_mask_S1 is not None:
        neighbor_mask_S1 = _contract.logicals(neighbor_mask_S1, "neighbor_mask_S1", ndim=2, order="F")
    if neighbor_mask_S2 is not None:
        neighbor_mask_S2 = _contract.logicals(neighbor_mask_S2, "neighbor_mask_S2", ndim=2, order="F")

    # what the inputs already say, rather than asking for it again
    n_reps_S1 = neighborhood_residuals_S1.shape[0]
    n_reps_S2 = neighborhood_residuals_S2.shape[0]
    n_neighbors = neighborhood_residuals_S1.shape[1]
    n_points = neighborhood_residuals_S1.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if neighborhood_residuals_S2.shape[1] != n_neighbors:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighborhood_residuals_S2' is declared with shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")
    if neighborhood_residuals_S2.shape[2] != n_points:
        raise _ToxArgumentError(f"'neighborhood_residuals_S2' has {neighborhood_residuals_S2.shape[2]} along axis 2, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighborhood_residuals_S2' is declared with shape (n_reps_S2, n_neighbors, n_points)", argument="neighborhood_residuals_S2")
    if neighbor_mask_S1 is not None and neighbor_mask_S1.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighbor_mask_S1' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S1")
    if neighbor_mask_S1 is not None and neighbor_mask_S1.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask_S1' has {neighbor_mask_S1.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighbor_mask_S1' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S1")
    if neighbor_mask_S2 is not None and neighbor_mask_S2.shape[0] != n_neighbors:
        raise _ToxArgumentError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[0]} along axis 0, but 'neighborhood_residuals_S1' implies n_neighbors == {n_neighbors}; 'neighbor_mask_S2' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S2")
    if neighbor_mask_S2 is not None and neighbor_mask_S2.shape[1] != n_points:
        raise _ToxArgumentError(f"'neighbor_mask_S2' has {neighbor_mask_S2.shape[1]} along axis 1, but 'neighborhood_residuals_S1' implies n_points == {n_points}; 'neighbor_mask_S2' is declared with shape (n_neighbors, n_points)", argument="neighbor_mask_S2")

    # outputs and work arrays, which the caller never sees
    jsd_null = _np.empty((n_permutations,), dtype=_np.float64, order='C')
    p_value = _ctypes.c_double(0)
    ierr = _ctypes.c_int(0)

    _lib.gjct_permutation_test_c(
        neighborhood_residuals_S1,
        neighborhood_residuals_S2,
        _ctypes.byref(_ctypes.c_int(n_reps_S1)),
        _ctypes.byref(_ctypes.c_int(n_reps_S2)),
        _ctypes.byref(_ctypes.c_int(n_neighbors)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_double(global_jsd_observed)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_double(shared_residual_range)),
        _ctypes.byref(_ctypes.c_int(n_permutations)),
        jsd_null,
        _ctypes.byref(p_value),
        None if random_seed is None else _ctypes.byref(_ctypes.c_int(random_seed)),
        neighbor_mask_S1,
        neighbor_mask_S2,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_ARGUMENTS, _GJCT_PERMUTATION_TEST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    jsd_null.flags.writeable = False

    return {
        "jsd_null": jsd_null,
        "p_value": p_value.value,
    }
