r"""tox_data_integration_stats

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Bootstrap

The bootstrap
:func:`tensor_omics.bootstrap_histogram` gives each
study's weighted global JSD an interval, by drawing each study's histograms again with
replacement from the study's own counts: a robustness interval, how much the observed value
would move if each study's residuals were drawn again. `calc_js_comp_test_n_top_k_jsds` sizes
its top-/bottom-k heaps for a caller that allocates its own. The significance test of the
JS-Comp-Test is the study-label Permutation-Test,
:func:`tensor_omics.run_js_comp_test_permutation`
and its adaptive counterpart.

Python binding, generated from tox_data_integration_stats. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.bootstrap_histogram_c.restype = None
_lib.bootstrap_histogram_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=3, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BOOTSTRAP_HISTOGRAM_ARGUMENTS = ("n_bootstraps", "n_bins", "n_points", "n_studies", "counts", "mean_pmf_included_n_reps", "included_n_reps", "confidence_interval", "two_sided_bootstrapping_significance_level", "random_seed", "divergence_mode", "ierr",)
#: For a derived argument, the one the caller passed it in
_BOOTSTRAP_HISTOGRAM_ARGUMENT_SOURCES = (None, "counts", "counts", "counts", None, None, None, None, None, None, None, None,)

def bootstrap_histogram(
        n_bootstraps,
        counts,
        mean_pmf_included_n_reps,
        included_n_reps,
        confidence_interval,
        two_sided_bootstrapping_significance_level=2.5,
        random_seed=42,
        divergence_mode='divergence_pairwise_k2',
):
    r"""Bootstrap an interval for each study's global JSD by resampling the studies' histograms

    In every bootstrap, each study's histogram at each reference point is drawn again with
    replacement (multinomial), with as many residuals as the study observed
    (`included_n_reps`); the consensus pmf, each study's per-point JSD and its weighted global
    JSD are recomputed from the draw, with the divergence `divergence_mode` selects (by default
    for exactly two studies their pairwise JSD, otherwise against the consensus). Each study is
    drawn from its own counts (`counts`): the interval shows how much the observed JSD would move
    if each study's own residuals were drawn again -- a robustness interval of the observed
    value. Drawing again adds sampling noise, so the bootstrapped JSDs are larger than the
    observed one on average: when the studies hardly differ, the interval can lie entirely above
    the observed JSD.

    `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
    slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
    values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
    smallest of the n_bootstrapping_top_k_jsds largest values]`.

    Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
    draw, when at some reference point the studies' `included_n_reps` do not add up to
    `mean_pmf_included_n_reps` (as
    :func:`tensor_omics.create_mean_pmf` produces
    them), or a study's `counts` at some reference point do not add up to its
    `included_n_reps`; `ERR_ALLOC_FAIL` when the random number generator cannot be
    created. The count checks, together with the wrappers' bounds on the counts, establish every
    draw's preconditions, so no draw can fail later on.

    Parameters
    ----------
    n_bootstraps : int
        Number of bootstrap resamples to perform
        The minimum valid value is `1`.
    counts : np.ndarray[np.int32] of shape (n_points, n_bins, n_studies,), column-major (order='F')
        Every study's own residual counts per bin: what each study is drawn from; at each
        reference point they must add up to the study's `included_n_reps` there
        The minimum valid value is `0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates per reference point of all studies together (the
        consensus's count, the denominator of the weights); must be the sum of the studies'
        `included_n_reps` there
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates per reference point, per study: the size of each study's
        draw. Every study needs at least one at every reference point
        The minimum valid value is `1`.
    confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), modified in place
        On entry the seed values of the heaps (see above); on return the bootstrapped
        `[lower, upper]` interval per study
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    two_sided_bootstrapping_significance_level : float, optional, default 2.5
        Only sizes n_bootstrapping_top_k_jsds (via calc_js_comp_test_n_top_k_jsds)
        The minimum valid value is `0.0`.
        The maximum valid value is `100.0`.
        The default value is `2.5`.
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.
    divergence_mode : str, one of 'divergence_symmetric' | 'divergence_pairwise_k2', optional, default 'divergence_pairwise_k2'
        How each study's per-point divergence of every bootstrap draw is computed, as in the
        routine whose interval this is

        The default value is `'divergence_pairwise_k2'`.
        Trailing blanks are not part of the value, as in Fortran: they are dropped.

    Returns
    -------
    None

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_stats::bootstrap_histogram`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_bootstraps = _contract.integer(n_bootstraps, "n_bootstraps")
    counts = _contract.integers(counts, "counts", ndim=3, order="F")
    mean_pmf_included_n_reps = _contract.integers(mean_pmf_included_n_reps, "mean_pmf_included_n_reps", ndim=1)
    included_n_reps = _contract.integers(included_n_reps, "included_n_reps", ndim=2, order="F")
    if not isinstance(confidence_interval, _np.ndarray) or confidence_interval.dtype != _np.float64:
        raise _ToxArgumentError("'confidence_interval' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="confidence_interval")
    if confidence_interval.ndim != 2:
        raise _ToxArgumentError(f"'confidence_interval' must have 2 dimensions, but has {confidence_interval.ndim}", argument="confidence_interval")
    if not confidence_interval.flags.f_contiguous:
        raise _ToxArgumentError("'confidence_interval' is modified in place, so it must already be column-major (order='F')", argument="confidence_interval")
    two_sided_bootstrapping_significance_level = _contract.real(two_sided_bootstrapping_significance_level, "two_sided_bootstrapping_significance_level")
    random_seed = _contract.integer(random_seed, "random_seed")
    divergence_mode = _np.array([_contract.mode(divergence_mode, "divergence_mode", ("divergence_symmetric", "divergence_pairwise_k2",)).ljust(22)], dtype="S22")

    # what the inputs already say, rather than asking for it again
    n_bins = counts.shape[1]
    n_points = counts.shape[0]
    n_studies = counts.shape[2]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if mean_pmf_included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'mean_pmf_included_n_reps' has {mean_pmf_included_n_reps.shape[0]} along axis 0, but 'counts' implies n_points == {n_points}; 'mean_pmf_included_n_reps' is declared with shape (n_points)", argument="mean_pmf_included_n_reps")
    if included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but 'counts' implies n_points == {n_points}; 'included_n_reps' is declared with shape (n_points, n_studies)", argument="included_n_reps")
    if included_n_reps.shape[1] != n_studies:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[1]} along axis 1, but 'counts' implies n_studies == {n_studies}; 'included_n_reps' is declared with shape (n_points, n_studies)", argument="included_n_reps")
    if confidence_interval.shape[0] != 2:
        raise _ToxArgumentError(f"'confidence_interval' has {confidence_interval.shape[0]} along axis 0, but its declared shape fixes it at {2}; 'confidence_interval' is declared with shape (2, n_studies)", argument="confidence_interval")
    if confidence_interval.shape[1] != n_studies:
        raise _ToxArgumentError(f"'confidence_interval' has {confidence_interval.shape[1]} along axis 1, but 'counts' implies n_studies == {n_studies}; 'confidence_interval' is declared with shape (2, n_studies)", argument="confidence_interval")

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.bootstrap_histogram_c(
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_studies)),
        counts,
        mean_pmf_included_n_reps,
        included_n_reps,
        confidence_interval,
        _ctypes.byref(_ctypes.c_double(two_sided_bootstrapping_significance_level)),
        _ctypes.byref(_ctypes.c_int(random_seed)),
        divergence_mode,
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BOOTSTRAP_HISTOGRAM_ARGUMENTS, _BOOTSTRAP_HISTOGRAM_ARGUMENT_SOURCES)

    return None
