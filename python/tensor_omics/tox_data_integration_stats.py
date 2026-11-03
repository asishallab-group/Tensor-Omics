r"""tox_data_integration_stats

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Permutation Test and Bootstrap

The permutation test of the JS-Comp-Test:
:func:`tensor_omics.gjct_permutation_test` estimates an
empirical p-value for each of K studies' weighted global JSD against the consensus pmf, by
resampling every reference point's pooled histogram counts into the studies without
replacement. `p_values` is the primary result; `p_values_observed_consensus` uses the older
null (a fixed observed consensus) and is kept for comparison.

The bootstrap
:func:`tensor_omics.bootstrap_histogram` gives each
study's weighted global JSD a confidence interval, by drawing each study's histograms with
replacement from the pooled consensus counts; `calc_js_comp_test_n_top_k_jsds` sizes its top-/bottom-k
heaps for a caller that allocates its own.

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
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_GJCT_PERMUTATION_TEST_ARGUMENTS = ("n_permutations", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf", "mean_pmf_included_n_reps", "included_n_reps", "global_jsd_observed", "p_values", "p_values_observed_consensus", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_GJCT_PERMUTATION_TEST_ARGUMENT_SOURCES = (None, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", None, None, None, None, None, None, None, None, None,)

_lib.bootstrap_histogram_c.restype = None
_lib.bootstrap_histogram_c.argtypes = (
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_BOOTSTRAP_HISTOGRAM_ARGUMENTS = ("n_bootstraps", "n_bins", "n_points", "n_studies", "mean_pmf_counts", "mean_pmf_included_n_reps", "included_n_reps", "confidence_interval", "two_sided_bootstrapping_significance_level", "random_seed", "ierr",)
#: For a derived argument, the one the caller passed it in
_BOOTSTRAP_HISTOGRAM_ARGUMENT_SOURCES = (None, "mean_pmf_counts", "mean_pmf_counts", "included_n_reps", None, None, None, None, None, None, None,)

def gjct_permutation_test(
        n_permutations,
        mean_pmf_counts,
        mean_pmf,
        mean_pmf_included_n_reps,
        included_n_reps,
        global_jsd_observed,
        random_seed=42,
):
    r"""Estimate how likely each study's observed weighted global JSD is to occur by chance

    In every permutation, each reference point's pooled consensus counts (`mean_pmf_counts`) are
    drawn without replacement (multivariate hypergeometric) into the studies, each study taking
    as many residuals as it observed (`included_n_reps`). The permuted pmfs, their equal-weight
    consensus (built exactly as
    :func:`tensor_omics.create_mean_pmf` builds the
    observed one) and each study's weighted global JSD against that consensus are recomputed, and
    compared with the observed `global_jsd_observed`. `p_values_observed_consensus` comes from the
    same draws, taking each permuted study's JSD against the fixed observed `mean_pmf` instead.

    Both p-values are `(count + 1)/(n_permutations + 1)`, so never exactly 0. With
    `n_permutations = 0` no permutation runs and both are NaN ("not computed").

    Errors (`ierr`): `ERR_INVALID_INPUT`, before any draw, when at some reference point the
    studies' counts do not add up to the pool's total (`mean_pmf_included_n_reps`) or that total
    is not the sum of the pool's bin counts; `ERR_ALLOC_FAIL` when the random number generator
    cannot be created. The p-values are then NaN. The count check, together with the wrappers'
    bounds on the counts, establishes every draw's preconditions, so no draw can fail later on.

    Parameters
    ----------
    n_permutations : int
        Number of permutations; `0` skips the test (both p-value outputs are then NaN)
        The minimum valid value is `0`.
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Pooled residual counts per bin of the consensus pmf: the pool each permutation draws from
        The minimum valid value is `0`.
    mean_pmf : np.ndarray[np.float64] of shape (n_points, n_bins,), column-major (order='F')
        Observed consensus pmf; used only for `p_values_observed_consensus`
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates (included ones) per reference point for the consensus pmf
        The minimum valid value is `0`.
    included_n_reps : np.ndarray[np.int32] of shape (n_points, n_studies,), column-major (order='F')
        Count of non-NaN replicates per reference point, per study: the size of each study's
        draw. Every study needs at least one at every reference point
        The minimum valid value is `1`.
    global_jsd_observed : np.ndarray[np.float64] of shape (n_studies,)
        Observed weighted global JSD of each study against the consensus pmf
    random_seed : int, optional, default 42
        Seed for the GSL random number generator
        The default value is `42`.

    Returns
    -------
    dict
        with keys:

        p_values : np.ndarray[np.float64] of shape (n_studies,), read-only
            Empirical p-value per study, `(count + 1)/(n_permutations + 1)` over the permutations
            whose global JSD reached the observed one; NaN when `n_permutations` is 0
            A result is a value; call `.copy()` to obtain a modifiable array.
        p_values_observed_consensus : np.ndarray[np.float64] of shape (n_studies,), read-only
            Like `p_values`, from the same permutations, but against the fixed observed consensus
            `mean_pmf`; kept for comparison only. NaN when `n_permutations` is 0
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
    Generated from the Fortran procedure `tox_data_integration_stats::gjct_permutation_test`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    n_permutations = _contract.integer(n_permutations, "n_permutations")
    mean_pmf_counts = _contract.integers(mean_pmf_counts, "mean_pmf_counts", ndim=2, order="F")
    mean_pmf = _contract.reals(mean_pmf, "mean_pmf", ndim=2, order="F")
    mean_pmf_included_n_reps = _contract.integers(mean_pmf_included_n_reps, "mean_pmf_included_n_reps", ndim=1)
    included_n_reps = _contract.integers(included_n_reps, "included_n_reps", ndim=2, order="F")
    global_jsd_observed = _contract.reals(global_jsd_observed, "global_jsd_observed", ndim=1)
    random_seed = _contract.integer(random_seed, "random_seed")

    # what the inputs already say, rather than asking for it again
    n_bins = mean_pmf_counts.shape[1]
    n_points = mean_pmf_counts.shape[0]
    n_studies = included_n_reps.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if mean_pmf.shape[0] != n_points:
        raise _ToxArgumentError(f"'mean_pmf' has {mean_pmf.shape[0]} along axis 0, but 'mean_pmf_counts' implies n_points == {n_points}; 'mean_pmf' is declared with shape (n_points, n_bins)", argument="mean_pmf")
    if mean_pmf.shape[1] != n_bins:
        raise _ToxArgumentError(f"'mean_pmf' has {mean_pmf.shape[1]} along axis 1, but 'mean_pmf_counts' implies n_bins == {n_bins}; 'mean_pmf' is declared with shape (n_points, n_bins)", argument="mean_pmf")
    if mean_pmf_included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'mean_pmf_included_n_reps' has {mean_pmf_included_n_reps.shape[0]} along axis 0, but 'mean_pmf_counts' implies n_points == {n_points}; 'mean_pmf_included_n_reps' is declared with shape (n_points)", argument="mean_pmf_included_n_reps")
    if included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but 'mean_pmf_counts' implies n_points == {n_points}; 'included_n_reps' is declared with shape (n_points, n_studies)", argument="included_n_reps")
    if global_jsd_observed.shape[0] != n_studies:
        raise _ToxArgumentError(f"'global_jsd_observed' has {global_jsd_observed.shape[0]} along axis 0, but 'included_n_reps' implies n_studies == {n_studies}; 'global_jsd_observed' is declared with shape (n_studies)", argument="global_jsd_observed")

    # outputs and work arrays, which the caller never sees
    p_values = _np.empty((n_studies,), dtype=_np.float64, order='C')
    p_values_observed_consensus = _np.empty((n_studies,), dtype=_np.float64, order='C')
    ierr = _ctypes.c_int(0)

    _lib.gjct_permutation_test_c(
        _ctypes.byref(_ctypes.c_int(n_permutations)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_studies)),
        mean_pmf_counts,
        mean_pmf,
        mean_pmf_included_n_reps,
        included_n_reps,
        global_jsd_observed,
        p_values,
        p_values_observed_consensus,
        _ctypes.byref(_ctypes.c_int(random_seed)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _GJCT_PERMUTATION_TEST_ARGUMENTS, _GJCT_PERMUTATION_TEST_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    p_values.flags.writeable = False
    p_values_observed_consensus.flags.writeable = False

    return {
        "p_values": p_values,
        "p_values_observed_consensus": p_values_observed_consensus,
    }

def bootstrap_histogram(
        n_bootstraps,
        mean_pmf_counts,
        mean_pmf_included_n_reps,
        included_n_reps,
        confidence_interval,
        two_sided_bootstrapping_significance_level=2.5,
        random_seed=42,
):
    r"""Bootstrap a confidence interval for each study's global JSD by resampling the pooled consensus histogram

    In every bootstrap, each study's histogram at each reference point is drawn with replacement
    (multinomial) from the pooled consensus counts `mean_pmf_counts`, with as many residuals as
    the study observed; the consensus pmf, each study's JSD against it and its weighted global
    JSD are recomputed from the draw.

    `confidence_interval` is input and output: its incoming `[lower, upper]` values seed every
    slot of the bottom-k/top-k heaps, so they are only displaced by more extreme bootstrap
    values. On return it holds `[largest of the n_bootstrapping_top_k_jsds smallest values,
    smallest of the n_bootstrapping_top_k_jsds largest values]`.

    Errors (`ierr`), with `confidence_interval` left unchanged: `ERR_INVALID_INPUT`, before any
    draw, when at some reference point the studies' counts do not add up to
    `mean_pmf_included_n_reps` or that total is not the sum of `mean_pmf_counts` (as
    :func:`tensor_omics.create_mean_pmf` produces
    them); `ERR_ALLOC_FAIL` when the random number generator cannot be created. The count check,
    together with the wrappers' bounds on the counts, establishes every draw's preconditions, so
    no draw can fail later on.

    Parameters
    ----------
    n_bootstraps : int
        Number of bootstrap resamples to perform
        The minimum valid value is `1`.
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Pooled residual counts per bin of the consensus pmf, the pool every bootstrap draws from
        The minimum valid value is `0`.
    mean_pmf_included_n_reps : np.ndarray[np.int32] of shape (n_points,)
        Count of non-NaN replicates per reference point of the pool
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
    mean_pmf_counts = _contract.integers(mean_pmf_counts, "mean_pmf_counts", ndim=2, order="F")
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

    # what the inputs already say, rather than asking for it again
    n_bins = mean_pmf_counts.shape[1]
    n_points = mean_pmf_counts.shape[0]
    n_studies = included_n_reps.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if mean_pmf_included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'mean_pmf_included_n_reps' has {mean_pmf_included_n_reps.shape[0]} along axis 0, but 'mean_pmf_counts' implies n_points == {n_points}; 'mean_pmf_included_n_reps' is declared with shape (n_points)", argument="mean_pmf_included_n_reps")
    if included_n_reps.shape[0] != n_points:
        raise _ToxArgumentError(f"'included_n_reps' has {included_n_reps.shape[0]} along axis 0, but 'mean_pmf_counts' implies n_points == {n_points}; 'included_n_reps' is declared with shape (n_points, n_studies)", argument="included_n_reps")
    if confidence_interval.shape[0] != 2:
        raise _ToxArgumentError(f"'confidence_interval' has {confidence_interval.shape[0]} along axis 0, but its declared shape fixes it at {2}; 'confidence_interval' is declared with shape (2, n_studies)", argument="confidence_interval")
    if confidence_interval.shape[1] != n_studies:
        raise _ToxArgumentError(f"'confidence_interval' has {confidence_interval.shape[1]} along axis 1, but 'included_n_reps' implies n_studies == {n_studies}; 'confidence_interval' is declared with shape (2, n_studies)", argument="confidence_interval")

    # outputs and work arrays, which the caller never sees
    ierr = _ctypes.c_int(0)

    _lib.bootstrap_histogram_c(
        _ctypes.byref(_ctypes.c_int(n_bootstraps)),
        _ctypes.byref(_ctypes.c_int(n_bins)),
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(n_studies)),
        mean_pmf_counts,
        mean_pmf_included_n_reps,
        included_n_reps,
        confidence_interval,
        _ctypes.byref(_ctypes.c_double(two_sided_bootstrapping_significance_level)),
        _ctypes.byref(_ctypes.c_int(random_seed)),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _BOOTSTRAP_HISTOGRAM_ARGUMENTS, _BOOTSTRAP_HISTOGRAM_ARGUMENT_SOURCES)

    return None
