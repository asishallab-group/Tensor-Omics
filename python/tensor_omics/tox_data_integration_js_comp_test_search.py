r"""tox_data_integration_js_comp_test_search

# Jensen-Shannon-Divergence (JSD) Compatibility Test (gJCT) Parameter-Search Criteria

What the two parameter searches of the JS-Comp-Test have in common: the fixed-k search
:func:`tensor_omics.run_js_comp_test_parameter_search`
and the adaptive search
:func:`tensor_omics.run_js_comp_test_adaptive_parameter_search`.
Both accept a candidate when its neighborhoods overlap enough and every consensus bin holds
enough residuals, and both stop at the first such candidate whose bootstrapped JSDs have
reached a plateau.

The admissibility gates are
:func:`tensor_omics.check_neighborhood_overlaps`
(consecutive neighborhoods overlap by at least `min_neighbor_overlap`) and
:func:`tensor_omics.check_mean_pmf_min_counts`
(every consensus bin reaches the minimum count). The plateau checks are
:func:`tensor_omics.check_plateau_condition`
(the confidence-interval overlap with the best candidate so far, joined over the studies by
`join_method`) and
:func:`tensor_omics.check_effect_size_plateau_condition`
(the relative change of each study's observed JSD between successive candidates); a search's
`plateau_mode` selects which of them stops it. Both searches evaluate at most
``MAX_CANDIDATE_PAIRS`` candidates.

Python binding, generated from tox_data_integration_js_comp_test_search. Do not edit.
"""

import ctypes as _ctypes
import os as _os

import numpy as _np

from . import _contract as _contract
from .error_handling import ToxArgumentError as _ToxArgumentError, check_err_code as _check_err_code
from .library import load_library as _load_library, nullable as _nullable

_lib = _load_library()

_lib.check_neighborhood_overlaps_c.restype = None
_lib.check_neighborhood_overlaps_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENTS = ("neighborhood_range", "n_points", "min_neighbor_overlap", "all_have_min_neighbor_overlap", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENT_SOURCES = (None, "neighborhood_range", None, None, None,)

_lib.check_mean_pmf_min_counts_c.restype = None
_lib.check_mean_pmf_min_counts_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.int32, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENTS = ("mean_pmf_counts", "n_bins", "n_bins_per_point", "n_points", "min_count", "all_bins_have_min_count", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENT_SOURCES = (None, "mean_pmf_counts", None, "mean_pmf_counts", None, None, None,)

_lib.check_plateau_condition_c.restype = None
_lib.check_plateau_condition_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=2, flags='F_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(ndim=1),
    _ctypes.POINTER(_ctypes.c_double),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_PLATEAU_CONDITION_ARGUMENTS = ("confidence_interval", "best_candidate_pair_confidence_interval", "n_studies", "best_candidate_index", "best_exceeded_ci_overlap_count", "candidate_index", "join_method", "succeeding_ci_overlap", "ci_overlap", "plateau_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_PLATEAU_CONDITION_ARGUMENT_SOURCES = (None, None, "confidence_interval", None, None, None, None, None, None, None, None,)

_lib.check_effect_size_plateau_condition_c.restype = None
_lib.check_effect_size_plateau_condition_c.argtypes = (
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_int),
    _ctypes.POINTER(_ctypes.c_int),
    _np.ctypeslib.ndpointer(dtype=_np.float64, ndim=1, flags='C_CONTIGUOUS'),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_double),
    _ctypes.POINTER(_ctypes.c_bool),
    _ctypes.POINTER(_ctypes.c_int),
)

#: The wrapped procedure's arguments, so an error can name one
_CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENTS = ("global_js_divergence", "prev_global_js_divergence", "n_studies", "has_previous", "delta_median_threshold", "delta_max_threshold", "delta_epsilon", "delta_min_consecutive_transitions", "n_consecutive_ok", "delta", "delta_median", "delta_max", "plateau_found", "ierr",)
#: For a derived argument, the one the caller passed it in
_CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENT_SOURCES = (None, None, "global_js_divergence", None, None, None, None, None, None, None, None, None, None, None,)

def check_neighborhood_overlaps(
        neighborhood_range,
        min_neighbor_overlap,
):
    r"""Test whether every pair of consecutive neighborhoods overlaps by at least a minimum fraction

    The first admissibility gate of a parameter-search candidate, on the `[min_idx, max_idx]`
    neighborhood spans of
    :func:`tensor_omics.construct_neighborhoods_ranged`
    (or the adaptive `pooled_neighborhood_range`). The overlap of two consecutive spans is the
    number of positions they share divided by the length of the earlier span; a reversed span
    (`max < min`) counts as empty. Not named `test_*`: the R test harness would run every
    `test_`-prefixed export as a test case.

    Parameters
    ----------
    neighborhood_range : np.ndarray[np.int32] of shape (2, n_points,), column-major (order='F')
        For each reference point, the `[min_idx, max_idx]` neighborhood span
        The minimum valid value is `1`.
    min_neighbor_overlap : float
        Minimum fractional overlap two consecutive neighborhoods must have
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    all_have_min_neighbor_overlap : bool
        `True` if every pair of consecutive neighborhoods overlaps by at least
        `min_neighbor_overlap`

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_search::check_neighborhood_overlaps`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    neighborhood_range = _contract.integers(neighborhood_range, "neighborhood_range", ndim=2, order="F")
    min_neighbor_overlap = _contract.real(min_neighbor_overlap, "min_neighbor_overlap")

    # what the inputs already say, rather than asking for it again
    n_points = neighborhood_range.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if neighborhood_range.shape[0] != 2:
        raise _ToxArgumentError(f"'neighborhood_range' has {neighborhood_range.shape[0]} along axis 0, but its declared shape fixes it at {2}; 'neighborhood_range' is declared with shape (2, n_points)", argument="neighborhood_range")

    # outputs and work arrays, which the caller never sees
    all_have_min_neighbor_overlap = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.check_neighborhood_overlaps_c(
        neighborhood_range,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_double(min_neighbor_overlap)),
        _ctypes.byref(all_have_min_neighbor_overlap),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENTS, _CHECK_NEIGHBORHOOD_OVERLAPS_ARGUMENT_SOURCES)

    return all_have_min_neighbor_overlap.value

def check_mean_pmf_min_counts(
        mean_pmf_counts,
        n_bins_per_point,
        min_count,
):
    r"""Test whether every bin of a mean pmf reaches a minimum absolute count

    The second admissibility gate of a parameter-search candidate, checked after
    :func:`tensor_omics.check_neighborhood_overlaps`.
    Each point is checked over its own `n_bins_per_point` bins only, so the zero padding beyond
    them does not count as a failure.

    Parameters
    ----------
    mean_pmf_counts : np.ndarray[np.int32] of shape (n_points, n_bins,), column-major (order='F')
        Absolute counts of a residual per bin for the mean pmf
        The minimum valid value is `0`.
    n_bins_per_point : np.ndarray[np.int32] of shape (n_points,)
        Each point's own bin count; only `mean_pmf_counts(i_point, 1:n_bins_per_point(i_point))`
        is inspected
        The minimum valid value is `1`.
        The maximum valid value is `n_bins`.
    min_count : int
        Minimum count each bin of the mean pmf must reach
        The minimum valid value is `0`.

    Returns
    -------
    all_bins_have_min_count : bool
        `True` if every bin within each reference point's own `n_bins_per_point`, at every
        reference point, reaches at least `min_count`

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_search::check_mean_pmf_min_counts`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    mean_pmf_counts = _contract.integers(mean_pmf_counts, "mean_pmf_counts", ndim=2, order="F")
    n_bins_per_point = _contract.integers(n_bins_per_point, "n_bins_per_point", ndim=1)
    min_count = _contract.integer(min_count, "min_count")

    # what the inputs already say, rather than asking for it again
    n_bins = mean_pmf_counts.shape[1]
    n_points = mean_pmf_counts.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if n_bins_per_point.shape[0] != n_points:
        raise _ToxArgumentError(f"'n_bins_per_point' has {n_bins_per_point.shape[0]} along axis 0, but 'mean_pmf_counts' implies n_points == {n_points}; 'n_bins_per_point' is declared with shape (n_points)", argument="n_bins_per_point")

    # outputs and work arrays, which the caller never sees
    all_bins_have_min_count = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.check_mean_pmf_min_counts_c(
        mean_pmf_counts,
        _ctypes.byref(_ctypes.c_int(n_bins)),
        n_bins_per_point,
        _ctypes.byref(_ctypes.c_int(n_points)),
        _ctypes.byref(_ctypes.c_int(min_count)),
        _ctypes.byref(all_bins_have_min_count),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENTS, _CHECK_MEAN_PMF_MIN_COUNTS_ARGUMENT_SOURCES)

    return all_bins_have_min_count.value

def check_plateau_condition(
        confidence_interval,
        best_candidate_pair_confidence_interval,
        best_candidate_index,
        best_exceeded_ci_overlap_count,
        candidate_index,
        join_method,
        succeeding_ci_overlap,
):
    r"""Test one candidate pair's bootstrapped confidence intervals against the current best, and detect a JSD plateau

    For each study, the overlap of the new candidate's interval with the best one so far is
    measured as a fraction of the new interval's width, and counted when it reaches
    `succeeding_ci_overlap`. A plateau is found when fewer studies pass than for the current
    best (the best is then kept), or, after the new candidate became the best, when the count
    meets the condition `join_method` names.

    Parameters
    ----------
    confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F')
        Bootstrapped JSD confidence interval `[lower, upper]` of the candidate under test
    best_candidate_pair_confidence_interval : np.ndarray[np.float64] of shape (2, n_studies,), column-major (order='F'), modified in place
        JSD confidence intervals for the current best candidate pair; overwritten with
        `confidence_interval` unless the new candidate is worse
    best_candidate_index : int
        Candidate-grid index of the current best candidate pair; overwritten with
        `candidate_index` unless the new candidate is worse
    best_exceeded_ci_overlap_count : int
        Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
        best candidate pair; overwritten unless the new candidate is worse
        The minimum valid value is `0`.
    candidate_index : int
        Candidate-grid index of the candidate pair that produced `confidence_interval`
        The minimum valid value is `1`.
    join_method : str, one of 'join_min' | 'join_max' | 'join_median'
        How the studies' overlaps are joined into the plateau condition; the median method
        requires `count > (n_studies - 1) / 2`
        Trailing blanks are not part of the value, as in Fortran: they are dropped.
    succeeding_ci_overlap : float
        Minimum fractional overlap with the best candidate's interval for a study to pass
        The minimum valid value is `0.0`.
        The maximum valid value is `1.0`.

    Returns
    -------
    dict
        with keys:

        best_candidate_index : int
            Candidate-grid index of the current best candidate pair; overwritten with
            `candidate_index` unless the new candidate is worse
        best_exceeded_ci_overlap_count : int
            Number of studies whose overlap exceeded `succeeding_ci_overlap` for the current
            best candidate pair; overwritten unless the new candidate is worse
            The minimum valid value is `0`.
        ci_overlap : np.ndarray[np.float64] of shape (n_studies,), read-only
            Per-study overlap, in `[0, 1]`, of `confidence_interval` with the best interval as
            it was on entry, relative to the new interval's width; always filled. Meaningless
            for the very first candidate, whose best interval is only a seed
            A result is a value; call `.copy()` to obtain a modifiable array.
        plateau_found : bool
            `True` once the new candidate is no better than the previous best, or once
            `join_method`'s overlap condition is met by the new candidate

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_search::check_plateau_condition`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    confidence_interval = _contract.reals(confidence_interval, "confidence_interval", ndim=2, order="F")
    if not isinstance(best_candidate_pair_confidence_interval, _np.ndarray) or best_candidate_pair_confidence_interval.dtype != _np.float64:
        raise _ToxArgumentError("'best_candidate_pair_confidence_interval' is modified in place, so it must already be a numpy array of {}".format(_np.float64), argument="best_candidate_pair_confidence_interval")
    if best_candidate_pair_confidence_interval.ndim != 2:
        raise _ToxArgumentError(f"'best_candidate_pair_confidence_interval' must have 2 dimensions, but has {best_candidate_pair_confidence_interval.ndim}", argument="best_candidate_pair_confidence_interval")
    if not best_candidate_pair_confidence_interval.flags.f_contiguous:
        raise _ToxArgumentError("'best_candidate_pair_confidence_interval' is modified in place, so it must already be column-major (order='F')", argument="best_candidate_pair_confidence_interval")
    best_candidate_index = _ctypes.c_int(_contract.integer(best_candidate_index, "best_candidate_index"))
    best_exceeded_ci_overlap_count = _ctypes.c_int(_contract.integer(best_exceeded_ci_overlap_count, "best_exceeded_ci_overlap_count"))
    candidate_index = _contract.integer(candidate_index, "candidate_index")
    join_method = _np.array([_contract.mode(join_method, "join_method", ("join_min", "join_max", "join_median",)).ljust(11)], dtype="S11")
    succeeding_ci_overlap = _contract.real(succeeding_ci_overlap, "succeeding_ci_overlap")

    # what the inputs already say, rather than asking for it again
    n_studies = confidence_interval.shape[1]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if confidence_interval.shape[0] != 2:
        raise _ToxArgumentError(f"'confidence_interval' has {confidence_interval.shape[0]} along axis 0, but its declared shape fixes it at {2}; 'confidence_interval' is declared with shape (2, n_studies)", argument="confidence_interval")
    if best_candidate_pair_confidence_interval.shape[0] != 2:
        raise _ToxArgumentError(f"'best_candidate_pair_confidence_interval' has {best_candidate_pair_confidence_interval.shape[0]} along axis 0, but its declared shape fixes it at {2}; 'best_candidate_pair_confidence_interval' is declared with shape (2, n_studies)", argument="best_candidate_pair_confidence_interval")
    if best_candidate_pair_confidence_interval.shape[1] != n_studies:
        raise _ToxArgumentError(f"'best_candidate_pair_confidence_interval' has {best_candidate_pair_confidence_interval.shape[1]} along axis 1, but 'confidence_interval' implies n_studies == {n_studies}; 'best_candidate_pair_confidence_interval' is declared with shape (2, n_studies)", argument="best_candidate_pair_confidence_interval")

    # outputs and work arrays, which the caller never sees
    ci_overlap = _np.empty((n_studies,), dtype=_np.float64, order='C')
    plateau_found = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.check_plateau_condition_c(
        confidence_interval,
        best_candidate_pair_confidence_interval,
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(best_candidate_index),
        _ctypes.byref(best_exceeded_ci_overlap_count),
        _ctypes.byref(_ctypes.c_int(candidate_index)),
        join_method,
        _ctypes.byref(_ctypes.c_double(succeeding_ci_overlap)),
        ci_overlap,
        _ctypes.byref(plateau_found),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CHECK_PLATEAU_CONDITION_ARGUMENTS, _CHECK_PLATEAU_CONDITION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    ci_overlap.flags.writeable = False

    return {
        "best_candidate_index": best_candidate_index.value,
        "best_exceeded_ci_overlap_count": best_exceeded_ci_overlap_count.value,
        "ci_overlap": ci_overlap,
        "plateau_found": plateau_found.value,
    }

def check_effect_size_plateau_condition(
        global_js_divergence,
        prev_global_js_divergence,
        has_previous,
        delta_median_threshold,
        delta_max_threshold,
        delta_epsilon,
        delta_min_consecutive_transitions,
        n_consecutive_ok,
):
    r"""Test one candidate's per-study JSD against the previous admissible candidate for a relative-effect-size plateau

    The relative-effect-size alternative to
    :func:`tensor_omics.check_plateau_condition`:
    for each study `i`, the relative change of the observed global JSD between successive
    admissible candidates,
    `delta(i) = |global_js_divergence(i) - prev_global_js_divergence(i)| / max(prev_global_js_divergence(i),
    delta_epsilon)`, is summarized by its median (`delta_median`) and maximum (`delta_max`). A
    plateau is found once both stay below their thresholds for
    `delta_min_consecutive_transitions` consecutive transitions, counted across calls in
    `n_consecutive_ok`.

    For the first admissible candidate (`has_previous = False`) there is no transition:
    `delta`, `delta_median` and `delta_max` are `-1.0` ("not computed").

    Parameters
    ----------
    global_js_divergence : np.ndarray[np.float64] of shape (n_studies,)
        Current admissible candidate's observed global JSD per study
        The minimum valid value is `0.0`.
    prev_global_js_divergence : np.ndarray[np.float64] of shape (n_studies,)
        Previous admissible candidate's observed global JSD per study; ignored when
        `has_previous` is `False`
        The minimum valid value is `0.0`.
    has_previous : bool
        `False` for the first admissible candidate, which has no transition
    delta_median_threshold : float
        Upper bound on the median relative change for a transition to qualify
        The minimum valid value is `above(0.0)`.
    delta_max_threshold : float
        Upper bound on the largest relative change for a transition to qualify
        The minimum valid value is `above(0.0)`.
    delta_epsilon : float
        Small constant preventing division by zero when a study's previous JSD was zero
        The minimum valid value is `above(0.0)`.
    delta_min_consecutive_transitions : int
        Number of consecutive qualifying transitions required to declare a plateau
        The minimum valid value is `1`.
    n_consecutive_ok : int
        Running count of consecutive qualifying transitions; reset to zero when this one does
        not qualify or `has_previous` is `False`
        The minimum valid value is `0`.

    Returns
    -------
    dict
        with keys:

        n_consecutive_ok : int
            Running count of consecutive qualifying transitions; reset to zero when this one does
            not qualify or `has_previous` is `False`
            The minimum valid value is `0`.
        delta : np.ndarray[np.float64] of shape (n_studies,), read-only
            Per-study relative JSD change from the previous admissible candidate; `-1.0`
            throughout iff `.not. has_previous`
            A result is a value; call `.copy()` to obtain a modifiable array.
        delta_median : float
            Median of `delta` across studies; `-1.0` iff `.not. has_previous`
        delta_max : float
            Maximum of `delta` across studies; `-1.0` iff `.not. has_previous`
        plateau_found : bool
            `True` once `n_consecutive_ok` reaches `delta_min_consecutive_transitions`

    Raises
    ------
    ToxArgumentError
        If an argument is refused before the library is called: its shape, type
        or value. A ToxInputError without a code, and a TypeError and a ValueError.
    ToxError
        If the underlying Fortran reports an error.

    Notes
    -----
    Generated from the Fortran procedure `tox_data_integration_js_comp_test_search::check_effect_size_plateau_condition`, whose argument names are
    the ones an error message reports.
    """
    # check every input by the input contract, converting only when C needs it
    global_js_divergence = _contract.reals(global_js_divergence, "global_js_divergence", ndim=1)
    prev_global_js_divergence = _contract.reals(prev_global_js_divergence, "prev_global_js_divergence", ndim=1)
    has_previous = _contract.logical(has_previous, "has_previous")
    delta_median_threshold = _contract.real(delta_median_threshold, "delta_median_threshold")
    delta_max_threshold = _contract.real(delta_max_threshold, "delta_max_threshold")
    delta_epsilon = _contract.real(delta_epsilon, "delta_epsilon")
    delta_min_consecutive_transitions = _contract.integer(delta_min_consecutive_transitions, "delta_min_consecutive_transitions")
    n_consecutive_ok = _ctypes.c_int(_contract.integer(n_consecutive_ok, "n_consecutive_ok"))

    # what the inputs already say, rather than asking for it again
    n_studies = global_js_divergence.shape[0]

    # Fortran cannot check that an array has the shape it is declared with; this can
    if prev_global_js_divergence.shape[0] != n_studies:
        raise _ToxArgumentError(f"'prev_global_js_divergence' has {prev_global_js_divergence.shape[0]} along axis 0, but 'global_js_divergence' implies n_studies == {n_studies}; 'prev_global_js_divergence' is declared with shape (n_studies)", argument="prev_global_js_divergence")

    # outputs and work arrays, which the caller never sees
    delta = _np.empty((n_studies,), dtype=_np.float64, order='C')
    delta_median = _ctypes.c_double(0)
    delta_max = _ctypes.c_double(0)
    plateau_found = _ctypes.c_bool(0)
    ierr = _ctypes.c_int(0)

    _lib.check_effect_size_plateau_condition_c(
        global_js_divergence,
        prev_global_js_divergence,
        _ctypes.byref(_ctypes.c_int(n_studies)),
        _ctypes.byref(_ctypes.c_bool(has_previous)),
        _ctypes.byref(_ctypes.c_double(delta_median_threshold)),
        _ctypes.byref(_ctypes.c_double(delta_max_threshold)),
        _ctypes.byref(_ctypes.c_double(delta_epsilon)),
        _ctypes.byref(_ctypes.c_int(delta_min_consecutive_transitions)),
        _ctypes.byref(n_consecutive_ok),
        delta,
        _ctypes.byref(delta_median),
        _ctypes.byref(delta_max),
        _ctypes.byref(plateau_found),
        _ctypes.byref(ierr),
    )

    _check_err_code(ierr.value, _CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENTS, _CHECK_EFFECT_SIZE_PLATEAU_CONDITION_ARGUMENT_SOURCES)

    # a result is a value: modify a copy, not this
    delta.flags.writeable = False

    return {
        "n_consecutive_ok": n_consecutive_ok.value,
        "delta": delta,
        "delta_median": delta_median.value,
        "delta_max": delta_max.value,
        "plateau_found": plateau_found.value,
    }
